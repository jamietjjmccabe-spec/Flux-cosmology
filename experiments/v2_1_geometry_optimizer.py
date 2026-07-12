#!/usr/bin/env python3
"""
v2_1_geometry_optimizer.py

Prospective geometry optimizer for frozen Flux-MIP Benchmark V2.1.

The benchmark apparatus is deliberately explicit:

* Source: finite rectangular tungsten plate, density 19,300 kg m^-3.
* Plate dimensions parallel to the atomic trajectories:
    x width = 0.16 m
    z height = 0.30 m
* Two straight, parallel 87Rb trajectories run along z for 1.0 s through
  the central 0.20 m of the plate height.
* Arm 1 is at y=0; arm 2 is at y=0.040 m.
* The plate lies between the arms. Its near face is y_gap from arm 1 and
  its thickness is d_plate. A finite clearance is required to arm 2.
* The optimizer sweeps y_gap and d_plate and ranks candidates by

      X = 1/2 integral [F_1(t)-F_2(t)]^2 dt.

This is a geometry-sensitivity benchmark, not yet a full engineering model.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.constants import hbar, physical_constants
from scipy.integrate import trapezoid

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from flux_lindblad_mip_v2_1 import BenchmarkV21, MassVoxels, evaluate_points  # noqa: E402

TUNGSTEN_DENSITY_KG_M3 = 19300.0
RB87_MASS_KG = 86.9091805310 * physical_constants["atomic mass constant"][0]


@dataclass(frozen=True)
class OptimizerSpec:
    arm_separation_m: float = 0.040
    interrogation_time_s: float = 1.0
    trajectory_length_m: float = 0.20
    trajectory_samples: int = 201
    plate_width_x_m: float = 0.16
    plate_height_z_m: float = 0.30
    tungsten_density_kg_m3: float = TUNGSTEN_DENSITY_KG_M3
    xz_voxel_step_m: float = 0.004
    y_layer_step_m: float = 0.001
    minimum_surface_clearance_m: float = 0.001
    minimum_plate_thickness_m: float = 0.001
    maximum_plate_thickness_m: float = 0.030
    ranking_metric: str = "X_half_dephasing_exposure_s"
    version: str = "Flux-MIP-V2.1-Dedicated-Plate-Optimizer-Benchmark-1"

    def signature(self, benchmark_signature: str) -> str:
        payload = {"optimizer": asdict(self), "benchmark_signature_sha256": benchmark_signature}
        text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode("utf-8")).hexdigest()


def cell_centres(length_m: float, step_m: float) -> np.ndarray:
    n = max(1, int(round(length_m / step_m)))
    actual_step = length_m / n
    return -0.5 * length_m + (np.arange(n) + 0.5) * actual_step


def make_trajectories(spec: OptimizerSpec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = np.linspace(0.0, spec.interrogation_time_s, spec.trajectory_samples)
    z = np.linspace(-0.5 * spec.trajectory_length_m, 0.5 * spec.trajectory_length_m, spec.trajectory_samples)
    x1 = np.column_stack((np.zeros_like(t), np.zeros_like(t), z))
    x2 = np.column_stack((np.zeros_like(t), np.full_like(t, spec.arm_separation_m), z))
    return t, x1, x2


def layer_source(y_center_m: float, spec: OptimizerSpec, benchmark: BenchmarkV21) -> MassVoxels:
    xs = cell_centres(spec.plate_width_x_m, spec.xz_voxel_step_m)
    zs = cell_centres(spec.plate_height_z_m, spec.xz_voxel_step_m)
    dx = spec.plate_width_x_m / len(xs)
    dz = spec.plate_height_z_m / len(zs)
    dy = spec.y_layer_step_m
    X, Z = np.meshgrid(xs, zs, indexing="ij")
    positions = np.column_stack((X.ravel(), np.full(X.size, y_center_m), Z.ravel()))
    voxel_volume = dx * dy * dz
    voxel_mass = spec.tungsten_density_kg_m3 * voxel_volume
    normalization_mass = benchmark.rho_star_kg_m3 * benchmark.reference_volume_m3
    weight = benchmark.kappa_reg * benchmark.coherence_factor * voxel_mass / normalization_mass
    weights = np.full(positions.shape[0], weight, dtype=float)
    return MassVoxels(
        positions_m=positions,
        weights=weights,
        voxel_volume_m3=voxel_volume,
        represented_mass_kg=float(voxel_mass * positions.shape[0]),
        config=benchmark,
    )


def precompute_layers(spec: OptimizerSpec, benchmark: BenchmarkV21, x1: np.ndarray, x2: np.ndarray,
                      point_chunk: int, source_chunk: int) -> dict[str, np.ndarray]:
    dy = spec.y_layer_step_m
    n_layers = int(math.floor(spec.arm_separation_m / dy))
    y_centres = (np.arange(n_layers) + 0.5) * dy
    points = np.vstack((x1, x2))
    S_layers = np.empty((n_layers, points.shape[0]), dtype=float)
    D_layers = np.empty_like(S_layers)
    gradS_layers = np.empty((n_layers, points.shape[0], 3), dtype=float)
    for i, yc in enumerate(y_centres):
        field = evaluate_points(
            layer_source(float(yc), spec, benchmark), points,
            point_chunk=point_chunk, source_chunk=source_chunk,
        )
        D_layers[i] = field.D_reg
        S_layers[i] = field.S_reg
        gradS_layers[i] = field.grad_S_per_m
        print(f"Precomputed layer {i+1:02d}/{n_layers}: y={yc*1e3:.1f} mm", flush=True)
    return {"y_centres_m": y_centres, "D_layers": D_layers,
            "S_layers": S_layers, "gradS_layers": gradS_layers}


def candidate_metrics(start: int, count: int, pre: dict[str, np.ndarray], t: np.ndarray, n_path: int,
                      spec: OptimizerSpec, benchmark: BenchmarkV21,
                      gamma_reference_s_inv: float = 1.0) -> dict[str, Any]:
    stop = start + count
    D = np.sum(pre["D_layers"][start:stop], axis=0)
    S = np.sum(pre["S_layers"][start:stop], axis=0)
    gradS = np.sum(pre["gradS_layers"][start:stop], axis=0)
    F = np.tanh(S / benchmark.susceptibility_star)
    gradF = ((1.0 - F**2) / benchmark.susceptibility_star)[:, None] * gradS
    D1, D2 = D[:n_path], D[n_path:]
    S1, S2 = S[:n_path], S[n_path:]
    F1, F2 = F[:n_path], F[n_path:]
    g1, g2 = gradF[:n_path], gradF[n_path:]
    deltaF = F1 - F2
    I_phase = float(trapezoid(deltaF, t))
    I_deph = float(trapezoid(deltaF**2, t))
    X = 0.5 * I_deph
    J = np.empty((3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            integrand = 0.5 * (g1[:, i]*g1[:, j] + g2[:, i]*g2[:, j])
            J[i, j] = trapezoid(integrand, t)
    Jtrace = float(np.trace(J))
    delta_E_gamma1 = gamma_reference_s_inv * hbar**2 * Jtrace / (2.0 * RB87_MASS_KG)
    visibility_loss_gamma1 = 1.0 - math.exp(-gamma_reference_s_inv * X)
    dy = spec.y_layer_step_m
    gap = start * dy
    thickness = count * dy
    far_gap = spec.arm_separation_m - gap - thickness
    plate_mass = spec.tungsten_density_kg_m3 * spec.plate_width_x_m * spec.plate_height_z_m * thickness
    return {
        "y_gap_m": gap, "plate_thickness_m": thickness,
        "far_arm_clearance_m": far_gap, "plate_mass_kg": plate_mass,
        "phase_exposure_integral_s": I_phase,
        "dephasing_exposure_integral_s": I_deph,
        "X_half_dephasing_exposure_s": X,
        "gradient_heating_exposure_s_per_m2": Jtrace,
        "visibility_loss_at_gamma_1": visibility_loss_gamma1,
        "delta_E_Rb87_at_gamma_1_J": delta_E_gamma1,
        "F1_min": float(F1.min()), "F1_max": float(F1.max()), "F1_mean": float(F1.mean()),
        "F2_min": float(F2.min()), "F2_max": float(F2.max()), "F2_mean": float(F2.mean()),
        "delta_F_rms": float(np.sqrt(np.mean(deltaF**2))),
        "S1_mean": float(S1.mean()), "S2_mean": float(S2.mean()),
        "D1_mean": float(D1.mean()), "D2_mean": float(D2.mean()),
        "start_layer": int(start), "layer_count": int(count),
    }


def optimize(spec: OptimizerSpec, benchmark: BenchmarkV21, point_chunk: int = 128,
             source_chunk: int = 20000) -> tuple[pd.DataFrame, dict[str, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    t, x1, x2 = make_trajectories(spec)
    pre = precompute_layers(spec, benchmark, x1, x2, point_chunk, source_chunk)
    dy = spec.y_layer_step_m
    n_layers = len(pre["y_centres_m"])
    min_clear_layers = int(math.ceil(spec.minimum_surface_clearance_m / dy))
    min_thick_layers = int(math.ceil(spec.minimum_plate_thickness_m / dy))
    max_thick_layers = int(math.floor(spec.maximum_plate_thickness_m / dy))
    rows: list[dict[str, Any]] = []
    for count in range(min_thick_layers, max_thick_layers + 1):
        for start in range(min_clear_layers, n_layers - min_clear_layers - count + 1):
            rows.append(candidate_metrics(start, count, pre, t, len(t), spec, benchmark))
    df = pd.DataFrame(rows).sort_values(
        ["X_half_dephasing_exposure_s", "gradient_heating_exposure_s_per_m2"],
        ascending=[False, True],
    ).reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df)+1))
    return df, pre, (t, x1, x2)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-prefix", type=Path, default=Path("v2_1_plate_optimizer"))
    p.add_argument("--top", type=int, default=25)
    p.add_argument("--point-chunk", type=int, default=128)
    p.add_argument("--source-chunk", type=int, default=20000)
    args = p.parse_args()
    benchmark = BenchmarkV21()
    spec = OptimizerSpec()
    df, pre, paths = optimize(spec, benchmark, args.point_chunk, args.source_chunk)
    top = df.head(args.top).copy()
    best = df.iloc[0].to_dict()
    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    all_csv = prefix.with_name(prefix.name + "_all.csv")
    top_csv = prefix.with_name(prefix.name + "_top.csv")
    summary_json = prefix.with_name(prefix.name + "_summary.json")
    pre_npz = prefix.with_name(prefix.name + "_precomputed_layers.npz")
    trajectories_npz = prefix.with_name(prefix.name + "_trajectories.npz")
    df.to_csv(all_csv, index=False)
    top.to_csv(top_csv, index=False)
    np.savez_compressed(pre_npz, **pre)
    t, x1, x2 = paths
    np.savez_compressed(trajectories_npz, t=t, x1=x1, x2=x2)
    summary = {
        "benchmark": {**asdict(benchmark), "reference_volume_m3": benchmark.reference_volume_m3,
                      "signature_sha256": benchmark.signature()},
        "optimizer_specification": asdict(spec),
        "optimizer_signature_sha256": spec.signature(benchmark.signature()),
        "candidate_count": int(len(df)), "best_candidate": best,
        "top_candidates": top.to_dict(orient="records"),
    }
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(top.head(10).to_string(index=False))
    print(f"V2.1 benchmark: {benchmark.signature()}")
    print(f"optimizer: {spec.signature(benchmark.signature())}")


if __name__ == "__main__":
    main()
