#!/usr/bin/env python3
"""
flux_lindblad_mip_v2_1.py

Flux-MIP Benchmark V2.1: regularized inverse-square precursor with a fixed
curvature high-pass susceptibility rule.

The absolute registration precursor is retained from V2:

    G_ell(r) = ell^2 / (r^2 + ell^2)

    D_reg(x) = kappa_reg/(rho_* V_*) * integral rho_b(r') G_ell(|x-r'|) d^3r'

The laboratory Lindblad operator does NOT couple to the absolute D_reg.  It
couples to the local curvature/high-pass response

    S_reg(x) = -ell^2 Laplacian[D_reg(x)]

and the bounded Hermitian response field

    F_reg(x) = tanh(S_reg(x)/S_*)

Frozen V2.1 values:

    ell       = 1.0 cm
    p         = 2 (fixed)
    rho_*     = 1000 kg m^-3
    V_*       = (4 pi/3) ell^3
    S_*       = 10
    kappa_reg = 1
    C         = 1

Why this is the background rule
-------------------------------
* Constant D backgrounds vanish under the Laplacian.
* Locally affine backgrounds (constant plus uniform gradient) also vanish.
* Smooth terrestrial fields enter only through curvature/tidal structure.
* For G_ell ~ r^-2, -ell^2 Laplacian(G_ell) ~ r^-4, so the laboratory
  susceptibility is infrared integrable in three spatial dimensions.
* No adaptive D_* and no post-data environmental rescaling are allowed.

For p=2, define a = ell^2, s = r^2+a.  The analytic high-pass kernel is

    H_ell(r) = -ell^2 Laplacian G_ell(r)
             = 2 a^2 (3a-r^2) / (r^2+a)^3

and

    grad H_ell = 8 a^2 (r^2-5a) (x-r') / (r^2+a)^4.

The Lindblad observables are

    I_phase = integral [F1-F2] dt
    X       = 1/2 integral [F1-F2]^2 dt
    J_ij    = 1/2 integral [d_iF1 d_jF1 + d_iF2 d_jF2] dt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.constants import hbar
from scipy.integrate import trapezoid


@dataclass(frozen=True)
class BenchmarkV21:
    ell_core_m: float = 1.0e-2
    power_p: float = 2.0
    rho_star_kg_m3: float = 1000.0
    reference_volume_factor: float = 4.1887902047863905
    susceptibility_star: float = 10.0
    kappa_reg: float = 1.0
    coherence_factor: float = 1.0
    precursor_kernel: str = "G(r)=ell^2/(r^2+ell^2)"
    background_rule: str = "S=-ell^2 Laplacian(D); F=tanh(S/S_star)"
    temporal_high_pass: bool = False
    adaptive_threshold: bool = False
    version: str = "Flux-MIP-Benchmark-V2.1-RPL-p2-curvature-HP"

    @property
    def reference_volume_m3(self) -> float:
        return self.reference_volume_factor * self.ell_core_m**3

    def signature(self) -> str:
        payload = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass
class MassVoxels:
    positions_m: np.ndarray
    weights: np.ndarray
    voxel_volume_m3: float
    represented_mass_kg: float
    config: BenchmarkV21


@dataclass
class PointField:
    D_reg: np.ndarray
    S_reg: np.ndarray
    F_reg: np.ndarray
    grad_D_per_m: np.ndarray
    grad_S_per_m: np.ndarray
    grad_F_per_m: np.ndarray


def _uniform_spacing(axis: np.ndarray, name: str) -> float:
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 1 or axis.size < 2:
        raise ValueError(f"{name} must be one-dimensional with at least two samples.")
    delta = np.diff(axis)
    if np.any(delta <= 0):
        raise ValueError(f"{name} must be strictly increasing.")
    mean = float(np.mean(delta))
    if not np.allclose(delta, mean, rtol=1e-8, atol=max(1e-14, abs(mean)*1e-10)):
        raise ValueError(f"{name} must be uniformly spaced.")
    return mean


def mass_grid_to_voxels(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    *,
    config: BenchmarkV21 = BenchmarkV21(),
    density_threshold_kg_m3: float = 0.0,
) -> MassVoxels:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    rho = np.asarray(rho, dtype=float)

    dx = _uniform_spacing(x, "x")
    dy = _uniform_spacing(y, "y")
    dz = _uniform_spacing(z, "z")
    if rho.shape != (x.size, y.size, z.size):
        raise ValueError(f"rho shape {rho.shape} does not match {(x.size, y.size, z.size)}.")
    if not np.all(np.isfinite(rho)) or np.any(rho < 0):
        raise ValueError("rho must be finite and non-negative.")

    mask = rho > density_threshold_kg_m3
    indices = np.argwhere(mask)
    if indices.size == 0:
        raise ValueError("No source voxels above threshold.")

    positions = np.column_stack((x[indices[:, 0]], y[indices[:, 1]], z[indices[:, 2]]))
    rho_values = rho[mask]
    dV = dx * dy * dz
    masses = rho_values * dV
    normalization_mass = config.rho_star_kg_m3 * config.reference_volume_m3
    weights = config.kappa_reg * config.coherence_factor * masses / normalization_mass

    return MassVoxels(
        positions_m=positions,
        weights=weights,
        voxel_volume_m3=dV,
        represented_mass_kg=float(np.sum(masses)),
        config=config,
    )


def evaluate_points(
    source: MassVoxels,
    points_m: np.ndarray,
    *,
    point_chunk: int = 128,
    source_chunk: int = 20000,
) -> PointField:
    points = np.asarray(points_m, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points_m must have shape (N,3).")
    if not np.all(np.isfinite(points)):
        raise ValueError("points_m contains non-finite values.")

    cfg = source.config
    if cfg.power_p != 2.0:
        raise ValueError("V2.1 is frozen to p=2.")

    a = cfg.ell_core_m**2
    n = points.shape[0]
    D = np.zeros(n, dtype=float)
    S = np.zeros(n, dtype=float)
    grad_D = np.zeros((n, 3), dtype=float)
    grad_S = np.zeros((n, 3), dtype=float)

    src_pos = source.positions_m
    src_w = source.weights

    for p0 in range(0, n, point_chunk):
        p1 = min(n, p0 + point_chunk)
        pts = points[p0:p1]
        D_part = np.zeros(p1-p0, dtype=float)
        S_part = np.zeros(p1-p0, dtype=float)
        gradD_part = np.zeros((p1-p0, 3), dtype=float)
        gradS_part = np.zeros((p1-p0, 3), dtype=float)

        for s0 in range(0, src_pos.shape[0], source_chunk):
            s1 = min(src_pos.shape[0], s0 + source_chunk)
            pos = src_pos[s0:s1]
            w = src_w[s0:s1]

            diff = pts[:, None, :] - pos[None, :, :]
            r2 = np.einsum("nmk,nmk->nm", diff, diff)
            den = r2 + a

            G = a / den
            D_part += G @ w
            coeff_grad_G = -2.0 * a / den**2
            gradD_part += np.einsum("nm,nmk->nk", coeff_grad_G * w[None, :], diff)

            H = 2.0 * a**2 * (3.0*a - r2) / den**3
            S_part += H @ w
            coeff_grad_H = 8.0 * a**2 * (r2 - 5.0*a) / den**4
            gradS_part += np.einsum("nm,nmk->nk", coeff_grad_H * w[None, :], diff)

        D[p0:p1] = D_part
        S[p0:p1] = S_part
        grad_D[p0:p1] = gradD_part
        grad_S[p0:p1] = gradS_part

    scaled = S / cfg.susceptibility_star
    F = np.tanh(scaled)
    dF_dS = (1.0 - F**2) / cfg.susceptibility_star
    grad_F = dF_dS[:, None] * grad_S

    return PointField(
        D_reg=D,
        S_reg=S,
        F_reg=F,
        grad_D_per_m=grad_D,
        grad_S_per_m=grad_S,
        grad_F_per_m=grad_F,
    )


def analyze_trajectories_v21(
    source: MassVoxels,
    t: np.ndarray,
    x1: np.ndarray,
    x2: np.ndarray,
    *,
    gamma_mip_s_inv: float = 0.0,
    beta_mip_s_inv: float = 0.0,
    atom_mass_kg: float = 1.44316060e-25,
    point_chunk: int = 128,
    source_chunk: int = 20000,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    t = np.asarray(t, dtype=float)
    x1 = np.asarray(x1, dtype=float)
    x2 = np.asarray(x2, dtype=float)
    if t.ndim != 1 or t.size < 2 or np.any(np.diff(t) <= 0):
        raise ValueError("t must be strictly increasing with at least two samples.")
    if x1.shape != (t.size, 3) or x2.shape != (t.size, 3):
        raise ValueError("x1 and x2 must have shape (len(t),3).")
    if gamma_mip_s_inv < 0:
        raise ValueError("gamma_mip_s_inv must be non-negative.")
    if atom_mass_kg <= 0:
        raise ValueError("atom_mass_kg must be positive.")

    f1 = evaluate_points(source, x1, point_chunk=point_chunk, source_chunk=source_chunk)
    f2 = evaluate_points(source, x2, point_chunk=point_chunk, source_chunk=source_chunk)

    delta_F = f1.F_reg - f2.F_reg
    I_phase = float(trapezoid(delta_F, t))
    I_deph = float(trapezoid(delta_F**2, t))
    X = 0.5 * I_deph

    J = np.empty((3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            integrand = 0.5 * (
                f1.grad_F_per_m[:, i] * f1.grad_F_per_m[:, j]
                + f2.grad_F_per_m[:, i] * f2.grad_F_per_m[:, j]
            )
            J[i, j] = trapezoid(integrand, t)

    delta_p_cov = gamma_mip_s_inv * hbar**2 * J
    delta_E = float(np.trace(delta_p_cov) / (2.0 * atom_mass_kg))
    cfg = source.config

    results: dict[str, Any] = {
        "benchmark_version": cfg.version,
        "benchmark_signature_sha256": cfg.signature(),
        "configuration": {**asdict(cfg), "reference_volume_m3": cfg.reference_volume_m3},
        "represented_source_mass_kg": source.represented_mass_kg,
        "source_voxel_count": int(source.positions_m.shape[0]),
        "trajectory_duration_s": float(t[-1]-t[0]),
        "gamma_mip_s_inv": float(gamma_mip_s_inv),
        "beta_mip_s_inv": float(beta_mip_s_inv),
        "atom_mass_kg": float(atom_mass_kg),
        "phase_exposure_integral_s": I_phase,
        "dephasing_exposure_integral_s": I_deph,
        "X_half_dephasing_exposure_s": X,
        "predicted_mip_phase_shift_rad": float(beta_mip_s_inv * I_phase),
        "predicted_visibility_ratio": float(math.exp(-gamma_mip_s_inv * X)),
        "diffusion_tensor_exposure_s_per_m2": J.tolist(),
        "gradient_heating_exposure_s_per_m2": float(np.trace(J)),
        "predicted_delta_p_covariance_SI": delta_p_cov.tolist(),
        "predicted_delta_kinetic_energy_j": delta_E,
        "path_statistics": {
            "D1_min": float(f1.D_reg.min()),
            "D1_max": float(f1.D_reg.max()),
            "D1_mean": float(f1.D_reg.mean()),
            "D2_min": float(f2.D_reg.min()),
            "D2_max": float(f2.D_reg.max()),
            "D2_mean": float(f2.D_reg.mean()),
            "S1_min": float(f1.S_reg.min()),
            "S1_max": float(f1.S_reg.max()),
            "S1_mean": float(f1.S_reg.mean()),
            "S2_min": float(f2.S_reg.min()),
            "S2_max": float(f2.S_reg.max()),
            "S2_mean": float(f2.S_reg.mean()),
            "F1_min": float(f1.F_reg.min()),
            "F1_max": float(f1.F_reg.max()),
            "F1_mean": float(f1.F_reg.mean()),
            "F2_min": float(f2.F_reg.min()),
            "F2_max": float(f2.F_reg.max()),
            "F2_mean": float(f2.F_reg.mean()),
            "delta_F_rms": float(np.sqrt(np.mean(delta_F**2))),
            "gradF1_rms_per_m": float(np.sqrt(np.mean(np.sum(f1.grad_F_per_m**2, axis=1)))),
            "gradF2_rms_per_m": float(np.sqrt(np.mean(np.sum(f2.grad_F_per_m**2, axis=1)))),
        },
    }

    sampled = {
        "t": t,
        "x1": x1,
        "x2": x2,
        "D1": f1.D_reg,
        "D2": f2.D_reg,
        "S1": f1.S_reg,
        "S2": f2.S_reg,
        "F1": f1.F_reg,
        "F2": f2.F_reg,
        "delta_F": delta_F,
        "grad_D1": f1.grad_D_per_m,
        "grad_D2": f2.grad_D_per_m,
        "grad_S1": f1.grad_S_per_m,
        "grad_S2": f2.grad_S_per_m,
        "grad_F1": f1.grad_F_per_m,
        "grad_F2": f2.grad_F_per_m,
    }
    return results, sampled


def load_mass_npz(path: Path):
    with np.load(path) as data:
        return data["x"], data["y"], data["z"], data["rho"]


def load_trajectory_npz(path: Path):
    with np.load(path) as data:
        return data["t"], data["x1"], data["x2"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate frozen Flux-MIP Benchmark V2.1.")
    parser.add_argument("--mass", type=Path, required=True)
    parser.add_argument("--trajectories", type=Path, required=True)
    parser.add_argument("--gamma", type=float, default=0.0)
    parser.add_argument("--beta", type=float, default=0.0)
    parser.add_argument("--atom-mass", type=float, default=1.44316060e-25)
    parser.add_argument("--point-chunk", type=int, default=128)
    parser.add_argument("--source-chunk", type=int, default=20000)
    parser.add_argument("--output-prefix", type=Path, default=Path("flux_mip_v2_1"))
    args = parser.parse_args()

    x, y, z, rho = load_mass_npz(args.mass)
    t, x1, x2 = load_trajectory_npz(args.trajectories)
    config = BenchmarkV21()
    source = mass_grid_to_voxels(x, y, z, rho, config=config)
    results, sampled = analyze_trajectories_v21(
        source, t, x1, x2,
        gamma_mip_s_inv=args.gamma,
        beta_mip_s_inv=args.beta,
        atom_mass_kg=args.atom_mass,
        point_chunk=args.point_chunk,
        source_chunk=args.source_chunk,
    )

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    results_path = prefix.with_name(prefix.name + "_results.json")
    sampled_path = prefix.with_name(prefix.name + "_sampled_fields.npz")
    results_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    np.savez_compressed(sampled_path, **sampled)
    print(json.dumps(results, indent=2))
    print(f"\nWrote: {results_path.resolve()}")
    print(f"Wrote: {sampled_path.resolve()}")


if __name__ == "__main__":
    main()
