#!/usr/bin/env python3
"""
supercluster.py

Flux Cosmology toy model for supercluster-scale coherence sinks.

Purpose
-------
Test the idea that structures such as the Great Attractor are not single hidden
objects, but large-scale flux/coherence basins visible through galaxy peculiar
velocity convergence.

Core interpretation
-------------------
Standard language:
    galaxies have peculiar velocities caused by nearby mass overdensities.

Flux language:
    galaxy envelopes drift along gradients of a supercluster flux field.
    A high-level sink is detected as negative velocity divergence:

        S_flux = ∫ max(0, -∇·v_pec) dV

Outputs
-------
Creates ./supercluster_outputs/ with:
    - supercluster_summary.txt
    - supercluster_galaxies.csv
    - flux_sink_grid.csv
    - velocity_convergence_map.png
    - flux_potential_map.png
    - basin_assignment_map.png
    - residual_flow_map.png
    - filament_channel_map.png

Dependencies
------------
    numpy, matplotlib

Run
---
    python supercluster.py

Optional examples
-----------------
    python supercluster.py --n-galaxies 1800 --grid 180 --seed 42
    python supercluster.py --no-plots
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Data containers
# ============================================================

@dataclass(frozen=True)
class Attractor:
    """A regional supercluster sink / basin centre."""
    name: str
    x: float          # Mpc
    y: float          # Mpc
    mass: float       # arbitrary mass / sink strength scale
    flux_gain: float  # extra coherence sink strength beyond mass-only model
    core: float       # softening/core scale in Mpc


@dataclass
class SimulationConfig:
    n_galaxies: int = 1200
    box_mpc: float = 260.0
    grid_n: int = 150
    seed: int = 7

    # Flow physics
    hubble_like_background: float = 0.0  # kept off; peculiar field only
    mobility: float = 0.085             # converts gradient into peculiar velocity
    noise_sigma: float = 18.0           # km/s-like random peculiar scatter
    filament_strength: float = 0.35
    residual_flux_boost: float = 0.45   # how strongly flux deviates from mass-only

    # Numerical controls
    eps: float = 1e-9
    output_dir: str = "supercluster_outputs"
    make_plots: bool = True


# ============================================================
# Default toy universe
# ============================================================

def default_attractors() -> List[Attractor]:
    """
    Three-basin toy model.

    Great Attractor is represented as a strong local basin.
    Shapley is farther and deeper.
    Vela-like hidden basin is offset behind the zone of avoidance analogue.
    """
    return [
        Attractor("Great_Attractor", x=-35.0, y=-10.0, mass=4.5, flux_gain=1.15, core=18.0),
        Attractor("Shapley",         x= 85.0, y= 45.0, mass=8.0, flux_gain=1.35, core=28.0),
        Attractor("Vela_Basin",      x= 15.0, y=-85.0, mass=6.2, flux_gain=1.55, core=24.0),
    ]


# ============================================================
# Field definitions
# ============================================================

def softened_potential(dx: np.ndarray, dy: np.ndarray, strength: float, core: float) -> np.ndarray:
    """
    Negative softened potential.

    Phi = -strength / sqrt(r^2 + core^2)
    """
    r2 = dx * dx + dy * dy + core * core
    return -strength / np.sqrt(r2)


def potential_field(
    x: np.ndarray,
    y: np.ndarray,
    attractors: List[Attractor],
    use_flux: bool = True,
) -> np.ndarray:
    """
    Total potential field.

    mass-only mode uses attractor.mass.
    flux mode boosts sink depth by flux_gain, modelling unactualized envelope draw.
    """
    phi = np.zeros_like(x, dtype=float)
    for a in attractors:
        dx = x - a.x
        dy = y - a.y
        strength = a.mass * (a.flux_gain if use_flux else 1.0)
        phi += softened_potential(dx, dy, strength=strength, core=a.core)
    return phi


def filament_channel_field(x: np.ndarray, y: np.ndarray, attractors: List[Attractor], strength: float) -> np.ndarray:
    """
    Add flux highways between basins.

    Filaments are modelled as broad negative channels along lines connecting
    attractors. This is not fitted cosmology; it is a toy test for whether
    residual flow aligns with large-scale channels.
    """
    channel = np.zeros_like(x, dtype=float)
    pairs = [(0, 1), (0, 2), (1, 2)]

    for i, j in pairs:
        a = attractors[i]
        b = attractors[j]

        ax, ay = a.x, a.y
        bx, by = b.x, b.y
        vx, vy = bx - ax, by - ay
        length2 = vx * vx + vy * vy + 1e-12

        # Projection parameter along segment, clipped to segment endpoints.
        t = ((x - ax) * vx + (y - ay) * vy) / length2
        t = np.clip(t, 0.0, 1.0)

        px = ax + t * vx
        py = ay + t * vy
        d2 = (x - px) ** 2 + (y - py) ** 2

        width = 18.0 + 0.04 * math.sqrt(length2)
        channel += -strength * np.exp(-d2 / (2.0 * width * width))

    return channel


def gradient_velocity(
    phi: np.ndarray,
    dx: float,
    dy: float,
    mobility: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Peculiar velocity follows downhill into sinks:

        v = -mu ∇Phi

    Since Phi is negative near attractors, this points inward.
    """
    dphi_dy, dphi_dx = np.gradient(phi, dy, dx)
    vx = -mobility * dphi_dx
    vy = -mobility * dphi_dy
    return vx, vy


def divergence(vx: np.ndarray, vy: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """2D divergence ∂vx/∂x + ∂vy/∂y."""
    _dvy_dy, dvy_dx_unused = np.gradient(vy, dy, dx)
    dvy_dy = _dvy_dy
    dvx_dy_unused, dvx_dx = np.gradient(vx, dy, dx)
    return dvx_dx + dvy_dy


def bilinear_sample(grid: np.ndarray, xs: np.ndarray, ys: np.ndarray, x_axis: np.ndarray, y_axis: np.ndarray) -> np.ndarray:
    """Sample 2D grid defined on y,x axes at arbitrary xs, ys."""
    nx = len(x_axis)
    ny = len(y_axis)
    dx = x_axis[1] - x_axis[0]
    dy = y_axis[1] - y_axis[0]

    fx = (xs - x_axis[0]) / dx
    fy = (ys - y_axis[0]) / dy

    x0 = np.floor(fx).astype(int)
    y0 = np.floor(fy).astype(int)
    x0 = np.clip(x0, 0, nx - 2)
    y0 = np.clip(y0, 0, ny - 2)
    x1 = x0 + 1
    y1 = y0 + 1

    wx = fx - x0
    wy = fy - y0

    q00 = grid[y0, x0]
    q10 = grid[y0, x1]
    q01 = grid[y1, x0]
    q11 = grid[y1, x1]

    return (
        q00 * (1 - wx) * (1 - wy)
        + q10 * wx * (1 - wy)
        + q01 * (1 - wx) * wy
        + q11 * wx * wy
    )


# ============================================================
# Galaxy generation and basin assignment
# ============================================================

def generate_galaxies(cfg: SimulationConfig, attractors: List[Attractor]) -> Dict[str, np.ndarray]:
    """
    Generate a toy supercluster sample.

    Mixes:
        - background field galaxies
        - overdensity clouds around attractors
        - filament galaxies between attractors
    """
    rng = np.random.default_rng(cfg.seed)
    half = cfg.box_mpc / 2.0

    n_bg = int(cfg.n_galaxies * 0.45)
    n_attr = int(cfg.n_galaxies * 0.35)
    n_fil = cfg.n_galaxies - n_bg - n_attr

    # Background.
    x_bg = rng.uniform(-half, half, n_bg)
    y_bg = rng.uniform(-half, half, n_bg)

    # Attractor clouds.
    xs_attr = []
    ys_attr = []
    weights = np.array([a.mass * a.flux_gain for a in attractors], dtype=float)
    weights /= weights.sum()
    counts = rng.multinomial(n_attr, weights)
    for count, a in zip(counts, attractors):
        xs_attr.append(rng.normal(a.x, a.core * 1.4, count))
        ys_attr.append(rng.normal(a.y, a.core * 1.4, count))
    x_attr = np.concatenate(xs_attr) if xs_attr else np.array([])
    y_attr = np.concatenate(ys_attr) if ys_attr else np.array([])

    # Filaments.
    pairs = [(0, 1), (0, 2), (1, 2)]
    pair_weights = np.array([1.0, 0.8, 0.7])
    pair_weights /= pair_weights.sum()
    pair_counts = rng.multinomial(n_fil, pair_weights)
    xs_fil = []
    ys_fil = []
    for count, (i, j) in zip(pair_counts, pairs):
        a = attractors[i]
        b = attractors[j]
        t = rng.uniform(0.0, 1.0, count)
        line_x = a.x + t * (b.x - a.x)
        line_y = a.y + t * (b.y - a.y)
        xs_fil.append(line_x + rng.normal(0.0, 11.0, count))
        ys_fil.append(line_y + rng.normal(0.0, 11.0, count))
    x_fil = np.concatenate(xs_fil) if xs_fil else np.array([])
    y_fil = np.concatenate(ys_fil) if ys_fil else np.array([])

    x = np.concatenate([x_bg, x_attr, x_fil])
    y = np.concatenate([y_bg, y_attr, y_fil])
    x = np.clip(x, -half, half)
    y = np.clip(y, -half, half)

    # Simple galaxy mass/luminosity proxy.
    stellar_mass = 10 ** rng.normal(10.3, 0.45, len(x))

    return {"x": x, "y": y, "stellar_mass": stellar_mass}


def assign_basin(xs: np.ndarray, ys: np.ndarray, attractors: List[Attractor]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Assign each galaxy to the attractor with strongest flux pull proxy:

        pull ~ mass*flux_gain / (r^2 + core^2)
    """
    pulls = []
    for a in attractors:
        r2 = (xs - a.x) ** 2 + (ys - a.y) ** 2 + a.core ** 2
        pulls.append(a.mass * a.flux_gain / r2)
    pulls = np.vstack(pulls)
    labels = np.argmax(pulls, axis=0)
    max_pull = np.max(pulls, axis=0)
    return labels, max_pull


# ============================================================
# Simulation
# ============================================================

def run_supercluster(cfg: SimulationConfig) -> Dict[str, object]:
    attractors = default_attractors()
    outdir = Path(cfg.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    half = cfg.box_mpc / 2.0
    x_axis = np.linspace(-half, half, cfg.grid_n)
    y_axis = np.linspace(-half, half, cfg.grid_n)
    dx = x_axis[1] - x_axis[0]
    dy = y_axis[1] - y_axis[0]
    X, Y = np.meshgrid(x_axis, y_axis)

    phi_mass = potential_field(X, Y, attractors, use_flux=False)
    phi_flux = potential_field(X, Y, attractors, use_flux=True)
    channels = filament_channel_field(X, Y, attractors, strength=cfg.filament_strength)
    phi_total = phi_flux + channels

    vx_flux, vy_flux = gradient_velocity(phi_total, dx=dx, dy=dy, mobility=cfg.mobility)
    vx_mass, vy_mass = gradient_velocity(phi_mass, dx=dx, dy=dy, mobility=cfg.mobility)

    div_flux = divergence(vx_flux, vy_flux, dx=dx, dy=dy)
    div_mass = divergence(vx_mass, vy_mass, dx=dx, dy=dy)

    convergence = np.maximum(0.0, -div_flux)
    convergence_mass = np.maximum(0.0, -div_mass)
    residual_convergence = convergence - convergence_mass

    cell_area = dx * dy
    s_flux = float(np.sum(convergence) * cell_area)
    s_mass = float(np.sum(convergence_mass) * cell_area)
    s_residual = float(np.sum(np.maximum(0.0, residual_convergence)) * cell_area)

    galaxies = generate_galaxies(cfg, attractors)
    basin_labels, basin_pull = assign_basin(galaxies["x"], galaxies["y"], attractors)

    rng = np.random.default_rng(cfg.seed + 10_000)
    gx = galaxies["x"]
    gy = galaxies["y"]
    gvx = bilinear_sample(vx_flux, gx, gy, x_axis, y_axis) + rng.normal(0.0, cfg.noise_sigma, len(gx))
    gvy = bilinear_sample(vy_flux, gx, gy, x_axis, y_axis) + rng.normal(0.0, cfg.noise_sigma, len(gx))
    gvx_mass = bilinear_sample(vx_mass, gx, gy, x_axis, y_axis)
    gvy_mass = bilinear_sample(vy_mass, gx, gy, x_axis, y_axis)

    residual_speed = np.sqrt((gvx - gvx_mass) ** 2 + (gvy - gvy_mass) ** 2)
    flow_speed = np.sqrt(gvx ** 2 + gvy ** 2)

    # Basin-level diagnostics.
    basin_stats = []
    for i, a in enumerate(attractors):
        mask = basin_labels == i
        if np.any(mask):
            basin_stats.append({
                "name": a.name,
                "n": int(mask.sum()),
                "mean_flow_speed": float(np.mean(flow_speed[mask])),
                "mean_residual_speed": float(np.mean(residual_speed[mask])),
                "stellar_mass_sum": float(np.sum(galaxies["stellar_mass"][mask])),
                "sink_strength": float(a.mass * a.flux_gain),
            })
        else:
            basin_stats.append({
                "name": a.name,
                "n": 0,
                "mean_flow_speed": 0.0,
                "mean_residual_speed": 0.0,
                "stellar_mass_sum": 0.0,
                "sink_strength": float(a.mass * a.flux_gain),
            })

    results = {
        "cfg": cfg,
        "attractors": attractors,
        "x_axis": x_axis,
        "y_axis": y_axis,
        "X": X,
        "Y": Y,
        "phi_mass": phi_mass,
        "phi_flux": phi_flux,
        "channels": channels,
        "phi_total": phi_total,
        "vx_flux": vx_flux,
        "vy_flux": vy_flux,
        "vx_mass": vx_mass,
        "vy_mass": vy_mass,
        "div_flux": div_flux,
        "div_mass": div_mass,
        "convergence": convergence,
        "residual_convergence": residual_convergence,
        "galaxies": galaxies,
        "basin_labels": basin_labels,
        "basin_pull": basin_pull,
        "gvx": gvx,
        "gvy": gvy,
        "gvx_mass": gvx_mass,
        "gvy_mass": gvy_mass,
        "residual_speed": residual_speed,
        "flow_speed": flow_speed,
        "s_flux": s_flux,
        "s_mass": s_mass,
        "s_residual": s_residual,
        "basin_stats": basin_stats,
    }

    write_outputs(results)
    if cfg.make_plots:
        make_plots(results)

    return results


# ============================================================
# Output writers
# ============================================================

def write_outputs(results: Dict[str, object]) -> None:
    cfg: SimulationConfig = results["cfg"]  # type: ignore[assignment]
    outdir = Path(cfg.output_dir)
    galaxies: Dict[str, np.ndarray] = results["galaxies"]  # type: ignore[assignment]
    attractors: List[Attractor] = results["attractors"]  # type: ignore[assignment]
    basin_labels: np.ndarray = results["basin_labels"]  # type: ignore[assignment]

    # Summary.
    summary_path = outdir / "supercluster_summary.txt"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write("Supercluster Flux Sink Scan\n")
        f.write("===========================\n")
        f.write(f"N galaxies                 : {cfg.n_galaxies}\n")
        f.write(f"Grid                       : {cfg.grid_n} x {cfg.grid_n}\n")
        f.write(f"Box size                   : {cfg.box_mpc:.1f} Mpc\n")
        f.write(f"S_flux convergence index   : {results['s_flux']:.6e}\n")
        f.write(f"S_mass-only index          : {results['s_mass']:.6e}\n")
        f.write(f"S_residual flux index      : {results['s_residual']:.6e}\n")
        f.write("\nBasin diagnostics\n")
        f.write("-----------------\n")
        for row in results["basin_stats"]:  # type: ignore[index]
            f.write(
                f"{row['name']:16s} | n={row['n']:4d} | "
                f"sink={row['sink_strength']:.3f} | "
                f"mean_flow={row['mean_flow_speed']:.3f} | "
                f"mean_residual={row['mean_residual_speed']:.3f} | "
                f"stellar_mass={row['stellar_mass_sum']:.3e}\n"
            )
        f.write("\nInterpretation\n")
        f.write("--------------\n")
        f.write("High S_flux marks a supercluster-scale coherence sink.\n")
        f.write("Positive S_residual means the flux model predicts convergence beyond a mass-only field.\n")
        f.write("Filament-aligned residuals support the active-channel reading of cosmic-web structure.\n")

    # Galaxy catalogue.
    gal_path = outdir / "supercluster_galaxies.csv"
    with gal_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "x_mpc", "y_mpc", "stellar_mass_proxy", "basin", "basin_pull",
            "vx_flux_observed", "vy_flux_observed", "vx_mass_only", "vy_mass_only",
            "flow_speed", "residual_speed",
        ])
        for i in range(len(galaxies["x"])):
            label = int(basin_labels[i])
            writer.writerow([
                galaxies["x"][i],
                galaxies["y"][i],
                galaxies["stellar_mass"][i],
                attractors[label].name,
                results["basin_pull"][i],  # type: ignore[index]
                results["gvx"][i],         # type: ignore[index]
                results["gvy"][i],         # type: ignore[index]
                results["gvx_mass"][i],    # type: ignore[index]
                results["gvy_mass"][i],    # type: ignore[index]
                results["flow_speed"][i],  # type: ignore[index]
                results["residual_speed"][i],  # type: ignore[index]
            ])

    # Grid catalogue, coarse but complete.
    grid_path = outdir / "flux_sink_grid.csv"
    x_axis: np.ndarray = results["x_axis"]  # type: ignore[assignment]
    y_axis: np.ndarray = results["y_axis"]  # type: ignore[assignment]
    with grid_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "x_mpc", "y_mpc", "phi_total", "phi_mass", "filament_channel",
            "vx_flux", "vy_flux", "vx_mass", "vy_mass", "div_flux", "convergence", "residual_convergence",
        ])
        for iy, y in enumerate(y_axis):
            for ix, x in enumerate(x_axis):
                writer.writerow([
                    x, y,
                    results["phi_total"][iy, ix],        # type: ignore[index]
                    results["phi_mass"][iy, ix],         # type: ignore[index]
                    results["channels"][iy, ix],         # type: ignore[index]
                    results["vx_flux"][iy, ix],          # type: ignore[index]
                    results["vy_flux"][iy, ix],          # type: ignore[index]
                    results["vx_mass"][iy, ix],          # type: ignore[index]
                    results["vy_mass"][iy, ix],          # type: ignore[index]
                    results["div_flux"][iy, ix],         # type: ignore[index]
                    results["convergence"][iy, ix],      # type: ignore[index]
                    results["residual_convergence"][iy, ix],  # type: ignore[index]
                ])


# ============================================================
# Plots
# ============================================================

def _save_imshow(path: Path, field: np.ndarray, extent: List[float], title: str, cbar: str) -> None:
    plt.figure(figsize=(8, 7))
    plt.imshow(field, origin="lower", extent=extent, aspect="equal")
    plt.colorbar(label=cbar)
    plt.xlabel("x [Mpc]")
    plt.ylabel("y [Mpc]")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def make_plots(results: Dict[str, object]) -> None:
    cfg: SimulationConfig = results["cfg"]  # type: ignore[assignment]
    outdir = Path(cfg.output_dir)
    half = cfg.box_mpc / 2.0
    extent = [-half, half, -half, half]

    galaxies: Dict[str, np.ndarray] = results["galaxies"]  # type: ignore[assignment]
    attractors: List[Attractor] = results["attractors"]  # type: ignore[assignment]
    labels: np.ndarray = results["basin_labels"]  # type: ignore[assignment]

    _save_imshow(
        outdir / "velocity_convergence_map.png",
        results["convergence"],  # type: ignore[arg-type]
        extent,
        "Supercluster velocity convergence: max(0, -div v)",
        "convergence index",
    )

    _save_imshow(
        outdir / "flux_potential_map.png",
        results["phi_total"],  # type: ignore[arg-type]
        extent,
        "Flux potential with supercluster sinks and filament channels",
        "Phi_flux",
    )

    _save_imshow(
        outdir / "filament_channel_map.png",
        results["channels"],  # type: ignore[arg-type]
        extent,
        "Filament-channel flux field",
        "channel potential",
    )

    _save_imshow(
        outdir / "residual_flow_map.png",
        results["residual_convergence"],  # type: ignore[arg-type]
        extent,
        "Residual convergence: flux field minus mass-only field",
        "residual convergence",
    )

    # Basin assignment scatter.
    plt.figure(figsize=(8, 7))
    plt.scatter(galaxies["x"], galaxies["y"], c=labels, s=6, alpha=0.65)
    for i, a in enumerate(attractors):
        plt.scatter([a.x], [a.y], marker="x", s=130)
        plt.text(a.x + 3, a.y + 3, a.name.replace("_", " "), fontsize=9)
    plt.xlabel("x [Mpc]")
    plt.ylabel("y [Mpc]")
    plt.title("Basin assignment map: galaxy envelopes drain to nearest flux sink")
    plt.xlim(-half, half)
    plt.ylim(-half, half)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(outdir / "basin_assignment_map.png", dpi=180)
    plt.close()

    # Flow arrows over potential.
    X: np.ndarray = results["X"]  # type: ignore[assignment]
    Y: np.ndarray = results["Y"]  # type: ignore[assignment]
    vx: np.ndarray = results["vx_flux"]  # type: ignore[assignment]
    vy: np.ndarray = results["vy_flux"]  # type: ignore[assignment]
    phi: np.ndarray = results["phi_total"]  # type: ignore[assignment]

    step = max(1, cfg.grid_n // 28)
    plt.figure(figsize=(8, 7))
    plt.imshow(phi, origin="lower", extent=extent, aspect="equal")
    plt.quiver(X[::step, ::step], Y[::step, ::step], vx[::step, ::step], vy[::step, ::step], scale=0.22)
    for a in attractors:
        plt.scatter([a.x], [a.y], marker="x", s=120)
    plt.xlabel("x [Mpc]")
    plt.ylabel("y [Mpc]")
    plt.title("Peculiar flow into supercluster flux basins")
    plt.tight_layout()
    plt.savefig(outdir / "peculiar_flow_vectors.png", dpi=180)
    plt.close()


# ============================================================
# CLI
# ============================================================

def parse_args() -> SimulationConfig:
    parser = argparse.ArgumentParser(description="Supercluster-scale flux sink toy simulation")
    parser.add_argument("--n-galaxies", type=int, default=1200)
    parser.add_argument("--box", type=float, default=260.0, help="box width in Mpc")
    parser.add_argument("--grid", type=int, default=150, help="grid resolution per axis")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--mobility", type=float, default=0.085)
    parser.add_argument("--noise", type=float, default=18.0)
    parser.add_argument("--filament-strength", type=float, default=0.35)
    parser.add_argument("--output-dir", type=str, default="supercluster_outputs")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    return SimulationConfig(
        n_galaxies=args.n_galaxies,
        box_mpc=args.box,
        grid_n=args.grid,
        seed=args.seed,
        mobility=args.mobility,
        noise_sigma=args.noise,
        filament_strength=args.filament_strength,
        output_dir=args.output_dir,
        make_plots=not args.no_plots,
    )


def main() -> None:
    cfg = parse_args()
    results = run_supercluster(cfg)

    print("Supercluster flux scan complete")
    print("===============================")
    print(f"N galaxies               : {cfg.n_galaxies}")
    print(f"Grid                     : {cfg.grid_n} x {cfg.grid_n}")
    print(f"S_flux convergence index : {results['s_flux']:.6e}")
    print(f"S_mass-only index        : {results['s_mass']:.6e}")
    print(f"S_residual flux index    : {results['s_residual']:.6e}")
    print(f"Outputs                  : {Path(cfg.output_dir).resolve()}")

    print("\nBasin diagnostics:")
    for row in results["basin_stats"]:  # type: ignore[index]
        print(
            f"  {row['name']:16s} n={row['n']:4d} "
            f"sink={row['sink_strength']:.3f} "
            f"mean_flow={row['mean_flow_speed']:.3f} "
            f"mean_residual={row['mean_residual_speed']:.3f}"
        )


if __name__ == "__main__":
    main()
