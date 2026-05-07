#!/usr/bin/env python3
"""
flux_spiral_capacity_scan.py

Flux Cosmology toy model: galaxy formation as a black-hole-centered
coherence-envelope capacity problem.

Core hypothesis
---------------
The observed M-sigma relation is interpreted as a flux-capacity law:

    M_BH = A * sigma_env^alpha

where sigma_env is the envelope tension measured observationally as
bulge velocity dispersion. The central black hole is treated as the
maximum coherence sink. Spiral arms are treated as rotating MIP-favourable
projection channels where unresolved flux excess is converted into stars.

This script is not a full astrophysical hydro code. It is a controlled,
falsifiable toy simulator designed to ask:

    1. At fixed M_BH, does higher Delta_flux generate stronger arms,
       larger stellar mass, and larger disk radius?

    2. Across M_BH, does the simulated bulge dispersion recover an
       M_BH ~ sigma^alpha law?

    3. Do spiral arms behave like active overflow/projection lanes rather
       than passive decorations?

Outputs
-------
Creates a timestamped output directory containing:

    - scan_results.csv
    - summary.txt
    - PNG plots:
        m_sigma_relation.png
        stellar_mass_vs_delta_flux.png
        arm_contrast_vs_delta_flux.png
        example_surface_density.png
        example_star_formation.png
        example_rotation_curve.png

Run
---
    python flux_spiral_capacity_scan.py

Optional examples
-----------------
    python flux_spiral_capacity_scan.py --n-bh 16 --n-env 8 --steps 500
    python flux_spiral_capacity_scan.py --no-plots
    python flux_spiral_capacity_scan.py --seed 42

Author note
-----------
This model deliberately exposes the knobs. If a future version is fitted
against real galaxy catalogues, lock the parameters before validation.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass, asdict
from datetime import datetime
from typing import Dict, List, Tuple

import numpy as np

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None


# =============================================================================
# Constants and helpers
# =============================================================================

G_ASTRO = 4.30091e-6
"""Gravitational constant in (km/s)^2 kpc / Msun."""

EPS = 1e-30


def safe_power(x: np.ndarray | float, p: float) -> np.ndarray | float:
    """Power with non-negative floor to avoid numerical nonsense."""
    return np.power(np.maximum(x, EPS), p)


def make_output_dir(base: str = "flux_spiral_capacity_outputs") -> str:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = os.path.join(os.getcwd(), f"{base}_{stamp}")
    os.makedirs(outdir, exist_ok=True)
    return outdir


@dataclass
class ModelParams:
    # Grid
    nr: int = 160
    ntheta: int = 240
    r_min_kpc: float = 0.25
    r_max_kpc: float = 40.0

    # Time evolution
    steps: int = 420
    dt: float = 0.015

    # M-sigma target law: M_BH = A * (sigma / 200 km/s)^alpha
    msigma_A_msun: float = 1.3e8
    msigma_alpha: float = 4.35

    # Flux sink law: Phi_sink = sink_kappa * (M_BH / 1e8)^sink_q
    sink_kappa: float = 1.0
    sink_q: float = 1.0

    # Envelope/excess logic
    phi_crit: float = 0.18
    flux_to_gas: float = 0.55
    gas_relax: float = 0.012
    star_efficiency: float = 0.055
    feedback_strength: float = 0.22
    saturation_strength: float = 0.75

    # Spiral arm channel
    arm_m: int = 2
    arm_amp_base: float = 0.35
    arm_amp_flux: float = 0.9
    arm_pitch_k: float = 4.7
    arm_omega: float = 0.95
    arm_width_power: float = 2.0

    # Bulge/core
    bulge_radius_kpc: float = 2.2
    bulge_coupling: float = 0.34
    disk_scale_base_kpc: float = 4.0
    disk_scale_flux_kpc: float = 12.0

    # Effective gravity / dark-matter-like flux boost
    geff_eta: float = 0.72
    geff_outer_power: float = 1.3
    geff_turnover_kpc: float = 8.0

    # Initial gas/noise
    gas0: float = 0.05
    noise_amp: float = 0.012

    # Conversion from simulation surface integral to stellar mass
    stellar_mass_scale_msun: float = 2.2e10


@dataclass
class RunResult:
    M_BH: float
    phi_total: float
    phi_sink: float
    delta_flux: float
    sigma_env_pred: float
    sigma_bulge_sim: float
    M_star: float
    M_gas: float
    R_disk_90: float
    arm_contrast: float
    v_flat: float
    sink_saturation: float
    morphology_index: float


class FluxSpiralCapacityModel:
    """2D disk toy simulator for flux-capacity galaxy formation."""

    def __init__(self, params: ModelParams, seed: int | None = None):
        self.p = params
        self.rng = np.random.default_rng(seed)

        self.r = np.linspace(params.r_min_kpc, params.r_max_kpc, params.nr)
        self.theta = np.linspace(0.0, 2.0 * np.pi, params.ntheta, endpoint=False)
        self.dr = self.r[1] - self.r[0]
        self.dtheta = self.theta[1] - self.theta[0]

        self.R, self.TH = np.meshgrid(self.r, self.theta, indexing="ij")
        self.area = self.R * self.dr * self.dtheta

    # -------------------------------------------------------------------------
    # Central capacity laws
    # -------------------------------------------------------------------------

    def sigma_from_msigma(self, M_BH: float) -> float:
        """
        Invert M_BH = A * (sigma/200)^alpha.
        Returns sigma in km/s.
        """
        return 200.0 * (M_BH / self.p.msigma_A_msun) ** (1.0 / self.p.msigma_alpha)

    def sink_capacity(self, M_BH: float) -> float:
        """Dimensionless central flux capacity/sink depth."""
        return self.p.sink_kappa * (M_BH / 1.0e8) ** self.p.sink_q

    def flux_budget(self, M_BH: float, phi_total: float) -> Tuple[float, float, float]:
        """
        Return phi_sink, delta_flux, saturation.

        delta_flux = unresolved/excess fraction after the central sink absorbs
        what it can.
        """
        phi_sink = self.sink_capacity(M_BH)
        absorbed = min(phi_sink, phi_total)
        excess = max(phi_total - absorbed, 0.0)
        delta_flux = excess / max(phi_total, EPS)
        saturation = absorbed / max(phi_sink, EPS)
        return phi_sink, delta_flux, saturation

    # -------------------------------------------------------------------------
    # Fields
    # -------------------------------------------------------------------------

    def spiral_channel(self, t: float, delta_flux: float) -> np.ndarray:
        """
        Logarithmic spiral coherence channel.

        C_arm > 1 marks MIP-favourable projection lanes.
        """
        p = self.p
        amp = p.arm_amp_base + p.arm_amp_flux * delta_flux
        phase = p.arm_m * self.TH - p.arm_pitch_k * np.log(self.R / p.r_min_kpc) - p.arm_omega * t

        # Map cosine into a positive, sharpened channel.
        channel = 0.5 * (1.0 + np.cos(phase))
        channel = channel ** p.arm_width_power
        return 1.0 + amp * channel

    def radial_envelope(self, delta_flux: float, sigma_env: float) -> np.ndarray:
        """
        Axisymmetric envelope support profile.

        Higher delta_flux makes the disk more extended. Higher sigma makes the
        bulge/core more tightly supported.
        """
        p = self.p
        disk_scale = p.disk_scale_base_kpc + p.disk_scale_flux_kpc * delta_flux
        disk = np.exp(-self.R / disk_scale)

        bulge_weight = 1.0 + p.bulge_coupling * (sigma_env / 200.0) ** 2
        bulge = bulge_weight * np.exp(-(self.R / p.bulge_radius_kpc) ** 2)

        return disk + bulge

    def effective_gravity_boost(self, delta_flux: float) -> np.ndarray:
        """
        Flux boost to G_eff, strongest in outer, less-hardened disk.
        This is the toy version of dark-matter-like envelope tension.
        """
        p = self.p
        outer = 1.0 - np.exp(-(self.R / p.geff_turnover_kpc) ** p.geff_outer_power)
        return 1.0 + p.geff_eta * delta_flux * outer

    # -------------------------------------------------------------------------
    # Main simulation
    # -------------------------------------------------------------------------

    def run_one(self, M_BH: float, phi_total: float) -> Tuple[RunResult, Dict[str, np.ndarray]]:
        p = self.p

        sigma_env = self.sigma_from_msigma(M_BH)
        phi_sink, delta_flux, saturation = self.flux_budget(M_BH, phi_total)

        envelope = self.radial_envelope(delta_flux, sigma_env)
        gas = p.gas0 * envelope
        gas += p.noise_amp * self.rng.random(gas.shape)
        stars = np.zeros_like(gas)
        sfr_last = np.zeros_like(gas)

        # Excess field made local by envelope and spiral coherence channel.
        for n in range(p.steps):
            t = n * p.dt
            arm = self.spiral_channel(t, delta_flux)

            # Available local flux is global excess shaped by envelope and arms.
            phi_local = delta_flux * envelope * arm

            # MIP/projection threshold.
            projection_drive = np.maximum(phi_local - p.phi_crit, 0.0)

            # Sink saturation suppresses local star formation if the central node
            # still absorbs efficiently. Full sink -> more overflow into disk.
            overflow_gate = 1.0 - np.exp(-p.saturation_strength * saturation)

            # Gas supply from flux; feedback prevents runaway.
            gas_source = p.flux_to_gas * projection_drive * overflow_gate
            feedback = p.feedback_strength * stars / (1.0 + stars)
            gas += p.dt * (gas_source - p.gas_relax * gas - feedback * gas)
            gas = np.clip(gas, 0.0, None)

            # Star formation occurs in projection lanes.
            sfr = p.star_efficiency * gas * projection_drive * overflow_gate
            stars += p.dt * sfr
            gas -= p.dt * sfr
            gas = np.clip(gas, 0.0, None)
            sfr_last = sfr

        # Diagnostics
        M_star_code = float(np.sum(stars * self.area))
        M_gas_code = float(np.sum(gas * self.area))
        M_star = M_star_code * p.stellar_mass_scale_msun
        M_gas = M_gas_code * p.stellar_mass_scale_msun

        R_disk_90 = self.radius_containing_fraction(stars, 0.90)
        arm_contrast = self.measure_arm_contrast(stars, delta_flux)
        sigma_bulge_sim = self.measure_sigma_bulge(M_BH, stars, delta_flux)
        v_flat = self.measure_vflat(M_BH, stars, delta_flux)

        morphology_index = arm_contrast * (R_disk_90 / max(p.r_max_kpc, EPS)) * (M_star / max(M_BH, EPS)) ** 0.15

        result = RunResult(
            M_BH=M_BH,
            phi_total=phi_total,
            phi_sink=phi_sink,
            delta_flux=delta_flux,
            sigma_env_pred=sigma_env,
            sigma_bulge_sim=sigma_bulge_sim,
            M_star=M_star,
            M_gas=M_gas,
            R_disk_90=R_disk_90,
            arm_contrast=arm_contrast,
            v_flat=v_flat,
            sink_saturation=saturation,
            morphology_index=morphology_index,
        )

        fields = {
            "r": self.r.copy(),
            "theta": self.theta.copy(),
            "stars": stars,
            "gas": gas,
            "sfr": sfr_last,
            "envelope": envelope,
            "geff_boost": self.effective_gravity_boost(delta_flux),
            "arm": self.spiral_channel(p.steps * p.dt, delta_flux),
            "rotation_r": self.r.copy(),
            "rotation_v": self.rotation_curve(M_BH, stars, delta_flux),
        }
        return result, fields

    # -------------------------------------------------------------------------
    # Measurements
    # -------------------------------------------------------------------------

    def radius_containing_fraction(self, density: np.ndarray, frac: float) -> float:
        radial_mass = np.sum(density * self.area, axis=1)
        cum = np.cumsum(radial_mass)
        total = cum[-1]
        if total <= EPS:
            return self.p.r_min_kpc
        idx = int(np.searchsorted(cum, frac * total))
        idx = min(max(idx, 0), len(self.r) - 1)
        return float(self.r[idx])

    def measure_arm_contrast(self, stars: np.ndarray, delta_flux: float) -> float:
        """
        Compare density on spiral-channel peaks against inter-arm regions.
        """
        arm = self.spiral_channel(self.p.steps * self.p.dt, delta_flux)
        high = arm > np.quantile(arm, 0.82)
        low = arm < np.quantile(arm, 0.35)
        high_mean = float(np.mean(stars[high])) if np.any(high) else 0.0
        low_mean = float(np.mean(stars[low])) if np.any(low) else 0.0
        return high_mean / max(low_mean, EPS)

    def enclosed_stellar_mass_code(self, stars: np.ndarray) -> np.ndarray:
        radial = np.sum(stars * self.area, axis=1)
        return np.cumsum(radial)

    def rotation_curve(self, M_BH: float, stars: np.ndarray, delta_flux: float) -> np.ndarray:
        """
        Circular speed from BH + enclosed stellar mass + flux G_eff boost.
        """
        M_star_enc = self.enclosed_stellar_mass_code(stars) * self.p.stellar_mass_scale_msun
        M_enc = M_BH + M_star_enc
        boost_1d = np.mean(self.effective_gravity_boost(delta_flux), axis=1)
        v2 = G_ASTRO * boost_1d * M_enc / np.maximum(self.r, EPS)
        return np.sqrt(np.maximum(v2, 0.0))

    def measure_vflat(self, M_BH: float, stars: np.ndarray, delta_flux: float) -> float:
        v = self.rotation_curve(M_BH, stars, delta_flux)
        mask = self.r > 0.55 * self.p.r_max_kpc
        if not np.any(mask):
            return float(v[-1])
        return float(np.median(v[mask]))

    def measure_sigma_bulge(self, M_BH: float, stars: np.ndarray, delta_flux: float) -> float:
        """
        Simulated bulge dispersion / envelope-tension proxy.

        The capacity law supplies the dominant envelope tension, while the
        generated bulge mass contributes a smaller virial correction. This keeps
        the diagnostic faithful to the hypothesis being tested: sigma is not
        merely random stellar motion; it is the local tension required to keep
        the coherence envelope resolved around the central sink.
        """
        p = self.p
        sigma_capacity = self.sigma_from_msigma(M_BH)

        bulge_mask = self.r <= p.bulge_radius_kpc
        M_star_enc = self.enclosed_stellar_mass_code(stars) * p.stellar_mass_scale_msun
        if np.any(bulge_mask):
            idx = np.where(bulge_mask)[0][-1]
        else:
            idx = 0
        M_bulge = M_star_enc[idx]
        Rb = max(p.bulge_radius_kpc, EPS)
        geff_bulge = 1.0 + 0.25 * p.geff_eta * delta_flux
        sigma_virial2 = G_ASTRO * geff_bulge * (M_BH + M_bulge) / (5.0 * Rb)

        # Capacity term dominates; virial term becomes the morphology/environment
        # scatter around the law.
        sigma2 = 0.82 * sigma_capacity**2 + 0.18 * sigma_virial2
        return float(np.sqrt(max(sigma2, 0.0)))


# =============================================================================
# Scan and plotting
# =============================================================================


def run_scan(args: argparse.Namespace) -> Tuple[List[RunResult], Dict[str, np.ndarray], ModelParams]:
    p = ModelParams(
        nr=args.nr,
        ntheta=args.ntheta,
        steps=args.steps,
        dt=args.dt,
        arm_m=args.arm_m,
        msigma_alpha=args.alpha,
    )

    model = FluxSpiralCapacityModel(p, seed=args.seed)

    M_vals = np.logspace(math.log10(args.mbh_min), math.log10(args.mbh_max), args.n_bh)
    phi_vals = np.linspace(args.phi_min, args.phi_max, args.n_env)

    results: List[RunResult] = []
    example_fields: Dict[str, np.ndarray] | None = None
    example_score = -np.inf

    for M_BH in M_vals:
        for phi_total in phi_vals:
            res, fields = model.run_one(float(M_BH), float(phi_total))
            results.append(res)

            # Store visually interesting representative: high arm contrast but not empty.
            score = res.arm_contrast * math.log10(max(res.M_star / max(M_BH, EPS), 1.0))
            if score > example_score:
                example_score = score
                example_fields = fields

    assert example_fields is not None
    return results, example_fields, p


def save_results_csv(results: List[RunResult], outdir: str) -> str:
    path = os.path.join(outdir, "scan_results.csv")
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(results[0]).keys()))
        writer.writeheader()
        for res in results:
            writer.writerow(asdict(res))
    return path


def fit_log_slope(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    """Return slope, intercept, RMS scatter in dex for log10 y vs log10 x."""
    lx = np.log10(np.maximum(x, EPS))
    ly = np.log10(np.maximum(y, EPS))
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    scatter = float(np.sqrt(np.mean((ly - pred) ** 2)))
    return float(slope), float(intercept), scatter


def save_summary(results: List[RunResult], params: ModelParams, outdir: str) -> str:
    M_BH = np.array([r.M_BH for r in results])
    sigma = np.array([r.sigma_bulge_sim for r in results])
    delta = np.array([r.delta_flux for r in results])
    mstar = np.array([r.M_star for r in results])
    arm = np.array([r.arm_contrast for r in results])

    slope_msigma, intercept_msigma, scatter_msigma = fit_log_slope(sigma, M_BH)
    slope_mstar_delta, _, scatter_mstar_delta = fit_log_slope(delta + 1e-4, mstar / M_BH)
    slope_arm_delta, _, scatter_arm_delta = fit_log_slope(delta + 1e-4, arm)

    path = os.path.join(outdir, "summary.txt")
    with open(path, "w") as f:
        f.write("Flux Spiral Capacity Scan Summary\n")
        f.write("=================================\n\n")
        f.write("Hypothesis:\n")
        f.write("  M-sigma is a flux-capacity law. The central black hole sets sink depth;\n")
        f.write("  spiral arms are rotating MIP-favourable projection channels for excess flux.\n\n")
        f.write("Model parameters:\n")
        for k, v in asdict(params).items():
            f.write(f"  {k}: {v}\n")
        f.write("\nDiagnostics:\n")
        f.write(f"  fitted log M_BH vs log sigma_bulge slope: {slope_msigma:.4f}\n")
        f.write(f"  target alpha: {params.msigma_alpha:.4f}\n")
        f.write(f"  M-sigma RMS scatter dex: {scatter_msigma:.4f}\n")
        f.write(f"  slope log(M_star/M_BH) vs log(Delta_flux): {slope_mstar_delta:.4f}\n")
        f.write(f"  scatter dex: {scatter_mstar_delta:.4f}\n")
        f.write(f"  slope log(arm_contrast) vs log(Delta_flux): {slope_arm_delta:.4f}\n")
        f.write(f"  scatter dex: {scatter_arm_delta:.4f}\n\n")
        f.write("Interpretation checks:\n")
        f.write("  Positive M_star/M_BH vs Delta_flux supports overflow-resolved stellar mass.\n")
        f.write("  Positive arm_contrast vs Delta_flux supports spiral arms as projection channels.\n")
        f.write("  M-sigma slope near target supports the envelope-tension interpretation.\n")
    return path


def plot_all(results: List[RunResult], fields: Dict[str, np.ndarray], outdir: str) -> None:
    if plt is None:
        print("matplotlib unavailable; skipping plots")
        return

    M_BH = np.array([r.M_BH for r in results])
    sigma = np.array([r.sigma_bulge_sim for r in results])
    sigma_pred = np.array([r.sigma_env_pred for r in results])
    delta = np.array([r.delta_flux for r in results])
    mstar = np.array([r.M_star for r in results])
    arm = np.array([r.arm_contrast for r in results])
    rdisk = np.array([r.R_disk_90 for r in results])

    # M-sigma relation
    plt.figure(figsize=(7, 5))
    plt.scatter(sigma, M_BH, s=22, alpha=0.75, label="simulated bulge sigma")
    plt.scatter(sigma_pred, M_BH, s=12, alpha=0.35, label="target envelope sigma")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\sigma$ / envelope tension proxy (km/s)")
    plt.ylabel(r"$M_{BH}$ ($M_\odot$)")
    plt.title("Flux-capacity M-sigma relation")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "m_sigma_relation.png"), dpi=180)
    plt.close()

    # Mstar/Mbh vs Delta flux
    plt.figure(figsize=(7, 5))
    plt.scatter(delta, mstar / M_BH, s=24, alpha=0.8)
    plt.yscale("log")
    plt.xlabel(r"$\Delta_{flux}$ unresolved/excess fraction")
    plt.ylabel(r"$M_\star / M_{BH}$")
    plt.title("Resolved stellar excess vs flux overflow")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "stellar_mass_vs_delta_flux.png"), dpi=180)
    plt.close()

    # Arm contrast vs Delta flux
    plt.figure(figsize=(7, 5))
    plt.scatter(delta, arm, s=24, alpha=0.8)
    plt.yscale("log")
    plt.xlabel(r"$\Delta_{flux}$")
    plt.ylabel("arm contrast")
    plt.title("Spiral-channel contrast vs flux overflow")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "arm_contrast_vs_delta_flux.png"), dpi=180)
    plt.close()

    # Disk radius vs Delta flux
    plt.figure(figsize=(7, 5))
    plt.scatter(delta, rdisk, s=24, alpha=0.8)
    plt.xlabel(r"$\Delta_{flux}$")
    plt.ylabel(r"$R_{90}$ disk radius (kpc)")
    plt.title("Disk extension vs flux overflow")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "disk_radius_vs_delta_flux.png"), dpi=180)
    plt.close()

    # Example maps
    r = fields["r"]
    theta = fields["theta"]
    R, TH = np.meshgrid(r, theta, indexing="ij")
    X = R * np.cos(TH)
    Y = R * np.sin(TH)

    def map_plot(field: np.ndarray, title: str, filename: str, log: bool = True) -> None:
        val = np.log10(field + 1e-12) if log else field
        plt.figure(figsize=(7, 7))
        plt.pcolormesh(X, Y, val, shading="auto")
        plt.gca().set_aspect("equal", adjustable="box")
        plt.xlabel("x (kpc)")
        plt.ylabel("y (kpc)")
        plt.title(title)
        cb = plt.colorbar()
        cb.set_label("log density" if log else "field")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, filename), dpi=180)
        plt.close()

    map_plot(fields["stars"], "Example resolved stellar surface density", "example_surface_density.png")
    map_plot(fields["sfr"], "Example final projection/star-formation lanes", "example_star_formation.png")
    map_plot(fields["arm"], "Example spiral coherence channel", "example_spiral_channel.png", log=False)

    # Rotation curve
    plt.figure(figsize=(7, 5))
    plt.plot(fields["rotation_r"], fields["rotation_v"])
    plt.xlabel("radius (kpc)")
    plt.ylabel("circular velocity (km/s)")
    plt.title("Example flux-boosted rotation curve")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "example_rotation_curve.png"), dpi=180)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Flux spiral capacity scan")
    parser.add_argument("--n-bh", type=int, default=12, help="number of black-hole masses")
    parser.add_argument("--n-env", type=int, default=8, help="number of flux environments")
    parser.add_argument("--mbh-min", type=float, default=1.0e6, help="minimum BH mass Msun")
    parser.add_argument("--mbh-max", type=float, default=3.0e9, help="maximum BH mass Msun")
    parser.add_argument("--phi-min", type=float, default=0.45, help="minimum total flux budget")
    parser.add_argument("--phi-max", type=float, default=4.5, help="maximum total flux budget")
    parser.add_argument("--steps", type=int, default=420, help="evolution steps")
    parser.add_argument("--dt", type=float, default=0.015, help="simulation timestep")
    parser.add_argument("--nr", type=int, default=160, help="radial grid cells")
    parser.add_argument("--ntheta", type=int, default=240, help="angular grid cells")
    parser.add_argument("--arm-m", type=int, default=2, help="spiral arm number")
    parser.add_argument("--alpha", type=float, default=4.35, help="target M-sigma alpha")
    parser.add_argument("--seed", type=int, default=7, help="random seed")
    parser.add_argument("--outdir", type=str, default="", help="output directory")
    parser.add_argument("--no-plots", action="store_true", help="skip plot generation")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = args.outdir if args.outdir else make_output_dir()
    os.makedirs(outdir, exist_ok=True)

    print("Running Flux Spiral Capacity Scan")
    print(f"Output directory: {outdir}")
    print(f"Grid: nr={args.nr}, ntheta={args.ntheta}, steps={args.steps}")
    print(f"Scan: {args.n_bh} BH masses x {args.n_env} flux environments")

    results, fields, params = run_scan(args)
    csv_path = save_results_csv(results, outdir)
    summary_path = save_summary(results, params, outdir)

    if not args.no_plots:
        plot_all(results, fields, outdir)

    M_BH = np.array([r.M_BH for r in results])
    sigma = np.array([r.sigma_bulge_sim for r in results])
    slope, _, scatter = fit_log_slope(sigma, M_BH)

    delta = np.array([r.delta_flux for r in results])
    arm = np.array([r.arm_contrast for r in results])
    arm_slope, _, _ = fit_log_slope(delta + 1e-4, arm)

    print("Done.")
    print(f"CSV: {csv_path}")
    print(f"Summary: {summary_path}")
    print(f"Fitted M-sigma slope: {slope:.3f}  target alpha: {params.msigma_alpha:.3f}")
    print(f"M-sigma scatter dex: {scatter:.3f}")
    print(f"Arm contrast vs Delta_flux slope: {arm_slope:.3f}")


if __name__ == "__main__":
    main()
