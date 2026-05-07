#!/usr/bin/env python3
"""
flux_sound_horizon_scan.py

Flux Cosmology CMB boundary-condition scanner.

Purpose
-------
This script turns the CMB sound horizon into a hard validator for the
flux-scalar backend. It scans flux-potential parameters and reports whether
the first large-scale render boundary survives:

    r_s(a*)   = integral_0^a* c_s / (a^2 H(a)) da
    chi_*     = integral_a*^1 1 / (a^2 H(a)) da
    theta_*   = r_s / chi_*
    ell_A     = pi / theta_*

Flux interpretation
-------------------
    r_s      : maximum coherent communication length before recombination
    chi_*    : projection distance from last scattering to today
    theta_*  : observed angular size of the first coherence cell
    ell_A    : acoustic angular scale

The code is intentionally standalone and uses the same structural backend as
cmbr.py: a sigma flux/coherence scalar with potential

    V(sigma) = V0 * (1 - exp(-sigma / mu))^2

Dimensionless units are used with 8*pi*G/3 = 1. The absolute physical value of
r_s is therefore not the primary validator here; the dimensionless ratio
r_s / chi_* is the boundary condition.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import csv
import math
from typing import Dict, Iterable, List, Tuple

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# =============================================================================
# Baseline cosmology / target boundary
# =============================================================================

OMEGA_R0_BASE = 9.0e-5
OMEGA_B0_BASE = 0.05
Z_STAR = 1100.0
A_STAR = 1.0 / (1.0 + Z_STAR)

# Planck-like acoustic angular scale: theta_* ~ 0.01041 rad, ell_A ~ 301.5
TARGET_THETA = 0.01041
TARGET_ELL_A = math.pi / TARGET_THETA

N_START = -10.0
N_END = 0.0
N_GRID = 2500

OUT_DIR = Path("/mnt/data/flux_sound_horizon_outputs")
CSV_PATH = OUT_DIR / "flux_sound_horizon_scan.csv"


@dataclass(frozen=True)
class FluxParams:
    """Parameter set for one flux-scalar CMB run."""

    Omega_m0: float = 0.30
    Omega_r0: float = OMEGA_R0_BASE
    Omega_b0: float = OMEGA_B0_BASE
    V0: float = 0.70
    mu: float = 2.0
    sigma_init: float = 5.0
    v_init: float = 0.0


@dataclass
class SoundHorizonResult:
    """CMB boundary diagnostics for one run."""

    Omega_m0: float
    Omega_r0: float
    Omega_b0: float
    V0: float
    mu: float
    sigma_init: float
    v_init: float
    success: bool
    r_s: float
    chi_star: float
    theta_star: float
    ell_A: float
    theta_frac_error: float
    ell_frac_error: float
    H0_model: float
    Omega_sigma_0: float
    w_sigma_0: float
    sigma_today: float
    v_today: float
    early_flux_fraction: float
    notes: str = ""


# =============================================================================
# Flux scalar backend
# =============================================================================


def V_sigma(sigma: np.ndarray | float, V0: float, mu: float) -> np.ndarray | float:
    """Flux/coherence scalar potential."""
    return V0 * (1.0 - np.exp(-np.asarray(sigma) / mu)) ** 2


def dV_dsigma(sigma: float, V0: float, mu: float) -> float:
    """Derivative of the flux potential."""
    exp_term = math.exp(-sigma / mu)
    return 2.0 * V0 * (1.0 - exp_term) * exp_term / mu


def H_of_N(N: np.ndarray | float, sigma: np.ndarray | float, v: np.ndarray | float, p: FluxParams) -> np.ndarray | float:
    """Dimensionless Hubble rate H(N)."""
    a = np.exp(N)
    rho_r = p.Omega_r0 * a ** -4
    rho_m = p.Omega_m0 * a ** -3
    rho_sigma = 0.5 * np.asarray(v) ** 2 + V_sigma(sigma, p.V0, p.mu)
    H2 = rho_r + rho_m + rho_sigma
    return np.sqrt(np.maximum(H2, 0.0))


def flux_scalar_system(N: float, y: np.ndarray, p: FluxParams) -> List[float]:
    """ODE system for sigma and v=d sigma/dt using N=ln(a)."""
    sigma, v = float(y[0]), float(y[1])
    H = float(H_of_N(N, sigma, v, p))
    H_safe = max(H, 1e-14)
    dsigma_dN = v / H_safe
    dv_dN = -3.0 * v - dV_dsigma(sigma, p.V0, p.mu) / H_safe
    return [dsigma_dN, dv_dN]


def integrate_flux_background(p: FluxParams) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Integrate sigma from early universe to today."""
    sol = solve_ivp(
        lambda N, y: flux_scalar_system(N, y, p),
        t_span=(N_START, N_END),
        y0=[p.sigma_init, p.v_init],
        dense_output=True,
        atol=1e-8,
        rtol=1e-8,
        max_step=0.05,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    N_vals = np.linspace(N_START, N_END, N_GRID)
    sigma_vals, v_vals = sol.sol(N_vals)
    H_vals = H_of_N(N_vals, sigma_vals, v_vals, p)
    return N_vals, sigma_vals, v_vals, H_vals


# =============================================================================
# CMB sound-horizon diagnostics
# =============================================================================


def compute_sound_horizon(p: FluxParams) -> SoundHorizonResult:
    """Compute r_s, chi_*, theta_*, and ell_A for one parameter set."""
    try:
        N_vals, sigma_vals, v_vals, H_vals = integrate_flux_background(p)
        a_vals = np.exp(N_vals)

        if not (a_vals[0] < A_STAR < a_vals[-1]):
            raise ValueError("A_STAR lies outside integration range")

        # Radiation in cmbr.py approximates photons as total radiation for c_s.
        rho_b = p.Omega_b0 * a_vals ** -3
        rho_gamma = p.Omega_r0 * a_vals ** -4
        cs2 = 1.0 / (3.0 * (1.0 + 3.0 * rho_b / (4.0 * rho_gamma)))
        c_s = np.sqrt(np.maximum(cs2, 0.0))

        # Interpolate onto dense a-grid that includes a_star exactly.
        a_dense = np.linspace(a_vals[0], 1.0, N_GRID * 2)
        H_dense = np.interp(a_dense, a_vals, H_vals)
        cs_dense = np.interp(a_dense, a_vals, c_s)

        early_mask = a_dense <= A_STAR
        los_mask = a_dense >= A_STAR

        a_early = np.concatenate([a_dense[early_mask], [A_STAR]])
        a_early = np.unique(np.sort(a_early))
        H_early = np.interp(a_early, a_dense, H_dense)
        cs_early = np.interp(a_early, a_dense, cs_dense)

        a_los = np.concatenate([[A_STAR], a_dense[los_mask]])
        a_los = np.unique(np.sort(a_los))
        H_los = np.interp(a_los, a_dense, H_dense)

        r_s = float(np.trapezoid(cs_early / (a_early ** 2 * H_early), a_early))
        chi_star = float(np.trapezoid(1.0 / (a_los ** 2 * H_los), a_los))
        theta_star = r_s / chi_star if chi_star > 0 else float("nan")
        ell_A = math.pi / theta_star if theta_star > 0 else float("nan")

        rho_sigma = 0.5 * v_vals ** 2 + V_sigma(sigma_vals, p.V0, p.mu)
        p_sigma = 0.5 * v_vals ** 2 - V_sigma(sigma_vals, p.V0, p.mu)
        rho_tot = p.Omega_r0 * a_vals ** -4 + p.Omega_m0 * a_vals ** -3 + rho_sigma

        H0_model = float(H_vals[-1])
        Omega_sigma_0 = float(rho_sigma[-1] / rho_tot[-1])
        w_sigma_0 = float(p_sigma[-1] / rho_sigma[-1]) if rho_sigma[-1] != 0 else float("nan")
        early_flux_fraction = float(rho_sigma[np.argmin(np.abs(a_vals - A_STAR))] / rho_tot[np.argmin(np.abs(a_vals - A_STAR))])

        theta_frac_error = (theta_star - TARGET_THETA) / TARGET_THETA
        ell_frac_error = (ell_A - TARGET_ELL_A) / TARGET_ELL_A

        return SoundHorizonResult(
            **asdict(p),
            success=True,
            r_s=r_s,
            chi_star=chi_star,
            theta_star=theta_star,
            ell_A=ell_A,
            theta_frac_error=theta_frac_error,
            ell_frac_error=ell_frac_error,
            H0_model=H0_model,
            Omega_sigma_0=Omega_sigma_0,
            w_sigma_0=w_sigma_0,
            sigma_today=float(sigma_vals[-1]),
            v_today=float(v_vals[-1]),
            early_flux_fraction=early_flux_fraction,
            notes="ok",
        )
    except Exception as exc:
        return SoundHorizonResult(
            **asdict(p),
            success=False,
            r_s=float("nan"),
            chi_star=float("nan"),
            theta_star=float("nan"),
            ell_A=float("nan"),
            theta_frac_error=float("nan"),
            ell_frac_error=float("nan"),
            H0_model=float("nan"),
            Omega_sigma_0=float("nan"),
            w_sigma_0=float("nan"),
            sigma_today=float("nan"),
            v_today=float("nan"),
            early_flux_fraction=float("nan"),
            notes=str(exc),
        )


# =============================================================================
# Scan / plotting
# =============================================================================


def parameter_grid() -> Iterable[FluxParams]:
    """Moderate scan around the current cmbr.py backend."""
    Omega_m_vals = [0.26, 0.28, 0.30, 0.32, 0.34]
    V0_vals = [0.55, 0.65, 0.75, 0.85, 0.95]
    mu_vals = [1.5, 2.0, 2.5, 3.0]
    sigma_vals = [3.0, 4.0, 5.0, 6.0]

    for Om in Omega_m_vals:
        for V0 in V0_vals:
            for mu in mu_vals:
                for sigma_init in sigma_vals:
                    yield FluxParams(
                        Omega_m0=Om,
                        Omega_r0=OMEGA_R0_BASE,
                        Omega_b0=OMEGA_B0_BASE,
                        V0=V0,
                        mu=mu,
                        sigma_init=sigma_init,
                        v_init=0.0,
                    )


def write_csv(results: List[SoundHorizonResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(asdict(results[0]).keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in results:
            writer.writerow(asdict(r))


def plot_results(results: List[SoundHorizonResult]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    good = [r for r in results if r.success and math.isfinite(r.theta_star)]
    if not good:
        return

    theta_err = np.array([100.0 * r.theta_frac_error for r in good])
    ell = np.array([r.ell_A for r in good])
    omega_sigma = np.array([r.Omega_sigma_0 for r in good])
    early_flux = np.array([r.early_flux_fraction for r in good])
    V0 = np.array([r.V0 for r in good])
    mu = np.array([r.mu for r in good])
    sigma_init = np.array([r.sigma_init for r in good])

    plt.figure(figsize=(8, 5))
    plt.scatter(omega_sigma, theta_err, s=18)
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("Omega_sigma(z=0)")
    plt.ylabel("theta_* fractional error (%)")
    plt.title("CMB boundary preservation vs flux density today")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "theta_error_vs_omega_sigma.png", dpi=160)

    plt.figure(figsize=(8, 5))
    plt.scatter(early_flux, ell, s=18)
    plt.axhline(TARGET_ELL_A, linestyle="--", linewidth=1)
    plt.xlabel("Flux scalar fraction at recombination")
    plt.ylabel("ell_A = pi / theta_*")
    plt.title("Acoustic scale vs early flux contamination")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "ellA_vs_early_flux_fraction.png", dpi=160)

    plt.figure(figsize=(8, 5))
    sc = plt.scatter(V0, mu, c=np.abs(theta_err), s=40)
    plt.colorbar(sc, label="|theta_* error| (%)")
    plt.xlabel("V0")
    plt.ylabel("mu")
    plt.title("Flux potential scan: CMB boundary error")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "potential_parameter_error_map.png", dpi=160)

    plt.figure(figsize=(8, 5))
    plt.scatter(sigma_init, theta_err, s=18)
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("sigma_init")
    plt.ylabel("theta_* fractional error (%)")
    plt.title("Initial coherence plateau vs sound-horizon boundary")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "theta_error_vs_sigma_init.png", dpi=160)



def print_summary(results: List[SoundHorizonResult]) -> None:
    good = [r for r in results if r.success and math.isfinite(r.theta_star)]
    ranked = sorted(good, key=lambda r: abs(r.theta_frac_error))
    accepted_1pct = [r for r in ranked if abs(r.theta_frac_error) <= 0.01]
    accepted_5pct = [r for r in ranked if abs(r.theta_frac_error) <= 0.05]

    print("Flux sound horizon scan complete")
    print("================================")
    print(f"Runs attempted              : {len(results)}")
    print(f"Successful runs             : {len(good)}")
    print(f"Target theta_*              : {TARGET_THETA:.6f} rad")
    print(f"Target ell_A                : {TARGET_ELL_A:.3f}")
    print(f"Accepted <= 1% theta error  : {len(accepted_1pct)}")
    print(f"Accepted <= 5% theta error  : {len(accepted_5pct)}")

    if ranked:
        b = ranked[0]
        print("\nBest boundary-preserving parameter set")
        print("--------------------------------------")
        print(f"Omega_m0       : {b.Omega_m0:.3f}")
        print(f"V0             : {b.V0:.3f}")
        print(f"mu             : {b.mu:.3f}")
        print(f"sigma_init     : {b.sigma_init:.3f}")
        print(f"r_s            : {b.r_s:.6e} H0^-1 units")
        print(f"chi_star       : {b.chi_star:.6e} H0^-1 units")
        print(f"theta_*        : {b.theta_star:.6f} rad")
        print(f"ell_A          : {b.ell_A:.3f}")
        print(f"theta error    : {100*b.theta_frac_error:+.3f}%")
        print(f"Omega_sigma_0  : {b.Omega_sigma_0:.3f}")
        print(f"w_sigma_0      : {b.w_sigma_0:.3f}")
        print(f"early flux frac: {b.early_flux_fraction:.3e}")



def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    params = list(parameter_grid())
    results = [compute_sound_horizon(p) for p in params]
    write_csv(results, CSV_PATH)
    plot_results(results)
    print_summary(results)
    print(f"\nCSV written to: {CSV_PATH}")
    print(f"Plots written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
