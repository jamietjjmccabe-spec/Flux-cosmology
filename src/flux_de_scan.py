#!/usr/bin/env python3
"""
flux_de_scan.py

Flux Cosmology dark-energy throughput scanner.

Purpose
-------
Search for parameter islands where a leakage-gated, structure-history source can
produce a late-time dark-energy-like sector while remaining invisible at CMB time.

Core model
----------
    rho_DE(t) = rho_floor + lambda0 * ∫ G_esc(t') G_leak(t') Sdot_struct(t') W(t,t') dt'

with:
    G_esc  = (1 - C_sigma)^n
    G_leak = 1 / [1 + ((1+z)/(1+z_leak))^m]

The script auto-normalizes lambda0 for each parameter point so that
Omega_DE(z=0) ~= target_omega_de0, then tests:
    - Omega_DE(z=0)
    - Omega_DE(z=1090)
    - w0
    - structure-driven w(z) hump in z=1..3

Notes
-----
This is a toy-to-intermediate validation module, not a final Boltzmann solver.
It uses a background time base from cmbr_universe_backend.py when available,
then builds a normalized throughput history on that time grid.

Run:
    python flux_de_scan.py
    python flux_de_scan.py --plot --top 20
    python flux_de_scan.py --z-leak 6 8 10 12 --m-values 4 6 8 --n-values 1 2 3 4

Outputs:
    flux_de_scan_results.csv
    flux_de_best_model.csv
    flux_de_best_model.png   (with --plot)
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Dict, List, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False

try:
    from scipy.integrate import cumulative_trapezoid
except Exception:  # pragma: no cover
    cumulative_trapezoid = None

EPS = 1e-30

# Baseline density units: H0^2 units with 8πG/3 = 1, matching the existing toy background.
OMEGA_M0 = 0.30
OMEGA_R0 = 9.0e-5
TARGET_OMEGA_DE0 = 0.70
Z_CMB = 1090.0


@dataclass
class ScanParams:
    z_leak: float
    m: float
    n: float
    rho_floor_frac: float
    lambda0: float


# -----------------------------------------------------------------------------
# Background time base
# -----------------------------------------------------------------------------

def load_background(n_steps: int = 8000, n_start: float = -10.0, n_end: float = 0.0) -> pd.DataFrame:
    """
    Load the project background if available. Falls back to an internal ΛCDM-like
    background with the same density convention.

    Returns a DataFrame sorted by increasing cosmic time, i.e. early -> late.
    Columns: t, a, z, H, sigma, C_sigma_base, rho_m, rho_r
    """
    try:
        # Works when this script is run in the repo beside cmbr_universe_backend.py.
        from cmbr_universe_backend import integrate_background
        bg = integrate_background(N_start=n_start, N_end=n_end, n_steps=n_steps)
        a = np.asarray(bg["a"], dtype=float)
        H = np.asarray(bg["H"], dtype=float)
        t = np.asarray(bg["t"], dtype=float)
        sigma = np.asarray(bg.get("sigma", np.ones_like(a)), dtype=float)
    except Exception:
        # Internal fallback: integrate t(N)=∫dN/H for a matter/radiation/Λ background.
        N = np.linspace(n_start, n_end, n_steps)
        a = np.exp(N)
        H = np.sqrt(OMEGA_R0 * a**-4 + OMEGA_M0 * a**-3 + TARGET_OMEGA_DE0)
        t = np.zeros_like(a)
        for i in range(1, len(a)):
            dN = N[i] - N[i - 1]
            t[i] = t[i - 1] + dN / max(0.5 * (H[i] + H[i - 1]), EPS)
        sigma = np.ones_like(a)

    z = 1.0 / np.maximum(a, EPS) - 1.0
    rho_m = OMEGA_M0 * a**-3
    rho_r = OMEGA_R0 * a**-4

    # Normalize sigma into a bounded coherence proxy.
    # This is intentionally conservative: high scalar order -> higher retention.
    smin = float(np.nanmin(sigma))
    smax = float(np.nanmax(sigma))
    if abs(smax - smin) < 1e-12:
        C_sigma_base = np.full_like(sigma, 0.5)
    else:
        C_sigma_base = (sigma - smin) / (smax - smin)
    C_sigma_base = np.clip(C_sigma_base, 0.0, 1.0)

    out = pd.DataFrame({
        "t": t,
        "a": a,
        "z": z,
        "H": H,
        "sigma": sigma,
        "C_sigma_base": C_sigma_base,
        "rho_m": rho_m,
        "rho_r": rho_r,
    })
    return out.sort_values("t").reset_index(drop=True)


# -----------------------------------------------------------------------------
# Source and gates
# -----------------------------------------------------------------------------

def madau_dickinson_sfr(z: np.ndarray) -> np.ndarray:
    """Cosmic star-formation-rate density proxy from Madau-Dickinson form."""
    z = np.asarray(z, dtype=float)
    return 0.015 * ((1.0 + z) ** 2.7) / (1.0 + ((1.0 + z) / 2.9) ** 5.6)


def structural_entropy_source(z: np.ndarray, bh_weight: float = 0.20, shock_weight: float = 0.15) -> np.ndarray:
    """
    Toy structural entropy source.

    SFR is the backbone. Extra terms are deliberately simple proxies:
      - BH/AGN channel: delayed/smoothed high-SFR contribution.
      - Shock/merger channel: nonlinear enhancement around peak structure activity.
    """
    sfr = madau_dickinson_sfr(z)
    sfr_norm = sfr / max(float(np.nanmax(sfr)), EPS)
    bh_proxy = bh_weight * sfr_norm * np.exp(-((z - 2.0) / 1.4) ** 2)
    shock_proxy = shock_weight * sfr_norm ** 1.25
    source = sfr_norm + bh_proxy + shock_proxy
    return np.clip(source, 0.0, None)


def leakage_gate(z: np.ndarray, z_leak: float = 8.0, m: float = 6.0) -> np.ndarray:
    """
    Smooth activation of large-scale entropy leakage.

    z_leak is a structure-maturity redshift, observationally anchored near
    reionization but not identical to ionization itself.
    """
    z = np.asarray(z, dtype=float)
    return 1.0 / (1.0 + ((1.0 + z) / (1.0 + z_leak)) ** m)


def coherence_profile(bg: pd.DataFrame, mode: str = "backend", C_floor: float = 0.05, C_ceiling: float = 0.95) -> np.ndarray:
    """
    Build C_sigma(z), bounded in [0,1].

    mode='backend' uses the normalized sigma from the existing background.
    mode='igm_low' treats the unbound substrate as low coherence at late times.
    mode='hybrid' blends both.
    """
    z = bg["z"].to_numpy(float)
    base = bg["C_sigma_base"].to_numpy(float)

    if mode == "backend":
        C = base
    elif mode == "igm_low":
        # Conservative toy: global unbound medium becomes lower-coherence after structure matures.
        C = C_floor + (C_ceiling - C_floor) * (1.0 / (1.0 + ((1.0 + 8.0) / (1.0 + z)) ** 3.0))
    elif mode == "hybrid":
        igm = C_floor + (C_ceiling - C_floor) * (1.0 / (1.0 + ((1.0 + 8.0) / (1.0 + z)) ** 3.0))
        C = 0.5 * base + 0.5 * igm
    else:
        raise ValueError(f"Unknown coherence mode: {mode}")

    return np.clip(C, 0.0, 1.0)


def escape_gate(C_sigma: np.ndarray, n: float = 2.0) -> np.ndarray:
    return np.clip((1.0 - np.asarray(C_sigma, dtype=float)) ** n, 0.0, 1.0)


# -----------------------------------------------------------------------------
# Throughput integration and diagnostics
# -----------------------------------------------------------------------------

def cumulative_time_integral(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Cumulative ∫ y dt on an increasing-time grid."""
    if cumulative_trapezoid is not None:
        return cumulative_trapezoid(y, t, initial=0.0)
    out = np.zeros_like(y, dtype=float)
    for i in range(1, len(y)):
        out[i] = out[i - 1] + 0.5 * (y[i] + y[i - 1]) * (t[i] - t[i - 1])
    return out


def compute_model(
    bg: pd.DataFrame,
    z_leak: float,
    m: float,
    n: float,
    rho_floor_frac: float,
    target_omega_de0: float = TARGET_OMEGA_DE0,
    coherence_mode: str = "hybrid",
    memory_tau: float | None = None,
) -> Tuple[pd.DataFrame, ScanParams]:
    """
    Compute one gated-throughput DE model.

    rho_floor_frac is a fraction of the target present-day physical DE density.
    lambda0 is solved so rho_DE(a=1) hits the target Omega_DE0 budget.
    """
    df = bg.copy()
    t = df["t"].to_numpy(float)
    z = df["z"].to_numpy(float)
    a = df["a"].to_numpy(float)

    C = coherence_profile(df, mode=coherence_mode)
    Gesc = escape_gate(C, n=n)
    Gleak = leakage_gate(z, z_leak=z_leak, m=m)
    Sdot = structural_entropy_source(z)
    raw_source = Gesc * Gleak * Sdot

    # Optional finite memory kernel. The default is permanent accumulation, matching the ledger concept.
    if memory_tau is None or memory_tau <= 0:
        throughput_raw = cumulative_time_integral(t, raw_source)
    else:
        throughput_raw = np.zeros_like(raw_source)
        for i in range(1, len(raw_source)):
            dt = t[i] - t[i - 1]
            decay = math.exp(-dt / memory_tau)
            throughput_raw[i] = throughput_raw[i - 1] * decay + 0.5 * (raw_source[i] + raw_source[i - 1]) * dt

    # Target physical rho_DE0 such that Omega_DE0 = target against present matter+radiation.
    rho_m0 = float(df["rho_m"].iloc[-1])
    rho_r0 = float(df["rho_r"].iloc[-1])
    rho_de0_target = target_omega_de0 / max(1.0 - target_omega_de0, EPS) * (rho_m0 + rho_r0)
    rho_floor = rho_floor_frac * rho_de0_target

    final_throughput = float(throughput_raw[-1])
    if final_throughput <= EPS:
        lambda0 = np.nan
        rho_de = np.full_like(throughput_raw, np.nan)
    else:
        lambda0 = max((rho_de0_target - rho_floor) / final_throughput, 0.0)
        rho_de = rho_floor + lambda0 * throughput_raw

    rho_tot = df["rho_m"].to_numpy(float) + df["rho_r"].to_numpy(float) + rho_de
    Omega_de = rho_de / np.maximum(rho_tot, EPS)

    # Effective w from conservation: dlnrho/dlna = -3(1+w).
    ln_a = np.log(np.maximum(a, EPS))
    ln_rho = np.log(np.maximum(rho_de, EPS))
    dlnrho_dlna = np.gradient(ln_rho, ln_a)
    w_eff = -1.0 - (1.0 / 3.0) * dlnrho_dlna

    df["C_sigma"] = C
    df["G_esc"] = Gesc
    df["G_leak"] = Gleak
    df["Sdot_struct"] = Sdot
    df["source_gated"] = raw_source
    df["throughput_raw"] = throughput_raw
    df["rho_DE_flux"] = rho_de
    df["Omega_DE_flux"] = Omega_de
    df["w_flux"] = w_eff

    params = ScanParams(z_leak=z_leak, m=m, n=n, rho_floor_frac=rho_floor_frac, lambda0=lambda0)
    return df, params


def interp_at_z(df: pd.DataFrame, col: str, z_target: float) -> float:
    # df is early->late, z descends. Sort ascending z for interpolation.
    tmp = df[["z", col]].dropna().sort_values("z")
    return float(np.interp(z_target, tmp["z"], tmp[col]))


def summarize_model(df: pd.DataFrame, params: ScanParams) -> Dict[str, float]:
    w0 = float(df["w_flux"].iloc[-1])
    omega0 = float(df["Omega_DE_flux"].iloc[-1])
    omega_cmb = interp_at_z(df, "Omega_DE_flux", Z_CMB)
    w_z1 = interp_at_z(df, "w_flux", 1.0)
    w_z2 = interp_at_z(df, "w_flux", 2.0)
    w_z3 = interp_at_z(df, "w_flux", 3.0)
    # Hump = max absolute deviation from -1 in the fossil window, excluding pathological edges.
    window = df[(df["z"] >= 1.0) & (df["z"] <= 3.0)]
    hump_abs = float(np.nanmax(np.abs(window["w_flux"].to_numpy(float) + 1.0))) if len(window) else np.nan
    hump_signed_at_z2 = w_z2 + 1.0

    # Scoring: low is better, with a small reward for measurable-but-not-huge hump.
    omega0_err = abs(omega0 - TARGET_OMEGA_DE0)
    cmb_penalty = max(math.log10(max(omega_cmb, EPS)) - math.log10(1e-8), 0.0)
    w0_err = abs(w0 + 1.0)
    # Prefer subtle hump roughly 0.005..0.08; penalize none or too much.
    if hump_abs < 0.002:
        hump_penalty = 1.0
    elif hump_abs > 0.12:
        hump_penalty = 5.0 * (hump_abs - 0.12)
    else:
        hump_penalty = 0.0

    score = 10.0 * omega0_err + 2.0 * cmb_penalty + 10.0 * w0_err + hump_penalty

    pass_basic = (
        abs(omega0 - TARGET_OMEGA_DE0) <= 0.03
        and omega_cmb <= 1e-8
        and abs(w0 + 1.0) <= 0.03
        and 0.002 <= hump_abs <= 0.12
    )

    return {
        "z_leak": params.z_leak,
        "m": params.m,
        "n": params.n,
        "rho_floor_frac": params.rho_floor_frac,
        "lambda0": params.lambda0,
        "Omega_DE_0": omega0,
        "Omega_DE_1090": omega_cmb,
        "w0": w0,
        "w_z1": w_z1,
        "w_z2": w_z2,
        "w_z3": w_z3,
        "hump_abs_z1_3": hump_abs,
        "hump_signed_z2": hump_signed_at_z2,
        "score": score,
        "pass_basic": bool(pass_basic),
    }


def run_scan(
    bg: pd.DataFrame,
    z_leak_values: Iterable[float],
    m_values: Iterable[float],
    n_values: Iterable[float],
    rho_floor_fracs: Iterable[float],
    coherence_mode: str,
    memory_tau: float | None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows: List[Dict[str, float]] = []
    best_df = None
    best_score = float("inf")

    for z_leak in z_leak_values:
        for m in m_values:
            for n in n_values:
                for floor in rho_floor_fracs:
                    model, params = compute_model(
                        bg,
                        z_leak=float(z_leak),
                        m=float(m),
                        n=float(n),
                        rho_floor_frac=float(floor),
                        coherence_mode=coherence_mode,
                        memory_tau=memory_tau,
                    )
                    summary = summarize_model(model, params)
                    rows.append(summary)
                    if np.isfinite(summary["score"]) and summary["score"] < best_score:
                        best_score = summary["score"]
                        best_df = model.copy()

    results = pd.DataFrame(rows).sort_values(["pass_basic", "score"], ascending=[False, True]).reset_index(drop=True)
    if best_df is None:
        best_df = pd.DataFrame()
    return results, best_df


# -----------------------------------------------------------------------------
# Plotting and CLI
# -----------------------------------------------------------------------------

def plot_best(best: pd.DataFrame, path: Path) -> None:
    if not HAS_MPL or best.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(best["z"], best["G_leak"], label="G_leak")
    ax.plot(best["z"], best["G_esc"], label="G_esc")
    ax.set_xlim(12, 0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("Gate value")
    ax.set_title("Leakage and escape gates")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(best["z"], best["Sdot_struct"], label="Sdot_struct")
    ax.plot(best["z"], best["source_gated"], label="gated source")
    ax.set_xlim(12, 0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("Normalized source")
    ax.set_title("Structural entropy source")
    ax.legend()

    ax = axes[1, 0]
    ax.semilogy(best["z"], np.maximum(best["Omega_DE_flux"], EPS))
    ax.axhline(1e-8, linestyle="--", label="CMB invisibility target")
    ax.set_xlim(1090, 0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("Omega_DE_flux")
    ax.set_title("DE fraction from CMB to today")
    ax.legend()

    ax = axes[1, 1]
    view = best[best["z"] <= 5]
    ax.plot(view["z"], view["w_flux"], label="w_flux")
    ax.axhline(-1.0, linestyle="--", label="Lambda")
    ax.set_xlim(5, 0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("w(z)")
    ax.set_title("Equation-of-state fossil hump")
    ax.legend()

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Scan leakage-gated Flux dark-energy throughput models.")
    p.add_argument("--output", default="flux_de_scan_results.csv", help="Scan CSV output path.")
    p.add_argument("--best-output", default="flux_de_best_model.csv", help="Best model time-series CSV output path.")
    p.add_argument("--plot-output", default="flux_de_best_model.png", help="Best model plot path.")
    p.add_argument("--plot", action="store_true", help="Write diagnostic plot for best model.")
    p.add_argument("--top", type=int, default=15, help="Rows to print.")
    p.add_argument("--z-leak", nargs="*", type=float, default=[6, 7, 8, 9, 10, 11, 12])
    p.add_argument("--m-values", nargs="*", type=float, default=[4, 5, 6, 7, 8])
    p.add_argument("--n-values", nargs="*", type=float, default=[1, 2, 3, 4])
    p.add_argument("--floor-fracs", nargs="*", type=float, default=[0.0, 0.02, 0.05, 0.10, 0.20, 0.40, 0.70, 0.90])
    p.add_argument("--coherence-mode", choices=["backend", "igm_low", "hybrid"], default="hybrid")
    p.add_argument("--memory-tau", type=float, default=0.0, help="Optional finite memory decay time in backend units. 0=permanent ledger.")
    p.add_argument("--n-steps", type=int, default=8000)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    bg = load_background(n_steps=args.n_steps, n_start=-10.0, n_end=0.0)
    memory_tau = None if args.memory_tau <= 0 else args.memory_tau

    results, best = run_scan(
        bg,
        z_leak_values=args.z_leak,
        m_values=args.m_values,
        n_values=args.n_values,
        rho_floor_fracs=args.floor_fracs,
        coherence_mode=args.coherence_mode,
        memory_tau=memory_tau,
    )

    out = Path(args.output)
    best_out = Path(args.best_output)
    results.to_csv(out, index=False)
    best.to_csv(best_out, index=False)

    if args.plot:
        plot_best(best, Path(args.plot_output))

    cols = [
        "pass_basic", "score", "z_leak", "m", "n", "rho_floor_frac", "lambda0",
        "Omega_DE_0", "Omega_DE_1090", "w0", "w_z1", "w_z2", "w_z3", "hump_abs_z1_3"
    ]
    print("\nFlux-DE leakage-gated scan")
    print("=" * 96)
    print(results[cols].head(args.top).to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print(f"\nWrote: {out.resolve()}")
    print(f"Wrote: {best_out.resolve()}")
    if args.plot:
        print(f"Wrote: {Path(args.plot_output).resolve()}")


if __name__ == "__main__":
    main()
