#!/usr/bin/env python3
"""
flux_de_scalar_bridge.py

Bridge diagnostic for the Flux dark-energy sector.

Purpose
-------
Compare two descriptions of the same proposed physical sector:

1. Scalar-sector background:
   The existing sigma/flux background from cmbr_universe_backend.py, interpreted as
   residual substrate tension / scalar flux energy.

2. Throughput-sector background:
   The leakage-gated dark-energy model from flux_de_best_model.csv, interpreted as
   coherence-gated escaped entropic throughput of structure history.

The test is not a final Boltzmann/CAMB/CLASS validation. It is a repo-level
consistency bridge asking:

    Is the scalar sector behaving like the smooth carrier/floor of the same
    dark-energy budget that the throughput ledger perturbs?

Run
---
    python flux_de_scalar_bridge.py
    python flux_de_scalar_bridge.py --throughput flux_de_best_model.csv --plot

Outputs
-------
    flux_de_scalar_bridge_comparison.csv
    flux_de_scalar_bridge_summary.csv
    flux_de_scalar_bridge.png
"""

from __future__ import annotations

import argparse
import contextlib
import io
import math
import os
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

# Avoid GUI backends when cmbr.py is imported by the backend.
os.environ.setdefault("MPLBACKEND", "Agg")

try:
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:  # pragma: no cover
    HAS_MPL = False

EPS = 1e-30
OMEGA_M0 = 0.30
OMEGA_R0 = 9.0e-5
TARGET_OMEGA_DE0 = 0.70
Z_CMB = 1090.0


def safe_gradient(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    if len(y) < 3:
        return np.zeros_like(y)
    return np.gradient(y, x)


def effective_w_from_rho(a: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Conservation definition: d ln rho / d ln a = -3(1+w)."""
    ln_a = np.log(np.maximum(np.asarray(a, dtype=float), EPS))
    ln_rho = np.log(np.maximum(np.asarray(rho, dtype=float), EPS))
    dlnrho_dlna = safe_gradient(ln_rho, ln_a)
    return -1.0 - (1.0 / 3.0) * dlnrho_dlna


def interp_on_z(source: pd.DataFrame, col: str, z_target: float) -> float:
    tmp = source[["z", col]].replace([np.inf, -np.inf], np.nan).dropna().sort_values("z")
    if tmp.empty:
        return float("nan")
    return float(np.interp(z_target, tmp["z"].to_numpy(float), tmp[col].to_numpy(float)))


def load_scalar_background(n_steps: int = 8000, n_start: float = -10.0, n_end: float = 0.0) -> pd.DataFrame:
    """
    Load scalar background using cmbr_universe_backend.py.

    The backend returns flux_frac as Omega_sigma. We reconstruct physical
    rho_sigma from H^2 * Omega_sigma and infer w_sigma via energy conservation.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))

    try:
        # cmbr.py has top-level print/plot code in some repo states; suppress import noise.
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            from cmbr_universe_backend import integrate_background
            bg = integrate_background(N_start=n_start, N_end=n_end, n_steps=n_steps)

        a = np.asarray(bg["a"], dtype=float)
        t = np.asarray(bg["t"], dtype=float)
        H = np.asarray(bg["H"], dtype=float)
        sigma = np.asarray(bg.get("sigma", np.full_like(a, np.nan)), dtype=float)
        Omega_sigma = np.asarray(bg.get("flux_frac", np.full_like(a, np.nan)), dtype=float)
        source = "cmbr_universe_backend.py"
    except Exception:
        # Fallback should rarely be used, but keeps the script runnable in a clean clone.
        N = np.linspace(n_start, n_end, n_steps)
        a = np.exp(N)
        H = np.sqrt(OMEGA_R0 * a**-4 + OMEGA_M0 * a**-3 + TARGET_OMEGA_DE0)
        t = np.zeros_like(a)
        for i in range(1, len(a)):
            dN = N[i] - N[i - 1]
            t[i] = t[i - 1] + dN / max(0.5 * (H[i] + H[i - 1]), EPS)
        sigma = np.full_like(a, np.nan)
        rho_sigma_fallback = np.full_like(a, TARGET_OMEGA_DE0)
        rho_tot = OMEGA_R0 * a**-4 + OMEGA_M0 * a**-3 + rho_sigma_fallback
        Omega_sigma = rho_sigma_fallback / rho_tot
        source = "internal_fallback_LCDM_like"

    z = 1.0 / np.maximum(a, EPS) - 1.0
    rho_m = OMEGA_M0 * a**-3
    rho_r = OMEGA_R0 * a**-4

    # In the backend convention H^2 = rho_total, so rho_sigma = Omega_sigma * H^2.
    rho_sigma = np.clip(Omega_sigma, 0.0, None) * H**2
    w_sigma = effective_w_from_rho(a, rho_sigma)

    # Normalize scalar physical density to the same present-day value as throughput.
    rho_sigma_norm = rho_sigma / max(float(rho_sigma[-1]), EPS)
    rho_sigma_scaled = rho_sigma_norm * (TARGET_OMEGA_DE0 / max(1.0 - TARGET_OMEGA_DE0, EPS)) * (rho_m[-1] + rho_r[-1])
    rho_tot_scaled = rho_m + rho_r + rho_sigma_scaled
    Omega_sigma_scaled = rho_sigma_scaled / np.maximum(rho_tot_scaled, EPS)

    out = pd.DataFrame({
        "t": t,
        "a": a,
        "z": z,
        "H_scalar": H,
        "sigma": sigma,
        "rho_m": rho_m,
        "rho_r": rho_r,
        "rho_sigma_raw": rho_sigma,
        "Omega_sigma_raw": Omega_sigma,
        "rho_sigma_scaled": rho_sigma_scaled,
        "Omega_sigma_scaled": Omega_sigma_scaled,
        "w_sigma": w_sigma,
    })
    out.attrs["source"] = source
    return out.sort_values("t").reset_index(drop=True)


def load_throughput(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"t", "a", "z", "rho_DE_flux", "Omega_DE_flux", "w_flux"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Throughput CSV missing required columns: {sorted(missing)}")
    return df.sort_values("t").reset_index(drop=True)


def build_comparison(scalar: pd.DataFrame, throughput: pd.DataFrame) -> pd.DataFrame:
    """Interpolate scalar-sector columns onto throughput redshift grid."""
    # Work in ascending z for interpolation.
    s = scalar.sort_values("z")
    zt = throughput["z"].to_numpy(float)

    comp = throughput.copy()
    for col in ["rho_sigma_scaled", "Omega_sigma_scaled", "w_sigma", "sigma", "Omega_sigma_raw", "rho_sigma_raw"]:
        comp[col] = np.interp(
            zt,
            s["z"].to_numpy(float),
            s[col].replace([np.inf, -np.inf], np.nan).interpolate().bfill().ffill().to_numpy(float),
        )

    # Useful ratios and residuals.
    comp["rho_flux_over_sigma"] = comp["rho_DE_flux"] / np.maximum(comp["rho_sigma_scaled"], EPS)
    comp["Omega_flux_minus_sigma"] = comp["Omega_DE_flux"] - comp["Omega_sigma_scaled"]
    comp["w_flux_minus_sigma"] = comp["w_flux"] - comp["w_sigma"]

    # Present-day floor/exhaust decomposition if available.
    if "throughput_raw" in comp.columns and "rho_DE_flux" in comp.columns:
        # Infer floor as first rho_DE value only if early throughput starts near zero.
        # Safer estimate from high-z plateau: median at z>1000 if present.
        high = comp[comp["z"] > 1000]
        if len(high):
            rho_floor_est = float(np.nanmedian(high["rho_DE_flux"]))
        else:
            rho_floor_est = float(np.nanmin(comp["rho_DE_flux"]))
        comp["rho_floor_est"] = rho_floor_est
        comp["rho_exhaust_est"] = np.maximum(comp["rho_DE_flux"] - rho_floor_est, 0.0)
        comp["exhaust_fraction_of_DE"] = comp["rho_exhaust_est"] / np.maximum(comp["rho_DE_flux"], EPS)

    return comp


def summarize(comp: pd.DataFrame, scalar_source: str) -> pd.DataFrame:
    # Compare on the observationally relevant window, excluding numerical edges.
    late = comp[(comp["z"] >= 0.0) & (comp["z"] <= 5.0)].copy()
    fossil = comp[(comp["z"] >= 1.0) & (comp["z"] <= 3.0)].copy()

    def rms(x: np.ndarray) -> float:
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        if len(x) == 0:
            return float("nan")
        return float(np.sqrt(np.mean(x**2)))

    corr = float("nan")
    if len(late) > 3:
        x = late["rho_DE_flux"].to_numpy(float)
        y = late["rho_sigma_scaled"].to_numpy(float)
        if np.nanstd(x) > 0 and np.nanstd(y) > 0:
            corr = float(np.corrcoef(x, y)[0, 1])

    rows: Dict[str, float | str] = {
        "scalar_source": scalar_source,
        "Omega_flux_0": interp_on_z(comp, "Omega_DE_flux", 0.0),
        "Omega_sigma_scaled_0": interp_on_z(comp, "Omega_sigma_scaled", 0.0),
        "Omega_flux_1090": interp_on_z(comp, "Omega_DE_flux", Z_CMB),
        "Omega_sigma_scaled_1090": interp_on_z(comp, "Omega_sigma_scaled", Z_CMB),
        "w_flux_0": interp_on_z(comp, "w_flux", 0.0),
        "w_sigma_0": interp_on_z(comp, "w_sigma", 0.0),
        "w_flux_z1": interp_on_z(comp, "w_flux", 1.0),
        "w_sigma_z1": interp_on_z(comp, "w_sigma", 1.0),
        "w_flux_z2": interp_on_z(comp, "w_flux", 2.0),
        "w_sigma_z2": interp_on_z(comp, "w_sigma", 2.0),
        "w_flux_z3": interp_on_z(comp, "w_flux", 3.0),
        "w_sigma_z3": interp_on_z(comp, "w_sigma", 3.0),
        "rms_Omega_difference_z0_5": rms(late["Omega_flux_minus_sigma"].to_numpy(float)) if len(late) else float("nan"),
        "rms_w_difference_z0_5": rms(late["w_flux_minus_sigma"].to_numpy(float)) if len(late) else float("nan"),
        "rho_correlation_z0_5": corr,
        "fossil_hump_flux_abs_z1_3": float(np.nanmax(np.abs(fossil["w_flux"].to_numpy(float) + 1.0))) if len(fossil) else float("nan"),
        "fossil_hump_sigma_abs_z1_3": float(np.nanmax(np.abs(fossil["w_sigma"].to_numpy(float) + 1.0))) if len(fossil) else float("nan"),
    }

    if "exhaust_fraction_of_DE" in comp.columns:
        rows["exhaust_fraction_DE_0_est"] = interp_on_z(comp, "exhaust_fraction_of_DE", 0.0)
        rows["exhaust_fraction_DE_z2_est"] = interp_on_z(comp, "exhaust_fraction_of_DE", 2.0)

    # A simple bridge verdict.
    omega_cmb_safe = rows["Omega_flux_1090"] <= 1e-8
    near_lcdm_today = abs(float(rows["w_flux_0"]) + 1.0) <= 0.03
    has_hump = 0.002 <= float(rows["fossil_hump_flux_abs_z1_3"]) <= 0.12
    scalar_smooth = float(rows["fossil_hump_sigma_abs_z1_3"]) < float(rows["fossil_hump_flux_abs_z1_3"]) if np.isfinite(float(rows["fossil_hump_sigma_abs_z1_3"])) else False
    rows["bridge_verdict"] = (
        "PASS_TOY_BRIDGE: throughput behaves as a small fossil correction on a smoother scalar floor"
        if omega_cmb_safe and near_lcdm_today and has_hump and scalar_smooth
        else "CHECK: bridge not yet clean; inspect scalar/throughput mismatch"
    )

    return pd.DataFrame([rows])


def plot_bridge(comp: pd.DataFrame, summary: pd.DataFrame, path: Path) -> None:
    if not HAS_MPL:
        return

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    view = comp[comp["z"] <= 12]
    ax.plot(view["z"], view["Omega_DE_flux"], label="Throughput Ω_DE")
    ax.plot(view["z"], view["Omega_sigma_scaled"], label="Scalar Ω_sigma scaled", linestyle="--")
    ax.set_xlim(12, 0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("Ω")
    ax.set_title("Scalar sector vs throughput sector")
    ax.legend()

    ax = axes[0, 1]
    ax.semilogy(comp["z"], np.maximum(comp["Omega_DE_flux"], EPS), label="Throughput Ω_DE")
    ax.semilogy(comp["z"], np.maximum(comp["Omega_sigma_scaled"], EPS), label="Scalar Ω_sigma scaled", linestyle="--")
    ax.axhline(1e-8, linestyle=":", label="CMB invisibility target")
    ax.set_xlim(1090, 0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("Ω, log scale")
    ax.set_title("CMB safety window")
    ax.legend()

    ax = axes[1, 0]
    view = comp[comp["z"] <= 5]
    ax.plot(view["z"], view["w_flux"], label="w_throughput")
    ax.plot(view["z"], view["w_sigma"], label="w_sigma", linestyle="--")
    ax.axhline(-1.0, linestyle=":", label="Λ")
    ax.set_xlim(5, 0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("w(z)")
    ax.set_title("Equation-of-state comparison")
    ax.legend()

    ax = axes[1, 1]
    view = comp[comp["z"] <= 5]
    ax.plot(view["z"], view["rho_flux_over_sigma"], label="ρ_throughput / ρ_sigma_scaled")
    if "exhaust_fraction_of_DE" in view.columns:
        ax.plot(view["z"], view["exhaust_fraction_of_DE"], label="estimated exhaust fraction")
    ax.axhline(1.0, linestyle=":")
    ax.set_xlim(5, 0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("ratio / fraction")
    ax.set_title("Bridge residuals")
    ax.legend()

    verdict = str(summary["bridge_verdict"].iloc[0]) if "bridge_verdict" in summary.columns else ""
    fig.suptitle(f"Flux DE Scalar Bridge — {verdict}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare Flux scalar DE background against leakage-throughput DE history.")
    p.add_argument("--throughput", default="flux_de_best_model.csv", help="Input throughput best-model CSV.")
    p.add_argument("--comparison-output", default="flux_de_scalar_bridge_comparison.csv")
    p.add_argument("--summary-output", default="flux_de_scalar_bridge_summary.csv")
    p.add_argument("--plot-output", default="flux_de_scalar_bridge.png")
    p.add_argument("--plot", action="store_true", help="Write diagnostic plot.")
    p.add_argument("--n-steps", type=int, default=8000)
    p.add_argument("--top-print", action="store_true", help="Print compact summary.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    here = Path(__file__).resolve().parent
    throughput_path = Path(args.throughput)
    if not throughput_path.is_absolute():
        throughput_path = here / throughput_path

    scalar = load_scalar_background(n_steps=args.n_steps, n_start=-10.0, n_end=0.0)
    throughput = load_throughput(throughput_path)
    comp = build_comparison(scalar, throughput)
    summ = summarize(comp, scalar.attrs.get("source", "unknown"))

    comp_out = Path(args.comparison_output)
    summ_out = Path(args.summary_output)
    plot_out = Path(args.plot_output)
    if not comp_out.is_absolute():
        comp_out = here / comp_out
    if not summ_out.is_absolute():
        summ_out = here / summ_out
    if not plot_out.is_absolute():
        plot_out = here / plot_out

    comp.to_csv(comp_out, index=False)
    summ.to_csv(summ_out, index=False)
    if args.plot:
        plot_bridge(comp, summ, plot_out)

    print("\nFlux-DE scalar bridge")
    print("=" * 96)
    show_cols = [
        "Omega_flux_0", "Omega_flux_1090", "w_flux_0", "w_flux_z2",
        "w_sigma_0", "w_sigma_z2", "rms_w_difference_z0_5",
        "fossil_hump_flux_abs_z1_3", "fossil_hump_sigma_abs_z1_3", "bridge_verdict"
    ]
    print(summ[show_cols].to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print(f"\nWrote: {comp_out.resolve()}")
    print(f"Wrote: {summ_out.resolve()}")
    if args.plot:
        print(f"Wrote: {plot_out.resolve()}")


if __name__ == "__main__":
    main()
