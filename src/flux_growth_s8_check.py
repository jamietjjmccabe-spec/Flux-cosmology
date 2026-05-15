#!/usr/bin/env python3
"""
flux_growth_s8_check.py

Growth-of-structure diagnostic for the Flux-DE Gated Ledger model.

Inputs
------
- flux_de_best_model.csv produced by flux_de_scan.py.

Outputs
-------
- flux_growth_s8_results.csv
- flux_growth_s8_summary.csv
- flux_growth_s8.png

Cases
-----
A. LCDM baseline: constant DE density, mu_eff = 1.
B. Flux-DE drag: throughput-derived rho_DE(z), mu_eff = 1.
C. Full Flux Ledger: throughput-derived rho_DE(z), plus a conservative
   retained-memory mu_eff(z) = 1 + mu_amp * G_ret(z) * G_mu_late(z).

Notes
-----
This is a linear-growth toy diagnostic, not a Boltzmann-code replacement.
Its purpose is to test whether the Flux-DE fossil hump moves sigma8/S8
in the observed low-S8 direction before running CLASS/CAMB-level checks.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False

EPS = 1e-30


@dataclass
class GrowthCase:
    name: str
    label: str
    N: np.ndarray
    a: np.ndarray
    z: np.ndarray
    H: np.ndarray
    Omega_m: np.ndarray
    mu_eff: np.ndarray
    D: np.ndarray
    dD_dN: np.ndarray
    f: np.ndarray
    gamma: np.ndarray
    D_norm: np.ndarray
    sigma8_z: np.ndarray
    fsigma8: np.ndarray
    sigma8_today: float
    S8_today: float


def load_flux_background(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")
    df = pd.read_csv(path)
    required = {"a", "z", "rho_m", "rho_r", "rho_DE_flux", "Omega_DE_flux", "w_flux"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input missing required columns: {sorted(missing)}")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["a", "z", "rho_m", "rho_r", "rho_DE_flux"])
    df = df.sort_values("a").drop_duplicates(subset=["a"], keep="last").reset_index(drop=True)
    return df


def make_grid(df: pd.DataFrame, z_init: float = 100.0, n_grid: int = 2500):
    a_init = 1.0 / (1.0 + z_init)
    a_min = max(a_init, float(df["a"].min()))
    a_max = 1.0
    N = np.linspace(np.log(a_min), np.log(a_max), n_grid)
    a = np.exp(N)
    z = 1.0 / a - 1.0
    return N, a, z


def interp_column(df: pd.DataFrame, a: np.ndarray, col: str) -> np.ndarray:
    xp = df["a"].to_numpy()
    fp = df[col].to_numpy()
    return np.interp(a, xp, fp)


def dlnH_dN(N: np.ndarray, H: np.ndarray) -> np.ndarray:
    lnH = np.log(np.maximum(H, EPS))
    return np.gradient(lnH, N)


def solve_growth_case(
    name: str,
    label: str,
    N: np.ndarray,
    a: np.ndarray,
    z: np.ndarray,
    rho_m: np.ndarray,
    rho_r: np.ndarray,
    rho_de: np.ndarray,
    mu_eff: np.ndarray,
    sigma8_anchor: float,
    sigma8_ratio_to_lcdm: float | None = None,
    delta_initial: float | None = None,
):
    rho_tot = rho_m + rho_r + rho_de
    H = np.sqrt(np.maximum(rho_tot, EPS))
    Om = rho_m / np.maximum(rho_tot, EPS)
    dlnH = dlnH_dN(N, H)

    if delta_initial is None:
        delta_initial = a[0]
    y0 = [delta_initial, delta_initial]  # matter-era growing mode: D ~ a, D' ~ a

    def rhs(Nx, y):
        D, Dp = y
        dlnH_x = np.interp(Nx, N, dlnH)
        Om_x = np.interp(Nx, N, Om)
        mu_x = np.interp(Nx, N, mu_eff)
        Dpp = -(2.0 + dlnH_x) * Dp + 1.5 * mu_x * Om_x * D
        return [Dp, Dpp]

    sol = solve_ivp(rhs, (N[0], N[-1]), y0, t_eval=N, atol=1e-10, rtol=1e-9)
    if not sol.success:
        raise RuntimeError(f"Growth solve failed for {name}: {sol.message}")

    D = sol.y[0]
    Dp = sol.y[1]
    f = Dp / np.maximum(D, EPS)
    valid = (f > 0) & (Om > 1e-4) & (Om < 0.9999)
    gamma = np.full_like(f, np.nan)
    gamma[valid] = np.log(f[valid]) / np.log(Om[valid])
    D_norm = D / D[-1]

    if sigma8_ratio_to_lcdm is None:
        sigma8_today = sigma8_anchor
    else:
        sigma8_today = sigma8_anchor * sigma8_ratio_to_lcdm
    Omega_m0 = Om[-1]
    S8_today = sigma8_today * np.sqrt(Omega_m0 / 0.3)
    sigma8_z = sigma8_today * D_norm
    fsigma8 = f * sigma8_z

    return GrowthCase(
        name=name,
        label=label,
        N=N,
        a=a,
        z=z,
        H=H,
        Omega_m=Om,
        mu_eff=mu_eff,
        D=D,
        dD_dN=Dp,
        f=f,
        gamma=gamma,
        D_norm=D_norm,
        sigma8_z=sigma8_z,
        fsigma8=fsigma8,
        sigma8_today=float(sigma8_today),
        S8_today=float(S8_today),
    )


def value_at_z(case: GrowthCase, arr: np.ndarray, z_target: float) -> float:
    idx = int(np.argmin(np.abs(case.z - z_target)))
    return float(arr[idx])


def build_mu_eff_full(df: pd.DataFrame, a: np.ndarray, mu_amp: float, mu_z_damp: float) -> np.ndarray:
    """
    Conservative retained-memory growth modifier.

    G_ret is the complement of the escape gate. The optional damping term prevents
    this toy halo-memory correction from acting as a large homogeneous early-universe
    modification. This keeps Case C as a bounded late-time diagnostic rather than a
    replacement for a full nonlinear halo model.
    """
    if "G_esc" in df.columns:
        G_esc = interp_column(df, a, "G_esc")
        G_ret = np.clip(1.0 - G_esc, 0.0, 1.0)
    elif "C_sigma" in df.columns:
        G_ret = np.clip(interp_column(df, a, "C_sigma"), 0.0, 1.0)
    else:
        G_ret = np.zeros_like(a)
    z = 1.0 / a - 1.0
    G_mu_late = 1.0 / (1.0 + (z / max(mu_z_damp, EPS))**4)
    return 1.0 + mu_amp * G_ret * G_mu_late


def plot_cases(cases: list[GrowthCase], outpath: Path):
    if not HAS_MPL:
        return
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    for c in cases:
        ax.plot(c.z, c.D_norm, label=c.label)
    ax.set_xlim(5, 0)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("D(z) / D(0)")
    ax.set_title("Normalized linear growth")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[0, 1]
    for c in cases:
        ax.plot(c.z, c.fsigma8, label=c.label)
    ax.set_xlim(3, 0)
    ax.set_xlabel("redshift z")
    ax.set_ylabel(r"$f\sigma_8(z)$")
    ax.set_title("RSD observable")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[1, 0]
    for c in cases:
        ax.plot(c.z, c.gamma, label=c.label)
    ax.set_xlim(3, 0)
    ax.set_ylim(0.35, 0.75)
    ax.set_xlabel("redshift z")
    ax.set_ylabel(r"growth index $\gamma$")
    ax.set_title(r"$f \approx \Omega_m^\gamma$")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[1, 1]
    labels = [c.name for c in cases]
    s8 = [c.S8_today for c in cases]
    ax.bar(labels, s8)
    ax.axhspan(0.77, 0.80, alpha=0.15, label="target low-S8 band")
    ax.set_ylabel(r"$S_8$ today")
    ax.set_title("CMB-anchored S8 comparison")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()

    fig.suptitle("Flux Gated-Ledger Growth Diagnostic", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="Linear growth / S8 diagnostic for Flux-DE Gated Ledger.")
    p.add_argument("--input", default="/mnt/data/flux_de_best_model.csv")
    p.add_argument("--outdir", default="/mnt/data")
    p.add_argument("--z-init", type=float, default=100.0)
    p.add_argument("--n-grid", type=int, default=2500)
    p.add_argument("--sigma8-anchor", type=float, default=0.811, help="CMB-anchored LCDM sigma8 reference.")
    p.add_argument("--mu-amp", type=float, default=0.03, help="Conservative retained-memory mu_eff amplitude for Case C.")
    p.add_argument("--mu-z-damp", type=float, default=2.0, help="Late-time damping scale for Case C toy mu_eff.")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = load_flux_background(Path(args.input))
    N, a, z = make_grid(df, z_init=args.z_init, n_grid=args.n_grid)

    rho_m = interp_column(df, a, "rho_m")
    rho_r = interp_column(df, a, "rho_r")
    rho_de_flux = interp_column(df, a, "rho_DE_flux")
    rho_de0 = float(rho_de_flux[-1])
    rho_de_lcdm = np.full_like(a, rho_de0)

    mu_A = np.ones_like(a)
    mu_B = np.ones_like(a)
    mu_C = build_mu_eff_full(df, a, mu_amp=args.mu_amp, mu_z_damp=args.mu_z_damp)

    # Solve A first, then anchor B/C to same initial amplitude and compare final D.
    case_A = solve_growth_case("A", "A: LCDM baseline", N, a, z, rho_m, rho_r, rho_de_lcdm, mu_A, args.sigma8_anchor)
    D_A_today = case_A.D[-1]
    case_B_raw = solve_growth_case("B_raw", "B raw", N, a, z, rho_m, rho_r, rho_de_flux, mu_B, args.sigma8_anchor)
    ratio_B = case_B_raw.D[-1] / D_A_today
    case_B = solve_growth_case("B", "B: Flux-DE drag", N, a, z, rho_m, rho_r, rho_de_flux, mu_B, args.sigma8_anchor, sigma8_ratio_to_lcdm=ratio_B)
    case_C_raw = solve_growth_case("C_raw", "C raw", N, a, z, rho_m, rho_r, rho_de_flux, mu_C, args.sigma8_anchor)
    ratio_C = case_C_raw.D[-1] / D_A_today
    case_C = solve_growth_case("C", "C: Full Flux ledger", N, a, z, rho_m, rho_r, rho_de_flux, mu_C, args.sigma8_anchor, sigma8_ratio_to_lcdm=ratio_C)
    cases = [case_A, case_B, case_C]

    rows = []
    for c in cases:
        rows.append(pd.DataFrame({
            "case": c.name,
            "label": c.label,
            "N": c.N,
            "a": c.a,
            "z": c.z,
            "H": c.H,
            "Omega_m": c.Omega_m,
            "mu_eff": c.mu_eff,
            "D_raw": c.D,
            "D_norm": c.D_norm,
            "f": c.f,
            "gamma": c.gamma,
            "sigma8_z": c.sigma8_z,
            "fsigma8": c.fsigma8,
        }))
    results = pd.concat(rows, ignore_index=True)
    results_path = outdir / "flux_growth_s8_results.csv"
    results.to_csv(results_path, index=False)

    summary_rows = []
    for c in cases:
        summary_rows.append({
            "case": c.name,
            "label": c.label,
            "sigma8_today": c.sigma8_today,
            "S8_today": c.S8_today,
            "delta_S8_vs_LCDM_percent": 100.0 * (c.S8_today / case_A.S8_today - 1.0),
            "D_today_raw_ratio_vs_LCDM": c.D[-1] / case_A.D[-1],
            "Omega_m0": c.Omega_m[-1],
            "mu_eff_z0": value_at_z(c, c.mu_eff, 0.0),
            "mu_eff_z1": value_at_z(c, c.mu_eff, 1.0),
            "mu_eff_z2": value_at_z(c, c.mu_eff, 2.0),
            "f_z0": value_at_z(c, c.f, 0.0),
            "f_z0p5": value_at_z(c, c.f, 0.5),
            "f_z1": value_at_z(c, c.f, 1.0),
            "fsigma8_z0": value_at_z(c, c.fsigma8, 0.0),
            "fsigma8_z0p5": value_at_z(c, c.fsigma8, 0.5),
            "fsigma8_z1": value_at_z(c, c.fsigma8, 1.0),
            "gamma_z0": value_at_z(c, c.gamma, 0.0),
            "gamma_z0p5": value_at_z(c, c.gamma, 0.5),
            "gamma_z1": value_at_z(c, c.gamma, 1.0),
        })
    summary = pd.DataFrame(summary_rows)

    # Add a simple verdict row as metadata columns repeated in a sidecar-friendly way.
    B_supp = float(summary.loc[summary.case == "B", "delta_S8_vs_LCDM_percent"].iloc[0])
    C_supp = float(summary.loc[summary.case == "C", "delta_S8_vs_LCDM_percent"].iloc[0])
    if B_supp < 0 and C_supp < 0:
        verdict = "PASS_DIRECTION: Flux-DE drag suppresses S8; full ledger remains below LCDM."
    elif B_supp < 0 <= C_supp:
        verdict = "MIXED: Flux-DE drag suppresses S8, but retained-memory mu_eff over-restores clustering."
    else:
        verdict = "FAIL_DIRECTION: background drag does not lower S8 in this toy setup."
    summary["diagnostic_verdict"] = verdict
    summary["input_file"] = str(Path(args.input))
    summary["z_init"] = args.z_init
    summary["sigma8_anchor_LCDM"] = args.sigma8_anchor
    summary["case_C_mu_amp"] = args.mu_amp
    summary["case_C_mu_z_damp"] = args.mu_z_damp

    summary_path = outdir / "flux_growth_s8_summary.csv"
    summary.to_csv(summary_path, index=False)

    plot_path = outdir / "flux_growth_s8.png"
    plot_cases(cases, plot_path)

    with pd.option_context("display.max_columns", None, "display.width", 220):
        print("\nFlux growth / S8 diagnostic")
        print("=" * 88)
        print(summary[["case", "sigma8_today", "S8_today", "delta_S8_vs_LCDM_percent", "fsigma8_z0p5", "fsigma8_z1", "gamma_z0", "gamma_z1", "diagnostic_verdict"]].to_string(index=False))
        print(f"\nWrote: {results_path}")
        print(f"Wrote: {summary_path}")
        if HAS_MPL:
            print(f"Wrote: {plot_path}")


if __name__ == "__main__":
    main()
