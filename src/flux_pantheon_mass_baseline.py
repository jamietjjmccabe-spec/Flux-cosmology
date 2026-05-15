#!/usr/bin/env python3
"""
flux_pantheon_mass_baseline.py

Pantheon+/SH0ES host-mass baseline audit for Flux/MIP Level-3A.

Purpose
-------
Before the full Flux capacity-overflow kill test can be run, we need a control
baseline using the columns already present in Pantheon+SH0ES.dat:

    CID / SNID
    zHD
    MU_SH0ES
    HOST_LOGMASS

This script tests the standard host-mass relationship:

    Delta_mu ~ HOST_LOGMASS

and builds the mass-step baseline that the later Flux F_fill / G_sat model must beat.

Outputs
-------
flux_pantheon_mass_baseline_data.csv
flux_pantheon_mass_baseline_results.csv
flux_pantheon_mass_baseline_summary.csv
flux_pantheon_mass_baseline_binned.csv
flux_pantheon_mass_baseline.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EPS = 1e-12
C_LIGHT = 299792.458


def read_pantheon_dat(path: str | Path) -> pd.DataFrame:
    """
    Read Pantheon+SH0ES.dat style whitespace table.
    """
    return pd.read_csv(path, delim_whitespace=True, comment="#")


def lcdm_mu(z, H0=67.4, omega_m=0.3):
    """
    Flat LCDM distance modulus relative to H0 and Omega_m.
    Good enough for a host-mass baseline residual construction.

    D_L = (1+z)c/H0 int dz/E(z)
    mu = 5log10(D_L/Mpc)+25
    """
    z = np.asarray(z, dtype=float)
    out = np.full_like(z, np.nan, dtype=float)

    for i, zz in enumerate(z):
        if not np.isfinite(zz) or zz <= 0:
            continue
        # adaptive enough for low-z to z~2 Pantheon range
        n = max(500, int(2000 * min(zz, 2.5)))
        grid = np.linspace(0, zz, n)
        E = np.sqrt(omega_m * (1.0 + grid)**3 + (1.0 - omega_m))
        dc = (C_LIGHT / H0) * np.trapz(1.0 / np.maximum(E, EPS), grid)
        dl = (1.0 + zz) * dc
        out[i] = 5.0 * np.log10(max(dl, EPS)) + 25.0
    return out


def prepare_data(df: pd.DataFrame, h0=67.4, omega_m=0.3) -> pd.DataFrame:
    """
    Normalize columns and create residuals.
    """
    out = df.copy()

    # Flexible expected columns.
    if "CID" in out.columns:
        out["SNID"] = out["CID"].astype(str)
    elif "SNID" not in out.columns:
        out["SNID"] = np.arange(len(out)).astype(str)

    # Host mass
    if "HOST_LOGMASS" in out.columns:
        out["logMstar"] = pd.to_numeric(out["HOST_LOGMASS"], errors="coerce")
    elif "MHOST" in out.columns:
        out["logMstar"] = pd.to_numeric(out["MHOST"], errors="coerce")
    else:
        raise ValueError("No HOST_LOGMASS/MHOST column found.")

    # Redshift
    zcol = "zHD" if "zHD" in out.columns else ("zCMB" if "zCMB" in out.columns else "z")
    out["z_use"] = pd.to_numeric(out[zcol], errors="coerce")

    # Distance modulus
    if "MU_SH0ES" in out.columns:
        out["mu_obs"] = pd.to_numeric(out["MU_SH0ES"], errors="coerce")
    elif "MU" in out.columns:
        out["mu_obs"] = pd.to_numeric(out["MU"], errors="coerce")
    else:
        raise ValueError("No MU_SH0ES/MU column found.")

    # Remove common sentinel values.
    for col in ["logMstar", "z_use", "mu_obs"]:
        out.loc[out[col] < -90, col] = np.nan

    out["mu_lcdm_67p4"] = lcdm_mu(out["z_use"].to_numpy(), H0=h0, omega_m=omega_m)
    out["mu_resid_67p4"] = out["mu_obs"] - out["mu_lcdm_67p4"]

    # Also create a mean-centered residual to remove absolute intercept/calibration offset.
    valid_mu = out["mu_resid_67p4"].replace([np.inf, -np.inf], np.nan)
    out["mu_resid_centered"] = valid_mu - valid_mu.mean(skipna=True)

    # H0 proxy from centered residual: negative residual => closer/brighter => higher inferred H0.
    out["H0_proxy_centered"] = h0 * 10.0 ** (-out["mu_resid_centered"] / 5.0)

    # Flags where available.
    for flag in ["IS_CALIBRATOR", "USED_IN_SH0ES_HF"]:
        if flag in out.columns:
            out[flag] = pd.to_numeric(out[flag], errors="coerce")

    # Usable rows
    out = out.dropna(subset=["logMstar", "z_use", "mu_obs", "mu_resid_centered"]).copy()

    # Keep plausible host masses.
    out = out[(out["logMstar"] > 5.0) & (out["logMstar"] < 13.0)].copy()

    return out


def design_matrix(df: pd.DataFrame, predictors: List[str]):
    X = [np.ones(len(df))]
    names = ["intercept"]
    for p in predictors:
        vals = pd.to_numeric(df[p], errors="coerce").to_numpy(float)
        mu, sd = np.nanmean(vals), np.nanstd(vals)
        if sd <= 0 or not np.isfinite(sd):
            sd = 1.0
        X.append((vals - mu) / sd)
        names.append(p)
    return np.vstack(X).T, names


def fit_ols(df: pd.DataFrame, model: str, predictors: List[str], ycol="mu_resid_centered"):
    d = df.dropna(subset=[ycol] + predictors).copy()
    y = d[ycol].to_numpy(float)
    X, names = design_matrix(d, predictors)
    mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y = y[mask]
    X = X[mask]

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    n, k = len(y), X.shape[1]
    rss = float(np.sum(resid**2))
    tss = float(np.sum((y - y.mean())**2))
    rmse = float(np.sqrt(rss / max(n, 1)))
    r2 = 1.0 - rss / max(tss, EPS)
    aic = n * np.log(max(rss / max(n, 1), EPS)) + 2 * k
    bic = n * np.log(max(rss / max(n, 1), EPS)) + k * np.log(max(n, 1))

    # Standard errors
    try:
        sigma2 = rss / max(n - k, 1)
        cov = sigma2 * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(cov))
    except Exception:
        se = np.full(k, np.nan)

    out = {
        "model": model,
        "predictors": ",".join(predictors) if predictors else "intercept_only",
        "N": n,
        "k_params": k,
        "RMSE": rmse,
        "R2": r2,
        "AIC": aic,
        "BIC": bic,
    }
    for nm, b, s in zip(names, beta, se):
        out[f"coef_{nm}"] = b
        out[f"se_{nm}"] = s
        out[f"t_{nm}"] = b / s if np.isfinite(s) and s > 0 else np.nan
    return out


def run_models(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.append(fit_ols(df, "M0_intercept_only", []))
    rows.append(fit_ols(df, "M1_linear_logMstar", ["logMstar"]))

    # Mass step model: binary high/low split at 10.
    tmp = df.copy()
    tmp["mass_high_10"] = (tmp["logMstar"] >= 10.0).astype(float)
    rows.append(fit_ols(tmp, "M2_mass_step_10", ["mass_high_10"]))

    # Redshift-control versions, to test whether mass trend is just redshift/sample mix.
    rows.append(fit_ols(df, "M3_logMstar_plus_z", ["logMstar", "z_use"]))
    tmp["logz"] = np.log10(np.maximum(tmp["z_use"], EPS))
    rows.append(fit_ols(tmp, "M4_logMstar_plus_logz", ["logMstar", "logz"]))

    return pd.DataFrame(rows)


def binned_mass_table(df: pd.DataFrame) -> pd.DataFrame:
    bins = [5, 8.5, 9.5, 10.0, 10.5, 11.0, 13.0]
    labels = ["<8.5", "8.5-9.5", "9.5-10", "10-10.5", "10.5-11", "11+"]
    tmp = df.copy()
    tmp["mass_bin"] = pd.cut(tmp["logMstar"], bins=bins, labels=labels, include_lowest=True)
    b = tmp.groupby("mass_bin", observed=False).agg(
        N=("SNID", "size"),
        mean_logMstar=("logMstar", "mean"),
        mean_z=("z_use", "mean"),
        mean_mu_resid=("mu_resid_centered", "mean"),
        std_mu_resid=("mu_resid_centered", "std"),
        median_mu_resid=("mu_resid_centered", "median"),
        mean_H0_proxy=("H0_proxy_centered", "mean"),
    ).reset_index()
    b["sem_mu_resid"] = b["std_mu_resid"] / np.sqrt(np.maximum(b["N"], 1))
    return b


def split_summary(df: pd.DataFrame) -> Dict[str, float]:
    low = df[df["logMstar"] < 10.0]
    high = df[df["logMstar"] >= 10.0]

    out = {
        "N_total": len(df),
        "N_low_mass": len(low),
        "N_high_mass": len(high),
        "mean_mu_low_mass": float(low["mu_resid_centered"].mean()),
        "mean_mu_high_mass": float(high["mu_resid_centered"].mean()),
        "mass_step_mu_high_minus_low": float(high["mu_resid_centered"].mean() - low["mu_resid_centered"].mean()),
        "mean_H0_low_mass": float(low["H0_proxy_centered"].mean()),
        "mean_H0_high_mass": float(high["H0_proxy_centered"].mean()),
        "H0_high_minus_low": float(high["H0_proxy_centered"].mean() - low["H0_proxy_centered"].mean()),
        "pearson_r_mu_logM": float(np.corrcoef(df["mu_resid_centered"], df["logMstar"])[0,1]),
    }

    if "IS_CALIBRATOR" in df.columns:
        cal = df[df["IS_CALIBRATOR"] == 1]
        non = df[df["IS_CALIBRATOR"] != 1]
        out.update({
            "N_calibrator": len(cal),
            "N_noncalibrator": len(non),
            "mean_mu_calibrator": float(cal["mu_resid_centered"].mean()) if len(cal) else np.nan,
            "mean_mu_noncalibrator": float(non["mu_resid_centered"].mean()) if len(non) else np.nan,
        })

    if "USED_IN_SH0ES_HF" in df.columns:
        hf = df[df["USED_IN_SH0ES_HF"] == 1]
        out.update({
            "N_used_in_shoes_hf": len(hf),
            "mean_mu_used_in_shoes_hf": float(hf["mu_resid_centered"].mean()) if len(hf) else np.nan,
        })

    return out


def save_plot(df: pd.DataFrame, results: pd.DataFrame, binned: pd.DataFrame, outdir: Path):
    plt.figure(figsize=(14, 10))

    plt.subplot(2, 2, 1)
    plt.scatter(df["logMstar"], df["mu_resid_centered"], s=9, alpha=0.35)
    m, b = np.polyfit(df["logMstar"], df["mu_resid_centered"], 1)
    xs = np.linspace(df["logMstar"].min(), df["logMstar"].max(), 100)
    plt.plot(xs, m*xs + b, color="black")
    plt.axvline(10.0, linestyle="--", color="black")
    plt.xlabel("HOST_LOGMASS")
    plt.ylabel("centered residual Δμ vs LCDM H0=67.4")
    plt.title("Pantheon+ mass-only baseline")

    plt.subplot(2, 2, 2)
    plt.errorbar(
        np.arange(len(binned)),
        binned["mean_mu_resid"],
        yerr=binned["sem_mu_resid"],
        marker="o",
    )
    plt.xticks(np.arange(len(binned)), binned["mass_bin"], rotation=25)
    plt.axhline(0, linestyle="--", color="black")
    plt.xlabel("host mass bin")
    plt.ylabel("mean centered Δμ")
    plt.title("Binned host-mass step")

    plt.subplot(2, 2, 3)
    plt.bar(results["model"], results["BIC"])
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("BIC")
    plt.title("Mass baseline model comparison")

    plt.subplot(2, 2, 4)
    plt.errorbar(
        np.arange(len(binned)),
        binned["mean_H0_proxy"],
        yerr=None,
        marker="o",
    )
    plt.xticks(np.arange(len(binned)), binned["mass_bin"], rotation=25)
    plt.axhline(67.4, linestyle="--", color="black", label="global/CMB ref")
    plt.xlabel("host mass bin")
    plt.ylabel("mean H0 proxy")
    plt.title("H0 proxy by host mass")
    plt.legend()

    plt.tight_layout()
    plt.savefig(outdir / "flux_pantheon_mass_baseline.png", dpi=180)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Pantheon+SH0ES.dat")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--h0", type=float, default=67.4)
    ap.add_argument("--omega-m", type=float, default=0.3)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = read_pantheon_dat(args.input)
    df = prepare_data(raw, h0=args.h0, omega_m=args.omega_m)
    results = run_models(df)
    bins = binned_mass_table(df)
    split = split_summary(df)

    best_model = str(results.sort_values("BIC").iloc[0]["model"])
    bic_intercept = float(results.loc[results["model"] == "M0_intercept_only", "BIC"].iloc[0])
    bic_mass = float(results.loc[results["model"] == "M1_linear_logMstar", "BIC"].iloc[0])
    bic_step = float(results.loc[results["model"] == "M2_mass_step_10", "BIC"].iloc[0])

    summary = pd.DataFrame([{
        "verdict": "PANTHEON_MASS_BASELINE_COMPLETE",
        "input_rows": len(raw),
        "usable_rows": len(df),
        "best_BIC_model": best_model,
        "BIC_intercept_only": bic_intercept,
        "BIC_linear_mass": bic_mass,
        "BIC_mass_step_10": bic_step,
        "delta_BIC_linear_mass_minus_intercept": bic_mass - bic_intercept,
        "delta_BIC_mass_step_minus_intercept": bic_step - bic_intercept,
        **split,
        "interpretation": (
            "This establishes the host-mass control baseline. The later Flux capacity model must beat this mass-only BIC and survive fixed-mass tests."
        ),
    }])

    df.to_csv(outdir / "flux_pantheon_mass_baseline_data.csv", index=False)
    results.to_csv(outdir / "flux_pantheon_mass_baseline_results.csv", index=False)
    summary.to_csv(outdir / "flux_pantheon_mass_baseline_summary.csv", index=False)
    bins.to_csv(outdir / "flux_pantheon_mass_baseline_binned.csv", index=False)
    save_plot(df, results, bins, outdir)

    print("\nPantheon+ host-mass baseline audit")
    print("=" * 88)
    print(summary.to_string(index=False))
    print("\nModel comparison")
    print(results[["model", "N", "k_params", "RMSE", "R2", "AIC", "BIC"]].sort_values("BIC").to_string(index=False))
    print("\nMass bins")
    print(bins.to_string(index=False))


if __name__ == "__main__":
    main()
