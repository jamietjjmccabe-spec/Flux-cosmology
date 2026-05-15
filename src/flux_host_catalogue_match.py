#!/usr/bin/env python3
"""
flux_host_catalogue_match.py

Flux/MIP host-environment forensic audit for Type Ia Supernova Hubble residuals.

Purpose
-------
This module prepares the observational audit phase for the Flux saturation branch.

Standard SN Ia cosmology often models residual host dependence as a "host mass step".
Flux/MIP predicts the deeper variable is registration depth:

    D_reg,host ~ (M_star / R_eff)^alpha * (t_stellar / t0)^beta

Mass alone is only a blunt proxy. If Flux saturation is real, Hubble residuals
should correlate better with registration-depth proxies than with logMstar alone.

Modes
-----
1. synthetic:
    Generates a Pantheon+/SH0ES-like mock host catalogue with:
        SN name, z, mu residual, host mass, effective radius, age, sSFR,
        metallicity, morphology, and synthetic true Flux residual.

2. catalogue:
    Loads a user-supplied CSV with comparable columns and performs the same audit.

Required/accepted input columns for catalogue mode
-------------------------------------------------
Flexible aliases are accepted where possible:

    sn_name / CID / SNID / name
    z / zHD / redshift
    mu_resid / residual_mu / Hubble_residual / HR
    logMstar / logM / host_logmass
    Reff_kpc / Reff / host_reff_kpc
    stellar_age_Gyr / age_Gyr / host_age
    sSFR / logsSFR / host_sSFR
    metallicity / host_metallicity
    morphology / host_type

Outputs
-------
flux_host_catalogue_match_data.csv
flux_host_catalogue_match_results.csv
flux_host_catalogue_match_summary.csv
flux_host_catalogue_match_binned.csv
flux_host_catalogue_match.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EPS = 1e-12
AGE0 = 5.0
H0_GLOBAL = 67.4
H0_LOCAL = 73.0


ALIASES = {
    "sn_name": ["sn_name", "SN", "CID", "SNID", "name"],
    "z": ["z", "zHD", "redshift", "z_helio", "zcmb"],
    "mu_resid": ["mu_resid", "residual_mu", "Hubble_residual", "HR", "mu_residual", "resid"],
    "logMstar": ["logMstar", "logM", "host_logmass", "log_mass", "HOST_LOGMASS"],
    "Reff_kpc": ["Reff_kpc", "Reff", "host_reff_kpc", "R_eff", "radius_kpc"],
    "stellar_age_Gyr": ["stellar_age_Gyr", "age_Gyr", "host_age", "stellar_age", "age"],
    "sSFR": ["sSFR", "logsSFR", "host_sSFR", "ssfr"],
    "metallicity": ["metallicity", "host_metallicity", "Z", "OH12"],
    "morphology": ["morphology", "host_type", "type", "host_morphology"],
}


def rng(seed=42):
    return np.random.default_rng(seed)


def alias_get(row: pd.Series, key: str, default=np.nan):
    for col in ALIASES.get(key, [key]):
        if col in row.index and pd.notna(row[col]):
            return row[col]
    return default


def make_synthetic_catalogue(n=2500, seed=42) -> pd.DataFrame:
    """
    Generate a mock Pantheon+/SH0ES-like host catalogue.

    The synthetic catalogue deliberately includes:
    - ordinary host mass effect
    - Flux registration-depth effect
    - scatter/noise
    - morphology/environment fields
    """
    R = rng(seed)

    morph = R.choice(["dwarf", "spiral", "elliptical"], size=n, p=[0.22, 0.55, 0.23])

    logM = np.zeros(n)
    Reff = np.zeros(n)
    age = np.zeros(n)
    logsSFR = np.zeros(n)
    metallicity = np.zeros(n)

    for i, m in enumerate(morph):
        if m == "dwarf":
            logM[i] = R.normal(8.8, 0.45)
            Reff[i] = 10 ** R.normal(0.30, 0.18)
            age[i] = np.clip(R.normal(3.2, 1.4), 0.5, 9.5)
            logsSFR[i] = R.normal(-9.4, 0.45)
            metallicity[i] = R.normal(8.35, 0.18)
        elif m == "spiral":
            logM[i] = R.normal(10.25, 0.40)
            Reff[i] = 10 ** R.normal(0.70, 0.20)
            age[i] = np.clip(R.normal(6.0, 2.0), 1.0, 11.5)
            logsSFR[i] = R.normal(-10.0, 0.45)
            metallicity[i] = R.normal(8.75, 0.15)
        else:
            logM[i] = R.normal(11.0, 0.35)
            Reff[i] = 10 ** R.normal(0.55, 0.20)
            age[i] = np.clip(R.normal(9.0, 1.5), 3.0, 13.0)
            logsSFR[i] = R.normal(-11.45, 0.35)
            metallicity[i] = R.normal(8.95, 0.12)

    logM = np.clip(logM, 7.2, 12.2)
    Reff = np.clip(Reff, 0.35, 25.0)

    # Redshift mix: many low-z ladder/calibrator-like objects plus Hubble-flow sample.
    z = np.concatenate([
        R.uniform(0.005, 0.05, int(0.35 * n)),
        R.uniform(0.05, 0.35, n - int(0.35 * n)),
    ])
    R.shuffle(z)

    df = pd.DataFrame({
        "sn_name": [f"SN_SYN_{i:05d}" for i in range(n)],
        "z": z,
        "logMstar": logM,
        "Reff_kpc": Reff,
        "stellar_age_Gyr": age,
        "logsSFR": logsSFR,
        "sSFR": 10 ** logsSFR,
        "metallicity": metallicity,
        "morphology": morph,
    })

    # Synthetic Flux residual construction.
    df = add_registration_depth(df)

    # Use a saturation-like residual but keep amplitude closer to observed mass-step scale
    # than the full 8% H0 synthetic bridge. This tests detectability in residuals.
    G = df["G_sat_proxy"].to_numpy()
    mass_step_component = -0.035 * (df["logMstar"].to_numpy() >= 10.0).astype(float)
    flux_component = -0.075 * G
    sfr_component = 0.012 * (df["logsSFR"].to_numpy() + 10.0)
    metallicity_component = -0.010 * (df["metallicity"].to_numpy() - 8.7)
    noise = R.normal(0.0, 0.10, n)

    df["mu_resid"] = mass_step_component + flux_component + sfr_component + metallicity_component + noise
    df["mu_resid_true_flux_component"] = flux_component
    df["mu_resid_true_mass_step_component"] = mass_step_component

    return df


def normalize_catalogue(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        out = {}
        for key in ["sn_name", "z", "mu_resid", "logMstar", "Reff_kpc", "stellar_age_Gyr", "sSFR", "metallicity", "morphology"]:
            out[key] = alias_get(row, key, default=np.nan)
        rows.append(out)
    out = pd.DataFrame(rows)

    out["sn_name"] = out["sn_name"].astype(str)
    for col in ["z", "mu_resid", "logMstar", "Reff_kpc", "stellar_age_Gyr", "sSFR", "metallicity"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    # If sSFR was supplied as log, detect and convert into linear where plausible.
    # Negative values in an sSFR column usually mean log10(sSFR).
    if out["sSFR"].notna().any() and out["sSFR"].median(skipna=True) < 0:
        out["logsSFR"] = out["sSFR"]
        out["sSFR"] = 10 ** out["logsSFR"]
    else:
        out["logsSFR"] = np.log10(np.maximum(out["sSFR"], EPS))

    # Fill missing optional host values conservatively.
    out["Reff_kpc"] = out["Reff_kpc"].fillna(out["Reff_kpc"].median(skipna=True) if out["Reff_kpc"].notna().any() else 4.0)
    out["stellar_age_Gyr"] = out["stellar_age_Gyr"].fillna(out["stellar_age_Gyr"].median(skipna=True) if out["stellar_age_Gyr"].notna().any() else 5.0)
    out["metallicity"] = out["metallicity"].fillna(out["metallicity"].median(skipna=True) if out["metallicity"].notna().any() else 8.7)
    out["morphology"] = out["morphology"].fillna("unknown").astype(str)

    return out


def add_registration_depth(df: pd.DataFrame, alpha=1.0, beta=0.8, Dsat=1.25, q=4.0) -> pd.DataFrame:
    out = df.copy()
    Mstar = 10 ** out["logMstar"].to_numpy()
    Reff = np.maximum(out["Reff_kpc"].to_numpy(), EPS)
    age = np.maximum(out["stellar_age_Gyr"].to_numpy(), EPS)

    raw = ((Mstar / Reff) ** alpha) * ((age / AGE0) ** beta)
    norm = np.nanmedian(raw)
    Dreg = raw / max(norm, EPS)

    Gsat = (Dreg ** q) / (Dreg ** q + Dsat ** q)

    out["D_reg_proxy"] = Dreg
    out["logD_reg_proxy"] = np.log10(np.maximum(Dreg, EPS))
    out["G_sat_proxy"] = Gsat
    out["potential_proxy_logM_over_Reff"] = out["logMstar"] - np.log10(np.maximum(out["Reff_kpc"], EPS))
    out["age_factor"] = out["stellar_age_Gyr"] / AGE0
    return out


def design_matrix(df: pd.DataFrame, model: str) -> Tuple[np.ndarray, List[str]]:
    """
    Build regression design matrix.
    Includes intercept automatically.
    """
    if model == "A_mass":
        cols = ["logMstar"]
    elif model == "B_mass_sSFR":
        cols = ["logMstar", "logsSFR"]
    elif model == "C_multivariate":
        cols = ["logMstar", "logsSFR", "Reff_kpc", "stellar_age_Gyr", "metallicity"]
    elif model == "D_flux_depth":
        cols = ["logD_reg_proxy"]
    elif model == "E_flux_depth_plus_controls":
        cols = ["logD_reg_proxy", "logsSFR", "metallicity"]
    else:
        raise ValueError(f"Unknown model {model}")

    Xcols = []
    names = ["intercept"]
    Xcols.append(np.ones(len(df)))

    for col in cols:
        vals = df[col].to_numpy(dtype=float)
        if col == "Reff_kpc":
            vals = np.log10(np.maximum(vals, EPS))
            name = "logReff_kpc"
        else:
            name = col

        # standardize non-intercept columns for comparable coefficients
        mu = np.nanmean(vals)
        sig = np.nanstd(vals)
        sig = sig if sig > 0 else 1.0
        Xcols.append((vals - mu) / sig)
        names.append(name)

    X = np.vstack(Xcols).T
    return X, names


def fit_ols(df: pd.DataFrame, model: str) -> Dict[str, float]:
    d = df.dropna(subset=["mu_resid", "logMstar", "logD_reg_proxy", "logsSFR", "Reff_kpc", "stellar_age_Gyr", "metallicity"]).copy()
    y = d["mu_resid"].to_numpy(dtype=float)

    X, names = design_matrix(d, model)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    n = len(y)
    k = X.shape[1]

    rss = float(np.sum(resid ** 2))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - rss / max(tss, EPS)
    sigma2 = rss / max(n, 1)
    aic = n * np.log(max(rss / n, EPS)) + 2 * k
    bic = n * np.log(max(rss / n, EPS)) + k * np.log(max(n, 1))

    # standard errors
    try:
        cov = sigma2 * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(cov))
    except Exception:
        se = np.full(k, np.nan)

    out = {
        "model": model,
        "N": n,
        "k_params": k,
        "RSS": rss,
        "RMSE": float(np.sqrt(rss / n)),
        "R2": r2,
        "AIC": aic,
        "BIC": bic,
    }

    for name, b, s in zip(names, beta, se):
        out[f"coef_{name}"] = b
        out[f"se_{name}"] = s
        out[f"t_{name}"] = b / s if np.isfinite(s) and s > 0 else np.nan

    return out


def correlation_table(df: pd.DataFrame) -> pd.DataFrame:
    variables = [
        "logMstar",
        "potential_proxy_logM_over_Reff",
        "stellar_age_Gyr",
        "logD_reg_proxy",
        "G_sat_proxy",
        "logsSFR",
        "metallicity",
    ]
    rows = []
    for v in variables:
        sub = df[["mu_resid", v]].dropna()
        if len(sub) < 3:
            r = np.nan
        else:
            r = float(np.corrcoef(sub["mu_resid"], sub[v])[0, 1])
        rows.append({"variable": v, "pearson_r_with_mu_resid": r, "N": len(sub)})
    return pd.DataFrame(rows)


def binned_depth_table(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.copy()
    tmp["D_bin"] = pd.qcut(tmp["logD_reg_proxy"], q=6, duplicates="drop")
    return tmp.groupby("D_bin", observed=False).agg(
        N=("mu_resid", "size"),
        mean_logD=("logD_reg_proxy", "mean"),
        mean_logMstar=("logMstar", "mean"),
        mean_mu_resid=("mu_resid", "mean"),
        std_mu_resid=("mu_resid", "std"),
        mean_Gsat=("G_sat_proxy", "mean"),
    ).reset_index()


def inferred_h0_from_mu_residual(mu_resid: np.ndarray, H0_ref=H0_GLOBAL):
    """
    Approximate relation:
        D_L inferred ∝ 10^(mu/5)
        H0 inferred ∝ 1/D_L
    If residual is negative/brighter/closer, H0 inferred increases.
    """
    return H0_ref * 10 ** (-np.asarray(mu_resid) / 5.0)


def add_h0_trend(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["H0_inferred_proxy"] = inferred_h0_from_mu_residual(out["mu_resid"].to_numpy())
    return out


def save_plot(df: pd.DataFrame, results: pd.DataFrame, binned: pd.DataFrame, outdir: Path):
    plt.figure(figsize=(14, 10))

    plt.subplot(2, 2, 1)
    plt.scatter(df["logD_reg_proxy"], df["mu_resid"], s=9, alpha=0.35)
    m, b = np.polyfit(df["logD_reg_proxy"], df["mu_resid"], 1)
    xs = np.linspace(df["logD_reg_proxy"].min(), df["logD_reg_proxy"].max(), 100)
    plt.plot(xs, m * xs + b, color="black")
    plt.xlabel("log registration-depth proxy")
    plt.ylabel("SN Ia Hubble residual Δμ")
    plt.title("Flux fingerprint: residual vs D_reg")

    plt.subplot(2, 2, 2)
    plt.scatter(df["logMstar"], df["mu_resid"], s=9, alpha=0.35)
    m, b = np.polyfit(df["logMstar"], df["mu_resid"], 1)
    xs = np.linspace(df["logMstar"].min(), df["logMstar"].max(), 100)
    plt.plot(xs, m * xs + b, color="black")
    plt.xlabel("log host stellar mass")
    plt.ylabel("SN Ia Hubble residual Δμ")
    plt.title("Mainstream proxy: residual vs mass")

    plt.subplot(2, 2, 3)
    order = results.sort_values("BIC")
    plt.bar(order["model"], order["BIC"])
    plt.xticks(rotation=25, ha="right")
    plt.ylabel("BIC")
    plt.title("Regression battle: lower BIC wins")

    plt.subplot(2, 2, 4)
    plt.errorbar(
        np.arange(len(binned)),
        binned["mean_mu_resid"],
        yerr=binned["std_mu_resid"] / np.sqrt(np.maximum(binned["N"], 1)),
        marker="o",
    )
    plt.xticks(np.arange(len(binned)), [str(x) for x in binned["D_bin"]], rotation=30, ha="right")
    plt.axhline(0.0, linestyle="--", color="black")
    plt.xlabel("D_reg quantile bin")
    plt.ylabel("mean Δμ")
    plt.title("Saturation sequence")

    plt.tight_layout()
    plt.savefig(outdir / "flux_host_catalogue_match.png", dpi=180)
    plt.close()


def run_audit(df: pd.DataFrame, outdir: Path):
    df = normalize_catalogue(df)
    df = df.dropna(subset=["mu_resid", "logMstar"]).copy()
    df = add_registration_depth(df)
    df = add_h0_trend(df)

    models = ["A_mass", "B_mass_sSFR", "C_multivariate", "D_flux_depth", "E_flux_depth_plus_controls"]
    results = pd.DataFrame([fit_ols(df, m) for m in models])
    corr = correlation_table(df)
    binned = binned_depth_table(df)

    best_bic = results.sort_values("BIC").iloc[0]["model"]
    bic_mass = float(results.loc[results["model"] == "A_mass", "BIC"].iloc[0])
    bic_flux = float(results.loc[results["model"] == "D_flux_depth", "BIC"].iloc[0])
    delta_bic_flux_vs_mass = bic_flux - bic_mass

    highD = df[df["logD_reg_proxy"] >= df["logD_reg_proxy"].quantile(0.75)]
    lowD = df[df["logD_reg_proxy"] <= df["logD_reg_proxy"].quantile(0.25)]
    h0_highD = highD["H0_inferred_proxy"].mean()
    h0_lowD = lowD["H0_inferred_proxy"].mean()

    verdict = "FLUX_PROXY_WINS_SYNTHETIC" if best_bic in ["D_flux_depth", "E_flux_depth_plus_controls"] else "MASS_OR_CONTROLS_WIN"

    summary = pd.DataFrame([{
        "verdict": verdict,
        "best_BIC_model": best_bic,
        "BIC_mass_model_A": bic_mass,
        "BIC_flux_model_D": bic_flux,
        "delta_BIC_flux_minus_mass": delta_bic_flux_vs_mass,
        "H0_proxy_lowD_quartile": h0_lowD,
        "H0_proxy_highD_quartile": h0_highD,
        "H0_proxy_high_minus_low": h0_highD - h0_lowD,
        "mean_mu_lowD_quartile": lowD["mu_resid"].mean(),
        "mean_mu_highD_quartile": highD["mu_resid"].mean(),
        "interpretation": (
            "Registration-depth proxy wins or remains competitive against host mass, supporting the saturation-audit strategy."
            if verdict == "FLUX_PROXY_WINS_SYNTHETIC"
            else
            "Mass or standard controls beat the registration-depth proxy in this audit; saturation branch would need revision or real-data retest."
        ),
    }])

    df.to_csv(outdir / "flux_host_catalogue_match_data.csv", index=False)
    results.to_csv(outdir / "flux_host_catalogue_match_results.csv", index=False)
    summary.to_csv(outdir / "flux_host_catalogue_match_summary.csv", index=False)
    binned.to_csv(outdir / "flux_host_catalogue_match_binned.csv", index=False)
    corr.to_csv(outdir / "flux_host_catalogue_match_correlations.csv", index=False)

    save_plot(df, results, binned, outdir)

    return summary, results, binned, corr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["synthetic", "catalogue"], default="synthetic")
    parser.add_argument("--input", type=str, default=None)
    parser.add_argument("--outdir", type=str, default=".")
    parser.add_argument("--n", type=int, default=2500)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.mode == "synthetic":
        df = make_synthetic_catalogue(n=args.n, seed=args.seed)
    else:
        if not args.input:
            raise ValueError("--input is required in catalogue mode")
        df = pd.read_csv(args.input)

    summary, results, binned, corr = run_audit(df, outdir)

    print("\nFlux host-catalogue forensic audit")
    print("=" * 88)
    print(summary.to_string(index=False))
    print("\nRegression battle")
    cols = ["model", "N", "k_params", "RMSE", "R2", "AIC", "BIC"]
    print(results[cols].sort_values("BIC").to_string(index=False))
    print("\nCorrelation table")
    print(corr.to_string(index=False))
    print("\nBinned depth trend")
    print(binned.to_string(index=False))


if __name__ == "__main__":
    main()
