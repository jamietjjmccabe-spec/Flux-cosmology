#!/usr/bin/env python3
"""
flux_pantheon_sh0es_kill_test.py

Level-3 object-level falsification template for the Flux/MIP Capacity-Overflow
host-environment hypothesis.

Purpose
-------
Merge real Pantheon+/SH0ES Type Ia SN residual data with object-level host
properties from Kelsey/TITAN/DES/SDSS-like catalogues, then test whether
the Flux well-filling variable predicts Hubble residuals better than host mass.

Core hypothesis
---------------
The H0 tension is a local calibration offset caused by registration saturation
in mature host wells.

Young / underfilled hosts:
    low F_fill, absorbing regime, H0 proxy near global/CMB value

Near-capacity hosts:
    transition regime, partial H0 offset

Old / compact / high-dispersion hosts:
    overflowing regime, strong H0 offset and residual shift

Key object-level test
---------------------
At fixed host mass, do older / more compact / higher-sigma hosts show residuals
consistent with higher F_fill?

If yes:
    host mass is a shadow variable.
If no:
    the saturation branch fails or must be revised.

Inputs
------
1. Pantheon+/SH0ES SN table:
    Pantheon+SH0ES.dat or CSV converted equivalent.

2. Host-property table:
    object-level host catalogue containing any of:
        stellar mass
        effective radius
        stellar/progenitor age
        velocity dispersion
        sSFR
        metallicity
        morphology

Flexible aliases are supported.

Outputs
-------
flux_killtest_merged.csv
flux_killtest_results.csv
flux_killtest_summary.csv
flux_killtest_regimes.csv
flux_killtest_binned.csv
flux_killtest_correlations.csv
flux_pantheon_sh0es_kill_test.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EPS = 1e-12
AGE_UNIVERSE = 13.8
H0_GLOBAL = 67.4

ALIASES = {
    "sn_id": ["SNID", "CID", "sn_id", "SN", "name", "IAU_NAME", "CIDint"],
    "host_id": ["host_id", "HOSTID", "host_name", "HOST_NAME", "galaxy_id", "host"],
    "z": ["zHD", "zCMB", "z", "redshift", "zcmb"],
    "mu": ["MU_SH0ES", "MU", "mu", "distance_modulus"],
    "mu_err": ["MU_SH0ES_ERR", "MUERR", "mu_err", "distance_modulus_err"],
    "mu_resid": ["mu_resid", "HR", "Hubble_residual", "residual_mu", "mu_residual"],
    "mb_corr": ["m_b_corr", "mB_corr", "mb_corr"],
    "logMstar": ["MHOST", "logMstar", "logM", "host_logmass", "HOST_LOGMASS", "log_mass"],
    "Reff_kpc": ["Reff_kpc", "Reff", "R_eff", "host_reff_kpc", "Rh", "R50_kpc", "half_light_radius_kpc"],
    "stellar_age_Gyr": ["stellar_age_Gyr", "age_Gyr", "host_age", "progenitor_age", "TITAN_age", "age"],
    "sigma_v_km_s": ["sigma_v_km_s", "sigma_v", "velocity_dispersion", "VDISP", "sigma"],
    "sSFR": ["sSFR", "logsSFR", "host_sSFR", "ssfr", "log_sSFR"],
    "metallicity": ["metallicity", "host_metallicity", "Z", "OH12", "12logOH"],
    "morphology": ["morphology", "host_type", "type", "host_morphology", "TType"],
}


def find_col(df: pd.DataFrame, canonical: str):
    for c in ALIASES.get(canonical, [canonical]):
        if c in df.columns:
            return c
    return None


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in [".csv"]:
        return pd.read_csv(path)
    if path.suffix.lower() in [".tsv", ".txt"]:
        return pd.read_csv(path, sep=None, engine="python")
    # Pantheon .dat is whitespace-separated
    return pd.read_csv(path, delim_whitespace=True, comment="#")


def standardize_ids(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().str.replace(r"\s+", "", regex=True)


def normalize_sn_table(sn: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=sn.index)

    for key in ["sn_id", "host_id", "z", "mu", "mu_err", "mu_resid", "mb_corr", "logMstar"]:
        col = find_col(sn, key)
        if col is not None:
            out[key] = sn[col]
        else:
            out[key] = np.nan

    out["sn_id"] = standardize_ids(out["sn_id"])
    out["host_id"] = standardize_ids(out["host_id"]) if out["host_id"].notna().any() else np.nan

    for c in ["z", "mu", "mu_err", "mu_resid", "mb_corr", "logMstar"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    return out


def normalize_host_table(host: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=host.index)

    for key in ["sn_id", "host_id", "logMstar", "Reff_kpc", "stellar_age_Gyr",
                "sigma_v_km_s", "sSFR", "metallicity", "morphology"]:
        col = find_col(host, key)
        if col is not None:
            out[key] = host[col]
        else:
            out[key] = np.nan

    out["sn_id"] = standardize_ids(out["sn_id"]) if out["sn_id"].notna().any() else np.nan
    out["host_id"] = standardize_ids(out["host_id"]) if out["host_id"].notna().any() else np.nan

    for c in ["logMstar", "Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s", "sSFR", "metallicity"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    # If sSFR column is log-like, preserve logsSFR.
    if out["sSFR"].notna().any() and out["sSFR"].median(skipna=True) < 0:
        out["logsSFR"] = out["sSFR"]
        out["sSFR"] = 10 ** out["logsSFR"]
    else:
        out["logsSFR"] = np.log10(np.maximum(out["sSFR"], EPS))

    out["morphology"] = out["morphology"].astype(str)
    return out


def merge_sn_host(sn: pd.DataFrame, host: pd.DataFrame, merge_on: str = "auto") -> pd.DataFrame:
    sn_norm = normalize_sn_table(sn)
    host_norm = normalize_host_table(host)

    if merge_on == "auto":
        # Prefer SN ID because host names are often inconsistent.
        if sn_norm["sn_id"].notna().any() and host_norm["sn_id"].notna().any():
            merge_on = "sn_id"
        elif sn_norm["host_id"].notna().any() and host_norm["host_id"].notna().any():
            merge_on = "host_id"
        else:
            raise ValueError("No usable common key found. Provide SNID/CID or host_id in both tables.")

    merged = pd.merge(sn_norm, host_norm, on=merge_on, how="inner", suffixes=("_sn", "_host"))

    # Coalesce duplicated logMstar if present in both.
    if "logMstar_sn" in merged.columns or "logMstar_host" in merged.columns:
        a = merged["logMstar_sn"] if "logMstar_sn" in merged.columns else np.nan
        b = merged["logMstar_host"] if "logMstar_host" in merged.columns else np.nan
        merged["logMstar"] = pd.Series(b).fillna(a)
    elif "logMstar" not in merged.columns:
        merged["logMstar"] = np.nan

    return merged


def lcdm_mu(z, H0=H0_GLOBAL, omega_m=0.3):
    """
    Low-z safe flat LCDM distance modulus approximation using numerical integral.
    c in km/s. Good enough for residual construction if no residual column exists.
    """
    c = 299792.458
    z = np.asarray(z, dtype=float)
    out = np.full_like(z, np.nan, dtype=float)

    for i, zz in enumerate(z):
        if not np.isfinite(zz) or zz <= 0:
            continue
        grid = np.linspace(0, zz, 800)
        E = np.sqrt(omega_m * (1 + grid) ** 3 + (1 - omega_m))
        dc = (c / H0) * np.trapz(1 / E, grid)
        dl = (1 + zz) * dc  # Mpc
        out[i] = 5 * np.log10(max(dl, EPS)) + 25
    return out


def ensure_residuals(df: pd.DataFrame, residual_mode: str = "auto") -> pd.DataFrame:
    out = df.copy()

    # Existing residual is preferred.
    if residual_mode == "auto" and "mu_resid" in out.columns and out["mu_resid"].notna().any():
        return out

    if residual_mode in ["auto", "lcdm"] and "mu" in out.columns and out["mu"].notna().any():
        out["mu_lcdm_67p4"] = lcdm_mu(out["z"].to_numpy())
        out["mu_resid"] = out["mu"] - out["mu_lcdm_67p4"]
        return out

    raise ValueError(
        "No usable mu_resid column and cannot compute residuals. "
        "Provide mu_resid/HR or MU_SH0ES + z."
    )


def fill_missing_host_values(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    # These fills are only for preliminary audits. The summary will report missingness.
    if "Reff_kpc" not in out.columns:
        out["Reff_kpc"] = np.nan
    if "stellar_age_Gyr" not in out.columns:
        out["stellar_age_Gyr"] = np.nan
    if "sigma_v_km_s" not in out.columns:
        out["sigma_v_km_s"] = np.nan
    if "logsSFR" not in out.columns:
        out["logsSFR"] = np.nan
    if "metallicity" not in out.columns:
        out["metallicity"] = np.nan

    # Conservative imputation flags.
    for c in ["Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s", "logsSFR", "metallicity"]:
        out[f"{c}_was_missing"] = out[c].isna()

    out["Reff_kpc"] = out["Reff_kpc"].fillna(out["Reff_kpc"].median(skipna=True) if out["Reff_kpc"].notna().any() else 5.0)
    out["stellar_age_Gyr"] = out["stellar_age_Gyr"].fillna(out["stellar_age_Gyr"].median(skipna=True) if out["stellar_age_Gyr"].notna().any() else 4.0)
    out["sigma_v_km_s"] = out["sigma_v_km_s"].fillna(out["sigma_v_km_s"].median(skipna=True) if out["sigma_v_km_s"].notna().any() else 150.0)
    out["logsSFR"] = out["logsSFR"].fillna(out["logsSFR"].median(skipna=True) if out["logsSFR"].notna().any() else -10.0)
    out["metallicity"] = out["metallicity"].fillna(out["metallicity"].median(skipna=True) if out["metallicity"].notna().any() else 8.7)

    return out


def normalize_to_median(x):
    x = np.asarray(x, dtype=float)
    return x / max(np.nanmedian(x), EPS)


def add_flux_capacity_variables(df: pd.DataFrame, Dsat=1.15, q_sat=7.0, proxy="hybrid") -> pd.DataFrame:
    out = df.copy()

    M = 10 ** pd.to_numeric(out["logMstar"], errors="coerce").to_numpy()
    Re = np.maximum(pd.to_numeric(out["Reff_kpc"], errors="coerce").to_numpy(), EPS)
    age = np.maximum(pd.to_numeric(out["stellar_age_Gyr"], errors="coerce").to_numpy(), EPS)
    sig = np.maximum(pd.to_numeric(out["sigma_v_km_s"], errors="coerce").to_numpy(), EPS)

    D_MR_age = normalize_to_median((M / Re) * (age / AGE_UNIVERSE))
    D_sigma_age = normalize_to_median((sig / 150.0) ** 2 * (age / AGE_UNIVERSE))
    D_hybrid = normalize_to_median(np.sqrt((M / Re) * (sig / 150.0) ** 2) * (age / AGE_UNIVERSE))

    if proxy == "MR_age":
        D = D_MR_age
    elif proxy == "sigma_age":
        D = D_sigma_age
    elif proxy == "hybrid":
        D = D_hybrid
    else:
        raise ValueError(proxy)

    F_fill = D / max(Dsat, EPS)
    G_sat = F_fill ** q_sat / (1.0 + F_fill ** q_sat)
    excess = np.maximum(F_fill - 1.0, 0.0)
    G_overflow = (excess / (1.0 + excess)) * G_sat
    G_absorb = 1.0 - G_sat

    out["D_MR_age"] = D_MR_age
    out["D_sigma_age"] = D_sigma_age
    out["D_hybrid"] = D_hybrid
    out["D_reg"] = D
    out["logD_reg"] = np.log10(np.maximum(D, EPS))
    out["F_fill"] = F_fill
    out["G_absorb"] = G_absorb
    out["G_sat"] = G_sat
    out["G_overflow"] = G_overflow

    out["well_regime"] = np.select(
        [F_fill < 0.75, (F_fill >= 0.75) & (F_fill < 1.25), F_fill >= 1.25],
        ["filling/absorbing", "near-capacity/saturated", "overflowing"],
        default="unknown",
    )

    # H0 proxy from residual: negative residual -> closer/brighter -> higher H0 proxy.
    out["H0_proxy"] = H0_GLOBAL * 10 ** (-out["mu_resid"].to_numpy(dtype=float) / 5.0)
    return out


def design_matrix(df, predictors: List[str]):
    X = [np.ones(len(df))]
    names = ["intercept"]
    for p in predictors:
        vals = pd.to_numeric(df[p], errors="coerce").to_numpy(dtype=float)
        mu = np.nanmean(vals)
        sd = np.nanstd(vals)
        if not np.isfinite(sd) or sd <= 0:
            sd = 1.0
        X.append((vals - mu) / sd)
        names.append(p)
    return np.vstack(X).T, names


def fit_ols(df, model, predictors):
    d = df.dropna(subset=["mu_resid"] + predictors).copy()
    y = d["mu_resid"].to_numpy(float)
    X, names = design_matrix(d, predictors)

    # Drop rows with any nonfinite X
    mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y = y[mask]
    X = X[mask]

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    n, k = len(y), X.shape[1]
    rss = float(np.sum(resid ** 2))
    tss = float(np.sum((y - y.mean()) ** 2))
    aic = n * np.log(max(rss / n, EPS)) + 2 * k
    bic = n * np.log(max(rss / n, EPS)) + k * np.log(max(n, 1))

    out = {
        "model": model,
        "predictors": ",".join(predictors),
        "N": n,
        "k_params": k,
        "RMSE": float(np.sqrt(rss / n)),
        "R2": 1.0 - rss / max(tss, EPS),
        "AIC": aic,
        "BIC": bic,
    }
    for nm, b in zip(names, beta):
        out[f"coef_{nm}"] = b
    return out


def regression_battle(df: pd.DataFrame) -> pd.DataFrame:
    models = {
        "A_mass_only": ["logMstar"],
        "B_mass_sSFR": ["logMstar", "logsSFR"],
        "C_kitchen_sink": ["logMstar", "logsSFR", "metallicity", "stellar_age_Gyr", "Reff_kpc"],
        "D_logD_reg": ["logD_reg"],
        "E_Gsat": ["G_sat"],
        "F_capacity_state": ["G_absorb", "G_sat", "G_overflow"],
        "G_capacity_controls": ["G_sat", "G_overflow", "logsSFR", "metallicity"],
    }
    return pd.DataFrame([fit_ols(df, m, p) for m, p in models.items()])


def regime_summary(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby("well_regime", observed=False).agg(
        N=("mu_resid", "size"),
        mean_logMstar=("logMstar", "mean"),
        mean_age=("stellar_age_Gyr", "mean"),
        mean_Reff=("Reff_kpc", "mean"),
        mean_sigma=("sigma_v_km_s", "mean"),
        mean_Ffill=("F_fill", "mean"),
        mean_Gsat=("G_sat", "mean"),
        mean_Goverflow=("G_overflow", "mean"),
        mean_mu_resid=("mu_resid", "mean"),
        std_mu_resid=("mu_resid", "std"),
        mean_H0_proxy=("H0_proxy", "mean"),
    ).reset_index()


def binned_table(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.copy()
    tmp["fill_bin"] = pd.qcut(tmp["F_fill"], q=min(7, max(2, len(tmp)//20)), duplicates="drop")
    return tmp.groupby("fill_bin", observed=False).agg(
        N=("mu_resid", "size"),
        mean_Ffill=("F_fill", "mean"),
        mean_Gsat=("G_sat", "mean"),
        mean_Goverflow=("G_overflow", "mean"),
        mean_mu_resid=("mu_resid", "mean"),
        std_mu_resid=("mu_resid", "std"),
        mean_H0_proxy=("H0_proxy", "mean"),
    ).reset_index()


def correlation_table(df: pd.DataFrame) -> pd.DataFrame:
    variables = ["logMstar", "logD_reg", "F_fill", "G_sat", "G_overflow",
                 "stellar_age_Gyr", "Reff_kpc", "sigma_v_km_s", "logsSFR", "metallicity"]
    rows = []
    for v in variables:
        sub = df[["mu_resid", v]].dropna()
        if len(sub) > 3:
            r = float(np.corrcoef(sub["mu_resid"], sub[v])[0, 1])
        else:
            r = np.nan
        rows.append({"variable": v, "pearson_r_with_mu_resid": r, "N": len(sub)})
    return pd.DataFrame(rows)


def fixed_mass_test(df: pd.DataFrame) -> pd.DataFrame:
    """
    Critical Flux test:
    within narrow mass bins, do high-F_fill hosts have different residuals?
    """
    tmp = df.copy()
    tmp["mass_bin"] = pd.cut(tmp["logMstar"], bins=[8, 9.5, 10.0, 10.5, 11.0, 12.5], include_lowest=True)
    rows = []
    for mbin, sub in tmp.groupby("mass_bin", observed=False):
        if len(sub) < 20:
            continue
        med = sub["F_fill"].median()
        lo = sub[sub["F_fill"] <= med]
        hi = sub[sub["F_fill"] > med]
        rows.append({
            "mass_bin": str(mbin),
            "N": len(sub),
            "mean_logMstar": sub["logMstar"].mean(),
            "lowF_mean_mu": lo["mu_resid"].mean(),
            "highF_mean_mu": hi["mu_resid"].mean(),
            "high_minus_low_mu": hi["mu_resid"].mean() - lo["mu_resid"].mean(),
            "lowF_H0": lo["H0_proxy"].mean(),
            "highF_H0": hi["H0_proxy"].mean(),
            "high_minus_low_H0": hi["H0_proxy"].mean() - lo["H0_proxy"].mean(),
        })
    return pd.DataFrame(rows)


def save_plot(df, results, regimes, binned, fixed_mass, outdir: Path):
    plt.figure(figsize=(15, 11))

    plt.subplot(2, 2, 1)
    plt.scatter(df["F_fill"], df["mu_resid"], s=12, alpha=0.35)
    plt.axvline(1.0, linestyle="--", color="black")
    plt.xscale("log")
    plt.xlabel("F_fill = D_reg / D_sat")
    plt.ylabel("Hubble residual Δμ")
    plt.title("Object-level Flux filling test")

    plt.subplot(2, 2, 2)
    order = results.sort_values("BIC")
    plt.bar(order["model"], order["BIC"])
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("BIC")
    plt.title("Regression battle")

    plt.subplot(2, 2, 3)
    if len(regimes):
        plt.bar(regimes["well_regime"], regimes["mean_H0_proxy"])
        plt.axhline(H0_GLOBAL, linestyle="--", color="black")
        plt.xticks(rotation=20, ha="right")
    plt.ylabel("mean H0 proxy")
    plt.title("H0 proxy by well regime")

    plt.subplot(2, 2, 4)
    if len(fixed_mass):
        plt.bar(fixed_mass["mass_bin"], fixed_mass["high_minus_low_mu"])
        plt.axhline(0, linestyle="--", color="black")
        plt.xticks(rotation=30, ha="right")
    plt.ylabel("high-F minus low-F Δμ")
    plt.title("Fixed-mass kill test")

    plt.tight_layout()
    plt.savefig(outdir / "flux_pantheon_sh0es_kill_test.png", dpi=180)
    plt.close()


def run_kill_test(sn_path, host_path, outdir, merge_on="auto", proxy="hybrid", Dsat=1.15, q_sat=7.0, residual_mode="auto"):
    sn = read_table(sn_path)
    host = read_table(host_path)
    merged = merge_sn_host(sn, host, merge_on=merge_on)
    merged = ensure_residuals(merged, residual_mode=residual_mode)
    merged = fill_missing_host_values(merged)

    # Must have mass and residual. Host radius/age/sigma can be imputed, but flagged.
    merged = merged.dropna(subset=["mu_resid", "logMstar", "z"]).copy()
    merged = add_flux_capacity_variables(merged, Dsat=Dsat, q_sat=q_sat, proxy=proxy)

    results = regression_battle(merged)
    regimes = regime_summary(merged)
    bins = binned_table(merged)
    corr = correlation_table(merged)
    fixed = fixed_mass_test(merged)

    best_model = str(results.sort_values("BIC").iloc[0]["model"])
    mass_bic = float(results.loc[results["model"] == "A_mass_only", "BIC"].iloc[0])
    cap_bic = float(results.loc[results["model"] == "F_capacity_state", "BIC"].iloc[0])
    gsat_bic = float(results.loc[results["model"] == "E_Gsat", "BIC"].iloc[0])

    missing_rates = {
        c: float(merged[f"{c}_was_missing"].mean()) if f"{c}_was_missing" in merged.columns else np.nan
        for c in ["Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s", "logsSFR", "metallicity"]
    }

    verdict = (
        "OBJECT_LEVEL_FLUX_SUPPORT"
        if (gsat_bic < mass_bic and cap_bic < mass_bic)
        else "OBJECT_LEVEL_NO_SUPPORT_OR_INCOMPLETE"
    )

    summary = pd.DataFrame([{
        "verdict": verdict,
        "N_merged": len(merged),
        "proxy": proxy,
        "Dsat": Dsat,
        "q_sat": q_sat,
        "best_BIC_model": best_model,
        "BIC_mass_only": mass_bic,
        "BIC_Gsat": gsat_bic,
        "BIC_capacity_state": cap_bic,
        "delta_BIC_Gsat_minus_mass": gsat_bic - mass_bic,
        "delta_BIC_capacity_minus_mass": cap_bic - mass_bic,
        "Reff_missing_rate": missing_rates["Reff_kpc"],
        "age_missing_rate": missing_rates["stellar_age_Gyr"],
        "sigma_missing_rate": missing_rates["sigma_v_km_s"],
        "sSFR_missing_rate": missing_rates["logsSFR"],
        "metallicity_missing_rate": missing_rates["metallicity"],
        "interpretation": (
            "Capacity variables outperform mass-only in the object-level table. Check fixed-mass bins and missingness before making any physical claim."
            if verdict == "OBJECT_LEVEL_FLUX_SUPPORT"
            else
            "Object-level table does not prefer Flux capacity variables over mass-only, or key host variables are incomplete."
        ),
    }])

    merged.to_csv(outdir / "flux_killtest_merged.csv", index=False)
    results.to_csv(outdir / "flux_killtest_results.csv", index=False)
    summary.to_csv(outdir / "flux_killtest_summary.csv", index=False)
    regimes.to_csv(outdir / "flux_killtest_regimes.csv", index=False)
    bins.to_csv(outdir / "flux_killtest_binned.csv", index=False)
    corr.to_csv(outdir / "flux_killtest_correlations.csv", index=False)
    fixed.to_csv(outdir / "flux_killtest_fixed_mass.csv", index=False)
    save_plot(merged, results, regimes, bins, fixed, outdir)

    return summary, results, regimes, fixed


def make_schema_template(path: Path):
    df = pd.DataFrame({
        "SNID": ["SN_EXAMPLE_001"],
        "zHD": [0.023],
        "MU_SH0ES": [35.12],
        "mu_resid": [np.nan],
        "MHOST": [10.55],
        "Reff_kpc": [4.8],
        "stellar_age_Gyr": [5.6],
        "sigma_v_km_s": [180.0],
        "logsSFR": [-10.8],
        "metallicity": [8.85],
        "morphology": ["early/quiescent"],
    })
    df.to_csv(path, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sn", default=None, help="Pantheon+/SH0ES SN table")
    parser.add_argument("--host", default=None, help="Host-property table")
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--merge-on", default="auto", choices=["auto", "sn_id", "host_id"])
    parser.add_argument("--proxy", default="hybrid", choices=["hybrid", "MR_age", "sigma_age"])
    parser.add_argument("--Dsat", type=float, default=1.15)
    parser.add_argument("--q-sat", type=float, default=7.0)
    parser.add_argument("--residual-mode", default="auto", choices=["auto", "lcdm"])
    parser.add_argument("--write-schema-only", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    schema_path = outdir / "flux_pantheon_host_schema_template.csv"
    make_schema_template(schema_path)

    if args.write_schema_only:
        print(f"Wrote schema template: {schema_path}")
        return

    if not args.sn or not args.host:
        raise ValueError("Provide --sn and --host, or use --write-schema-only.")

    summary, results, regimes, fixed = run_kill_test(
        args.sn, args.host, outdir,
        merge_on=args.merge_on,
        proxy=args.proxy,
        Dsat=args.Dsat,
        q_sat=args.q_sat,
        residual_mode=args.residual_mode,
    )

    print("\nFlux Pantheon+/SH0ES object-level kill test")
    print("=" * 92)
    print(summary.to_string(index=False))
    print("\nRegression battle")
    print(results[["model", "N", "k_params", "RMSE", "R2", "AIC", "BIC"]].sort_values("BIC").to_string(index=False))
    print("\nRegime summary")
    print(regimes.to_string(index=False))
    print("\nFixed-mass test")
    print(fixed.to_string(index=False))


if __name__ == "__main__":
    main()
