#!/usr/bin/env python3
"""
flux_neighborhood_audit.py

Level-3B Flux/MIP Neighborhood Forensic Audit.

Purpose
-------
This script starts from the already-built Pantheon+ mass-baseline file:

    flux_pantheon_mass_baseline_data.csv

and merges it with an object-level host-structure/maturity proxy table containing
any combination of:

    Reff_kpc
    stellar_age_Gyr
    sigma_v_km_s
    logsSFR / sSFR
    metallicity
    morphology

It then computes the Flux well-filling variables:

    D_reg
    F_fill
    G_absorb
    G_sat
    G_overflow
    well_regime

and tests whether these beat the mass-only baseline.

Critical benchmark from current Pantheon+ baseline
--------------------------------------------------
    BIC intercept-only ≈ -5948.39
    BIC mass-only     ≈ -5941.02
    BIC mass-step     ≈ -5941.75

Flux kill condition
-------------------
A serious positive result requires:

1. G_sat or capacity-state variables beat mass-only by BIC.
2. Preferably beat the intercept-only/null BIC.
3. Fixed-mass test remains nonzero:
   within the same logMstar bin, high-F_fill hosts differ from low-F_fill hosts.
4. Missing proxy values do not dominate the result.

Outputs
-------
flux_neighborhood_merged.csv
flux_neighborhood_results.csv
flux_neighborhood_summary.csv
flux_neighborhood_regimes.csv
flux_neighborhood_binned.csv
flux_neighborhood_fixed_mass.csv
flux_neighborhood_correlations.csv
flux_neighborhood_audit.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EPS = 1e-12
AGE_UNIVERSE = 13.8
H0_GLOBAL = 67.4

# Official current baseline targets from flux_pantheon_mass_baseline.py run.
BIC_INTERCEPT_BASELINE = -5948.386225
BIC_MASS_ONLY_BASELINE = -5941.017977
BIC_MASS_STEP_BASELINE = -5941.746913

ALIASES = {
    "SNID": ["SNID", "CID", "sn_id", "SN", "name", "IAU_NAME"],
    "logMstar": ["logMstar", "HOST_LOGMASS", "MHOST", "logM", "host_logmass", "HOST_LOGMASS"],
    "Reff_kpc": ["Reff_kpc", "Reff", "R_eff", "host_reff_kpc", "Rh", "R50_kpc", "half_light_radius_kpc"],
    "stellar_age_Gyr": ["stellar_age_Gyr", "age_Gyr", "host_age", "progenitor_age", "TITAN_age", "age"],
    "sigma_v_km_s": ["sigma_v_km_s", "sigma_v", "velocity_dispersion", "VDISP", "sigma"],
    "sSFR": ["sSFR", "ssfr", "host_sSFR"],
    "logsSFR": ["logsSFR", "log_sSFR", "LOGSFR", "host_logsSFR"],
    "metallicity": ["metallicity", "host_metallicity", "Z", "OH12", "12logOH"],
    "morphology": ["morphology", "host_type", "type", "host_morphology", "TType"],
}


def find_col(df: pd.DataFrame, canonical: str):
    for c in ALIASES.get(canonical, [canonical]):
        if c in df.columns:
            return c
    return None


def standardize_id(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().str.replace(r"\s+", "", regex=True)


def read_csv_any(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def normalize_proxy_table(host: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=host.index)
    for key in ["SNID", "logMstar", "Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s",
                "sSFR", "logsSFR", "metallicity", "morphology"]:
        col = find_col(host, key)
        out[key] = host[col] if col is not None else np.nan

    out["SNID"] = standardize_id(out["SNID"])
    for c in ["logMstar", "Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s", "sSFR", "logsSFR", "metallicity"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    # If sSFR is provided but logsSFR missing, make logsSFR.
    if out["logsSFR"].isna().all() and out["sSFR"].notna().any():
        if out["sSFR"].median(skipna=True) < 0:
            out["logsSFR"] = out["sSFR"]
            out["sSFR"] = 10 ** out["logsSFR"]
        else:
            out["logsSFR"] = np.log10(np.maximum(out["sSFR"], EPS))
    elif out["sSFR"].isna().all() and out["logsSFR"].notna().any():
        out["sSFR"] = 10 ** out["logsSFR"]

    out["morphology"] = out["morphology"].astype(str)
    return out


def normalize_baseline_table(base: pd.DataFrame) -> pd.DataFrame:
    out = base.copy()
    if "SNID" not in out.columns:
        if "CID" in out.columns:
            out["SNID"] = out["CID"]
        else:
            raise ValueError("Baseline file must contain SNID or CID.")
    out["SNID"] = standardize_id(out["SNID"])
    if "logMstar" not in out.columns:
        if "HOST_LOGMASS" in out.columns:
            out["logMstar"] = out["HOST_LOGMASS"]
        else:
            raise ValueError("Baseline file must contain logMstar or HOST_LOGMASS.")
    if "mu_resid_centered" not in out.columns:
        raise ValueError("Baseline file must contain mu_resid_centered from flux_pantheon_mass_baseline.py.")
    return out


def merge_baseline_with_proxies(base: pd.DataFrame, proxies: pd.DataFrame) -> pd.DataFrame:
    base_n = normalize_baseline_table(base)
    prox_n = normalize_proxy_table(proxies)

    merged = pd.merge(base_n, prox_n, on="SNID", how="left", suffixes=("", "_proxy"))

    # Coalesce proxy mass if baseline missing.
    if "logMstar_proxy" in merged.columns:
        merged["logMstar"] = merged["logMstar"].fillna(merged["logMstar_proxy"])

    return merged


def impute_and_flag(df: pd.DataFrame, allow_impute: bool = True) -> pd.DataFrame:
    out = df.copy()
    required_proxy_cols = ["Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s", "logsSFR", "metallicity"]
    defaults = {
        "Reff_kpc": 5.0,
        "stellar_age_Gyr": 4.0,
        "sigma_v_km_s": 150.0,
        "logsSFR": -10.0,
        "metallicity": 8.7,
    }

    for c in required_proxy_cols:
        if c not in out.columns:
            out[c] = np.nan
        out[f"{c}_was_missing"] = out[c].isna()
        if allow_impute:
            fallback = out[c].median(skipna=True) if out[c].notna().any() else defaults[c]
            out[c] = out[c].fillna(fallback)

    if not allow_impute:
        out = out.dropna(subset=["Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s"]).copy()

    return out


def normalize_to_median(x):
    x = np.asarray(x, dtype=float)
    return x / max(np.nanmedian(x), EPS)


def add_capacity_variables(df: pd.DataFrame, proxy="hybrid", Dsat=1.15, q_sat=7.0) -> pd.DataFrame:
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

    F = D / max(Dsat, EPS)
    G_sat = F ** q_sat / (1.0 + F ** q_sat)
    excess = np.maximum(F - 1.0, 0.0)
    G_overflow = (excess / (1.0 + excess)) * G_sat
    G_absorb = 1.0 - G_sat

    out["D_MR_age"] = D_MR_age
    out["D_sigma_age"] = D_sigma_age
    out["D_hybrid"] = D_hybrid
    out["D_reg"] = D
    out["logD_reg"] = np.log10(np.maximum(D, EPS))
    out["F_fill"] = F
    out["G_absorb"] = G_absorb
    out["G_sat"] = G_sat
    out["G_overflow"] = G_overflow
    out["well_regime"] = np.select(
        [F < 0.75, (F >= 0.75) & (F < 1.25), F >= 1.25],
        ["filling/absorbing", "near-capacity/saturated", "overflowing"],
        default="unknown",
    )
    out["H0_proxy_centered"] = H0_GLOBAL * 10 ** (-out["mu_resid_centered"].to_numpy(dtype=float) / 5.0)
    return out


def design_matrix(df: pd.DataFrame, predictors: List[str]):
    X = [np.ones(len(df))]
    names = ["intercept"]
    for p in predictors:
        vals = pd.to_numeric(df[p], errors="coerce").to_numpy(float)
        mu, sd = np.nanmean(vals), np.nanstd(vals)
        if not np.isfinite(sd) or sd <= 0:
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
    rss = float(np.sum(resid ** 2))
    tss = float(np.sum((y - y.mean()) ** 2))
    rmse = float(np.sqrt(rss / max(n, 1)))
    r2 = 1.0 - rss / max(tss, EPS)
    aic = n * np.log(max(rss / max(n, 1), EPS)) + 2 * k
    bic = n * np.log(max(rss / max(n, 1), EPS)) + k * np.log(max(n, 1))

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
    for nm, b in zip(names, beta):
        out[f"coef_{nm}"] = b
    return out


def run_models(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.copy()
    tmp["mass_high_10"] = (tmp["logMstar"] >= 10.0).astype(float)
    models = {
        "M0_intercept_only": [],
        "M1_mass_only": ["logMstar"],
        "M2_mass_step_10": ["mass_high_10"],
        "M3_logD_reg": ["logD_reg"],
        "M4_Gsat": ["G_sat"],
        "M5_capacity_state": ["G_absorb", "G_sat", "G_overflow"],
        "M6_capacity_controls": ["G_sat", "G_overflow", "logsSFR", "metallicity"],
        "M7_kitchen_sink": ["logMstar", "Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s", "logsSFR", "metallicity"],
    }
    return pd.DataFrame([fit_ols(tmp, name, preds) for name, preds in models.items()])


def fixed_mass_test(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.copy()
    tmp["mass_bin"] = pd.cut(tmp["logMstar"], bins=[5, 8.5, 9.5, 10.0, 10.5, 11.0, 13.0], include_lowest=True)
    rows = []
    for mbin, sub in tmp.groupby("mass_bin", observed=False):
        if len(sub) < 20:
            continue
        med = sub["F_fill"].median()
        lowF = sub[sub["F_fill"] <= med]
        highF = sub[sub["F_fill"] > med]
        rows.append({
            "mass_bin": str(mbin),
            "N": len(sub),
            "mean_logMstar": sub["logMstar"].mean(),
            "lowF_mean_mu": lowF["mu_resid_centered"].mean(),
            "highF_mean_mu": highF["mu_resid_centered"].mean(),
            "high_minus_low_mu": highF["mu_resid_centered"].mean() - lowF["mu_resid_centered"].mean(),
            "lowF_H0": lowF["H0_proxy_centered"].mean(),
            "highF_H0": highF["H0_proxy_centered"].mean(),
            "high_minus_low_H0": highF["H0_proxy_centered"].mean() - lowF["H0_proxy_centered"].mean(),
        })
    return pd.DataFrame(rows)


def regime_summary(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby("well_regime", observed=False).agg(
        N=("SNID", "size"),
        mean_logMstar=("logMstar", "mean"),
        mean_Reff=("Reff_kpc", "mean"),
        mean_age=("stellar_age_Gyr", "mean"),
        mean_sigma=("sigma_v_km_s", "mean"),
        mean_Ffill=("F_fill", "mean"),
        mean_Gsat=("G_sat", "mean"),
        mean_Goverflow=("G_overflow", "mean"),
        mean_mu_resid=("mu_resid_centered", "mean"),
        std_mu_resid=("mu_resid_centered", "std"),
        mean_H0_proxy=("H0_proxy_centered", "mean"),
    ).reset_index()


def binned_table(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.copy()
    tmp["fill_bin"] = pd.qcut(tmp["F_fill"], q=min(7, max(2, len(tmp)//30)), duplicates="drop")
    return tmp.groupby("fill_bin", observed=False).agg(
        N=("SNID", "size"),
        mean_Ffill=("F_fill", "mean"),
        mean_Gsat=("G_sat", "mean"),
        mean_Goverflow=("G_overflow", "mean"),
        mean_mu_resid=("mu_resid_centered", "mean"),
        std_mu_resid=("mu_resid_centered", "std"),
        mean_H0_proxy=("H0_proxy_centered", "mean"),
    ).reset_index()


def correlations(df: pd.DataFrame) -> pd.DataFrame:
    vars_ = ["logMstar", "logD_reg", "F_fill", "G_sat", "G_overflow",
             "Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s", "logsSFR", "metallicity"]
    rows = []
    for v in vars_:
        sub = df[["mu_resid_centered", v]].dropna()
        r = float(np.corrcoef(sub["mu_resid_centered"], sub[v])[0, 1]) if len(sub) > 3 else np.nan
        rows.append({"variable": v, "pearson_r_with_mu_resid": r, "N": len(sub)})
    return pd.DataFrame(rows)


def save_plot(df, results, regimes, fixed, outdir):
    plt.figure(figsize=(15, 11))

    plt.subplot(2, 2, 1)
    plt.scatter(df["F_fill"], df["mu_resid_centered"], s=9, alpha=0.35)
    plt.axvline(1.0, linestyle="--", color="black")
    plt.xscale("log")
    plt.xlabel("F_fill = D_reg/Dsat")
    plt.ylabel("centered residual Δμ")
    plt.title("Neighborhood maturity audit")

    plt.subplot(2, 2, 2)
    order = results.sort_values("BIC")
    plt.bar(order["model"], order["BIC"])
    plt.axhline(BIC_INTERCEPT_BASELINE, linestyle="--", color="black", label="baseline null")
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("BIC")
    plt.title("BIC battle")
    plt.legend()

    plt.subplot(2, 2, 3)
    if len(regimes):
        plt.bar(regimes["well_regime"], regimes["mean_H0_proxy"])
        plt.axhline(H0_GLOBAL, linestyle="--", color="black")
        plt.xticks(rotation=20, ha="right")
    plt.ylabel("mean H0 proxy")
    plt.title("H0 proxy by well regime")

    plt.subplot(2, 2, 4)
    if len(fixed):
        plt.bar(fixed["mass_bin"], fixed["high_minus_low_mu"])
        plt.axhline(0, linestyle="--", color="black")
        plt.xticks(rotation=30, ha="right")
    plt.ylabel("high-F minus low-F Δμ")
    plt.title("Fixed-mass maturity test")

    plt.tight_layout()
    plt.savefig(outdir / "flux_neighborhood_audit.png", dpi=180)
    plt.close()


def make_schema(path: Path):
    schema = pd.DataFrame({
        "SNID": ["1990O", "1990af"],
        "Reff_kpc": [4.8, 7.2],
        "stellar_age_Gyr": [5.6, 2.1],
        "sigma_v_km_s": [180.0, 105.0],
        "logsSFR": [-10.8, -9.5],
        "metallicity": [8.85, 8.62],
        "morphology": ["early/quiescent", "late/star-forming"],
    })
    schema.to_csv(path, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=False, default=None, help="flux_pantheon_mass_baseline_data.csv")
    parser.add_argument("--proxies", required=False, default=None, help="host structure/maturity proxy CSV")
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--proxy", choices=["hybrid", "MR_age", "sigma_age"], default="hybrid")
    parser.add_argument("--Dsat", type=float, default=1.15)
    parser.add_argument("--q-sat", type=float, default=7.0)
    parser.add_argument("--no-impute", action="store_true")
    parser.add_argument("--write-schema-only", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    make_schema(outdir / "flux_neighborhood_host_proxy_schema.csv")

    if args.write_schema_only:
        print(f"Wrote schema to {outdir / 'flux_neighborhood_host_proxy_schema.csv'}")
        return

    if not args.baseline or not args.proxies:
        raise ValueError("Provide --baseline and --proxies, or use --write-schema-only.")

    base = read_csv_any(args.baseline)
    prox = read_csv_any(args.proxies)
    merged = merge_baseline_with_proxies(base, prox)
    merged = impute_and_flag(merged, allow_impute=(not args.no_impute))
    merged = merged.dropna(subset=["SNID", "logMstar", "mu_resid_centered"]).copy()
    merged = add_capacity_variables(merged, proxy=args.proxy, Dsat=args.Dsat, q_sat=args.q_sat)

    results = run_models(merged)
    fixed = fixed_mass_test(merged)
    regimes = regime_summary(merged)
    bins = binned_table(merged)
    corr = correlations(merged)

    bic_mass = float(results.loc[results["model"] == "M1_mass_only", "BIC"].iloc[0])
    bic_gsat = float(results.loc[results["model"] == "M4_Gsat", "BIC"].iloc[0])
    bic_cap = float(results.loc[results["model"] == "M5_capacity_state", "BIC"].iloc[0])
    best_model = str(results.sort_values("BIC").iloc[0]["model"])

    missing_rates = {
        c: float(merged[f"{c}_was_missing"].mean()) if f"{c}_was_missing" in merged.columns else np.nan
        for c in ["Reff_kpc", "stellar_age_Gyr", "sigma_v_km_s", "logsSFR", "metallicity"]
    }

    verdict = "NEIGHBORHOOD_FLUX_SUPPORT" if (
        (bic_gsat < bic_mass or bic_cap < bic_mass) and
        min(bic_gsat, bic_cap) < BIC_INTERCEPT_BASELINE
    ) else "NO_SUPPORT_OR_INCOMPLETE"

    summary = pd.DataFrame([{
        "verdict": verdict,
        "N_rows": len(merged),
        "proxy": args.proxy,
        "Dsat": args.Dsat,
        "q_sat": args.q_sat,
        "best_BIC_model": best_model,
        "baseline_BIC_intercept": BIC_INTERCEPT_BASELINE,
        "baseline_BIC_mass_only": BIC_MASS_ONLY_BASELINE,
        "current_BIC_mass_only": bic_mass,
        "current_BIC_Gsat": bic_gsat,
        "current_BIC_capacity_state": bic_cap,
        "delta_BIC_Gsat_minus_mass": bic_gsat - bic_mass,
        "delta_BIC_capacity_minus_mass": bic_cap - bic_mass,
        "beats_original_null": min(bic_gsat, bic_cap) < BIC_INTERCEPT_BASELINE,
        "Reff_missing_rate": missing_rates["Reff_kpc"],
        "age_missing_rate": missing_rates["stellar_age_Gyr"],
        "sigma_missing_rate": missing_rates["sigma_v_km_s"],
        "sSFR_missing_rate": missing_rates["logsSFR"],
        "metallicity_missing_rate": missing_rates["metallicity"],
        "interpretation": (
            "Capacity/maturity variables beat mass and the original null baseline. Check fixed-mass bins and missingness before claiming physical support."
            if verdict == "NEIGHBORHOOD_FLUX_SUPPORT"
            else
            "Capacity/maturity variables do not yet beat the required baselines, or proxy coverage is incomplete."
        ),
    }])

    merged.to_csv(outdir / "flux_neighborhood_merged.csv", index=False)
    results.to_csv(outdir / "flux_neighborhood_results.csv", index=False)
    summary.to_csv(outdir / "flux_neighborhood_summary.csv", index=False)
    regimes.to_csv(outdir / "flux_neighborhood_regimes.csv", index=False)
    bins.to_csv(outdir / "flux_neighborhood_binned.csv", index=False)
    fixed.to_csv(outdir / "flux_neighborhood_fixed_mass.csv", index=False)
    corr.to_csv(outdir / "flux_neighborhood_correlations.csv", index=False)
    save_plot(merged, results, regimes, fixed, outdir)

    print("\nFlux Neighborhood Audit")
    print("=" * 88)
    print(summary.to_string(index=False))
    print("\nRegression battle")
    print(results[["model", "N", "k_params", "RMSE", "R2", "AIC", "BIC"]].sort_values("BIC").to_string(index=False))
    print("\nRegime summary")
    print(regimes.to_string(index=False))
    print("\nFixed-mass test")
    print(fixed.to_string(index=False))
    print("\nCorrelations")
    print(corr.to_string(index=False))


if __name__ == "__main__":
    main()
