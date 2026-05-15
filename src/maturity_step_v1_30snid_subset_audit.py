#!/usr/bin/env python3
"""
maturity_step_v1_30snid_subset_audit.py

Strict populated-subset audit for Flux Neighborhood / Maturity Step v1.

This script uses only SNIDs with BOTH:
    Reff_kpc
    stellar_age_Gyr

It avoids the full-sample imputation dilution problem.
"""

import numpy as np
import pandas as pd

BASELINE = "flux_pantheon_mass_baseline_data.csv"
PROXIES = "host_structure_maturity_proxies.csv"

AGE_UNIVERSE = 13.8
H0_GLOBAL = 67.4
DSAT = 1.15
Q_SAT = 7.0


def fit_ols(df, name, predictors, ycol="mu_resid_centered"):
    d = df.dropna(subset=[ycol] + predictors).copy()
    y = d[ycol].to_numpy(float)

    X = [np.ones(len(d))]
    for col in predictors:
        v = pd.to_numeric(d[col], errors="coerce").to_numpy(float)
        mu = np.nanmean(v)
        sd = np.nanstd(v)
        if not np.isfinite(sd) or sd <= 0:
            sd = 1.0
        X.append((v - mu) / sd)

    X = np.vstack(X).T
    mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y = y[mask]
    X = X[mask]

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat

    n = len(y)
    k = X.shape[1]
    rss = float(np.sum(resid ** 2))
    tss = float(np.sum((y - y.mean()) ** 2))
    rmse = float(np.sqrt(rss / max(n, 1)))
    r2 = 1.0 - rss / max(tss, 1e-12)
    aic = n * np.log(max(rss / max(n, 1), 1e-12)) + 2 * k
    bic = n * np.log(max(rss / max(n, 1), 1e-12)) + k * np.log(max(n, 1))

    return {
        "model": name,
        "predictors": ",".join(predictors) if predictors else "intercept_only",
        "N": n,
        "k": k,
        "RMSE": rmse,
        "R2": r2,
        "AIC": aic,
        "BIC": bic,
    }


def main():
    baseline = pd.read_csv(BASELINE)
    proxy = pd.read_csv(PROXIES)

    p = proxy[proxy["Reff_kpc"].notna() & proxy["stellar_age_Gyr"].notna()].copy()

    cols = [
        "SNID",
        "host_name",
        "Reff_kpc",
        "stellar_age_Gyr",
        "sigma_v_km_s",
        "logsSFR",
        "metallicity",
        "morphology",
        "notes",
    ]
    cols = [c for c in cols if c in p.columns]

    m = baseline.merge(p[cols], on="SNID", how="inner")

    # Optional values: only fill for computation if missing.
    if "sigma_v_km_s" not in m.columns:
        m["sigma_v_km_s"] = np.nan
    if "logsSFR" not in m.columns:
        m["logsSFR"] = np.nan
    if "metallicity" not in m.columns:
        m["metallicity"] = np.nan

    m["sigma_v_km_s"] = pd.to_numeric(m["sigma_v_km_s"], errors="coerce").fillna(150.0)
    m["logsSFR"] = pd.to_numeric(m["logsSFR"], errors="coerce").fillna(-10.0)
    m["metallicity"] = pd.to_numeric(m["metallicity"], errors="coerce").fillna(8.7)

    M = 10 ** pd.to_numeric(m["logMstar"], errors="coerce").to_numpy(float)
    Re = pd.to_numeric(m["Reff_kpc"], errors="coerce").to_numpy(float)
    age = pd.to_numeric(m["stellar_age_Gyr"], errors="coerce").to_numpy(float)
    sig = pd.to_numeric(m["sigma_v_km_s"], errors="coerce").to_numpy(float)

    D_MR = (M / Re) * (age / AGE_UNIVERSE)
    D_sig = (sig / 150.0) ** 2 * (age / AGE_UNIVERSE)
    D_hybrid = np.sqrt(D_MR * D_sig)

    m["D_reg"] = D_hybrid / np.nanmedian(D_hybrid)
    m["logD_reg"] = np.log10(np.maximum(m["D_reg"], 1e-12))
    m["F_fill"] = m["D_reg"] / DSAT
    m["G_sat"] = m["F_fill"] ** Q_SAT / (1.0 + m["F_fill"] ** Q_SAT)
    m["G_overflow"] = (
        np.maximum(m["F_fill"] - 1.0, 0.0)
        / (1.0 + np.maximum(m["F_fill"] - 1.0, 0.0))
    ) * m["G_sat"]
    m["G_absorb"] = 1.0 - m["G_sat"]
    m["H0_proxy"] = H0_GLOBAL * 10 ** (-m["mu_resid_centered"] / 5.0)

    m["well_regime"] = np.select(
        [
            m["F_fill"] < 0.75,
            (m["F_fill"] >= 0.75) & (m["F_fill"] < 1.25),
            m["F_fill"] >= 1.25,
        ],
        ["filling/absorbing", "near-capacity/saturated", "overflowing"],
        default="unknown",
    )

    m["mass_high_10"] = (m["logMstar"] >= 10.0).astype(float)

    models = pd.DataFrame(
        [
            fit_ols(m, "intercept", []),
            fit_ols(m, "mass_only", ["logMstar"]),
            fit_ols(m, "mass_step_10", ["mass_high_10"]),
            fit_ols(m, "logD_reg", ["logD_reg"]),
            fit_ols(m, "G_sat", ["G_sat"]),
            fit_ols(m, "capacity_state", ["G_absorb", "G_sat", "G_overflow"]),
            fit_ols(m, "capacity_controls", ["G_sat", "G_overflow", "logsSFR", "metallicity"]),
            fit_ols(
                m,
                "kitchen_sink",
                [
                    "logMstar",
                    "Reff_kpc",
                    "stellar_age_Gyr",
                    "sigma_v_km_s",
                    "logsSFR",
                    "metallicity",
                ],
            ),
        ]
    ).sort_values("BIC")

    regimes = (
        m.groupby("well_regime", observed=False)
        .agg(
            N=("SNID", "size"),
            unique_SNIDs=("SNID", "nunique"),
            mean_logMstar=("logMstar", "mean"),
            mean_age=("stellar_age_Gyr", "mean"),
            mean_Reff=("Reff_kpc", "mean"),
            mean_sigma=("sigma_v_km_s", "mean"),
            mean_Ffill=("F_fill", "mean"),
            mean_Gsat=("G_sat", "mean"),
            mean_Goverflow=("G_overflow", "mean"),
            mean_mu=("mu_resid_centered", "mean"),
            std_mu=("mu_resid_centered", "std"),
            mean_H0=("H0_proxy", "mean"),
        )
        .reset_index()
    )

    snid = (
        m.groupby(["SNID", "host_name", "morphology"], dropna=False)
        .agg(
            n_entries=("SNID", "size"),
            logMstar=("logMstar", "mean"),
            Reff_kpc=("Reff_kpc", "mean"),
            stellar_age_Gyr=("stellar_age_Gyr", "mean"),
            sigma_v_km_s=("sigma_v_km_s", "mean"),
            F_fill=("F_fill", "mean"),
            G_sat=("G_sat", "mean"),
            G_overflow=("G_overflow", "mean"),
            mean_mu=("mu_resid_centered", "mean"),
            mean_H0=("H0_proxy", "mean"),
        )
        .reset_index()
        .sort_values("F_fill")
    )

    # Fixed-mass test
    m["mass_bin"] = pd.cut(
        m["logMstar"],
        bins=[5, 8.5, 9.5, 10.0, 10.5, 11.0, 13.0],
        include_lowest=True,
    )

    fixed_rows = []
    for mbin, sub in m.groupby("mass_bin", observed=False):
        if len(sub) < 3:
            continue
        med = sub["F_fill"].median()
        lo = sub[sub["F_fill"] <= med]
        hi = sub[sub["F_fill"] > med]
        fixed_rows.append(
            {
                "mass_bin": str(mbin),
                "N": len(sub),
                "unique_SNIDs": sub["SNID"].nunique(),
                "mean_logMstar": sub["logMstar"].mean(),
                "lowF_mu": lo["mu_resid_centered"].mean(),
                "highF_mu": hi["mu_resid_centered"].mean(),
                "high_minus_low_mu": hi["mu_resid_centered"].mean()
                - lo["mu_resid_centered"].mean(),
                "lowF_H0": lo["H0_proxy"].mean(),
                "highF_H0": hi["H0_proxy"].mean(),
                "high_minus_low_H0": hi["H0_proxy"].mean() - lo["H0_proxy"].mean(),
            }
        )

    fixed = pd.DataFrame(fixed_rows)

    summary = pd.DataFrame(
        [
            {
                "verdict": "MATURITY_STEP_V1_30SNID_STRICT_SUBSET",
                "unique_Reff_Age_SNIDs": p["SNID"].nunique(),
                "matched_Pantheon_rows": len(m),
                "best_BIC_model": models.iloc[0]["model"],
                "BIC_intercept": float(models.loc[models["model"] == "intercept", "BIC"].iloc[0]),
                "BIC_mass_only": float(models.loc[models["model"] == "mass_only", "BIC"].iloc[0]),
                "BIC_Gsat": float(models.loc[models["model"] == "G_sat", "BIC"].iloc[0]),
                "BIC_capacity_state": float(models.loc[models["model"] == "capacity_state", "BIC"].iloc[0]),
                "corr_mu_Ffill": float(np.corrcoef(m["mu_resid_centered"], m["F_fill"])[0, 1]),
                "corr_mu_Gsat": float(np.corrcoef(m["mu_resid_centered"], m["G_sat"])[0, 1]),
                "warning": "Strict subset only; mixed sourced and estimated host values.",
            }
        ]
    )

    m.to_csv("maturity_step_v1_30snid_matched_rows.csv", index=False)
    models.to_csv("maturity_step_v1_30snid_models.csv", index=False)
    regimes.to_csv("maturity_step_v1_30snid_regimes.csv", index=False)
    snid.to_csv("maturity_step_v1_30snid_snid_sequence.csv", index=False)
    fixed.to_csv("maturity_step_v1_30snid_fixed_mass.csv", index=False)
    summary.to_csv("maturity_step_v1_30snid_summary.csv", index=False)

    print("\n30-SNID populated subset audit")
    print("=" * 88)
    print(summary.to_string(index=False))

    print("\nModel BIC ranking:")
    print(models.to_string(index=False))

    print("\nRegime summary:")
    print(regimes.to_string(index=False))

    print("\nSNID-level sequence:")
    print(snid.to_string(index=False))

    print("\nFixed-mass test:")
    print(fixed.to_string(index=False))

    print("\nSaved:")
    print("maturity_step_v1_30snid_matched_rows.csv")
    print("maturity_step_v1_30snid_models.csv")
    print("maturity_step_v1_30snid_regimes.csv")
    print("maturity_step_v1_30snid_snid_sequence.csv")
    print("maturity_step_v1_30snid_fixed_mass.csv")
    print("maturity_step_v1_30snid_summary.csv")


if __name__ == "__main__":
    main()
