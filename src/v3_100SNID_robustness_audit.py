
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def bic_from_rss(n, rss, k):
    if n <= 0 or rss <= 0:
        return np.nan
    return n * np.log(rss / n) + k * np.log(n)


def fit_ols_bic(df, y_col, x_cols, model_name):
    cols = [y_col] + x_cols
    d = df[cols].replace([np.inf, -np.inf], np.nan).dropna()
    n = len(d)
    if n <= len(x_cols) + 1:
        return {
            "model": model_name, "predictors": ",".join(x_cols) if x_cols else "intercept_only",
            "N": n, "k": len(x_cols) + 1, "RMSE": np.nan, "R2": np.nan,
            "AIC": np.nan, "BIC": np.nan
        }

    y = d[y_col].to_numpy(float)
    X = np.ones((n, len(x_cols) + 1), dtype=float)
    for i, col in enumerate(x_cols, start=1):
        X[:, i] = d[col].to_numpy(float)

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    resid = y - pred
    rss = float(np.sum(resid ** 2))
    tss = float(np.sum((y - y.mean()) ** 2))
    k = X.shape[1]
    rmse = float(np.sqrt(rss / n))
    r2 = float(1 - rss / tss) if tss > 0 else np.nan
    aic = float(n * np.log(rss / n) + 2 * k) if rss > 0 else np.nan
    bic = float(bic_from_rss(n, rss, k))
    return {
        "model": model_name, "predictors": ",".join(x_cols) if x_cols else "intercept_only",
        "N": n, "k": k, "RMSE": rmse, "R2": r2, "AIC": aic, "BIC": bic
    }


def classify_regime(f):
    if f < 0.5:
        return "filling/absorbing"
    if f < 1.5:
        return "near-capacity/saturated"
    return "overflowing"


def add_flux_columns(df, dsat=1.15, qsat=7.0):
    d = df.copy()

    # Safety defaults
    for col in ["sigma_v_km_s", "logsSFR", "metallicity"]:
        if col not in d.columns:
            d[col] = np.nan

    # Use robust imputation for sigma only inside the active subset.
    sigma_med = d["sigma_v_km_s"].replace([np.inf, -np.inf], np.nan).median()
    if not np.isfinite(sigma_med):
        sigma_med = 110.0
    d["sigma_eff"] = d["sigma_v_km_s"].fillna(sigma_med)

    # D_MR = (10^logM / Reff) * (Age/13.8), normalized later.
    # Use log space for numerical stability.
    d["logD_MR_raw"] = d["logMstar"] - np.log10(d["Reff_kpc"]) + np.log10(d["stellar_age_Gyr"] / 13.8)

    # D_sigma = (sigma/150)^2 * (Age/13.8)
    d["logD_sigma_raw"] = 2 * np.log10(d["sigma_eff"] / 150.0) + np.log10(d["stellar_age_Gyr"] / 13.8)

    # Hybrid geometric mean
    d["logD_hybrid_raw"] = 0.5 * (d["logD_MR_raw"] + d["logD_sigma_raw"])

    # Normalize by median hybrid depth inside the current active audit sample.
    median_logD = d["logD_hybrid_raw"].replace([np.inf, -np.inf], np.nan).median()
    d["D_reg"] = 10 ** (d["logD_hybrid_raw"] - median_logD)
    d["F_fill"] = d["D_reg"] / dsat

    # Gates
    d["G_sat"] = (d["F_fill"] ** qsat) / (1.0 + d["F_fill"] ** qsat)
    d["G_overflow"] = np.clip(d["F_fill"] - 1.0, 0.0, None) / (1.0 + np.clip(d["F_fill"] - 1.0, 0.0, None))
    d["G_absorb"] = 1.0 - d["G_sat"]

    d["well_regime"] = d["F_fill"].apply(classify_regime)
    d["mass_step_10"] = (d["logMstar"] >= 10.0).astype(int)

    return d


def fixed_mass_table(df):
    d = df.copy()
    bins = [0, 8.5, 9.5, 10.0, 10.5, 11.0, 13.5]
    d["mass_bin"] = pd.cut(d["logMstar"], bins=bins, include_lowest=True)
    out = []
    for mb, g in d.groupby("mass_bin", observed=False):
        if len(g) < 2:
            continue
        med = g["F_fill"].median()
        low = g[g["F_fill"] <= med]
        high = g[g["F_fill"] > med]
        out.append({
            "mass_bin": str(mb),
            "N": len(g),
            "unique_SNIDs": g["SNID"].astype(str).nunique(),
            "mean_logMstar": g["logMstar"].mean(),
            "lowF_mean_mu": low["mean_mu_resid"].mean(),
            "highF_mean_mu": high["mean_mu_resid"].mean(),
            "high_minus_low_mu": high["mean_mu_resid"].mean() - low["mean_mu_resid"].mean(),
            "lowF_H0": low["mean_H0_proxy"].mean(),
            "highF_H0": high["mean_H0_proxy"].mean(),
            "high_minus_low_H0": high["mean_H0_proxy"].mean() - low["mean_H0_proxy"].mean(),
        })
    return pd.DataFrame(out)


def run_audit(input_path, outdir, dsat=1.15, qsat=7.0):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(input_path)
    active = raw.dropna(subset=["Reff_kpc", "stellar_age_Gyr", "logMstar", "mean_mu_resid", "mean_H0_proxy"]).copy()
    active = active.drop_duplicates(subset=["SNID"], keep="first").copy()

    active = add_flux_columns(active, dsat=dsat, qsat=qsat)

    # Model battle
    candidate_models = [
        ("intercept", []),
        ("mass_only", ["logMstar"]),
        ("mass_step_10", ["mass_step_10"]),
        ("logD_hybrid", ["logD_hybrid_raw"]),
        ("G_sat", ["G_sat"]),
        ("capacity_state", ["G_absorb", "G_sat", "G_overflow"]),
        ("capacity_controls", ["G_sat", "G_overflow", "logsSFR", "metallicity"]),
        ("kitchen_sink", ["logMstar", "Reff_kpc", "stellar_age_Gyr", "sigma_eff", "logsSFR", "metallicity"]),
    ]

    model_rows = [fit_ols_bic(active, "mean_mu_resid", xs, name) for name, xs in candidate_models]
    models = pd.DataFrame(model_rows).sort_values("BIC", ascending=True)

    regimes = (
        active.groupby("well_regime", observed=False)
        .agg(
            N=("SNID", "size"),
            mean_logMstar=("logMstar", "mean"),
            mean_age=("stellar_age_Gyr", "mean"),
            mean_Reff=("Reff_kpc", "mean"),
            mean_sigma=("sigma_eff", "mean"),
            mean_Ffill=("F_fill", "mean"),
            mean_Gsat=("G_sat", "mean"),
            mean_Goverflow=("G_overflow", "mean"),
            mean_mu=("mean_mu_resid", "mean"),
            std_mu=("mean_mu_resid", "std"),
            mean_H0=("mean_H0_proxy", "mean"),
        )
        .reset_index()
    )

    fixed_mass = fixed_mass_table(active)

    corr_vars = ["logMstar", "logD_hybrid_raw", "D_reg", "F_fill", "G_sat", "G_overflow",
                 "Reff_kpc", "stellar_age_Gyr", "sigma_eff", "logsSFR", "metallicity"]
    corr_rows = []
    for v in corr_vars:
        sub = active[["mean_mu_resid", v]].replace([np.inf, -np.inf], np.nan).dropna()
        corr_rows.append({
            "variable": v,
            "pearson_r_with_mu_resid": sub["mean_mu_resid"].corr(sub[v]) if len(sub) > 2 else np.nan,
            "N": len(sub),
        })
    corrs = pd.DataFrame(corr_rows)

    # Robustness filters
    filters = {
        "all_active": active,
        "outlier_clipped_45_90": active[(active["mean_H0_proxy"] > 45) & (active["mean_H0_proxy"] < 90)],
        "high_confidence_non_estimate": active[active.get("value_confidence", "").astype(str).str.lower() != "estimate"],
        "calibrators_only": active[active.get("is_calibrator", 0) == 1],
        "hubble_flow_only": active[active.get("is_calibrator", 0) == 0],
    }

    robustness_rows = []
    for name, sub in filters.items():
        if len(sub) == 0:
            continue
        rsum = sub.groupby("well_regime", observed=False)["mean_H0_proxy"].mean().to_dict()
        fm = fixed_mass_table(sub)
        positive_bins = int((fm["high_minus_low_H0"] > 0).sum()) if len(fm) else 0
        tested_bins = int(fm["high_minus_low_H0"].notna().sum()) if len(fm) else 0
        robustness_rows.append({
            "filter": name,
            "N": len(sub),
            "unique_SNIDs": sub["SNID"].astype(str).nunique(),
            "filling_H0": rsum.get("filling/absorbing", np.nan),
            "saturated_H0": rsum.get("near-capacity/saturated", np.nan),
            "overflowing_H0": rsum.get("overflowing", np.nan),
            "adult_minus_child_H0": rsum.get("overflowing", np.nan) - rsum.get("filling/absorbing", np.nan),
            "fixed_mass_positive_bins": positive_bins,
            "fixed_mass_tested_bins": tested_bins,
        })
    robustness = pd.DataFrame(robustness_rows)

    # Summary/verdict
    best = models.iloc[0].to_dict()
    bics = {r["model"]: r["BIC"] for _, r in models.iterrows()}
    delta_gsat_mass = bics.get("G_sat", np.nan) - bics.get("mass_only", np.nan)
    delta_hybrid_mass = bics.get("logD_hybrid", np.nan) - bics.get("mass_only", np.nan)
    adult_minus_child = (
        regimes.set_index("well_regime").loc["overflowing", "mean_H0"]
        - regimes.set_index("well_regime").loc["filling/absorbing", "mean_H0"]
        if set(["overflowing", "filling/absorbing"]).issubset(set(regimes["well_regime"]))
        else np.nan
    )

    fixed_positive = int((fixed_mass["high_minus_low_H0"] > 0).sum()) if len(fixed_mass) else 0
    fixed_tested = int(fixed_mass["high_minus_low_H0"].notna().sum()) if len(fixed_mass) else 0

    verdict = "MIXED_SUPPORT"
    if adult_minus_child > 0 and fixed_positive >= max(1, fixed_tested // 2):
        verdict = "DIRECTIONAL_SUPPORT"
    if (best["model"] in ["G_sat", "logD_hybrid", "capacity_state", "capacity_controls", "kitchen_sink"]) and adult_minus_child > 0:
        verdict = "MODEL_SUPPORT"
    if (delta_gsat_mass < -10 or delta_hybrid_mass < -10) and adult_minus_child > 0:
        verdict = "DECISIVE_PROXY_SUPPORT"

    summary = pd.DataFrame([{
        "verdict": verdict,
        "input_file": str(input_path),
        "unique_active_SNIDs": active["SNID"].astype(str).nunique(),
        "Dsat": dsat,
        "q_sat": qsat,
        "best_BIC_model": best["model"],
        "best_BIC": best["BIC"],
        "BIC_intercept": bics.get("intercept", np.nan),
        "BIC_mass_only": bics.get("mass_only", np.nan),
        "BIC_logD_hybrid": bics.get("logD_hybrid", np.nan),
        "BIC_Gsat": bics.get("G_sat", np.nan),
        "BIC_capacity_state": bics.get("capacity_state", np.nan),
        "delta_BIC_Gsat_minus_mass": delta_gsat_mass,
        "delta_BIC_logD_minus_mass": delta_hybrid_mass,
        "adult_minus_child_H0": adult_minus_child,
        "fixed_mass_positive_bins": fixed_positive,
        "fixed_mass_tested_bins": fixed_tested,
    }])

    # Outputs
    active.to_csv(outdir / "v3_100snid_active_with_flux_columns.csv", index=False)
    models.to_csv(outdir / "v3_100snid_model_battle.csv", index=False)
    regimes.to_csv(outdir / "v3_100snid_regime_summary.csv", index=False)
    fixed_mass.to_csv(outdir / "v3_100snid_fixed_mass.csv", index=False)
    corrs.to_csv(outdir / "v3_100snid_correlations.csv", index=False)
    robustness.to_csv(outdir / "v3_100snid_robustness_filters.csv", index=False)
    summary.to_csv(outdir / "v3_100snid_summary.csv", index=False)

    print("\nFlux Neighborhood Audit v3 / 100-SNID")
    print("=" * 88)
    print(summary.to_string(index=False))

    print("\nModel BIC ranking:")
    print(models.to_string(index=False))

    print("\nRegime summary:")
    print(regimes.to_string(index=False))

    print("\nFixed-mass test:")
    print(fixed_mass.to_string(index=False))

    print("\nCorrelations:")
    print(corrs.to_string(index=False))

    print("\nRobustness filters:")
    print(robustness.to_string(index=False))

    print("\nSaved outputs to:", outdir)
    return summary, models, regimes, fixed_mass, corrs, robustness


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--proxies", default="host_structure_maturity_proxies_v3_100SNID_FINAL.csv")
    parser.add_argument("--outdir", default="neighborhood_results_v3_100SNID")
    parser.add_argument("--Dsat", type=float, default=1.15)
    parser.add_argument("--q_sat", type=float, default=7.0)
    args = parser.parse_args()

    run_audit(args.proxies, args.outdir, dsat=args.Dsat, qsat=args.q_sat)
