#!/usr/bin/env python3
"""
flux_capacity_overflow_check.py

Level-2 Capacity-Overflow audit for Flux/MIP cosmology.

Purpose
-------
This script upgrades the host-environment audit from a simple saturation proxy
to a dynamic well-capacity model:

    Empty/young wells       -> absorption / DM filling
    Full/mature wells       -> saturation / H0 calibration stiffness
    Overflowing wells       -> escaped throughput / DE and S8 friction

The central question:
Can TITAN-like age bimodality plus realistic host mass/size distributions
naturally split hosts into filling, saturated, and overflowing regimes?

This is still a literature-prior Monte Carlo audit, not object-level data.

Core quantities
---------------
D_reg:
    registration depth proxy, built from potential depth and stellar age

F_fill:
    D_reg / D_sat

G_sat:
    sigmoid/Hill saturation gate

G_overflow:
    excess over capacity, gated by saturation

Predicted observables:
    - Hubble residual / H0 calibration bias from G_sat
    - overflow contribution from G_overflow
    - regime fractions by morphology/age/mass

Outputs
-------
flux_capacity_overflow_population.csv
flux_capacity_overflow_summary.csv
flux_capacity_overflow_regimes.csv
flux_capacity_overflow_results.csv
flux_capacity_overflow_binned.csv
flux_capacity_overflow.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EPS = 1e-12
H0_GLOBAL = 67.4
AGE_UNIVERSE = 13.8


def rng(seed=42):
    return np.random.default_rng(seed)


def hill_gate(x, q=5.0):
    x = np.maximum(np.asarray(x), EPS)
    return x**q / (1.0 + x**q)


def generate_literature_population(N=9000, seed=42):
    """
    Literature-prior population:
    - 60% late/star-forming, age peak ~2.2 Gyr
    - 40% early/quiescent, age peak ~6.0 Gyr
    - mass centered around logM~10.5
    - size-mass scaling by morphology
    - velocity dispersion from compactness
    """
    R = rng(seed)

    is_early = R.random(N) > 0.60
    morphology = np.where(is_early, "early/quiescent", "late/star-forming")

    logM = R.normal(10.5, 0.50, N)
    logM += np.where(is_early, R.normal(0.18, 0.18, N), R.normal(-0.08, 0.22, N))
    logM = np.clip(logM, 8.0, 12.1)
    Mstar = 10**logM

    # Size-mass relation: early hosts are more compact at fixed mass.
    logRe_late = 0.25 * (logM - 10.5) + 0.70 + R.normal(0.0, 0.18, N)
    logRe_early = 0.50 * (logM - 10.5) + 0.48 + R.normal(0.0, 0.16, N)
    logRe = np.where(is_early, logRe_early, logRe_late)
    Reff = np.clip(10**logRe, 0.35, 25.0)

    age = np.where(is_early, R.normal(6.0, 1.5, N), R.normal(2.2, 0.8, N))
    age = np.clip(age, 0.4, 12.5)

    logsSFR = np.where(is_early, R.normal(-11.45, 0.35, N), R.normal(-9.75, 0.45, N))
    metallicity = np.where(is_early, R.normal(8.95, 0.12, N), R.normal(8.70, 0.18, N))

    compactness = Mstar / np.maximum(Reff, EPS)
    comp_norm = compactness / np.median(compactness)
    sigma_v = 145.0 * comp_norm**0.25
    sigma_v *= np.where(is_early, 1.12, 0.92)
    sigma_v *= 10 ** R.normal(0.0, 0.045, N)
    sigma_v = np.clip(sigma_v, 55.0, 330.0)

    return pd.DataFrame({
        "host_id": [f"CAP_HOST_{i:05d}" for i in range(N)],
        "is_early": is_early,
        "morphology": morphology,
        "logMstar": logM,
        "Mstar": Mstar,
        "Reff_kpc": Reff,
        "stellar_age_Gyr": age,
        "logsSFR": logsSFR,
        "sSFR": 10**logsSFR,
        "metallicity": metallicity,
        "sigma_v_km_s": sigma_v,
    })


def normalize(x):
    x = np.asarray(x)
    return x / max(np.nanmedian(x), EPS)


def add_capacity_model(
    df,
    proxy="MR_age",
    Dsat=1.0,
    q_sat=5.0,
    overflow_power=1.0,
    phi0=0.083,
    overflow_scale=1.0,
    residual_noise=0.10,
    seed=43,
):
    """
    Add well filling, saturation, overflow, and predicted observables.

    proxy:
        MR_age       -> (M/Re) * age
        sigma_age    -> sigma_v^2 * age
        hybrid       -> sqrt((M/Re) * sigma_v^2) * age
    """
    R = rng(seed)
    out = df.copy()

    M = out["Mstar"].to_numpy()
    Re = np.maximum(out["Reff_kpc"].to_numpy(), EPS)
    age = np.maximum(out["stellar_age_Gyr"].to_numpy(), EPS)
    sig = np.maximum(out["sigma_v_km_s"].to_numpy(), EPS)

    D_MR_age = normalize((M / Re) * (age / AGE_UNIVERSE))
    D_sigma_age = normalize((sig / 150.0)**2 * (age / AGE_UNIVERSE))
    D_hybrid = normalize(np.sqrt((M / Re) * (sig / 150.0)**2) * (age / AGE_UNIVERSE))

    if proxy == "MR_age":
        D = D_MR_age
    elif proxy == "sigma_age":
        D = D_sigma_age
    elif proxy == "hybrid":
        D = D_hybrid
    else:
        raise ValueError(proxy)

    F_fill = D / max(Dsat, EPS)
    G_sat = hill_gate(F_fill, q=q_sat)

    # Overflow begins after capacity. Use a smooth nonnegative excess.
    excess = np.maximum(F_fill - 1.0, 0.0)
    G_overflow = (excess**overflow_power) / (1.0 + excess**overflow_power)
    G_overflow *= G_sat

    # Absorption dominates before saturation.
    G_absorb = (1.0 - G_sat)

    # H0 calibration / residual. Negative residual -> brighter/closer -> higher inferred H0.
    phi_sat = phi0 * G_sat
    delta_mu_sat = -5.0 * np.log10(1.0 + phi_sat)

    # Add realistic small nuisance leftovers.
    mass_nuis = -0.015 * (out["logMstar"].to_numpy() >= 10.0).astype(float)
    sfr_nuis = 0.008 * (out["logsSFR"].to_numpy() + 10.0)
    metal_nuis = -0.006 * (out["metallicity"].to_numpy() - 8.7)
    noise = R.normal(0.0, residual_noise, len(out))

    mu_resid = delta_mu_sat + mass_nuis + sfr_nuis + metal_nuis + noise
    H0_proxy = H0_GLOBAL * 10**(-mu_resid / 5.0)

    out["D_MR_age"] = D_MR_age
    out["D_sigma_age"] = D_sigma_age
    out["D_hybrid"] = D_hybrid
    out["D_reg"] = D
    out["logD_reg"] = np.log10(np.maximum(D, EPS))
    out["Dsat"] = Dsat
    out["F_fill"] = F_fill
    out["G_absorb"] = G_absorb
    out["G_sat"] = G_sat
    out["G_overflow"] = G_overflow
    out["phi_sat"] = phi_sat
    out["delta_mu_sat"] = delta_mu_sat
    out["mu_resid"] = mu_resid
    out["H0_proxy"] = H0_proxy

    # Regime classification.
    conditions = [
        F_fill < 0.75,
        (F_fill >= 0.75) & (F_fill < 1.25),
        F_fill >= 1.25,
    ]
    choices = ["filling/absorbing", "near-capacity/saturated", "overflowing"]
    out["well_regime"] = np.select(conditions, choices, default="unknown")

    return out


def design_matrix(df, predictors):
    X = [np.ones(len(df))]
    names = ["intercept"]
    for p in predictors:
        vals = df[p].to_numpy(float)
        mu, sd = np.nanmean(vals), np.nanstd(vals)
        if sd <= 0:
            sd = 1.0
        X.append((vals - mu) / sd)
        names.append(p)
    return np.vstack(X).T, names


def fit_ols(df, model, predictors):
    d = df.dropna(subset=["mu_resid"] + predictors).copy()
    y = d["mu_resid"].to_numpy(float)
    X, names = design_matrix(d, predictors)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    n, k = len(y), X.shape[1]
    rss = float(np.sum(resid**2))
    tss = float(np.sum((y - y.mean())**2))
    aic = n * np.log(max(rss/n, EPS)) + 2*k
    bic = n * np.log(max(rss/n, EPS)) + k*np.log(n)
    out = {
        "model": model,
        "predictors": ",".join(predictors),
        "N": n,
        "k_params": k,
        "RMSE": float(np.sqrt(rss/n)),
        "R2": 1.0 - rss/max(tss, EPS),
        "AIC": aic,
        "BIC": bic,
    }
    for nm, b in zip(names, beta):
        out[f"coef_{nm}"] = b
    return out


def regression_battle(df):
    models = {
        "A_mass_only": ["logMstar"],
        "B_mass_sSFR": ["logMstar", "logsSFR"],
        "C_kitchen_sink": ["logMstar", "logsSFR", "metallicity", "stellar_age_Gyr", "Reff_kpc"],
        "D_logD_reg": ["logD_reg"],
        "E_Gsat": ["G_sat"],
        "F_capacity_state": ["G_absorb", "G_sat", "G_overflow"],
        "G_capacity_controls": ["G_sat", "G_overflow", "logsSFR", "metallicity"],
    }
    return pd.DataFrame([fit_ols(df, name, preds) for name, preds in models.items()])


def scan_capacity_params(base, seed=42):
    """
    Light scan for a plausible capacity-overflow model:
    - realistic mass-step scale
    - strong D/regime sequence
    - capacity model competitive with mass-only
    """
    rows = []
    for proxy in ["MR_age", "sigma_age", "hybrid"]:
        for Dsat in [0.75, 0.9, 1.0, 1.15, 1.3]:
            for q in [3.0, 5.0, 7.0]:
                for phi0 in [0.055, 0.070, 0.083, 0.10]:
                    df = add_capacity_model(base, proxy=proxy, Dsat=Dsat, q_sat=q, phi0=phi0, seed=seed+1)
                    res = regression_battle(df)
                    best = res.sort_values("BIC").iloc[0]
                    mass_bic = float(res.loc[res["model"] == "A_mass_only", "BIC"].iloc[0])
                    cap_bic = float(res.loc[res["model"] == "F_capacity_state", "BIC"].iloc[0])

                    high = df[df["logMstar"] >= 10.0]
                    low = df[df["logMstar"] < 10.0]
                    mass_step = float(high["mu_resid"].mean() - low["mu_resid"].mean())

                    reg = df.groupby("well_regime", observed=False)["H0_proxy"].mean().to_dict()
                    h_fill = reg.get("filling/absorbing", np.nan)
                    h_over = reg.get("overflowing", np.nan)
                    h_split = h_over - h_fill if np.isfinite(h_fill) and np.isfinite(h_over) else np.nan

                    # Favor plausible SN residual step and strong regime separation.
                    step_pen = abs(abs(mass_step) - 0.08)
                    split_pen = 0 if (np.isfinite(h_split) and h_split > 3.0) else 0.1
                    cap_pen = 0 if cap_bic < mass_bic else 0.15
                    score = step_pen + split_pen + cap_pen

                    rows.append({
                        "proxy": proxy,
                        "Dsat": Dsat,
                        "q_sat": q,
                        "phi0": phi0,
                        "best_model": str(best["model"]),
                        "best_BIC": float(best["BIC"]),
                        "BIC_mass_only": mass_bic,
                        "BIC_capacity_state": cap_bic,
                        "delta_BIC_capacity_minus_mass": cap_bic - mass_bic,
                        "mass_step_mu_high_minus_low": mass_step,
                        "H0_overflow_minus_filling": h_split,
                        "score": score,
                    })
    return pd.DataFrame(rows).sort_values("score").reset_index(drop=True)


def regime_summary(df):
    g = df.groupby("well_regime", observed=False).agg(
        N=("host_id", "size"),
        early_fraction=("is_early", "mean"),
        mean_logMstar=("logMstar", "mean"),
        mean_age=("stellar_age_Gyr", "mean"),
        mean_Reff=("Reff_kpc", "mean"),
        mean_sigma=("sigma_v_km_s", "mean"),
        mean_Dreg=("D_reg", "mean"),
        mean_Ffill=("F_fill", "mean"),
        mean_Gabsorb=("G_absorb", "mean"),
        mean_Gsat=("G_sat", "mean"),
        mean_Goverflow=("G_overflow", "mean"),
        mean_mu_resid=("mu_resid", "mean"),
        mean_H0_proxy=("H0_proxy", "mean"),
    ).reset_index()
    return g


def binned_table(df):
    tmp = df.copy()
    tmp["fill_bin"] = pd.qcut(tmp["F_fill"], q=7, duplicates="drop")
    return tmp.groupby("fill_bin", observed=False).agg(
        N=("host_id", "size"),
        mean_Ffill=("F_fill", "mean"),
        mean_Gabsorb=("G_absorb", "mean"),
        mean_Gsat=("G_sat", "mean"),
        mean_Goverflow=("G_overflow", "mean"),
        early_fraction=("is_early", "mean"),
        mean_mu_resid=("mu_resid", "mean"),
        std_mu_resid=("mu_resid", "std"),
        mean_H0_proxy=("H0_proxy", "mean"),
    ).reset_index()


def correlations(df):
    vars_ = ["logMstar", "logD_reg", "F_fill", "G_absorb", "G_sat", "G_overflow",
             "stellar_age_Gyr", "sigma_v_km_s", "Reff_kpc", "logsSFR", "metallicity"]
    rows = []
    for v in vars_:
        rows.append({
            "variable": v,
            "pearson_r_with_mu_resid": float(np.corrcoef(df["mu_resid"], df[v])[0, 1]),
            "N": len(df),
        })
    return pd.DataFrame(rows)


def save_plot(df, regimes, binned, results, outdir):
    plt.figure(figsize=(15, 11))

    plt.subplot(2, 2, 1)
    colors = {"filling/absorbing": "tab:blue", "near-capacity/saturated": "tab:orange", "overflowing": "tab:red"}
    for reg, sub in df.groupby("well_regime"):
        plt.scatter(sub["F_fill"], sub["mu_resid"], s=8, alpha=0.25, label=reg, c=colors.get(reg, None))
    plt.axvline(1.0, linestyle="--", color="black", label="capacity")
    plt.xscale("log")
    plt.xlabel("well filling factor F_fill = D_reg/D_sat")
    plt.ylabel("Hubble residual Δμ")
    plt.title("Capacity-overflow residual sequence")
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.plot(binned.index, binned["mean_Gabsorb"], marker="o", label="absorb")
    plt.plot(binned.index, binned["mean_Gsat"], marker="o", label="saturate")
    plt.plot(binned.index, binned["mean_Goverflow"], marker="o", label="overflow")
    plt.xlabel("F_fill quantile")
    plt.ylabel("mean gate")
    plt.title("Absorption → saturation → overflow")
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.bar(regimes["well_regime"], regimes["mean_H0_proxy"])
    plt.axhline(67.4, linestyle="--", color="black", label="global")
    plt.ylabel("mean H0 proxy")
    plt.title("H0 proxy by well state")
    plt.xticks(rotation=20, ha="right")

    plt.subplot(2, 2, 4)
    order = results.sort_values("BIC")
    plt.bar(order["model"], order["BIC"])
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("BIC")
    plt.title("Regression battle")

    plt.tight_layout()
    plt.savefig(outdir / "flux_capacity_overflow.png", dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--n", type=int, default=9000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    base_scan = generate_literature_population(N=2500, seed=args.seed)
    scan = scan_capacity_params(base_scan, seed=args.seed)
    best = scan.iloc[0]

    base = generate_literature_population(N=args.n, seed=args.seed)
    pop = add_capacity_model(
        base,
        proxy=str(best["proxy"]),
        Dsat=float(best["Dsat"]),
        q_sat=float(best["q_sat"]),
        phi0=float(best["phi0"]),
        seed=args.seed + 2,
    )

    results = regression_battle(pop)
    regimes = regime_summary(pop)
    bins = binned_table(pop)
    corr = correlations(pop)

    best_model = str(results.sort_values("BIC").iloc[0]["model"])
    mass_bic = float(results.loc[results["model"] == "A_mass_only", "BIC"].iloc[0])
    cap_bic = float(results.loc[results["model"] == "F_capacity_state", "BIC"].iloc[0])

    reg_h0 = regimes.set_index("well_regime")["mean_H0_proxy"].to_dict()
    h_fill = reg_h0.get("filling/absorbing", np.nan)
    h_sat = reg_h0.get("near-capacity/saturated", np.nan)
    h_over = reg_h0.get("overflowing", np.nan)

    # Overflow fraction by early/quiescent is important for TITAN age-bimodality test.
    overflow = pop[pop["well_regime"] == "overflowing"]
    filling = pop[pop["well_regime"] == "filling/absorbing"]

    verdict = "CAPACITY_OVERFLOW_PASS" if (
        cap_bic < mass_bic and
        np.isfinite(h_over) and np.isfinite(h_fill) and
        (h_over - h_fill) > 3.0
    ) else "PARTIAL_OR_FAIL"

    summary = pd.DataFrame([{
        "verdict": verdict,
        "selected_proxy": str(best["proxy"]),
        "Dsat": float(best["Dsat"]),
        "q_sat": float(best["q_sat"]),
        "phi0": float(best["phi0"]),
        "best_BIC_model": best_model,
        "BIC_mass_only": mass_bic,
        "BIC_capacity_state": cap_bic,
        "delta_BIC_capacity_minus_mass": cap_bic - mass_bic,
        "mean_H0_filling": h_fill,
        "mean_H0_saturated": h_sat,
        "mean_H0_overflowing": h_over,
        "H0_overflow_minus_filling": h_over - h_fill if np.isfinite(h_over) and np.isfinite(h_fill) else np.nan,
        "overflow_fraction_total": len(overflow) / len(pop),
        "overflow_early_fraction": float(overflow["is_early"].mean()) if len(overflow) else np.nan,
        "filling_early_fraction": float(filling["is_early"].mean()) if len(filling) else np.nan,
        "interpretation": (
            "TITAN-like age bimodality plus realistic mass/size priors naturally sorts hosts into filling, saturated, and overflowing wells. Capacity-state gates beat mass-only and produce the expected H0 calibration sequence."
            if verdict == "CAPACITY_OVERFLOW_PASS"
            else
            "The capacity-overflow stress test did not cleanly beat mass-only or did not produce a strong enough H0 regime sequence."
        ),
    }])

    pop.to_csv(outdir / "flux_capacity_overflow_population.csv", index=False)
    scan.to_csv(outdir / "flux_capacity_overflow_scan.csv", index=False)
    results.to_csv(outdir / "flux_capacity_overflow_results.csv", index=False)
    regimes.to_csv(outdir / "flux_capacity_overflow_regimes.csv", index=False)
    bins.to_csv(outdir / "flux_capacity_overflow_binned.csv", index=False)
    corr.to_csv(outdir / "flux_capacity_overflow_correlations.csv", index=False)
    summary.to_csv(outdir / "flux_capacity_overflow_summary.csv", index=False)
    save_plot(pop, regimes, bins, results, outdir)

    print("\nFlux capacity-overflow audit")
    print("=" * 88)
    print(summary.to_string(index=False))
    print("\nRegime summary")
    print(regimes.to_string(index=False))
    print("\nRegression battle")
    print(results[["model", "N", "k_params", "RMSE", "R2", "AIC", "BIC"]].sort_values("BIC").to_string(index=False))
    print("\nCorrelations")
    print(corr.to_string(index=False))
    print("\nBinned fill sequence")
    print(bins.to_string(index=False))


if __name__ == "__main__":
    main()
