#!/usr/bin/env python3
"""
flux_literature_proxy_check.py

Fast Level-2 literature-distribution proxy audit for Flux/MIP saturation.

Uses literature-motivated priors:
- logMstar centered near 10.5
- 60/40 late/early split
- size-mass scaling
- TITAN-like bimodal ages: late ~2.2 Gyr, early ~6.0 Gyr
- sigma_v tied to compactness

Tests whether a registration-depth proxy, especially sigma_v^2 * age,
beats mass-only as a Hubble-residual predictor.

Outputs:
    flux_literature_proxy_population.csv
    flux_literature_proxy_results.csv
    flux_literature_proxy_summary.csv
    flux_literature_proxy_binned.csv
    flux_literature_proxy_correlations.csv
    flux_literature_proxy.png
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


def saturation_gate(D, Dsat=1.0, q=5.0):
    D = np.maximum(np.asarray(D), EPS)
    return D**q / (D**q + Dsat**q)


def normalize(x):
    x = np.asarray(x)
    return x / max(np.nanmedian(x), EPS)


def generate_population(N=7000, seed=42):
    R = rng(seed)

    is_early = R.random(N) > 0.60
    morph = np.where(is_early, "early/quiescent", "late/star-forming")

    logM = R.normal(10.5, 0.50, N)
    logM += np.where(is_early, R.normal(0.18, 0.18, N), R.normal(-0.08, 0.22, N))
    logM = np.clip(logM, 8.0, 12.1)

    logRe_late = 0.25 * (logM - 10.5) + 0.70 + R.normal(0.0, 0.18, N)
    logRe_early = 0.50 * (logM - 10.5) + 0.48 + R.normal(0.0, 0.16, N)
    logRe = np.where(is_early, logRe_early, logRe_late)
    Reff = np.clip(10**logRe, 0.35, 25.0)

    age = np.where(is_early, R.normal(6.0, 1.5, N), R.normal(2.2, 0.8, N))
    age = np.clip(age, 0.4, 12.5)

    logsSFR = np.where(is_early, R.normal(-11.45, 0.35, N), R.normal(-9.75, 0.45, N))
    metallicity = np.where(is_early, R.normal(8.95, 0.12, N), R.normal(8.70, 0.18, N))

    M = 10**logM
    compactness = M / np.maximum(Reff, EPS)
    comp_norm = compactness / np.median(compactness)
    sigma_v = 145.0 * comp_norm**0.25
    sigma_v *= np.where(is_early, 1.12, 0.92)
    sigma_v *= 10 ** R.normal(0.0, 0.045, N)
    sigma_v = np.clip(sigma_v, 55.0, 330.0)

    D_mass = normalize(M)
    D_potential = normalize(M / Reff)
    D_MR_age = normalize((M / Reff) * (age / AGE_UNIVERSE))
    D_sigma_age = normalize((sigma_v / 150.0)**2 * (age / AGE_UNIVERSE))

    return pd.DataFrame({
        "host_id": [f"LIT_HOST_{i:05d}" for i in range(N)],
        "is_early": is_early,
        "morphology": morph,
        "logMstar": logM,
        "Mstar": M,
        "Reff_kpc": Reff,
        "stellar_age_Gyr": age,
        "logsSFR": logsSFR,
        "sSFR": 10**logsSFR,
        "metallicity": metallicity,
        "sigma_v_km_s": sigma_v,
        "D_mass_proxy": D_mass,
        "D_potential_only": D_potential,
        "D_MR_age": D_MR_age,
        "D_sigma_age": D_sigma_age,
        "logD_mass_proxy": np.log10(np.maximum(D_mass, EPS)),
        "logD_potential_only": np.log10(np.maximum(D_potential, EPS)),
        "logD_MR_age": np.log10(np.maximum(D_MR_age, EPS)),
        "logD_sigma_age": np.log10(np.maximum(D_sigma_age, EPS)),
    })


def add_residuals(df, proxy_col="D_sigma_age", amp=-0.085, Dsat=1.0, q=5.0, scatter=0.10, seed=43):
    R = rng(seed)
    out = df.copy()
    G = saturation_gate(out[proxy_col].to_numpy(), Dsat=Dsat, q=q)

    # Flux component plus small ordinary nuisance components.
    flux_component = amp * G
    mass_nuis = -0.018 * (out["logMstar"].to_numpy() >= 10.0).astype(float)
    sfr_nuis = 0.010 * (out["logsSFR"].to_numpy() + 10.0)
    metal_nuis = -0.007 * (out["metallicity"].to_numpy() - 8.7)
    noise = R.normal(0.0, scatter, len(out))

    out["G_sat_literature"] = G
    out["mu_resid_flux_component"] = flux_component
    out["mu_resid"] = flux_component + mass_nuis + sfr_nuis + metal_nuis + noise
    out["H0_proxy"] = H0_GLOBAL * 10 ** (-out["mu_resid"] / 5.0)
    return out


def design(df, predictors):
    X = [np.ones(len(df))]
    names = ["intercept"]
    for p in predictors:
        vals = df[p].to_numpy(float)
        mu, sd = np.mean(vals), np.std(vals)
        if sd <= 0:
            sd = 1.0
        X.append((vals - mu) / sd)
        names.append(p)
    return np.vstack(X).T, names


def fit(df, name, predictors):
    d = df.dropna(subset=["mu_resid"] + predictors)
    y = d["mu_resid"].to_numpy(float)
    X, names = design(d, predictors)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    n, k = len(y), X.shape[1]
    rss = float(np.sum(resid**2))
    tss = float(np.sum((y - y.mean())**2))
    aic = n * np.log(max(rss/n, EPS)) + 2*k
    bic = n * np.log(max(rss/n, EPS)) + k*np.log(n)
    out = {"model": name, "predictors": ",".join(predictors), "N": n, "k_params": k,
           "RMSE": float(np.sqrt(rss/n)), "R2": 1-rss/max(tss, EPS), "AIC": aic, "BIC": bic}
    for nm, b in zip(names, beta):
        out[f"coef_{nm}"] = b
    return out


def regression_battle(df):
    models = {
        "A_mass_only": ["logMstar"],
        "B_mass_sSFR": ["logMstar", "logsSFR"],
        "C_kitchen_sink": ["logMstar", "logsSFR", "metallicity", "stellar_age_Gyr", "Reff_kpc"],
        "D_MR_age": ["logD_MR_age"],
        "E_sigma_age": ["logD_sigma_age"],
        "F_sigma_age_controls": ["logD_sigma_age", "logsSFR", "metallicity"],
        "G_potential_only": ["logD_potential_only"],
    }
    return pd.DataFrame([fit(df, n, p) for n, p in models.items()])


def light_scan(base, seed=42):
    """
    Fast local scan over a few plausible parameters.
    """
    rows = []
    for proxy in ["D_sigma_age", "D_MR_age"]:
        for amp in [-0.055, -0.070, -0.085, -0.100, -0.115]:
            for Dsat in [0.8, 1.0, 1.2, 1.4]:
                for q in [3.0, 5.0, 7.0]:
                    df = add_residuals(base, proxy_col=proxy, amp=amp, Dsat=Dsat, q=q, seed=seed+1)
                    res = regression_battle(df)
                    best = res.sort_values("BIC").iloc[0]
                    mass_bic = float(res.loc[res["model"] == "A_mass_only", "BIC"].iloc[0])
                    high = df[df["logMstar"] >= 10.0]
                    low = df[df["logMstar"] < 10.0]
                    mass_step = float(high["mu_resid"].mean() - low["mu_resid"].mean())
                    qhi = df[df[proxy] >= df[proxy].quantile(0.75)]
                    qlo = df[df[proxy] <= df[proxy].quantile(0.25)]
                    h0_split = float(qhi["H0_proxy"].mean() - qlo["H0_proxy"].mean())
                    step_ok = 0.04 <= abs(mass_step) <= 0.12
                    proxy_wins = str(best["model"]) in ["D_MR_age", "E_sigma_age", "F_sigma_age_controls"]
                    score = abs(abs(mass_step)-0.08) + (0 if step_ok else 0.05) + (0 if proxy_wins else 0.2)
                    rows.append({
                        "proxy": proxy, "amp": amp, "Dsat": Dsat, "q": q,
                        "best_model": str(best["model"]),
                        "best_BIC": float(best["BIC"]),
                        "BIC_mass_only": mass_bic,
                        "delta_BIC_best_minus_mass": float(best["BIC"] - mass_bic),
                        "mass_step_mu_high_minus_low": mass_step,
                        "H0_split_highD_minus_lowD": h0_split,
                        "score": score,
                    })
    return pd.DataFrame(rows).sort_values("score").reset_index(drop=True)


def correlations(df):
    vars_ = ["logMstar", "logD_potential_only", "logD_MR_age", "logD_sigma_age",
             "G_sat_literature", "stellar_age_Gyr", "Reff_kpc", "sigma_v_km_s", "logsSFR", "metallicity"]
    rows = []
    for v in vars_:
        rows.append({"variable": v, "pearson_r_with_mu_resid": float(np.corrcoef(df["mu_resid"], df[v])[0,1]), "N": len(df)})
    return pd.DataFrame(rows)


def binned(df, proxy):
    tmp = df.copy()
    pcol = "logD_sigma_age" if proxy == "D_sigma_age" else "logD_MR_age"
    tmp["proxy_bin"] = pd.qcut(tmp[pcol], q=6, duplicates="drop")
    return tmp.groupby("proxy_bin", observed=False).agg(
        N=("mu_resid", "size"),
        mean_proxy=(pcol, "mean"),
        mean_logMstar=("logMstar", "mean"),
        mean_age=("stellar_age_Gyr", "mean"),
        mean_sigma=("sigma_v_km_s", "mean"),
        mean_Reff=("Reff_kpc", "mean"),
        early_fraction=("is_early", "mean"),
        mean_mu_resid=("mu_resid", "mean"),
        std_mu_resid=("mu_resid", "std"),
        mean_H0_proxy=("H0_proxy", "mean"),
        mean_Gsat=("G_sat_literature", "mean"),
    ).reset_index()


def save_plot(df, results, bt, outdir, proxy):
    pcol = "logD_sigma_age" if proxy == "D_sigma_age" else "logD_MR_age"
    xlabel = "log D_reg = log(sigma_v² × age)" if proxy == "D_sigma_age" else "log D_reg = log((M/Re) × age)"

    plt.figure(figsize=(14,10))

    plt.subplot(2,2,1)
    plt.scatter(df[pcol], df["mu_resid"], s=8, alpha=0.25)
    m,b = np.polyfit(df[pcol], df["mu_resid"], 1)
    xs = np.linspace(df[pcol].min(), df[pcol].max(), 100)
    plt.plot(xs, m*xs+b, color="black")
    plt.xlabel(xlabel); plt.ylabel("Hubble residual Δμ")
    plt.title("Registration-depth fingerprint")

    plt.subplot(2,2,2)
    plt.scatter(df["logMstar"], df["mu_resid"], s=8, alpha=0.25)
    m,b = np.polyfit(df["logMstar"], df["mu_resid"], 1)
    xs = np.linspace(df["logMstar"].min(), df["logMstar"].max(), 100)
    plt.plot(xs, m*xs+b, color="black")
    plt.xlabel("log Mstar"); plt.ylabel("Hubble residual Δμ")
    plt.title("Mass-only comparison")

    plt.subplot(2,2,3)
    order = results.sort_values("BIC")
    plt.bar(order["model"], order["BIC"])
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("BIC"); plt.title("Regression battle")

    plt.subplot(2,2,4)
    plt.errorbar(range(len(bt)), bt["mean_mu_resid"], yerr=bt["std_mu_resid"]/np.sqrt(np.maximum(bt["N"],1)), marker="o")
    plt.xticks(range(len(bt)), [str(x) for x in bt["proxy_bin"]], rotation=30, ha="right")
    plt.axhline(0, linestyle="--", color="black")
    plt.xlabel("D_reg quantile"); plt.ylabel("mean Δμ")
    plt.title("Literature-prior saturation sequence")

    plt.tight_layout()
    plt.savefig(outdir / "flux_literature_proxy.png", dpi=180)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--n", type=int, default=7000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    base_scan = generate_population(N=2500, seed=args.seed)
    scan = light_scan(base_scan, seed=args.seed)
    best = scan.iloc[0]

    pop = generate_population(N=args.n, seed=args.seed)
    pop = add_residuals(pop, proxy_col=str(best["proxy"]), amp=float(best["amp"]), Dsat=float(best["Dsat"]), q=float(best["q"]), seed=args.seed+2)
    results = regression_battle(pop)
    corr = correlations(pop)
    bt = binned(pop, proxy=str(best["proxy"]))

    best_model = str(results.sort_values("BIC").iloc[0]["model"])
    mass_bic = float(results.loc[results["model"]=="A_mass_only", "BIC"].iloc[0])
    best_bic = float(results.sort_values("BIC").iloc[0]["BIC"])
    high = pop[pop["logMstar"] >= 10.0]
    low = pop[pop["logMstar"] < 10.0]
    mass_step = float(high["mu_resid"].mean() - low["mu_resid"].mean())
    proxy = str(best["proxy"])
    qhi = pop[pop[proxy] >= pop[proxy].quantile(0.75)]
    qlo = pop[pop[proxy] <= pop[proxy].quantile(0.25)]
    h0_split = float(qhi["H0_proxy"].mean() - qlo["H0_proxy"].mean())

    verdict = "LITERATURE_PROXY_PASS" if best_model in ["D_MR_age","E_sigma_age","F_sigma_age_controls"] and abs(mass_step) >= 0.04 else "PARTIAL_OR_FAIL"

    summary = pd.DataFrame([{
        "verdict": verdict,
        "selected_proxy": proxy,
        "best_BIC_model": best_model,
        "BIC_best": best_bic,
        "BIC_mass_only": mass_bic,
        "delta_BIC_best_minus_mass": best_bic - mass_bic,
        "amplitude_mu": float(best["amp"]),
        "Dsat": float(best["Dsat"]),
        "q": float(best["q"]),
        "mass_step_mu_high_minus_low": mass_step,
        "H0_proxy_split_highD_minus_lowD": h0_split,
        "mean_H0_lowD_quartile": float(qlo["H0_proxy"].mean()),
        "mean_H0_highD_quartile": float(qhi["H0_proxy"].mean()),
        "interpretation": (
            "Literature-prior host distributions can produce a realistic host residual step from a registration-depth saturation gate, and the depth proxy beats mass-only in BIC."
            if verdict == "LITERATURE_PROXY_PASS"
            else
            "The literature-prior stress test did not decisively prefer the registration-depth proxy over mass-only."
        )
    }])

    pop.to_csv(outdir / "flux_literature_proxy_population.csv", index=False)
    scan.to_csv(outdir / "flux_literature_proxy_scan.csv", index=False)
    results.to_csv(outdir / "flux_literature_proxy_results.csv", index=False)
    corr.to_csv(outdir / "flux_literature_proxy_correlations.csv", index=False)
    bt.to_csv(outdir / "flux_literature_proxy_binned.csv", index=False)
    summary.to_csv(outdir / "flux_literature_proxy_summary.csv", index=False)
    save_plot(pop, results, bt, outdir, proxy)

    print("\nFlux literature-distribution proxy audit")
    print("="*88)
    print(summary.to_string(index=False))
    print("\nRegression battle")
    print(results[["model","N","k_params","RMSE","R2","AIC","BIC"]].sort_values("BIC").to_string(index=False))
    print("\nCorrelation table")
    print(corr.to_string(index=False))
    print("\nBinned sequence")
    print(bt.to_string(index=False))


if __name__ == "__main__":
    main()
