#!/usr/bin/env python3
"""
flux_saturation_calibration.py

Synthetic host-population scan for the Flux/MIP H0 saturation-calibration hypothesis.

Purpose
-------
The late-time phantom-lite DE background did not raise H0. The next Flux-native
possibility is that H0 tension is a local-to-global calibration mismatch:
standard candles sit inside mature registration wells, while the CMB ruler is
set by the early unsaturated substrate.

This script builds a synthetic host population and tests whether an environment-
dependent saturation gate can produce an ~8.3% inferred H0 offset without using
ordinary GR gravitational redshift.

Core equations
--------------
D_reg proxy:
    D_reg ∝ (M_star / R_eff)^alpha * (age / age0)^beta

Saturation gate:
    G_sat = D_reg^q / (D_reg^q + D_sat^q)

Calibration bias:
    phi_sat = phi0 * G_sat

Local inferred H0:
    H0_inferred = H0_global * (1 + phi_sat)

Distance modulus shift:
    Delta_mu = -5 log10(1 + phi_sat)

Interpretation
--------------
This is an empirical standard-candle calibration bias, not conventional
gravitational time dilation.

Outputs
-------
flux_saturation_hosts.csv
flux_saturation_summary.csv
flux_saturation_binned.csv
flux_saturation_calibration.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

H0_GLOBAL = 67.4
H0_LOCAL_TARGET = 73.0
PHI_REQUIRED = H0_LOCAL_TARGET / H0_GLOBAL - 1.0
AGE0 = 5.0
EPS = 1e-12


def rng_from_seed(seed: int):
    return np.random.default_rng(seed)


def make_synthetic_hosts(n=5000, seed=42):
    """
    Synthetic SN Ia host population.
    Produces a mix of dwarfs, spirals, and ellipticals.
    """
    rng = rng_from_seed(seed)

    # Type mixture
    types = rng.choice(["dwarf", "spiral", "elliptical"], size=n, p=[0.25, 0.50, 0.25])

    logM = np.zeros(n)
    Reff = np.zeros(n)
    age = np.zeros(n)
    sSFR = np.zeros(n)

    for i, t in enumerate(types):
        if t == "dwarf":
            logM[i] = rng.normal(8.8, 0.45)
            Reff[i] = 10 ** rng.normal(0.30, 0.18)  # kpc
            age[i] = np.clip(rng.normal(3.0, 1.5), 0.5, 9.0)
            sSFR[i] = 10 ** rng.normal(-9.4, 0.45)
        elif t == "spiral":
            logM[i] = rng.normal(10.3, 0.38)
            Reff[i] = 10 ** rng.normal(0.70, 0.18)
            age[i] = np.clip(rng.normal(6.0, 2.0), 1.0, 11.0)
            sSFR[i] = 10 ** rng.normal(-10.0, 0.45)
        else:
            logM[i] = rng.normal(11.0, 0.35)
            Reff[i] = 10 ** rng.normal(0.55, 0.20)
            age[i] = np.clip(rng.normal(9.0, 1.5), 3.0, 13.0)
            sSFR[i] = 10 ** rng.normal(-11.5, 0.35)

    # Clip to plausible synthetic range
    logM = np.clip(logM, 7.2, 12.1)
    Mstar = 10 ** logM
    Reff = np.clip(Reff, 0.4, 20.0)

    return pd.DataFrame({
        "host_type": types,
        "logMstar": logM,
        "Mstar": Mstar,
        "Reff_kpc": Reff,
        "stellar_age_Gyr": age,
        "sSFR": sSFR,
    })


def compute_registration_depth(df, alpha=1.0, beta=0.7):
    """
    Dimensionless registration-depth proxy.

    Uses Mstar/Reff as potential-depth proxy and age as persistence proxy.
    Normalized to the median spiral-like host scale.
    """
    potential_proxy = (df["Mstar"].to_numpy() / np.maximum(df["Reff_kpc"].to_numpy(), EPS))
    age_proxy = (df["stellar_age_Gyr"].to_numpy() / AGE0)

    raw = (potential_proxy ** alpha) * (age_proxy ** beta)

    # Normalize by median to produce order-unity depths.
    norm = np.nanmedian(raw)
    return raw / max(norm, EPS)


def saturation_gate(Dreg, Dsat=1.0, q=3.0):
    D = np.maximum(np.asarray(Dreg), EPS)
    return (D ** q) / (D ** q + Dsat ** q)


def apply_flux_saturation(df, phi0=PHI_REQUIRED, Dsat=1.0, q=3.0, alpha=1.0, beta=0.7, noise_mu=0.08, seed=43):
    rng = rng_from_seed(seed)
    out = df.copy()
    Dreg = compute_registration_depth(out, alpha=alpha, beta=beta)
    Gsat = saturation_gate(Dreg, Dsat=Dsat, q=q)
    phi = phi0 * Gsat
    H0_inf = H0_GLOBAL * (1.0 + phi)
    dmu = -5.0 * np.log10(1.0 + phi)

    # Synthetic observed Hubble residual with SN scatter
    residual_mu_obs = dmu + rng.normal(0.0, noise_mu, len(out))

    out["D_reg"] = Dreg
    out["G_sat"] = Gsat
    out["phi_sat"] = phi
    out["H0_inferred"] = H0_inf
    out["delta_mu_flux"] = dmu
    out["hubble_residual_mu_obs"] = residual_mu_obs
    return out


def summarize(df):
    rows = []
    for label, sub in [("all", df)] + [(t, df[df["host_type"] == t]) for t in ["dwarf", "spiral", "elliptical"]]:
        rows.append({
            "sample": label,
            "N": len(sub),
            "mean_logMstar": sub["logMstar"].mean(),
            "mean_D_reg": sub["D_reg"].mean(),
            "median_G_sat": sub["G_sat"].median(),
            "mean_phi_sat": sub["phi_sat"].mean(),
            "mean_H0_inferred": sub["H0_inferred"].mean(),
            "median_H0_inferred": sub["H0_inferred"].median(),
            "mean_delta_mu_flux": sub["delta_mu_flux"].mean(),
            "std_delta_mu_flux": sub["delta_mu_flux"].std(),
        })

    # Mass-step style split
    high = df[df["logMstar"] >= 10.0]
    low = df[df["logMstar"] < 10.0]
    mass_step_mu = high["delta_mu_flux"].mean() - low["delta_mu_flux"].mean()
    mass_step_H0 = high["H0_inferred"].mean() - low["H0_inferred"].mean()

    extra = {
        "sample": "mass_step_high_minus_low",
        "N": len(df),
        "mean_logMstar": np.nan,
        "mean_D_reg": np.nan,
        "median_G_sat": np.nan,
        "mean_phi_sat": np.nan,
        "mean_H0_inferred": mass_step_H0,
        "median_H0_inferred": np.nan,
        "mean_delta_mu_flux": mass_step_mu,
        "std_delta_mu_flux": np.nan,
    }
    rows.append(extra)

    return pd.DataFrame(rows)


def binned_stats(df):
    bins = [7.0, 8.5, 9.5, 10.0, 10.5, 11.0, 12.3]
    labels = ["7-8.5", "8.5-9.5", "9.5-10", "10-10.5", "10.5-11", "11+"]
    tmp = df.copy()
    tmp["mass_bin"] = pd.cut(tmp["logMstar"], bins=bins, labels=labels, include_lowest=True)
    b = tmp.groupby("mass_bin", observed=False).agg(
        N=("logMstar", "size"),
        mean_logMstar=("logMstar", "mean"),
        mean_D_reg=("D_reg", "mean"),
        mean_G_sat=("G_sat", "mean"),
        mean_phi_sat=("phi_sat", "mean"),
        mean_H0_inferred=("H0_inferred", "mean"),
        mean_delta_mu_flux=("delta_mu_flux", "mean"),
        mean_residual_obs=("hubble_residual_mu_obs", "mean"),
    ).reset_index()
    return b


def scan_params(n=4000, seed=42):
    """
    Scan saturation parameters for a plausible bridge.
    Need:
      - all-host mean H0 can approach local target if calibrators are high-depth.
      - low-depth hosts stay closer to global.
      - mass step in magnitude is not absurdly huge.
    """
    base = make_synthetic_hosts(n=n, seed=seed)
    rows = []
    for phi0 in np.linspace(0.04, 0.12, 17):
        for Dsat in np.linspace(0.6, 2.2, 17):
            for q in [1.5, 2.0, 3.0, 4.0, 6.0]:
                df = apply_flux_saturation(base, phi0=phi0, Dsat=Dsat, q=q, noise_mu=0.08, seed=seed+1)
                high = df[df["logMstar"] >= 10.0]
                low = df[df["logMstar"] < 10.0]
                ell = df[df["host_type"] == "elliptical"]
                dwarf = df[df["host_type"] == "dwarf"]

                mass_step_mu = high["delta_mu_flux"].mean() - low["delta_mu_flux"].mean()
                mass_step_abs = abs(mass_step_mu)
                h0_high = high["H0_inferred"].mean()
                h0_low = low["H0_inferred"].mean()
                h0_ell = ell["H0_inferred"].mean()
                h0_dwarf = dwarf["H0_inferred"].mean()

                # Score wants saturated massive hosts near 73, dwarfs closer to 67,
                # and mass-step not insanely larger than known SN environmental effects.
                score = (
                    abs(h0_high - H0_LOCAL_TARGET)
                    + 0.4 * abs(h0_low - H0_GLOBAL)
                    + 10.0 * max(mass_step_abs - 0.18, 0.0) # avoid huge mag step
                )

                rows.append({
                    "phi0": phi0,
                    "Dsat": Dsat,
                    "q": q,
                    "h0_high_mass": h0_high,
                    "h0_low_mass": h0_low,
                    "h0_elliptical": h0_ell,
                    "h0_dwarf": h0_dwarf,
                    "h0_split_high_low": h0_high - h0_low,
                    "mass_step_mu_high_minus_low": mass_step_mu,
                    "mass_step_abs": mass_step_abs,
                    "score": score,
                })
    return pd.DataFrame(rows).sort_values("score").reset_index(drop=True)


def save_plot(df, binned, summary, outdir):
    plt.figure(figsize=(13, 10))

    plt.subplot(2, 2, 1)
    colors = {"dwarf": "tab:blue", "spiral": "tab:green", "elliptical": "tab:red"}
    for t, sub in df.groupby("host_type"):
        plt.scatter(sub["logMstar"], sub["H0_inferred"], s=8, alpha=0.25, label=t)
    plt.axhline(H0_GLOBAL, linestyle="--", color="k", label="global/CMB")
    plt.axhline(H0_LOCAL_TARGET, linestyle=":", color="k", label="local target")
    plt.xlabel("log10 stellar mass")
    plt.ylabel("inferred H0")
    plt.title("Host mass vs inferred H0")
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.scatter(df["D_reg"], df["G_sat"], s=8, alpha=0.25)
    plt.xscale("log")
    plt.xlabel("registration depth proxy D_reg")
    plt.ylabel("G_sat")
    plt.title("Saturation gate")

    plt.subplot(2, 2, 3)
    plt.errorbar(
        np.arange(len(binned)),
        binned["mean_H0_inferred"],
        yerr=None,
        marker="o",
    )
    plt.xticks(np.arange(len(binned)), binned["mass_bin"], rotation=25)
    plt.axhline(H0_GLOBAL, linestyle="--", color="k")
    plt.axhline(H0_LOCAL_TARGET, linestyle=":", color="k")
    plt.xlabel("host stellar-mass bin")
    plt.ylabel("mean inferred H0")
    plt.title("Binned H0 calibration bias")

    plt.subplot(2, 2, 4)
    plt.errorbar(
        np.arange(len(binned)),
        binned["mean_delta_mu_flux"],
        marker="o",
    )
    plt.xticks(np.arange(len(binned)), binned["mass_bin"], rotation=25)
    plt.axhline(0, linestyle="--", color="k")
    plt.xlabel("host stellar-mass bin")
    plt.ylabel("mean Delta mu")
    plt.title("Flux mass-step-like magnitude shift")

    plt.tight_layout()
    plt.savefig(outdir / "flux_saturation_calibration.png", dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--n-hosts", type=int, default=7000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    scan = scan_params(n=3000, seed=args.seed)
    best = scan.iloc[0]

    hosts0 = make_synthetic_hosts(n=args.n_hosts, seed=args.seed)
    hosts = apply_flux_saturation(
        hosts0,
        phi0=float(best["phi0"]),
        Dsat=float(best["Dsat"]),
        q=float(best["q"]),
        noise_mu=0.08,
        seed=args.seed + 10,
    )
    summ = summarize(hosts)
    binned = binned_stats(hosts)

    high_h0 = float(summ.loc[summ["sample"] == "mass_step_high_minus_low", "mean_H0_inferred"].iloc[0])
    high_mass_mean = float(hosts[hosts["logMstar"] >= 10]["H0_inferred"].mean())
    low_mass_mean = float(hosts[hosts["logMstar"] < 10]["H0_inferred"].mean())

    verdict = "SATURATION_BRIDGE_TOY_PASS" if high_mass_mean > 71.5 and low_mass_mean < 70.5 else "PARTIAL_OR_FAIL"

    summary_meta = pd.DataFrame([{
        "verdict": verdict,
        "H0_global": H0_GLOBAL,
        "H0_local_target": H0_LOCAL_TARGET,
        "phi_required": PHI_REQUIRED,
        "best_phi0": float(best["phi0"]),
        "best_Dsat": float(best["Dsat"]),
        "best_q": float(best["q"]),
        "mean_H0_high_mass": high_mass_mean,
        "mean_H0_low_mass": low_mass_mean,
        "H0_split_high_minus_low": high_mass_mean - low_mass_mean,
        "mass_step_mu_high_minus_low": float(best["mass_step_mu_high_minus_low"]),
        "interpretation": (
            "Synthetic saturated hosts can produce a local-calibration H0 elevation while low-depth hosts remain closer to the global value. This is a testable host-environment bias, not a conventional gravitational redshift."
            if verdict == "SATURATION_BRIDGE_TOY_PASS"
            else
            "The synthetic saturation bridge did not cleanly separate high-depth and low-depth hosts enough in this scan."
        ),
    }])

    hosts.to_csv(outdir / "flux_saturation_hosts.csv", index=False)
    summ.to_csv(outdir / "flux_saturation_population_summary.csv", index=False)
    binned.to_csv(outdir / "flux_saturation_binned.csv", index=False)
    scan.to_csv(outdir / "flux_saturation_scan_results.csv", index=False)
    summary_meta.to_csv(outdir / "flux_saturation_summary.csv", index=False)
    save_plot(hosts, binned, summ, outdir)

    print("\nFlux saturation calibration scan")
    print("=" * 80)
    print(summary_meta.to_string(index=False))
    print("\nPopulation summary")
    print(summ.to_string(index=False))
    print("\nBinned summary")
    print(binned.to_string(index=False))


if __name__ == "__main__":
    main()
