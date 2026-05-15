#!/usr/bin/env python3
"""
mip_sound_horizon_bridge.py

Fast toy MIP / flux-coherence sound-horizon bridge for the H0 tension.

Question:
Can an early Minimum Interface Point / flux-coherence modulation shrink the
sound horizon r_s enough to map H0_CMB ~ 67.4 toward H0_local ~ 73?

This is a controlled r_s sensitivity test, not a full Boltzmann solver.

Outputs:
    mip_sound_horizon_scan_results.csv
    mip_sound_horizon_summary.csv
    mip_sound_horizon_best_curve.csv
    mip_sound_horizon_bridge.png
"""

from __future__ import annotations
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

H0_CMB = 67.4
H0_LOCAL = 73.0
TARGET_RATIO = H0_CMB / H0_LOCAL

OMEGA_M0 = 0.315
OMEGA_B0 = 0.049
OMEGA_GAMMA0 = 5.38e-5
OMEGA_R0 = 9.2e-5
OMEGA_L0 = 1.0 - OMEGA_M0 - OMEGA_R0
Z_STAR = 1090.0
EPS = 1e-12


def E_lcdm_z(z):
    z = np.asarray(z)
    return np.sqrt(OMEGA_R0 * (1+z)**4 + OMEGA_M0 * (1+z)**3 + OMEGA_L0)


def cs_lcdm_z(z):
    z = np.asarray(z)
    a = 1.0 / (1.0 + z)
    R = (3.0 * OMEGA_B0 * a**-3) / (4.0 * OMEGA_GAMMA0 * a**-4)
    return 1.0 / np.sqrt(3.0 * (1.0 + R))


def mip_profile(z, sigma_i, width_ln):
    x = np.log((1.0 + z) / (1.0 + Z_STAR))
    return np.clip(sigma_i * np.exp(-0.5 * (x / width_ln)**2), 0.0, 1.0)


def build_grid(zmax=1e7, n=1600):
    z = np.geomspace(Z_STAR, zmax, n)
    z[0] = Z_STAR
    E0 = E_lcdm_z(z)
    cs0 = cs_lcdm_z(z)
    base = cs0 / np.maximum(E0, EPS)
    rs0 = np.trapz(base, z)
    return z, E0, cs0, base, rs0


def run_scan():
    z, E0, cs0, base, rs0 = build_grid()

    sigma_vals = np.linspace(0.0, 1.0, 26)
    AH_vals = np.linspace(0.0, 0.25, 26)
    Acs_vals = np.linspace(0.0, 0.12, 21)
    width_vals = [0.45, 0.65, 0.85]

    rows = []
    for width in width_vals:
        profile_unit = mip_profile(z, 1.0, width)
        for sigma in sigma_vals:
            C = sigma * profile_unit
            Cstar = float(sigma)
            for AH in AH_vals:
                Hfac = np.sqrt(1.0 + AH * C)
                ede_like = (AH * Cstar) / (1.0 + AH * Cstar) if AH > 0 else 0.0
                for Acs in Acs_vals:
                    cs_factor = np.clip(1.0 - Acs * C, 0.05, 1.0)
                    rs = np.trapz(base * cs_factor / Hfac, z)
                    ratio = rs / rs0
                    h0_eff = H0_CMB / max(ratio, EPS)
                    cs_ratio_star = max(1.0 - Acs * Cstar, 0.05)
                    cmb_mild = (ede_like < 0.10) and (cs_ratio_star > 0.88)
                    rows.append({
                        "sigma_i": sigma,
                        "A_H": AH,
                        "A_cs": Acs,
                        "width_ln": width,
                        "rs_dimless": rs,
                        "rs_ratio": ratio,
                        "target_ratio": TARGET_RATIO,
                        "H0_eff": h0_eff,
                        "delta_H0": h0_eff - H0_CMB,
                        "target_distance": abs(ratio - TARGET_RATIO),
                        "C_star": Cstar,
                        "Hboost_star": float(np.sqrt(1.0 + AH * Cstar)),
                        "cs_ratio_star": cs_ratio_star,
                        "ede_like_frac_star": ede_like,
                        "cmb_mild": cmb_mild,
                    })

    df = pd.DataFrame(rows)
    df["score"] = df["target_distance"] + np.where(df["cmb_mild"], 0.0, 10.0)
    df = df.sort_values("score").reset_index(drop=True)
    return df, rs0


def make_curve(best):
    z, E0, cs0, base, rs0 = build_grid(n=3000)
    C = mip_profile(z, float(best["sigma_i"]), float(best["width_ln"]))
    Hfac = np.sqrt(1.0 + float(best["A_H"]) * C)
    csfac = np.clip(1.0 - float(best["A_cs"]) * C, 0.05, 1.0)
    return pd.DataFrame({
        "z": z,
        "C_mip": C,
        "E_LCDM": E0,
        "E_Flux": E0 * Hfac,
        "H_boost": Hfac,
        "cs_LCDM": cs0,
        "cs_Flux": cs0 * csfac,
        "cs_ratio": csfac,
        "rs_integrand_LCDM": base,
        "rs_integrand_Flux": base * csfac / Hfac,
    })


def save_plot(best, curve, outdir):
    plt.figure(figsize=(12, 9))

    plt.subplot(2, 2, 1)
    plt.plot(curve["z"], curve["C_mip"])
    plt.xscale("log")
    plt.xlabel("redshift z")
    plt.ylabel("C_MIP")
    plt.title("MIP coherence profile")

    plt.subplot(2, 2, 2)
    plt.plot(curve["z"], curve["H_boost"], label="H boost")
    plt.plot(curve["z"], curve["cs_ratio"], label="sound-speed ratio")
    plt.xscale("log")
    plt.xlabel("redshift z")
    plt.ylabel("ratio")
    plt.title("Early ruler modification")
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(curve["z"], curve["rs_integrand_LCDM"], label="LCDM")
    plt.plot(curve["z"], curve["rs_integrand_Flux"], label="Flux/MIP")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("redshift z")
    plt.ylabel("c_s / E(z)")
    plt.title("Sound-horizon integrand")
    plt.legend()

    plt.subplot(2, 2, 4)
    plt.bar(["H0 CMB", "H0 via r_s bridge", "H0 local"], [H0_CMB, float(best["H0_eff"]), H0_LOCAL])
    plt.ylabel("H0 [km/s/Mpc]")
    plt.title("Effective H0 from sound-horizon shrinkage")

    plt.tight_layout()
    plt.savefig(outdir / "mip_sound_horizon_bridge.png", dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df, rs0 = run_scan()
    best = df.iloc[0]
    curve = make_curve(best)

    verdict = "RS_BRIDGE_WINDOW" if abs(float(best["rs_ratio"]) - TARGET_RATIO) < 0.003 and bool(best["cmb_mild"]) else "NO_RS_BRIDGE"

    summary = pd.DataFrame([{
        "verdict": verdict,
        "rs_lcdm_dimless": rs0,
        "target_rs_ratio": TARGET_RATIO,
        "best_rs_ratio": float(best["rs_ratio"]),
        "H0_CMB": H0_CMB,
        "H0_local_target": H0_LOCAL,
        "best_H0_eff": float(best["H0_eff"]),
        "best_sigma_i": float(best["sigma_i"]),
        "best_A_H": float(best["A_H"]),
        "best_A_cs": float(best["A_cs"]),
        "best_width_ln": float(best["width_ln"]),
        "best_C_star": float(best["C_star"]),
        "best_Hboost_star": float(best["Hboost_star"]),
        "best_cs_ratio_star": float(best["cs_ratio_star"]),
        "best_ede_like_frac_star": float(best["ede_like_frac_star"]),
        "interpretation": (
            "A mild recombination-era MIP/coherence modulation can shrink the sound horizon enough to map CMB H0 toward local H0 in this toy sensitivity model. Full Boltzmann/CMB peak testing remains required."
            if verdict == "RS_BRIDGE_WINDOW"
            else
            "The scanned mild MIP/coherence modulation did not hit the required sound-horizon shrinkage window."
        ),
    }])

    df.to_csv(outdir / "mip_sound_horizon_scan_results.csv", index=False)
    summary.to_csv(outdir / "mip_sound_horizon_summary.csv", index=False)
    curve.to_csv(outdir / "mip_sound_horizon_best_curve.csv", index=False)
    save_plot(best, curve, outdir)

    print("\nMIP sound-horizon bridge scan")
    print("=" * 72)
    print(summary.to_string(index=False))
    print("\nTop 10")
    cols = ["sigma_i", "A_H", "A_cs", "width_ln", "rs_ratio", "H0_eff", "Hboost_star", "cs_ratio_star", "ede_like_frac_star", "cmb_mild"]
    print(df[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
