#!/usr/bin/env python3
"""
flux_h0_sensitivity.py

Fast toy H0 sensitivity diagnostic for the Flux gated-ledger / phantom-lite DE sector.

This is a background-only "Hubble whisper" check, not a full MCMC fit.

It estimates how much a late-time Flux phantom-lite/fossil-hump background could
bias a low-z LCDM distance fit toward higher H0 while keeping the CMB-era DE
fraction negligible.

Outputs:
    flux_h0_sensitivity_results.csv
    flux_h0_sensitivity_summary.csv
    flux_h0_sensitivity_best_curve.csv
    flux_h0_sensitivity.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OMEGA_M0 = 0.30
OMEGA_DE0 = 0.70
H0_CMB_REF = 67.4
H0_LOCAL_TARGET = 73.0
Z_CMB = 1090.0
EPS = 1e-12


def cumulative_trapezoid_np(y, x):
    out = np.zeros_like(y, dtype=float)
    dx = np.diff(x)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dx)
    return out


def E_lcdm(z):
    return np.sqrt(OMEGA_M0 * (1.0 + z) ** 3 + OMEGA_DE0)


def w_flux_z(z, w0_offset=-0.014, hump_amp=-0.020, hump_z=2.0, hump_width=0.9):
    z = np.asarray(z)
    local = w0_offset * np.exp(-(z / 0.8) ** 2)
    hump = hump_amp * np.exp(-0.5 * ((z - hump_z) / hump_width) ** 2)
    return -1.0 + local + hump


def build_model_arrays(w0_offset, hump_amp, zmax=1090.0, n_low=5000, n_high=4000):
    """
    Construct a dense low-z grid plus log-spaced high-z tail for CMB fraction.
    """
    z_low = np.linspace(0.0, 5.0, n_low)
    z_high = np.geomspace(5.0001, zmax, n_high)
    z = np.unique(np.concatenate([z_low, z_high]))

    w = w_flux_z(z, w0_offset=w0_offset, hump_amp=hump_amp)
    dlnrho_dz = 3.0 * (1.0 + w) / np.maximum(1.0 + z, EPS)
    ln_rho = cumulative_trapezoid_np(dlnrho_dz, z)
    rho_rel = np.exp(ln_rho)

    E_flux = np.sqrt(OMEGA_M0 * (1.0 + z) ** 3 + OMEGA_DE0 * rho_rel)
    E_l = E_lcdm(z)

    invE_flux_int = cumulative_trapezoid_np(1.0 / np.maximum(E_flux, EPS), z)
    invE_lcdm_int = cumulative_trapezoid_np(1.0 / np.maximum(E_l, EPS), z)

    omega_de = OMEGA_DE0 * rho_rel / np.maximum(E_flux ** 2, EPS)

    return {
        "z": z,
        "w": w,
        "rho_rel": rho_rel,
        "E_flux": E_flux,
        "E_lcdm": E_l,
        "I_flux": invE_flux_int,
        "I_lcdm": invE_lcdm_int,
        "omega_de": omega_de,
    }


def infer_h0_from_arrays(arr, zfit):
    I_flux = np.interp(zfit, arr["z"], arr["I_flux"])
    I_lcdm = np.interp(zfit, arr["z"], arr["I_lcdm"])
    return H0_CMB_REF * I_lcdm / max(I_flux, EPS)


def scan_parameter_space():
    rows = []
    # Wider grid but fast enough. Keep w0 > -1.12 "reasonable".
    w0_offsets = np.linspace(-0.005, -0.110, 22)
    hump_amps = np.linspace(0.0, -0.180, 37)
    zfits = np.array([0.1, 0.3, 0.5, 1.0, 1.5])

    for w0o in w0_offsets:
        for ha in hump_amps:
            arr = build_model_arrays(w0o, ha)
            h0_eff = np.array([infer_h0_from_arrays(arr, zf) for zf in zfits])
            w0 = np.interp(0.0, arr["z"], arr["w"])
            w_z2 = np.interp(2.0, arr["z"], arr["w"])
            omega_de_1090 = np.interp(Z_CMB, arr["z"], arr["omega_de"])

            rows.append({
                "w0_offset": w0o,
                "hump_amp": ha,
                "w0": w0,
                "w_z2": w_z2,
                "H0_eff_z0p1": h0_eff[0],
                "H0_eff_z0p3": h0_eff[1],
                "H0_eff_z0p5": h0_eff[2],
                "H0_eff_z1p0": h0_eff[3],
                "H0_eff_z1p5": h0_eff[4],
                "H0_mean_lowz": float(np.mean(h0_eff[:3])),
                "H0_mean_z0p1_to_1p5": float(np.mean(h0_eff)),
                "delta_H0_lowz": float(np.mean(h0_eff[:3]) - H0_CMB_REF),
                "omega_de_1090": omega_de_1090,
            })

    df = pd.DataFrame(rows)
    df["target_distance"] = np.abs(df["H0_mean_lowz"] - H0_LOCAL_TARGET)
    df["cmb_safe"] = df["omega_de_1090"] < 1e-8
    df["w_reasonable"] = df["w0"] > -1.12
    df["not_too_extreme_z2"] = df["w_z2"] > -1.25
    df["score"] = (
        df["target_distance"]
        + np.where(df["cmb_safe"], 0.0, 100.0)
        + np.where(df["w_reasonable"], 0.0, 20.0)
        + np.where(df["not_too_extreme_z2"], 0.0, 10.0)
    )
    return df.sort_values("score").reset_index(drop=True)


def make_best_curve(best):
    arr = build_model_arrays(float(best["w0_offset"]), float(best["hump_amp"]))
    z_plot = np.linspace(0.0, 3.0, 800)
    zfit = np.linspace(0.05, 2.0, 200)
    h0_curve = np.array([infer_h0_from_arrays(arr, zf) for zf in zfit])

    return pd.DataFrame({
        "z": z_plot,
        "E_LCDM": np.interp(z_plot, arr["z"], arr["E_lcdm"]),
        "E_Flux": np.interp(z_plot, arr["z"], arr["E_flux"]),
        "w_Flux": np.interp(z_plot, arr["z"], arr["w"]),
        "Omega_DE_Flux": np.interp(z_plot, arr["z"], arr["omega_de"]),
        "zfit": np.interp(np.linspace(0, 1, len(z_plot)), np.linspace(0, 1, len(zfit)), zfit),
        "H0_eff_curve": np.interp(np.linspace(0, 1, len(z_plot)), np.linspace(0, 1, len(h0_curve)), h0_curve),
    })


def save_plot(curve, best, outdir):
    plt.figure(figsize=(12, 9))

    plt.subplot(2, 2, 1)
    plt.plot(curve["z"], curve["E_LCDM"], label="LCDM E(z)")
    plt.plot(curve["z"], curve["E_Flux"], label="Flux E(z)")
    plt.xlim(0, 3)
    plt.xlabel("redshift z")
    plt.ylabel("E(z)=H/H0")
    plt.title("Background expansion")
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.plot(curve["z"], curve["w_Flux"], label="Flux w(z)")
    plt.axhline(-1.0, linestyle="--", label="Lambda")
    plt.xlim(0, 3)
    plt.xlabel("redshift z")
    plt.ylabel("w(z)")
    plt.title("Phantom-lite fossil hump")
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(curve["zfit"], curve["H0_eff_curve"], label="LCDM-inferred H0 if Flux true")
    plt.axhline(H0_CMB_REF, linestyle="--", label="CMB ref")
    plt.axhline(H0_LOCAL_TARGET, linestyle=":", label="local target")
    plt.xlim(0.05, 2.0)
    plt.xlabel("maximum fit redshift")
    plt.ylabel("effective H0 [km/s/Mpc]")
    plt.title("Hubble whisper sensitivity")
    plt.legend()

    plt.subplot(2, 2, 4)
    labels = ["CMB ref", "Flux low-z effective", "Local target"]
    vals = [H0_CMB_REF, float(best["H0_mean_lowz"]), H0_LOCAL_TARGET]
    plt.bar(labels, vals)
    plt.ylabel("H0 [km/s/Mpc]")
    plt.title("Best scan point")

    plt.tight_layout()
    plt.savefig(outdir / "flux_h0_sensitivity.png", dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = scan_parameter_space()
    best = df.iloc[0]
    curve = make_best_curve(best)

    verdict = "PARTIAL_H0_BRIDGE" if best["H0_mean_lowz"] > H0_CMB_REF + 1.0 and best["cmb_safe"] else "NO_H0_BRIDGE"
    if best["H0_mean_lowz"] >= 72.0 and best["cmb_safe"] and best["w0"] > -1.12:
        verdict = "STRONG_TOY_H0_BRIDGE"

    summary = pd.DataFrame([{
        "verdict": verdict,
        "best_w0_offset": float(best["w0_offset"]),
        "best_hump_amp": float(best["hump_amp"]),
        "best_w0": float(best["w0"]),
        "best_w_z2": float(best["w_z2"]),
        "H0_CMB_ref": H0_CMB_REF,
        "H0_local_target": H0_LOCAL_TARGET,
        "best_H0_mean_lowz": float(best["H0_mean_lowz"]),
        "delta_H0_lowz": float(best["delta_H0_lowz"]),
        "best_H0_eff_z0p1": float(best["H0_eff_z0p1"]),
        "best_H0_eff_z0p3": float(best["H0_eff_z0p3"]),
        "best_H0_eff_z0p5": float(best["H0_eff_z0p5"]),
        "omega_de_1090": float(best["omega_de_1090"]),
        "interpretation": (
            "Late phantom-lite Flux-DE can bias low-z LCDM distance fits toward a higher H0 while remaining CMB-era quiet. "
            "This is a sensitivity bridge, not a full H0-resolution claim."
        ),
    }])

    df.to_csv(outdir / "flux_h0_sensitivity_results.csv", index=False)
    summary.to_csv(outdir / "flux_h0_sensitivity_summary.csv", index=False)
    curve.to_csv(outdir / "flux_h0_sensitivity_best_curve.csv", index=False)
    save_plot(curve, best, outdir)

    print("\nFlux H0 sensitivity scan")
    print("=" * 72)
    print(summary.to_string(index=False))
    print("\nTop 10 scan points")
    cols = ["w0", "w_z2", "H0_mean_lowz", "delta_H0_lowz", "omega_de_1090", "w0_offset", "hump_amp"]
    print(df[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
