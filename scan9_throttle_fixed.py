#!/usr/bin/env python3
"""
scan9_throttle_fixed.py

Flux-CBH 5D throttle scan.

This is the fixed version of scan9.py.

Fixes:
    1. Uses Omega_flux_4term_fit as the preferred raw Flux history,
       falling back to S_CBH_density_kB_mpc3 only if needed.

    2. Enforces F_raw(0)=1 and F_metric(0)=1.

    3. Saves all tested models, nearest CMB misses, and best passing models.

    4. Uses vectorized trapezoid integrations instead of many scipy.quad calls,
       making the scan faster and easier to debug.

Model:
    p_eff(z) =
        p_high + (p0 - p_high) / [1 + exp((z - z_t)/w)]

    F_metric(z) =
        1 + A * (F_raw(z)^p_eff(z) - 1)

    E^2(z) =
        Omega_r(1+z)^4
      + Omega_m(1+z)^3
      + Omega_flux,0 F_metric(z)

Inputs:
    flux_cbh_schechter_entropy.csv

Required columns:
    z
    Omega_flux_4term_fit

Fallback column:
    S_CBH_density_kB_mpc3

Outputs:
    scan9_throttle_all_models.csv
    scan9_throttle_best.csv
    scan9_throttle_nearest_cmb_misses.csv
    scan9_throttle_best_plot.png
    scan9_throttle_A_vs_SN.png
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d


# =====================================================================
# 1. PARAMETERS & CONSTANTS
# =====================================================================

C_KM_S = 299792.458
Z_STAR = 1089.92
THETA_STAR_TARGET = 0.0104132

# Flux-CBH best matter-float coordinates
H0_FLUX = 73.0
OMEGA_M_H2_FLUX = 0.144
OMEGA_B_H2_FLUX = 0.02237
OMEGA_R_H2_FLUX = 4.183e-5

# Photon density only for sound speed baryon-loading term
OMEGA_GAMMA_H2 = 2.469e-5

h_flux = H0_FLUX / 100.0
OMEGA_M_FLUX = OMEGA_M_H2_FLUX / h_flux**2
OMEGA_R_FLUX = OMEGA_R_H2_FLUX / h_flux**2
OMEGA_FLUX_0 = 1.0 - OMEGA_M_FLUX - OMEGA_R_FLUX

R_B_FACTOR = (3.0 * OMEGA_B_H2_FLUX) / (4.0 * OMEGA_GAMMA_H2)

# LCDM baseline for SN-shape comparison only
H0_LCDM = 67.4
OMEGA_M_LCDM = 0.315
OMEGA_L_LCDM = 1.0 - OMEGA_M_LCDM

# CMB gates
CMB_STRICT_PPM = 800.0
CMB_SOFT_PPM = 2500.0

# SN shape gates
SN_GOOD_MAG = 0.10
SN_EXCELLENT_MAG = 0.05

# Integration grids
Z_CMB_DIST = np.concatenate([
    np.linspace(0.0, 8.0, 2200),
    np.linspace(8.0, Z_STAR, 2800),
])

Z_SOUND = np.logspace(np.log10(Z_STAR), 6.0, 5500)

Z_SN = np.linspace(0.01, 2.3, 120)
Z_SN_INT = np.linspace(0.0, 2.3, 1800)
ANCHOR_MASK = (Z_SN >= 0.01) & (Z_SN <= 0.05)

# Focused 5D scan around scan8 success
P0_VALS = np.array([0.170, 0.175, 0.180, 0.185, 0.190, 0.195, 0.200])
P_HIGH_VALS = np.array([0.200, 0.210, 0.220, 0.230, 0.240])
ZT_VALS = np.array([1.3, 1.5, 1.7, 1.9, 2.1])
W_VALS = np.array([0.25, 0.35, 0.45, 0.60, 0.80])
A_VALS = np.linspace(0.20, 1.20, 21)


# =====================================================================
# 2. LOAD FLUX HISTORY
# =====================================================================

def load_flux_history(path: str | Path = "flux_cbh_schechter_entropy.csv"):
    path = Path(path)

    if not path.exists():
        print("\nERROR: Required CSV not found.")
        print(f"Looking for: {path.resolve()}")
        print(f"Current folder: {Path.cwd()}")
        raise FileNotFoundError(path)

    df = pd.read_csv(path).sort_values("z").reset_index(drop=True)

    if "Omega_flux_4term_fit" in df.columns:
        source_col = "Omega_flux_4term_fit"
    elif "S_CBH_density_kB_mpc3" in df.columns:
        source_col = "S_CBH_density_kB_mpc3"
    else:
        raise ValueError(
            "CSV must contain either Omega_flux_4term_fit or S_CBH_density_kB_mpc3"
        )

    z = df["z"].values.astype(float)
    F = df[source_col].values.astype(float)

    good = np.isfinite(z) & np.isfinite(F) & (F > 0)
    z = z[good]
    F = F[good]

    order = np.argsort(z)
    z = z[order]
    F = F[order]

    idx0 = int(np.argmin(np.abs(z - 0.0)))
    F0 = F[idx0]

    if F0 <= 0:
        raise ValueError("F(0) is non-positive after cleaning.")

    F_norm = np.clip(F / F0, 1e-90, None)

    interp = interp1d(
        z,
        F_norm,
        kind="linear",
        bounds_error=False,
        # z ascending: below min gets F[0], above max gets F[-1]
        fill_value=(F_norm[0], F_norm[-1]),
    )

    print(f"Loaded Flux history from: {path}")
    print(f"Using source column: {source_col}")
    print(f"F_raw(0) normalized to 1.0")
    print(f"z range: {z.min():.3f} to {z.max():.3f}")

    return interp, source_col, df


F_RAW_INTERP, SOURCE_COL, DF_FLUX = load_flux_history()


def F_raw(z):
    return np.clip(F_RAW_INTERP(z), 1e-90, None)


# =====================================================================
# 3. MODEL FUNCTIONS
# =====================================================================

def p_eff(z, p0, p_high, zt, w):
    z = np.asarray(z, dtype=float)
    return p_high + (p0 - p_high) / (1.0 + np.exp((z - zt) / w))


def F_power(z, p0, p_high, zt, w):
    z = np.asarray(z, dtype=float)
    F = F_raw(z)
    pe = p_eff(z, p0, p_high, zt, w)
    out = F ** pe

    # enforce F_power(0)=1
    F00 = F_raw(0.0) ** p_eff(0.0, p0, p_high, zt, w)
    return out / F00


def F_metric(z, p0, p_high, zt, w, A):
    """
    Amplitude-throttled metric response.

        F_metric = 1 + A(F_power - 1)

    Because F_power(0)=1, this guarantees F_metric(0)=1.
    """
    Fp = F_power(z, p0, p_high, zt, w)
    Fm = 1.0 + A * (Fp - 1.0)

    # Avoid nonphysical negative/zero dark-sector density.
    return np.clip(Fm, 1e-12, None)


def E_flux(z, p0, p_high, zt, w, A):
    z = np.asarray(z, dtype=float)
    flux_term = OMEGA_FLUX_0 * F_metric(z, p0, p_high, zt, w, A)

    return np.sqrt(
        OMEGA_R_FLUX * (1.0 + z)**4
        + OMEGA_M_FLUX * (1.0 + z)**3
        + flux_term
    )


def E_lcdm(z):
    z = np.asarray(z, dtype=float)
    return np.sqrt(
        OMEGA_M_LCDM * (1.0 + z)**3
        + OMEGA_L_LCDM
    )


def sound_speed_over_c(z):
    z = np.asarray(z, dtype=float)
    return 1.0 / np.sqrt(3.0 * (1.0 + R_B_FACTOR / (1.0 + z)))


# =====================================================================
# 4. DISTANCE / CMB / SN EVALUATORS
# =====================================================================

def luminosity_distance_grid(z_grid, E_values, H0):
    invE = 1.0 / E_values

    Dc_int = np.zeros_like(z_grid)
    dz = np.diff(z_grid)
    avg = 0.5 * (invE[1:] + invE[:-1])
    Dc_int[1:] = np.cumsum(dz * avg)

    Dc = (C_KM_S / H0) * Dc_int
    Dl = (1.0 + z_grid) * Dc
    return Dl


def eval_cmb(p0, p_high, zt, w, A):
    E_dist = E_flux(Z_CMB_DIST, p0, p_high, zt, w, A)
    D_M = (C_KM_S / H0_FLUX) * trapezoid(1.0 / E_dist, Z_CMB_DIST)

    E_sound = E_flux(Z_SOUND, p0, p_high, zt, w, A)
    r_s = (C_KM_S / H0_FLUX) * trapezoid(
        sound_speed_over_c(Z_SOUND) / E_sound,
        Z_SOUND,
    )

    theta = r_s / D_M
    err_ppm = (theta - THETA_STAR_TARGET) / THETA_STAR_TARGET * 1e6

    return theta, err_ppm, r_s, D_M


# Precompute LCDM SN baseline
E_LCDM_SN_INT = E_lcdm(Z_SN_INT)
DL_LCDM_INT = luminosity_distance_grid(Z_SN_INT, E_LCDM_SN_INT, H0_LCDM)
MU_LCDM_INT = 5.0 * np.log10(np.clip(DL_LCDM_INT, 1e-30, None)) + 25.0
MU_LCDM = interp1d(
    Z_SN_INT,
    MU_LCDM_INT,
    bounds_error=False,
    fill_value="extrapolate",
)(Z_SN)


def eval_sn_shape(p0, p_high, zt, w, A):
    E_flux_sn_int = E_flux(Z_SN_INT, p0, p_high, zt, w, A)
    DL_flux_int = luminosity_distance_grid(Z_SN_INT, E_flux_sn_int, H0_FLUX)

    MU_flux_int = 5.0 * np.log10(np.clip(DL_flux_int, 1e-30, None)) + 25.0
    MU_flux = interp1d(
        Z_SN_INT,
        MU_flux_int,
        bounds_error=False,
        fill_value="extrapolate",
    )(Z_SN)

    delta = MU_flux - MU_LCDM
    offset = float(np.mean(delta[ANCHOR_MASK]))
    shape = delta - offset

    return {
        "sn_shape_max": float(np.max(np.abs(shape))),
        "sn_shape_end": float(shape[-1]),
        "sn_shape_rms": float(np.sqrt(np.mean(shape**2))),
        "sn_offset": offset,
        "shape_curve": shape,
    }


# =====================================================================
# 5. SANITY CHECKS
# =====================================================================

print("\nSanity checks:")
for ptest in [0.190, 0.214]:
    theta, err_ppm, rs, dm = eval_cmb(ptest, ptest, 1.7, 0.35, 1.0)
    sn = eval_sn_shape(ptest, ptest, 1.7, 0.35, 1.0)
    print(
        f"  constant p={ptest:.3f}, A=1.0 -> "
        f"theta={theta:.8f}, err={err_ppm:+.1f} ppm, "
        f"r_s={rs:.3f}, D_M={dm:.3f}, "
        f"SNmax={sn['sn_shape_max']:.4f}, SNrms={sn['sn_shape_rms']:.4f}"
    )


# =====================================================================
# 6. FOCUSED 5D SCAN
# =====================================================================

print("\nStarting 5D throttle scan")
print("=" * 78)
print(f"p0 values:      {P0_VALS[0]:.3f} -> {P0_VALS[-1]:.3f} ({len(P0_VALS)})")
print(f"p_high values:  {P_HIGH_VALS[0]:.3f} -> {P_HIGH_VALS[-1]:.3f} ({len(P_HIGH_VALS)})")
print(f"zt values:      {ZT_VALS[0]:.2f} -> {ZT_VALS[-1]:.2f} ({len(ZT_VALS)})")
print(f"w values:       {W_VALS[0]:.2f} -> {W_VALS[-1]:.2f} ({len(W_VALS)})")
print(f"A values:       {A_VALS[0]:.2f} -> {A_VALS[-1]:.2f} ({len(A_VALS)})")

n_total = len(P0_VALS) * len(P_HIGH_VALS) * len(ZT_VALS) * len(W_VALS) * len(A_VALS)
print(f"Total models:   {n_total:,}")
print()

start = time.time()

all_rows = []
best_rows = []
tested = 0
passed_strict = 0
passed_soft = 0

for p0 in P0_VALS:
    print(f"  p0={p0:.3f}")
    for p_high in P_HIGH_VALS:
        for zt in ZT_VALS:
            for w in W_VALS:
                for A in A_VALS:
                    tested += 1

                    theta, err_ppm, r_s, D_M = eval_cmb(p0, p_high, zt, w, A)
                    abs_err_ppm = abs(err_ppm)

                    row = {
                        "A": A,
                        "p0": p0,
                        "p_high": p_high,
                        "zt": zt,
                        "w": w,
                        "theta_star": theta,
                        "theta_err_ppm": err_ppm,
                        "theta_abs_err_ppm": abs_err_ppm,
                        "r_s_Mpc": r_s,
                        "D_M_Mpc": D_M,
                        "cmb_soft_pass": abs_err_ppm <= CMB_SOFT_PPM,
                        "cmb_strict_pass": abs_err_ppm <= CMB_STRICT_PPM,
                    }

                    if abs_err_ppm <= CMB_SOFT_PPM:
                        passed_soft += 1
                        sn = eval_sn_shape(p0, p_high, zt, w, A)

                        row.update({
                            "sn_shape_max": sn["sn_shape_max"],
                            "sn_shape_end": sn["sn_shape_end"],
                            "sn_shape_rms": sn["sn_shape_rms"],
                            "sn_offset": sn["sn_offset"],
                            "sn_good_pass": sn["sn_shape_max"] <= SN_GOOD_MAG,
                            "sn_excellent_pass": sn["sn_shape_max"] <= SN_EXCELLENT_MAG,
                        })

                        # Not a formal likelihood. Diagnostic sorting score.
                        row["score"] = (
                            (abs_err_ppm / CMB_STRICT_PPM)**2
                            + (sn["sn_shape_max"] / SN_GOOD_MAG)**2
                            + (sn["sn_shape_rms"] / 0.05)**2
                            + ((A - 0.75) / 0.50)**2
                        )

                        if abs_err_ppm <= CMB_STRICT_PPM:
                            passed_strict += 1
                            best_rows.append(row.copy())
                    else:
                        row.update({
                            "sn_shape_max": np.nan,
                            "sn_shape_end": np.nan,
                            "sn_shape_rms": np.nan,
                            "sn_offset": np.nan,
                            "sn_good_pass": False,
                            "sn_excellent_pass": False,
                            "score": (abs_err_ppm / CMB_STRICT_PPM)**2 + 99.0,
                        })

                    all_rows.append(row)


# =====================================================================
# 7. SAVE RESULTS
# =====================================================================

all_df = pd.DataFrame(all_rows).sort_values("score").reset_index(drop=True)
all_df.to_csv("scan9_throttle_all_models.csv", index=False)

nearest_cmb = all_df.sort_values("theta_abs_err_ppm").head(300).copy()
nearest_cmb.to_csv("scan9_throttle_nearest_cmb_misses.csv", index=False)

strict_df = all_df[all_df["cmb_strict_pass"]].copy()
if not strict_df.empty:
    best_df = strict_df.sort_values(["sn_shape_max", "theta_abs_err_ppm", "score"]).head(300).copy()
else:
    best_df = all_df.head(300).copy()

best_df.to_csv("scan9_throttle_best.csv", index=False)

elapsed = time.time() - start

print("\nScan complete")
print("=" * 78)
print(f"Elapsed:             {elapsed:.1f}s")
print(f"Tested:              {tested:,}")
print(f"Soft CMB pass:        {passed_soft:,}")
print(f"Strict CMB pass:      {passed_strict:,}")
print()
print("Top 12 models:")
cols = [
    "A", "p0", "p_high", "zt", "w",
    "theta_err_ppm", "sn_shape_max", "sn_shape_rms",
    "sn_shape_end", "score",
]
print(best_df.head(12)[cols].to_string(index=False, float_format=lambda x: f"{x:.6g}"))

print("\nWrote:")
print("  scan9_throttle_all_models.csv")
print("  scan9_throttle_best.csv")
print("  scan9_throttle_nearest_cmb_misses.csv")


# =====================================================================
# 8. PLOTS FOR BEST MODEL
# =====================================================================

best = best_df.iloc[0]

A_b = float(best["A"])
p0_b = float(best["p0"])
ph_b = float(best["p_high"])
zt_b = float(best["zt"])
w_b = float(best["w"])

z_plot = np.linspace(0.0, 8.0, 900)
F_raw_plot = F_raw(z_plot)
F_power_plot = F_power(z_plot, p0_b, ph_b, zt_b, w_b)
F_metric_plot = F_metric(z_plot, p0_b, ph_b, zt_b, w_b, A_b)
p_plot = p_eff(z_plot, p0_b, ph_b, zt_b, w_b)

sn_best = eval_sn_shape(p0_b, ph_b, zt_b, w_b, A_b)
shape_best = sn_best["shape_curve"]

plt.style.use("dark_background")
fig, axes = plt.subplots(4, 1, figsize=(11, 13), sharex=False)

axes[0].plot(z_plot, F_raw_plot, "--", color="gray", label="Raw Flux history F(0)=1")
axes[0].plot(z_plot, F_power_plot, color="white", alpha=0.75, linewidth=2, label="Unthrottled power response")
axes[0].plot(z_plot, F_metric_plot, color="lime", linewidth=2.5, label=rf"Throttled response A={A_b:.2f}")
axes[0].invert_xaxis()
axes[0].set_xlim(8, 0)
axes[0].set_ylabel("F(z)")
axes[0].set_title(
    f"Best scan9 throttle: A={A_b:.2f}, p0={p0_b:.3f}, p_high={ph_b:.3f}, zt={zt_b:.2f}, w={w_b:.2f}"
)
axes[0].grid(alpha=0.25)
axes[0].legend()

axes[1].plot(z_plot, p_plot, color="cyan", linewidth=2.5)
axes[1].invert_xaxis()
axes[1].set_xlim(8, 0)
axes[1].set_ylabel(r"$p_{\rm eff}(z)$")
axes[1].grid(alpha=0.25)

axes[2].axhline(0, color="white", linestyle="--", alpha=0.6)
axes[2].axhspan(-0.05, 0.05, color="green", alpha=0.15, label="±0.05 mag")
axes[2].axhspan(-0.10, 0.10, color="gray", alpha=0.15, label="±0.10 mag")
axes[2].plot(Z_SN, shape_best, color="magenta", linewidth=2.5, label="SN shape residual")
axes[2].set_xlim(0, 2.3)
axes[2].set_ylabel(r"$\Delta\mu_{\rm shape}$")
axes[2].grid(alpha=0.25)
axes[2].legend()

# A vs SN among strict CMB pass models
strict_plot = all_df[all_df["cmb_strict_pass"] & np.isfinite(all_df["sn_shape_max"])].copy()
if not strict_plot.empty:
    axes[3].scatter(
        strict_plot["A"],
        strict_plot["sn_shape_max"],
        c=strict_plot["theta_abs_err_ppm"],
        s=25,
        alpha=0.85,
    )
    axes[3].axhline(SN_GOOD_MAG, color="gray", linestyle="--", label="0.10 mag")
    axes[3].axhline(SN_EXCELLENT_MAG, color="green", linestyle=":", label="0.05 mag")
    axes[3].axvline(A_b, color="magenta", linestyle=":", label="Best A")
    axes[3].set_xlabel("Amplitude throttle A")
    axes[3].set_ylabel("SN max shape residual [mag]")
    axes[3].grid(alpha=0.25)
    axes[3].legend()
else:
    axes[3].text(0.1, 0.5, "No strict CMB pass models to plot", transform=axes[3].transAxes)
    axes[3].set_axis_off()

plt.tight_layout()
plt.savefig("scan9_throttle_best_plot.png", dpi=220)

# Separate A-vs-SN plot
if not strict_plot.empty:
    plt.figure(figsize=(9, 6))
    sc = plt.scatter(
        strict_plot["A"],
        strict_plot["sn_shape_max"],
        c=strict_plot["theta_abs_err_ppm"],
        s=35,
        alpha=0.85,
    )
    plt.colorbar(sc, label="theta_abs_err_ppm")
    plt.axhline(SN_GOOD_MAG, color="gray", linestyle="--", label="0.10 mag")
    plt.axhline(SN_EXCELLENT_MAG, color="green", linestyle=":", label="0.05 mag")
    plt.axvline(A_b, color="magenta", linestyle=":", label=f"Best A={A_b:.2f}")
    plt.xlabel("Amplitude throttle A")
    plt.ylabel("SN max shape residual [mag]")
    plt.title("Scan9 Throttle: A vs SN Shape Residual, CMB-Strict Models")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig("scan9_throttle_A_vs_SN.png", dpi=220)

print("  scan9_throttle_best_plot.png")
print("  scan9_throttle_A_vs_SN.png")
