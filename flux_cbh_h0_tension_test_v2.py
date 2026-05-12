#!/usr/bin/env python3
"""
flux_cbh_h0_tension_test_v2.py

Updated Flux-CBH H0 tension diagnostic.

Fixes / upgrades:
    1. Fixes the plotting normalization bug:
        F(0) is always normalized to 1 using the row closest to z=0.

    2. Extends H0 scan:
        H0 range now 55 -> 100 km/s/Mpc.

    3. Adds nonlinear/compressed Flux response:
        F_metric(z) = F_flux(z)^p

       This tests whether the metric responds to a compressed/screened
       CBH entropy register rather than linearly to raw entropy.

    4. Scans several FLUX_POWER values and reports which one gives an
       acoustic-angle inferred H0 closest to a local target, default 73.

Input:
    flux_cbh_schechter_entropy.csv

Required column:
    Omega_flux_4term_fit

Outputs:
    flux_cbh_h0_scan_v2.csv
    flux_cbh_flux_power_summary.csv
    flux_cbh_dynamic_density_history_v2.png
    flux_cbh_h0_theta_scan_v2.png
    flux_cbh_h0_distance_response_v2.png
    flux_cbh_power_vs_h0_v2.png

Important:
    This is still an approximate acoustic-angle diagnostic.
    It is not a replacement for CLASS/CAMB.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar


# ============================================================
# 1. Constants
# ============================================================

C_KM_S = 299792.458

# Approx Planck acoustic angle target
THETA_STAR_TARGET = 0.0104132
Z_STAR = 1089.92

# Physical densities, Planck-like baseline
OMEGA_B_H2 = 0.02237
OMEGA_M_H2 = 0.1430

# Radiation physical density, photons + relativistic neutrino correction
OMEGA_R_H2 = 4.18e-5

# H0 scan range
H0_MIN = 55.0
H0_MAX = 100.0
N_H0 = 1000

# Local H0 target for the power summary only
LOCAL_H0_TARGET = 73.0

# Flux powers to test.
# p=1.0 means raw linear response.
# p<1.0 means compressed/screened metric response.
FLUX_POWER_GRID = [
    0.05,
    0.075,
    0.10,
    0.125,
    0.15,
    0.175,
    0.20,
    0.25,
    0.30,
    0.40,
    0.50,
    0.75,
    1.00,
]


# ============================================================
# 2. Loading and normalizing Flux history
# ============================================================

def normalize_flux_history(z, F_raw, flux_power=1.0):
    """
    Build metric response history:

        F_metric(z) = [F_raw(z) / F_raw(0)]^p

    then renormalize again so F_metric(0)=1.

    This makes the present-day dark sector density fixed while the
    redshift history changes shape.
    """
    z = np.asarray(z, dtype=float)
    F = np.asarray(F_raw, dtype=float)

    if len(z) != len(F):
        raise ValueError("z and F_raw must have same length.")

    idx0 = np.argmin(np.abs(z - 0.0))
    F0 = F[idx0]

    if not np.isfinite(F0) or F0 <= 0:
        raise ValueError(f"Flux history has invalid F(0): {F0}")

    # Normalize raw history to today
    F = F / F0

    # Prevent negative values from breaking fractional powers
    F = np.clip(F, 0.0, None)

    # Compressed/screened response
    F = F ** flux_power

    # Enforce F(0)=1 exactly after power transform
    F0_metric = F[idx0]
    if F0_metric <= 0:
        raise ValueError("Metric flux history has non-positive F_metric(0).")

    F = F / F0_metric

    return F


def load_flux_history(path="flux_cbh_schechter_entropy.csv", flux_power=1.0):
    path = Path(path)

    if not path.exists():
        print("\nERROR: Required CSV not found.")
        print(f"Looking for: {path.resolve()}")
        print(f"Current working folder: {Path.cwd()}")
        print("\nFix:")
        print("  1. Copy flux_cbh_schechter_entropy.csv into this folder, or")
        print("  2. Run flux_cbh_entropy_schechter.py first, or")
        print("  3. Use a full Windows path in load_flux_history(...).")
        raise FileNotFoundError(path)

    df = pd.read_csv(path)

    if "z" not in df.columns:
        raise ValueError("CSV must contain column: z")

    if "Omega_flux_4term_fit" not in df.columns:
        raise ValueError("CSV must contain column: Omega_flux_4term_fit")

    # Sort ascending z for interpolation
    df = df.sort_values("z").reset_index(drop=True)

    z = df["z"].values
    F_raw = df["Omega_flux_4term_fit"].values.astype(float)

    F_metric = normalize_flux_history(z, F_raw, flux_power=flux_power)

    df["F_flux_raw_norm_today"] = normalize_flux_history(z, F_raw, flux_power=1.0)
    df[f"F_metric_power_{flux_power:g}"] = F_metric

    interp = interp1d(
        z,
        F_metric,
        kind="linear",
        bounds_error=False,
        # Above max z, hold earliest value.
        # Below min z, hold today value.
        fill_value=(F_metric[0], F_metric[-1]),
    )

    return df, interp


# ============================================================
# 3. Background expansion models
# ============================================================

def params_from_H0(H0):
    h = H0 / 100.0

    Om_b = OMEGA_B_H2 / h**2
    Om_m = OMEGA_M_H2 / h**2
    Om_r = OMEGA_R_H2 / h**2

    # Flat closure
    Om_de = 1.0 - Om_m - Om_r

    if Om_de <= 0:
        return None

    return {
        "h": h,
        "H0": H0,
        "Om_b": Om_b,
        "Om_m": Om_m,
        "Om_r": Om_r,
        "Om_de": Om_de,
    }


def E_lcdm(z, p):
    z = np.asarray(z, dtype=float)
    return np.sqrt(
        p["Om_r"] * (1.0 + z)**4
        + p["Om_m"] * (1.0 + z)**3
        + p["Om_de"]
    )


def E_flux(z, p, flux_interp):
    z = np.asarray(z, dtype=float)
    F = flux_interp(z)
    F = np.clip(F, 0.0, None)

    return np.sqrt(
        p["Om_r"] * (1.0 + z)**4
        + p["Om_m"] * (1.0 + z)**3
        + p["Om_de"] * F
    )


# ============================================================
# 4. Sound horizon and CMB distance
# ============================================================

def sound_speed_over_c(z):
    """
    Photon-baryon sound speed / c.

    c_s/c = 1 / sqrt[3(1 + R)]

    Approx:
        R(z) = 31500 * omega_b / (1+z)
    """
    z = np.asarray(z, dtype=float)
    R = 31500.0 * OMEGA_B_H2 / (1.0 + z)
    return 1.0 / np.sqrt(3.0 * (1.0 + R))


def compute_theta_star(H0, model="lcdm", flux_interp=None):
    p = params_from_H0(H0)
    if p is None:
        return np.nan, np.nan, np.nan

    # Distance to last scattering: 0 -> z*
    z_dist = np.linspace(0.0, Z_STAR, 14000)

    if model == "lcdm":
        E_dist = E_lcdm(z_dist, p)
    elif model == "flux":
        if flux_interp is None:
            raise ValueError("flux_interp required for model='flux'")
        E_dist = E_flux(z_dist, p, flux_interp)
    else:
        raise ValueError("model must be 'lcdm' or 'flux'")

    D_M = (C_KM_S / H0) * trapezoid(1.0 / E_dist, z_dist)

    # Sound horizon: z* -> high z
    z_hi = np.logspace(np.log10(Z_STAR), 6.0, 18000)

    if model == "lcdm":
        E_hi = E_lcdm(z_hi, p)
    else:
        E_hi = E_flux(z_hi, p, flux_interp)

    cs = sound_speed_over_c(z_hi)
    r_s = (C_KM_S / H0) * trapezoid(cs / E_hi, z_hi)

    theta = r_s / D_M

    return theta, r_s, D_M


def find_best_H0(model="lcdm", flux_interp=None):
    def loss(H0):
        theta, _, _ = compute_theta_star(H0, model=model, flux_interp=flux_interp)
        if not np.isfinite(theta):
            return 1e99
        return (theta - THETA_STAR_TARGET)**2

    result = minimize_scalar(loss, bounds=(H0_MIN, H0_MAX), method="bounded")

    best_H0 = result.x
    theta, rs, DM = compute_theta_star(best_H0, model=model, flux_interp=flux_interp)

    return {
        "H0": best_H0,
        "theta": theta,
        "r_s_Mpc": rs,
        "D_M_Mpc": DM,
        "loss": result.fun,
    }


# ============================================================
# 5. Scans
# ============================================================

def scan_H0(model="lcdm", flux_interp=None, flux_power=np.nan):
    H0_grid = np.linspace(H0_MIN, H0_MAX, N_H0)

    theta_vals = []
    rs_vals = []
    DM_vals = []

    for H0 in H0_grid:
        theta, rs, DM = compute_theta_star(H0, model=model, flux_interp=flux_interp)
        theta_vals.append(theta)
        rs_vals.append(rs)
        DM_vals.append(DM)

    return pd.DataFrame({
        "H0": H0_grid,
        "theta_star": theta_vals,
        "r_s_Mpc": rs_vals,
        "D_M_Mpc": DM_vals,
        "model": model,
        "flux_power": flux_power,
    })


def scan_flux_powers(csv_path="flux_cbh_schechter_entropy.csv"):
    rows = []

    best_lcdm = find_best_H0(model="lcdm")

    for power in FLUX_POWER_GRID:
        _, flux_interp = load_flux_history(csv_path, flux_power=power)
        best_flux = find_best_H0(model="flux", flux_interp=flux_interp)

        rows.append({
            "flux_power": power,
            "H0_flux": best_flux["H0"],
            "theta_flux": best_flux["theta"],
            "r_s_flux_Mpc": best_flux["r_s_Mpc"],
            "D_M_flux_Mpc": best_flux["D_M_Mpc"],
            "H0_lcdm": best_lcdm["H0"],
            "delta_H0_vs_LCDM": best_flux["H0"] - best_lcdm["H0"],
            "distance_to_local_H0_target": abs(best_flux["H0"] - LOCAL_H0_TARGET),
            "loss": best_flux["loss"],
        })

    summary = pd.DataFrame(rows)
    summary = summary.sort_values("distance_to_local_H0_target").reset_index(drop=True)
    return best_lcdm, summary


# ============================================================
# 6. Plotting
# ============================================================

def make_plots(
    df_flux_history,
    scan_lcdm,
    scan_flux_best_power,
    best_lcdm,
    best_flux,
    power_summary,
    best_power,
    prefix="flux_cbh",
):
    plt.style.use("dark_background")

    z = df_flux_history["z"].values
    F_raw_norm = df_flux_history["F_flux_raw_norm_today"].values
    F_metric = df_flux_history[f"F_metric_power_{best_power:g}"].values

    # --------------------------------------------------------
    # Flux history, fixed normalization
    # --------------------------------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(
        z,
        F_raw_norm,
        color="gray",
        linewidth=1.8,
        linestyle="--",
        label=r"Raw 4-term Flux history, normalized $F(0)=1$"
    )
    plt.plot(
        z,
        F_metric,
        color="white",
        linewidth=2.5,
        label=rf"Metric response $F(z)^p$, $p={best_power:g}$"
    )
    plt.axhline(1.0, color="white", alpha=0.25)
    plt.gca().invert_xaxis()
    plt.xlim(8, 0)
    plt.xlabel("Redshift z")
    plt.ylabel("F(z), normalized to today")
    plt.title("Dynamic CBH-Flux Dark Sector History")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{prefix}_dynamic_density_history_v2.png", dpi=180)

    # --------------------------------------------------------
    # H0 theta scan
    # --------------------------------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(scan_lcdm["H0"], scan_lcdm["theta_star"], color="cyan", label="ΛCDM")
    plt.plot(
        scan_flux_best_power["H0"],
        scan_flux_best_power["theta_star"],
        color="magenta",
        label=rf"Flux-CBH, $p={best_power:g}$"
    )
    plt.axhline(THETA_STAR_TARGET, color="white", linestyle="--", label="Target θ*")
    plt.axvline(best_lcdm["H0"], color="cyan", linestyle=":")
    plt.axvline(best_flux["H0"], color="magenta", linestyle=":")
    plt.xlabel(r"$H_0$ [km/s/Mpc]")
    plt.ylabel(r"$\theta_\star = r_s / D_M$")
    plt.title("CMB Acoustic Angle H0 Scan")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{prefix}_h0_theta_scan_v2.png", dpi=180)

    # --------------------------------------------------------
    # Distance response
    # --------------------------------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(scan_lcdm["H0"], scan_lcdm["r_s_Mpc"], color="cyan", linestyle="--", label=r"ΛCDM $r_s$")
    plt.plot(
        scan_flux_best_power["H0"],
        scan_flux_best_power["r_s_Mpc"],
        color="magenta",
        linestyle="--",
        label=rf"Flux-CBH $r_s$, $p={best_power:g}$"
    )
    plt.plot(scan_lcdm["H0"], scan_lcdm["D_M_Mpc"] / 100.0, color="cyan", label=r"ΛCDM $D_M/100$")
    plt.plot(
        scan_flux_best_power["H0"],
        scan_flux_best_power["D_M_Mpc"] / 100.0,
        color="magenta",
        label=rf"Flux-CBH $D_M/100$, $p={best_power:g}$"
    )
    plt.xlabel(r"$H_0$ [km/s/Mpc]")
    plt.ylabel("Mpc scale")
    plt.title("Sound Horizon and CMB Distance Response")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{prefix}_h0_distance_response_v2.png", dpi=180)

    # --------------------------------------------------------
    # Power vs inferred H0
    # --------------------------------------------------------
    sorted_power = power_summary.sort_values("flux_power")
    plt.figure(figsize=(10, 6))
    plt.plot(
        sorted_power["flux_power"],
        sorted_power["H0_flux"],
        marker="o",
        linewidth=2,
        color="lime",
        label="Flux-CBH inferred H0"
    )
    plt.axhline(best_lcdm["H0"], color="cyan", linestyle="--", label=f"ΛCDM H0={best_lcdm['H0']:.2f}")
    plt.axhline(LOCAL_H0_TARGET, color="white", linestyle=":", label=f"Local target H0={LOCAL_H0_TARGET:.1f}")
    plt.axvline(best_power, color="magenta", linestyle=":", label=f"Best p={best_power:g}")
    plt.xlabel(r"Flux compression power $p$")
    plt.ylabel(r"Inferred $H_0$ [km/s/Mpc]")
    plt.title("Flux-CBH Compression Power vs Acoustic-Angle H0")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{prefix}_power_vs_h0_v2.png", dpi=180)

    plt.close("all")


# ============================================================
# 7. Main
# ============================================================

def main():
    csv_path = "flux_cbh_schechter_entropy.csv"

    # Power scan summary
    best_lcdm, power_summary = scan_flux_powers(csv_path)
    power_summary.to_csv("flux_cbh_flux_power_summary.csv", index=False)

    best_row = power_summary.iloc[0]
    best_power = float(best_row["flux_power"])

    # Build best-power history and scans
    df_flux_history, flux_interp = load_flux_history(csv_path, flux_power=best_power)
    best_flux = find_best_H0(model="flux", flux_interp=flux_interp)

    scan_lcdm = scan_H0(model="lcdm")
    scan_flux_best = scan_H0(model="flux", flux_interp=flux_interp, flux_power=best_power)

    output = pd.concat([scan_lcdm, scan_flux_best], ignore_index=True)
    output.to_csv("flux_cbh_h0_scan_v2.csv", index=False)

    make_plots(
        df_flux_history=df_flux_history,
        scan_lcdm=scan_lcdm,
        scan_flux_best_power=scan_flux_best,
        best_lcdm=best_lcdm,
        best_flux=best_flux,
        power_summary=power_summary,
        best_power=best_power,
    )

    print("\nFlux-CBH H0 Tension Diagnostic v2")
    print("=" * 82)
    print("Approximate acoustic-angle scan, not a CLASS/CAMB replacement.")
    print()
    print("Best-fit LCDM from theta_star:")
    print(f"  H0     = {best_lcdm['H0']:.3f} km/s/Mpc")
    print(f"  theta* = {best_lcdm['theta']:.8f}")
    print(f"  r_s    = {best_lcdm['r_s_Mpc']:.3f} Mpc")
    print(f"  D_M    = {best_lcdm['D_M_Mpc']:.3f} Mpc")
    print()
    print("Flux power scan, closest to local H0 target:")
    print(power_summary[[
        "flux_power",
        "H0_flux",
        "delta_H0_vs_LCDM",
        "distance_to_local_H0_target",
        "r_s_flux_Mpc",
        "D_M_flux_Mpc",
    ]].head(8).to_string(index=False, float_format=lambda x: f"{x:.6f}"))
    print()
    print(f"Selected best compression power p = {best_power:g}")
    print()
    print("Best selected Flux-CBH model:")
    print(f"  H0     = {best_flux['H0']:.3f} km/s/Mpc")
    print(f"  theta* = {best_flux['theta']:.8f}")
    print(f"  r_s    = {best_flux['r_s_Mpc']:.3f} Mpc")
    print(f"  D_M    = {best_flux['D_M_Mpc']:.3f} Mpc")
    print(f"  Delta H0 = {best_flux['H0'] - best_lcdm['H0']:.3f} km/s/Mpc")
    print()
    print("Wrote:")
    print("  flux_cbh_h0_scan_v2.csv")
    print("  flux_cbh_flux_power_summary.csv")
    print("  flux_cbh_dynamic_density_history_v2.png")
    print("  flux_cbh_h0_theta_scan_v2.png")
    print("  flux_cbh_h0_distance_response_v2.png")
    print("  flux_cbh_power_vs_h0_v2.png")


if __name__ == "__main__":
    main()
