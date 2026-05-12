#!/usr/bin/env python3
"""
flux_cbh_h0_tension_test.py

Purpose:
    Test whether the CBH-driven Flux expansion history changes the CMB-inferred H0.

Input:
    flux_cbh_schechter_entropy.csv

This file should contain:
    z
    Omega_flux_4term_fit
    S_CBH_density_kB_mpc3
    rho_star_msun_mpc3
    sfrd_msun_yr_mpc3
    dS_CBH_dt_kB_mpc3_Gyr

Method:
    Replace constant Λ density evolution with a normalized dynamic Flux density:

        rho_DE(z) = rho_DE0 * F_flux(z)

    where F_flux(0) = 1.

    Then scan H0 while holding physical matter densities fixed:

        omega_b = Ω_b h^2
        omega_m = Ω_m h^2
        omega_r = Ω_r h^2

    and compute:

        theta_star = r_s(z*) / D_M(z*)

    Compare dynamic Flux-CBH model against vanilla ΛCDM.
"""

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

# CMB target acoustic angle, approximate Planck value.
# This is enough for a toy diagnostic, not a CLASS/CAMB replacement.
THETA_STAR_TARGET = 0.0104132

Z_STAR = 1089.92

# Physical densities, Planck-like baseline.
# omega_x = Omega_x h^2
OMEGA_B_H2 = 0.02237
OMEGA_M_H2 = 0.1430

# Radiation physical density including photons + relativistic neutrino correction.
# Approx value for Tcmb ~2.7255 and Neff ~3.046.
OMEGA_R_H2 = 4.18e-5

# Scan range
H0_MIN = 55.0
H0_MAX = 80.0
N_H0 = 600


# ============================================================
# 2. Load Flux-CBH history
# ============================================================

def load_flux_history(path="flux_cbh_schechter_entropy.csv"):
    df = pd.read_csv(path)

    if "Omega_flux_4term_fit" not in df.columns:
        raise ValueError("CSV must contain Omega_flux_4term_fit")

    # Sort ascending z for interpolation
    df = df.sort_values("z").reset_index(drop=True)

    z = df["z"].values
    F = df["Omega_flux_4term_fit"].values.astype(float)

    # Normalize to today
    F0 = F[np.argmin(np.abs(z - 0.0))]
    if F0 <= 0:
        raise ValueError("Flux history has non-positive F(0).")
    F = F / F0

    # Force non-negative
    F = np.clip(F, 0.0, None)

    # Interpolation.
    # For z above simulation max, hold the earliest high-z value.
    # Usually this should be tiny.
    interp = interp1d(
        z,
        F,
        kind="linear",
        bounds_error=False,
        fill_value=(F[0], F[-1]),
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

    # Flat universe closure
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
    return np.sqrt(
        p["Om_r"] * (1.0 + z)**4
        + p["Om_m"] * (1.0 + z)**3
        + p["Om_de"]
    )


def E_flux(z, p, flux_interp):
    F = flux_interp(z)
    F = np.clip(F, 0.0, None)

    return np.sqrt(
        p["Om_r"] * (1.0 + z)**4
        + p["Om_m"] * (1.0 + z)**3
        + p["Om_de"] * F
    )


# ============================================================
# 4. Sound horizon and distance
# ============================================================

def sound_speed_over_c(z):
    """
    Photon-baryon sound speed / c.

    c_s/c = 1 / sqrt[3(1 + R)]
    R = 3 rho_b / 4 rho_gamma

    Approx:
        R(z) = 31500 * omega_b * (Tcmb/2.7)^-4 / (1+z)
    """
    R = 31500.0 * OMEGA_B_H2 / (1.0 + z)
    return 1.0 / np.sqrt(3.0 * (1.0 + R))


def compute_theta_star(H0, model="lcdm", flux_interp=None):
    p = params_from_H0(H0)
    if p is None:
        return np.nan, np.nan, np.nan

    # Distance integral: 0 -> z*
    z_dist = np.linspace(0.0, Z_STAR, 12000)

    if model == "lcdm":
        E_dist = E_lcdm(z_dist, p)
    elif model == "flux":
        E_dist = E_flux(z_dist, p, flux_interp)
    else:
        raise ValueError("model must be 'lcdm' or 'flux'")

    D_M = (C_KM_S / H0) * trapezoid(1.0 / E_dist, z_dist)

    # Sound horizon integral: z* -> high z
    # Use logarithmic grid for early universe.
    z_hi = np.logspace(np.log10(Z_STAR), 6.0, 16000)

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
# 5. Scan and plots
# ============================================================

def scan_H0(model="lcdm", flux_interp=None):
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
    })


def make_plots(df_flux_history, scan_lcdm, scan_flux, best_lcdm, best_flux):
    plt.style.use("dark_background")

    # Flux history
    plt.figure(figsize=(10, 6))
    plt.plot(
        df_flux_history["z"],
        df_flux_history["Omega_flux_4term_fit"]
        / df_flux_history["Omega_flux_4term_fit"].iloc[-1],
        color="white",
        linewidth=2.5,
        label=r"Normalized Flux-CBH density history $F_{\rm flux}(z)$"
    )
    plt.gca().invert_xaxis()
    plt.xlim(8, 0)
    plt.xlabel("Redshift z")
    plt.ylabel("F(z), normalized to today")
    plt.title("Dynamic CBH-Flux Dark Sector History")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig("flux_cbh_dynamic_density_history.png", dpi=180)

    # theta scan
    plt.figure(figsize=(10, 6))
    plt.plot(scan_lcdm["H0"], scan_lcdm["theta_star"], color="cyan", label="ΛCDM")
    plt.plot(scan_flux["H0"], scan_flux["theta_star"], color="magenta", label="Flux-CBH")
    plt.axhline(THETA_STAR_TARGET, color="white", linestyle="--", label="Target θ*")
    plt.axvline(best_lcdm["H0"], color="cyan", linestyle=":")
    plt.axvline(best_flux["H0"], color="magenta", linestyle=":")
    plt.xlabel(r"$H_0$ [km/s/Mpc]")
    plt.ylabel(r"$\theta_\star = r_s / D_M$")
    plt.title("CMB Acoustic Angle H0 Scan")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig("flux_cbh_h0_theta_scan.png", dpi=180)

    # distance and sound horizon
    plt.figure(figsize=(10, 6))
    plt.plot(scan_lcdm["H0"], scan_lcdm["r_s_Mpc"], color="cyan", linestyle="--", label=r"ΛCDM $r_s$")
    plt.plot(scan_flux["H0"], scan_flux["r_s_Mpc"], color="magenta", linestyle="--", label=r"Flux-CBH $r_s$")
    plt.plot(scan_lcdm["H0"], scan_lcdm["D_M_Mpc"] / 100.0, color="cyan", label=r"ΛCDM $D_M/100$")
    plt.plot(scan_flux["H0"], scan_flux["D_M_Mpc"] / 100.0, color="magenta", label=r"Flux-CBH $D_M/100$")
    plt.xlabel(r"$H_0$ [km/s/Mpc]")
    plt.ylabel("Mpc scale")
    plt.title("Sound Horizon and CMB Distance Response")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig("flux_cbh_h0_distance_response.png", dpi=180)

    plt.close("all")


# ============================================================
# 6. Main
# ============================================================

def main():
    df_flux_history, flux_interp = load_flux_history("flux_cbh_schechter_entropy.csv")

    best_lcdm = find_best_H0(model="lcdm")
    best_flux = find_best_H0(model="flux", flux_interp=flux_interp)

    scan_lcdm = scan_H0(model="lcdm")
    scan_flux = scan_H0(model="flux", flux_interp=flux_interp)

    output = pd.concat([scan_lcdm, scan_flux], ignore_index=True)
    output.to_csv("flux_cbh_h0_scan.csv", index=False)

    make_plots(df_flux_history, scan_lcdm, scan_flux, best_lcdm, best_flux)

    print("\nFlux-CBH H0 Tension Diagnostic")
    print("=" * 78)
    print("Approximate acoustic-angle scan, not a CLASS/CAMB replacement.")
    print()
    print("Best-fit H0 from theta_star:")
    print(f"  ΛCDM:     H0 = {best_lcdm['H0']:.3f} km/s/Mpc")
    print(f"            theta* = {best_lcdm['theta']:.8f}")
    print(f"            r_s    = {best_lcdm['r_s_Mpc']:.3f} Mpc")
    print(f"            D_M    = {best_lcdm['D_M_Mpc']:.3f} Mpc")
    print()
    print(f"  Flux-CBH: H0 = {best_flux['H0']:.3f} km/s/Mpc")
    print(f"            theta* = {best_flux['theta']:.8f}")
    print(f"            r_s    = {best_flux['r_s_Mpc']:.3f} Mpc")
    print(f"            D_M    = {best_flux['D_M_Mpc']:.3f} Mpc")
    print()
    print(f"  ΔH0 = H0_Flux - H0_LCDM = {best_flux['H0'] - best_lcdm['H0']:.3f} km/s/Mpc")
    print()
    print("Wrote:")
    print("  flux_cbh_h0_scan.csv")
    print("  flux_cbh_dynamic_density_history.png")
    print("  flux_cbh_h0_theta_scan.png")
    print("  flux_cbh_h0_distance_response.png")


if __name__ == "__main__":
    main()