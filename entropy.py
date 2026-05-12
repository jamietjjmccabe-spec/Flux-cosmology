#!/usr/bin/env python3
"""
madau_dickinson_stellar_entropy.py

Integrates the Madau-Dickinson cosmic star-formation-rate density curve:

    SFRD(z) = 0.015 * (1+z)^2.7 / [1 + ((1+z)/2.9)^5.6]
              Msun yr^-1 Mpc^-3

Then estimates:
    - accumulated stellar mass density rho_star(z)
    - stellar luminosity density proxy
    - cosmic fusion-chain rate density
    - stellar entropy-production rate density

Flux-use:
    This gives a numerical stellar-registration source term that can be
    compared against H(z), dH/dt, q(z), or a Flux entropy term dS/dt = gamma Phi.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

# ============================================================
# 1. Constants
# ============================================================

# Cosmology
H0_km_s_Mpc = 67.4
Omega_m0 = 0.315
Omega_L0 = 0.685

# Unit conversions
MPC_M = 3.0856775814913673e22
YEAR_S = 365.25 * 24 * 3600
MSUN_KG = 1.98847e30
L_SUN_W = 3.828e26
K_B = 1.380649e-23

# Fusion energy: 4p -> He chain, approx 26.7 MeV
MEV_J = 1.602176634e-13
Q_PP_CHAIN_J = 26.7 * MEV_J

# Approx stellar mass-to-light ratio.
# This is a crude population-averaged value.
# Lower = brighter young stellar population.
MASS_TO_LIGHT_SOLAR = 3.0  # Msun / Lsun

# Recycled mass fraction: fraction of stellar mass returned to gas
RETURN_FRACTION = 0.27

# Entropy sink temperature.
# 2.725*(1+z) for CMB floor.
T_CMB0 = 2.7255

# Effective stellar emission temperature
T_STAR_EFF = 5772.0


# ============================================================
# 2. Cosmology
# ============================================================

def H_z_SI(z):
    """
    H(z) in s^-1.
    """
    H0_SI = H0_km_s_Mpc * 1000.0 / MPC_M
    Ez = np.sqrt(Omega_m0 * (1.0 + z)**3 + Omega_L0)
    return H0_SI * Ez


def dt_dz(z):
    """
    dt/dz in seconds per redshift.
    Negative because time increases as z decreases.
    """
    return -1.0 / ((1.0 + z) * H_z_SI(z))


# ============================================================
# 3. Madau-Dickinson star-formation curve
# ============================================================

def sfrd_madau_dickinson(z):
    """
    Cosmic star formation rate density.

    Units:
        Msun yr^-1 Mpc^-3
    """
    return 0.015 * (1.0 + z)**2.7 / (1.0 + ((1.0 + z) / 2.9)**5.6)


# ============================================================
# 4. Derived stellar-registration terms
# ============================================================

def compute_history(z_max=20.0, n=6000):
    """
    Integrate from high redshift to today.

    Returns dataframe ordered from high z -> low z.
    """

    # High-z to low-z grid
    z = np.linspace(z_max, 0.0, n)

    sfrd = sfrd_madau_dickinson(z)  # Msun yr^-1 Mpc^-3

    # Convert dt/dz to positive elapsed time along high-z -> low-z direction
    abs_dt_dz_years = np.abs(dt_dz(z)) / YEAR_S

    # Stellar mass density formed:
    # d rho_star = SFRD * dt
    # include surviving locked mass fraction (1 - R)
    integrand_mass = sfrd * abs_dt_dz_years * (1.0 - RETURN_FRACTION)

    rho_star = cumulative_trapezoid(
        integrand_mass,
        z,
        initial=0.0
    )

    # Since z grid decreases, cumulative_trapezoid returns negative accumulation.
    # Flip sign to get positive accumulated stellar mass density.
    rho_star = -rho_star

    # Luminosity density proxy:
    # rho_star / M/L = Lsun Mpc^-3
    lum_density_lsun_mpc3 = rho_star / MASS_TO_LIGHT_SOLAR
    lum_density_w_mpc3 = lum_density_lsun_mpc3 * L_SUN_W

    # Fusion-chain rate density:
    # L / Q
    fusion_chain_rate_mpc3_s = lum_density_w_mpc3 / Q_PP_CHAIN_J

    # Proton processing rate density:
    proton_processing_rate_mpc3_s = 4.0 * fusion_chain_rate_mpc3_s

    # Entropy production from radiation degrading from star temp to sink temp.
    # Use CMB as minimum sink. You can later swap this for dust/IGM temp.
    T_sink = T_CMB0 * (1.0 + z)

    entropy_rate_J_K_mpc3_s = lum_density_w_mpc3 * (
        (1.0 / T_sink) - (1.0 / T_STAR_EFF)
    )

    entropy_rate_kB_mpc3_s = entropy_rate_J_K_mpc3_s / K_B

    # Cosmic age since z_max, useful as integration clock
    elapsed_years = cumulative_trapezoid(
        abs_dt_dz_years,
        z,
        initial=0.0
    )
    elapsed_years = -elapsed_years

    df = pd.DataFrame({
        "z": z,
        "elapsed_years_since_zmax": elapsed_years,
        "sfrd_msun_yr_mpc3": sfrd,
        "rho_star_msun_mpc3": rho_star,
        "lum_density_lsun_mpc3": lum_density_lsun_mpc3,
        "lum_density_w_mpc3": lum_density_w_mpc3,
        "fusion_chain_rate_mpc3_s": fusion_chain_rate_mpc3_s,
        "proton_processing_rate_mpc3_s": proton_processing_rate_mpc3_s,
        "T_sink_K": T_sink,
        "entropy_rate_J_K_mpc3_s": entropy_rate_J_K_mpc3_s,
        "entropy_rate_kB_mpc3_s": entropy_rate_kB_mpc3_s,
    })

    return df


# ============================================================
# 5. Plotting
# ============================================================

def make_plots(df, prefix="madau_dickinson"):
    z = df["z"].values

    plt.figure(figsize=(8, 5))
    plt.plot(z, df["sfrd_msun_yr_mpc3"])
    plt.gca().invert_xaxis()
    plt.xlabel("Redshift z")
    plt.ylabel(r"SFRD [$M_\odot$ yr$^{-1}$ Mpc$^{-3}$]")
    plt.title("Madau-Dickinson cosmic star-formation-rate density")
    plt.tight_layout()
    plt.savefig(f"{prefix}_sfrd.png", dpi=180)

    plt.figure(figsize=(8, 5))
    plt.plot(z, df["rho_star_msun_mpc3"])
    plt.gca().invert_xaxis()
    plt.xlabel("Redshift z")
    plt.ylabel(r"Accumulated stellar mass density [$M_\odot$ Mpc$^{-3}$]")
    plt.title("Integrated stellar mass density proxy")
    plt.tight_layout()
    plt.savefig(f"{prefix}_rho_star.png", dpi=180)

    plt.figure(figsize=(8, 5))
    plt.semilogy(z, df["fusion_chain_rate_mpc3_s"])
    plt.gca().invert_xaxis()
    plt.xlabel("Redshift z")
    plt.ylabel(r"Fusion chains s$^{-1}$ Mpc$^{-3}$")
    plt.title("Estimated stellar fusion-chain throughput density")
    plt.tight_layout()
    plt.savefig(f"{prefix}_fusion_rate.png", dpi=180)

    plt.figure(figsize=(8, 5))
    plt.semilogy(z, df["entropy_rate_kB_mpc3_s"])
    plt.gca().invert_xaxis()
    plt.xlabel("Redshift z")
    plt.ylabel(r"Entropy production [$k_B$ s$^{-1}$ Mpc$^{-3}$]")
    plt.title("Estimated stellar entropy-production rate density")
    plt.tight_layout()
    plt.savefig(f"{prefix}_entropy_rate.png", dpi=180)

    plt.close("all")


# ============================================================
# 6. Main
# ============================================================

def main():
    df = compute_history(z_max=20.0, n=6000)

    output_csv = "madau_dickinson_stellar_entropy.csv"
    df.to_csv(output_csv, index=False)

    make_plots(df)

    # Useful diagnostics
    peak_idx = df["sfrd_msun_yr_mpc3"].idxmax()
    today = df.iloc[-1]

    print("\nMadau-Dickinson stellar-registration integration")
    print("=" * 72)
    print(f"Peak SFRD redshift: z = {df.loc[peak_idx, 'z']:.3f}")
    print(f"Peak SFRD:          {df.loc[peak_idx, 'sfrd_msun_yr_mpc3']:.6e} Msun yr^-1 Mpc^-3")
    print()
    print("Today, z = 0:")
    print(f"  Stellar mass density proxy:      {today['rho_star_msun_mpc3']:.6e} Msun Mpc^-3")
    print(f"  Luminosity density proxy:        {today['lum_density_lsun_mpc3']:.6e} Lsun Mpc^-3")
    print(f"  Fusion chain rate density:       {today['fusion_chain_rate_mpc3_s']:.6e} chains s^-1 Mpc^-3")
    print(f"  Proton processing rate density:  {today['proton_processing_rate_mpc3_s']:.6e} protons s^-1 Mpc^-3")
    print(f"  Entropy production density:      {today['entropy_rate_kB_mpc3_s']:.6e} kB s^-1 Mpc^-3")
    print()
    print(f"Wrote CSV: {output_csv}")
    print("Wrote plots:")
    print("  madau_dickinson_sfrd.png")
    print("  madau_dickinson_rho_star.png")
    print("  madau_dickinson_fusion_rate.png")
    print("  madau_dickinson_entropy_rate.png")


if __name__ == "__main__":
    main()