"""
clock_drift_hubble_scan.py

Flux Cosmology / MIP clock-drift test for the Hubble tension.

Purpose
-------
This file keeps the CMB geometry from cmbr.py fixed, especially:

    theta_* = r_s(a_*) / chi_*

Then it adds a separate environmental clock-rate law for MIP/actualization.
The aim is to test whether the same background can yield:

    H0_CMB   ~= 67.4 km/s/Mpc    from the high-sigma CMB calibration
    H0_local ~= 73.0 km/s/Mpc    from low-sigma local distance-ladder clocks

without moving the CMB acoustic angle theta_*.

Interpretation
--------------
Planck time remains the hard floor. The zeptosecond/MIP interval is treated as
an emergent coarse-grained actualization cadence. In a low-sigma environment,
the MIP cadence can run at a slightly different rate relative to the high-sigma
CMB reference environment.

This script does NOT claim the model is correct. It asks a clean numerical
question:

    How large must the environmental clock drift be to map Planck-like H0 to
    SH0ES-like H0 while leaving theta_* unchanged?

If the required drift is small, smooth, and observationally testable, it becomes
worth adding to the main solver. If it requires huge or pathological drift, the
idea fails cleanly.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ============================================================
# 1. Baseline cosmology parameters from cmbr.py
# ============================================================
# Dimensionless units with 8*pi*G/3 = 1.
Omega_m0 = 0.3
Omega_r0 = 9e-5
Omega_b0 = 0.05

V0 = 0.7
mu = 2.0

# Reference H0 values in km/s/Mpc
H0_CMB_REF = 67.4
H0_LOCAL_TARGET = 73.0

# Recombination redshift
z_star = 1100.0

# Integration range in N = ln(a)
N_start = -10.0
N_end = 0.0
N_samples = 4000

# Initial flux/coherence scalar
sigma_init = 5.0
v_init = 0.0


# ============================================================
# 2. Flux scalar background
# ============================================================
def V_sigma(sigma):
    """Flux scalar potential."""
    return V0 * (1.0 - np.exp(-sigma / mu)) ** 2


def dV_dsigma(sigma):
    """Derivative of the flux scalar potential."""
    exp_term = np.exp(-sigma / mu)
    return 2.0 * V0 * (1.0 - exp_term) * exp_term / mu


def H_of_N(N, sigma, v):
    """Dimensionless geometric Hubble rate H(N)."""
    a = np.exp(N)
    rho_r = Omega_r0 * a ** -4
    rho_m = Omega_m0 * a ** -3
    rho_sigma = 0.5 * v ** 2 + V_sigma(sigma)
    return np.sqrt(rho_r + rho_m + rho_sigma)


def flux_scalar_system(N, y):
    """ODE system for sigma(N), v(N), where v = d sigma / dt."""
    sigma, v = y
    H = H_of_N(N, sigma, v)
    return [v / H, -3.0 * v - dV_dsigma(sigma) / H]


def solve_background():
    """Solve the baseline cmbr.py-style background."""
    sol = solve_ivp(
        flux_scalar_system,
        t_span=(N_start, N_end),
        y0=[sigma_init, v_init],
        dense_output=True,
        atol=1e-9,
        rtol=1e-9,
    )

    N_vals = np.linspace(N_start, N_end, N_samples)
    sigma_vals, v_vals = sol.sol(N_vals)
    a_vals = np.exp(N_vals)
    z_vals = 1.0 / a_vals - 1.0
    H_vals = H_of_N(N_vals, sigma_vals, v_vals)

    return {
        "N": N_vals,
        "a": a_vals,
        "z": z_vals,
        "sigma": sigma_vals,
        "v": v_vals,
        "H": H_vals,
    }


# ============================================================
# 3. CMB geometry: sound horizon and last-scattering distance
# ============================================================
def cmb_geometry(bg):
    """Compute r_s, chi_*, and theta_* from the baseline geometry."""
    a_vals = bg["a"]
    H_vals = bg["H"]

    a_star = 1.0 / (1.0 + z_star)

    # Sound horizon: r_s = int_0^a* c_s da / (a^2 H)
    mask_star = a_vals <= a_star
    a_early = a_vals[mask_star]
    H_early = H_vals[mask_star]

    rho_b_early = Omega_b0 * a_early ** -3
    rho_gamma_early = Omega_r0 * a_early ** -4

    cs2 = 1.0 / (3.0 * (1.0 + 3.0 * rho_b_early / (4.0 * rho_gamma_early)))
    cs = np.sqrt(cs2)

    r_s = np.trapezoid(cs / (a_early ** 2 * H_early), a_early)

    # Comoving distance: chi_* = int_a*^1 da / (a^2 H)
    mask_los = a_vals >= a_star
    a_los = a_vals[mask_los]
    H_los = H_vals[mask_los]
    chi_star = np.trapezoid(1.0 / (a_los ** 2 * H_los), a_los)

    theta_star = r_s / chi_star

    return {
        "a_star": a_star,
        "r_s": r_s,
        "chi_star": chi_star,
        "theta_star": theta_star,
        "ell_A_proxy": np.pi / theta_star,
    }


# ============================================================
# 4. Environmental MIP clock-drift law
# ============================================================
def mip_clock_rate_factor(sigma_env_norm, beta, power=1.0):
    """
    Environmental MIP clock-rate factor.

    sigma_env_norm:
        Normalized environmental coherence, 0 <= sigma <= 1.
        sigma = 1 means high-coherence/high-sigma reference environment.
        sigma = 0 means maximally low-coherence environment.

    beta:
        Maximum fractional clock-rate lift in the lowest-sigma environment.

    power:
        Shape parameter. power=1 is linear; power>1 makes the drift weaker
        until very low sigma.

    Returns
    -------
    rate_factor:
        Multiplicative factor on the inferred local H0.

        H0_inferred(env) = H0_CMB_reference * rate_factor

    MIP interval factor:
        tau_env / tau_ref = 1 / rate_factor.
        So a rate_factor > 1 means shorter coarse-grained MIP interval.
    """
    sigma_env_norm = np.clip(sigma_env_norm, 0.0, 1.0)
    deficit = 1.0 - sigma_env_norm
    return 1.0 + beta * deficit ** power


def inferred_H0(H0_cmb, sigma_env_norm, beta, power=1.0):
    """Observed H0 inferred under the environmental clock-rate law."""
    return H0_cmb * mip_clock_rate_factor(sigma_env_norm, beta, power)


def beta_required_for_target(H0_cmb, H0_target, sigma_env_norm, power=1.0):
    """Analytic beta needed to map H0_cmb to H0_target for a given environment."""
    sigma_env_norm = np.clip(sigma_env_norm, 0.0, 1.0)
    deficit = 1.0 - sigma_env_norm
    if deficit <= 0.0:
        return np.inf
    return (H0_target / H0_cmb - 1.0) / (deficit ** power)


# ============================================================
# 5. Scan and diagnostics
# ============================================================
def run_scan():
    bg = solve_background()
    geom = cmb_geometry(bg)

    # Local environments to test.
    # 1.0 = CMB reference; lower values = lower local coherence.
    sigma_env_grid = np.linspace(1.0, 0.05, 200)
    power_values = [1.0, 1.5, 2.0]

    rows = []
    for power in power_values:
        beta_req = np.array([
            beta_required_for_target(H0_CMB_REF, H0_LOCAL_TARGET, s, power)
            for s in sigma_env_grid
        ])
        rows.append((power, beta_req))

    # Fiducial example: if the local ladder samples a strongly low-sigma
    # environment, sigma_env_norm = 0 means the needed beta is just the H0 ratio - 1.
    sigma_local_fid = 0.0
    beta_fid = beta_required_for_target(
        H0_CMB_REF, H0_LOCAL_TARGET, sigma_local_fid, power=1.0
    )
    H0_local_fid = inferred_H0(H0_CMB_REF, sigma_local_fid, beta_fid, power=1.0)
    tau_interval_factor = 1.0 / mip_clock_rate_factor(sigma_local_fid, beta_fid, power=1.0)

    print("================ Flux Clock-Drift Hubble Scan ================")
    print(f"Baseline H0_CMB reference        : {H0_CMB_REF:.3f} km/s/Mpc")
    print(f"Target local H0                 : {H0_LOCAL_TARGET:.3f} km/s/Mpc")
    print(f"Required full-deficit beta      : {beta_fid:.6f}")
    print(f"Clock-rate lift needed          : {(H0_LOCAL_TARGET/H0_CMB_REF - 1.0)*100:.3f}%")
    print(f"MIP interval factor tau/tau_ref : {tau_interval_factor:.6f}")
    print(f"Recovered H0_local              : {H0_local_fid:.3f} km/s/Mpc")
    print()
    print("CMB geometry is unchanged by this first-order clock law:")
    print(f"r_s                             : {geom['r_s']:.8e} H0^-1")
    print(f"chi_*                           : {geom['chi_star']:.8e} H0^-1")
    print(f"theta_* = r_s/chi_*             : {geom['theta_star']:.8e}")
    print(f"ell_A proxy = pi/theta_*        : {geom['ell_A_proxy']:.3f}")
    print()
    print("Interpretation:")
    print("A beta of ~0.083 means the local actualization clock would need to")
    print("run about 8.3% faster than the high-sigma CMB reference clock, or")
    print("equivalently the coarse MIP interval would be about 7.7% shorter.")

    # Plot 1: beta required versus assumed local sigma.
    plt.figure(figsize=(8, 5))
    for power, beta_req in rows:
        plt.plot(sigma_env_grid, beta_req, label=f"power={power}")
    plt.axhline(beta_fid, linestyle="--", label="full-deficit beta")
    plt.xlabel("local normalized sigma environment")
    plt.ylabel("beta required to reach H0 = 73")
    plt.title("Clock-drift strength required by local sigma environment")
    plt.legend()
    plt.tight_layout()
    plt.savefig("clock_drift_beta_required.png", dpi=180)

    # Plot 2: inferred H0 for a few beta values.
    plt.figure(figsize=(8, 5))
    for beta in [0.02, 0.05, beta_fid, 0.12]:
        H0_curve = inferred_H0(H0_CMB_REF, sigma_env_grid, beta, power=1.0)
        plt.plot(sigma_env_grid, H0_curve, label=f"beta={beta:.3f}")
    plt.axhline(H0_CMB_REF, linestyle="--", label="CMB reference")
    plt.axhline(H0_LOCAL_TARGET, linestyle="--", label="local target")
    plt.xlabel("local normalized sigma environment")
    plt.ylabel("inferred H0 [km/s/Mpc]")
    plt.title("Apparent H0 under environmental MIP clock drift")
    plt.legend()
    plt.tight_layout()
    plt.savefig("clock_drift_h0_scan.png", dpi=180)

    return bg, geom


if __name__ == "__main__":
    run_scan()
