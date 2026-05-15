# cmbr_universe_backend.py
import numpy as np
from cmbr import V_sigma, dV_dsigma, flux_scalar_system, H_of_N  # adjust imports to match your file

def integrate_background(N_start=-10.0, N_end=1.5, n_steps=4000):
    """
    Integrate flux cosmology and return t, a(t), H(t), sigma(t), Omegas etc.
    N = ln a; we map to a 'fake' time using t(N) = ∫ dN / H(N) for visualization.
    """
    from scipy.integrate import solve_ivp

    sigma_init = 5.0
    v_init = 0.0
    y0 = [sigma_init, v_init]

    sol = solve_ivp(
        flux_scalar_system,
        t_span=(N_start, N_end),
        y0=y0,
        dense_output=True,
        atol=1e-8,
        rtol=1e-8,
    )

    N_vals = np.linspace(N_start, N_end, n_steps)
    sigma_vals, v_vals = sol.sol(N_vals)
    a_vals = np.exp(N_vals)
    H_vals = H_of_N(N_vals, sigma_vals, v_vals)

    # Build an effective cosmic time t(N) just for animation pacing:
    t_vals = np.zeros_like(N_vals)
    for i in range(1, len(N_vals)):
        dN = N_vals[i] - N_vals[i-1]
        H_mid = 0.5 * (H_vals[i] + H_vals[i-1])
        t_vals[i] = t_vals[i-1] + dN / max(H_mid, 1e-8)

    # flux scalar energy density
    rho_r = 9e-5 * a_vals**-4
    rho_m = 0.3 * a_vals**-3
    rho_sigma = 0.5 * v_vals**2 + V_sigma(sigma_vals)
    rho_tot = rho_r + rho_m + rho_sigma
    Omega_sigma = rho_sigma / rho_tot

    # G_eff template (you can later import g1 from your lunar file if you want exact match)
    g1 = -1e-6
    flux_frac = Omega_sigma  # as a simple proxy here
    G_eff = 1.0 * (1.0 + g1 * flux_frac)

    return {
        "t": t_vals,
        "N": N_vals,
        "a": a_vals,
        "H": H_vals,
        "sigma": sigma_vals,
        "flux_frac": flux_frac,
        "G_eff": G_eff,
    }
