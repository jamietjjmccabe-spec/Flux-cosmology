import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ======================================
# 1. Cosmological parameters
# ======================================
Omega_m0 = 0.3
Omega_r0 = 9e-5
Omega_b0 = 0.05

# ======================================
# 2. MIP / coherence-field parameters
# ======================================
sigma0 = 0.0       # background coherence value
m_sigma = 0.2      # coherence restoration scale
lam = 0.0          # matter-coherence coupling (start with 0 for GR-safe tests)

def V_sigma(sigma):
    """Quadratic restoring potential."""
    return 0.5 * m_sigma**2 * (sigma - sigma0)**2

def dV_dsigma(sigma):
    """Derivative of quadratic restoring potential."""
    return m_sigma**2 * (sigma - sigma0)

# ======================================
# 3. Background densities and pressure
# ======================================
def background_densities(a):
    """Matter and radiation only."""
    rho_m = Omega_m0 * a**-3
    rho_r = Omega_r0 * a**-4
    p_m = 0.0
    p_r = rho_r / 3.0
    return rho_m, rho_r, p_m, p_r

def sigma_density_pressure(v, sigma):
    """Homogeneous scalar-field energy density and pressure."""
    rho_sigma = 0.5 * v**2 + V_sigma(sigma)
    p_sigma = 0.5 * v**2 - V_sigma(sigma)
    return rho_sigma, p_sigma

# ======================================
# 4. Hubble function
# ======================================
def H_of_N(N, sigma, v):
    a = np.exp(N)

    rho_m, rho_r, _, _ = background_densities(a)
    rho_sigma, _ = sigma_density_pressure(v, sigma)

    rho_tot = rho_m + rho_r + rho_sigma
    return np.sqrt(rho_tot)

# ======================================
# 5. ODE system in N = ln a
# ======================================
def flux_scalar_system(N, y):
    """
    y = [sigma, v], where
      sigma = coherence scalar
      v     = d sigma / dt

    Evolution equations:
      d sigma / dN = v / H

      FLRW-corrected MIP equation:
      (1 - 2 λ ρ_bg) * dv/dt + 3H (1 + 2 λ p_bg) * v + dV/dsigma = 0

    so
      dv/dt = - [3H (1 + 2 λ p_bg) v + dV/dsigma] / (1 - 2 λ ρ_bg)

      dv/dN = (dv/dt) / H
    """
    sigma, v = y
    a = np.exp(N)

    rho_m, rho_r, p_m, p_r = background_densities(a)
    rho_bg = rho_m + rho_r
    p_bg = p_m + p_r

    H = H_of_N(N, sigma, v)

    inertia = 1.0 - 2.0 * lam * rho_bg
    if np.any(np.abs(inertia) < 1e-10):
        raise RuntimeError(
            f"Encountered near-singular effective inertia: 1 - 2 λ ρ_bg = {inertia}"
        )

    dsigma_dN = v / H

    dv_dt = -(
        3.0 * H * (1.0 + 2.0 * lam * p_bg) * v
        + dV_dsigma(sigma)
    ) / inertia

    dv_dN = dv_dt / H
    return [dsigma_dN, dv_dN]

# ======================================
# 6. Integrate from early universe to today
# ======================================
N_start = -10.0
N_end = 0.0

# Choose a displaced initial sigma if you want dynamics.
# For GR/ΛCDM recovery, start sigma_init = sigma0 and v_init = 0.
sigma_init = 0.3
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

N_vals = np.linspace(N_start, N_end, 2000)
sigma_vals, v_vals = sol.sol(N_vals)

a_vals = np.exp(N_vals)
z_vals = 1.0 / a_vals - 1.0

H_vals = H_of_N(N_vals, sigma_vals, v_vals)

rho_m_vals = Omega_m0 * a_vals**-3
rho_r_vals = Omega_r0 * a_vals**-4
p_r_vals = rho_r_vals / 3.0

rho_sigma_vals = 0.5 * v_vals**2 + V_sigma(sigma_vals)
p_sigma_vals = 0.5 * v_vals**2 - V_sigma(sigma_vals)

rho_tot_vals = rho_m_vals + rho_r_vals + rho_sigma_vals
p_tot_vals = p_r_vals + p_sigma_vals

# ======================================
# 7. CMB-relevant integrals
# ======================================
z_star = 1100.0
a_star = 1.0 / (1.0 + z_star)

mask_star = a_vals <= a_star
a_early = a_vals[mask_star]
H_early = H_vals[mask_star]

rho_b_early = Omega_b0 * a_early**-3
rho_gamma_early = Omega_r0 * a_early**-4

cs2 = 1.0 / (3.0 * (1.0 + 3.0 * rho_b_early / (4.0 * rho_gamma_early)))
cs = np.sqrt(cs2)

integrand_rs = cs / (a_early**2 * H_early)
r_s = np.trapz(integrand_rs, a_early)

mask_los = a_vals >= a_star
a_los = a_vals[mask_los]
H_los = H_vals[mask_los]

integrand_chi = 1.0 / (a_los**2 * H_los)
chi_star = np.trapz(integrand_chi, a_los)

print("Sound horizon r_s (in H0^-1 units):", r_s)
print("Comoving distance to last scattering χ_*:", chi_star)

# ======================================
# 8. Derived quantities
# ======================================
w_sigma_vals = p_sigma_vals / np.maximum(rho_sigma_vals, 1e-30)
Omega_sigma_vals = rho_sigma_vals / rho_tot_vals

print("w_sigma(z=0):", w_sigma_vals[-1])
print("Omega_sigma(z=0):", Omega_sigma_vals[-1])

# ======================================
# 9. Plots
# ======================================
plt.figure()
plt.loglog(a_vals, rho_r_vals, label="Radiation")
plt.loglog(a_vals, rho_m_vals, label="Matter")
plt.loglog(a_vals, rho_sigma_vals, label="Coherence field (σ)")
plt.xlabel("Scale factor a")
plt.ylabel("Energy density")
plt.legend()
plt.title("Background densities")

plt.figure()
plt.loglog(a_vals, H_vals)
plt.xlabel("Scale factor a")
plt.ylabel("H(a)")
plt.title("Hubble rate")

plt.figure()
plt.plot(a_vals, sigma_vals)
plt.xlabel("Scale factor a")
plt.ylabel("σ(a)")
plt.title("Coherence field evolution")

plt.tight_layout()
plt.show()