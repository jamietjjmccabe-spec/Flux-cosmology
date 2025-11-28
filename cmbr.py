import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ================================
# 1. Cosmological parameters
# ================================
# Dimensionless units with 8πG/3 = 1, H0 ~ O(1)
Omega_m0 = 0.3    # total matter today
Omega_r0 = 9e-5   # radiation today (photons+neutrinos approx)
Omega_b0 = 0.05   # baryons subset of matter

# ================================
# 2. Flux scalar potential
# ================================
# sigma = your flux/coherence scalar
# V(sigma) = "flux tension" / residual flux energy

V0 = 0.7   # overall scale of flux energy
mu = 2.0   # controls how fast flux exhausts as sigma changes

def V_sigma(sigma):
    """
    Flux potential:
    - For large sigma: approaches plateau ~ V0 (quasi-Λ, accelerated expansion).
    - Near sigma ~ 0: ~quadratic minimum.
    """
    return V0 * (1.0 - np.exp(-sigma/mu))**2

def dV_dsigma(sigma):
    """Derivative dV/dsigma."""
    exp_term = np.exp(-sigma/mu)
    return 2.0 * V0 * (1.0 - exp_term) * exp_term / mu

# ================================
# 3. Hubble parameter H(N)
# ================================
def H_of_N(N, sigma, v):
    """
    Hubble parameter H(N) given:
      N     = ln a
      sigma = flux scalar
      v     = d sigma / dt  (time derivative)
    Energy densities (dimensionless):
      rho_r     ~ Omega_r0 a^-4
      rho_m     ~ Omega_m0 a^-3
      rho_sigma = 0.5 v^2 + V(sigma)
    We set 8πG/3 = 1, so H^2 = rho_total.
    """
    a = np.exp(N)
    rho_r = Omega_r0 * a**-4
    rho_m = Omega_m0 * a**-3
    rho_sigma = 0.5 * v**2 + V_sigma(sigma)
    H2 = rho_r + rho_m + rho_sigma
    return np.sqrt(H2)

# ================================
# 4. ODE system for flux scalar
# ================================
def flux_scalar_system(N, y):
    """
    Evolve the flux scalar in terms of N = ln(a).
    
    y = [sigma, v], where:
      sigma(N) = flux scalar
      v(N)     = d sigma / dt

    Using:
      dN/dt = H  => d/dN = (1/H) d/dt

    Then:
      d sigma / dN = (d sigma / dt) / H = v / H
      d v / dN     = (1/H) d v / dt
                   = (1/H) [ -3 H v - dV/dsigma ]
                   = -3 v - (1/H) dV/dsigma
    """
    sigma, v = y
    H = H_of_N(N, sigma, v)
    dsigma_dN = v / H
    dv_dN = -3.0 * v - dV_dsigma(sigma) / H
    return [dsigma_dN, dv_dN]

# ================================
# 5. Integrate from early universe to today
# ================================
N_start = -10.0  # early: a ~ e^-10 ~ 4.5e-5
N_end   =  0.0   # today: a = 1

# Initial conditions for your scalar at early times:
sigma_init = 5.0   # high on plateau => almost constant early
v_init     = 0.0   # frozen / slow-roll start

y0 = [sigma_init, v_init]

sol = solve_ivp(
    flux_scalar_system,
    t_span=(N_start, N_end),
    y0=y0,
    dense_output=True,
    atol=1e-8,
    rtol=1e-8,
)

# Sample solution on a grid
N_vals = np.linspace(N_start, N_end, 2000)
sigma_vals, v_vals = sol.sol(N_vals)

a_vals = np.exp(N_vals)
z_vals = 1.0 / a_vals - 1.0

H_vals = H_of_N(N_vals, sigma_vals, v_vals)

rho_r_vals = Omega_r0 * a_vals**-4
rho_m_vals = Omega_m0 * a_vals**-3
rho_sigma_vals = 0.5 * v_vals**2 + V_sigma(sigma_vals)
rho_tot_vals = rho_r_vals + rho_m_vals + rho_sigma_vals

# ================================
# 6. CMB-relevant integrals
# ================================
# Recombination redshift (approx Planck value)
z_star = 1100.0
a_star = 1.0 / (1.0 + z_star)

# --- Sound horizon r_s(a_*) ---
mask_star = a_vals <= a_star
a_early = a_vals[mask_star]
H_early = H_vals[mask_star]

rho_b_early = Omega_b0 * a_early**-3
rho_gamma_early = Omega_r0 * a_early**-4  # approx: photons ≈ total rad early

# Photon-baryon sound speed c_s^2 = 1 / [3(1 + 3ρ_b / 4ρ_γ)]
cs2 = 1.0 / (3.0 * (1.0 + 3.0 * rho_b_early / (4.0 * rho_gamma_early)))
cs = np.sqrt(cs2)

# r_s = ∫_0^{a_*} c_s / (a'^2 H(a')) da'
integrand_rs = cs / (a_early**2 * H_early)
r_s = np.trapz(integrand_rs, a_early)

# --- Comoving distance to last scattering χ_* ---
mask_los = a_vals >= a_star
a_los = a_vals[mask_los]
H_los = H_vals[mask_los]

# χ_* = ∫_{a_*}^1 da / (a^2 H(a))
integrand_chi = 1.0 / (a_los**2 * H_los)
chi_star = np.trapz(integrand_chi, a_los)

print("Sound horizon r_s (in H0^-1 units)      :", r_s)
print("Comoving distance to last scattering χ_*:", chi_star)

# ================================
# 7. Some background plots
# ================================

# Energy densities vs scale factor
plt.figure()
plt.loglog(a_vals, rho_r_vals, label="Radiation")
plt.loglog(a_vals, rho_m_vals, label="Matter")
plt.loglog(a_vals, rho_sigma_vals, label="Flux scalar (σ)")
plt.xlabel("Scale factor a")
plt.ylabel("Energy density (arb. units)")
plt.legend()
plt.title("Background densities with flux scalar")

# Hubble vs scale factor
plt.figure()
plt.loglog(a_vals, H_vals)
plt.xlabel("Scale factor a")
plt.ylabel("H(a)")
plt.title("Hubble rate with flux-driven dark energy")

# Flux scalar evolution
plt.figure()
plt.plot(a_vals, sigma_vals)
plt.xlabel("Scale factor a")
plt.ylabel("σ(a)  (flux scalar)")
plt.title("Evolution of flux/coherence scalar σ")
# ======================================
# 8. Derived quantities: w_sigma(a), Omega_sigma(a)
# ======================================
# In our units 8πG/3 = 1, so:
# H^2 = rho_total
# rho_sigma = 0.5 v^2 + V(sigma)
# p_sigma   = 0.5 v^2 - V(sigma)

rho_sigma_vals = 0.5 * v_vals**2 + V_sigma(sigma_vals)
p_sigma_vals   = 0.5 * v_vals**2 - V_sigma(sigma_vals)

w_sigma_vals = p_sigma_vals / rho_sigma_vals

# total density we already had as rho_tot_vals
Omega_sigma_vals = rho_sigma_vals / rho_tot_vals

# Today's values (a = 1, N = 0)
w_sigma_0 = w_sigma_vals[-1]
Omega_sigma_0 = Omega_sigma_vals[-1]

print("w_sigma(z=0)        :", w_sigma_0)
print("Omega_sigma(z=0)    :", Omega_sigma_0)


plt.tight_layout()
plt.show()
