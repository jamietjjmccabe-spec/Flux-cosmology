import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 0. UNITS & RADIAL GRID
# ============================================================
# Length unit: 1 = 1 kpc
# Velocity unit: 1 = 220 km/s
G = 1.0

r_max = 25.0
n_r   = 400
r = np.linspace(0.1, r_max, n_r)  # avoid r=0

# ============================================================
# 1. INITIAL FLUX PROFILE φ(r, t=0)
# ============================================================
phi0    = 1.0   # central normalisation
R_phi   = 2.5   # kpc  (inner CMZ/bar)
phi_halo = 0.3  # outer floor
R_h      = 8.0  # kpc  (halo transition)

phi_init = phi0 * np.exp(-r / R_phi) + phi_halo * (1.0 - np.exp(-r / R_h))

# ============================================================
# 2. STELLAR BULGE + DISK (STARS) -- mostly static
# ============================================================
M_bulge_total = 1.0
R_b = 1.0

Sigma_b = (M_bulge_total / (2.0 * np.pi * R_b**2)) * np.exp(-r / R_b)

M_disk_star_initial = 3.0   # initial stellar disk mass (smaller than final)
R_d = 3.0

Sigma_d_star = (M_disk_star_initial / (2.0 * np.pi * R_d**2)) * np.exp(-r / R_d)

# ============================================================
# 3. GAS DISK (EVOLVES) + INITIAL SFR
# ============================================================
M_gas_total = 3.0     # initially gas-rich
R_g = 4.5

Sigma_g = (M_gas_total / (2.0 * np.pi * R_g**2)) * np.exp(-r / R_g)

# add a Gaussian gas ring (like MW molecular ring)
ring_radius = 8.0
ring_width  = 1.5
ring_boost  = np.exp(-0.5 * ((r - ring_radius) / ring_width)**2)
Sigma_g *= (1.0 + 1.0 * ring_boost)

# ============================================================
# 4. FLUX-HALO DENSITY FROM φ
# ============================================================
rho_h0     = 0.05   # overall strength of flux-halo
alpha_flux = 1.5    # extra gravity factor for flux-mass

def rho_flux_from_phi(phi_r):
    """Spherical flux-halo density profile at given φ(r)."""
    return rho_h0 * (phi_r / phi0) / (1.0 + (r / R_h)**2)

def cumulative_mass_from_density(r_grid, rho_grid):
    """M(<r) = ∫ 4π r'^2 rho(r') dr'."""
    M = np.zeros_like(r_grid)
    integrand = 4.0 * np.pi * (r_grid**2) * rho_grid
    for i in range(1, len(r_grid)):
        M[i] = M[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * (r_grid[i] - r_grid[i-1])
    return M

def cumulative_mass_from_surface_density(r_grid, Sigma):
    """M(<r) ≈ ∫ 2π r' Σ(r') dr' (thin axisymmetric disk)."""
    M = np.zeros_like(r_grid)
    integrand = 2.0 * np.pi * r_grid * Sigma
    for i in range(1, len(r_grid)):
        M[i] = M[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * (r_grid[i] - r_grid[i-1])
    return M

# ============================================================
# 5. TIME EVOLUTION PARAMETERS
# ============================================================
t_max = 10.0          # "Gyr" in toy units
Nt    = 200           # time steps
t = np.linspace(0.0, t_max, Nt)
dt = t[1] - t[0]

# Couplings
epsilon_SF  = 0.3     # efficiency of SFR from gas+flux
gamma_star  = 0.2     # how strongly SFR depletes φ
gamma_in    = 0.02    # inflow / fresh flux source
gamma_relax = 0.05    # relaxation to floor flux
phi_floor   = 0.15    # minimal flux level

# Arrays to store evolution
phi      = np.zeros((Nt, n_r))
Sigma_g_t = np.zeros((Nt, n_r))
Sigma_d_star_t = np.zeros((Nt, n_r))
SFR_t    = np.zeros((Nt, n_r))

phi[0, :]           = phi_init
Sigma_g_t[0, :]     = Sigma_g
Sigma_d_star_t[0,:] = Sigma_d_star

# radial profile where inflow is strongest (inside-out formation)
inflow_profile = np.exp(-r / 5.0)

# ============================================================
# 6. MAIN EVOLUTION LOOP
# ============================================================
for n in range(Nt - 1):
    phi_n   = phi[n, :]
    gas_n   = Sigma_g_t[n, :]
    stars_n = Sigma_d_star_t[n, :]

    # Star formation rate surface density
    # SFR ∝ ε * φ * gas
    SFR = epsilon_SF * phi_n * gas_n
    SFR_t[n, :] = SFR

    # Evolution of gas: consumed by SFR
    dSigma_g = -SFR

    # Stellar disk grows from SFR
    dSigma_star = SFR

    # Evolution of flux:
    #  - depleted where SFR is high
    #  - replenished by inflow (more at small r)
    #  - relaxes toward a floor value
    dphi = (
        -gamma_star * SFR
        + gamma_in * inflow_profile
        - gamma_relax * (phi_n - phi_floor)
    )

    # Euler update
    gas_next   = gas_n   + dt * dSigma_g
    star_next  = stars_n + dt * dSigma_star
    phi_next   = phi_n   + dt * dphi

    # Physical bounds
    gas_next  = np.clip(gas_next, 0.0, None)
    star_next = np.clip(star_next, 0.0, None)
    phi_next  = np.clip(phi_next, 0.01, 2.0)

    Sigma_g_t[n+1, :]     = gas_next
    Sigma_d_star_t[n+1, :] = star_next
    phi[n+1, :]           = phi_next

# compute SFR for last step
SFR_t[-1,:] = epsilon_SF * phi[-1,:] * Sigma_g_t[-1,:]

# ============================================================
# 7. ROTATION CURVES AT SELECTED TIMES
# ============================================================
snap_indices = [0, Nt//3, 2*Nt//3, Nt-1]
snap_labels  = [f"t = {t[i]:.1f}" for i in snap_indices]

plt.figure(figsize=(7,5))

for i, idx in enumerate(snap_indices):
    phi_snap   = phi[idx, :]
    Sigma_g_snap   = Sigma_g_t[idx, :]
    Sigma_star_snap = Sigma_b + Sigma_d_star_t[idx, :]  # bulge + disk stars

    # Masses
    M_b = cumulative_mass_from_surface_density(r, Sigma_b)
    M_star_disk = cumulative_mass_from_surface_density(r, Sigma_d_star_t[idx,:])
    M_star_total = M_b + M_star_disk

    rho_flux = rho_flux_from_phi(phi_snap)
    M_flux   = cumulative_mass_from_density(r, rho_flux)

    v_star2 = G * M_star_total / r
    v_flux2 = alpha_flux * G * M_flux / r
    v_tot   = np.sqrt(v_star2 + v_flux2)

    plt.plot(r, v_tot, label=snap_labels[i])

plt.xlabel("r [kpc]")
plt.ylabel("v_c(r) [v_unit ≈ 220 km/s]")
plt.title("Evolution of rotation curve in flux Milky Way")
plt.legend()
plt.grid(alpha=0.3)

# ============================================================
# 8. SFR & GAS EVOLUTION
# ============================================================
plt.figure(figsize=(7,5))
for i, idx in enumerate(snap_indices):
    plt.plot(r, SFR_t[idx,:], label=snap_labels[i])
plt.xlabel("r [kpc]")
plt.ylabel("SFR surface density (arb.)")
plt.title("Flux-driven star formation over time")
plt.legend()
plt.grid(alpha=0.3)

plt.figure(figsize=(7,5))
for i, idx in enumerate(snap_indices):
    plt.plot(r, Sigma_g_t[idx,:], label=snap_labels[i])
plt.xlabel("r [kpc]")
plt.ylabel("Gas surface density Σ_g (arb.)")
plt.title("Gas depletion and inside-out quenching")
plt.legend()
plt.grid(alpha=0.3)

# ============================================================
# 9. FLUX PROFILE EVOLUTION
# ============================================================
plt.figure(figsize=(7,5))
for i, idx in enumerate(snap_indices):
    plt.plot(r, phi[idx,:], label=snap_labels[i])
plt.xlabel("r [kpc]")
plt.ylabel(r"$\phi(r,t)$")
plt.title("Flux fraction profile evolving in time")
plt.legend()
plt.grid(alpha=0.3)

# ============================================================
# 10. 2D FLUX FIELD AT EARLY & LATE TIMES
# ============================================================
n_xy = 300
x = np.linspace(-r_max, r_max, n_xy)
y = np.linspace(-r_max, r_max, n_xy)
X, Y = np.meshgrid(x, y)
R_xy = np.sqrt(X**2 + Y**2)

def phi_2d_from_snapshot(idx):
    phi_snap = phi[idx,:]
    phi_2d = np.interp(R_xy, r, phi_snap, left=phi_snap[0], right=phi_snap[-1])
    phi_2d[R_xy > r_max] = np.nan
    return phi_2d

phi2d_early = phi_2d_from_snapshot(snap_indices[0])
phi2d_late  = phi_2d_from_snapshot(snap_indices[-1])

plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.imshow(phi2d_early, origin="lower",
           extent=[x.min(), x.max(), y.min(), y.max()],
           aspect="equal")
plt.colorbar(label=r"$\phi(x,y)$")
plt.title("Flux field (early)")
plt.xlabel("x [kpc]")
plt.ylabel("y [kpc]")

plt.subplot(1,2,2)
plt.imshow(phi2d_late, origin="lower",
           extent=[x.min(), x.max(), y.min(), y.max()],
           aspect="equal")
plt.colorbar(label=r"$\phi(x,y)$")
plt.title("Flux field (late)")
plt.xlabel("x [kpc]")
plt.ylabel("y [kpc]")

plt.tight_layout()
plt.show()
