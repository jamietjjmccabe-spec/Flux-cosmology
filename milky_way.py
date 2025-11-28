import numpy as np
import matplotlib.pyplot as plt

# ============================================
# 1. UNITS & BASIC CONSTANTS
# ============================================

# Galaxy units:
#  - Length unit: 1 = 1 kpc
#  - Velocity unit: 1 = 220 km/s (roughly local MW rotation speed)
#  - G is chosen so that v^2 = G M(<r) / r in these units.

G = 1.0  # toy gravitational constant

# Radial grid (0 to 25 kpc)
r_max = 25.0
n_r   = 500
r = np.linspace(0.1, r_max, n_r)  # avoid r=0 exactly

# ============================================
# 2. FLUX PROFILE φ(r)
# ============================================

# Inner "bar/CMZ" flux scale
phi0    = 1.0   # normalisation
R_phi   = 2.5   # kpc, inner flux scale

# Outer flux-halo floor
phi_halo = 0.3  # background halo flux fraction
R_h      = 8.0  # kpc, halo transition scale

phi_r = phi0 * np.exp(-r / R_phi) + phi_halo * (1.0 - np.exp(-r / R_h))

# ============================================
# 3. MATTER PROFILES: BULGE + DISK
# ============================================

# Simple exponential bulge
M_bulge_total = 1.0   # arbitrary units
R_b           = 1.0   # kpc

Sigma_b = (M_bulge_total / (2.0 * np.pi * R_b**2)) * np.exp(-r / R_b)

# Exponential stellar disk
M_disk_total = 5.0    # arbitrary units
R_d          = 3.0    # kpc

Sigma_d = (M_disk_total / (2.0 * np.pi * R_d**2)) * np.exp(-r / R_d)

# ============================================
# 4. FLUX-HALO EFFECTIVE DENSITY
# ============================================

# Overall strength of the flux-halo
rho_h0      = 0.05     # sets how much "pseudo mass" sits in flux
alpha_flux  = 1.5      # scaling of flux-mass to gravity (like a fifth-force factor)

# We approximate the halo as spherical with density:
# rho_flux(r) ∝ (phi(r)/phi0) * [1 / (1 + (r/R_h)^2)]
rho_flux = rho_h0 * (phi_r / phi0) / (1.0 + (r / R_h)**2)

# To get M_flux(<r), we integrate 4π r'^2 rho_flux(r') dr'
# We'll do this numerically.
def cumulative_mass_from_density(r_grid, rho_grid):
    """
    Spherical cumulative mass: M(<r) = ∫ 4π r'^2 rho(r') dr'
    Using cumulative trapezoidal integration.
    """
    M = np.zeros_like(r_grid)
    integrand = 4.0 * np.pi * (r_grid**2) * rho_grid
    for i in range(1, len(r_grid)):
        M[i] = M[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * (r_grid[i] - r_grid[i-1])
    return M

M_flux = cumulative_mass_from_density(r, rho_flux)

# ============================================
# 5. BULGE & DISK CUMULATIVE MASS
# ============================================

# For a thin axisymmetric disk, we approximate cumulative mass as:
# M(<r) ≈ ∫_0^r 2π r' Σ(r') dr'
def cumulative_mass_from_surface_density(r_grid, Sigma):
    M = np.zeros_like(r_grid)
    integrand = 2.0 * np.pi * r_grid * Sigma
    for i in range(1, len(r_grid)):
        M[i] = M[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * (r_grid[i] - r_grid[i-1])
    return M

M_b = cumulative_mass_from_surface_density(r, Sigma_b)
M_d = cumulative_mass_from_surface_density(r, Sigma_d)

# ============================================
# 6. ROTATION CURVE v_c(r)
# ============================================

# Each component v^2 = G M(<r) / r
v_b2    = G * M_b    / r
v_d2    = G * M_d    / r
v_flux2 = alpha_flux * G * M_flux / r  # flux-enhanced halo gravity

v_c2 = v_b2 + v_d2 + v_flux2
v_c  = np.sqrt(v_c2)

# ============================================
# 7. FLUX-DRIVEN STAR FORMATION PROFILE
# ============================================

# Simple gas disk: a bit more extended than stellar disk
M_gas_total = 1.5
R_g         = 4.5

Sigma_g = (M_gas_total / (2.0 * np.pi * R_g**2)) * np.exp(-r / R_g)

# Add a Gaussian "ring" at ~8 kpc (Solar circle style)
ring_radius = 8.0
ring_width  = 1.5

ring_boost = np.exp(-0.5 * ((r - ring_radius)/ring_width)**2)
Sigma_g *= (1.0 + 1.0 * ring_boost)

# Star formation rate surface density ~ phi(r) * gas surface density
SFR = phi_r * Sigma_g

# ============================================
# 8. 2D FLUX FIELD Φ(x, y) FOR VISUALIZATION
# ============================================

# Make a 2D grid in the galactic plane
n_xy = 300
x = np.linspace(-r_max, r_max, n_xy)
y = np.linspace(-r_max, r_max, n_xy)
X, Y = np.meshgrid(x, y)
R_xy = np.sqrt(X**2 + Y**2)

# Interpolate phi(r) onto R_xy
phi_2d = np.interp(R_xy, r, phi_r, left=phi_r[0], right=phi_r[-1])

# Optional: mask radii beyond r_max
phi_2d[R_xy > r_max] = np.nan

# ============================================
# 9. PLOTTING
# ============================================

plt.figure(figsize=(6,4))
plt.plot(r, phi_r)
plt.xlabel("r [kpc]")
plt.ylabel(r"$\phi(r)$ (flux fraction)")
plt.title("Flux fraction profile in the Milky Way toy model")
plt.grid(alpha=0.3)

plt.figure(figsize=(6,4))
plt.plot(r, Sigma_b, label="Bulge Σ_b")
plt.plot(r, Sigma_d, label="Disk Σ_d")
plt.plot(r, rho_flux, label="Flux-halo ρ_flux (arb.)")
plt.yscale("log")
plt.xlabel("r [kpc]")
plt.ylabel("Density / Surface density (arb.)")
plt.title("Matter & flux-halo components")
plt.legend()
plt.grid(alpha=0.3)

plt.figure(figsize=(6,4))
plt.plot(r, np.sqrt(v_b2),    label="Bulge")
plt.plot(r, np.sqrt(v_d2),    label="Disk")
plt.plot(r, np.sqrt(v_flux2), label="Flux-halo")
plt.plot(r, v_c,              label="Total", linewidth=2)
plt.xlabel("r [kpc]")
plt.ylabel("v_c(r) [v_unit ≈ 220 km/s]")
plt.title("Rotation curve with flux-halo")
plt.legend()
plt.grid(alpha=0.3)

plt.figure(figsize=(6,4))
plt.plot(r, SFR)
plt.xlabel("r [kpc]")
plt.ylabel("SFR surface density (arb.)")
plt.title("Flux-driven star formation profile")
plt.grid(alpha=0.3)

plt.figure(figsize=(6,5))
plt.imshow(phi_2d, origin="lower",
           extent=[x.min(), x.max(), y.min(), y.max()],
           aspect="equal")
plt.colorbar(label=r"$\phi(x,y)$")
plt.xlabel("x [kpc]")
plt.ylabel("y [kpc]")
plt.title("2D flux field in the galactic plane")

plt.tight_layout()
plt.show()
