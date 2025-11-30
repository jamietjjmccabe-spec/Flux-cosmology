import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# GLOBAL GRID
# ==========================================

NR = 1000
NZ = 1000
RMAX = 2.0
ZMAX = 4.0

r = np.linspace(-RMAX, RMAX, NR)
z = np.linspace(-ZMAX, ZMAX, NZ)
R, Z = np.meshgrid(r, z)

# ==========================================
# BASE BIPOLAR GEOMETRY
# ==========================================

v_lobe = 1.6
t = 2.4
z_extent = v_lobe * t

waist_radius = 0.35
curvature    = 1.7
shell_thick  = 0.035

# Smooth bipolar shell
r_base = waist_radius + (np.abs(Z) / z_extent) ** curvature

# ==========================================
# PRECESSION (S-SHAPE)
# ==========================================

A_p     = 0.35   # amplitude of precession
omega_p = 0.50   # frequency

# sinusoidal S-bend + small linear bias
x_shift = A_p * np.sin(omega_p * Z) + 0.12 * (Z / z_extent)

# Use precessed radius
R_pre = R - x_shift

# ==========================================
# FRACTAL TURBULENCE (CHAOTIC WALLS)
# ==========================================

def fractal_noise(Z, octaves=6):
    noise = np.zeros_like(Z)
    amp   = 0.18
    freq  = 1.3
    for _ in range(octaves):
        phase1 = 2*np.pi*np.random.rand()
        phase2 = 2*np.pi*np.random.rand()
        noise += amp * np.sin(freq * np.abs(Z)       + phase1)
        noise += amp * 0.7 * np.cos(freq * np.abs(Z)**1.15 + phase2)
        amp  *= 0.55
        freq *= 1.75
    return noise

fractal   = fractal_noise(Z, octaves=6)
wrinkles  = 0.04 * np.sin(18 * Z + 2*np.pi*np.random.rand())
turb      = fractal + wrinkles

# strongest turbulence near waist, fades outward
turb *= (1.3 * np.exp(-np.abs(Z)/3.0) * (1.0 + np.exp(-(Z**2)/0.7)))

# final shell radius
r_shell = r_base + turb

# shell mask
shell = np.abs(np.abs(R_pre) - r_shell) < shell_thick
shell &= (np.abs(Z) < z_extent)

# ==========================================
# INNER HOT CAVITY
# ==========================================

inner = (R_pre**2 + (0.55*Z)**2) < (0.55**2)

# ==========================================
# DUST LANE (EQUATORIAL BELT)
# ==========================================

dust_depth = 0.8   # how dark the lane is
dust_sigma = 0.25  # vertical thickness
dust_rad   = 0.9   # radial extent

dust_mask = np.exp(-(Z**2) / (2*dust_sigma**2)) * (np.abs(R_pre) < dust_rad)
dust_attenuation = 1.0 - dust_depth * dust_mask   # multiply emissivity by this

# ==========================================
# JET KNOTS / BULLETS
# ==========================================

knots = np.array([0.8, 1.6, 2.4, 3.2,
                 -0.8, -1.6, -2.4, -3.2])

knot_rad_sigma = 0.09
knot_z_sigma   = 0.13

emiss_knots = np.zeros_like(R)
for z0 in knots:
    gauss = (
        np.exp(-(R_pre**2) / (2*knot_rad_sigma**2)) *
        np.exp(-((Z - z0)**2) / (2*knot_z_sigma**2))
    )
    emiss_knots += gauss

# normalize knot brightness
emiss_knots /= emiss_knots.max() + 1e-12
emiss_knots *= 0.7   # relative strength vs shell

# ==========================================
# EMISSION MAP (COMBINED)
# ==========================================

emiss = np.zeros_like(R)
emiss[shell] = 1.0
emiss[inner] = 0.4
emiss += emiss_knots

# apply dust attenuation
emiss *= dust_attenuation

# rescale
emiss /= emiss.max() + 1e-12

# ==========================================
# VELOCITY FIELD (DOPPLER MAP)
# ==========================================

# simple bipolar flow: away from center along Z
v0 = 1.0
v_z = v0 * np.tanh(Z / 1.0)    # saturates outward
# small precession contribution (tilt)
v_z += 0.2 * v0 * np.sin(omega_p * Z)

# normalize to [-1,1] for plotting
v_norm = v_z / (np.max(np.abs(v_z)) + 1e-12)

# mask velocities where emission is weak
vel_mask = emiss > 0.05
v_plot   = np.zeros_like(v_norm)
v_plot[vel_mask] = v_norm[vel_mask]

# ==========================================
# PLOTS
# ==========================================

# 1) Emission (for textures / brightness)
plt.figure(figsize=(6,10))
plt.imshow(
    emiss,
    origin="lower",
    extent=[-RMAX, RMAX, -ZMAX, ZMAX],
    cmap="inferno",
    aspect="auto"
)
plt.colorbar(label="Emission")
plt.title("Red Spider Nebula – Full Toy Model\n(fractal walls, precession, dust lane, jet knots)")
plt.xlabel("r (symmetric)")
plt.ylabel("z")
plt.tight_layout()

# 2) Doppler map (velocity field masked by emission)
plt.figure(figsize=(6,10))
plt.imshow(
    v_plot,
    origin="lower",
    extent=[-RMAX, RMAX, -ZMAX, ZMAX],
    cmap="seismic",        # blue = approaching, red = receding
    aspect="auto",
    vmin=-1, vmax=1
)
plt.colorbar(label="Normalized line-of-sight velocity")
plt.contour(
    emiss,
    levels=[0.2, 0.5, 0.8],
    colors="k",
    linewidths=0.4,
    origin="lower",
    extent=[-RMAX, RMAX, -ZMAX, ZMAX],
)
plt.title("Red Spider Nebula – Doppler Velocity Map (toy)")
plt.xlabel("r (symmetric)")
plt.ylabel("z")
plt.tight_layout()

plt.show()
