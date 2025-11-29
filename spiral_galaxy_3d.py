import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed for 3D)

# ============================================================
# Parameters
# ============================================================

N_STARS      = 20000       # number of stars in the disk
R_DISK_MAX   = 15.0        # kpc, outer radius of the visible disk
R_SCALE      = 4.0         # kpc, exponential scale length
Z_SCALE      = 0.3         # kpc, vertical thickness (Gaussian sigma)

N_ARMS       = 2           # number of major arms
ARM_SPREAD   = 0.35        # rad, how wide the initial arm sectors are

T_MAX_GYR    = 6.0         # total evolution time in Gyr
DT_GYR       = 0.5         # time step in Gyr between snapshots

# Rotation curve (flat) parameters (toy but roughly Milky Way-ish)
R0_SUN       = 8.0         # kpc, reference radius
T_ORBIT_SUN  = 0.23        # Gyr, orbital period at R0 ~ 220 km/s
OMEGA0       = 2.0 * np.pi / T_ORBIT_SUN  # rad / Gyr at R = R0

# Plot / export options
POINT_SIZE   = 1.0         # scatter point size
OUT_PREFIX   = "galaxy_t_" # filename prefix


# ============================================================
# Utility: sample radii from an exponential disk
# p(R) ~ R * exp(-R / R_SCALE) on [0, R_DISK_MAX]
# ============================================================

def sample_exponential_disk(n, r_scale, r_max):
    """
    Sample n radii from a truncated exponential disk:
      p(R) ∝ R * exp(-R / r_scale)  for 0 <= R <= r_max
    Using inverse transform on the cumulative with truncation.
    """
    # Use rejection sampling for simplicity (fine at these Ns)
    radii = []
    r_peak = r_scale  # where R * exp(-R/Rs) peaks
    p_peak = r_peak * np.exp(-r_peak / r_scale)

    while len(radii) < n:
        R_try = np.random.uniform(0.0, r_max)
        p_try = R_try * np.exp(-R_try / r_scale)
        if np.random.uniform(0.0, p_peak) < p_try:
            radii.append(R_try)

    return np.array(radii)


# ============================================================
# Initial conditions
# ============================================================

# Radii: exponential disk
R = sample_exponential_disk(N_STARS, R_SCALE, R_DISK_MAX)

# Vertical positions: Gaussian disk thickness
Z = np.random.normal(loc=0.0, scale=Z_SCALE, size=N_STARS)

# Arms: assign each star to one of N_ARMS arms
arm_index = np.random.randint(0, N_ARMS, size=N_STARS)
arm_angle_base = 2.0 * np.pi * arm_index / N_ARMS

# Initial azimuth: concentrated around arm_angle_base
phi0 = arm_angle_base + np.random.normal(loc=0.0, scale=ARM_SPREAD, size=N_STARS)

# Small random noise in R to avoid perfect rings
R += np.random.normal(loc=0.0, scale=0.1, size=N_STARS)
R = np.clip(R, 0.0, None)

# Pre-generate colours by arm, just for visual distinction
arm_colors = plt.cm.plasma(np.linspace(0.1, 0.9, N_ARMS))
colors = arm_colors[arm_index]


# ============================================================
# Rotation law
# ============================================================

def omega_of_R(R):
    """
    Angular speed Ω(R) for a flat rotation curve.
    We keep it simple: Ω(R) = Ω0 * (R0 / max(R, eps)).
    Units: rad / Gyr if R in kpc.
    """
    eps = 0.2
    return OMEGA0 * (R0_SUN / np.maximum(R, eps))


# ============================================================
# Snapshot generator
# ============================================================

def make_snapshot(t_gyr, show=False):
    """
    Generate and save a 3D snapshot of the galaxy at time t_gyr (in Gyr).
    The arms are formed by differential rotation shearing the initial sectors.
    """
    # Evolve angles
    Omega = omega_of_R(R)           # rad / Gyr
    phi_t = phi0 + Omega * t_gyr    # current angle

    # Cartesian coordinates
    x = R * np.cos(phi_t)
    y = R * np.sin(phi_t)
    z = Z

    # 3D plot
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z, s=POINT_SIZE, c=colors, alpha=0.8)

    # Aesthetic tweaks
    R_plot = R_DISK_MAX
    ax.set_xlim(-R_plot, R_plot)
    ax.set_ylim(-R_plot, R_plot)
    ax.set_zlim(-Z_SCALE * 5, Z_SCALE * 5)

    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")
    ax.set_zlabel("z [kpc]")
    ax.set_title(f"Spiral disk at t = {t_gyr:.1f} Gyr")

    # Turn off grid for cleaner look
    ax.grid(False)

    # Slight viewing angle tweak so you see the spiral + thickness
    ax.view_init(elev=20, azim=135)

    # Save
    fname = f"{OUT_PREFIX}{t_gyr:.1f}Gyr.png"
    plt.tight_layout()
    plt.savefig(fname, dpi=200)
    if show:
        plt.show()
    else:
        plt.close(fig)


# ============================================================
# Main loop: generate snapshots every DT_GYR
# ============================================================

if __name__ == "__main__":
    times = np.arange(0.0, T_MAX_GYR + 1e-6, DT_GYR)
    for t in times:
        print(f"Rendering t = {t:.1f} Gyr...")
        make_snapshot(t_gyr=t, show=False)

    print("Done. PNG frames written to disk.")
