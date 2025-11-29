import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ============================================================
# 1) FLUX COSMOLOGY MINI-MODEL (copied from your unified code)
# ============================================================

# Time in this block is in "ticks"; we will later map it to Gyr in the galaxy block.
T_MAX_COSMO = 10.0
DT_COSMO    = 0.01

# Flux–cosmology parameters
s0      = 1.0
gamma_S = 0.02
Omega_m = 0.3
Omega_r = 0.0
k_flux  = 0.05

G0 = 1.0
c0 = 1.0
H0 = 1.0

# Initial cosmology
a0 = 0.5
S0 = 0.0

def step_cosmology(a, S, dt):
    """
    Single step of the flux cosmology.
    Returns: a_new, S_new, H, flux_frac, G_eff, c_eff
    """
    if a <= 1e-8:
        a = 1e-8

    S_max = s0 * a * a
    Phi   = max(S_max - S, 0.0)
    flux_frac = 0.0 if S_max == 0 else Phi / S_max

    # entropy production
    dS_dt = gamma_S * Phi
    S_new = S + dS_dt * dt

    # Hubble parameter
    a_safe = max(a, 1e-4)
    H_sq = Omega_m / (a_safe**3) + Omega_r / (a_safe**4) + k_flux * flux_frac
    H_sq = max(H_sq, 0.0)
    H = math.sqrt(H_sq)

    # scale factor evolution
    a_new = a + (H * a) * dt

    # Effective gravity and c:
    # as flux saturates, gravity and c soften slightly
    G_eff = G0 * (1.0 - 0.1 * flux_frac)
    if H0 != 0:
        c_eff = c0 * (1.0 - 0.05 * flux_frac)
    else:
        c_eff = c0

    return a_new, S_new, H, flux_frac, G_eff, c_eff


def precompute_flux_cosmology():
    """
    Integrate the flux cosmology once and tabulate:
      t_cosmo, a(t), flux_frac(t), G_eff(t)
    """
    n_steps = int(T_MAX_COSMO / DT_COSMO)
    a = a0
    S = S0
    t = 0.0

    t_list         = []
    a_list         = []
    flux_frac_list = []
    G_eff_list     = []

    for i in range(n_steps + 1):
        t_list.append(t)
        a_list.append(a)
        # step_cosmology gives next-step quantities based on current a, S
        a_next, S_next, H, flux_frac, G_eff, c_eff = step_cosmology(a, S, DT_COSMO)
        flux_frac_list.append(flux_frac)
        G_eff_list.append(G_eff)

        a = a_next
        S = S_next
        t += DT_COSMO

    return np.array(t_list), np.array(a_list), np.array(flux_frac_list), np.array(G_eff_list)


# Precompute cosmology once at import
T_COSMO, A_COSMO, FLUX_FRAC_COSMO, G_EFF_COSMO = precompute_flux_cosmology()

# Normalise some helpful quantities
FLUX_FRAC_NORM = (FLUX_FRAC_COSMO - FLUX_FRAC_COSMO.min()) / max(
    FLUX_FRAC_COSMO.max() - FLUX_FRAC_COSMO.min(), 1e-8
)
A_NORM = (A_COSMO - A_COSMO.min()) / max(A_COSMO.max() - A_COSMO.min(), 1e-8)


def flux_quantities_at_t_cosmo(t_cosmo):
    """
    Simple linear interpolation to get a(t), flux_frac(t), G_eff(t)
    at an arbitrary cosmological time t_cosmo in [0, T_MAX_COSMO].
    """
    t_clamped = np.clip(t_cosmo, T_COSMO[0], T_COSMO[-1])
    a = np.interp(t_clamped, T_COSMO, A_COSMO)
    flux = np.interp(t_clamped, T_COSMO, FLUX_FRAC_COSMO)
    G_eff = np.interp(t_clamped, T_COSMO, G_EFF_COSMO)
    return a, flux, G_eff


# ============================================================
# 2) GALAXY MODEL PARAMETERS (now coupled to flux cosmology)
# ============================================================

# Galaxy times in Gyr; we will map 0..T_MAX_GAL_GYR -> 0..T_MAX_COSMO
T_MAX_GAL_GYR = 6.0
DT_GAL_GYR    = 0.5

N_STARS       = 25000        # total number of tracer stars
R_DISK_MAX    = 15.0         # kpc
R_SCALE       = 4.0          # kpc
Z_SCALE       = 0.3          # kpc

N_ARMS        = 2
ARM_SPREAD    = 0.35         # basic angular spread for arm stars

# Fraction of stars that start in coherent arms vs diffuse field
F_COHERENT    = 0.6

# Rotation curve parameters (base, before flux modulation)
R0_SUN        = 8.0
T_ORBIT_SUN   = 0.23         # Gyr
OMEGA0        = 2.0 * math.pi / T_ORBIT_SUN

# How strongly flux affects pattern speed and arm strength
PATTERN_SPEED_BOOST = 0.6    # multiplies around 1 ± 0.3 over the run
ARM_STRENGTH_BOOST  = 0.9    # controls relative brightness of arms vs field

POINT_SIZE     = 1.0
OUT_PREFIX     = "galaxy_flux_t_"


# ============================================================
# 3) SAMPLERS
# ============================================================

def sample_exponential_disk(n, r_scale, r_max):
    """
    Sample n radii from a truncated exponential disk:
      p(R) ∝ R * exp(-R / r_scale)  for 0 <= R <= r_max
    via rejection sampling.
    """
    radii = []
    r_peak = r_scale
    p_peak = r_peak * math.exp(-r_peak / r_scale)

    while len(radii) < n:
        R_try = np.random.uniform(0.0, r_max)
        p_try = R_try * math.exp(-R_try / r_scale)
        if np.random.uniform(0.0, p_peak) < p_try:
            radii.append(R_try)

    return np.array(radii)


# ============================================================
# 4) INITIAL GALAXY CONDITIONS (with coherent vs field stars)
# ============================================================

# Radii
R_all = sample_exponential_disk(N_STARS, R_SCALE, R_DISK_MAX)
R_all += np.random.normal(loc=0.0, scale=0.1, size=N_STARS)
R_all = np.clip(R_all, 0.0, None)

# Vertical
Z_all = np.random.normal(loc=0.0, scale=Z_SCALE, size=N_STARS)

# Coherent vs field population
is_coherent = np.random.uniform(0.0, 1.0, size=N_STARS) < F_COHERENT

# Arm assignment only matters for coherent stars
arm_index = np.random.randint(0, N_ARMS, size=N_STARS)
arm_angle_base = 2.0 * math.pi * arm_index / N_ARMS

# Initial azimuth:
# - coherent stars: clustered around arm_angle_base
# - field stars: random 0..2π
phi0 = np.empty(N_STARS)
phi0[is_coherent] = arm_angle_base[is_coherent] + np.random.normal(
    loc=0.0, scale=ARM_SPREAD, size=is_coherent.sum()
)
phi0[~is_coherent] = np.random.uniform(0.0, 2.0 * math.pi, size=(~is_coherent).sum())

# A simple proxy for "stellar mass at birth": we will normalise using a(t)
# We assign each star a formation time t_form in galaxy Gyr based on the
# cosmological a(t) – higher a(t) corresponds to later-forming stars.
# For a toy mapping, draw a uniform random x in [0,1], map to t_cosmo via
# inverse CDF of A_NORM (monotonic), then project to t_gal.
cum_A = np.cumsum(A_NORM)
cum_A /= cum_A[-1]

def invert_cdf(u):
    idx = np.searchsorted(cum_A, u, side="right")
    idx = np.clip(idx, 0, len(T_COSMO) - 1)
    return T_COSMO[idx]

u_birth = np.random.uniform(0.0, 1.0, size=N_STARS)
t_cosmo_birth = np.array([invert_cdf(u) for u in u_birth])

# --- SIMPLE MODE: all stars exist from the start ---
t_gal_birth = np.zeros(N_STARS)


# ============================================================
# 5) ROTATION LAW WITH FLUX-DEPENDENT PATTERN SPEED
# ============================================================

def omega_of_R_base(R):
    """Flat rotation curve baseline Ω(R)."""
    eps = 0.2
    return OMEGA0 * (R0_SUN / np.maximum(R, eps))


def flux_modulation_factors(t_gyr):
    """
    Given a galaxy time in Gyr, compute:
      - pattern_speed_factor(t)
      - arm_brightness_coherent(t)
      - arm_brightness_field(t)
    using the flux fraction and scale factor from the cosmology.
    """
    # Map galaxy time 0..T_MAX_GAL_GYR -> cosmology time 0..T_MAX_COSMO
    t_cosmo = (t_gyr / T_MAX_GAL_GYR) * T_MAX_COSMO
    a, flux, G_eff = flux_quantities_at_t_cosmo(t_cosmo)

    # Normalised proxies
    flux_n = (flux - FLUX_FRAC_COSMO.min()) / max(
        FLUX_FRAC_COSMO.max() - FLUX_FRAC_COSMO.min(), 1e-8
    )
    flux_n = float(np.clip(flux_n, 0.0, 1.0))

    # Pattern speed:
    # high flux (strong coherence) -> faster, more rigid spiral pattern
    pattern_speed_factor = 1.0 + PATTERN_SPEED_BOOST * (flux_n - 0.5)

    # Arm strength:
    # high flux -> bright, coherent arms; low flux -> arms fade, field dominates
    arm_brightness_coherent = 0.2 + ARM_STRENGTH_BOOST * flux_n
    arm_brightness_field    = 0.2 + ARM_STRENGTH_BOOST * (1.0 - flux_n)

    return pattern_speed_factor, arm_brightness_coherent, arm_brightness_field


# ============================================================
# 6) SNAPSHOT GENERATOR
# ============================================================

def make_snapshot_flux_coupled(t_gyr, show=False):
    """
    Snapshot of the galaxy at time t_gyr, with:
      - pattern speed scaled by flux cosmology
      - coherent vs field brightness controlled by flux_frac
      - birth times tied loosely to a(t) (later a -> later star formation)
    """
    pattern_factor, bright_coh, bright_field = flux_modulation_factors(t_gyr)

    # Time since each star's birth (if negative, star has not formed yet)
    dt_orbit = t_gyr - t_gal_birth
    alive_mask = dt_orbit > 0.0
    if not np.any(alive_mask):
        # nothing has formed yet – just skip plotting
        return

    # Only orbit stars that have formed
    R = R_all[alive_mask]
    Z = Z_all[alive_mask]
    phi_init = phi0[alive_mask]
    is_coh_alive = is_coherent[alive_mask]
    dt_alive = dt_orbit[alive_mask]

    Omega_base = omega_of_R_base(R)
    phi_t = phi_init + Omega_base * pattern_factor * dt_alive

    x = R * np.cos(phi_t)
    y = R * np.sin(phi_t)
    z = Z

    # Build per-point alpha via coherent/field brightness
    alpha_vals = np.where(is_coh_alive, bright_coh, bright_field)

    # Rescale to 0..1 for plotting
    alpha_vals = alpha_vals / alpha_vals.max()

    # 3D plot
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Use a single colormap but with alpha encoding coherence / field dominance
    colors = plt.cm.plasma(0.7 * np.ones_like(alpha_vals))

    ax.scatter(x, y, z, s=POINT_SIZE, c=colors, alpha=alpha_vals)

    R_plot = R_DISK_MAX
    ax.set_xlim(-R_plot, R_plot)
    ax.set_ylim(-R_plot, R_plot)
    ax.set_zlim(-Z_SCALE * 5, Z_SCALE * 5)

    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")
    ax.set_zlabel("z [kpc]")
    ax.set_title(f"Flux-coupled spiral disk at t = {t_gyr:.1f} Gyr")

    ax.grid(False)
    ax.view_init(elev=20, azim=135)

    fname = f"{OUT_PREFIX}{t_gyr:.1f}Gyr.png"
    plt.tight_layout()
    plt.savefig(fname, dpi=200)
    if show:
        plt.show()
    else:
        plt.close(fig)


# ============================================================
# 7) MAIN LOOP
# ============================================================

if __name__ == "__main__":
    times = np.arange(0.0, T_MAX_GAL_GYR + 1e-6, DT_GAL_GYR)
    for t in times:
        print(f"Rendering flux-coupled galaxy at t = {t:.1f} Gyr...")
        make_snapshot_flux_coupled(t_gyr=t, show=False)

    print("Done.")
