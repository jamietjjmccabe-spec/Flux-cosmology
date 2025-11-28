import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# ============================================================
# 0. BASIC SETUP: UNITS & GRIDS
# ============================================================

# Units:
#   length unit  = 1 kpc
#   velocity     = arbitrary (we don't integrate orbits here)
#   time unit    = Gyr-like toy

# Radial grid
r_max = 25.0
n_r   = 200
r = np.linspace(0.3, r_max, n_r)  # avoid r=0

# Angular grid
n_theta = 256
theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)

# 2D polar grids
R, TH = np.meshgrid(r, theta, indexing="ij")

# ============================================================
# 1. INITIAL RADIAL FLUX & GAS
# ============================================================

# Flux profile: inner peak + outer floor
phi0    = 1.0
R_phi   = 2.5
phi_halo = 0.3
R_h      = 8.0

phi_init = phi0 * np.exp(-r / R_phi) + phi_halo * (1.0 - np.exp(-r / R_h))

# Gas disk: exponential + molecular ring
M_gas_total = 3.0
R_g         = 4.5

Sigma_g_init = (M_gas_total / (2.0 * np.pi * R_g**2)) * np.exp(-r / R_g)
ring_radius  = 8.0
ring_width   = 1.5
ring_boost   = np.exp(-0.5 * ((r - ring_radius) / ring_width)**2)
Sigma_g_init *= (1.0 + 1.0 * ring_boost)

# Optional extra central gas (CMZ–like)
Sigma_g_init += 0.004 * np.exp(-r / 1.5)

# ============================================================
# 1b. ARMS / PHYSICS TUNING
# ============================================================

# Spiral structure
m_arm      = 2          # 2-armed grand design
A_spiral   = 0.55       # strong flux spiral amplitude
B_spiral   = 0.85       # strong gas spiral response

# Open spiral pitch angle
pitch_angle_deg = 18.0
pitch_angle     = np.deg2rad(pitch_angle_deg)
k_spiral        = 1.0 / np.tan(pitch_angle)

# Boost initial gas for visible star-forming arms
Sigma_g_init *= 3.0

# Couplings (used everywhere below)
epsilon_SF  = 0.15   # efficiency: SFR ∝ epsilon * φ * gas
gamma_star  = 0.10   # SFR depletion of φ
gamma_in    = 0.03   # inflow creating φ
gamma_relax = 0.03   # relaxation to floor
phi_floor   = 0.15

# ============================================================
# 2. TIME EVOLUTION (RADIAL, 1D)
# ============================================================

t_max = 10.0
Nt    = 200
t_arr = np.linspace(0.0, t_max, Nt)
dt    = t_arr[1] - t_arr[0]

# Inflow stronger toward center (inside-out growth)
inflow_profile = np.exp(-r / 5.0)

phi_rad   = np.zeros((Nt, n_r))
Sigma_g   = np.zeros((Nt, n_r))
SFR_rad   = np.zeros((Nt, n_r))
Sigma_star = np.zeros((Nt, n_r))
Sigma_star[0, :] = 0.0

phi_rad[0, :] = phi_init
Sigma_g[0, :] = Sigma_g_init

for n in range(Nt - 1):
    phi_n = phi_rad[n, :]
    gas_n = Sigma_g[n, :]

    # Star formation rate surface density (radial)
    SFR = epsilon_SF * phi_n * gas_n
    SFR_rad[n, :] = SFR

    # Gas evolution
    dSigma_g = -SFR

    # Flux evolution
    dphi = (
        -gamma_star * SFR
        + gamma_in * inflow_profile
        - gamma_relax * (phi_n - phi_floor)
    )

    # Build up stellar mass / light from past SFR
    # (simple running integral; could add exponential fading if you like)
    Sigma_star[n+1, :] = Sigma_star[n, :] + SFR * dt


    phi_next = phi_n + dt * dphi
    gas_next = gas_n + dt * dSigma_g

    phi_next = np.clip(phi_next, 0.01, 2.0)
    gas_next = np.clip(gas_next, 0.0, None)

    phi_rad[n+1, :] = phi_next
    Sigma_g[n+1, :] = gas_next

# Last-step SFR
SFR_rad[-1, :] = epsilon_SF * phi_rad[-1, :] * Sigma_g[-1, :]

# ============================================================
# 3. BAR & SPIRAL GEOMETRY
# ============================================================

# Bar parameters
A_bar     = 0.25         # bar amplitude in φ
R_bar     = 4.0          # bar extent
Omega_bar = 0.8          # bar pattern speed (rad / time-unit)

# Spiral pattern speed (linked to bar)
Omega_p = Omega_bar * 0.75

def bar_factor(R, TH, t):
    """Rotating m=2 bar, strongest inside R_bar."""
    phase = 2.0 * (TH - Omega_bar * t)
    radial_env = np.exp(-R / R_bar)
    return 1.0 + A_bar * radial_env * np.cos(phase)

def spiral_phase(R, TH, t):
    """
    Logarithmic spiral: m * (theta - Omega_p * t - k ln(r/R_ref))
    """
    R_ref = 3.0  # reference radius where arm crosses theta=0
    return m_arm * (TH - Omega_p * t - k_spiral * np.log(R / R_ref))

def spiral_factor_phi(R, TH, t):
    """Spiral modulation for flux."""
    phase = spiral_phase(R, TH, t)
    radial_env = 1.0 / (1.0 + (R / r_max)**2)  # gently weaker at very large R
    return 1.0 + A_spiral * radial_env * np.cos(phase)

def spiral_factor_gas(R, TH, t):
    """Spiral modulation for gas."""
    phase = spiral_phase(R, TH, t)
    radial_env = np.exp(-((R - ring_radius) / (2.0 * ring_width))**2)
    return 1.0 + B_spiral * radial_env * np.cos(phase)

# ============================================================
# 4. FUNCTIONS TO BUILD FULL 2D FIELDS AT ANY TIME (POLAR)
# ============================================================

def snapshot_2d(idx):
    """
    Build 2D φ(R,θ), gas(R,θ), and SFR(R,θ) for time index idx.
    Uses radial φ(r,t) and Σ_g(r,t) plus bar and spiral geometry.
    """
    t = t_arr[idx]
    phi_base = phi_rad[idx, :]     # radial profile
    gas_base = Sigma_g[idx, :]

    # Lift radial profiles to 2D via broadcasting
    phi_base_2d = phi_base[:, None] * np.ones_like(TH)
    gas_base_2d = gas_base[:, None] * np.ones_like(TH)

    # Apply bar and spiral modulations
    bar_mod        = bar_factor(R, TH, t)
    spiral_mod_phi = spiral_factor_phi(R, TH, t)
    spiral_mod_gas = spiral_factor_gas(R, TH, t)

    phi_2d = phi_base_2d * bar_mod * spiral_mod_phi
    gas_2d = gas_base_2d * spiral_mod_gas

    # Physical bounds
    phi_2d = np.clip(phi_2d, 0.01, 2.0)
    gas_2d = np.clip(gas_2d, 0.0, None)

    # SFR in 2D
    SFR_2d = epsilon_SF * phi_2d * gas_2d

    return phi_2d, gas_2d, SFR_2d

# ============================================================
# 5. POLAR SNAPSHOT PLOTS (OPTIONAL)
# ============================================================

snap_indices = [0, Nt//3, 2*Nt//3, Nt-1]
snap_titles  = [f"t = {t_arr[i]:.1f}" for i in snap_indices]

plt.figure(figsize=(7,5))
for i, idx in enumerate(snap_indices):
    plt.plot(r, phi_rad[idx,:], label=snap_titles[i])
plt.xlabel("r [kpc]")
plt.ylabel(r"$\phi(r,t)$")
plt.title("Radial flux evolution")
plt.legend()
plt.grid(alpha=0.3)

plt.figure(figsize=(7,5))
for i, idx in enumerate(snap_indices):
    plt.plot(r, Sigma_g[idx,:], label=snap_titles[i])
plt.xlabel("r [kpc]")
plt.ylabel(r"$\Sigma_g(r,t)$")
plt.title("Radial gas evolution")
plt.legend()
plt.grid(alpha=0.3)

# ============================================================
# 6. (x,y) CARTESIAN RENDERER
# ============================================================

# Cartesian grid
n_xy = 400
x = np.linspace(-r_max, r_max, n_xy)
y = np.linspace(-r_max, r_max, n_xy)
X, Y = np.meshgrid(x, y)

R_xy  = np.sqrt(X**2 + Y**2)
TH_xy = np.mod(np.arctan2(Y, X), 2.0 * np.pi)  # 0..2π

def polar_field_to_cartesian(field_polar):
    """
    Bilinear interpolation from polar grid (r,theta) to Cartesian (x,y).
    field_polar has shape (n_r, n_theta).
    """
    # Normalised indices in polar grid
    r_min, r_max_val = r[0], r[-1]
    r_idx = (R_xy - r_min) / (r_max_val - r_min) * (n_r - 1)
    th_idx = TH_xy / (2.0 * np.pi) * n_theta

    # Clip inside bounds
    r_idx = np.clip(r_idx, 0.0, n_r - 1.001)
    th_idx = np.clip(th_idx, 0.0, n_theta - 1.001)

    i0 = np.floor(r_idx).astype(int)
    j0 = np.floor(th_idx).astype(int)
    i1 = np.clip(i0 + 1, 0, n_r - 1)
    j1 = (j0 + 1) % n_theta   # wrap in theta

    fr = r_idx - i0
    ft = th_idx - j0

    # Bilinear interpolation
    f00 = field_polar[i0, j0]
    f10 = field_polar[i1, j0]
    f01 = field_polar[i0, j1]
    f11 = field_polar[i1, j1]

    f0 = f00 * (1 - fr) + f10 * fr
    f1 = f01 * (1 - fr) + f11 * fr
    f  = f0 * (1 - ft) + f1 * ft

    # Mask outside galaxy radius
    f[R_xy > r_max_val] = np.nan
    return f

# ============================================================
# 7. (x,y) SNAPSHOTS (TOP-DOWN GALAXY IMAGES)
# ============================================================

idx_early = snap_indices[0]
idx_late  = snap_indices[-1]

phi2d_e, gas2d_e, SFR2d_e = snapshot_2d(idx_early)
phi2d_l, gas2d_l, SFR2d_l = snapshot_2d(idx_late)

phi_xy_e  = polar_field_to_cartesian(phi2d_e)
phi_xy_l  = polar_field_to_cartesian(phi2d_l)
gas_xy_e  = polar_field_to_cartesian(gas2d_e)
gas_xy_l  = polar_field_to_cartesian(gas2d_l)
SFR_xy_e  = polar_field_to_cartesian(SFR2d_e)**0.5  # gamma for visibility
SFR_xy_l  = polar_field_to_cartesian(SFR2d_l)**0.5

extent_xy = [x.min(), x.max(), y.min(), y.max()]

plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.imshow(SFR_xy_e, origin="lower", extent=extent_xy, aspect="equal")
plt.colorbar(label="SFR (early)")
plt.xlabel("x [kpc]")
plt.ylabel("y [kpc]")
plt.title("Star formation (top-down, early)")

plt.subplot(1,2,2)
plt.imshow(SFR_xy_l, origin="lower", extent=extent_xy, aspect="equal")
plt.colorbar(label="SFR (late)")
plt.xlabel("x [kpc]")
plt.ylabel("y [kpc]")
plt.title("Star formation (top-down, late)")

plt.tight_layout()

# ============================================================
# 8. CINEMATIC (x,y) ANIMATION – ROTATING GALAXY
# ============================================================

fig_xy, ax_xy = plt.subplots(figsize=(6, 6))

# Camera rotation speed (visual only, in rotations per t_max)
camera_turns = 0.5          # half a full orbit over the whole simulation
camera_omega = 2.0 * np.pi * camera_turns / t_max  # rad / time-unit

def polar_field_to_cartesian_rot(field_polar, rot_angle):
    """
    Bilinear interpolation from polar grid (r,theta) to Cartesian (x,y),
    with an extra camera rotation angle applied.
    field_polar has shape (n_r, n_theta).
    """
    r_min, r_max_val = r[0], r[-1]

    # rotate camera: equivalent to shifting theta by -rot_angle
    TH_xy_rot = np.mod(TH_xy - rot_angle, 2.0 * np.pi)

    r_idx = (R_xy - r_min) / (r_max_val - r_min) * (n_r - 1)
    th_idx = TH_xy_rot / (2.0 * np.pi) * n_theta

    r_idx = np.clip(r_idx, 0.0, n_r - 1.001)
    th_idx = np.clip(th_idx, 0.0, n_theta - 1.001)

    i0 = np.floor(r_idx).astype(int)
    j0 = np.floor(th_idx).astype(int)
    i1 = np.clip(i0 + 1, 0, n_r - 1)
    j1 = (j0 + 1) % n_theta

    fr = r_idx - i0
    ft = th_idx - j0

    f00 = field_polar[i0, j0]
    f10 = field_polar[i1, j0]
    f01 = field_polar[i0, j1]
    f11 = field_polar[i1, j1]

    f0 = f00 * (1 - fr) + f10 * fr
    f1 = f01 * (1 - fr) + f11 * fr
    f  = f0 * (1 - ft) + f1 * ft

    f[R_xy > r_max_val] = np.nan
    return f

# Initial frame
phi2d_0, gas2d_0, SFR2d_0 = snapshot_2d(0)
SFR_xy_0 = polar_field_to_cartesian_rot(SFR2d_0, rot_angle=0.0)**0.5

im_xy = ax_xy.imshow(SFR_xy_0, origin="lower", extent=extent_xy,
                     aspect="equal")
cbar_xy = fig_xy.colorbar(im_xy, ax=ax_xy)
cbar_xy.set_label("SFR (arb.)")

ax_xy.set_xlabel("x [kpc]")
ax_xy.set_ylabel("y [kpc]")
ax_xy.set_title(f"Flux-spiral galaxy SFR  (t = {t_arr[0]:.2f})")

def update_xy(frame):
    idx = frame
    t_now = t_arr[idx]

    # physical pattern rotation from Omega_p / Omega_bar is already
    # encoded inside snapshot_2d via the bar/spiral factors
    _, _, SFR2d = snapshot_2d(idx)

    # add smooth camera orbit around z-axis
    rot_angle = camera_omega * t_now
    SFR_xy = polar_field_to_cartesian_rot(SFR2d, rot_angle=rot_angle)**0.5

    im_xy.set_data(SFR_xy)
    ax_xy.set_title(f"Flux-spiral galaxy SFR  (t = {t_now:.2f})")
    return [im_xy]

anim_xy = animation.FuncAnimation(
    fig_xy, update_xy, frames=Nt, interval=40, blit=True
)

plt.show()

# To export a cinematic mp4 (uncomment this):
# anim_xy.save("flux_spiral_galaxy_cinematic.mp4", fps=25, dpi=180)

# ============================================================
# 9. FULL RGB CINEMATIC GALAXY (x,y)
# ============================================================

fig_rgb, ax_rgb = plt.subplots(figsize=(6, 6))

# Camera rotation (same idea as before)
camera_turns = 0.5                  # half an orbit over t_max
camera_omega = 2.0 * np.pi * camera_turns / t_max  # rad / time-unit

def polar_field_to_cartesian_rot(field_polar, rot_angle):
    """
    Bilinear interpolation from polar grid (r,theta) to Cartesian (x,y),
    with an extra camera rotation angle applied.
    field_polar has shape (n_r, n_theta).
    """
    r_min, r_max_val = r[0], r[-1]

    # rotate camera: equivalent to shifting theta by -rot_angle
    TH_xy_rot = np.mod(TH_xy - rot_angle, 2.0 * np.pi)

    r_idx = (R_xy - r_min) / (r_max_val - r_min) * (n_r - 1)
    th_idx = TH_xy_rot / (2.0 * np.pi) * n_theta

    r_idx = np.clip(r_idx, 0.0, n_r - 1.001)
    th_idx = np.clip(th_idx, 0.0, n_theta - 1.001)

    i0 = np.floor(r_idx).astype(int)
    j0 = np.floor(th_idx).astype(int)
    i1 = np.clip(i0 + 1, 0, n_r - 1)
    j1 = (j0 + 1) % n_theta   # wrap in theta

    fr = r_idx - i0
    ft = th_idx - j0

    f00 = field_polar[i0, j0]
    f10 = field_polar[i1, j0]
    f01 = field_polar[i0, j1]
    f11 = field_polar[i1, j1]

    f0 = f00 * (1 - fr) + f10 * fr
    f1 = f01 * (1 - fr) + f11 * fr
    f  = f0 * (1 - ft) + f1 * ft

    f[R_xy > r_max_val] = np.nan
    return f

def norm_channel(C):
    """Normalise channel to [0,1] with 99th percentile clipping."""
    C = np.nan_to_num(C, nan=0.0)
    if np.all(C <= 0):
        return C
    vmax = np.percentile(C[C > 0], 99.0)
    if vmax <= 0:
        vmax = C.max()
    if vmax <= 0:
        vmax = 1.0
    return np.clip(C / vmax, 0.0, 1.0)

# Precompute a convenience radial grid for stars -> 2D
ONES_TH = np.ones_like(TH)

# Initial frame for RGB
idx0   = 0
t0     = t_arr[idx0]
phi2d0, gas2d0, SFR2d0 = snapshot_2d(idx0)

# Stellar light: radial history at t0 projected with bar+spiral shape
star_base0 = Sigma_star[idx0, :]
star_base0_2d = star_base0[:, None] * ONES_TH
shape0 = bar_factor(R, TH, t0) * spiral_factor_phi(R, TH, t0)
star2d0 = star_base0_2d * shape0**0.5   # softer shaping than gas/SFR

rot0 = camera_omega * t0

R_xy0 = polar_field_to_cartesian_rot(star2d0, rot0)
G_xy0 = polar_field_to_cartesian_rot(gas2d0,  rot0)
B_xy0 = polar_field_to_cartesian_rot(SFR2d0,  rot0)**0.5  # gamma

R0 = norm_channel(R_xy0)
G0 = norm_channel(G_xy0)
B0 = norm_channel(B_xy0)

RGB0 = np.dstack([R0, G0, B0])

im_rgb = ax_rgb.imshow(RGB0, origin="lower", extent=extent_xy,
                       aspect="equal")

ax_rgb.set_xlabel("x [kpc]")
ax_rgb.set_ylabel("y [kpc]")
ax_rgb.set_title(f"Flux-spiral galaxy RGB  (t = {t0:.2f})")

def update_rgb(frame):
    idx = frame
    t_now = t_arr[idx]

    phi2d, gas2d, SFR2d = snapshot_2d(idx)

    # Stellar light from integrated SFR history (radial)
    star_base = Sigma_star[idx, :]
    star_base_2d = star_base[:, None] * ONES_TH
    shape = bar_factor(R, TH, t_now) * spiral_factor_phi(R, TH, t_now)
    star2d = star_base_2d * shape**0.5

    rot_angle = camera_omega * t_now

    R_xy = polar_field_to_cartesian_rot(star2d, rot_angle)
    G_xy = polar_field_to_cartesian_rot(gas2d,  rot_angle)
    B_xy = polar_field_to_cartesian_rot(SFR2d,  rot_angle)**0.5

    Rn = norm_channel(R_xy)
    Gn = norm_channel(G_xy)
    Bn = norm_channel(B_xy)

    RGB = np.dstack([Rn, Gn, Bn])
    im_rgb.set_data(RGB)

    ax_rgb.set_title(f"Flux-spiral galaxy RGB  (t = {t_now:.2f})")
    return [im_rgb]

anim_rgb = animation.FuncAnimation(
    fig_rgb, update_rgb, frames=Nt, interval=40, blit=True
)

plt.show()

# To save the cinematic RGB movie:
# anim_rgb.save("flux_spiral_galaxy_RGB.mp4", fps=25, dpi=180)

