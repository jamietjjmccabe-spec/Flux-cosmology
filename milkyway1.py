
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# COMMON GRIDS (polar + Cartesian)
# ============================================================

r_max = 25.0
n_r   = 200
r = np.linspace(0.3, r_max, n_r)

n_theta = 256
theta = np.linspace(0.0, 2.0*np.pi, n_theta, endpoint=False)

R, TH = np.meshgrid(r, theta, indexing='ij')

# Cartesian grid for rendering
n_xy = 256
x = np.linspace(-r_max, r_max, n_xy)
y = np.linspace(-r_max, r_max, n_xy)
X, Y = np.meshgrid(x, y)
R_xy = np.sqrt(X**2 + Y**2)
TH_xy = np.mod(np.arctan2(Y, X), 2.0*np.pi)

extent_xy = [x.min(), x.max(), y.min(), y.max()]

ONES_TH = np.ones_like(TH)


def polar_field_to_cartesian(field_polar):
    """
    Bilinear interpolation from polar (r,theta) grid to Cartesian (x,y) grid.
    field_polar has shape (n_r, n_theta).
    """
    r_min, r_max_val = r[0], r[-1]
    nr, nth = field_polar.shape

    r_idx = (R_xy - r_min) / (r_max_val - r_min) * (nr - 1)
    th_idx = TH_xy / (2.0*np.pi) * nth

    r_idx = np.clip(r_idx, 0.0, nr - 1.001)
    th_idx = np.clip(th_idx, 0.0, nth - 1.001)

    i0 = np.floor(r_idx).astype(int)
    j0 = np.floor(th_idx).astype(int)
    i1 = np.clip(i0 + 1, 0, nr - 1)
    j1 = (j0 + 1) % nth

    fr = r_idx - i0
    ft = th_idx - j0

    f00 = field_polar[i0, j0]
    f10 = field_polar[i1, j0]
    f01 = field_polar[i0, j1]
    f11 = field_polar[i1, j1]

    f0 = f00 * (1 - fr) + f10 * fr
    f1 = f01 * (1 - fr) + f11 * fr
    f  = f0 * (1 - ft) + f1 * ft

    # mask outside disk
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


# ============================================================
# MAIN GALAXY SIMULATOR
# ============================================================

def simulate_flux_galaxy(params):
    """
    Run a 1D radial flux+gas evolution, then imprint a bar+spiral pattern
    and produce an RGB (x,y) image.
    """
    # ----- unpack params -----
    m_arm      = params.get('m_arm', 2)
    A_spiral   = params.get('A_spiral', 0.4)
    B_spiral   = params.get('B_spiral', 0.6)
    pitch_deg  = params.get('pitch_deg', 18.0)
    gas_boost  = params.get('gas_boost', 3.0)
    cmz_boost  = params.get('cmz_boost', 0.004)

    epsilon_SF = params.get('epsilon_SF', 0.15)
    gamma_star = params.get('gamma_star', 0.10)
    gamma_in   = params.get('gamma_in', 0.03)
    gamma_rel  = params.get('gamma_rel', 0.03)
    phi_floor  = params.get('phi_floor', 0.15)

    A_bar      = params.get('A_bar', 0.25)
    R_bar      = params.get('R_bar', 4.0)
    Omega_bar  = params.get('Omega_bar', 0.8)
    Omega_p    = params.get('Omega_p', 0.8 * 0.75)

    t_max      = params.get('t_max', 10.0)
    Nt         = params.get('Nt', 200)
    t_snap_frac = params.get('t_snap_frac', 0.6)

    # ----- bulge params -----##
    bulge_strength = params.get('bulge_strength', 0.0)  # relative to disk
    bulge_Re       = params.get('bulge_Re', 1.0)        # kpc, half-light radius
    bulge_n        = params.get('bulge_n', 2.0)         # Sérsic index (2–4 typical)

    # ----- halo params -----##

    halo_strength = params.get('halo_strength', 0.4)   # flux halo brightness

    # pitch
    if pitch_deg <= 0:
        pitch_deg = 5.0
    pitch = np.deg2rad(pitch_deg)
    k_spiral = 1.0 / np.tan(pitch)

    # ----- initial radial flux profile -----
    phi0     = 1.0
    R_phi    = 2.5
    phi_halo = 0.3
    R_h      = 8.0
    phi_init = phi0 * np.exp(-r / R_phi) + phi_halo * (1.0 - np.exp(-r / R_h))

    # ----- initial gas profile -----
    M_gas_total = 3.0
    R_g         = 4.5
    Sigma_g_init = (M_gas_total / (2.0 * np.pi * R_g**2)) * np.exp(-r / R_g)

    ring_radius = params.get('ring_radius', 8.0)
    ring_width  = params.get('ring_width', 1.5)
    ring_boost  = np.exp(-0.5 * ((r - ring_radius) / ring_width)**2)
    Sigma_g_init *= (1.0 + 1.0 * ring_boost)

    if cmz_boost > 0:
        Sigma_g_init += cmz_boost * np.exp(-r / 1.5)

    Sigma_g_init *= gas_boost

    # ----- time grid -----
    t_arr = np.linspace(0.0, t_max, Nt)
    dt = t_arr[1] - t_arr[0]

    inflow_profile = np.exp(-r / 5.0)

    phi_rad    = np.zeros((Nt, n_r))
    Sigma_g    = np.zeros((Nt, n_r))
    SFR_rad    = np.zeros((Nt, n_r))
    Sigma_star = np.zeros((Nt, n_r))

    phi_rad[0, :]    = phi_init
    Sigma_g[0, :]    = Sigma_g_init
    Sigma_star[0, :] = 0.0

    # ----- radial evolution -----
    for n in range(Nt - 1):
        phi_n = phi_rad[n, :]
        gas_n = Sigma_g[n, :]

        SFR = epsilon_SF * phi_n * gas_n
        SFR_rad[n, :] = SFR

        dSigma_g = -SFR
        dphi = (
            -gamma_star * SFR
            + gamma_in * inflow_profile
            - gamma_rel * (phi_n - phi_floor)
        )

        phi_next = phi_n + dt * dphi
        gas_next = gas_n + dt * dSigma_g

        phi_next = np.clip(phi_next, 0.01, 2.0)
        gas_next = np.clip(gas_next, 0.0, None)

        phi_rad[n+1, :] = phi_next
        Sigma_g[n+1, :] = gas_next
        Sigma_star[n+1, :] = Sigma_star[n, :] + SFR * dt

    SFR_rad[-1, :] = epsilon_SF * phi_rad[-1, :] * Sigma_g[-1, :]

    # ----- bar & spiral geometry -----
    def bar_factor(Rloc, THloc, t):
        phase = 2.0 * (THloc - Omega_bar * t)
        radial_env = np.exp(-Rloc / R_bar)
        return 1.0 + A_bar * radial_env * np.cos(phase)

    def spiral_phase(Rloc, THloc, t):
        if m_arm == 0 or A_spiral == 0:
            return np.zeros_like(Rloc)
        R_ref = 3.0
        return m_arm * (THloc - Omega_p * t - k_spiral * np.log(Rloc / R_ref))

    def spiral_factor_phi(Rloc, THloc, t):
        if m_arm == 0 or A_spiral == 0:
            return np.ones_like(Rloc)
        phase = spiral_phase(Rloc, THloc, t)
        radial_env = 1.0 / (1.0 + (Rloc / r_max)**2)
        return 1.0 + A_spiral * radial_env * np.cos(phase)

    def spiral_factor_gas(Rloc, THloc, t):
        if m_arm == 0 or B_spiral == 0:
            return np.ones_like(Rloc)
        phase = spiral_phase(Rloc, THloc, t)
        radial_env = np.exp(-((Rloc - ring_radius) / (2.0 * ring_width))**2)
        return 1.0 + B_spiral * radial_env * np.cos(phase)

    def snapshot_2d(idx):
        t = t_arr[idx]
        phi_base = phi_rad[idx, :]
        gas_base = Sigma_g[idx, :]

        phi_base_2d = phi_base[:, None] * ONES_TH
        gas_base_2d = gas_base[:, None] * ONES_TH

        bar_mod        = bar_factor(R, TH, t)
        spiral_mod_phi = spiral_factor_phi(R, TH, t)
        spiral_mod_gas = spiral_factor_gas(R, TH, t)

        phi_2d = phi_base_2d * bar_mod * spiral_mod_phi
        gas_2d = gas_base_2d * spiral_mod_gas

        phi_2d = np.clip(phi_2d, 0.01, 2.0)
        gas_2d = np.clip(gas_2d, 0.0, None)

        SFR_2d = epsilon_SF * phi_2d * gas_2d
        return t, phi_2d, gas_2d, SFR_2d

    # choose snapshot time
    idx_snap = int(t_snap_frac * (Nt - 1))
    t_snap, phi2d, gas2d, SFR2d = snapshot_2d(idx_snap)

    # stellar light from integrated SFR
    star_base = Sigma_star[idx_snap, :]
    star_base_2d = star_base[:, None] * ONES_TH
    shape = bar_factor(R, TH, t_snap) * spiral_factor_phi(R, TH, t_snap)
    star2d = star_base_2d * (shape**0.5)

    # track bulge separately for reddening
    bulge2d = np.zeros_like(star2d)

    # ----- add Sérsic bulge component -----
    if bulge_strength > 0.0:
        r_eff = np.clip(r / bulge_Re, 1e-3, None)
        b_n = 2.0 * bulge_n - 0.327  # standard approximation
        bulge_radial = np.exp(-b_n * (r_eff**(1.0 / bulge_n)))
        bulge_radial /= bulge_radial.max()

        bulge2d = bulge_radial[:, None] * ONES_TH

        # add bulge light to stars
        star2d += bulge_strength * bulge2d



    # ----- map to (x,y) -----
    phi_xy   = polar_field_to_cartesian(phi2d)
    gas_xy   = polar_field_to_cartesian(gas2d)
    SFR_xy   = polar_field_to_cartesian(SFR2d)**0.5
    star_xy  = polar_field_to_cartesian(star2d)
    bulge_xy = polar_field_to_cartesian(bulge2d)

    # normalised bulge map for colour weighting
    bulge_mask = norm_channel(bulge_xy)  # 0 in disk, 1 in bulge core

    # ----- flux halo from phi_xy -----
    # turn φ into a smooth outer envelope, strongest at large radius
    phi_norm = norm_channel(phi_xy)

    # use global R_xy (in kpc) to weight outer region
    R_norm = np.clip(R_xy / r_max, 0.0, 1.0)

    # halo strongest at outer radii, suppressed in the inner disk
    halo_raw = phi_norm * (R_norm**1.5)

    # optional extra falloff beyond the edge to keep it soft
    halo_raw *= np.exp(-2.0 * (R_norm**2))

    halo_mask = norm_channel(halo_raw)   # 0–1 halo brightness map

    # =======================================================
    # MULTI-WAVELENGTH INTENSITY MAPS (scalar, before colour)
    # =======================================================

    # UV: dominated by recent massive star formation
    I_UV  = norm_channel(SFR_xy)

    # Optical (roughly V-band): mixed disk+bulge starlight + some young stars
    I_OPT = norm_channel(star_xy + 0.3 * SFR_xy)

    # NIR/IR: old stellar mass + bulge + a bit of flux halo
    I_NIR = norm_channel(star_xy + 1.5 * bulge_xy + 0.5 * phi_xy)



    # how strongly to redden the bulge
    alpha_R = 1.0   # boost red in bulge
    beta_G  = 0.3   # slightly suppress green
    gamma_B = 0.5   # suppress blue more

    # base channels
    R_base = star_xy
    G_base = gas_xy
    B_base = SFR_xy

    # apply bulge reddening
    R_col = R_base * (1.0 + alpha_R * bulge_mask)
    G_col = G_base * (1.0 - beta_G  * bulge_mask)
    B_col = B_base * (1.0 - gamma_B * bulge_mask)

    # ----- add flux halo glow (teal/cyan) -----
    if halo_strength > 0.0:
        # halo roughly teal: little red, more green/blue
        R_col += halo_strength * 0.2 * halo_mask
        G_col += halo_strength * 1.0 * halo_mask
        B_col += halo_strength * 1.2 * halo_mask


    # normalise channels to [0,1]
    Rch = norm_channel(R_col)
    Gch = norm_channel(G_col)
    Bch = norm_channel(B_col)

    RGB_OPT = np.dstack([Rch, Gch, Bch])

    # ===============================
    # FALSE-COLOUR PER BAND
    # ===============================

    # UV: blue/cyan
    UV_R = 0.15 * I_UV
    UV_G = 0.55 * I_UV
    UV_B = 1.00 * I_UV
    RGB_UV = np.dstack([UV_R, UV_G, UV_B])

    # NIR: red/orange
    NIR_R = 1.00 * I_NIR
    NIR_G = 0.55 * I_NIR
    NIR_B = 0.30 * I_NIR
    RGB_NIR = np.dstack([NIR_R, NIR_G, NIR_B])

    # collect bands
    bands = {
        "RGB": RGB_OPT,  # optical composite
        "UV":  RGB_UV,
        "NIR": RGB_NIR,
    }
    
    return t_snap,  bands



# ============================================================
# 3D TILT RENDERER FOR A GALAXY RGB IMAGE
# ============================================================

def tilt_RGB_image(RGB, inclination_deg=60, thickness=0.12):
    """
    Takes an RGB (x,y) galaxy image and renders it as an inclined disk.
    - inclination_deg = 0° (face-on) to 90° (edge-on)
    - thickness = vertical thickness of the disk relative to radius
    """

    inc = np.deg2rad(inclination_deg)
    cosi = np.cos(inc)

    # Input image size
    Ny, Nx, _ = RGB.shape

    # Normalized coordinate grid
    y_idx = np.linspace(-1, 1, Ny)
    x_idx = np.linspace(-1, 1, Nx)
    Xg, Yg = np.meshgrid(x_idx, y_idx)

    # Apply tilt: compress Y axis
    Yt = Yg * cosi

    # Add vertical disk thickness (Gaussian bulge)
    Zt = (Yg * np.sin(inc)) * thickness

    # Perspective projection (optional)
    perspective_strength = 0.35
    Zscale = 1 / (1 + perspective_strength * np.abs(Zt))

    Xp = Xg * Zscale
    Yp = Yt * Zscale

    # Convert projected [-1,1] space back to pixel grid
    xi = ((Xp + 1) * 0.5 * (Nx - 1)).astype(int)
    yi = ((Yp + 1) * 0.5 * (Ny - 1)).astype(int)

    # Valid bounds
    xi = np.clip(xi, 0, Nx - 1)
    yi = np.clip(yi, 0, Ny - 1)

    # Create new tilted image
    tilted = RGB[yi, xi]

    # Darken edges slightly for realism
    vignette = np.exp(-0.8 * (Xg**2 + (Yg / cosi)**2))
    tilted = tilted * vignette[..., None]

    return tilted




# ============================================================
# DEFINE GALAXY TYPES (FLUX HUBBLE SEQUENCE)
# ============================================================

galaxies = [
    ("Flux-E", {
        "m_arm": 0, "A_spiral": 0.0, "B_spiral": 0.0,
        "A_bar": 0.05, "gas_boost": 1.0,
        "epsilon_SF": 0.05, "gamma_star": 0.02,
        "gamma_in": 0.01, "gamma_rel": 0.02,
        "pitch_deg": 5.0, "t_snap_frac": 0.9,
        # big, concentrated bulge (almost entire galaxy)
        "bulge_strength": 3.0, "bulge_Re": 1.2, "bulge_n": 4.0
    }),

    ("Flux-S0", {
        "m_arm": 2, "A_spiral": 0.15, "B_spiral": 0.15,
        "A_bar": 0.15, "gas_boost": 1.5,
        "epsilon_SF": 0.08, "t_snap_frac": 0.8,
        # strong bulge, n~3
        "bulge_strength": 2.0, "bulge_Re": 1.5, "bulge_n": 3.0
    }),

    ("Flux-Sa (barred)", {
        "m_arm": 2, "A_spiral": 0.3, "B_spiral": 0.3,
        "pitch_deg": 10.0, "A_bar": 0.35,
        "gas_boost": 2.0, "t_snap_frac": 0.6,
        # big Sa bulge
        "bulge_strength": 1.5, "bulge_Re": 1.8, "bulge_n": 2.5
    }),

    ("Flux-Sb (barred)", {
        "m_arm": 2, "A_spiral": 0.55, "B_spiral": 0.85,
        "pitch_deg": 18.0, "A_bar": 0.25,
        "gas_boost": 3.0, "t_snap_frac": 0.5,
        # Milky Way: moderate bulge
        "bulge_strength": 1.0, "bulge_Re": 2.0, "bulge_n": 2.0
    }),

    ("Flux-Sc", {
        "m_arm": 3, "A_spiral": 0.6, "B_spiral": 0.9,
        "pitch_deg": 25.0, "A_bar": 0.1,
        "gas_boost": 3.0, "t_snap_frac": 0.4,
        # small bulge
        "bulge_strength": 0.5, "bulge_Re": 2.5, "bulge_n": 1.5
    }),

    ("Flux-Sd", {
        "m_arm": 4, "A_spiral": 0.7, "B_spiral": 1.0,
        "pitch_deg": 30.0, "A_bar": 0.05,
        "gas_boost": 2.5, "t_snap_frac": 0.3,
        # almost bulgeless
        "bulge_strength": 0.2, "bulge_Re": 3.0, "bulge_n": 1.5
    }),
]


# ============================================================
# 3D TILT RENDERER FOR A GALAXY RGB IMAGE (WITH DUST LANES)
# ============================================================

def tilt_RGB_image(RGB,
                   inclination_deg=60,
                   thickness=0.12,
                   dust_strength=0.7,
                   dust_scale_rad=0.6,
                   dust_scale_vert=0.08):
    """
    Takes an RGB (x,y) galaxy image and renders it as an inclined disk.
    - inclination_deg = 0° (face-on) to 90° (edge-on)
    - thickness       = vertical thickness of the disk relative to radius
    - dust_strength   = overall optical depth scale for the dust lane
    - dust_scale_rad  = radial scale (0–1) of dust along the major axis
    - dust_scale_vert = vertical scale (0–1) of dust thickness in projection
    """

    inc = np.deg2rad(inclination_deg)
    cosi = np.cos(inc)
    sini = np.sin(inc)

    # Input image size
    Ny, Nx, _ = RGB.shape

    # Normalized coordinate grid [-1,1] in both axes
    y_idx = np.linspace(-1, 1, Ny)
    x_idx = np.linspace(-1, 1, Nx)
    Xg, Yg = np.meshgrid(x_idx, y_idx)

    # --------------------------------------------------------
    # Geometric tilt + perspective
    # --------------------------------------------------------
    # Apply tilt: compress Y axis
    Yt = Yg * cosi

    # Add vertical disk thickness (Gaussian bulge)
    Zt = (Yg * sini) * thickness

    # Perspective projection
    perspective_strength = 0.35
    Zscale = 1.0 / (1.0 + perspective_strength * np.abs(Zt))

    Xp = Xg * Zscale
    Yp = Yt * Zscale

    # Convert projected [-1,1] space back to pixel grid
    xi = ((Xp + 1) * 0.5 * (Nx - 1)).astype(int)
    yi = ((Yp + 1) * 0.5 * (Ny - 1)).astype(int)

    xi = np.clip(xi, 0, Nx - 1)
    yi = np.clip(yi, 0, Ny - 1)

    # Base tilted image
    tilted = RGB[yi, xi]

    # Vignette for soft outer fading
    vignette = np.exp(-0.8 * (Xg**2 + (Yg / max(cosi, 1e-3))**2))
    tilted = tilted * vignette[..., None]

    # --------------------------------------------------------
    # Dust lane model (projected mid-plane attenuation)
    # --------------------------------------------------------
    # Only meaningful if moderately inclined
    if dust_strength > 0 and inclination_deg > 30:
        # Midplane around Y=0, strongest near center along X
        # Use normalized coords, scaled by chosen dust widths
        Xd = Xg / dust_scale_rad
        Yd = Yg / dust_scale_vert

        # Gaussian dust optical depth concentrated along the projected plane
        midplane = np.exp(-0.5 * (Xd**2)) * np.exp(-0.5 * (Yd**2))

        # Optical depth grows with inclination (edge-on → stronger)
        tau = dust_strength * (sini**1.5) * midplane

        # Apply more extinction to blue, less to red
        atten_B = np.exp(-tau)          # strong
        atten_G = np.exp(-0.7 * tau)    # medium
        atten_R = np.exp(-0.3 * tau)    # weak

        tilted[..., 0] *= atten_R  # R
        tilted[..., 1] *= atten_G  # G
        tilted[..., 2] *= atten_B  # B

    return tilted



# ============================================================
# RENDER THE FLUX HUBBLE SEQUENCE (TILTED + DUST)
# ============================================================

fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.ravel()

thickness_param = 0.18

for ax, (name, params) in zip(axes, galaxies):
    t_snap, bands = simulate_flux_galaxy(params)

    # choose inclination per type
    if "Flux-E" in name:
        inc = 20
    elif "Flux-S0" in name:
        inc = 40
    else:
        inc = 63  # spirals more inclined

    # Choose which band to render:
    img_faceon = bands["RGB"]   # or "UV", "NIR", etc.

    # Tilt it
    img_tilted = tilt_RGB_image(
        img_faceon,
        inclination_deg=inc,
        thickness=thickness_param,
        dust_strength=0.8,
        dust_scale_rad=0.7,
        dust_scale_vert=0.06
    )

    ax.imshow(img_tilted, origin="lower", extent=extent_xy, aspect="equal")
    ax.set_title(f"{name}\n(t={t_snap:.1f}, i={inc}°)")
    ax.set_xticks([])
    ax.set_yticks([])

plt.tight_layout()
plt.show()

