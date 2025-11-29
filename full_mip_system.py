import math
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
#  GLOBAL / COSMOLOGY PARAMETERS
# ============================================================

T_MAX = 10.0      # total time in "ticks" (can be thought of as years)
DT    = 0.01      # time step

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

# ============================================================
#  BLACK HOLE + DISK + JET PARAMETERS
# ============================================================

# Black hole (dimensionless)
M_BH = 10.0           # larger mass than earlier toy models
r0_infall = 10.0      # initial test particle radius (for context, optional)

# Disk geometry
R_IN  = 1.5           # inner disk radius ~ few r_s
R_OUT = 6.0           # outer disk radius
N_R   = 64            # radial grid points

# Jet geometry
Z_MAX = 20.0
N_Z   = 128           # vertical grid points

# Disk turbulence / luminosity parameters
phi_disk_bg   = 0.01  # background MIP level
disk_relax    = 0.5   # how quickly disk turbulence relaxes to background
disk_flux_coupling = 3.0   # strength of BH flux pumping the inner disk
disk_radial_index   = 2.0  # L(r) ~ r^{-q}

# Jet parameters
phi_jet_bg      = 0.4
jet_relax       = 0.3
jet_flux_coupling = 1.5    # how strongly disk feeds jet base flux
jet_accel_coeff   = 0.12   # dv/dt ~ jet_accel_coeff * phi_jet (flux thrust)
jet_drag          = 0.05   # simple linear drag

# Flux effect on BH horizon
flux_horizon_boost = 3.0   # how strongly flux deficit expands horizon

# ------------------------------------------------------------
#  MAGNETIC SECTOR: semi-decohered matter still carries B
# ------------------------------------------------------------

B0      = 1.0      # base magnetic field scale (tune)
k_Bmag  = 0.5      # converts B^2 gradient into acceleration
rho0    = 1.0      # dimensionless density
nu_damp = 0.05     # extra velocity damping (MHD turbulence, reconnection)


# ============================================================
#  COSMOLOGY STEPPER
# ============================================================

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

    # Effective gravity and c
    # as flux saturates, gravity and c soften slightly
    G_eff = G0 * (1.0 - 0.1 * flux_frac)
    if H0 != 0:
        c_eff = c0 * (1.0 - 0.05 * flux_frac)
    else:
        c_eff = c0

    return a_new, S_new, H, flux_frac, G_eff, c_eff


# ============================================================
#  MAIN UNIFIED SIMULATOR
# ============================================================

def run_unified():
    n_steps = int(T_MAX / DT)

    # time series
    t_list        = []
    flux_frac_ts  = []
    G_eff_ts      = []
    c_eff_ts      = []
    L_jet_ts      = []

    # BH horizons (for diagnostics)
    r_s_GR_ts     = []
    r_h_flux_ts   = []

    # --- spatial grids ---
    r_grid = np.linspace(R_IN, R_OUT, N_R)
    dr = r_grid[1] - r_grid[0]

    z_grid = np.linspace(0.0, Z_MAX, N_Z)
    dz = z_grid[1] - z_grid[0]

    # Disk variables
    phi_disk = np.ones(N_R) * phi_disk_bg  # MIP turbulence level

    # Jet variables
    phi_jet = np.ones(N_Z) * phi_jet_bg
    v_jet   = np.zeros(N_Z)               # vertical velocity along axis

    # Cosmology initial
    a = a0
    S = S0

    for n in range(n_steps):
        t = n * DT
        t_list.append(t)

        # ---------------------------------------------
        # Cosmology & effective constants
        # ---------------------------------------------
        a, S, H, flux_frac, G_eff, c_eff = step_cosmology(a, S, DT)

        flux_frac_ts.append(flux_frac)
        G_eff_ts.append(G_eff)
        c_eff_ts.append(c_eff)

        # ---------------------------------------------
        # Black hole horizons
        # ---------------------------------------------
        # Schwarzschild radius (GR) r_s = 2GM/c^2 (dimensionless units)
        r_s_GR = 2.0 * G_eff * M_BH / max(c_eff**2, 1e-6)

        # Flux-defined horizon grows when flux deficit (1 - flux_frac) is large
        flux_deficit = 1.0 - flux_frac
        r_h_flux = r_s_GR * (1.0 + flux_horizon_boost * flux_deficit)

        r_s_GR_ts.append(r_s_GR)
        r_h_flux_ts.append(r_h_flux)

        # ---------------------------------------------
        # Accretion disk MIP turbulence & luminosity
        # ---------------------------------------------
        # Source strongest at inner radii and when horizon expanded
        radial_weight = np.exp(-(r_grid - R_IN))
        radial_weight /= radial_weight.max()  # normalize

        # disk source ∝ (horizon expansion) * flux_deficit
        disk_source_strength = disk_flux_coupling * flux_deficit * (r_h_flux / max(r_s_GR, 1e-6))
        disk_source_profile = disk_source_strength * radial_weight

        # Relaxation toward background + source
        dphi_disk_dt = -disk_relax * (phi_disk - phi_disk_bg) + disk_source_profile
        phi_disk += dphi_disk_dt * DT
        phi_disk = np.clip(phi_disk, 0.0, None)

        # Disk local emissivity L(r) ~ phi_disk * r^{-q}
        emissivity = phi_disk * (r_grid ** (-disk_radial_index))
        L_disk = np.trapz(emissivity, r_grid)

        # ---------------------------------------------
        # Jet flux: feed from disk + global flux deficit
        # ---------------------------------------------
        base_feed = jet_flux_coupling * L_disk * flux_deficit

        # flux source concentrated near z=0
        z_weight = np.exp(-z_grid / (0.3 * Z_MAX))
        z_weight /= z_weight.max()
        jet_source_profile = base_feed * z_weight

        dphi_jet_dt = -jet_relax * (phi_jet - phi_jet_bg) + jet_source_profile
        phi_jet += dphi_jet_dt * DT
        phi_jet = np.clip(phi_jet, 0.0, None)

        # -------------------------------------------------
        # MAGNETIC + FLUX ACCELERATION (this is the new bit)
        # -------------------------------------------------
        # 1. Coherence fraction: semi-decohered matter still carrying B
        chi = np.clip(1.0 - phi_jet, 0.0, 1.0)     # 1 = fully coherent, 0 = fully hardened

        # 2. Magnetic field amplitude from semi-decohered matter
        B = B0 * chi * phi_jet                     # needs both flux + coherence

        # 3. Magnetic pressure
        P_B = k_Bmag * B**2

        # 4. Magnetic acceleration a_mag = -(1/rho) dP_B/dz
        a_mag = np.zeros_like(v_jet)
        a_mag[1:-1] = -(P_B[2:] - P_B[1:-1]) / (dz * rho0)
        a_mag[0]    = a_mag[1]     # simple reflective-ish base
        a_mag[-1]   = 0.0          # open boundary at top

        # 5. Flux-driven thrust + simple drag
        a_flux = jet_accel_coeff * phi_jet - jet_drag * v_jet

        # 6. Total acceleration and velocity update
        dv_dt = a_flux + a_mag - nu_damp * v_jet
        v_jet += dv_dt * DT

        # Cap at sub-relativistic speeds wrt c_eff to avoid numerical blow-up
        v_max = 0.9 * c_eff
        v_jet = np.clip(v_jet, 0.0, v_max)

        # Jet "luminosity proxy" ~ integral of kinetic energy flux
        # L_jet ∝ ∫ rho v^3 dz; take rho=1 for toy
        L_jet = np.trapz(v_jet**3, z_grid)
        L_jet_ts.append(L_jet)

    # Convert to arrays
    t_arr        = np.array(t_list)
    flux_frac_ts = np.array(flux_frac_ts)
    G_eff_ts     = np.array(G_eff_ts)
    c_eff_ts     = np.array(c_eff_ts)
    L_jet_ts     = np.array(L_jet_ts)
    r_s_GR_ts    = np.array(r_s_GR_ts)
    r_h_flux_ts  = np.array(r_h_flux_ts)

    # Final profiles / diagnostics
    # Lorentz factor along jet axis using final c_eff (last timestep)
    c_eff_final = c_eff_ts[-1]
    beta = v_jet / max(c_eff_final, 1e-6)
    beta = np.clip(beta, 0.0, 0.999999)
    gamma = 1.0 / np.sqrt(1.0 - beta**2)

    # Build a 2D flux field phi(z,r) as an outer product of vertical and radial
    phi_2d = np.outer(phi_jet, phi_disk)  # shape (N_Z, N_R)

    # Vertical flux profile averaged over radius
    phi_vertical_avg = phi_2d.mean(axis=1)

    results = {
        "t": t_arr,
        "flux_frac": flux_frac_ts,
        "G_eff": G_eff_ts,
        "c_eff": c_eff_ts,
        "L_jet": L_jet_ts,
        "r_s_GR": r_s_GR_ts,
        "r_h_flux": r_h_flux_ts,
        "r_grid": r_grid,
        "z_grid": z_grid,
        "phi_disk_final": phi_disk,
        "phi_jet_final": phi_jet,
        "phi_2d": phi_2d,
        "v_jet_final": v_jet,
        "gamma_final": gamma,
        "phi_vertical_avg": phi_vertical_avg,
    }
    return results


# ============================================================
#  PLOTTING
# ============================================================

if __name__ == "__main__":
    data = run_unified()

    t    = data["t"]
    z    = data["z_grid"]
    r    = data["r_grid"]

    # 1) Jet luminosity vs time
    plt.figure()
    plt.plot(t, data["L_jet"])
    plt.xlabel("t (ticks)")
    plt.ylabel("Jet luminosity proxy L(t)")
    plt.title("MIP-driven jet luminosity vs time")

    # 2) Final vertical velocity along jet axis
    plt.figure()
    plt.plot(z, data["v_jet_final"])
    plt.xlabel("z (height)")
    plt.ylabel("v(z)")
    plt.title("Final vertical velocity along jet axis")

    # 3) Final flux field phi(z,r)
    plt.figure()
    Z, R = np.meshgrid(z, r, indexing="ij")
    plt.imshow(
        data["phi_2d"],
        origin="lower",
        aspect="auto",
        extent=[r.min(), r.max(), z.min(), z.max()],
    )
    plt.colorbar(label="phi (flux / MIP density)")
    plt.xlabel("r")
    plt.ylabel("z")
    plt.title("Final flux field phi(z,r)")

    # 4) Final Lorentz factor along jet axis
    plt.figure()
    plt.plot(z, data["gamma_final"])
    plt.xlabel("z (height)")
    plt.ylabel("gamma(z)")
    plt.title("Final Lorentz factor along jet axis")

    # 5) Global flux fraction vs time
    plt.figure()
    plt.plot(t, data["flux_frac"])
    plt.xlabel("t")
    plt.ylabel("Phi / S_max")
    plt.title("Global flux fraction vs time")

    # 6) G_eff and c_eff vs time
    plt.figure()
    plt.plot(t, data["G_eff"], label="G_eff(t)")
    plt.plot(t, data["c_eff"], label="c_eff(t)")
    plt.xlabel("t")
    plt.ylabel("Effective constants")
    plt.title("G_eff and c_eff from flux cosmology")
    plt.legend()

    # 7) Final vertical flux fraction profile (averaged over radius)
    plt.figure()
    plt.plot(z, data["phi_vertical_avg"])
    plt.xlabel("z (axis)")
    plt.ylabel("phi(z)")
    plt.title("Final vertical flux fraction profile (averaged over r)")

    plt.tight_layout()
    plt.show()
