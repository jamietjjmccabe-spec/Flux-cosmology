import math
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
#  GLOBAL / COSMOLOGY PARAMETERS (unchanged physics)
# ============================================================

T_MAX = 30.0      # total time in "ticks" – now think ~ Gyr-scale
DT    = 0.05      # bigger timestep: 0.05 tick ~ 50–100 Myr if you like

# Flux–cosmology parameters (same structure, just slower evolution)
s0      = 1.0
gamma_S = 0.01       # slightly smaller: entropy production per tick
Omega_m = 0.3
Omega_r = 0.0
k_flux  = 0.05

G0 = 1.0
c0 = 1.0
H0 = 1.0

# Initial cosmology (embed the supercluster somewhere around a ~ 0.5–0.7)
a0 = 0.6
S0 = 0.0

# ============================================================
#  SUPERCLUSTER CORE + FILAMENT PARAMETERS
# ============================================================

# Effective "BH" = supercluster core mass (dimensionless).
M_BH      = 200.0      # was 10.0 – now a massive supercluster core
r0_infall = 5.0        # characteristic radius where infall gets serious

# "Disk" geometry now = supercluster core + virial radius
R_IN  = 0.3            # inner ~ 0.3 Mpc core
R_OUT = 10.0           # outer ~ 10 Mpc supercluster / filament radius
N_R   = 128            # more radial resolution

# "Jet" geometry now = large-scale filament spine
Z_MAX = 40.0           # ~ 40 Mpc filament height
N_Z   = 256            # more vertical resolution

# Disk turbulence / luminosity parameters
phi_disk_bg        = 0.02   # background MIP level in the core
disk_relax         = 0.2    # slower relaxation (cluster memory over Gyr)
disk_flux_coupling = 5.0    # supercluster strongly pumped by flux deficit
disk_radial_index  = 1.5    # L(r) ~ r^{-q}, flatter profile for cluster

# Jet/filament parameters (flux-only part)
phi_jet_bg        = 0.2
jet_relax         = 0.15
jet_flux_coupling = 2.0
jet_accel_coeff   = 0.005   # 10x weaker flux thrust
jet_drag          = 0.05    # stronger drag against IGM

# Flux effect on "horizon" (cluster core boundary)
flux_horizon_boost = 2.0    # how strongly flux deficit expands core radius

# ------------------------------------------------------------
#  MAGNETIC SECTOR: semi-decohered matter still carries B
# ------------------------------------------------------------

B0      = 1.0
k_Bmag  = 0.1       # softer magnetic pressure gradient
rho0    = 1.0
nu_damp = 0.05      # stronger generic damping

# vertical scale where infalling-matter B is important (fraction of Z_MAX)
B_infall_z_factor = 0.3

# ------------------------------------------------------------
#  INFALLING MATTER RESERVOIR / ACCRETION PARAMETERS
# ------------------------------------------------------------

Mdot0        = 5.0      # MUCH larger supply than single AGN
tau_decoh    = 5.0      # Gyr-scale "processing" time of infalling matter
f_disk       = 0.6      # most processed matter ends up in core region
f_jet        = 0.4      # some is ejected along filaments
k_mass_disk  = 3.0      # mass feed strongly boosts core activity
k_mass_jet   = 1.5      # mass feed also boosts filament flux

# ------------------------------------------------------------
#  JET COLLIMATION / EXTRA DYNAMICS
# ------------------------------------------------------------

coll_z_factor = 0.6   # collimation over most of the filament
k_coll        = 0.02  # weak "pinching" – just shapes the filament profile


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

    # Effective gravity and c:
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

    # Infalling / coherent matter reservoirs
    M_incoh_ts    = []
    M_coh_ts      = []
    Mdot_in_ts    = []
    Mdot_decoh_ts = []
    Mdot_acc_ts   = []

    # --- spatial grids ---
    r_grid = np.linspace(R_IN, R_OUT, N_R)
    z_grid = np.linspace(0.0, Z_MAX, N_Z)
    dz = z_grid[1] - z_grid[0]

    Z_B    = B_infall_z_factor * Z_MAX
    Z_coll = coll_z_factor * Z_MAX

    # Disk variables
    phi_disk = np.ones(N_R) * phi_disk_bg  # MIP turbulence level

    # Jet variables
    phi_jet = np.ones(N_Z) * phi_jet_bg
    v_jet   = np.zeros(N_Z)               # vertical velocity along axis

    # Cosmology initial
    a = a0
    S = S0

    # Mass reservoirs
    M_incoh = 0.0   # matter still incoherent / not yet decohered
    M_coh   = 0.0   # coherent/hardened matter that has reached the inner region

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
        r_s_GR = 2.0 * G_eff * M_BH / max(c_eff**2, 1e-6)

        flux_deficit = 1.0 - flux_frac
        r_h_flux = r_s_GR * (1.0 + flux_horizon_boost * flux_deficit)

        r_s_GR_ts.append(r_s_GR)
        r_h_flux_ts.append(r_h_flux)

        # ---------------------------------------------
        # Infalling matter reservoir & decoherence
        # ---------------------------------------------
        Mdot_in = Mdot0 * (0.2 + flux_deficit) * (r_h_flux / max(r_s_GR, 1e-6))
        Mdot_decoh = M_incoh / max(tau_decoh, 1e-6)

        dM_incoh = (Mdot_in - Mdot_decoh) * DT
        dM_coh   = Mdot_decoh * DT

        M_incoh += dM_incoh
        M_coh   += dM_coh

        M_incoh = max(M_incoh, 0.0)
        M_coh   = max(M_coh,   0.0)

        Mdot_acc = Mdot_decoh

        M_incoh_ts.append(M_incoh)
        M_coh_ts.append(M_coh)
        Mdot_in_ts.append(Mdot_in)
        Mdot_decoh_ts.append(Mdot_decoh)
        Mdot_acc_ts.append(Mdot_acc)

        # Mass feed split between disk and jet
        disk_feed_from_mass = f_disk * Mdot_decoh
        jet_feed_from_mass  = f_jet  * Mdot_decoh

        # ---------------------------------------------
        # Accretion disk MIP turbulence & luminosity
        # ---------------------------------------------
        radial_weight = np.exp(-(r_grid - R_IN))
        radial_weight /= radial_weight.max()

        disk_source_strength = (
            disk_flux_coupling * flux_deficit * (r_h_flux / max(r_s_GR, 1e-6))
            + k_mass_disk * disk_feed_from_mass
        )
        disk_source_profile = disk_source_strength * radial_weight

        dphi_disk_dt = -disk_relax * (phi_disk - phi_disk_bg) + disk_source_profile
        phi_disk += dphi_disk_dt * DT
        phi_disk = np.clip(phi_disk, 0.0, None)

        emissivity = phi_disk * (r_grid ** (-disk_radial_index))
        L_disk = np.trapz(emissivity, r_grid)

        # ---------------------------------------------
        # Jet flux: feed from disk + global flux deficit + mass feed
        # ---------------------------------------------
        base_feed = (
            jet_flux_coupling * L_disk * flux_deficit
            + k_mass_jet * jet_feed_from_mass
        )

        z_weight = np.exp(-z_grid / (0.3 * Z_MAX))
        z_weight /= z_weight.max()
        jet_source_profile = base_feed * z_weight

        dphi_jet_dt = -jet_relax * (phi_jet - phi_jet_bg) + jet_source_profile
        phi_jet += dphi_jet_dt * DT
        phi_jet = np.clip(phi_jet, 0.0, None)

        # -------------------------------------------------
        # MAGNETIC + FLUX ACCELERATION + GRAVITY + COLLIMATION
        # -------------------------------------------------
        chi = np.clip(1.0 - phi_jet, 0.0, 1.0)

        b_infall_amp     = M_incoh / (1.0 + M_incoh)
        b_infall_profile = b_infall_amp * np.exp(-z_grid / max(Z_B, 1e-6))

        B = B0 * (chi * phi_jet + b_infall_profile)

        P_B = k_Bmag * B**2

        a_mag = np.zeros_like(v_jet)
        a_mag[1:-1] = -(P_B[2:] - P_B[1:-1]) / (dz * rho0)
        a_mag[0]    = a_mag[1]
        a_mag[-1]   = 0.0

        a_flux = jet_accel_coeff * phi_jet - jet_drag * v_jet

        a_grav = -G_eff * M_BH / ((r0_infall + z_grid)**2)
        a_grav *= chi

        coll_profile = np.exp(-z_grid / max(Z_coll, 1e-6))
        a_coll = k_coll * (phi_jet - phi_jet_bg) * coll_profile

        dv_dt = a_flux + a_mag + a_coll + a_grav - nu_damp * v_jet
        v_jet += dv_dt * DT

        # cluster-like cap
        v_max = 0.2 * c_eff
        v_jet = np.clip(v_jet, 0.0, v_max)

        L_jet = np.trapz(v_jet**3, z_grid)
        L_jet_ts.append(L_jet)

    # ---- pack results ----
    t_arr        = np.array(t_list)
    flux_frac_ts = np.array(flux_frac_ts)
    G_eff_ts     = np.array(G_eff_ts)
    c_eff_ts     = np.array(c_eff_ts)
    L_jet_ts     = np.array(L_jet_ts)
    r_s_GR_ts    = np.array(r_s_GR_ts)
    r_h_flux_ts  = np.array(r_h_flux_ts)
    M_incoh_ts   = np.array(M_incoh_ts)
    M_coh_ts     = np.array(M_coh_ts)
    Mdot_in_ts   = np.array(Mdot_in_ts)
    Mdot_decoh_ts = np.array(Mdot_decoh_ts)
    Mdot_acc_ts  = np.array(Mdot_acc_ts)

    c_eff_final = c_eff_ts[-1]
    beta = v_jet / max(c_eff_final, 1e-6)
    beta = np.clip(beta, 0.0, 0.999999)
    gamma = 1.0 / np.sqrt(1.0 - beta**2)

    coll_r_profile = np.exp(-((r_grid - R_IN) / 0.8)**2)
    coll_r_profile /= coll_r_profile.max()

    phi_2d = np.outer(phi_jet, phi_disk * coll_r_profile)
    phi_vertical_avg = phi_2d.mean(axis=1)

    return {
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
        "M_incoh": M_incoh_ts,
        "M_coh": M_coh_ts,
        "Mdot_in": Mdot_in_ts,
        "Mdot_decoh": Mdot_decoh_ts,
        "Mdot_acc": Mdot_acc_ts,
    }


# ============================================================
#  PLOTTING
# ============================================================

if __name__ == "__main__":
    data = run_unified()

    t = data["t"]
    z = data["z_grid"]
    r = data["r_grid"]

    plt.figure()
    plt.plot(t, data["L_jet"])
    plt.xlabel("t (ticks ~ Gyr)")
    plt.ylabel("Filament power proxy L(t)")
    plt.title("Supercluster filament power vs time")

    plt.figure()
    plt.imshow(
        data["phi_2d"],
        origin="lower",
        aspect="auto",
        extent=[r.min(), r.max(), z.min(), z.max()],
    )
    plt.colorbar(label="phi (flux / MIP density)")
    plt.xlabel("r (Mpc)")
    plt.ylabel("z (Mpc)")
    plt.title("Supercluster node + filament flux field phi(z,r)")

    plt.figure()
    plt.plot(z, data["gamma_final"])
    plt.xlabel("z (height)")
    plt.ylabel("gamma(z)")
    plt.title("Final Lorentz factor along jet axis")

    plt.figure()
    plt.plot(t, data["flux_frac"])
    plt.xlabel("t")
    plt.ylabel("Phi / S_max")
    plt.title("Global flux fraction vs time")

    plt.figure()
    plt.plot(t, data["G_eff"], label="G_eff(t)")
    plt.plot(t, data["c_eff"], label="c_eff(t)")
    plt.xlabel("t")
    plt.ylabel("Effective constants")
    plt.title("G_eff and c_eff from flux cosmology")
    plt.legend()

    plt.figure()
    plt.plot(z, data["phi_vertical_avg"])
    plt.xlabel("z (axis)")
    plt.ylabel("phi(z)")
    plt.title("Final vertical flux fraction profile (averaged over r)")

    plt.figure()
    plt.plot(t, data["M_incoh"], label="M_incoh(t)  (infalling / incoherent)")
    plt.plot(t, data["M_coh"],   label="M_coh(t)   (coherent / hardened)")
    plt.xlabel("t")
    plt.ylabel("Mass (arb. units)")
    plt.title("Infalling matter & coherent reservoir")
    plt.legend()

    plt.figure()
    plt.plot(t, data["Mdot_acc"], label="Mdot_acc(t)")
    plt.xlabel("t")
    plt.ylabel("Accretion rate (arb. units)")
    plt.title("Accretion rate from decoherence")
    plt.legend()

    plt.figure()
    plt.plot(t, data["r_s_GR"],   label="r_s^GR(t)")
    plt.plot(t, data["r_h_flux"], label="r_h^flux(t)")
    plt.xlabel("t")
    plt.ylabel("radius")
    plt.title("GR vs flux-defined horizon radii")
    plt.legend()

    plt.tight_layout()
    plt.show()
