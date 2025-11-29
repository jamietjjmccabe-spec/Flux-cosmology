import math
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
#  GLOBAL / COSMOLOGY PARAMETERS
# ============================================================

T_MAX = 15.0      # total time in "ticks" (can be thought of as years)
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
M_BH      = 10.0   # larger mass than earlier toy models
r0_infall = 10.0   # characteristic radius for infall / gravity scale

# Disk geometry
R_IN  = 1.5        # inner disk radius ~ few r_s
R_OUT = 6.0        # outer disk radius
N_R   = 64         # radial grid points

# Jet geometry
Z_MAX = 20.0
N_Z   = 128        # vertical grid points

# Disk turbulence / luminosity parameters
phi_disk_bg        = 0.01   # background MIP level
disk_relax         = 0.5    # how quickly disk turbulence relaxes to background
disk_flux_coupling = 3.0    # strength of BH flux pumping the inner disk
disk_radial_index  = 2.0    # L(r) ~ r^{-q}

# Jet parameters (flux-only part)
phi_jet_bg        = 0.4
jet_relax         = 0.3
jet_flux_coupling = 1.5     # how strongly disk feeds jet base flux
jet_accel_coeff   = 0.45    # dv/dt ~ jet_accel_coeff * phi_jet (flux thrust)
jet_drag          = 0.01    # simple linear drag

# Flux effect on BH horizon
flux_horizon_boost = 3.0    # how strongly flux deficit expands horizon

# ------------------------------------------------------------
#  MAGNETIC SECTOR: semi-decohered matter still carries B
# ------------------------------------------------------------

B0      = 2.0      # base magnetic field scale (tune)
k_Bmag  = 1.0      # converts B^2 gradient into acceleration
rho0    = 1.0      # dimensionless density
nu_damp = 0.01     # extra velocity damping (MHD turbulence, reconnection)

# vertical scale over which infalling-matter B-field is important
B_infall_z_factor = 0.2     # Z_B = B_infall_z_factor * Z_MAX

# ------------------------------------------------------------
#  INFALLING MATTER RESERVOIR / ACCRETION PARAMETERS
# ------------------------------------------------------------

Mdot0        = 0.2    # base supply rate of incoherent matter
tau_decoh    = 1.5    # timescale for incoherent -> coherent / hardened matter
f_disk       = 0.7    # fraction of decohered matter feeding disk
f_jet        = 0.3    # fraction feeding jet
k_mass_disk  = 2.0    # how strongly mass feed boosts disk flux source
k_mass_jet   = 1.0    # how strongly mass feed boosts jet base flux

# ------------------------------------------------------------
#  JET COLLIMATION / EXTRA DYNAMICS
# ------------------------------------------------------------

# vertical collimation scale (in units of Z_MAX)
coll_z_factor = 0.25   # Z_coll = coll_z_factor * Z_MAX
k_coll        = 0.10  # strength of collimation acceleration

# ============================================================
#  COSMOLOGY STEPPER
# ============================================================

def step_cosmology(a, S, dt):
    # Single step of the flux cosmology.
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
    L_disk_ts     = []   # NEW: disk luminosity vs time

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
        # Schwarzschild radius (GR) r_s = 2GM/c^2 (dimensionless units)
        r_s_GR = 2.0 * G_eff * M_BH / max(c_eff**2, 1e-6)

        # Flux-defined horizon grows when flux deficit (1 - flux_frac) is large
        flux_deficit = 1.0 - flux_frac
        r_h_flux = r_s_GR * (1.0 + flux_horizon_boost * flux_deficit)

        r_s_GR_ts.append(r_s_GR)
        r_h_flux_ts.append(r_h_flux)

        # ---------------------------------------------
        # Infalling matter reservoir & decoherence
        # ---------------------------------------------
        # Base inflow rate boosted when horizon is large & flux deficit high
        Mdot_in = Mdot0 * (0.2 + flux_deficit) * (r_h_flux / max(r_s_GR, 1e-6))

        # Decoherence rate from incoherent reservoir to coherent / accreted
        Mdot_decoh = M_incoh / max(tau_decoh, 1e-6)

        # Update reservoirs
        dM_incoh = (Mdot_in - Mdot_decoh) * DT
        dM_coh   = Mdot_decoh * DT

        M_incoh += dM_incoh
        M_coh   += dM_coh

        # enforce non-negative
        M_incoh = max(M_incoh, 0.0)
        M_coh   = max(M_coh,   0.0)

        # For this toy model, take accretion rate as the decoherence rate
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
        # Source strongest at inner radii and when horizon expanded
        radial_weight = np.exp(-(r_grid - R_IN))
        radial_weight /= radial_weight.max()  # normalize

        # disk source ∝ (horizon expansion) * flux_deficit + mass feed
        disk_source_strength = (
            disk_flux_coupling * flux_deficit * (r_h_flux / max(r_s_GR, 1e-6))
            + k_mass_disk * disk_feed_from_mass
        )
        disk_source_profile = disk_source_strength * radial_weight

        # Relaxation toward background + source
        dphi_disk_dt = -disk_relax * (phi_disk - phi_disk_bg) + disk_source_profile
        phi_disk += dphi_disk_dt * DT
        phi_disk = np.clip(phi_disk, 0.0, None)

        # Disk local emissivity L(r) ~ phi_disk * r^{-q}
        emissivity = phi_disk * (r_grid ** (-disk_radial_index))
        L_disk = np.trapz(emissivity, r_grid)
        L_disk_ts.append(L_disk)

        # ---------------------------------------------
        # Jet flux: feed from disk + global flux deficit + mass feed
        # ---------------------------------------------
        base_feed = (
            jet_flux_coupling * L_disk * flux_deficit
            + k_mass_jet * jet_feed_from_mass
        )

        # flux source concentrated near z=0
        z_weight = np.exp(-z_grid / (0.3 * Z_MAX))
        z_weight /= z_weight.max()
        jet_source_profile = base_feed * z_weight

        dphi_jet_dt = -jet_relax * (phi_jet - phi_jet_bg) + jet_source_profile
        phi_jet += dphi_jet_dt * DT
        phi_jet = np.clip(phi_jet, 0.0, None)

        # -------------------------------------------------
        # MAGNETIC + FLUX ACCELERATION + GRAVITY + COLLIMATION
        # -------------------------------------------------

        # 1. Coherence fraction: semi-decohered matter still carrying B
        chi = np.clip(1.0 - phi_jet, 0.0, 1.0)  # 1 = fully coherent, 0 = fully hardened

        # 2. Magnetic field amplitude from semi-decohered matter
        #    plus additional B from infalling matter near the base
        b_infall_amp     = M_incoh / (1.0 + M_incoh)
        b_infall_profile = b_infall_amp * np.exp(-z_grid / max(Z_B, 1e-6))

        B = B0 * (chi * phi_jet + b_infall_profile)   # needs both flux + coherence

        # 3. Magnetic pressure
        P_B = k_Bmag * B**2

        # 4. Magnetic acceleration a_mag = -(1/rho) dP_B/dz
        a_mag = np.zeros_like(v_jet)
        a_mag[1:-1] = -(P_B[2:] - P_B[1:-1]) / (dz * rho0)
        a_mag[0]    = a_mag[1]     # simple reflective-ish base
        a_mag[-1]   = 0.0          # open boundary at top

        # 5. Flux-driven thrust + simple drag
        a_flux = jet_accel_coeff * phi_jet - jet_drag * v_jet

        # 6. Gravitational compression along the jet:
        #    only the coherent part (chi) feels it strongly.
        a_grav = -G_eff * M_BH / ((r0_infall + z_grid)**2)
        a_grav *= chi

        # 7. Effective collimation term:
        #    near-base booster that favours high-phi regions.
        coll_profile = np.exp(-z_grid / max(Z_coll, 1e-6))
        a_coll = k_coll * (phi_jet - phi_jet_bg) * coll_profile

        # 8. Total acceleration and velocity update
        dv_dt = a_flux + a_mag + a_coll + a_grav - nu_damp * v_jet
        v_jet += dv_dt * DT

        # Cap at sub-relativistic speeds wrt c_eff to avoid numerical blow-up
        v_max = 0.995 * c_eff
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
    L_disk_ts    = np.array(L_disk_ts)
    r_s_GR_ts    = np.array(r_s_GR_ts)
    r_h_flux_ts  = np.array(r_h_flux_ts)
    M_incoh_ts   = np.array(M_incoh_ts)
    M_coh_ts     = np.array(M_coh_ts)
    Mdot_in_ts   = np.array(Mdot_in_ts)
    Mdot_decoh_ts = np.array(Mdot_decoh_ts)
    Mdot_acc_ts  = np.array(Mdot_acc_ts)

    # Final profiles / diagnostics
    # Lorentz factor along jet axis using final c_eff (last timestep)
    c_eff_final = c_eff_ts[-1]
    beta = v_jet / max(c_eff_final, 1e-6)
    beta = np.clip(beta, 0.0, 0.999999)
    gamma = 1.0 / np.sqrt(1.0 - beta**2)

    # Build a 2D flux field phi(z,r) with radial collimation
    coll_r_profile = np.exp(-((r_grid - R_IN) / 0.8)**2)
    coll_r_profile /= coll_r_profile.max()

    phi_2d = np.outer(phi_jet, phi_disk * coll_r_profile)  # shape (N_Z, N_R)

    # Vertical flux profile averaged over radius
    phi_vertical_avg = phi_2d.mean(axis=1)

    results = {
        "t": t_arr,
        "flux_frac": flux_frac_ts,
        "G_eff": G_eff_ts,
        "c_eff": c_eff_ts,
        "L_jet": L_jet_ts,
        "L_disk": L_disk_ts,
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

    return results


# ============================================================
#  QUASAR OBSERVABLES
# ============================================================

def build_quasar_observables(
    results,
    theta_view_deg=10.0,   # viewing angle to jet axis
    D_L=1.0,               # luminosity distance in arbitrary units
    p_beam=3.0,            # beaming index (≈3 for continuous jet)
):
    """
    Turn your unified BH+disk+jet toy into a 'quasar' with observables.
    """
    t       = results["t"]
    L_disk  = results["L_disk"]
    L_jet   = results["L_jet"]
    gamma_z = results["gamma_final"]

    # 1) Choose a characteristic gamma for the part of the jet that dominates the emission.
    gamma_char = float(np.max(gamma_z))
    if gamma_char < 1.0:
        gamma_char = 1.0

    beta_char = math.sqrt(1.0 - 1.0 / gamma_char**2)

    # 2) Doppler factor δ for a given viewing angle
    theta = math.radians(theta_view_deg)
    denom = gamma_char * (1.0 - beta_char * math.cos(theta))
    denom = max(denom, 1e-6)
    delta = 1.0 / denom

    # 3) Beamed jet luminosity as seen by the observer
    L_jet_obs = (delta**p_beam) * L_jet

    # 4) Total bolometric luminosity (disk + jet)
    L_bol_int = L_disk + L_jet       # intrinsic (toy units)
    L_bol_obs = L_disk + L_jet_obs   # observed (disk ~ isotropic)

    # 5) Convert to observed fluxes (F = L / 4π D_L^2)
    four_pi_D2 = 4.0 * math.pi * (D_L**2)
    F_disk_obs = L_disk    / four_pi_D2
    F_jet_obs  = L_jet_obs / four_pi_D2
    F_tot_obs  = L_bol_obs / four_pi_D2

    return {
        "t": t,
        "L_disk_int": L_disk,
        "L_jet_int": L_jet,
        "L_bol_int": L_bol_int,
        "L_jet_obs": L_jet_obs,
        "L_bol_obs": L_bol_obs,
        "F_disk_obs": F_disk_obs,
        "F_jet_obs": F_jet_obs,
        "F_tot_obs": F_tot_obs,
        "theta_view_deg": theta_view_deg,
        "delta": delta,
        "gamma_char": gamma_char,
    }


# ============================================================
#  PLOTTING
# ============================================================

if __name__ == "__main__":
    data = run_unified()

    t = data["t"]
    z = data["z_grid"]
    r = data["r_grid"]

    # -------------------------------------------
    # QUASAR DEMO: build observed light curve
    # -------------------------------------------
    quasar = build_quasar_observables(
        data,
        theta_view_deg=5.0,   # small angle => bright blazar-like thing
        D_L=10.0,             # arbitrary distance
        p_beam=3.0,
    )

    # Quasar light curve
    plt.figure()
    plt.plot(quasar["t"], quasar["F_tot_obs"], label="Total F_obs(t)")
    plt.plot(quasar["t"], quasar["F_disk_obs"], label="Disk only")
    plt.plot(quasar["t"], quasar["F_jet_obs"],  label="Jet only")
    plt.xlabel("t (ticks)")
    plt.ylabel("Observed flux (arb. units)")
    plt.title(f"Quasar light curve (θ={quasar['theta_view_deg']}°, δ≈{quasar['delta']:.1f})")
    plt.legend()

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
    plt.imshow(
        data["phi_2d"],
        origin="lower",
        aspect="auto",
        extent=[r.min(), r.max(), z.min(), z.max()],
    )
    plt.colorbar(label="phi (flux / MIP density)")
    plt.xlabel("r")
    plt.ylabel("z")
    plt.title("Final flux field phi(z,r) with radial collimation")

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

    # 8) Infalling vs coherent matter reservoirs
    plt.figure()
    plt.plot(t, data["M_incoh"], label="M_incoh(t)  (infalling / incoherent)")
    plt.plot(t, data["M_coh"],   label="M_coh(t)   (coherent / hardened)")
    plt.xlabel("t")
    plt.ylabel("Mass (arb. units)")
    plt.title("Infalling matter & coherent reservoir")
    plt.legend()

    # 9) Accretion rate vs time
    plt.figure()
    plt.plot(t, data["Mdot_acc"], label="Mdot_acc(t)")
    plt.xlabel("t")
    plt.ylabel("Accretion rate (arb. units)")
    plt.title("Accretion rate from decoherence")
    plt.legend()

    # 10) GR vs flux-defined horizons
    plt.figure()
    plt.plot(t, data["r_s_GR"],   label="r_s^GR(t)")
    plt.plot(t, data["r_h_flux"], label="r_h^flux(t)")
    plt.xlabel("t")
    plt.ylabel("radius")
    plt.title("GR vs flux-defined horizon radii")
    plt.legend()

    plt.tight_layout()
    plt.show()
