import math
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
#  GLOBAL / COSMOLOGY PARAMETERS (copied-style from your flux models)
# ============================================================

T_MAX = 5.0      # total time in "Gyr units" (toy)
DT    = 0.001    # timestep

# Flux–cosmology parameters (closed quantum flux style)
s0       = 1.0
gamma_S  = 0.02
Omega_m  = 0.3
Omega_r  = 0.0
k_flux   = 0.05

G0 = 1.0
c0 = 1.0
H0 = 1.0

# Initial cosmology
a0 = 0.5       # scale factor
S0 = 0.0       # hardened / entropy variable


# ============================================================
#  COSMOLOGY STEPPER (very close to your BH + jet step_cosmology)
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
    flux_frac = 0.0 if S_max == 0 else Phi / S_max  # Phi / S_max

    # entropy / hardening evolution
    dS_dt = gamma_S * Phi
    S_new = S + dS_dt * dt

    # Hubble parameter
    a_safe = max(a, 1e-4)
    H_sq = Omega_m / (a_safe**3) + Omega_r / (a_safe**4) + k_flux * flux_frac
    H_sq = max(H_sq, 0.0)
    H = math.sqrt(H_sq)

    # scale factor evolution
    a_new = a + (H * a) * dt

    # Effective gravity and c (softening as flux saturates)
    G_eff = G0 * (1.0 - 0.1 * flux_frac)
    if H0 != 0:
        c_eff = c0 * (1.0 - 0.05 * flux_frac)
    else:
        c_eff = c0

    return a_new, S_new, H, flux_frac, G_eff, c_eff


# ============================================================
#  RED GIANT TOY MODEL (FLUX-COUPLED)
# ============================================================

# Stellar parameters (dimensionless)
M_star = 1.0       # ~ 1 M_sun
M_core_init = 0.20 # initial core mass at subgiant onset
M_env_init  = M_star - M_core_init

R_init   = 1.0     # 1 unit ~ R_sun (entering subgiant phase)
L_init   = 1.0

# Flux / MIP parameters for core & envelope
phi_core_bg = 0.1
phi_env_bg  = 0.2

core_relax = 0.5      # how fast core flux relaxes to bg
env_relax  = 0.3      # how fast envelope flux relaxes to bg

k_core_flux = 3.0     # how strongly cosmological flux_deficit pumps core
k_core_grav = 2.5     # gravitational compression pumping core flux

k_env_from_core = 1.5 # coupling: core flux -> envelope flux
k_env_flux      = 0.8 # cosmological flux_deficit pumping envelope

# Shell burning / mass transfer parameters
k_burn = 0.3           # stronger burning
min_env_mass = 0.01

# Envelope flux coupling
k_env_from_core = 1.5
k_env_flux      = 2.5

# Core flux feedback
k_core_flux = 3.0


# Luminosity scaling
L0        = 1.0
alpha_env = 0.5   # weight of envelope mass in luminosity

# Radius floor
R_min = 0.2       # minimum radius (pre-subgiant)


# Radius evolution parameters
k_R_expand = 1.5
k_R_grav   = 0.3

# Initial masses
M_core_init = 0.20

# ============================================================
#  MAIN RED GIANT SIMULATOR
# ============================================================

def run_red_giant():
    n_steps = int(T_MAX / DT)

    # Time series arrays
    t_list        = []
    a_list        = []
    flux_frac_ts  = []
    G_eff_ts      = []
    c_eff_ts      = []

    M_core_ts     = []
    M_env_ts      = []
    R_ts          = []
    L_ts          = []
    T_eff_ts      = []
    phi_core_ts   = []
    phi_env_ts    = []

    # Initial cosmology
    a = a0
    S = S0

    # Initial star
    M_core = M_core_init
    M_env  = M_env_init
    R      = R_init

    phi_core = phi_core_bg
    phi_env  = phi_env_bg

    for n in range(n_steps):
        t = n * DT
        t_list.append(t)

        # -------------------------
        # Cosmology / flux background
        # -------------------------
        a, S, H, flux_frac, G_eff, c_eff = step_cosmology(a, S, DT)
        flux_deficit = 1.0 - flux_frac

        a_list.append(a)
        flux_frac_ts.append(flux_frac)
        G_eff_ts.append(G_eff)
        c_eff_ts.append(c_eff)

        # -------------------------
        # Core flux evolution
        # -------------------------
        grav_comp_term = M_core / max(R**2, 1e-6)

        dphi_core_dt = (
            -core_relax * (phi_core - phi_core_bg)
            + k_core_flux * flux_deficit
            + k_core_grav * grav_comp_term
        )
        phi_core += dphi_core_dt * DT
        phi_core = max(phi_core, 0.0)

        # -------------------------
        # Envelope flux evolution
        # -------------------------
        # coupling from core flux, weaker for very extended envelopes
        core_coupling = k_env_from_core * phi_core / (1.0 + R)
        dphi_env_dt = (
            -env_relax * (phi_env - phi_env_bg)
            + core_coupling
            + k_env_flux * flux_deficit
        )
        phi_env += dphi_env_dt * DT
        phi_env = max(phi_env, 0.0)

        # -------------------------
        # Shell burning & mass transfer
        # -------------------------
        # Burn rate ∝ core flux * core mass; clamp by available envelope mass
        burn_rate = k_burn * phi_core * M_core
        max_burn  = max((M_env - min_env_mass) / DT, 0.0)
        burn_rate = min(burn_rate, max_burn)

        dM_core = burn_rate * DT
        dM_env  = -burn_rate * DT

        M_core += dM_core
        M_env  += dM_env

        M_core = max(M_core, 0.0)
        M_env  = max(M_env, min_env_mass)

        # -------------------------
        # Radius evolution (red giant swelling)
        # -------------------------
        grav_term = k_R_grav * M_core / max(R**2, 1e-6)
        dR_dt = k_R_expand * phi_env - grav_term
        R += dR_dt * DT
        R = max(R, R_min)

        # -------------------------
        # Luminosity & effective temperature (toy scaling)
        # -------------------------

        M_tot_effective = M_core + alpha_env * M_env
        L = L0 * phi_env**2 * M_tot_effective
        # Stefan-Boltzmann scaling: L ~ R^2 T^4 => T ~ (L/R^2)^(1/4)
        T_eff = (L / max(R**2, 1e-6))**0.25

        # Record
        M_core_ts.append(M_core)
        M_env_ts.append(M_env)
        R_ts.append(R)
        L_ts.append(L)
        T_eff_ts.append(T_eff)
        phi_core_ts.append(phi_core)
        phi_env_ts.append(phi_env)

    # Pack results
    results = {
        "t": np.array(t_list),
        "a": np.array(a_list),
        "flux_frac": np.array(flux_frac_ts),
        "G_eff": np.array(G_eff_ts),
        "c_eff": np.array(c_eff_ts),
        "M_core": np.array(M_core_ts),
        "M_env": np.array(M_env_ts),
        "R": np.array(R_ts),
        "L": np.array(L_ts),
        "T_eff": np.array(T_eff_ts),
        "phi_core": np.array(phi_core_ts),
        "phi_env": np.array(phi_env_ts),
    }
    return results


# ============================================================
#  PLOTTING
# ============================================================

if __name__ == "__main__":
    data = run_red_giant()

    t = data["t"]
    R = data["R"]
    L = data["L"]
    T_eff = data["T_eff"]
    M_core = data["M_core"]
    M_env  = data["M_env"]
    phi_core = data["phi_core"]
    phi_env  = data["phi_env"]
    flux_frac = data["flux_frac"]

    # 1) Radius vs time (red giant swelling)
    plt.figure()
    plt.plot(t, R)
    plt.xlabel("t (Gyr units, toy)")
    plt.ylabel("R / R_sun")
    plt.title("Flux-driven red giant radius evolution")

    # 2) Luminosity vs time
    plt.figure()
    plt.plot(t, L)
    plt.xlabel("t")
    plt.ylabel("L (arb. units)")
    plt.title("Red giant luminosity (flux-coupled)")

    # 3) Core / envelope masses
    plt.figure()
    plt.plot(t, M_core, label="M_core")
    plt.plot(t, M_env,  label="M_env")
    plt.xlabel("t")
    plt.ylabel("Mass (M_star units)")
    plt.title("Core growth & envelope depletion")
    plt.legend()

    # 4) Effective temperature vs time
    plt.figure()
    plt.plot(t, T_eff)
    plt.xlabel("t")
    plt.ylabel("T_eff (arb.)")
    plt.title("Effective temperature evolution")

    # 5) Core & envelope flux variables
    plt.figure()
    plt.plot(t, phi_core, label="phi_core")
    plt.plot(t, phi_env,  label="phi_env")
    plt.xlabel("t")
    plt.ylabel("Flux / MIP level")
    plt.title("Core & envelope flux evolution")
    plt.legend()

    # 6) Global flux fraction vs time
    plt.figure()
    plt.plot(t, flux_frac)
    plt.xlabel("t")
    plt.ylabel("Phi / S_max (cosmic)")
    plt.title("Global flux fraction from cosmology")

    plt.tight_layout()
    plt.show()
