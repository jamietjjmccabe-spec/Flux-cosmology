import math
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
#  SIMPLE FLUX COSMOLOGY BACKGROUND
# ============================================================

def step_cosmology(a, S, dt):
    s0      = 1.0
    gamma_S = 0.02
    Omega_m = 0.3
    Omega_r = 0.0
    k_flux  = 0.05

    G0 = 1.0
    c0 = 1.0

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

    # mild G, c running with flux
    G_eff = G0 * (1.0 - 0.05 * flux_frac)
    c_eff = c0 * (1.0 - 0.02 * flux_frac)

    return a_new, S_new, H, flux_frac, G_eff, c_eff


# ============================================================
#  FLUX-DRIVEN CORE-COLLAPSE SUPERNOVA TOY MODEL
# ============================================================

def run_supernova():
    T_MAX = 5.0
    DT    = 0.001
    n_steps = int(T_MAX / DT)

    # time series
    t_arr        = np.linspace(0.0, T_MAX, n_steps)
    R_core_arr   = np.zeros(n_steps)
    v_core_arr   = np.zeros(n_steps)
    flux_core_arr = np.zeros(n_steps)
    fuel_arr     = np.zeros(n_steps)
    R_shock_arr  = np.zeros(n_steps)
    v_shock_arr  = np.zeros(n_steps)
    flux_env_arr = np.zeros(n_steps)
    H_arr        = np.zeros(n_steps)
    E_expl_arr   = np.zeros(n_steps)

    # background
    a    = 1.0
    Senv = 0.0

    # core parameters
    M_core        = 1.4       # core mass (Chandra-ish)
    R_core        = 1.0       # initial core radius
    v_core        = 0.0
    S_core        = 0.0
    s0_core       = 1.5       # sets S_max_core ~ s0_core * R^3
    gamma_core    = 0.25      # how fast core hardens when flux deficit is large
    k_flux_support = 0.9      # outward support from live flux

    R_ns          = 0.15      # neutron-star radius floor
    k_bounce      = 40.0      # stiffness when R < R_ns
    core_damp     = 0.5       # velocity damping

    fuel_frac     = 1.0
    fuel_burn_rate = 0.05     # secular fuel usage

    # envelope / shock
    M_env         = 8.0
    eta_expl      = 0.8       # fraction of binding energy → explosion
    k_shock_drag  = 0.08
    Hubble_drag_factor = 1.0

    collapse_happened = False
    E_explosion       = 0.0
    R_shock           = R_core
    v_shock           = 0.0

    for i in range(n_steps):
        t = t_arr[i]

        # ---- cosmology / environment ----
        a, Senv, H, flux_env, G_eff, c_eff = step_cosmology(a, Senv, DT)
        flux_env_arr[i] = flux_env
        H_arr[i]        = H

        # ---- core flux + entropy ----
        R_safe = max(R_core, 0.03)
        S_max_core = s0_core * (R_safe**3)
        Phi_core   = max(S_max_core - S_core, 0.0)
        flux_frac_core = 0.0 if S_max_core == 0 else Phi_core / S_max_core

        dS_core_dt = gamma_core * Phi_core
        S_core += dS_core_dt * DT

        # simple fuel burning that speeds up when flux support is strong
        dfuel_dt = -fuel_burn_rate * fuel_frac * (0.7 + 0.6 * flux_frac_core)
        fuel_frac = max(fuel_frac + dfuel_dt * DT, 0.0)

        # ---- core equation of motion ----
        a_grav   = -G_eff * M_core / (R_safe**2)
        a_flux   =  k_flux_support * flux_frac_core / (R_safe**2)
        a_bounce =  k_bounce * (R_ns - R_core) if R_core < R_ns else 0.0

        a_net_core = a_grav + a_flux + a_bounce - core_damp * v_core
        v_core    += a_net_core * DT
        R_core    += v_core * DT

        # enforce neutron-star floor
        if R_core < R_ns:
            R_core = R_ns
            v_core = 0.0

        # ---- detect collapse once (flux-support death) ----
        if (not collapse_happened) and (flux_frac_core < 0.15):
            collapse_happened = True

            R_collapse  = R_core
            E_bind      = G_eff * (M_core**2) / R_collapse
            E_explosion = eta_expl * E_bind

            v_shock = math.sqrt(2.0 * E_explosion / max(M_env, 1e-6))
            R_shock = R_core

        # ---- shock / ejecta evolution every step after collapse ----
        if collapse_happened:
            dv_shock_dt = -k_shock_drag * v_shock * abs(v_shock) - Hubble_drag_factor * H * v_shock
            v_shock    += dv_shock_dt * DT
            v_shock     = max(v_shock, 0.0)
            R_shock    += v_shock * DT

        # ---- store ----
        R_core_arr[i]   = R_core
        v_core_arr[i]   = v_core
        flux_core_arr[i] = flux_frac_core
        fuel_arr[i]     = fuel_frac
        R_shock_arr[i]  = R_shock
        v_shock_arr[i]  = v_shock
        E_expl_arr[i]   = E_explosion

    return {
        "t": t_arr,
        "R_core": R_core_arr,
        "v_core": v_core_arr,
        "flux_core": flux_core_arr,
        "fuel": fuel_arr,
        "R_shock": R_shock_arr,
        "v_shock": v_shock_arr,
        "flux_env": flux_env_arr,
        "H": H_arr,
        "E_explosion": E_expl_arr,
    }


if __name__ == "__main__":
    data = run_supernova()
    t = data["t"]

    # 1) Core radius vs shock radius
    plt.figure()
    plt.plot(t, data["R_core"], label="Core radius R_core(t)")
    plt.plot(t, data["R_shock"], label="Shock radius R_shock(t)")
    plt.xlabel("t (toy seconds)")
    plt.ylabel("Radius (arb. units)")
    plt.title("Flux-driven core collapse and shock propagation")
    plt.legend()

    # 2) Core flux fraction
    plt.figure()
    plt.plot(t, data["flux_core"], label="Core flux fraction")
    plt.xlabel("t")
    plt.ylabel("Flux fraction Φ_core / S_max_core")
    plt.title("Core flux support decay")
    plt.legend()

    # 3) Fuel and core radial velocity
    plt.figure()
    plt.plot(t, data["fuel"], label="Fuel fraction")
    plt.plot(t, data["v_core"], label="Core radial velocity")
    plt.xlabel("t")
    plt.ylabel("Fuel, v_core")
    plt.title("Fuel exhaustion and core infall")
    plt.legend()

    # 4) Shock velocity
    plt.figure()
    plt.plot(t, data["v_shock"])
    plt.xlabel("t")
    plt.ylabel("v_shock")
    plt.title("Shock velocity vs time")

    plt.tight_layout()
    plt.show()
