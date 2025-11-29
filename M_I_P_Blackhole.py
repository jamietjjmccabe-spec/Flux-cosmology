import math
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Global cosmology parameters
# -----------------------------

T_MAX = 50.0       # total time (ticks / years / whatever)
DT    = 0.01       # time step

s0    = 1.0
gamma = 0.01
Omega_m = 0.3
Omega_r = 0.0
k_flux  = 0.05

G0 = 1.0
c0 = 1.0
H0 = 1.0

# initial cosmology
a0 = 0.1
S0 = 0.0


# ======================================================================
#  FLUX COSMOLOGY STEPPER  (same structure as Moon model)
# ======================================================================

def step_cosmology(a, S, dt, g1):
    if a <= 1e-8:
        a = 1e-8

    S_max = s0 * a * a
    Phi   = max(S_max - S, 0.0)
    flux_fraction = 0.0 if S_max == 0 else Phi / S_max

    # entropy evolution
    dS_dt = gamma * Phi
    S_new = S + dS_dt * dt

    # Hubble parameter
    a_safe = max(a, 1e-4)
    H_sq = Omega_m/(a_safe**3) + Omega_r/(a_safe**4) + k_flux * flux_fraction
    H_sq = max(H_sq, 0.0)
    H = math.sqrt(H_sq)

    # scale factor
    a_new = a + (H * a) * dt

    # variable G, c from flux
    G_eff = G0 * (1.0 + g1 * flux_fraction)

    if H0 != 0:
        c_eff = c0 * (1.0 + 0.3 * (H/H0 - 1.0))
    else:
        c_eff = c0

    return a_new, S_new, H, flux_fraction, G_eff, c_eff


# ======================================================================
#  BLACK HOLE RADIAL INFALL WITH DYNAMIC FLUX HORIZON
# ======================================================================

def run_blackhole_sim(
    g1=-1e-6,
    M_BH=1.0,
    r0=10.0,
    v0=0.0,
    alpha_flux=2.0,
):
    """
    Toy model: test particle falling radially into a black hole
    in your flux cosmology.

    Parameters
    ----------
    g1 : float
        Flux → G_eff coupling (same meaning as in lunar model).
    M_BH : float
        Black hole mass (dimensionless).
    r0 : float
        Initial radius.
    v0 : float
        Initial radial velocity (negative = inward).
    alpha_flux : float
        Controls how strongly the horizon responds to flux:

        r_h(t) = r_s_GR(t) * [1 + alpha_flux * flux_fraction(t)]

        flux_fraction ~ 1  (flux-rich era)  -> horizon larger than GR
        flux_fraction ~ 0  (saturated era)  -> horizon ~ GR
    """

    n_steps = int(T_MAX / DT)

    # arrays for output
    t_vals   = []
    r_vals   = []
    rs_GR    = []
    r_flux   = []
    G_vals   = []
    c_vals   = []
    flux_vals = []
    dil_vals = []      # gravitational time dilation dτ/dt
    mip_stack = []     # "MIP density" = global ticks / local ticks

    # initial cosmology state
    a = a0
    S = S0

    # initial particle state
    r = r0
    v = v0

    for i in range(n_steps):
        t = i * DT
        t_vals.append(t)

        # ---- update cosmology ----
        a, S, H, flux_frac, G_eff, c_eff = step_cosmology(a, S, DT, g1)

        flux_vals.append(flux_frac)
        G_vals.append(G_eff)
        c_vals.append(c_eff)

        # GR Schwarzschild radius with current G,c
        rs = 2.0 * G_eff * M_BH / (c_eff**2 + 1e-12)
        rs_GR.append(rs)

        # FLUX-DEFINED HORIZON
        # r_h = rs * (1 + alpha * flux_fraction)
        r_h = rs * (1.0 + alpha_flux * flux_frac)
        r_flux.append(r_h)

        # If we've effectively crossed the flux horizon, clamp and extend
        if r <= 1.01 * r_h:
            r = 1.01 * r_h
            r_vals.append(r)
            dil_vals.append(0.0)
            mip_stack.append(float("inf"))

            # freeze everything from here on
            for j in range(i+1, n_steps):
                t_vals.append(j * DT)
                r_vals.append(r)
                rs_GR.append(rs)
                r_flux.append(r_h)
                G_vals.append(G_eff)
                c_vals.append(c_eff)
                flux_vals.append(flux_frac)
                dil_vals.append(0.0)
                mip_stack.append(float("inf"))
            break

        # ---- radial infall dynamics ----
        # Newtonian-ish radial acceleration, damped near horizon
        grav_factor = - G_eff * M_BH / (r * r)
        correction  = max(1.0 - r_h / r, 0.0)
        a_r = grav_factor * correction

        v  += a_r * DT
        r  += v * DT

        r_vals.append(r)

        # ---- gravitational time dilation & "MIP stacking" ----
        # Use flux horizon:
        # dτ = dt * sqrt(1 - r_h / r)
        if r > r_h:
            factor = math.sqrt(max(1.0 - r_h/r, 0.0))
        else:
            factor = 0.0

        dil_vals.append(factor)                   # local_tick / global_tick
        mip_stack.append(1.0 / max(factor, 1e-9)) # global ticks per local

    # ensure all arrays same length
    n = len(r_vals)
    t_vals    = np.array(t_vals[:n])
    rs_GR     = np.array(rs_GR[:n])
    r_flux    = np.array(r_flux[:n])
    G_vals    = np.array(G_vals[:n])
    c_vals    = np.array(c_vals[:n])
    flux_vals = np.array(flux_vals[:n])
    dil_vals  = np.array(dil_vals[:n])
    mip_stack = np.array(mip_stack[:n])
    r_vals    = np.array(r_vals[:n])

    return {
        "t": t_vals,
        "r": r_vals,
        "rs_GR": rs_GR,
        "r_flux": r_flux,
        "G": G_vals,
        "c": c_vals,
        "flux": flux_vals,
        "dil": dil_vals,
        "mip": mip_stack,
    }


# ======================================================================
#  MAIN: RUN A DYNAMIC-HORIZON BLACK HOLE SIM
# ======================================================================

if __name__ == "__main__":
    # parameters you can play with
    g1          = -1e-6   # flux → G coupling
    M_BH        = 1.0     # BH mass
    r0          = 10.0    # starting radius
    v0          = 0.0     # dropped from rest
    alpha_flux  = 2.0     # strength of flux-horizon coupling

    data = run_blackhole_sim(
        g1=g1,
        M_BH=M_BH,
        r0=r0,
        v0=v0,
        alpha_flux=alpha_flux,
    )

    t      = data["t"]
    r      = data["r"]
    rs_GR  = data["rs_GR"]
    r_flux = data["r_flux"]

    # 1) radius vs time with BOTH horizons
    plt.figure()
    plt.plot(t, r, label="r(t)")
    plt.plot(t, rs_GR, "--", label="r_s^GR(t)")
    plt.plot(t, r_flux, ":", label="r_h^flux(t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("radius")
    plt.title("Radial infall with GR vs flux horizon")
    plt.legend()

    # 2) Local tick vs global tick
    plt.figure()
    plt.plot(t, data["dil"])
    plt.xlabel("t (ticks)")
    plt.ylabel("dτ/dt")
    plt.title("Local tick vs global tick (flux horizon)")

    # 3) MIP stacking (log scale)
    plt.figure()
    plt.plot(t, data["mip"])
    plt.yscale("log")
    plt.xlabel("t (ticks)")
    plt.ylabel("MIP stacking factor")
    plt.title("MIP density near black hole (log scale)")

    # 4) Effective G and c
    plt.figure()
    plt.plot(t, data["G"], label="G_eff(t)")
    plt.plot(t, data["c"], label="c_eff(t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("G_eff, c_eff")
    plt.title("Effective G and c from flux")
    plt.legend()

    # 5) Flux fraction vs time (for intuition)
    plt.figure()
    plt.plot(t, data["flux"])
    plt.xlabel("t (ticks)")
    plt.ylabel("flux_fraction(t)")
    plt.title("Global flux fraction vs time")

    plt.show()
