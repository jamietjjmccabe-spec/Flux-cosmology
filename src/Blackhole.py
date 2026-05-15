import math
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Global cosmology parameters
# -----------------------------

T_MAX = 50.0       # total time (arbitrary units, say years or seconds)
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
#  FLUX COSMOLOGY STEPPER  (same logic as your lunar model)
# ======================================================================

def step_cosmology(a, S, dt, g1):
    if a <= 1e-8:
        a = 1e-8

    S_max = s0 * a * a
    Phi   = max(S_max - S, 0.0)
    flux_fraction = 0 if S_max == 0 else Phi / S_max

    # entropy evolution
    dS_dt = gamma * Phi
    S_new = S + dS_dt * dt

    # Hubble parameter
    a_safe = max(a, 1e-4)
    H_sq = Omega_m/(a_safe**3) + Omega_r/(a_safe**4) + k_flux*flux_fraction
    H_sq = max(H_sq, 0.0)
    H = math.sqrt(H_sq)

    # scale factor
    a_new = a + (H * a) * dt

    # variable G, c from flux
    G_eff = G0 * (1.0 + g1 * flux_fraction)

    if H0 != 0:
        c_eff = c0 * (1 + 0.3*(H/H0 - 1))
    else:
        c_eff = c0

    return a_new, S_new, H, flux_fraction, G_eff, c_eff


# ======================================================================
#  BLACK HOLE RADIAL INFALL MODEL
# ======================================================================

def run_blackhole_sim(g1=-1e-6,
                      M_BH=1.0,
                      r0=10.0,
                      v0=0.0):
    """
    Toy model: test particle falling radially into a black hole
    in your flux cosmology.

    - g1 : flux → G_eff coupling
    - M_BH : black hole mass (dimensionless)
    - r0 : initial radius (in units where r_s ~ O(1–10))
    - v0 : initial radial velocity (negative = inward)
    """

    n_steps = int(T_MAX / DT)

    # arrays for output
    t_vals = []
    r_vals = []
    rs_vals = []
    G_vals = []
    c_vals = []
    flux_vals = []
    dil_vals = []      # gravitational time dilation factor
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

        # Schwarzschild radius for current G,c
        rs = 2.0 * G_eff * M_BH / (c_eff**2 + 1e-12)
        rs_vals.append(rs)

        # stop if we've crossed the horizon (or extremely close)
        if r <= 1.01 * rs:
            # clamp and replicate last values
            r_vals.append(r)
            dil_vals.append(0.0)
            mip_stack.append(float('inf'))
            # fill remaining steps with flat values and break
            for j in range(i+1, n_steps):
                t_vals.append(j*DT)
                r_vals.append(r)
                rs_vals.append(rs)
                G_vals.append(G_eff)
                c_vals.append(c_eff)
                flux_vals.append(flux_frac)
                dil_vals.append(0.0)
                mip_stack.append(float('inf'))
            break

        # ---- radial infall dynamics ----
        # Newtonian-ish radial acceleration with a GR-ish correction factor
        # (1 - rs/r) damps the force as you get close to the horizon.
        grav_factor = - G_eff * M_BH / (r*r)
        correction  = max(1.0 - rs/r, 0.0)
        a_r = grav_factor * correction

        v  += a_r * DT
        r  += v   * DT

        r_vals.append(r)

        # ---- gravitational time dilation & "MIP stacking" ----
        # GR: dτ = dt * sqrt(1 - rs/r)
        # interpret: local tick = global tick * sqrt(1 - rs/r)
        # your language: more MIP events per global tick as rs/r → 1
        if r > rs:
            factor = math.sqrt(1.0 - rs/r)
        else:
            factor = 0.0

        dil_vals.append(factor)                 # local_tick / global_tick
        mip_stack.append(1.0 / max(factor, 1e-9))  # how many global ticks per local

    # trim everything to same length (if we broke early)
    n = len(r_vals)
    t_vals      = t_vals[:n]
    rs_vals     = rs_vals[:n]
    G_vals      = G_vals[:n]
    c_vals      = c_vals[:n]
    flux_vals   = flux_vals[:n]
    dil_vals    = dil_vals[:n]
    mip_stack   = mip_stack[:n]

    return {
        "t": np.array(t_vals),
        "r": np.array(r_vals),
        "rs": np.array(rs_vals),
        "G": np.array(G_vals),
        "c": np.array(c_vals),
        "flux": np.array(flux_vals),
        "dil": np.array(dil_vals),
        "mip": np.array(mip_stack),
    }


# ======================================================================
#  MAIN: RUN A BLACK HOLE FLUX SIM
# ======================================================================

if __name__ == "__main__":
    # choose parameters
    g1 = -1e-6      # same sign as in lunar model: flux hardening → stronger G
    M_BH = 1.0      # mass of the black hole
    r0   = 10.0     # start 10 "units" away from center
    v0   = 0.0      # dropped from rest

    data = run_blackhole_sim(g1=g1, M_BH=M_BH, r0=r0, v0=v0)

    t   = data["t"]
    r   = data["r"]
    rs  = data["rs"]
    dil = data["dil"]
    mip = data["mip"]

    # 1) radius vs time with horizon
    plt.figure()
    plt.plot(t, r, label="r(t)")
    plt.plot(t, rs, "--", label="r_s(t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("radius")
    plt.title("Radial infall toward black hole in flux cosmology")
    plt.legend()

    # 2) gravitational time dilation (local tick / global tick)
    plt.figure()
    plt.plot(t, dil)
    plt.xlabel("t (ticks)")
    plt.ylabel("dτ/dt")
    plt.title("Local tick vs global tick (time dilation)")

    # 3) MIP stacking = how many global ticks per local tick
    plt.figure()
    plt.plot(t, mip)
    plt.yscale("log")
    plt.xlabel("t (ticks)")
    plt.ylabel("MIP stacking factor")
    plt.title("MIP density near black hole (log scale)")

    # 4) G_eff & c_eff evolution (from flux field)
    plt.figure()
    plt.plot(t, data["G"], label="G_eff(t)")
    plt.plot(t, data["c"], label="c_eff(t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("G_eff, c_eff")
    plt.title("Effective G and c from flux")
    plt.legend()

    plt.show()
