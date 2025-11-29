import math
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
#  GLOBAL COSMOLOGY / FLUX PARAMETERS (same family as before)
# =========================================================

T_MAX = 50.0       # total time (arbitrary units)
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

# black hole properties
M_BH = 1.0


# =========================================================
#  FLUX COSMOLOGY STEPPER (same logic as lunar + BH models)
# =========================================================

def step_cosmology(a, S, dt, g1):
    """
    Single global tick of your flux cosmology:
    - updates scale factor a(t)
    - updates entropy S(t)
    - computes flux fraction Phi/S_max
    - returns effective G_eff(t), c_eff(t)
    """
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


# =========================================================
#  DYNAMIC FLUX HORIZON (no particles yet)
# =========================================================

def compute_horizons(flux_frac, G_eff, c_eff, alpha_flux=2.0):
    """
    Given current flux fraction + G_eff + c_eff, return:
      rs_GR   : GR Schwarzschild radius
      r_flux  : flux-defined "MIP horizon"
    """
    rs_GR = 2.0 * G_eff * M_BH / (c_eff**2 + 1e-12)
    r_flux = rs_GR * (1.0 + alpha_flux * flux_frac)
    return rs_GR, r_flux


# =========================================================
#  ACCRETION DISK MODEL (MIP turbulence)
# =========================================================

def run_disk_sim(
    g1=-1e-6,
    alpha_flux=2.0,
    N_rings=64,
    r_in_factor=1.3,   # initial inner edge in units of r_flux(0)
    r_out=15.0,
    nu_visc=0.02,      # viscous drift strength
    qpo_amp=0.3,       # amplitude of periodic modulation of viscosity
    qpo_period=5.0,    # period of modulation (in time units)
):
    """
    MIP accretion disk model around a flux-horizon black hole.

    - Disk represented by N_rings circular rings.
    - Each ring has radius r_i(t), azimuthal angle phi_i(t).
    - Keplerian shear Omega_i ~ (G_eff M / r^3)^{1/2}.
    - Radial drift from viscosity pulls rings slowly toward flux horizon.
    - MIP turbulence at each ring:
          MIP_i ∝ |dΩ/dr| * flux_fraction * (r_flux / r_i)^2
      (shear × global flux × compression factor)
    - "Luminosity" L(t) ~ sum of inner-ring MIP_i.
    """

    n_steps = int(T_MAX / DT)

    # --- arrays to record global behaviour ---
    t_vals   = np.zeros(n_steps)
    flux_vals = np.zeros(n_steps)
    G_vals   = np.zeros(n_steps)
    c_vals   = np.zeros(n_steps)
    rs_vals  = np.zeros(n_steps)
    rh_vals  = np.zeros(n_steps)
    L_vals   = np.zeros(n_steps)      # light curve from inner disk

    # --- cosmology state ---
    a = a0
    S = S0

    # --- ring radii & angles ---
    # start inner edge slightly outside initial flux horizon
    # (we'll set this properly after first cosmo step)
    r_rings = None
    phi_rings = None

    # for diagnostics: store final MIP profile
    final_MIP = None
    final_radii = None

    for n in range(n_steps):
        t = n * DT
        t_vals[n] = t

        # ---- global flux step ----
        a, S, H, flux_frac, G_eff, c_eff = step_cosmology(a, S, DT, g1)
        flux_vals[n] = flux_frac
        G_vals[n]    = G_eff
        c_vals[n]    = c_eff

        rs_GR, r_flux = compute_horizons(flux_frac, G_eff, c_eff, alpha_flux=alpha_flux)
        rs_vals[n] = rs_GR
        rh_vals[n] = r_flux

        # ---- initialise disk once we know first r_flux ----
        if r_rings is None:
            r_in0 = r_in_factor * r_flux
            # log-spaced rings from inner edge to outer edge
            r_rings = np.exp(np.linspace(math.log(r_in0), math.log(r_out), N_rings))
            phi_rings = 2.0 * math.pi * np.random.rand(N_rings)

        # ---- compute local orbital frequencies & shear ----
        # Omega_i ~= sqrt(G M / r^3)
        Omega = np.sqrt(G_eff * M_BH / (r_rings**3 + 1e-12))

        # numerical derivative dOmega/dr for each ring
        # central differences internally, forward/backward at edges
        dOdr = np.zeros_like(Omega)
        dr   = np.gradient(r_rings)
        dOdr[1:-1] = (Omega[2:] - Omega[:-2]) / (r_rings[2:] - r_rings[:-2])
        dOdr[0]    = (Omega[1]  - Omega[0])  / (r_rings[1]  - r_rings[0])
        dOdr[-1]   = (Omega[-1] - Omega[-2]) / (r_rings[-1] - r_rings[-2])

        # ---- viscous radial drift ----
        # base viscosity
        nu = nu_visc

        # optional QPO-like modulation of viscosity
        # (drives quasi-periodic oscillations in L(t))
        if qpo_amp != 0.0 and qpo_period > 0.0:
            mod = 1.0 + qpo_amp * math.sin(2.0 * math.pi * t / qpo_period)
            nu *= mod

        # simple inward drift law: pull rings toward r_flux
        # more strongly if they are close to inner edge
        dr_dt = -nu * (r_rings - r_flux) / (r_rings + 1e-6)

        # never let rings cross the flux horizon
        r_new = r_rings + dr_dt * DT
        r_new = np.maximum(r_new, 1.05 * r_flux)

        # ---- update azimuthal angle (Keplerian rotation) ----
        phi_rings = phi_rings + Omega * DT
        phi_rings = np.mod(phi_rings, 2.0 * math.pi)

        r_rings = r_new

        # ---- MIP turbulence field ----
        # Shear magnitude * flux fraction * compression factor
        shear_mag = np.abs(dOdr)
        compression = (r_flux / r_rings)**2
        MIP_turb = shear_mag * flux_frac * compression

        # store final profile for later plotting
        if n == n_steps - 1:
            final_MIP   = MIP_turb.copy()
            final_radii = r_rings.copy()

        # "Luminosity": sum of MIP_turb in inner half of disk
        inner_mask = r_rings < (r_flux + 3.0 * (r_out - r_flux) / 5.0)
        L_vals[n] = np.sum(MIP_turb[inner_mask])

    return {
        "t": t_vals,
        "flux": flux_vals,
        "G": G_vals,
        "c": c_vals,
        "rs": rs_vals,
        "rh": rh_vals,
        "L": L_vals,
        "final_radii": final_radii,
        "final_MIP": final_MIP,
    }


# =========================================================
#  MAIN: RUN DISK + MIP TURBULENCE MODEL
# =========================================================

if __name__ == "__main__":
    # parameters you can tweak
    g1         = -1e-6   # flux → G coupling
    alpha_flux = 2.0     # how much flux inflates the horizon
    N_rings    = 64
    r_out      = 15.0

    result = run_disk_sim(
        g1=g1,
        alpha_flux=alpha_flux,
        N_rings=N_rings,
        r_out=r_out,
        nu_visc=0.02,
        qpo_amp=0.4,    # try 0 to turn off, 0.4 for chunky QPO
        qpo_period=5.0,
    )

    t   = result["t"]
    L   = result["L"]

    # ----------------- PLOTS -----------------

    # 1) Light curve (this is where QPOs live)
    plt.figure()
    plt.plot(t, L)
    plt.xlabel("t (ticks)")
    plt.ylabel("Relative luminosity L(t)")
    plt.title("MIP-driven disk luminosity (inner disk)")

    # 2) Final MIP turbulence vs radius
    plt.figure()
    plt.plot(result["final_radii"], result["final_MIP"])
    plt.xlabel("radius r")
    plt.ylabel("MIP turbulence level")
    plt.title("Final radial MIP turbulence profile")
    plt.axvline(result["rh"][-1], color="k", linestyle="--", label="flux horizon")
    plt.legend()

    # 3) Flux fraction vs time
    plt.figure()
    plt.plot(t, result["flux"])
    plt.xlabel("t (ticks)")
    plt.ylabel("flux_fraction(t)")
    plt.title("Global flux fraction vs time")

    # 4) Flux horizon vs GR horizon vs time
    plt.figure()
    plt.plot(t, result["rs"], label="r_s^GR(t)")
    plt.plot(t, result["rh"], label="r_h^flux(t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("radius")
    plt.title("GR vs flux-defined horizon radii")
    plt.legend()

    # 5) Effective G and c
    plt.figure()
    plt.plot(t, result["G"], label="G_eff(t)")
    plt.plot(t, result["c"], label="c_eff(t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("G_eff, c_eff")
    plt.title("Effective G and c from flux")
    plt.legend()

    plt.show()
