import math
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
#  GLOBAL COSMOLOGY / FLUX PARAMETERS  (same family as before)
# =========================================================

T_MAX = 20.0       # total time (ticks)
DT    = 0.01
N_STEPS = int(T_MAX / DT)

s0      = 1.0
gamma   = 0.01
Omega_m = 0.3
Omega_r = 0.0
k_flux  = 0.05

G0 = 1.0
c0 = 1.0
H0 = 1.0

# initial cosmology
a0 = 0.1
S0 = 0.0

# black hole mass (dimensionless, toy units)
M_BH = 1.0


# =========================================================
#  FLUX COSMOLOGY STEPPER (your standard one)
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
#  HORIZON RADII (GR vs flux horizon)
# =========================================================

def compute_horizons(flux_frac, G_eff, c_eff, alpha_flux=2.0):
    """
    Return:
      r_s^GR : GR Schwarzschild radius
      r_h    : flux-defined "MIP horizon"
    """
    rs_GR = 2.0 * G_eff * M_BH / (c_eff**2 + 1e-12)
    r_flux = rs_GR * (1.0 + alpha_flux * flux_frac)
    return rs_GR, r_flux


# =========================================================
#  MIP-DRIVEN JET SIMULATOR (1D vertical column)
# =========================================================

def run_jet_sim(
    g1=-1e-6,
    alpha_flux=2.0,
    z_max=20.0,
    N_z=200,
    horizon_strength=0.7,  # how strongly horizon compression suppresses phi at base
    z_chimney=4.0,         # vertical scale for flux recovery (chimney height scale)
    tau_phi=3.0,           # relaxation timescale for phi -> equilibrium
    accel_strength=10.0,   # K in dv/dt = K * phi_z / phi^2  (MIP pressure)
    v_damp=0.05            # small damping to keep v finite
):
    """
    Simulate a 1D vertical jet along the rotation axis of the black hole.

    Grid: z in [0, z_max], with N_z points.
    State fields:
      phi(z,t): vertical flux fraction (MIP availability)
      v(z,t):   jet velocity along +z (dimensionless, capped by c_eff)

    Physics:
      - Global flux fraction phi_global(t) from your cosmology.
      - Horizon compression suppresses phi at base:
            phi_base(t) = phi_global(t) * (1 - horizon_strength * f_comp(r_h))
      - Flux recovers with height to phi_global(t) on scale z_chimney.
      - Jet acceleration:
            dv/dt = K * (dphi/dz) / phi^2  - v_damp * v
      - v is capped at 0.99 * c_eff(t).
    """

    # vertical grid
    z = np.linspace(0.0, z_max, N_z)
    dz = z[1] - z[0]

    # histories
    t_vals    = np.zeros(N_STEPS)
    flux_vals = np.zeros(N_STEPS)
    G_vals    = np.zeros(N_STEPS)
    c_vals    = np.zeros(N_STEPS)
    rs_vals   = np.zeros(N_STEPS)
    rh_vals   = np.zeros(N_STEPS)
    L_vals    = np.zeros(N_STEPS)      # jet luminosity proxy

    phi_hist  = np.zeros((N_STEPS, N_z))
    v_hist    = np.zeros((N_STEPS, N_z))

    # initial cosmology
    a = a0
    S = S0

    # initial fields
    phi = np.ones(N_z)         # start with uniform flux
    v   = np.zeros(N_z)        # no jet initially

    for n in range(N_STEPS):
        t = n * DT
        t_vals[n] = t

        # ---- global cosmology tick ----
        a, S, H, flux_global, G_eff, c_eff = step_cosmology(a, S, DT, g1)
        flux_vals[n] = flux_global
        G_vals[n]    = G_eff
        c_vals[n]    = c_eff

        rs_GR, r_flux = compute_horizons(flux_global, G_eff, c_eff, alpha_flux=alpha_flux)
        rs_vals[n] = rs_GR
        rh_vals[n] = r_flux

        # ---- horizon compression → base flux deficit ----
        # compression factor grows with r_h (saturates as r_h >> 1)
        comp_factor = r_flux / (r_flux + 1.0)      # in (0,1)
        # base flux at z=0 is reduced relative to global flux
        phi_base = flux_global * (1.0 - horizon_strength * comp_factor)
        phi_base = max(phi_base, 0.03)            # don't let it hit 0

        # flux at large height tends to global value
        phi_far = max(flux_global, phi_base + 0.01)

        # equilibrium vertical profile (MIP chimney):
        # low flux at base, rising toward phi_far with height
        phi_eq = phi_base + (phi_far - phi_base) * (1.0 - np.exp(-z / z_chimney))

        # relax phi toward this profile
        dphi_dt = (phi_eq - phi) / tau_phi
        phi = phi + dphi_dt * DT
        phi = np.clip(phi, 0.01, 1.5)

        # ---- jet acceleration from flux gradient ----
        phi_z = np.gradient(phi, dz)
        dv_dt = accel_strength * phi_z / (phi**2 + 1e-12)

        # damping
        dv_dt -= v_damp * v

        v = v + dv_dt * DT
        v = np.maximum(v, 0.0)     # no downward jet

        # local speed-of-light cap from c_eff(t)
        v_max_allowed = 0.99 * c_eff
        v = np.minimum(v, v_max_allowed)

        # ---- jet luminosity proxy ----
        # take L ~ ∫ v^3 dz (up to a constant factor)
        L_vals[n] = np.sum(v**3) * dz

        # store histories
        phi_hist[n, :] = phi
        v_hist[n, :]   = v

    return {
        "t": t_vals,
        "z": z,
        "phi_hist": phi_hist,
        "v_hist": v_hist,
        "flux": flux_vals,
        "G": G_vals,
        "c": c_vals,
        "rs": rs_vals,
        "rh": rh_vals,
        "L": L_vals,
    }


# =========================================================
#  MAIN: RUN JET SIM + PLOTS
# =========================================================

if __name__ == "__main__":
    g1         = -1e-6
    alpha_flux = 2.0

    result = run_jet_sim(
        g1=g1,
        alpha_flux=alpha_flux,
        z_max=20.0,
        N_z=200,
        horizon_strength=0.7,
        z_chimney=4.0,
        tau_phi=3.0,
        accel_strength=10.0,
        v_damp=0.05,
    )

    t  = result["t"]
    z  = result["z"]
    L  = result["L"]
    phi_hist = result["phi_hist"]
    v_hist   = result["v_hist"]

    # --------- 1) Jet luminosity L(t) ----------
    plt.figure()
    plt.plot(t, L)
    plt.xlabel("t (ticks)")
    plt.ylabel("Jet luminosity proxy L(t)")
    plt.title("MIP-driven jet luminosity")

    # --------- 2) Final jet velocity profile ----------
    plt.figure()
    plt.plot(z, v_hist[-1, :])
    plt.xlabel("z (height above BH)")
    plt.ylabel("v(z)")
    plt.title("Final jet velocity profile")

    # --------- 3) Final vertical flux fraction profile ----------
    plt.figure()
    plt.plot(z, phi_hist[-1, :])
    plt.xlabel("z (height above BH)")
    plt.ylabel("phi(z)")
    plt.title("Final vertical flux fraction profile")

    # --------- 4) Jet velocity history v(z,t) ----------
    plt.figure()
    plt.imshow(
        v_hist.T,
        origin="lower",
        aspect="auto",
        extent=[t[0], t[-1], z[0], z[-1]],
    )
    plt.colorbar(label="v(z,t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("z (height)")
    plt.title("Jet velocity history v(z,t)")

    # --------- 5) Effective G_eff and c_eff vs time ----------
    plt.figure()
    plt.plot(t, result["G"], label="G_eff(t)")
    plt.plot(t, result["c"], label="c_eff(t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("G_eff, c_eff")
    plt.title("Effective G and c from flux cosmology")
    plt.legend()

    plt.tight_layout()
    plt.show()

