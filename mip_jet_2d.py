"""
mip_jet_2d.py

High-resolution 2D axisymmetric jet model driven by your MIP / flux cosmology idea.

- Global flux field evolves via simple entropy -> flux_fraction(t) model
- G_eff(t) and c_eff(t) are derived from flux_fraction(t)
- A vertical jet is launched from the z=0 boundary by a flux pressure gradient
- We evolve v_z(z,r,t), v_r(z,r,t), and local flux fraction phi(z,r,t)
- Outputs:
    * Jet luminosity vs time
    * Final v_z(z,r) map (jet structure)
    * Final phi(z,r) map (flux / MIP density)
    * Final gamma(z) along axis (Lorentz factor profile)
"""

import math
import numpy as np
import matplotlib.pyplot as plt

# ================================================================
#  GLOBAL PARAMETERS
# ================================================================

# --- Grid / resolution ---
NZ = 400          # number of cells in z (height)
NR = 200          # number of cells in r (radius)

Z_MAX = 20.0      # total jet height (code units)
R_MAX = 6.0       # maximum radius simulated

dz = Z_MAX / (NZ - 1)
dr = R_MAX / (NR - 1)

# --- Time stepping ---
T_MAX = 10.0      # total simulation time (code units ~ "ticks")
DT    = 0.003     # time step; CFL-safe for v ~ c_eff ~ 1

N_STEPS = int(T_MAX / DT)

# --- Flux cosmology parameters (global) ---
s0      = 1.0
gamma_S = 0.02        # entropy production rate
Omega_m = 0.3
Omega_r = 0.0
k_flux  = 0.05

G0 = 1.0
c0 = 1.0
H0 = 1.0

a0 = 0.1
S0 = 0.0

# Coupling of flux to gravity / light speed
g1_G = -5.0    # G_eff = G0 * (1 + g1_G * (1 - flux_fraction))
g1_c = -3.0    # c_eff = c0 * (1 + g1_c * (1 - flux_fraction))

# --- Jet / MIP parameters ---
PHI_BASE      = 0.8     # base flux fraction injected in core
PHI_AMBIENT   = 0.4     # ambient flux fraction outside jet
K_P           = 2.0     # pressure coefficient: P_MIP = K_P * phi
JET_RADIUS    = 1.0     # radius of jet nozzle at z=0
JET_ACCEL_GAIN = 3.0    # how strongly local flux gradient accelerates jet

# Small floor values for stability
EPS = 1e-8


# ================================================================
#  COSMOLOGY STEPPER (GLOBAL FLUX FRACTION)
# ================================================================

def step_cosmology(a, S, dt):
    """
    Minimal FRW-like + flux model to give us a(t), flux_fraction(t),
    and thus G_eff(t), c_eff(t).
    """
    # avoid zero scale factor
    if a <= 1e-8:
        a = 1e-8

    S_max = s0 * a * a
    Phi   = max(S_max - S, 0.0)
    flux_fraction = 0.0 if S_max == 0 else Phi / S_max

    # entropy production
    dS_dt = gamma_S * Phi
    S_new = S + dS_dt * dt

    # Hubble parameter
    a_safe = max(a, 1e-4)
    H_sq = Omega_m / (a_safe ** 3) + Omega_r / (a_safe ** 4) + k_flux * flux_fraction
    H_sq = max(H_sq, 0.0)
    H    = math.sqrt(H_sq)

    # scale factor evolution
    a_new = a + H * a * dt

    # effective constants
    # flux_fraction ≈ 1 at early times, then gently decreases
    f = flux_fraction
    # we define (1 - f) as "flux hardness"; when flux is exhausted, f→0, hardness→1
    hardness = 1.0 - f

    G_eff = G0 * (1.0 + g1_G * hardness)
    c_eff = c0 * (1.0 + g1_c * hardness)

    # keep them positive
    G_eff = max(G_eff, 0.1 * G0)
    c_eff = max(c_eff, 0.3 * c0)

    return a_new, S_new, flux_fraction, G_eff, c_eff, H


# ================================================================
#  UTILITY: INITIAL CONDITIONS
# ================================================================

def initialize_fields():
    """
    Allocate and initialize v_z, v_r, phi.

    Convention:
      - z index: 0 ... NZ-1
      - r index: 0 ... NR-1
    """
    v_z  = np.zeros((NZ, NR), dtype=np.float64)
    v_r  = np.zeros((NZ, NR), dtype=np.float64)

    # flux field starts with ambient + small bump near the base
    phi  = PHI_AMBIENT * np.ones((NZ, NR), dtype=np.float64)

    # Add a smooth flux core near the axis for z ~ 0..1
    z_grid = np.linspace(0.0, Z_MAX, NZ)
    r_grid = np.linspace(0.0, R_MAX, NR)
    Z, R = np.meshgrid(z_grid, r_grid, indexing="ij")

    core_mask = (Z < 1.0) & (R < JET_RADIUS)
    phi[core_mask] = PHI_BASE

    return v_z, v_r, phi


# ================================================================
#  SIMPLE FINITE DIFFERENCE OPERATORS
# ================================================================

def grad_z(field, dz):
    """
    ∂field/∂z using central differences in the interior,
    one-sided at boundaries.
    """
    d = np.zeros_like(field)
    d[1:-1, :] = (field[2:, :] - field[:-2, :]) / (2.0 * dz)
    d[0,  :]   = (field[1,  :] - field[0,  :]) / dz
    d[-1, :]   = (field[-1, :] - field[-2, :]) / dz
    return d

def grad_r(field, dr):
    """
    ∂field/∂r using central differences in the interior,
    one-sided at boundaries.
    Axis (r=0) uses symmetric boundary condition.
    """
    d = np.zeros_like(field)
    # interior
    d[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / (2.0 * dr)
    # axis r=0: symmetric ∂/∂r ≈ (f1 - f0)/dr but enforce v_r=0 at axis
    d[:, 0] = (field[:, 1] - field[:, 0]) / dr
    # outer boundary
    d[:, -1] = (field[:, -1] - field[:, -2]) / dr
    return d

def divergence(v_z, v_r, dz, dr):
    """
    ∇·v in cylindrical coordinates (approx):

    ∇·v ≈ ∂v_z/∂z + ∂v_r/∂r  (we ignore the explicit (1/r)*v_r term
    for simplicity; can be added later if needed)
    """
    dvz_dz = grad_z(v_z, dz)
    dvr_dr = grad_r(v_r, dr)
    return dvz_dz + dvr_dr


# ================================================================
#  MAIN SIMULATOR
# ================================================================

def run_simulation():
    # allocate fields
    v_z, v_r, phi = initialize_fields()

    # global cosmology
    a = a0
    S = S0

    # preallocate diagnostics
    t_arr    = np.zeros(N_STEPS)
    lum_arr  = np.zeros(N_STEPS)
    f_arr    = np.zeros(N_STEPS)  # global flux_fraction
    G_arr    = np.zeros(N_STEPS)
    c_arr    = np.zeros(N_STEPS)

    # coordinates
    z_grid = np.linspace(0.0, Z_MAX, NZ)
    r_grid = np.linspace(0.0, R_MAX, NR)
    Z, R = np.meshgrid(z_grid, r_grid, indexing="ij")

    # jet core mask at base (z=0)
    jet_core = (R[0, :] <= JET_RADIUS)

    for n in range(N_STEPS):
        t = n * DT
        t_arr[n] = t

        # --- 1. cosmology step (global flux & constants) ---
        a, S, flux_frac, G_eff, c_eff, H = step_cosmology(a, S, DT)
        f_arr[n] = flux_frac
        G_arr[n] = G_eff
        c_arr[n] = c_eff

        # --- 2. local pressure from flux field ---
        P = K_P * phi

        # gradients
        dP_dz = grad_z(P, dz)
        dP_dr = grad_r(P, dr)

        # --- 3. update velocities (explicit Euler) ---
        # dv/dt = - grad(P) + MIP acceleration term along z
        # accel_z ~ JET_ACCEL_GAIN * (local hardness) near axis
        hardness_local = 1.0 - phi   # cold / hardened flux
        accel_MIP = JET_ACCEL_GAIN * hardness_local

        v_z += DT * (-dP_dz + accel_MIP)
        v_r += DT * (-dP_dr)

        # --- 4. enforce relativistic speed limit with local c_eff ---
        v2 = v_z**2 + v_r**2
        vmax2 = (0.99 * c_eff)**2
        mask_superluminal = v2 > vmax2
        if np.any(mask_superluminal):
            factor = np.sqrt(vmax2 / (v2[mask_superluminal] + EPS))
            v_z[mask_superluminal] *= factor
            v_r[mask_superluminal] *= factor

        # --- 5. update flux field phi with advection + compression ---
        # dphi/dt = - (v · ∇phi) - phi (∇·v)
        dphi_dz = grad_z(phi, dz)
        dphi_dr = grad_r(phi, dr)
        div_v   = divergence(v_z, v_r, dz, dr)

        adv_term = v_z * dphi_dz + v_r * dphi_dr
        comp_term = phi * div_v

        phi += DT * ( -adv_term - comp_term )

        # keep phi in [0, 1.5] range for stability (can be hardened > 1 slightly)
        phi = np.clip(phi, 0.0, 1.5)

        # --- 6. boundary conditions ---

        # axis r = 0: no radial flow, symmetric fields
        v_r[:, 0] = 0.0

        # outer radius: zero-gradient (open)
        v_z[:, -1] = v_z[:, -2]
        v_r[:, -1] = v_r[:, -2]
        phi[:, -1] = phi[:, -2]

        # top boundary (z = Z_MAX): open
        v_z[-1, :] = v_z[-2, :]
        v_r[-1, :] = v_r[-2, :]
        phi[-1, :] = phi[-2, :]

        # base (z=0): inject jet inside JET_RADIUS
        # flux gets refreshed from global flux fraction & base value
        phi_base_now = PHI_BASE * flux_frac
        phi[0, :] = PHI_AMBIENT
        phi[0, jet_core] = phi_base_now

        # vertical velocity injection close to c_eff
        v_inj = 0.8 * c_eff
        v_z[0, jet_core] = v_inj
        # no inflow from outside region; let that be ambient
        v_z[0, ~jet_core] *= 0.5

        # --- 7. compute jet luminosity proxy ---
        # L(t) ~ integral over area of phi * v_z at some sampling height
        # take a slice at z_slice ~ 0.75 * Z_MAX
        z_slice_index = int(0.75 * (NZ - 1))
        phi_slice = phi[z_slice_index, :]
        vz_slice  = v_z[z_slice_index, :]
        # cylindrical ring integration: sum phi * v_z * 2π r dr
        L = 0.0
        for j in range(NR):
            r_mid = r_grid[j]
            area_ring = 2.0 * math.pi * r_mid * dr
            L += phi_slice[j] * vz_slice[j] * area_ring
        lum_arr[n] = L

        # (optional) basic progress print every few hundred steps
        if (n % 500) == 0 or n == N_STEPS - 1:
            print(f"Step {n+1}/{N_STEPS}, t={t:.3f}, flux_frac={flux_frac:.3f}, G_eff={G_eff:.3f}, c_eff={c_eff:.3f}")

    return {
        "t": t_arr,
        "lum": lum_arr,
        "flux_frac": f_arr,
        "G_eff": G_arr,
        "c_eff": c_arr,
        "v_z": v_z,
        "v_r": v_r,
        "phi": phi,
        "z_grid": z_grid,
        "r_grid": r_grid,
    }


# ================================================================
#  POST-PROCESSING / PLOTTING
# ================================================================

def make_plots(result):
    t       = result["t"]
    lum     = result["lum"]
    f_arr   = result["flux_frac"]
    G_arr   = result["G_eff"]
    c_arr   = result["c_eff"]
    v_z     = result["v_z"]
    v_r     = result["v_r"]
    phi     = result["phi"]
    z_grid  = result["z_grid"]
    r_grid  = result["r_grid"]

    Z, R = np.meshgrid(z_grid, r_grid, indexing="ij")

    # --- 1. Jet luminosity vs time ---
    plt.figure(figsize=(6, 4))
    plt.plot(t, lum)
    plt.xlabel("t (code units)")
    plt.ylabel("Jet luminosity proxy L(t)")
    plt.title("MIP-driven jet luminosity vs time")
    plt.grid(True)

    # --- 2. Final v_z map ---
    plt.figure(figsize=(6, 8))
    plt.imshow(
        v_z,
        origin="lower",
        extent=[0, R_MAX, 0, Z_MAX],
        aspect="auto",
    )
    plt.colorbar(label="v_z")
    plt.xlabel("r")
    plt.ylabel("z")
    plt.title("Final vertical velocity v_z(z,r)")

    # --- 3. Final phi map (flux field) ---
    plt.figure(figsize=(6, 8))
    plt.imshow(
        phi,
        origin="lower",
        extent=[0, R_MAX, 0, Z_MAX],
        aspect="auto",
    )
    plt.colorbar(label="phi (flux / MIP density)")
    plt.xlabel("r")
    plt.ylabel("z")
    plt.title("Final flux field phi(z,r)")

    # --- 4. Final gamma along axis ---
    # Take r ≈ 0 column
    v_axis = v_z[:, 0]   # near axis
    # Use final c_eff as local light speed for gamma
    c_eff_final = c_arr[-1]
    beta_axis = v_axis / (c_eff_final + EPS)
    beta_axis = np.clip(beta_axis, -0.999, 0.999)
    gamma_axis = 1.0 / np.sqrt(1.0 - beta_axis**2)

    plt.figure(figsize=(6, 4))
    plt.plot(z_grid, gamma_axis)
    plt.xlabel("z")
    plt.ylabel("gamma(z)")
    plt.title("Final Lorentz factor along jet axis")
    plt.grid(True)

    # --- 5. Global flux fraction and constants vs time ---
    plt.figure(figsize=(6, 4))
    plt.plot(t, f_arr, label="flux_fraction(t)")
    plt.xlabel("t")
    plt.ylabel("flux_fraction")
    plt.title("Global flux fraction vs time")
    plt.grid(True)

    plt.figure(figsize=(6, 4))
    plt.plot(t, G_arr, label="G_eff(t)")
    plt.plot(t, c_arr, label="c_eff(t)")
    plt.xlabel("t")
    plt.ylabel("Effective constants")
    plt.legend()
    plt.title("G_eff and c_eff from flux cosmology")
    plt.grid(True)

    plt.tight_layout()
    plt.show()


# ================================================================
#  MAIN
# ================================================================

if __name__ == "__main__":
    result = run_simulation()
    make_plots(result)
