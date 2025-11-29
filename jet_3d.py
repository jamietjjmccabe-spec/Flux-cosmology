import math
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
#  FLUX COSMOLOGY PARAMETERS
# ============================================================

T_MAX = 10.0       # total "cosmic" time (arbitrary units, e.g. years)
DT    = 0.05       # cosmology / jet time step
N_STEPS = int(T_MAX / DT)

s0      = 1.0
gamma_S = 0.01
Omega_m = 0.3
Omega_r = 0.0
k_flux  = 0.05

G0 = 1.0
c0 = 1.0
H0 = 1.0

a0 = 0.1   # initial scale factor
S0 = 0.0   # initial entropy


def step_cosmology(a, S, dt, g1):
    """
    Simple FRW + flux exhaustion + variable G and c.
    Returns: a_new, S_new, H, flux_fraction, G_eff, c_eff
    """
    if a <= 1e-8:
        a = 1e-8

    S_max = s0 * a * a
    Phi   = max(S_max - S, 0.0)
    flux_fraction = 0 if S_max == 0 else Phi / S_max

    # entropy
    dS_dt = gamma_S * Phi
    S_new = S + dS_dt * dt

    # Hubble parameter
    a_safe = max(a, 1e-4)
    H_sq = Omega_m / (a_safe**3) + Omega_r / (a_safe**4) + k_flux * flux_fraction
    H_sq = max(H_sq, 0.0)
    H    = math.sqrt(H_sq)

    # scale factor
    a_new = a + (H * a) * dt

    # variable gravity
    G_eff = G0 * (1.0 + g1 * flux_fraction)

    # variable c
    if H0 != 0:
        c_eff = c0 * (1.0 + 0.3 * (H / H0 - 1.0))
    else:
        c_eff = c0

    return a_new, S_new, H, flux_fraction, G_eff, c_eff


# ============================================================
#  3D JET GRID SETUP
# ============================================================

# Grid dimensions (keep modest so it runs quickly)
NX, NY, NZ = 32, 32, 48

# Physical spacing (arbitrary units)
DX = DY = DZ = 1.0

# 3D coordinate grids for convenience (not strictly needed everywhere)
x = (np.arange(NX) - (NX-1)/2) * DX
y = (np.arange(NY) - (NY-1)/2) * DY
z = np.arange(NZ) * DZ

X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

# Jet axis is at (x=0, y=0)
R = np.sqrt(X**2 + Y**2)


# ============================================================
#  INITIAL CONDITIONS (phi and vz in 3D)
# ============================================================

def initialize_jet_fields():
    """
    Initial flux density phi(x,y,z) and vertical velocity vz(x,y,z).
    Start with a concentrated "nozzle" at the bottom around the axis.
    """
    # Flux concentrated near base and axis
    r0 = 2.0
    z0 = 1.0

    phi0 = np.exp(-(R / r0)**2) * np.exp(-(Z / z0)**2)

    # Initial vertical velocity: small upward flow near base
    vz0 = 0.05 * np.exp(-(R / r0)**2) * np.exp(-(Z / z0)**2)

    return phi0, vz0


# ============================================================
#  3D JET EVOLUTION
# ============================================================

def evolve_jet_3d(g1=-1e-6):
    """
    Full 3D jet + cosmology evolution.

    g1 < 0   → G_eff increases as flux hardens (your "mass from flux" behaviour).
    """
    # Allocate cosmology history
    times      = np.zeros(N_STEPS)
    flux_hist  = np.zeros(N_STEPS)
    G_hist     = np.zeros(N_STEPS)
    c_hist     = np.zeros(N_STEPS)
    H_hist     = np.zeros(N_STEPS)
    L_hist     = np.zeros(N_STEPS)   # jet luminosity proxy

    # Initial cosmology
    a = a0
    S = S0

    # Jet fields
    phi, vz = initialize_jet_fields()

    # Constants for jet evolution
    D_phi       = 0.02     # diffusion coefficient for phi (MIP diffusion)
    decay_phi   = 0.05     # decay rate of phi
    accel_coeff = 0.8      # how strongly grad(phi) accelerates vz
    drag_v      = 0.2      # velocity damping
    injection_radius = 3.0 # radius at base where fresh flux is injected

    # Mask for base injection region
    base_mask = (Z == 0) & (R <= injection_radius)

    for n in range(N_STEPS):
        t = n * DT
        times[n] = t

        # ---- Cosmology step ----
        a, S, H, flux_frac, G_eff, c_eff = step_cosmology(a, S, DT, g1)

        flux_hist[n] = flux_frac
        G_hist[n]    = G_eff
        c_hist[n]    = c_eff
        H_hist[n]    = H

        # ---- 3D phi evolution: diffusion + decay + injection at the base ----
        # 3D Laplacian using periodic rolls, then we'll "fix" boundaries to Neumann
        lap_phi = (
            np.roll(phi, +1, axis=0) + np.roll(phi, -1, axis=0) +
            np.roll(phi, +1, axis=1) + np.roll(phi, -1, axis=1) +
            np.roll(phi, +1, axis=2) + np.roll(phi, -1, axis=2) -
            6.0 * phi
        ) / (DX * DY * DZ)  # overall scale isn't super important here

        # Source term: flux injection scales with global flux fraction
        source = 0.5 * flux_frac * base_mask.astype(float)

        phi = phi + DT * (D_phi * lap_phi - decay_phi * phi + source)

        # Enforce non-negative phi
        phi = np.clip(phi, 0.0, None)

        # Neumann-like boundary: copy near-boundary values
        phi[0, :, :]   = phi[1, :, :]
        phi[-1, :, :]  = phi[-2, :, :]
        phi[:, 0, :]   = phi[:, 1, :]
        phi[:, -1, :]  = phi[:, -2, :]
        phi[:, :, 0]   = phi[:, :, 1]
        phi[:, :, -1]  = phi[:, :, -2]

        # ---- vz evolution: acceleration from vertical gradient of phi ----
        # grad_z(phi) ~ forward difference
        grad_phi_z = (np.roll(phi, -1, axis=2) - phi) / DZ

        vz = vz + DT * (accel_coeff * grad_phi_z - drag_v * vz)

        # No downward flow in this simple model
        vz = np.clip(vz, 0.0, None)

        # Limit to subluminal speed relative to current c_eff
        vmax = 0.95 * c_eff
        vz = np.clip(vz, 0.0, vmax)

        # ---- Jet luminosity proxy L(t) ----
        # Simple proxy: integral of phi * vz over the column
        # Only count upward flow
        dV = DX * DY * DZ
        L_hist[n] = np.sum(phi * vz) * dV

    # Final state diagnostics
    # Axis index (closest to x=0,y=0)
    ix0 = NX // 2
    iy0 = NY // 2

    vz_axis    = vz[ix0, iy0, :]
    phi_axis   = phi[ix0, iy0, :]

    # Lorentz factor along axis using final c_eff
    c_final = c_hist[-1]
    beta    = np.clip(vz_axis / max(c_final, 1e-6), 0.0, 0.99)
    gamma   = 1.0 / np.sqrt(1.0 - beta**2)

    # 2D slice of phi through y = 0 plane
    phi_slice = phi[:, iy0, :]

    # Vertical profile: average phi over x,y at each z
    phi_z_profile = np.mean(phi, axis=(0, 1))

    results = {
        "t": times,
        "flux_frac": flux_hist,
        "G_eff": G_hist,
        "c_eff": c_hist,
        "H": H_hist,
        "L": L_hist,
        "z": z,
        "vz_axis": vz_axis,
        "phi_axis": phi_axis,
        "gamma_axis": gamma,
        "phi_slice": phi_slice,
        "phi_z_profile": phi_z_profile,
        "phi_3d": phi,
        "vz_3d": vz,
    }

    return results


# ============================================================
#  MAIN: RUN AND PLOT
# ============================================================

if __name__ == "__main__":
    # Tune g1 if you want stronger/weaker coupling between flux and gravity
    g1 = -1e-6

    data = evolve_jet_3d(g1=g1)

    t  = data["t"]
    z  = data["z"]

    # ---------- Figure 1: Jet luminosity vs time ----------
    plt.figure()
    plt.plot(t, data["L"])
    plt.xlabel("t (ticks)")
    plt.ylabel("Jet luminosity proxy L(t)")
    plt.title("MIP-driven jet luminosity vs time")

    # ---------- Figure 2: Final vertical velocity along jet axis ----------
    plt.figure()
    plt.plot(z, data["vz_axis"])
    plt.xlabel("z (height)")
    plt.ylabel("v_z(z)")
    plt.title("Final vertical velocity along jet axis")

    # ---------- Figure 3: Final flux field slice phi(x,z) at y=0 ----------
    plt.figure()
    plt.imshow(
        data["phi_slice"].T,
        origin="lower",
        extent=[x[0], x[-1], z[0], z[-1]],
        aspect="auto"
    )
    plt.colorbar(label="phi (flux / MIP density)")
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title("Final flux field φ(x,z) slice at y=0")

    # ---------- Figure 4: Lorentz factor along axis ----------
    plt.figure()
    plt.plot(z, data["gamma_axis"])
    plt.xlabel("z (height)")
    plt.ylabel("γ(z)")
    plt.title("Final Lorentz factor along jet axis")

    # ---------- Figure 5: G_eff and c_eff vs time ----------
    plt.figure()
    plt.plot(t, data["G_eff"], label="G_eff(t)")
    plt.plot(t, data["c_eff"], label="c_eff(t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("Effective constants")
    plt.title("G_eff and c_eff from flux cosmology")
    plt.legend()

    # ---------- Figure 6: Global flux fraction vs time ----------
    plt.figure()
    plt.plot(t, data["flux_frac"])
    plt.xlabel("t (ticks)")
    plt.ylabel("Φ / S_max")
    plt.title("Global flux fraction vs time")

    # ---------- Figure 7: Vertical flux fraction profile ----------
    plt.figure()
    plt.plot(z, data["phi_z_profile"])
    plt.xlabel("z (height)")
    plt.ylabel("⟨φ⟩(z)")
    plt.title("Final vertical flux profile (averaged over x,y)")

    plt.tight_layout()
    plt.show()
