import math
import numpy as np
import matplotlib.pyplot as plt

# ===============================================================
#  GLOBAL / COSMOLOGY PARAMETERS
# ===============================================================

T_MAX = 10.0       # total time (ticks / code-years)
DT    = 0.01       # time step

# simple closed-FRW + flux field toy model
s0      = 1.0
gamma_S = 0.02
Omega_m = 0.3
Omega_r = 0.0
k_flux  = 0.10

G0 = 1.0
c0 = 1.0
H0 = 1.0

a0 = 0.5          # initial scale factor
S0 = 0.0          # initial entropy

# gravity / c dependence on global flux fraction Φ
g1 = -0.10        # negative: more flux → slightly stronger gravity
c1 = -0.05        # more flux → slightly lower c (toy effect)

# ===============================================================
#  JET / GRID PARAMETERS
# ===============================================================

Z_MAX = 20.0      # vertical extent of jet
R_MAX = 6.0       # radial extent
NZ    = 160       # vertical grid points
NR    = 80        # radial grid points

z = np.linspace(0.0, Z_MAX, NZ)
r = np.linspace(0.0, R_MAX, NR)
dz = z[1] - z[0]
dr = r[1] - r[0]

Z, R = np.meshgrid(z, r, indexing="ij")  # shape (NZ, NR)

# flux field preferences
phi_core   = 0.95     # central high-flux column
phi_outer  = 0.40     # background flux level
R_core     = 1.0      # size of inner column
Z_base     = 2.0      # vertical scale for strong base injection

lambda_phi = 1.5      # relaxation rate of phi → phi_target
turb_amp   = 0.03     # tiny noise to break symmetry

# velocity evolution parameters
A_z      = 0.6        # strength of vertical acceleration from phi
A_r      = 0.4        # radial response to grad_r(phi)
A_phi    = 0.5        # spin-up from flux near axis
damp_z   = 0.30
damp_r   = 0.30
damp_phi = 0.30

vmax_frac = 0.95      # cap at 0.95 * c_eff


# ===============================================================
#  COSMOLOGY STEPPER
# ===============================================================

def step_cosmology(a, S, dt):
    """
    Minimal FRW + flux model:
    - S_max(a) = s0 * a^2
    - Φ = (S_max - S) / S_max
    - dS/dt = gamma_S * (S_max - S)
    - H^2 = Ω_m/a^3 + Ω_r/a^4 + k_flux * Φ
    """
    if a <= 1e-8:
        a = 1e-8

    S_max = s0 * a * a
    Phi   = max(S_max - S, 0.0)
    flux_fraction = 0.0 if S_max == 0.0 else Phi / S_max

    # entropy production
    dS_dt = gamma_S * Phi
    S_new = S + dS_dt * dt

    # Hubble parameter
    a_safe = max(a, 1e-4)
    H_sq = Omega_m / (a_safe**3) + Omega_r / (a_safe**4) + k_flux * flux_fraction
    H_sq = max(H_sq, 0.0)
    H = math.sqrt(H_sq)

    # scale factor evolution
    a_new = a + (H * a) * dt

    # flux-dependent effective G and c
    G_eff = G0 * (1.0 + g1 * flux_fraction)
    c_eff = c0 * (1.0 + c1 * flux_fraction)

    return a_new, S_new, flux_fraction, G_eff, c_eff, H


# ===============================================================
#  BUILD TARGET FLUX GEOMETRY
# ===============================================================

def phi_target_field(global_flux_fraction):
    """
    Desired quasi-static φ(r,z) profile for this global tick.
    High near base & axis, decaying with r and z, scaled by global Φ.
    """
    # base central column
    core_shape = np.exp(- (R / R_core)**2 - (Z / Z_base)**2)
    phi_target = phi_outer + (phi_core - phi_outer) * core_shape

    # modulate by global flux fraction (Φ → 0 shuts everything down)
    phi_target *= global_flux_fraction

    # hard bounds
    phi_target = np.clip(phi_target, 0.0, 1.0)
    return phi_target


# ===============================================================
#  MAIN JET SIMULATION
# ===============================================================

def run_full_jet_sim():
    n_steps = int(T_MAX / DT)

    # --- storage for diagnostics ---
    times        = []
    jet_L        = []   # luminosity proxy vs time
    flux_global  = []
    G_eff_hist   = []
    c_eff_hist   = []
    H_hist       = []

    # --- cosmology state ---
    a = a0
    S = S0

    # --- fields on (z, r) grid ---
    phi   = np.zeros((NZ, NR), dtype=float)
    v_z   = np.zeros_like(phi)
    v_r   = np.zeros_like(phi)
    v_phi = np.zeros_like(phi)

    # --- main time loop ---
    for step in range(n_steps):
        t = step * DT
        times.append(t)

        # ---- cosmology update ----
        a, S, flux_frac, G_eff, c_eff, H = step_cosmology(a, S, DT)

        flux_global.append(flux_frac)
        G_eff_hist.append(G_eff)
        c_eff_hist.append(c_eff)
        H_hist.append(H)

        # ---- update flux field φ(z,r) ----
        phi_tgt = phi_target_field(flux_frac)

        # simple relaxation + tiny noise near base for "turbulence"
        noise = turb_amp * np.random.normal(size=phi.shape) * np.exp(-Z / (0.5 * Z_base))
        dphi_dt = -lambda_phi * (phi - phi_tgt) + noise
        phi += dphi_dt * DT
        phi = np.clip(phi, 0.0, 1.0)

        # ---- gradients of φ (for forces) ----
        # gradient returns [∂/∂z, ∂/∂r] for axis 0=z,1=r
        dphi_dz, dphi_dr = np.gradient(phi, dz, dr, edge_order=2)

        # ---- update velocities ----
        # vertical acceleration: push from φ, plus slight response to grad_z
        a_z = A_z * phi - 0.2 * dphi_dz - damp_z * v_z

        # radial acceleration: respond to grad_r and centrifugal term
        eps_r = 0.05
        a_r = -A_r * dphi_dr - damp_r * v_r + (v_phi**2) / (R + eps_r)

        # azimuthal acceleration: spin up near axis from φ, damp with shear
        a_phi = (A_phi * phi / (R + eps_r)) - damp_phi * v_phi - (v_r * v_phi) / (R + eps_r)

        # integrate velocities
        v_z   += a_z   * DT
        v_r   += a_r   * DT
        v_phi += a_phi * DT

        # ---- cap speeds at < c_eff ----
        # c_eff is global so we can do one cap per grid.
        v_tot_sq = v_z**2 + v_r**2 + v_phi**2
        v_max_sq = (vmax_frac * c_eff)**2
        mask = v_tot_sq > v_max_sq
        if np.any(mask):
            scale = np.sqrt(v_max_sq / v_tot_sq[mask])
            v_z[mask]   *= scale
            v_r[mask]   *= scale
            v_phi[mask] *= scale

        # ---- luminosity proxy: base vertical energy flux ----
        # integrate over first few z cells
        z_base_cells = 3
        base_slice = slice(0, z_base_cells)
        L_t = np.sum(phi[base_slice, :] * np.abs(v_z[base_slice, :])) * dz * dr
        jet_L.append(L_t)

    # package useful final fields
    result = {
        "t":           np.array(times),
        "L":           np.array(jet_L),
        "flux_global": np.array(flux_global),
        "G_eff":       np.array(G_eff_hist),
        "c_eff":       np.array(c_eff_hist),
        "H":           np.array(H_hist),
        "phi_final":   phi,
        "v_z_final":   v_z,
        "v_r_final":   v_r,
        "v_phi_final": v_phi,
    }
    return result


# ===============================================================
#  DIAGNOSTIC PLOTS
# ===============================================================

def make_plots(result):
    t    = result["t"]
    L    = result["L"]
    phi  = result["phi_final"]
    v_z  = result["v_z_final"]
    v_r  = result["v_r_final"]
    v_ph = result["v_phi_final"]

    flux_global = result["flux_global"]
    G_eff_hist  = result["G_eff"]
    c_eff_hist  = result["c_eff"]

    # choose axis index near r=0 (avoid exact 0 for division)
    axis_idx = 1 if NR > 1 else 0

    v_z_axis  = v_z[:, axis_idx]
    v_r_axis  = v_r[:, axis_idx]
    v_phi_axis = v_ph[:, axis_idx]
    phi_axis  = phi[:, axis_idx]

    # use final c_eff for gamma (slowly varying anyway)
    c_eff_final = c_eff_hist[-1]
    v_axis_tot_sq = v_z_axis**2 + v_r_axis**2 + v_phi_axis**2
    beta_sq = np.clip(v_axis_tot_sq / (c_eff_final**2 + 1e-12), 0.0, 0.999999)
    gamma_axis = 1.0 / np.sqrt(1.0 - beta_sq)

    # ---- Figure 1: jet luminosity vs time ----
    plt.figure(figsize=(5, 4))
    plt.plot(t, L)
    plt.xlabel("t (ticks)")
    plt.ylabel("Jet luminosity proxy L(t)")
    plt.title("MIP-driven jet luminosity vs time")
    plt.grid(True)

    # ---- Figure 2: final vertical velocity profile along axis ----
    plt.figure(figsize=(5, 4))
    plt.plot(z, v_z_axis)
    plt.xlabel("z (height)")
    plt.ylabel("v_z(z)")
    plt.title("Final vertical velocity along jet axis")
    plt.grid(True)

    # ---- Figure 3: final flux field φ(z,r) as 2D map ----
    plt.figure(figsize=(5, 4))
    im = plt.imshow(
        phi,
        origin="lower",
        extent=[0.0, R_MAX, 0.0, Z_MAX],
        aspect="auto"
    )
    plt.colorbar(im, label="φ(z,r) (MIP density)")
    plt.xlabel("r")
    plt.ylabel("z")
    plt.title("Final flux field φ(z,r)")

    # ---- Figure 4: Lorentz factor along jet axis ----
    plt.figure(figsize=(5, 4))
    plt.plot(z, gamma_axis)
    plt.xlabel("z (height)")
    plt.ylabel("γ(z)")
    plt.title("Final Lorentz factor along jet axis")
    plt.grid(True)

    # ---- Figure 5: global flux fraction vs time ----
    plt.figure(figsize=(5, 4))
    plt.plot(t, flux_global)
    plt.xlabel("t (ticks)")
    plt.ylabel("Global flux fraction Φ(t)")
    plt.title("Global flux fraction vs time")
    plt.grid(True)

    # ---- Figure 6: effective G and c vs time ----
    plt.figure(figsize=(5, 4))
    plt.plot(t, G_eff_hist, label="G_eff(t)")
    plt.plot(t, c_eff_hist, label="c_eff(t)")
    plt.xlabel("t (ticks)")
    plt.ylabel("Effective constants")
    plt.title("G_eff and c_eff from flux cosmology")
    plt.legend()
    plt.grid(True)

    # ---- Figure 7: final vertical φ profile along axis ----
    plt.figure(figsize=(5, 4))
    plt.plot(z, phi_axis)
    plt.xlabel("z (height)")
    plt.ylabel("φ(z, axis)")
    plt.title("Final vertical flux fraction profile")
    plt.grid(True)

    plt.tight_layout()
    plt.show()


# ===============================================================
#  MAIN
# ===============================================================

if __name__ == "__main__":
    np.random.seed(1)  # for repeatable "turbulence" patterns

    result = run_full_jet_sim()
    make_plots(result)
