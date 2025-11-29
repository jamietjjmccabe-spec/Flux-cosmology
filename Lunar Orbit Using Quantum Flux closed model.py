import math
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Cosmology parameters
# -----------------------------

T_MAX = 50.0       # total time (years)
DT    = 0.01       # time step (years)

s0    = 1.0
gamma = 0.01
Omega_m = 0.3
Omega_r = 0.0
k_flux  = 0.05

G0 = 1.0
c0 = 1.0
H0 = 1.0

# Initial cosmological conditions
a0 = 0.1
S0 = 0.0

# -----------------------------
# Earth–Moon dimensionless parameters
# -----------------------------

M_earth = 1.0
M_moon  = 0.0123      # Moon / Earth mass ratio
R0      = 1.0         # 1 unit = current Moon distance
V0      = math.sqrt(G0 * M_earth / R0)


# ======================================================================
#  COSMOLOGY STEPPER
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
    H_sq = Omega_m / (a_safe**3) + Omega_r/(a_safe**4) + k_flux*flux_fraction
    H_sq = max(H_sq, 0.0)
    H = math.sqrt(H_sq)

    # scale factor evolution
    a_new = a + (H * a) * dt

    # variable gravity
    G_eff = G0 * (1.0 + g1 * flux_fraction)

    # effective speed of light
    if H0 != 0:
        c_eff = c0 * (1 + 0.3*(H/H0 - 1))
    else:
        c_eff = c0

    return a_new, S_new, H, flux_fraction, G_eff, c_eff


# ======================================================================
#  ORBITAL SIMULATOR
# ======================================================================

def run_sim(g1, tidal_torque_factor=1e-4):
    """Simulate Moon orbit for 50 years with flux gravity and tidal torques."""
    n_steps = int(T_MAX / DT)
    times = []

    a_vals, S_vals, flux_vals = [], [], []
    H_vals, G_vals = [], []
    R_vals, x_vals, y_vals = [], [], []

    # initial cosmology
    a = a0
    S = S0

    # initial orbit (use G0 here; G_eff will evolve later)
    x, y  = R0, 0.0
    vx = 0.0
    vy = V0

    for i in range(n_steps):
        t = i * DT
        times.append(t)

        # ---- cosmology update ----
        a_vals.append(a)
        S_vals.append(S)

        a, S, H, flux_frac, G_eff, c_eff = step_cosmology(a, S, DT, g1)

        flux_vals.append(flux_frac)
        H_vals.append(H)
        G_vals.append(G_eff)

        # ---- Newtonian gravity ----
        R = math.sqrt(x*x + y*y)
        R_safe = max(R, 1e-6)

        ax = -G_eff * M_earth * x / (R_safe**3)
        ay = -G_eff * M_earth * y / (R_safe**3)

        vx += ax * DT
        vy += ay * DT
        x  += vx * DT
        y  += vy * DT

        # ======================================================================
        #  TIDAL TORQUE MODEL (actual cause of Moon recession)
        # ======================================================================

        # Earth-Moon tidal angular momentum gain
        L_dot = tidal_torque_factor * (M_earth * M_moon / R_safe**2) * G_eff

        # current orbital angular momentum per unit mass
        h = math.sqrt(G_eff * M_earth * R_safe)

        # update with torque
        h += (L_dot * DT) / M_moon

        # new semi-major axis (circular orbit approx)
        R_new = h*h / (G_eff * M_earth)

        # scale orbit outward
        scale = R_new / R_safe
        x  *= scale
        y  *= scale
        vx *= math.sqrt(scale)
        vy *= math.sqrt(scale)

        # ======================================================================

        R = math.sqrt(x*x + y*y)
        R_vals.append(R)
        x_vals.append(x)
        y_vals.append(y)

    return {
        "t": times,
        "a": a_vals,
        "S": S_vals,
        "flux_frac": flux_vals,
        "H": H_vals,
        "G_eff": G_vals,
        "R": R_vals,
        "x": x_vals,
        "y": y_vals,
    }


# ======================================================================
#  MEASURE FRACTIONAL RECESSION RATE
# ======================================================================

def fractional_recession_rate(data):
    """Return fractional rate (1/R)(dR/dt) per year."""
    t = data["t"]
    R = data["R"]
    return (R[-1] - R[0]) / (R[0] * (t[-1] - t[0]))


# ======================================================================
#  MAIN: AUTO-CALIBRATE TIDAL TORQUE TO 3.8 cm/year
# ======================================================================

if __name__ == "__main__":
    target_rate = 9.8855e-11  # real Moon recession 3.8 cm/year
    print("Target fractional recession:", target_rate)

    # Step 1 — Run once with strong torque (1e-4) to measure rate
    g1_fixed = -1e-6
    baseline_torque = 1e-4

    base_data = run_sim(g1_fixed, tidal_torque_factor=baseline_torque)
    base_rate = fractional_recession_rate(base_data)

    print(f"Baseline torque {baseline_torque:.3e} → rate {base_rate:.3e}")

    # Step 2 — Calculate correction factor
    scale = target_rate / base_rate
    torque_calib = baseline_torque * scale

    print(f"Calibrated tidal_torque_factor = {torque_calib:.3e}")

    # Step 3 — Re-run with calibrated torque
    data = run_sim(g1_fixed, tidal_torque_factor=torque_calib)
    rate = fractional_recession_rate(data)

    print(f"Check run rate = {rate:.3e}  (target {target_rate:.3e})")

    # ==================================================================
    #  PLOTS
    # ==================================================================

    t = data["t"]

    # raw radius
    plt.figure()
    plt.plot(t, data["R"])
    plt.xlabel("Time (years)")
    plt.ylabel("R / R0")
    plt.title(f"Orbital radius (torque={torque_calib:.2e})")

    # orbital path
    plt.figure()
    plt.plot(data["x"], data["y"])
    plt.axis("equal")
    plt.title("Orbital path")

    # G_eff
    plt.figure()
    plt.plot(t, data["G_eff"])
    plt.xlabel("Time (years)")
    plt.ylabel("G_eff(t)")
    plt.title(f"G_eff vs time (g1={g1_fixed})")

    # ---------- smoothed radius plot ----------
    R_arr = np.array(data["R"])
    window = 500            # ~5-year smoothing for DT=0.01

    def moving_average(x, w):
        x = np.array(x)
        out = np.zeros_like(x)
        for i in range(len(x)):
            start = max(0, i-w)
            out[i] = np.mean(x[start:i+1])
        return out

    smooth_R = moving_average(R_arr, window)

    plt.figure()
    plt.plot(t, smooth_R, label="Smoothed R(t)")
    plt.plot(t, R_arr, alpha=0.3, label="Raw R(t)")
    plt.xlabel("Time (years)")
    plt.ylabel("R / R0")
    plt.title("Smoothed Orbital Radius (Tidal Signal Visible)")
    plt.legend()

    plt.show()
