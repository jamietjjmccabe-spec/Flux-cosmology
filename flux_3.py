import math
import matplotlib.pyplot as plt

# -----------------------------
# Cosmology parameters
# -----------------------------

T_MAX = 50.0       # total time (interpret as years)
DT    = 0.01       # step in "years" (so 100 steps per year)

s0    = 1.0        # S_max(a) = s0 * a^2
gamma = 0.01       # slower entropy production
Omega_m = 0.3
Omega_r = 0.0
k_flux  = 0.05     # weaker flux-driven expansion

G0 = 1.0
c0 = 1.0
H0 = 1.0

# Initial cosmological conditions
a0 = 0.1
S0 = 0.0

# -----------------------------
# Earth–Moon toy parameters (dimensionless)
# -----------------------------

M_earth = 1.0      # central mass
R0      = 1.0      # current mean Earth–Moon distance
V0      = math.sqrt(G0 * M_earth / R0)   # circular orbit in static space


# -----------------------------
# Cosmology stepper
# -----------------------------

def step_cosmology(a, S, dt, g1):
    """
    One integration step for the cosmology.
    Returns: a_new, S_new, H, flux_fraction, G_eff, c_eff
    """
    if a <= 1e-8:
        a = 1e-8

    S_max = s0 * a * a
    Phi   = max(S_max - S, 0.0)
    flux_fraction = 0.0 if S_max == 0 else Phi / S_max

    # entropy
    dS_dt = gamma * Phi
    S_new = S + dS_dt * dt

    # Hubble parameter
    a_safe = max(a, 1e-4)
    H_sq = Omega_m / (a_safe ** 3) + Omega_r / (a_safe ** 4) + k_flux * flux_fraction
    H_sq = max(H_sq, 0.0)
    H    = math.sqrt(H_sq)

    # scale factor evolution
    da_dt = H * a
    a_new = a + da_dt * dt

    # effective gravity and c
    G_eff = G0 * (1.0 + g1 * flux_fraction)
    if H0 != 0:
        c_eff = c0 * (1.0 + 0.3 * (H / H0 - 1.0))
    else:
        c_eff = c0

    return a_new, S_new, H, flux_fraction, G_eff, c_eff

# --- 1. TIDAL TORQUE MODEL (added outward push) ---

# adjustable factor to calibrate to 3.8 cm/year
tidal_torque_factor = 1e-4   # tune this later

# current radius
R_moon = math.sqrt(x*x + y*y)

# 1. Angular momentum gain from tidal torque
#    L_dot has dimensions of angular momentum per unit time
L_dot = tidal_torque_factor * (M_earth * M_moon / R_moon**2) * G_eff

# 2. Convert L_dot → specific angular momentum change (h)
#    h = |r × v| for circular-ish orbit ≈ sqrt(G M R)
h = math.sqrt(G_eff * M_earth * R_moon)

# update specific angular momentum
h += L_dot * DT / M_moon

# 3. Update semi-major axis based on new angular momentum
#    For a circular orbit: h^2 = G M a   →  a = h^2 / (G M)
R_moon_new = h*h / (G_eff * M_earth)

# scale the orbit outward by adjusting x,y proportionally
scale = R_moon_new / R_moon
x *= scale
y *= scale

# velocity must also be scaled to keep centrifugal balance
vx *= math.sqrt(scale)
vy *= math.sqrt(scale)


# -----------------------------
# Main simulator for a given g1
# -----------------------------

def run_sim(g1):
    n_steps = int(T_MAX / DT)
    times   = []

    a_vals      = []
    S_vals      = []
    flux_vals   = []
    H_vals      = []
    G_vals      = []

    R_vals = []
    x_vals = []
    y_vals = []

    # cosmology
    a = a0
    S = S0

    # Moon orbit (physical radius R, not comoving)
    x  = R0
    y  = 0.0
    vx = 0.0
    vy = V0

    for i in range(n_steps):
        t = i * DT
        times.append(t)

        # --- cosmology ---
        a_vals.append(a)
        S_vals.append(S)

        a, S, H, flux_frac, G_eff, c_eff = step_cosmology(a, S, DT, g1)

        flux_vals.append(flux_frac)
        H_vals.append(H)
        G_vals.append(G_eff)

        # --- Earth–Moon orbit with NO local expansion term ---
        R = math.sqrt(x * x + y * y)
        R_safe = max(R, 1e-6)

        # a_vec = -G_eff * M * r_vec / R^3
        grav_factor = - G_eff * M_earth / (R_safe ** 3)
        ax = grav_factor * x
        ay = grav_factor * y

        vx += ax * DT
        vy += ay * DT
        x  += vx * DT
        y  += vy * DT

        R_vals.append(R_safe)
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


# -----------------------------
# Calibration scan for g1
# -----------------------------

def fractional_recession_rate(data):
    """
    Returns <(1/R) dR/dt> over whole run (per time unit).
    """
    t = data["t"]
    R = data["R"]
    R0 = R[0]
    R_end = R[-1]
    T = t[-1] - t[0]
    return (R_end - R0) / (R0 * T)


if __name__ == "__main__":
    # target fractional rate: 3.8 cm / 384400 km per year
    target_rate = 9.8855e-11

    print("Target fractional recession rate per year:", target_rate)

    # coarse scan over g1 values (all negative: flux hardening -> stronger gravity)
    candidates = [-1e-3, -5e-4, -1e-4, -5e-5, -1e-5, -5e-6, -1e-6]

    best_g1 = None
    best_err = None

    for g in candidates:
        data = run_sim(g)
        rate = fractional_recession_rate(data)
        err  = abs(rate - target_rate)

        print(f"g1 = {g: .1e}  ->  frac rate = {rate: .3e}  (err = {err: .3e})")

        if best_err is None or err < best_err:
            best_err = err
            best_g1  = g

    print("\nBest coarse g1 ≈", best_g1, "with error", best_err)

    # rerun with best g1 and show plots
    data = run_sim(best_g1)
    t = data["t"]

    plt.figure()
    plt.plot(t, data["G_eff"])
    plt.xlabel("t (years)")
    plt.ylabel("G_eff(t)")
    plt.title(f"Effective gravity vs time (g1 = {best_g1})")

    plt.figure()
    plt.plot(t, data["R"])
    plt.xlabel("t (years)")
    plt.ylabel("R(t) [units of current Moon distance]")
    plt.title("Orbital radius vs time")

    plt.figure()
    plt.plot(data["x"], data["y"])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.title("Orbital path of Moon (dimensionless)")

    plt.show()
