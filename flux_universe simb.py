import math
import matplotlib.pyplot as plt

# -----------------------------
# Parameters (tweak these)
# -----------------------------

# Time settings
T_MAX = 50.0       # total dimensionless time
DT    = 0.01       # time step (global tick size)

# Holographic bound scaling
s0 = 1.0           # S_max(a) = s0 * a^2

# Entropy production rate
gamma = 0.2        # how fast flux hardens (dS/dt = gamma * (S_max - S))

# FRW-like density parameters (dimensionless)
Omega_m = 0.3      # "matter"
Omega_r = 0.0      # "radiation" (set small or zero for simplicity)
k_flux   = 0.7     # weight of flux term in H^2

# Reference scales for G and c (arbitrary units)
G0 = 1.0
c0 = 1.0
H0 = 1.0           # reference H for c_eff coupling

# Coupling strengths to flux / expansion
#   g1 > 0 : more flux -> stronger gravity
#   g1 < 0 : more flux -> weaker gravity (gravity grows as flux hardens)
g1 = -0.5          # try -0.5 or +0.5 and compare
c1 = 0.3           # controls how strongly H alters c_eff

# Initial cosmological conditions
a0 = 0.1           # small universe
S0 = 0.0           # zero initial entropy

# -----------------------------
# Earth–Moon toy model params
# -----------------------------

M_earth = 1.0      # central mass (dimensionless)
R0 = 1.0           # initial physical orbital radius

# For a circular orbit in *static* space with G0:
V0 = math.sqrt(G0 * M_earth / R0)


# -----------------------------
# Cosmology stepper
# -----------------------------

def step_cosmology(a, S, dt):
    """
    One integration step for the cosmology using simple Euler scheme.
    Returns updated (a_new, S_new, H, flux_fraction, G_eff, c_eff).
    """
    # Prevent degenerate scale factor
    if a <= 1e-8:
        a = 1e-8

    # Holographic bound and flux
    S_max = s0 * a * a
    Phi = max(S_max - S, 0.0)
    flux_fraction = 0.0 if S_max == 0 else Phi / S_max  # between 0 and 1

    # Entropy evolution
    dS_dt = gamma * Phi
    S_new = S + dS_dt * dt

    # Expansion rate H from toy Friedmann
    a_safe = max(a, 1e-4)
    H_sq = Omega_m / (a_safe ** 3) + Omega_r / (a_safe ** 4) + k_flux * flux_fraction
    H_sq = max(H_sq, 0.0)
    H = math.sqrt(H_sq)

    # Scale factor evolution
    da_dt = H * a
    a_new = a + da_dt * dt

    # Effective gravity and speed of light
    G_eff = G0 * (1.0 + g1 * flux_fraction)
    if H0 != 0:
        c_eff = c0 * (1.0 + c1 * (H / H0 - 1.0))
    else:
        c_eff = c0

    return a_new, S_new, H, flux_fraction, G_eff, c_eff


def run_sim():
    # Time grid
    n_steps = int(T_MAX / DT)
    times = []

    # Cosmology histories
    a_vals = []
    S_vals = []
    flux_vals = []
    H_vals = []
    G_vals = []
    c_vals = []
    addota_vals = []   # (ddot a / a) history

    # Earth–Moon histories
    R_vals = []     # orbital radius vs time
    x_vals = []     # x position (for orbit plot)
    y_vals = []     # y position

    # Initial cosmological state
    a = a0
    S = S0

    # Initial orbital state (Moon at (R0, 0), velocity (0, V0))
    x = R0
    y = 0.0
    vx = 0.0
    vy = V0

    prev_H = None

    for i in range(n_steps):
        t = i * DT
        times.append(t)

        # ---- cosmology step ----
        a_vals.append(a)
        S_vals.append(S)

        a, S, H, flux_frac, G_eff, c_eff = step_cosmology(a, S, DT)

        flux_vals.append(flux_frac)
        H_vals.append(H)
        G_vals.append(G_eff)
        c_vals.append(c_eff)

        # Compute (ddot a / a) ≈ dH/dt + H^2   (McVittie / FRW relation)
        if prev_H is None:
            dH_dt = 0.0
        else:
            dH_dt = (H - prev_H) / DT
        prev_H = H

        a_ddot_over_a = 1e-20

        addota_vals.append(a_ddot_over_a)

        # ---- Earth–Moon step (McVittie-style) ----
        # Physical position relative to Earth (at origin)
        R = math.sqrt(x * x + y * y)
        R_safe = max(R, 1e-6)   # avoid divide-by-zero

        # Gravitational + expansion acceleration:
        # a_vec = [ -G_eff * M / R^3  + (ddot a/a) ] * R_vec
        grav_factor = - G_eff * M_earth / (R_safe ** 3)
        exp_factor  = a_ddot_over_a    # multiplies R_vec linearly

        ax = grav_factor * x + exp_factor * x
        ay = grav_factor * y + exp_factor * y

        # Integrate orbit
        vx += ax * DT
        vy += ay * DT
        x += vx * DT
        y += vy * DT

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
        "c_eff": c_vals,
        "addota": addota_vals,
        "R": R_vals,
        "x": x_vals,
        "y": y_vals,
    }


if __name__ == "__main__":
    data = run_sim()

    t = data["t"]

    # 1. Scale factor
    plt.figure()
    plt.plot(t, data["a"])
    plt.xlabel("t (ticks)")
    plt.ylabel("a(t)")
    plt.title("Scale factor a(t)")

    # 2. Entropy and flux fraction
    plt.figure()
    plt.plot(t, data["S"], label="S(t)")
    plt.plot(t, data["flux_frac"], label="Flux fraction Φ/S_max")
    plt.xlabel("t (ticks)")
    plt.ylabel("Entropy / Flux")
    plt.title("Entropy and flux fraction")
    plt.legend()

    # 3. Effective gravity
    plt.figure()
    plt.plot(t, data["G_eff"])
    plt.xlabel("t (ticks)")
    plt.ylabel("G_eff(t)")
    plt.title("Effective gravity vs time")

    # 4. Effective c
    plt.figure()
    plt.plot(t, data["c_eff"])
    plt.xlabel("t (ticks)")
    plt.ylabel("c_eff(t)")
    plt.title("Effective c vs time")

    # 5. (ddot a / a)
    plt.figure()
    plt.plot(t[:len(data["addota"])], data["addota"])
    plt.xlabel("t (ticks)")
    plt.ylabel("ddot a / a")
    plt.title("Cosmic acceleration (ddot a / a)")

    # 6. Orbital radius vs time
    plt.figure()
    plt.plot(t, data["R"])
    plt.xlabel("t (ticks)")
    plt.ylabel("R(t)")
    plt.title("Orbital radius (Moon) vs time (McVittie toy)")

    # 7. Orbital trajectory
    plt.figure()
    plt.plot(data["x"], data["y"])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.title("Orbital path of Moon around Earth (toy McVittie orbit)")

    plt.show()
