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
g1 = -0.5
.5           # controls how strongly flux alters G_eff
c1 = 0.3           # controls how strongly H alters c_eff

# Initial conditions
a0 = 0.1           # small universe
S0 = 0.0           # zero initial entropy


# -----------------------------
# Integration setup
# -----------------------------

def step(a, S, dt):
    """
    One integration step using simple Euler scheme.
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
    # Avoid division blow-ups for very small a
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
    a_vals = []
    S_vals = []
    flux_vals = []
    H_vals = []
    G_vals = []
    c_vals = []

    a = a0
    S = S0

    for i in range(n_steps):
        t = i * DT
        times.append(t)
        a_vals.append(a)
        S_vals.append(S)

        a, S, H, flux_frac, G_eff, c_eff = step(a, S, DT)

        flux_vals.append(flux_frac)
        H_vals.append(H)
        G_vals.append(G_eff)
        c_vals.append(c_eff)

    return {
        "t": times,
        "a": a_vals,
        "S": S_vals,
        "flux_frac": flux_vals,
        "H": H_vals,
        "G_eff": G_vals,
        "c_eff": c_vals,
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

    # 4. Effective speed of light
    plt.figure()
    plt.plot(t, data["c_eff"])
    plt.xlabel("t (ticks)")
    plt.ylabel("c_eff(t)")
    plt.title("Effective c vs time")

    plt.show()
