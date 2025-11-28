import numpy as np

# Reuse the same structures from the last execution (assuming state persists)
planet_orbits = {
    "Mercury": (0.39, 440),
    "Venus": (0.72, 730),
    "Earth": (1.0, 288),
    "Mars": (1.52, 210),
    "Jupiter": (5.2, 125),
    "Saturn": (9.5, 95),
    "Uranus": (19.2, 60),
    "Neptune": (30.1, 60),
    "Pluto": (39.5, 44)
}

# We already have x_AU, y_AU, Phi_planetary_2d_scaled from previous cell
Phi_all = {}
for name, (r_au, T) in planet_orbits.items():
    ix = np.argmin(np.abs(x_AU - r_au))
    iy = np.argmin(np.abs(y_AU - 0.0))
    Phi_all[name] = Phi_planetary_2d_scaled[iy, ix]

# Rank by Φ descending
ranked = sorted(Phi_all.items(), key=lambda kv: kv[1], reverse=True)
ranked

