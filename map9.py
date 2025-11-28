import numpy as np

# We'll recompute moon Φ more directly: the puff peak at the moon center is just strength * scale_factor
# Using previous planet_orbits, moons, scale_factor, Phi_solar_2d, Phi_planetary_2d_scaled, x_AU, y_AU

moons = {
    "Io":      {"host": "Jupiter", "d_host_AU": 0.0028, "T": 110},
    "Europa":  {"host": "Jupiter", "d_host_AU": 0.0045, "T": 100},
    "Titan":   {"host": "Saturn",  "d_host_AU": 0.0082, "T": 95},
}

results = {}

for moon, props in moons.items():
    host = props["host"]
    d_host = props["d_host_AU"]
    T_moon = props["T"]
    r_host_au, T_host = planet_orbits[host]
    
    x_moon = r_host_au
    y_moon = d_host
    
    ix = np.argmin(np.abs(x_AU - x_moon))
    iy = np.argmin(np.abs(y_AU - y_moon))
    
    Phi_sun_moon = Phi_solar_2d[iy, ix]
    Phi_host_IR_at_moon = Phi_planetary_2d_scaled[iy, ix]
    
    # Peak Φ for moon IR puff at its own center
    moon_strength = (T_moon**4)
    Phi_moon_peak = moon_strength * scale_factor
    
    results[moon] = {
        "Phi_moon": Phi_moon_peak,
        "Phi_sun": Phi_sun_moon,
        "Phi_host_IR": Phi_host_IR_at_moon,
        "Phi_moon_over_sun": Phi_moon_peak / Phi_sun_moon,
        "Phi_moon_over_hostIR": Phi_moon_peak / Phi_host_IR_at_moon if Phi_host_IR_at_moon > 0 else np.inf
    }

results
