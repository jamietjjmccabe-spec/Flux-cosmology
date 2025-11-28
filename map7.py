import numpy as np

# Reuse solar + planetary setup from previous runs

L_sun = 3.828e26        # W
c = 3e8                 # m/s
eV_to_J = 1.602176634e-19

E_eff_eV = 2.0
E_eff_J = E_eff_eV * eV_to_J

E0_eV = 1.0
p = 1.0
phi_max = 1.0

def phi_E(E_eV, E0_eV=1.0, p=1.0, phi_max=1.0):
    E = E_eV
    return phi_max * (E**p) / (E**p + E0_eV**p)

phi_solar = phi_E(E_eff_eV, E0_eV, p, phi_max)

n_cmb = 4.11e8  # m^-3
E_cmb_eV = 6e-4
phi_cmb = phi_E(E_cmb_eV, E0_eV, p, phi_max)
Phi_cmb = n_cmb * phi_cmb

AU = 1.496e11
max_r_AU = 50
n_grid = 500

x_AU = np.linspace(-max_r_AU, max_r_AU, n_grid)
y_AU = np.linspace(-max_r_AU, max_r_AU, n_grid)
X_AU, Y_AU = np.meshgrid(x_AU, y_AU)

R_AU = np.sqrt(X_AU**2 + Y_AU**2)
R_m = R_AU * AU
R_m_safe = np.where(R_m == 0, 1e7, R_m)

n_solar_2d = L_sun / (4 * np.pi * R_m_safe**2 * c * E_eff_J)
Phi_solar_2d = n_solar_2d * phi_solar

# Planet orbits + temps
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

# Rebuild planetary IR field
Phi_planetary_2d = np.zeros_like(Phi_solar_2d)
sigma_AU = 0.5

for name, (r_au, T) in planet_orbits.items():
    strength = (T**4)
    dx = X_AU - r_au
    dy = Y_AU - 0.0
    R_local2 = dx**2 + dy**2
    gaussian = np.exp(-R_local2 / (2 * sigma_AU**2))
    Phi_planetary_2d += strength * gaussian

# Normalise planetary IR as before
idx_1au_x = np.argmin(np.abs(x_AU - 1.0))
idx_1au_y = np.argmin(np.abs(y_AU - 0.0))
Phi_solar_1au = Phi_solar_2d[idx_1au_y, idx_1au_x]

scale_factor = Phi_solar_1au * 1e-2 / Phi_planetary_2d.max()
Phi_planetary_2d_scaled = Phi_planetary_2d * scale_factor

# Define some moons: Io, Europa (around Jupiter), Titan (around Saturn)
# We'll give them approximate effective temperatures (very rough)
moons = {
    "Io":      {"host": "Jupiter", "d_host_AU": 0.0028, "T": 110},  # Io's orbital radius ~ 421,700 km ≈ 0.0028 AU
    "Europa":  {"host": "Jupiter", "d_host_AU": 0.0045, "T": 100},  # ~ 671,000 km ≈ 0.0045 AU
    "Titan":   {"host": "Saturn",  "d_host_AU": 0.0082, "T": 95},   # ~ 1,221,000 km ≈ 0.0082 AU
}

# We'll compute each moon's IR Φ at its location and compare:
# - to solar Φ at that radius
# - to host planet's IR Φ at same location

results = {}

for moon, props in moons.items():
    host = props["host"]
    d_host = props["d_host_AU"]
    T_moon = props["T"]
    r_host_au, T_host = planet_orbits[host]
    
    # Place host on +x axis; put moon slightly offset in y by d_host for simplicity
    x_moon = r_host_au
    y_moon = d_host  # small vertical offset
    
    ix = np.argmin(np.abs(x_AU - x_moon))
    iy = np.argmin(np.abs(y_AU - y_moon))
    
    # Solar Φ at moon location
    Phi_sun_moon = Phi_solar_2d[iy, ix]
    
    # Host planetary IR Φ at moon location (from the already-built planetary field)
    Phi_host_IR_at_moon = Phi_planetary_2d_scaled[iy, ix]
    
    # Now add a small Gaussian puff for the moon itself and sample its peak
    sigma_moon_AU = 0.0005  # tighter puff than planets
    dxm = X_AU - x_moon
    dym = Y_AU - y_moon
    Rm2 = dxm**2 + dym**2
    moon_strength = (T_moon**4)
    moon_gauss = np.exp(-Rm2 / (2 * sigma_moon_AU**2))
    Phi_moon_field = moon_strength * moon_gauss
    
    # Scale moons using the same overall scale factor (for consistency)
    Phi_moon_field_scaled = Phi_moon_field * scale_factor
    
    Phi_moon_peak = Phi_moon_field_scaled[iy, ix]
    
    results[moon] = {
        "Phi_moon": Phi_moon_peak,
        "Phi_sun": Phi_sun_moon,
        "Phi_host_IR": Phi_host_IR_at_moon,
        "Phi_moon_over_sun": Phi_moon_peak / Phi_sun_moon,
        "Phi_moon_over_hostIR": Phi_moon_peak / Phi_host_IR_at_moon if Phi_host_IR_at_moon > 0 else np.inf
    }

results
