import numpy as np

# Rebuild the same grid and planetary IR layer used before
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
Phi_total_2d = Phi_solar_2d + Phi_cmb

Phi_plot = Phi_total_2d.copy()
Phi_plot[R_AU < 0.1] = Phi_plot[R_AU >= 0.1].max()

# Planet orbits and approximate equilibrium temperatures (K)
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

Phi_planetary_2d = np.zeros_like(Phi_total_2d)
sigma_AU = 0.5

for name, (r_au, T) in planet_orbits.items():
    strength = (T**4)
    dx = X_AU - r_au
    dy = Y_AU - 0.0
    R_local2 = dx**2 + dy**2
    gaussian = np.exp(-R_local2 / (2 * sigma_AU**2))
    Phi_planetary_2d += strength * gaussian

# Normalise as before
idx_1au_x = np.argmin(np.abs(x_AU - 1.0))
idx_1au_y = np.argmin(np.abs(y_AU - 0.0))
Phi_solar_1au = Phi_solar_2d[idx_1au_y, idx_1au_x]

scale_factor = Phi_solar_1au * 1e-2 / Phi_planetary_2d.max()
Phi_planetary_2d_scaled = Phi_planetary_2d * scale_factor

# Now sample Φ_planetary at Venus, Earth, Jupiter
targets = ["Venus", "Earth", "Jupiter"]
Phi_vals = {}

for name in targets:
    r_au, T = planet_orbits[name]
    ix = np.argmin(np.abs(x_AU - r_au))
    iy = np.argmin(np.abs(y_AU - 0.0))
    Phi_vals[name] = Phi_planetary_2d_scaled[iy, ix]

Phi_vals
