import numpy as np
import matplotlib.pyplot as plt

# Physical constants
L_sun = 3.828e26        # W
c = 3e8                 # m/s
eV_to_J = 1.602176634e-19

# Photon energy assumptions
E_eff_eV = 2.0          # effective solar photon energy (visible-ish)
E_eff_J = E_eff_eV * eV_to_J

# Phi model params
E0_eV = 1.0
p = 1.0
phi_max = 1.0

def phi_E(E_eV, E0_eV=1.0, p=1.0, phi_max=1.0):
    E = E_eV
    return phi_max * (E**p) / (E**p + E0_eV**p)

phi_solar = phi_E(E_eff_eV, E0_eV, p, phi_max)

# CMB floor
T_cmb = 2.725  # K (not used exactly, but for context)
n_cmb = 4.11e8  # m^-3
E_cmb_eV = 6e-4
phi_cmb = phi_E(E_cmb_eV, E0_eV, p, phi_max)
Phi_cmb = n_cmb * phi_cmb

# Radial grid: 0.1 to 50 AU
AU = 1.496e11
r_AU = np.linspace(0.1, 50, 500)
r_m = r_AU * AU

# Solar photon number density at radius r:
# n(r) = L / (4π r^2 c E_ph)
n_solar = L_sun / (4 * np.pi * r_m**2 * c * E_eff_J)

Phi_solar = n_solar * phi_solar
Phi_total = Phi_solar + Phi_cmb

# Key planet positions
planet_positions = {
    "Mercury": 0.39,
    "Venus": 0.72,
    "Earth": 1.0,
    "Mars": 1.52,
    "Jupiter": 5.2,
    "Saturn": 9.5,
    "Uranus": 19.2,
    "Neptune": 30.1,
    "Pluto": 39.5
}

planet_data = {}
for name, r_au in planet_positions.items():
    r = r_au * AU
    n = L_sun / (4 * np.pi * r**2 * c * E_eff_J)
    Phi_s = n * phi_solar
    Phi_t = Phi_s + Phi_cmb
    planet_data[name] = {
        "r_AU": r_au,
        "n_solar": n,
        "Phi_solar": Phi_s,
        "Phi_total": Phi_t
    }

planet_data

# Make a plot of Phi_total vs radius in AU on log scale

plt.figure()
plt.plot(r_AU, Phi_total)
plt.yscale("log")
plt.xlabel("Radius r [AU]")
plt.ylabel("Flux potential density Φ_tot [arb units]")
plt.title("Toy Flux-Potential Profile of the Solar System (Sun + CMB floor)")

# Mark planet positions
for name, r_au in planet_positions.items():
    r_index = np.argmin(np.abs(r_AU - r_au))
    plt.scatter(r_au, Phi_total[r_index])
    plt.text(r_au, Phi_total[r_index]*1.1, name, fontsize=8, ha='center')

plt.tight_layout()
plt.show()
