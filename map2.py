
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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
n_cmb = 4.11e8  # m^-3
E_cmb_eV = 6e-4
phi_cmb = phi_E(E_cmb_eV, E0_eV, p, phi_max)
Phi_cmb = n_cmb * phi_cmb

AU = 1.496e11

# 2D grid in x,y (in AU), centered on the Sun
max_r_AU = 50
n_grid = 500

x_AU = np.linspace(-max_r_AU, max_r_AU, n_grid)
y_AU = np.linspace(-max_r_AU, max_r_AU, n_grid)
X_AU, Y_AU = np.meshgrid(x_AU, y_AU)

R_AU = np.sqrt(X_AU**2 + Y_AU**2)
R_m = R_AU * AU

# Avoid division by zero at the origin by setting a minimum radius
R_m_safe = np.where(R_m == 0, 1e7, R_m)  # 10,000 km as arbitrary tiny radius

# Solar photon number density and flux potential
n_solar_2d = L_sun / (4 * np.pi * R_m_safe**2 * c * E_eff_J)
Phi_solar_2d = n_solar_2d * phi_solar

Phi_total_2d = Phi_solar_2d + Phi_cmb

# Mask the very central pixel area for plotting clarity (optional)
Phi_plot = Phi_total_2d.copy()
Phi_plot[R_AU < 0.1] = Phi_plot[R_AU >= 0.1].max()

# Plot heatmap
plt.figure(figsize=(8, 7))
im = plt.imshow(
    Phi_plot,
    origin="lower",
    extent=[-max_r_AU, max_r_AU, -max_r_AU, max_r_AU],
    norm=LogNorm(),
)
plt.colorbar(im, label="Flux potential density Φ_tot [arb units]")
plt.scatter(0, 0, color="white", s=40, marker="*")
plt.text(0, 0.5, "Sun", color="white", ha="center", va="bottom", fontsize=9)

plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
plt.title("2D Flux-Web Heatmap of the Solar System (Sun + CMB floor)")
plt.tight_layout()
plt.show()
