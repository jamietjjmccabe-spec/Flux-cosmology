import numpy as np

# We'll build a local Earth–Moon IR Φ model in 2D around Earth
AU = 1.496e11  # m
km_to_m = 1000.0

# Earth–Moon distance and a small shift
d0_km = 384400.0     # current mean distance
d0_AU = (d0_km * km_to_m) / AU

# We'll use a larger delta (1000 km) to compute a derivative, then scale down to 1 cm
dd_km = 1000.0
dd_AU = (dd_km * km_to_m) / AU

# Temperatures (rough averages)
T_earth = 288.0  # K
T_moon = 220.0   # K

# Relative IR "strength" ~ T^4 (no absolute scaling needed)
S_earth = T_earth**4
S_moon = T_moon**4

# Gaussian puff widths in AU (chosen so blobs are local but overlapping)
sigma_earth_AU = 0.001   # ~150,000 km
sigma_moon_AU = 0.001

# 2D grid around Earth: +/- 0.01 AU (~1.5 million km)
extent_AU = 0.01
n_grid = 600
x_AU = np.linspace(-extent_AU, extent_AU, n_grid)
y_AU = np.linspace(-extent_AU, extent_AU, n_grid)
X_AU, Y_AU = np.meshgrid(x_AU, y_AU)

def phi_field(d_AU):
    """
    Build Φ-field from Earth at (0,0) and Moon at (d_AU,0)
    using Gaussian IR puffs with strengths S_earth, S_moon.
    """
    # Earth Gaussian
    R2_earth = X_AU**2 + Y_AU**2
    phi_earth = S_earth * np.exp(-R2_earth / (2 * sigma_earth_AU**2))
    # Moon Gaussian
    R2_moon = (X_AU - d_AU)**2 + Y_AU**2
    phi_moon = S_moon * np.exp(-R2_moon / (2 * sigma_moon_AU**2))
    return phi_earth + phi_moon

# Compute Φ-fields at d0 and d0 + dd
Phi0 = phi_field(d0_AU)
Phi1 = phi_field(d0_AU + dd_AU)

# Define a threshold as a fraction of the peak Φ at current configuration
Phi0_max = Phi0.max()
threshold = 0.1 * Phi0_max  # 10% of max as "environment" contour

# Compute "environment area" where Φ >= threshold
dx = (2 * extent_AU) / (n_grid - 1)
dy = dx
cell_area_AU2 = dx * dy

A0 = np.sum(Phi0 >= threshold) * cell_area_AU2
A1 = np.sum(Phi1 >= threshold) * cell_area_AU2

dA_dAU = (A1 - A0) / dd_AU       # derivative of area wrt separation, in AU^2 / AU = AU
# Now scale that to area change per 1 cm of lunar recession
cm_to_AU = (0.01 * km_to_m) / AU  # 1 cm in AU
dA_per_cm_AU2 = dA_dAU * cm_to_AU

A0_km2 = A0 * (AU/1000.0)**2
dA_per_cm_km2 = dA_per_cm_AU2 * (AU/1000.0)**2

(A0, A0_km2, dA_per_cm_AU2, dA_per_cm_km2)
