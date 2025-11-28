import numpy as np
import matplotlib.pyplot as plt

# Spatial grid
L = 10.0
N = 300
x = np.linspace(-L/2, L/2, N)
y = np.linspace(-L/2, L/2, N)
X, Y = np.meshgrid(x, y)
R2 = X**2 + Y**2

dx = x[1] - x[0]
dy = y[1] - y[0]
cell_area = dx * dy

Phi_c = 0.4  # threshold for gravitating mass

# Time evolution for sigma(t): smooth collapse from 3.0 to 1.0
T = 100
t_arr = np.linspace(0, 1, T)
sigma_start = 3.0
sigma_end = 1.0
sigma_t = sigma_start + (sigma_end - sigma_start) * t_arr  # linear in time; could choose nonlinear if desired

total_mass = []

for sigma in sigma_t:
    Phi = np.exp(-R2 / (2 * sigma**2))
    rho = np.where(Phi > Phi_c, Phi, 0.0)
    M = rho.sum() * cell_area
    total_mass.append(M)

total_mass = np.array(total_mass)

# Normalise mass to initial and final for easier comparison
M0 = total_mass[0]
Mf = total_mass[-1]

plt.figure(figsize=(7,4))
plt.plot(t_arr, total_mass, label="Total gravitating mass M(t)")
plt.axhline(M0, linestyle='--', label=f"M(0)={M0:.3f}")
plt.axhline(Mf, linestyle=':', label=f"M(1)={Mf:.3f}")
plt.xlabel("t (normalized time)")
plt.ylabel("M(t) (arbitrary units)")
plt.title("Total gravitating mass vs collapse (σ: 3.0 → 1.0)")
plt.legend()
plt.tight_layout()
plt.show()

M0, Mf
