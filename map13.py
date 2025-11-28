import numpy as np
import matplotlib.pyplot as plt

# Base grid
L = 10.0
N = 400
x = np.linspace(-L/2, L/2, N)
y = np.linspace(-L/2, L/2, N)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2)

# Base Φ well
Phi0 = 1.0
sigma = 1.0
Phi = Phi0 * np.exp(-R**2 / (2 * sigma**2))

# Helper to solve Poisson for a given rho
def solve_potential(rho):
    rho_k = np.fft.fft2(rho)
    kx = np.fft.fftfreq(N, d=(x[1]-x[0])) * 2*np.pi
    ky = np.fft.fftfreq(N, d=(y[1]-y[0])) * 2*np.pi
    KX, KY = np.meshgrid(kx, ky)
    K2 = KX**2 + KY**2
    K2[0, 0] = 1.0
    V_k = -4*np.pi * rho_k / K2
    V_k[0, 0] = 0.0
    V = np.fft.ifft2(V_k).real
    return V

# Linear: rho ∝ Φ
rho_lin = Phi
V_lin = solve_potential(rho_lin)

# Nonlinear: rho ∝ Φ^2
rho_sq = Phi**2
V_sq = solve_potential(rho_sq)

# Also compute gravity magnitude along x-axis (y=0) for both
ix0 = np.argmin(np.abs(y - 0.0))

dVdx_lin, dVdy_lin = np.gradient(V_lin, x[1]-x[0], y[1]-y[0])
gX_lin = -dVdx_lin
g_mag_lin = np.sqrt(gX_lin**2 + dVdy_lin**2)
g_line_lin = g_mag_lin[ix0, :]

dVdx_sq, dVdy_sq = np.gradient(V_sq, x[1]-x[0], y[1]-y[0])
gX_sq = -dVdx_sq
g_mag_sq = np.sqrt(gX_sq**2 + dVdy_sq**2)
g_line_sq = g_mag_sq[ix0, :]

# Plot potentials and 1D gravity cuts
plt.figure(figsize=(12,5))

plt.subplot(1,3,1)
plt.imshow(V_lin, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
plt.colorbar(label="V_lin")
plt.title("Linear mapping: ρ ∝ Φ")

plt.subplot(1,3,2)
plt.imshow(V_sq, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
plt.colorbar(label="V_sq")
plt.title("Nonlinear mapping: ρ ∝ Φ²")

plt.subplot(1,3,3)
plt.plot(x, g_line_lin, label="ρ ∝ Φ")
plt.plot(x, g_line_sq, label="ρ ∝ Φ²", linestyle='--')
plt.xlabel("x (y=0 slice)")
plt.ylabel("|g|")
plt.title("Gravity strength along x-axis")
plt.legend()

plt.tight_layout()
plt.show()
