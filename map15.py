import numpy as np
import matplotlib.pyplot as plt

# Base grid (reuse idea)
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

# Mappings
rho_lin = Phi
rho_sq = Phi**2
rho_cub = Phi**3

# Thresholded: only regions with Φ > 0.5 contribute, and linearly
threshold = 0.5
rho_thr = np.where(Phi > threshold, Phi, 0.0)

V_lin = solve_potential(rho_lin)
V_sq = solve_potential(rho_sq)
V_cub = solve_potential(rho_cub)
V_thr = solve_potential(rho_thr)

# Gravity along x-axis (y=0)
ix0 = np.argmin(np.abs(y - 0.0))

def gravity_slice(V):
    dVdx, dVdy = np.gradient(V, x[1]-x[0], y[1]-y[0])
    gX = -dVdx
    gY = -dVdy
    g_mag = np.sqrt(gX**2 + gY**2)
    return g_mag[ix0, :]

g_lin = gravity_slice(V_lin)
g_sq = gravity_slice(V_sq)
g_cub = gravity_slice(V_cub)
g_thr = gravity_slice(V_thr)

# Plot potentials
plt.figure(figsize=(12,6))

plt.subplot(2,3,1)
plt.imshow(V_lin, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
plt.title("V (ρ ∝ Φ)")
plt.colorbar()

plt.subplot(2,3,2)
plt.imshow(V_cub, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
plt.title("V (ρ ∝ Φ³)")
plt.colorbar()

plt.subplot(2,3,3)
plt.imshow(V_thr, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
plt.title(f"V (thresholded Φ>{threshold})")
plt.colorbar()

# Gravity slices
plt.subplot(2,1,2)
plt.plot(x, g_lin, label="ρ ∝ Φ")
plt.plot(x, g_sq, label="ρ ∝ Φ²", linestyle='--')
plt.plot(x, g_cub, label="ρ ∝ Φ³", linestyle='-.')
plt.plot(x, g_thr, label=f"threshold Φ>{threshold}", linestyle=':')
plt.xlabel("x (y=0 slice)")
plt.ylabel("|g|")
plt.title("Gravity strength along x-axis for different ρ(Φ) mappings")
plt.legend()

plt.tight_layout()
plt.show()
