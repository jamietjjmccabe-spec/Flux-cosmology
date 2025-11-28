import numpy as np
import matplotlib.pyplot as plt

# Toy model: gravity emerging from Φ-density in 2D

# 1. Define a 2D grid
L = 10.0   # box size in arbitrary units
N = 400
x = np.linspace(-L/2, L/2, N)
y = np.linspace(-L/2, L/2, N)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2)

# 2. Define a Φ-density field: a central Gaussian "flux well"
Phi0 = 1.0
sigma = 1.0
Phi = Phi0 * np.exp(-R**2 / (2 * sigma**2))

# 3. Define an effective mass density from Φ
# Here we simply take rho_eff ∝ Φ (toy identification)
rho_eff = Phi

# 4. Solve a discrete Poisson equation ∇²V = 4πG * rho_eff (G=1 in toy units)
# We'll do it in Fourier space for speed: V(k) = -4πG * rho(k) / k^2
G_eff = 1.0

# FFT of rho
rho_k = np.fft.fft2(rho_eff)

# k-space grid
kx = np.fft.fftfreq(N, d=(x[1]-x[0])) * 2*np.pi
ky = np.fft.fftfreq(N, d=(y[1]-y[0])) * 2*np.pi
KX, KY = np.meshgrid(kx, ky)
K2 = KX**2 + KY**2

# Avoid division by zero at k=0 (set mean potential to zero)
K2[0, 0] = 1.0

V_k = -4 * np.pi * G_eff * rho_k / K2

# Set k=0 mode to zero (no arbitrary constant offset)
V_k[0, 0] = 0.0

# Inverse FFT to get potential in real space
V = np.fft.ifft2(V_k).real

# 5. Compute "gravitational field" g = -∇V
dVdx, dVdy = np.gradient(V, x[1]-x[0], y[1]-y[0])
gX = -dVdx
gY = -dVdy
g_mag = np.sqrt(gX**2 + gY**2)

# 6. Plot Φ, V, and |g|
plt.figure(figsize=(14,4))

plt.subplot(1,3,1)
plt.imshow(Phi, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
plt.colorbar(label="Φ-density")
plt.title("Φ-density (flux well)")
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(1,3,2)
plt.imshow(V, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
plt.colorbar(label="V (toy gravitational potential)")
plt.title("Emergent potential from Φ")
plt.xlabel("x")

plt.subplot(1,3,3)
plt.imshow(g_mag, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
plt.colorbar(label="|g| (toy gravity strength)")
plt.title("Gravity magnitude from ∇V")
plt.xlabel("x")

plt.tight_layout()
plt.show()
