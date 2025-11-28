import numpy as np
import matplotlib.pyplot as plt

# 2-cloud time-dependent Φ model

# Spatial grid
L = 10.0
N = 300
x = np.linspace(-L/2, L/2, N)
y = np.linspace(-L/2, L/2, N)
X, Y = np.meshgrid(x, y)

dx = x[1] - x[0]
dy = y[1] - y[0]
cell_area = dx * dy

# Cloud positions
x1, y1 = -2.0, 0.0
x2, y2 =  2.0, 0.0
sigma = 1.0
Phi_c = 0.4

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

def gravity_mag(V):
    dVdx, dVdy = np.gradient(V, x[1]-x[0], y[1]-y[0])
    gX = -dVdx
    gY = -dVdy
    return np.sqrt(gX**2 + gY**2)

# Time evolution of peak Φ0(t)
T = 80
t_arr = np.linspace(0, 1, T)
Phi0_start = 0.1
Phi0_end = 1.0
Phi0_t = Phi0_start + (Phi0_end - Phi0_start) * t_arr

# Choose 4 representative times to show snapshots
snapshot_indices = [0, int(T*0.3), int(T*0.6), T-1]

snapshots = []

for idx, Phi0 in enumerate(Phi0_t):
    # Two Gaussian Φ clouds
    R1_2 = (X - x1)**2 + (Y - y1)**2
    R2_2 = (X - x2)**2 + (Y - y2)**2
    Phi = Phi0 * np.exp(-R1_2/(2*sigma**2)) + Phi0 * np.exp(-R2_2/(2*sigma**2))
    
    # Thresholded mapping
    rho = np.where(Phi > Phi_c, Phi, 0.0)
    M = rho.sum() * cell_area
    
    # Potential and gravity only for selected snapshots
    if idx in snapshot_indices:
        if rho.max() > 0:
            V = solve_potential(rho)
            gmag = gravity_mag(V)
        else:
            V = np.zeros_like(Phi)
            gmag = np.zeros_like(Phi)
        snapshots.append((t_arr[idx], Phi.copy(), rho.copy(), V.copy(), gmag.copy(), M))

# Plot snapshots: Φ, ρ, and V for each time
n_snaps = len(snapshots)
fig, axes = plt.subplots(n_snaps, 3, figsize=(11, 3*n_snaps))

for i, (t, Phi_s, rho_s, V_s, gmag_s, M_s) in enumerate(snapshots):
    ax0 = axes[i, 0]
    im0 = ax0.imshow(Phi_s, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
    ax0.set_title(f"t={t:.2f}: Φ (M≈{M_s:.2f})")
    fig.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)
    
    ax1 = axes[i, 1]
    im1 = ax1.imshow(rho_s, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
    ax1.set_title("ρ(Φ) (gravitating region)")
    fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    
    ax2 = axes[i, 2]
    im2 = ax2.imshow(V_s, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
    ax2.set_title("V (gravitational potential)")
    fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

for j in range(3):
    axes[-1, j].set_xlabel("x")

plt.tight_layout()
plt.show()
