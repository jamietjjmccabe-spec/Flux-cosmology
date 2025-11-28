import numpy as np
import matplotlib.pyplot as plt

# Time-dependent toy model: collapsing Φ cloud -> turning on gravity

# Spatial grid
L = 10.0
N = 300
x = np.linspace(-L/2, L/2, N)
y = np.linspace(-L/2, L/2, N)
X, Y = np.meshgrid(x, y)
R2 = X**2 + Y**2

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

# Time steps: simulate collapse by decreasing sigma (cloud size) over time
times = [0, 1, 2, 3]
sigmas = [3.0, 2.0, 1.5, 1.0]  # cloud shrinks, peak Φ rises
Phi_c = 0.4                    # threshold for "turning on" gravitating mass

fig, axes = plt.subplots(len(times), 3, figsize=(10, 12))

for i, (t, sigma) in enumerate(zip(times, sigmas)):
    # Low-Φ cloud collapsing into high-Φ core
    Phi = np.exp(-R2 / (2 * sigma**2))  # peak fixed at 1, width shrinks
    
    # Thresholded mapping: only Φ > Phi_c gravitates
    rho = np.where(Phi > Phi_c, Phi, 0.0)
    
    # Solve for potential and gravity
    if rho.max() <= 0:
        V = np.zeros_like(Phi)
        gmag = np.zeros_like(Phi)
    else:
        V = solve_potential(rho)
        gmag = gravity_mag(V)
    
    # Plot Φ
    ax0 = axes[i, 0]
    im0 = ax0.imshow(Phi, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
    ax0.set_title(f"t={t}: Φ (σ={sigma})")
    ax0.set_ylabel("y")
    fig.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)
    
    # Plot rho (where gravity is sourced)
    ax1 = axes[i, 1]
    im1 = ax1.imshow(rho, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
    ax1.set_title("ρ(Φ) (gravitating core)")
    fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    
    # Plot |g|
    ax2 = axes[i, 2]
    im2 = ax2.imshow(gmag, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower')
    ax2.set_title("|g| (gravity strength)")
    fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

for j in range(3):
    axes[-1, j].set_xlabel("x")

plt.tight_layout()
plt.show()
