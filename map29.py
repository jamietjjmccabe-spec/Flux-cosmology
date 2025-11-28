import numpy as np
import matplotlib.pyplot as plt

# ============================================================
#  CONFIG
# ============================================================

# Grid size
NX, NY = 256, 256

# Physical size of the box (arbitrary units)
Lx, Ly = 10.0, 10.0

# Photon deficit parameters (two pits = binary)
A_deficit = 1.0   # amplitude of photon deficit
R_core    = 0.7   # "size" of each deficit pit
d_sep     = 2.5   # separation between the two wells (centers at ±d_sep/2)

# Flux and matter coupling constants (tunable)
kappa_gamma = 1.0   # how strongly photon deficit sources Φ
alpha_m     = 1.0   # how strongly curvature of Δn_gamma gives matter

# ============================================================
#  GRID
# ============================================================

x = np.linspace(-Lx/2, Lx/2, NX, endpoint=False)
y = np.linspace(-Ly/2, Ly/2, NY, endpoint=False)
X, Y = np.meshgrid(x, y, indexing="ij")

dx = x[1] - x[0]
dy = y[1] - y[0]

# ============================================================
#  PHOTON DEFICIT FIELD Δn_gamma (binary pits)
# ============================================================

# Centers of the two deficits
x1, y1 = -d_sep/2.0, 0.0
x2, y2 =  d_sep/2.0, 0.0

# Gaussian photon deficits
def gaussian_deficit(X, Y, xc, yc, R, A):
    return A * np.exp(-((X - xc)**2 + (Y - yc)**2) / (2.0 * R**2))

Delta_n1 = gaussian_deficit(X, Y, x1, y1, R_core, A_deficit)
Delta_n2 = gaussian_deficit(X, Y, x2, y2, R_core, A_deficit)

Delta_n = Delta_n1 + Delta_n2   # total photon deficit field (two wells)

# ============================================================
#  FFT-BASED POISSON SOLVER FOR Φ:
#      ∇^2 Φ = - kappa_gamma * Δn
#  (periodic BCs in this toy model)
# ============================================================

# Wave numbers for FFT (2π n / L)
kx = 2.0 * np.pi * np.fft.fftfreq(NX, d=dx)
ky = 2.0 * np.pi * np.fft.fftfreq(NY, d=dy)
KX, KY = np.meshgrid(kx, ky, indexing="ij")

k2 = KX**2 + KY**2

Delta_n_k = np.fft.fftn(Delta_n)

Phi_k = np.zeros_like(Delta_n_k, dtype=np.complex128)

# Avoid division by zero at k=0 (set mean of Φ to zero)
mask_nonzero = (k2 != 0)
Phi_k[mask_nonzero] = kappa_gamma * Delta_n_k[mask_nonzero] / k2[mask_nonzero]

Phi = np.fft.ifftn(Phi_k).real

# ============================================================
#  LAPLACIAN OPERATOR (FINITE DIFFERENCE) FOR MATTER DENSITY
# ============================================================

def laplacian(field, dx, dy):
    """
    2D Laplacian with periodic boundary conditions.
    """
    f = field
    lap = (
        (np.roll(f, +1, axis=0) - 2.0*f + np.roll(f, -1, axis=0)) / dx**2 +
        (np.roll(f, +1, axis=1) - 2.0*f + np.roll(f, -1, axis=1)) / dy**2
    )
    return lap

# Matter density from negative curvature of photon deficit
lap_Delta_n = laplacian(Delta_n, dx, dy)

rho_m = alpha_m * np.maximum(-lap_Delta_n, 0.0)

# Normalise for nicer plotting (optional)
rho_m /= rho_m.max() + 1e-12
Phi   -= Phi.min()  # shift so minimum is zero for visualization

# ============================================================
#  PLOTTING
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# 1) Photon deficit field Δn_gamma
im0 = axes[0].imshow(
    Delta_n.T,
    origin="lower",
    extent=[x.min(), x.max(), y.min(), y.max()],
    aspect="equal"
)
axes[0].set_title("Photon deficit Δnγ(x,y)\n(two pits = binary)")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")
fig.colorbar(im0, ax=axes[0], shrink=0.8, label="Δnγ")

# 2) Flux potential Φ
im1 = axes[1].imshow(
    Phi.T,
    origin="lower",
    extent=[x.min(), x.max(), y.min(), y.max()],
    aspect="equal"
)
axes[1].set_title("Flux potential Φ(x,y)\nfrom ∇²Φ = -κγ Δnγ")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")
fig.colorbar(im1, ax=axes[1], shrink=0.8, label="Φ")

# 3) Matter density ρ_m
im2 = axes[2].imshow(
    rho_m.T,
    origin="lower",
    extent=[x.min(), x.max(), y.min(), y.max()],
    aspect="equal",
    vmin=0.0, vmax=1.0
)
axes[2].set_title("Matter density ρ_m(x,y)\n~ [-∇²Δnγ]_+")
axes[2].set_xlabel("x")
axes[2].set_ylabel("y")
fig.colorbar(im2, ax=axes[2], shrink=0.8, label="ρ_m (norm)")

plt.tight_layout()
plt.show()
