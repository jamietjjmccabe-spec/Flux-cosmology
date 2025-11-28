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
sigma = 2.0  # fixed cloud size

# Time evolution for peak Φ0(t): from 0.1 to 1.0
T = 80
t_arr = np.linspace(0, 1, T)
Phi0_start = 0.1
Phi0_end = 1.0
Phi0_t = Phi0_start + (Phi0_end - Phi0_start) * t_arr

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

total_mass = []
g_center = []

for Phi0 in Phi0_t:
    Phi = Phi0 * np.exp(-R2 / (2 * sigma**2))
    rho = np.where(Phi > Phi_c, Phi, 0.0)
    M = rho.sum() * cell_area
    total_mass.append(M)
    
    if rho.max() <= 0:
        g_center.append(0.0)
    else:
        V = solve_potential(rho)
        gmag = gravity_mag(V)
        # take gravity magnitude at center
        ix0 = np.argmin(np.abs(x - 0.0))
        iy0 = np.argmin(np.abs(y - 0.0))
        g_center.append(gmag[iy0, ix0])

total_mass = np.array(total_mass)
g_center = np.array(g_center)

# Find switch-on point: first index where M>0
on_idx = np.argmax(total_mass > 0)

plt.figure(figsize=(8,5))
plt.plot(t_arr, total_mass, label="Total gravitating mass M(t)")
plt.plot(t_arr, g_center / g_center.max() * total_mass.max(), 
         label="Central gravity |g(0)| (rescaled)", linestyle='--')
if on_idx > 0:
    plt.axvline(t_arr[on_idx], color='red', linestyle=':', 
                label=f'Gravity switches on at t≈{t_arr[on_idx]:.2f}')
plt.xlabel("t (normalized time)")
plt.ylabel("Arbitrary units")
plt.title("Gravity switching on as Φ-peak rises across threshold")
plt.legend()
plt.tight_layout()
plt.show()

total_mass[0], total_mass[-1], g_center[0], g_center[-1], t_arr[on_idx]
