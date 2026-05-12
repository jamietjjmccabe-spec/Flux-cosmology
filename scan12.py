import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# =====================================================================
# 1. PARAMETERS & CONSTANTS
# =====================================================================
# Planck 2018 LCDM Baseline
H0_lcdm = 67.4
Om_lcdm = 0.315
Or_lcdm = 4.183e-5 / (H0_lcdm/100.0)**2
sigma8_lcdm = 0.811
S8_lcdm = sigma8_lcdm * np.sqrt(Om_lcdm / 0.3)

# Flux-CBH Best Joint Fit
H0_flux = 73.0
omega_m_flux = 0.140
h_flux = H0_flux / 100.0
Om_flux = omega_m_flux / h_flux**2
Or_flux = 4.183e-5 / h_flux**2
Oflux_0 = 1.0 - Om_flux - Or_flux

p0 = 0.18
p_high = 0.23
zt = 1.70
w = 0.35

# =====================================================================
# 2. LOAD 4-TERM FLUX DATA
# =====================================================================
df = pd.read_csv('flux_cbh_schechter_entropy.csv').sort_values('z')
source_col = "Omega_flux_4term_fit"

F_interp = interp1d(
    df["z"].values,
    df[source_col].values,
    fill_value=(df[source_col].values[0], df[source_col].values[-1]),
    bounds_error=False
)
F0 = F_interp(0.0)

def F_raw(z):
    return np.maximum(F_interp(z) / F0, 1e-90)

# =====================================================================
# 3. BACKGROUND KINEMATICS
# =====================================================================
def E_lcdm(z):
    return np.sqrt(Or_lcdm*(1+z)**4 + Om_lcdm*(1+z)**3 + (1.0 - Om_lcdm - Or_lcdm))

def p_eff(z):
    return p_high + (p0 - p_high) / (1.0 + np.exp((z - zt) / w))

def E_flux(z):
    flux_term = Oflux_0 * (F_raw(z)**p_eff(z))
    return np.sqrt(Or_flux*(1+z)**4 + Om_flux*(1+z)**3 + flux_term)

# Derivatives using central difference
def dlnH_dN_lcdm(z):
    dz = 1e-4
    dE_dz = (E_lcdm(z+dz) - E_lcdm(z-dz)) / (2*dz)
    return - (1+z) / E_lcdm(z) * dE_dz

def dlnH_dN_flux(z):
    dz = 1e-4
    dE_dz = (E_flux(z+dz) - E_flux(z-dz)) / (2*dz)
    return - (1+z) / E_flux(z) * dE_dz

# =====================================================================
# 4. LINEAR GROWTH ODE
# =====================================================================
# Equation: D'' + (2 + dlnH/dN) D' - 1.5 * Om(a) * D = 0
# State vector Y = [D, dD/dN]

def ode_lcdm(N, Y):
    y, v = Y
    z = np.exp(-N) - 1.0
    Om_z = Om_lcdm * (1+z)**3 / E_lcdm(z)**2
    dv_dN = - (2.0 + dlnH_dN_lcdm(z)) * v + 1.5 * Om_z * y
    dy_dN = v
    return [dy_dN, dv_dN]

def ode_flux(N, Y):
    y, v = Y
    z = np.exp(-N) - 1.0
    Om_z = Om_flux * (1+z)**3 / E_flux(z)**2
    # Setting mu_eff = 1 for now (pure background effect)
    dv_dN = - (2.0 + dlnH_dN_flux(z)) * v + 1.5 * Om_z * y
    dy_dN = v
    return [dy_dN, dv_dN]

# =====================================================================
# 5. INTEGRATION
# =====================================================================
z_start = 100.0  # Deep matter domination
N_start = -np.log(1.0 + z_start)
N_end = 0.0

# Initial conditions: D ~ a => D = e^N, dD/dN = e^N
Y0 = [np.exp(N_start), np.exp(N_start)]
N_eval = np.linspace(N_start, N_end, 1000)
z_eval = np.exp(-N_eval) - 1.0

sol_lcdm = solve_ivp(ode_lcdm, [N_start, N_end], Y0, t_eval=N_eval, rtol=1e-8, atol=1e-8)
sol_flux = solve_ivp(ode_flux, [N_start, N_end], Y0, t_eval=N_eval, rtol=1e-8, atol=1e-8)

D_lcdm = sol_lcdm.y[0]
v_lcdm = sol_lcdm.y[1]
D_flux = sol_flux.y[0]
v_flux = sol_flux.y[1]

# Growth rate f = dlnD / dln(a) = v / D
f_lcdm = v_lcdm / D_lcdm
f_flux = v_flux / D_flux

# =====================================================================
# 6. S8 AND SIGMA_8 CALCULATION
# =====================================================================
# We anchor the primordial power to match LCDM. Thus, sigma_8 is scaled 
# strictly by the ratio of the linear growth factors today.
ratio_D0 = D_flux[-1] / D_lcdm[-1]
sigma8_flux = sigma8_lcdm * ratio_D0
S8_flux = sigma8_flux * np.sqrt(Om_flux / 0.3)

print("=== GROWTH OF STRUCTURE RESULTS ===")
print(f"{'Model':<15} | {'Omega_m':<10} | {'D(z=0)':<10} | {'sigma_8':<10} | {'S_8':<10}")
print("-" * 65)
print(f"{'Planck LCDM':<15} | {Om_lcdm:<10.4f} | {D_lcdm[-1]:<10.4f} | {sigma8_lcdm:<10.4f} | {S8_lcdm:<10.4f}")
print(f"{'Flux-CBH':<15} | {Om_flux:<10.4f} | {D_flux[-1]:<10.4f} | {sigma8_flux:<10.4f} | {S8_flux:<10.4f}")
print(f"Shift in S_8 : {(S8_flux - S8_lcdm):.4f}")

# =====================================================================
# 7. VISUALIZATION
# =====================================================================
fsig8_lcdm = f_lcdm * sigma8_lcdm * (D_lcdm / D_lcdm[-1])
fsig8_flux = f_flux * sigma8_flux * (D_flux / D_flux[-1])

plt.style.use('dark_background')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Panel 1: Normalized Growth Factor D(z)/a
ax1.plot(z_eval, D_lcdm / np.exp(N_eval), color='white', linestyle='--', linewidth=2, label=r'$\Lambda$CDM')
ax1.plot(z_eval, D_flux / np.exp(N_eval), color='cyan', linewidth=2.5, label='Flux-CBH')
ax1.set_xlabel('Redshift (z)', fontsize=12)
ax1.set_ylabel(r'$D(z) / a$', fontsize=12)
ax1.set_title('Normalized Linear Growth Factor', fontsize=14)
ax1.set_xlim(0, 10)
ax1.legend()
ax1.grid(alpha=0.2)

# Panel 2: Growth Rate f*sigma_8
ax2.plot(z_eval, fsig8_lcdm, color='white', linestyle='--', linewidth=2, label=r'$\Lambda$CDM')
ax2.plot(z_eval, fsig8_flux, color='magenta', linewidth=2.5, label='Flux-CBH')

# Approximate typical data scatter bounds for f*sigma8
ax2.fill_between(z_eval, fsig8_lcdm*0.9, fsig8_lcdm*1.1, color='grey', alpha=0.15, label=r'Typical $\Lambda$CDM scatter')

ax2.set_xlabel('Redshift (z)', fontsize=12)
ax2.set_ylabel(r'$f(z) \sigma_8(z)$', fontsize=12)
ax2.set_title('Observable Growth Rate', fontsize=14)
ax2.set_xlim(0, 2)
ax2.legend()
ax2.grid(alpha=0.2)

plt.tight_layout()
plt.savefig('scan12_growth_history.png', dpi=300)
print("\nSaved 'scan12_growth_history.png'.")