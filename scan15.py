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

# Flux-CBH BAO-Rescued Best Fit
H0_flux = 71.0
omega_m_flux = 0.145
h_flux = H0_flux / 100.0
Om_flux = omega_m_flux / h_flux**2
Or_flux = 4.183e-5 / h_flux**2
Oflux_0 = 1.0 - Om_flux - Or_flux

p0 = 0.19
p_high = 0.25
zt = 1.70
w = 0.35

# =====================================================================
# 2. RSD DATA (Consensus Compilation)
# =====================================================================
rsd_data = np.array([
    [0.067, 0.423, 0.055], [0.17,  0.510, 0.060], [0.22,  0.420, 0.070],
    [0.25,  0.351, 0.058], [0.32,  0.384, 0.095], [0.38,  0.440, 0.060],
    [0.41,  0.450, 0.040], [0.57,  0.427, 0.066], [0.60,  0.430, 0.040],
    [0.73,  0.437, 0.072], [0.86,  0.480, 0.100], [1.40,  0.482, 0.116],
    [1.52,  0.420, 0.076]
])
z_obs = rsd_data[:, 0]
fsigma8_obs = rsd_data[:, 1]
err_obs = rsd_data[:, 2]

# =====================================================================
# 3. LOAD 4-TERM FLUX DATA
# =====================================================================
df = pd.read_csv('flux_cbh_schechter_entropy.csv').sort_values('z')
source_col = "Omega_flux_4term_fit"

F_interp = interp1d(
    df["z"].values, df[source_col].values, 
    fill_value=(df[source_col].values[0], df[source_col].values[-1]), 
    bounds_error=False
)
F0 = F_interp(0.0)

def F_raw(z):
    return np.maximum(F_interp(z) / F0, 1e-90)

# =====================================================================
# 4. BACKGROUND KINEMATICS
# =====================================================================
def E_lcdm(z):
    return np.sqrt(Or_lcdm*(1+z)**4 + Om_lcdm*(1+z)**3 + (1.0 - Om_lcdm - Or_lcdm))

def p_eff(z):
    return p_high + (p0 - p_high) / (1.0 + np.exp((z - zt) / w))

def E_flux(z):
    return np.sqrt(Or_flux*(1+z)**4 + Om_flux*(1+z)**3 + Oflux_0 * (F_raw(z)**p_eff(z)))

def dlnH_dN_lcdm(z):
    dz = 1e-4
    return - (1+z) / E_lcdm(z) * ((E_lcdm(z+dz) - E_lcdm(z-dz)) / (2*dz))

def dlnH_dN_flux(z):
    dz = 1e-4
    return - (1+z) / E_flux(z) * ((E_flux(z+dz) - E_flux(z-dz)) / (2*dz))

# =====================================================================
# 5. LINEAR GROWTH ODE
# =====================================================================
def ode_lcdm(N, Y):
    y, v = Y; z = np.exp(-N) - 1.0
    return [v, - (2.0 + dlnH_dN_lcdm(z)) * v + 1.5 * (Om_lcdm * (1+z)**3 / E_lcdm(z)**2) * y]

def ode_flux(N, Y):
    y, v = Y; z = np.exp(-N) - 1.0
    return [v, - (2.0 + dlnH_dN_flux(z)) * v + 1.5 * (Om_flux * (1+z)**3 / E_flux(z)**2) * y]

z_start = 100.0
N_start, N_end = -np.log(1.0 + z_start), 0.0
Y0 = [np.exp(N_start), np.exp(N_start)]
N_eval = np.linspace(N_start, N_end, 1000)
z_eval = np.exp(-N_eval) - 1.0

sol_lcdm = solve_ivp(ode_lcdm, [N_start, N_end], Y0, t_eval=N_eval, rtol=1e-8, atol=1e-8)
sol_flux = solve_ivp(ode_flux, [N_start, N_end], Y0, t_eval=N_eval, rtol=1e-8, atol=1e-8)

D_lcdm, v_lcdm = sol_lcdm.y
D_flux, v_flux = sol_flux.y

f_lcdm = v_lcdm / D_lcdm
f_flux = v_flux / D_flux

# Anchor sigma8
sigma8_flux = sigma8_lcdm * (D_flux[-1] / D_lcdm[-1])
S8_flux = sigma8_flux * np.sqrt(Om_flux / 0.3)

fsig8_lcdm = f_lcdm * sigma8_lcdm * (D_lcdm / D_lcdm[-1])
fsig8_flux = f_flux * sigma8_flux * (D_flux / D_flux[-1])

# Create sorted interpolators for chi^2
order = np.argsort(z_eval)
fsig8_lcdm_interp = interp1d(z_eval[order], fsig8_lcdm[order], kind='cubic', bounds_error=False, fill_value="extrapolate")
fsig8_flux_interp = interp1d(z_eval[order], fsig8_flux[order], kind='cubic', bounds_error=False, fill_value="extrapolate")

# =====================================================================
# 6. CHI-SQUARE CALCULATION
# =====================================================================
pred_lcdm = fsig8_lcdm_interp(z_obs)
pred_flux = fsig8_flux_interp(z_obs)

chi2_lcdm = np.sum(((pred_lcdm - fsigma8_obs) / err_obs)**2)
chi2_flux = np.sum(((pred_flux - fsigma8_obs) / err_obs)**2)

print("=== FINAL BAO-RESCUED GROWTH & RSD RESULTS ===")
print(f"\nModel Parameters: H0={H0_flux:.1f}, Om_m={Om_flux:.4f}, p0={p0}, p_high={p_high}")
print("-" * 65)
print(f"{'Model':<15} | {'Omega_m':<10} | {'D(z=0)':<10} | {'sigma_8':<10} | {'S_8':<10}")
print("-" * 65)
print(f"{'Planck LCDM':<15} | {Om_lcdm:<10.4f} | {D_lcdm[-1]:<10.4f} | {sigma8_lcdm:<10.4f} | {S8_lcdm:<10.4f}")
print(f"{'Flux-CBH':<15} | {Om_flux:<10.4f} | {D_flux[-1]:<10.4f} | {sigma8_flux:<10.4f} | {S8_flux:<10.4f}")
print(f"Shift in S_8 : {(S8_flux - S8_lcdm):.4f}")
print("-" * 65)
print(f"RSD Chi^2 (LCDM): {chi2_lcdm:.2f}")
print(f"RSD Chi^2 (Flux): {chi2_flux:.2f}")
print(f"Delta Chi^2     : {(chi2_flux - chi2_lcdm):.2f}")

# =====================================================================
# 7. VISUALIZATION
# =====================================================================
plt.style.use('dark_background')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Panel 1: S8/sigma8 text
ax1.plot(z_eval, D_lcdm / np.exp(N_eval), color='gray', linestyle='--', linewidth=2, label=r'$\Lambda$CDM')
ax1.plot(z_eval, D_flux / np.exp(N_eval), color='cyan', linewidth=2.5, label='Flux-CBH (BAO-Rescued)')
ax1.set_xlabel('Redshift (z)', fontsize=12)
ax1.set_ylabel(r'$D(z) / a$', fontsize=12)
ax1.set_title(r'Normalized Linear Growth Factor', fontsize=14)
ax1.set_xlim(0, 10)
ax1.legend()
ax1.grid(alpha=0.2)
ax1.text(0.5, 0.5, f"$\Delta S_8 = {(S8_flux - S8_lcdm):.3f}$", transform=ax1.transAxes, fontsize=14, color='white')

# Panel 2: RSD Fit
ax2.errorbar(z_obs, fsigma8_obs, yerr=err_obs, fmt='o', color='white', ecolor='gray', capsize=4, label='RSD Data')
ax2.plot(z_eval, fsig8_lcdm, color='gray', linestyle='--', linewidth=2, label=rf'$\Lambda$CDM ($\chi^2={chi2_lcdm:.1f}$)')
ax2.plot(z_eval, fsig8_flux, color='magenta', linewidth=3, label=rf'Flux-CBH ($\chi^2={chi2_flux:.1f}$)')

ax2.set_xlabel('Redshift (z)', fontsize=12)
ax2.set_ylabel(r'$f(z) \sigma_8(z)$', fontsize=12)
ax2.set_title(r'Redshift-Space Distortions: Growth Rate', fontsize=14)
ax2.set_xlim(0, 1.8)
ax2.legend(loc='lower left', fontsize=12)
ax2.grid(alpha=0.2)

plt.tight_layout()
plt.savefig('scan16_bao_rescued_growth.png', dpi=300)
print("\nSaved 'scan16_bao_rescued_growth.png'.")