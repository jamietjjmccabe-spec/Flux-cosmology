import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d

# =====================================================================
# 1. PARAMETERS & CONSTANTS
# =====================================================================
c = 299792.458 # km/s

# LambdaCDM (Planck 2018 baseline)
H0_lcdm = 67.4
Om_lcdm = 0.315
Ol_lcdm = 1.0 - Om_lcdm

# Flux-CBH (Best Matter-Float Fit)
H0_flux = 73.0
p_flux = 0.214
omega_m_flux = 0.144
omega_r_flux = 4.183e-5
# Convert physical densities to dimensionless fractions for Flux
h_flux = H0_flux / 100.0
Om_flux = omega_m_flux / (h_flux**2)
Or_flux = omega_r_flux / (h_flux**2)
Oflux_0 = 1.0 - Om_flux - Or_flux # Flatness

# =====================================================================
# 2. LOAD CBH ENTROPY DATA
# =====================================================================
df = pd.read_csv('flux_cbh_schechter_entropy.csv').sort_values('z')
z_data = df['z'].values
S_CBH = df['S_CBH_density_kB_mpc3'].values

S_CBH_interp = interp1d(z_data, S_CBH, kind='linear', fill_value=(S_CBH[0], S_CBH[-1]), bounds_error=False)
S_CBH_z0 = S_CBH_interp(0.0)

def F_CBH(z):
    return S_CBH_interp(z) / S_CBH_z0

# =====================================================================
# 3. KINEMATIC FUNCTIONS
# =====================================================================
def E_lcdm(z):
    return np.sqrt(Om_lcdm*(1+z)**3 + Ol_lcdm)

def E_flux(z):
    flux_term = Oflux_0 * (F_CBH(z)**p_flux)
    return np.sqrt(Or_flux*(1+z)**4 + Om_flux*(1+z)**3 + flux_term)

def integrand_Dc_lcdm(z):
    return 1.0 / E_lcdm(z)

def integrand_Dc_flux(z):
    return 1.0 / E_flux(z)

def distance_modulus(z, model='lcdm'):
    if z == 0:
        return 0
        
    if model == 'lcdm':
        val, _ = quad(integrand_Dc_lcdm, 0, z, limit=500)
        Dc = (c / H0_lcdm) * val
    else:
        val, _ = quad(integrand_Dc_flux, 0, z, limit=500)
        Dc = (c / H0_flux) * val
        
    Dl = (1 + z) * Dc # Luminosity distance
    
    # Distance modulus: mu = 5 log10(Dl_Mpc) + 25
    return 5.0 * np.log10(Dl) + 25.0

# =====================================================================
# 4. COMPUTE OVER SUPERNOVA DOMAIN
# =====================================================================
print(f"Computing distance moduli over SN domain (0 < z < 2.3)...")
z_sn = np.linspace(0.01, 2.3, 150) # Start slightly above 0 to avoid log(0)

mu_lcdm = np.array([distance_modulus(z, 'lcdm') for z in z_sn])
mu_flux = np.array([distance_modulus(z, 'flux') for z in z_sn])

# Residual: Delta mu
delta_mu = mu_flux - mu_lcdm

# =====================================================================
# 5. VISUALIZATION
# =====================================================================
plt.style.use('dark_background')
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

# Top panel: Absolute Distance Modulus
ax1.plot(z_sn, mu_lcdm, color='white', linestyle='--', linewidth=2, label=r'Planck $\Lambda$CDM ($H_0=67.4$)')
ax1.plot(z_sn, mu_flux, color='cyan', linewidth=2, label=r'Flux-CBH ($H_0=73.0, p=0.214$)')
ax1.set_ylabel(r'Distance Modulus $\mu(z)$', fontsize=12)
ax1.set_title(r'Supernova Domain Diagnostics: Flux-CBH vs $\Lambda$CDM', fontsize=14, pad=10)
ax1.legend(loc='lower right', fontsize=10)
ax1.grid(True, alpha=0.2)

# Bottom panel: Residual
ax2.axhline(0, color='white', linestyle='--', alpha=0.5)
ax2.plot(z_sn, delta_mu, color='magenta', linewidth=2.5, label=r'$\Delta\mu$ (Flux - $\Lambda$CDM)')

# Typical Pantheon observational scatter bounds (approximate reference limits)
ax2.axhspan(-0.1, 0.1, color='grey', alpha=0.2, label=r'Typical SN Scatter Bound ($\pm 0.1$ mag)')

ax2.set_xlabel(r'Redshift ($z$)', fontsize=12)
ax2.set_ylabel(r'$\Delta\mu$ [mag]', fontsize=12)
ax2.legend(loc='lower right', fontsize=10)
ax2.grid(True, alpha=0.2)

plt.xlim(0, 2.3)
plt.tight_layout()
plt.savefig('flux_cbh_sn_residual_scan.png', dpi=300)
print("Saved 'flux_cbh_sn_residual_scan.png'.")

# Output a quick summary
print(f"Max residual: {np.max(delta_mu):.3f} mag")
print(f"Min residual: {np.min(delta_mu):.3f} mag")