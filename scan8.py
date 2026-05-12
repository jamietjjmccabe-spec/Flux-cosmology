import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d

# =====================================================================
# 1. PARAMETERS & CONSTANTS
# =====================================================================
c = 299792.458 # km/s
z_star = 1089.92
theta_star_planck = 0.010411

# Planck Baseline for LCDM
H0_lcdm = 67.4
Om_lcdm = 0.315
Ol_lcdm = 1.0 - Om_lcdm

# Flux-CBH Base (from your best matter-float fit)
H0_flux = 73.0
omega_m_flux = 0.144
omega_b_flux = 0.02237
omega_r_flux = 4.183e-5

h_flux = H0_flux / 100.0
Om_flux = omega_m_flux / (h_flux**2)
Or_flux = omega_r_flux / (h_flux**2)
Oflux_0 = 1.0 - Om_flux - Or_flux
R_b_factor = (3.0 * omega_b_flux) / (4.0 * omega_gamma_flux) if 'omega_gamma_flux' in locals() else (3.0 * 0.02237) / (4.0 * 2.469e-5)

# =====================================================================
# 2. LOAD CBH ENTROPY DATA
# =====================================================================
print("Loading Schechter Entropy data...")
df = pd.read_csv('flux_cbh_schechter_entropy.csv').sort_values('z')
S_CBH_interp = interp1d(df['z'].values, df['S_CBH_density_kB_mpc3'].values, kind='linear', fill_value=(df['S_CBH_density_kB_mpc3'].values[0], 0.0), bounds_error=False)
S_CBH_z0 = S_CBH_interp(0.0)

def F_CBH(z):
    # Ensure it never drops exactly to zero to avoid log/power errors
    return max(S_CBH_interp(z) / S_CBH_z0, 1e-90)

# =====================================================================
# 3. KINEMATIC FUNCTIONS
# =====================================================================
def p_dynamic(z, p0, dp, zt, w=1.0):
    """Dynamic compression envelope"""
    return p0 + dp * (1.0 - np.exp(- (z / zt)**w ))

def E_flux(z, p0, dp, zt):
    p_z = p_dynamic(z, p0, dp, zt)
    flux_term = Oflux_0 * (F_CBH(z)**p_z)
    return np.sqrt(Or_flux*(1+z)**4 + Om_flux*(1+z)**3 + flux_term)

def sound_speed(z):
    return c / np.sqrt(3.0 * (1.0 + (R_b_factor / (1.0 + z))))

def integrand_rs(z, p0, dp, zt):
    return sound_speed(z) / (100.0 * h_flux * E_flux(z, p0, dp, zt))

def integrand_Dm(z, p0, dp, zt):
    return c / (100.0 * h_flux * E_flux(z, p0, dp, zt))

def evaluate_cmb(p0, dp, zt):
    D_M, _ = quad(integrand_Dm, 0, z_star, args=(p0, dp, zt), limit=300)
    r_s, _ = quad(integrand_rs, z_star, 1e5, args=(p0, dp, zt), limit=300)
    return r_s / D_M

# SN Pre-compute for baseline LCDM
z_sn_eval = np.linspace(0.01, 2.3, 25)
def E_lcdm(z): return np.sqrt(Om_lcdm*(1+z)**3 + Ol_lcdm)
mu_lcdm = np.array([5.0 * np.log10((1+z) * (c/H0_lcdm) * quad(lambda x: 1/E_lcdm(x), 0, z)[0]) + 25.0 for z in z_sn_eval])

def evaluate_sn_shape(p0, dp, zt):
    mu_flux = np.array([5.0 * np.log10((1+z) * (c/H0_flux) * quad(integrand_Dm, 0, z, args=(p0, dp, zt))[0] * (100.0*h_flux)/c) + 25.0 for z in z_sn_eval])
    delta_mu = mu_flux - mu_lcdm
    offset = np.mean(delta_mu[(z_sn_eval >= 0.01) & (z_sn_eval <= 0.05)])
    return np.max(np.abs(delta_mu - offset))

# =====================================================================
# 4. GRID SCAN WITH SHORT-CIRCUIT
# =====================================================================
print("Starting 3D Scan (p0, delta_p, z_t)...")
p0_vals = np.linspace(0.15, 0.35, 12)
dp_vals = np.linspace(-0.15, 0.15, 12)
zt_vals = np.linspace(0.3, 1.5, 10)

results = []

for p0 in p0_vals:
    for dp in dp_vals:
        for zt in zt_vals:
            # Step 1: Check CMB first
            theta = evaluate_cmb(p0, dp, zt)
            err_ppm = abs(theta - theta_star_planck) / theta_star_planck * 1e6
            
            # Short-circuit: Only check SN if CMB is strictly preserved (< 800 ppm)
            if err_ppm < 800:
                sn_shape_max = evaluate_sn_shape(p0, dp, zt)
                results.append({'p0': p0, 'dp': dp, 'zt': zt, 'theta_err_ppm': err_ppm, 'sn_shape_max': sn_shape_max})

# =====================================================================
# 5. RESULTS SUMMARY
# =====================================================================
df_res = pd.DataFrame(results)
if not df_res.empty:
    df_res = df_res.sort_values('sn_shape_max')
    print("\n=== TOP 5 GATED MODELS ===")
    print(df_res.head(5).to_string(index=False))
    df_res.to_csv('flux_cbh_gated_scan_results.csv', index=False)
else:
    print("No models passed the CMB constraint filter.")