import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
import time

# =====================================================================
# 1. PARAMETERS & CONSTANTS
# =====================================================================
c = 299792.458
z_star = 1089.92
z_d = 1059.39
theta_star_planck = 0.0104132

omega_b_flux = 0.02237
omega_gamma_flux = 2.469e-5
R_b_factor = (3.0 * omega_b_flux) / (4.0 * omega_gamma_flux)

zt_fixed = 1.70
w_fixed = 0.35

# =====================================================================
# 2. LOAD 4-TERM FLUX DATA
# =====================================================================
df = pd.read_csv('flux_cbh_schechter_entropy.csv').sort_values('z')
source_col = "Omega_flux_4term_fit"

F_interp = interp1d(df["z"].values, df[source_col].values, fill_value=(df[source_col].values[0], df[source_col].values[-1]), bounds_error=False)
F0 = F_interp(0.0)

if F0 <= 0 or not np.isfinite(F0):
    raise ValueError(f"Invalid F0={F0}")

def F_raw(z):
    return np.maximum(F_interp(z) / F0, 1e-90)

# =====================================================================
# 3. BAO CONSENSUS DATA
# =====================================================================
bao_data = [
    ['DV', 0.106, 2.976, 0.133, '6dFGS'],
    ['DV', 0.150, 4.466, 0.168, 'SDSS MGS'],
    ['DM', 0.380, 10.23, 0.17,  'BOSS DR12'],
    ['DH', 0.380, 25.00, 0.76,  'BOSS DR12'],
    ['DM', 0.510, 13.36, 0.21,  'BOSS DR12'],
    ['DH', 0.510, 22.33, 0.58,  'BOSS DR12'],
    ['DM', 0.610, 15.36, 0.27,  'BOSS DR12'],
    ['DH', 0.610, 20.53, 0.60,  'BOSS DR12'],
    ['DM', 0.698, 17.86, 0.33,  'eBOSS LRG'],
    ['DH', 0.698, 19.33, 0.53,  'eBOSS LRG'],
    ['DM', 1.480, 30.69, 0.80,  'eBOSS QSO'],
    ['DH', 1.480, 13.26, 0.55,  'eBOSS QSO'],
    ['DM', 2.330, 37.41, 1.86,  'eBOSS LyA'],
    ['DH', 2.330,  8.94, 0.22,  'eBOSS LyA']
]

# =====================================================================
# 4. PRE-COMPUTE LCDM SN BASELINES
# =====================================================================
z_sn_eval = np.linspace(0.01, 2.3, 30)
anchor_mask = (z_sn_eval >= 0.01) & (z_sn_eval <= 0.05)
lcdm_baselines = {}
Om_refs = np.linspace(0.24, 0.34, 11)
H0_baseline = 73.0

for Om_ref in Om_refs:
    Ol_ref = 1.0 - Om_ref
    def E_lcdm(z, Om=Om_ref, Ol=Ol_ref): return np.sqrt(Om*(1+z)**3 + Ol)
    mu = np.array([5.0 * np.log10((1+z) * (c/H0_baseline) * quad(lambda x: 1/E_lcdm(x), 0, z)[0]) + 25.0 for z in z_sn_eval])
    lcdm_baselines[Om_ref] = mu

# =====================================================================
# 5. KINEMATICS & EVALUATION
# =====================================================================
def E_flux(z, H0, omega_m, p0, p_high):
    h = H0 / 100.0
    Om = omega_m / h**2
    Or = 4.183e-5 / h**2
    Oflux_0 = 1.0 - Om - Or
    p_z = p_high + (p0 - p_high) / (1.0 + np.exp((z - zt_fixed) / w_fixed))
    return np.sqrt(Or*(1+z)**4 + Om*(1+z)**3 + Oflux_0 * (F_raw(z)**p_z))

def cs_flux(z):
    return c / np.sqrt(3.0 * (1.0 + (R_b_factor / (1.0 + z))))

def evaluate_model(H0, omega_m, p0, p_high):
    h = H0 / 100.0
    
    # Integrands
    int_Dm = lambda z: c / (100.0 * h * E_flux(z, H0, omega_m, p0, p_high))
    int_rs = lambda z: cs_flux(z) / (100.0 * h * E_flux(z, H0, omega_m, p0, p_high))
    
    # 1. CMB Gate
    D_M_star, _ = quad(int_Dm, 0, z_star, limit=300)
    r_s_star, _ = quad(int_rs, z_star, 1e6, limit=300)
    theta_star = r_s_star / D_M_star
    err_ppm = abs(theta_star - theta_star_planck) / theta_star_planck * 1e6
    
    if err_ppm > 2000:
        return err_ppm, np.nan, np.nan, np.nan

    # 2. BAO Gate
    r_d, _ = quad(int_rs, z_d, 1e6, limit=300)
    chi2_bao = 0.0
    dm_cache = {}
    
    for row in bao_data:
        obs_type, z, val_obs, err, survey = row
        if z not in dm_cache:
            Dm_z, _ = quad(int_Dm, 0, z, limit=300)
            Dh_z = c / (100.0 * h * E_flux(z, H0, omega_m, p0, p_high))
            Dv_z = (z * Dm_z**2 * Dh_z)**(1/3.0)
            dm_cache[z] = {'DM': Dm_z, 'DH': Dh_z, 'DV': Dv_z}
        
        val_model = dm_cache[z][obs_type] / r_d
        chi2_bao += ((val_model - val_obs) / err)**2
        
    if chi2_bao > 50.0:  # Relaxing the filter to see the gradient
        return err_ppm, chi2_bao, np.nan, np.nan

    # 3. SN Shape
    mu_flux = np.array([5.0 * np.log10((1+z) * (c/H0) * quad(int_Dm, 0, z)[0] * 100.0*h/c) + 25.0 for z in z_sn_eval])
    best_shape_max = 999.0
    best_Om_ref = None
    for Om_ref, mu_lcdm in lcdm_baselines.items():
        delta = mu_flux - mu_lcdm
        offset = np.mean(delta[anchor_mask])
        shape_max = np.max(np.abs(delta - offset))
        if shape_max < best_shape_max:
            best_shape_max = shape_max
            best_Om_ref = Om_ref
            
    return err_ppm, chi2_bao, best_shape_max, best_Om_ref

# =====================================================================
# 6. SCAN GRID
# =====================================================================
print("Starting Scan 15: Joint BAO Rescue...")
start = time.time()

H0_vals = np.linspace(69, 74, 6)
om_vals = np.linspace(0.135, 0.150, 4)
p0_vals = np.linspace(0.16, 0.22, 7)
ph_vals = np.linspace(0.20, 0.26, 7)

results = []
tested = 0

for H0 in H0_vals:
    for om in om_vals:
        for p0 in p0_vals:
            for ph in ph_vals:
                tested += 1
                err_ppm, chi2_bao, sn_max, best_Om_ref = evaluate_model(H0, om, p0, ph)
                
                if not np.isnan(chi2_bao):
                    results.append({
                        'H0': H0, 'omega_m': om, 'Omega_m': om/(H0/100.0)**2,
                        'p0': p0, 'p_high': ph, 'theta_err_ppm': err_ppm,
                        'chi2_bao': chi2_bao, 'sn_shape_max': sn_max, 'best_Om_ref': best_Om_ref
                    })

print(f"Scan complete in {time.time()-start:.1f}s.")
print(f"Tested {tested} models. Evaluated full stats for {len([r for r in results if not np.isnan(r['sn_shape_max'])])} models.")

df_res = pd.DataFrame(results)

if not df_res.empty:
    df_valid = df_res.dropna(subset=['sn_shape_max']).copy()
    if not df_valid.empty:
        # Score emphasizes BAO heavily since that is the broken constraint
        df_valid['score'] = df_valid['chi2_bao'] + (df_valid['theta_err_ppm']/200.0)**2 + (df_valid['sn_shape_max']*50)**2
        df_valid = df_valid.sort_values('score')
        
        print("\n=== TOP 10 MODELS (BAO Rescue) ===")
        print(df_valid.head(10).to_string(index=False))
        df_valid.to_csv('scan15_bao_rescue_results.csv', index=False)
    else:
        print("No models passed the relaxed BAO gate.")
else:
    print("No models evaluated.")