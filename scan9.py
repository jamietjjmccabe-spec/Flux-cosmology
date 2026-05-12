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
theta_star_planck = 0.0104132

H0_flux = 73.0
omega_m_flux = 0.144
omega_b_flux = 0.02237
omega_r_flux = 4.183e-5
omega_gamma_flux = 2.469e-5

h_flux = H0_flux / 100.0
Om_flux = omega_m_flux / h_flux**2
Or_flux = omega_r_flux / h_flux**2
Oflux_0 = 1.0 - Om_flux - Or_flux
R_b_factor = (3.0 * omega_b_flux) / (4.0 * omega_gamma_flux)

# LCDM Baseline for SN Comparison
H0_lcdm = 67.4
Om_lcdm = 0.315
Ol_lcdm = 1.0 - Om_lcdm

# =====================================================================
# 2. LOAD CBH DATA
# =====================================================================
df = pd.read_csv('flux_cbh_schechter_entropy.csv').sort_values('z')
S_CBH_interp = interp1d(df['z'].values, df['S_CBH_density_kB_mpc3'].values, fill_value=(df['S_CBH_density_kB_mpc3'].values[0], 0.0), bounds_error=False)
S_CBH_z0 = S_CBH_interp(0.0)

def F_raw(z):
    return max(S_CBH_interp(z) / S_CBH_z0, 1e-90)

# =====================================================================
# 3. KINEMATICS & THE THROTTLE
# =====================================================================
def p_eff(z, p0, p_high, zt, w):
    return p_high + (p0 - p_high) / (1.0 + np.exp((z - zt) / w))

def E_flux(z, p0, p_high, zt, w, A):
    p_z = p_eff(z, p0, p_high, zt, w)
    # The New Amplitude Throttle Equation
    F_metric = 1.0 + A * (F_raw(z)**p_z - 1.0)
    flux_term = Oflux_0 * F_metric
    return np.sqrt(Or_flux*(1+z)**4 + Om_flux*(1+z)**3 + flux_term)

def sound_speed(z):
    return c / np.sqrt(3.0 * (1.0 + (R_b_factor / (1.0 + z))))

def integrand_rs(z, p0, p_high, zt, w, A):
    return sound_speed(z) / (100.0 * h_flux * E_flux(z, p0, p_high, zt, w, A))

def integrand_Dm(z, p0, p_high, zt, w, A):
    return c / (100.0 * h_flux * E_flux(z, p0, p_high, zt, w, A))

def eval_cmb(p0, p_high, zt, w, A):
    D_M, _ = quad(integrand_Dm, 0, z_star, args=(p0, p_high, zt, w, A), limit=200)
    r_s, _ = quad(integrand_rs, z_star, 1e5, args=(p0, p_high, zt, w, A), limit=200)
    return r_s / D_M

# SN Pre-compute LCDM Baseline
z_sn_eval = np.linspace(0.01, 2.3, 30)
def E_lcdm(z): return np.sqrt(Om_lcdm*(1+z)**3 + Ol_lcdm)
mu_lcdm = np.array([5.0 * np.log10((1+z) * (c/H0_lcdm) * quad(lambda x: 1/E_lcdm(x), 0, z)[0]) + 25.0 for z in z_sn_eval])

def eval_sn(p0, p_high, zt, w, A):
    mu_flux = np.array([5.0 * np.log10((1+z) * (c/H0_flux) * quad(integrand_Dm, 0, z, args=(p0, p_high, zt, w, A))[0]) + 25.0 for z in z_sn_eval])
    delta = mu_flux - mu_lcdm
    offset = np.mean(delta[(z_sn_eval >= 0.01) & (z_sn_eval <= 0.05)])
    shape_resid = np.abs(delta - offset)
    return np.max(shape_resid), shape_resid[-1], offset

# =====================================================================
# 4. FOCUSED 5D SCAN
# =====================================================================
print("Starting 5D Throttle Scan (Short-circuiting on CMB)...")
start_time = time.time()

# Narrow ranges based on the successful Scan 8
p0_vals = np.array([0.18, 0.185, 0.19])
phigh_vals = np.array([0.21, 0.22, 0.23])
zt_vals = np.array([1.5, 1.7, 1.9])
w_vals = np.array([0.3, 0.4])
A_vals = np.linspace(0.4, 0.9, 10) # Testing the Throttle from 40% to 90%

results = []
tested = 0
passed_cmb = 0

for p0 in p0_vals:
    for ph in phigh_vals:
        for zt in zt_vals:
            for w in w_vals:
                for A in A_vals:
                    tested += 1
                    theta = eval_cmb(p0, ph, zt, w, A)
                    err_ppm = abs(theta - theta_star_planck) / theta_star_planck * 1e6
                    
                    if err_ppm < 800:
                        passed_cmb += 1
                        sn_max, sn_end, sn_off = eval_sn(p0, ph, zt, w, A)
                        # We want models that drive the SN max residual down
                        score = err_ppm/100.0 + sn_max*100.0
                        results.append({
                            'A': A, 'p0': p0, 'p_high': ph, 'zt': zt, 'w': w,
                            'theta_err_ppm': err_ppm, 'sn_shape_max': sn_max, 'score': score
                        })

print(f"\nScan complete in {time.time()-start_time:.1f}s.")
print(f"Tested {tested} models. {passed_cmb} passed the CMB gate.")

df_res = pd.DataFrame(results)
if not df_res.empty:
    df_res = df_res.sort_values('score')
    print("\n=== TOP 5 THROTTLED MODELS ===")
    print(df_res.head(5).to_string(index=False))
    df_res.to_csv('scan9_throttle_best.csv', index=False)
else:
    print("No models passed.")