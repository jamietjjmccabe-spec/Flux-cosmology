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

# Flux-CBH Base (from successful Matter-Float fit)
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

# =====================================================================
# 2. LOAD 4-TERM FLUX DATA
# =====================================================================
print("Loading Schechter Entropy data...")
df = pd.read_csv('flux_cbh_schechter_entropy.csv').sort_values('z')
source_col = "Omega_flux_4term_fit"

if source_col not in df.columns:
    print(f"Warning: {source_col} not found. Falling back to S_CBH_density_kB_mpc3")
    source_col = "S_CBH_density_kB_mpc3"

F_interp = interp1d(
    df["z"].values,
    df[source_col].values,
    fill_value=(df[source_col].values[0], df[source_col].values[-1]),
    bounds_error=False
)
F0 = F_interp(0.0)

def F_raw(z):
    return max(F_interp(z) / F0, 1e-90)

# =====================================================================
# 3. PRE-COMPUTE LCDM SN BASELINES (Floating Om_ref)
# =====================================================================
print("Pre-computing LCDM baselines for SN shape fitting...")
z_sn_eval = np.linspace(0.01, 2.3, 30)
anchor_mask = (z_sn_eval >= 0.01) & (z_sn_eval <= 0.05)
lcdm_baselines = {}
Om_refs = np.linspace(0.24, 0.34, 11)

for Om_ref in Om_refs:
    Ol_ref = 1.0 - Om_ref
    # We use H0_flux here to align the baseline conceptually;
    # the intercept offset will be strictly marginalized away anyway.
    def E_lcdm(z, Om=Om_ref, Ol=Ol_ref): return np.sqrt(Om*(1+z)**3 + Ol)
    mu = np.array([5.0 * np.log10((1+z) * (c/H0_flux) * quad(lambda x: 1/E_lcdm(x), 0, z)[0]) + 25.0 for z in z_sn_eval])
    lcdm_baselines[Om_ref] = mu

# =====================================================================
# 4. KINEMATICS & THE PIVOT
# =====================================================================
# Fixed transition for the compression exponent (from Scan 8)
zt_fixed = 1.70
w_fixed = 0.35

# Width for the Pivot step function
wp_fixed = 0.30

def p_eff(z, p0, p_high):
    return p_high + (p0 - p_high) / (1.0 + np.exp((z - zt_fixed) / w_fixed))

def G_z(z, z_p):
    """Smooth step function turning ON at z_p"""
    return 1.0 / (1.0 + np.exp(-(z - z_p) / wp_fixed))

def E_flux(z, p0, p_high, beta, z_p):
    p_z = p_eff(z, p0, p_high)
    
    # Normalize Pivot so F_metric(0) = 1 exactly
    G_0 = G_z(0, z_p)
    pivot_factor = (1.0 + beta * G_z(z, z_p)) / (1.0 + beta * G_0)
    
    F_metric = (F_raw(z)**p_z) * pivot_factor
    flux_term = Oflux_0 * F_metric
    return np.sqrt(Or_flux*(1+z)**4 + Om_flux*(1+z)**3 + flux_term)

def sound_speed(z):
    return c / np.sqrt(3.0 * (1.0 + (R_b_factor / (1.0 + z))))

def integrand_rs(z, p0, p_high, beta, z_p):
    return sound_speed(z) / (100.0 * h_flux * E_flux(z, p0, p_high, beta, z_p))

def integrand_Dm(z, p0, p_high, beta, z_p):
    return c / (100.0 * h_flux * E_flux(z, p0, p_high, beta, z_p))

def eval_cmb(p0, p_high, beta, z_p):
    D_M, _ = quad(integrand_Dm, 0, z_star, args=(p0, p_high, beta, z_p), limit=200)
    r_s, _ = quad(integrand_rs, z_star, 1e5, args=(p0, p_high, beta, z_p), limit=200)
    return r_s / D_M

def eval_sn(p0, p_high, beta, z_p):
    mu_flux = np.array([5.0 * np.log10((1+z) * (c/H0_flux) * quad(integrand_Dm, 0, z, args=(p0, p_high, beta, z_p))[0]) + 25.0 for z in z_sn_eval])
    
    best_shape_max = 999.0
    best_Om_ref = None
    best_shape_resid = None
    best_offset = None
    
    # Test against all LCDM baselines and keep the tightest shape match
    for Om_ref, mu_lcdm in lcdm_baselines.items():
        delta = mu_flux - mu_lcdm
        offset = np.mean(delta[anchor_mask])
        shape_resid = delta - offset
        shape_max = np.max(np.abs(shape_resid))
        
        if shape_max < best_shape_max:
            best_shape_max = shape_max
            best_Om_ref = Om_ref
            best_shape_resid = shape_resid
            best_offset = offset
            
    return best_shape_max, best_Om_ref, best_offset, best_shape_resid

# =====================================================================
# 5. FOCUSED 4D SCAN
# =====================================================================
print("\nStarting Scan 10: Redshift Pivot...")
start_time = time.time()

p0_vals = np.array([0.17, 0.185, 0.20])
phigh_vals = np.array([0.21, 0.22, 0.23])
beta_vals = np.linspace(0.0, 0.6, 7)
zp_vals = np.linspace(0.5, 2.0, 4)

results = []
tested = 0
passed_cmb_gate = 0

for p0 in p0_vals:
    for ph in phigh_vals:
        for beta in beta_vals:
            for zp in zp_vals:
                tested += 1
                theta = eval_cmb(p0, ph, beta, zp)
                err_ppm = abs(theta - theta_star_planck) / theta_star_planck * 1e6
                
                # CMB Filter: Only proceed if error is < 1500 ppm
                if err_ppm < 1500:
                    passed_cmb_gate += 1
                    sn_max, best_Om_ref, sn_off, _ = eval_sn(p0, ph, beta, zp)
                    
                    # Score penalizes both CMB error and SN residual
                    score = (err_ppm/500.0)**2 + (sn_max*100.0)**2
                    results.append({
                        'beta': beta, 'z_p': zp, 'p0': p0, 'p_high': ph,
                        'theta_err_ppm': err_ppm, 'sn_shape_max': sn_max, 
                        'best_Om_ref': best_Om_ref, 'score': score
                    })

print(f"\nScan complete in {time.time()-start_time:.1f}s.")
print(f"Tested {tested} models. {passed_cmb_gate} passed the loose CMB gate (< 1500 ppm).")

df_res = pd.DataFrame(results)
if not df_res.empty:
    df_res = df_res.sort_values('score')
    print("\n=== TOP 5 PIVOT MODELS ===")
    print(df_res.head(5).to_string(index=False))
    df_res.to_csv('scan10_pivot_results.csv', index=False)
    
    # Plot the best fit
    best = df_res.iloc[0]
    _, best_Om_ref, _, best_shape_resid = eval_sn(best['p0'], best['p_high'], best['beta'], best['z_p'])
    
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axhline(0, color='white', linestyle='--', alpha=0.5)
    ax.plot(z_sn_eval, best_shape_resid, color='magenta', linewidth=3, 
             label=f"beta={best['beta']:.2f}, zp={best['z_p']:.2f}, p0={best['p0']:.3f}")
    
    ax.axhspan(-0.1, 0.1, color='grey', alpha=0.2, label=r'SN Target Zone ($\pm 0.10$ mag)')
    ax.axhspan(-0.05, 0.05, color='green', alpha=0.2, label=r'Excellent ($\pm 0.05$ mag)')
    
    ax.set_title(f"Best Pivot Model vs LCDM (Om={best_Om_ref:.2f}) | Max |Δμ| = {best['sn_shape_max']:.3f}", fontsize=14, pad=15)
    ax.set_xlabel("Redshift (z)")
    ax.set_ylabel(r"$\Delta\mu_{shape}$ [mag]")
    ax.set_xlim(0, 2.3)
    ax.legend(loc='lower left')
    ax.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig('scan10_best_pivot_sn.png', dpi=300)
    print("\nSaved plot as scan10_best_pivot_sn.png")
else:
    print("No models passed.")