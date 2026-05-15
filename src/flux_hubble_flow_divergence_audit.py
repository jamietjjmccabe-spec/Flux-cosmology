
import pandas as pd
import numpy as np
import os
from pathlib import Path

def run_divergence_audit():
    df = pd.read_csv("host_structure_maturity_proxies_v3_100SNID_FINAL.csv")
    df = df.dropna(subset=['Reff_kpc', 'stellar_age_Gyr']).copy()
    
    # 1. Standard Maturity Calculation (D_hybrid)
    df['D_MR'] = (10**df['logMstar'] / df['Reff_kpc']) * (df['stellar_age_Gyr'] / 13.8)
    df['D_sigma'] = (df['sigma_v_km_s'] / 150)**2 * (df['stellar_age_Gyr'] / 13.8)
    df['D_hybrid'] = np.sqrt(df['D_MR'] * df['D_sigma'])
    
    median_D = df['D_hybrid'].median()
    df['D_reg'] = df['D_hybrid'] / median_D
    Dsat = 1.15
    df['F_fill'] = df['D_reg'] / Dsat
    
    def get_regime(f):
        if f < 0.5: return 'filling'
        if f <= 1.5: return 'saturated'
        return 'overflowing'
        
    df['regime'] = df['F_fill'].apply(get_regime)
    df['high_F'] = df['regime'] == 'overflowing'
    df['low_F']  = df['regime'] == 'filling'

    # Create output directory
    outdir = Path("v3.1_hubble_flow_divergence")
    outdir.mkdir(exist_ok=True)

    # Breakdown 1: Overall Calibrator vs Hubble-flow
    calibrators = df[df['is_calibrator'] == 1]
    hf = df[df['is_calibrator'] == 0]
    
    splits = []
    
    def analyze_subset(name, subset):
        if len(subset) == 0: return None
        n = len(subset)
        mean_H0 = subset['mean_H0_proxy'].mean()
        highF_subset = subset[subset['high_F']]
        lowF_subset = subset[subset['low_F']]
        
        hf_h0 = highF_subset['mean_H0_proxy'].mean() if len(highF_subset) > 0 else np.nan
        lf_h0 = lowF_subset['mean_H0_proxy'].mean() if len(lowF_subset) > 0 else np.nan
        delta = hf_h0 - lf_h0
        
        return {
            'subset': name,
            'N': n,
            'mean_H0': mean_H0,
            'highF_H0': hf_h0,
            'lowF_H0': lf_h0,
            'Adult_minus_Child_H0': delta
        }

    splits.append(analyze_subset('All Active', df))
    splits.append(analyze_subset('Calibrators Only', calibrators))
    splits.append(analyze_subset('Hubble-flow Only', hf))
    
    # Breakdown 2: Hubble-flow by Redshift (low, mid, high)
    hf_lowz = hf[hf['mean_z'] < 0.02]
    hf_midz = hf[(hf['mean_z'] >= 0.02) & (hf['mean_z'] < 0.05)]
    hf_highz = hf[hf['mean_z'] >= 0.05]
    
    splits.append(analyze_subset('HF: z < 0.02 (Local)', hf_lowz))
    splits.append(analyze_subset('HF: 0.02 <= z < 0.05 (Mid)', hf_midz))
    splits.append(analyze_subset('HF: z >= 0.05 (High)', hf_highz))
    
    # Breakdown 3: Hubble-flow by Outlier Status
    hf_core = hf[hf['mean_mu_resid'].abs() <= 0.2]
    hf_outlier = hf[hf['mean_mu_resid'].abs() > 0.2]
    
    splits.append(analyze_subset('HF: Core (|mu_resid| <= 0.2)', hf_core))
    splits.append(analyze_subset('HF: Outliers (|mu_resid| > 0.2)', hf_outlier))
    
    # Breakdown 4: Calibrators vs HF in Fixed-Mass Bins
    df['mass_bin'] = pd.cut(df['logMstar'], bins=[0, 9.5, 10.5, 20], labels=['Low (<9.5)', 'Mid (9.5-10.5)', 'High (>10.5)'])
    
    for mb in ['Low (<9.5)', 'Mid (9.5-10.5)', 'High (>10.5)']:
        splits.append(analyze_subset(f'Calibrators: Mass {mb}', calibrators[calibrators['mass_bin'] == mb]))
        splits.append(analyze_subset(f'HF: Mass {mb}', hf[hf['mass_bin'] == mb]))

    summary_df = pd.DataFrame([s for s in splits if s is not None])
    summary_df.to_csv(outdir / "hubble_flow_divergence_summary.csv", index=False)
    
    print("V3.1 Hubble-Flow Divergence Audit")
    print("=================================")
    print(summary_df.to_string(index=False))

if __name__ == '__main__':
    run_divergence_audit()
