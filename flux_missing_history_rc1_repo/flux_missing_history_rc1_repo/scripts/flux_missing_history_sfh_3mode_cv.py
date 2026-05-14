"""
Flux Cosmology: Missing History 3-Mode Constrained CV

Tests the Empirical SFH upgrade by running the constrained cross-validation 
across three distinct data subsets:
Mode A: All mixed (Empirical + Fallback)
Mode B: z0MGS-only (Empirical SFH only)
Mode C: Fallback-only (Delayed-Tau only)
"""

import argparse
import numpy as np
import pandas as pd
import json
from scipy.optimize import nnls

# ============================================================
# FROZEN FLUX BENCHMARK PARAMETERS
# ============================================================
K_CLOSURE = 18.0
GAMMA = 0.08
SIGMA_REF = 4.08e9
EPS = 1e-12
Q_VAL = 4.0  
K_FOLDS = 5  

def compute_rmse(y_true, y_pred):
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))

def fit_constrained_two_term(X_acc, X_sat, y_train):
    A_matrix = np.column_stack([X_acc, -X_sat])
    w, _ = nnls(A_matrix, y_train)
    return w[0], -w[1]

def run_cv_mode(df_subset, mode_name):
    N = len(df_subset)
    if N < K_FOLDS:
        return {"mode": mode_name, "N": N, "status": "Failed: Not enough samples for 5 folds"}
        
    df = df_subset.copy().reset_index(drop=True)
    
    indices = np.random.RandomState(42).permutation(N) 
    fold_sizes = np.full(K_FOLDS, N // K_FOLDS, dtype=int)
    fold_sizes[:N % K_FOLDS] += 1
    
    current = 0
    folds = []
    for fold_size in fold_sizes:
        start, stop = current, current + fold_size
        folds.append(indices[start:stop])
        current = stop

    results = {"baryon": [], "two_term": []}
    two_term_coefficients = []

    for i in range(K_FOLDS):
        test_idx = folds[i]
        train_idx = np.hstack([folds[j] for j in range(K_FOLDS) if j != i])
        
        train_df = df.iloc[train_idx].copy()
        test_df = df.iloc[test_idx].copy()
        
        y_train = train_df["Y_obs"].values
        y_test = test_df["Y_obs"].values
        
        # Hygiene: D_star is medians of training folds only
        D_star_train = np.nanmedian(train_df["D_hist_sfh"])
        
        for sub_df in [train_df, test_df]:
            sub_df["D_norm"] = np.log(1.0 + sub_df["D_hist_sfh"] / D_star_train)
            sub_df["G_bound"] = (sub_df["D_hist_sfh"]**Q_VAL) / (sub_df["D_hist_sfh"]**Q_VAL + D_star_train**Q_VAL)
            sub_df["X_acc"] = sub_df["D_norm"] * sub_df["G_geom"]
            sub_df["X_sat"] = sub_df["G_bound"] * sub_df["G_geom"]
        
        rmse_baryon = compute_rmse(y_test, np.zeros(len(y_test)))
        
        A_D, A_S = fit_constrained_two_term(train_df["X_acc"].values, train_df["X_sat"].values, y_train)
        y_pred_two = A_D * test_df["X_acc"].values + A_S * test_df["X_sat"].values
        rmse_two = compute_rmse(y_test, y_pred_two)
        
        results["baryon"].append(rmse_baryon)
        results["two_term"].append(rmse_two)
        two_term_coefficients.append((A_D, A_S))

    mean_baryon = np.mean(results["baryon"])
    mean_two = np.mean(results["two_term"])
    ad_folds = sum(1 for a_d, a_s in two_term_coefficients if a_d > 0)
    as_folds = sum(1 for a_d, a_s in two_term_coefficients if a_s < 0)

    return {
        "mode": mode_name,
        "N": N,
        "status": "Success",
        "rmse_baryon": mean_baryon,
        "rmse_two_term": mean_two,
        "ad_folds": ad_folds,
        "as_folds": as_folds
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, default="sparc_sfh_depth_catalog.csv")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input).dropna(
            subset=["V_obs", "V_b", "Rb", "Sigma_b", "D_hist_sfh", "D_hist_source"]
        ).copy()
    except FileNotFoundError:
        print(f"Error: {args.input} not found.")
        return

    # 1. Compute Base Constraints
    df["V_b_safe"] = np.maximum(df["V_b"], EPS)
    df["Y_obs"] = (df["V_obs"]**2 - df["V_b"]**2) / (df["V_b_safe"]**2)
    df = df[(df["Y_obs"] > -1) & (df["Y_obs"] < 100)].copy()
    
    df["Sigma_b_safe"] = np.maximum(df["Sigma_b"], EPS)
    df["L_chi"] = K_CLOSURE * df["Rb"] * (df["Sigma_b_safe"] / SIGMA_REF)**(-GAMMA)
    df["G_geom"] = df["Rb"] / (df["Rb"] + df["L_chi"])
    
    # 2. Segment by Data Source Mode
    mode_a = df
    mode_b = df[df["D_hist_source"] == "z0MGS_Empirical"]
    mode_c = df[df["D_hist_source"] == "Delayed_Tau_Fallback"]

    results = []
    results.append(run_cv_mode(mode_a, "All Mixed"))
    results.append(run_cv_mode(mode_b, "z0MGS-only (Empirical)"))
    results.append(run_cv_mode(mode_c, "Fallback-only"))

    # 3. Print Defensive Table
    print("\n" + "="*88)
    print("MISSING HISTORY: 3-MODE EMPIRICAL SFH CROSS-VALIDATION")
    print("="*88)
    print(f"{'Sample Mode':<25} | {'N':<4} | {'RMSE baryon':<12} | {'RMSE two-term':<14} | {'(A_D>0) folds':<13} | {'(A_S<0) folds'}")
    print("-" * 88)
    
    for r in results:
        if r["status"] == "Success":
            print(f"{r['mode']:<25} | {r['N']:<4} | {r['rmse_baryon']:<12.4f} | {r['rmse_two_term']:<14.4f} | {r['ad_folds']}/5           | {r['as_folds']}/5")
        else:
            print(f"{r['mode']:<25} | {r['N']:<4} | {r['status']}")
            
    print("="*88 + "\n")

    with open("flux_missing_history_3mode_cv_summary.json", "w") as f:
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    main()