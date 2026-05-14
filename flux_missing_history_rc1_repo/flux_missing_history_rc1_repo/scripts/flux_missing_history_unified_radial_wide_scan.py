"""
Flux Cosmology: Unified Radial Bounded Map - Wide Grid Scan

Performs a wide-grid search to calibrate the sigmoid bounded map linking
empirical Missing-History depth to radial residual amplitude.

This corrected version uses the same SPARC radial parser/baryon convention as
flux_missing_history_unified_radial.py:

    V_b^2 = sign(V_gas) * V_gas^2 + V_disk^2 + V_bulge^2

This is necessary because SPARC gas contributions can be signed.
"""

import os
import json
import numpy as np
import pandas as pd

# ============================================================
# FROZEN FLUX BENCHMARK PARAMETERS
# ============================================================
K_CLOSURE = 18.0
GAMMA = 0.08
SIGMA_REF = 4.08e9
Q_VAL = 4.0
EPS = 1e-12

# Use the same values as the successful unified radial benchmark.
# Replace later with empirical fold-mean coefficients if you export them.
A_D_FIXED = 0.50
A_S_FIXED = -9.50

# ============================================================
# WIDE GRID SEARCH PARAMETERS
# ============================================================
Y_MAX_GRID = np.linspace(1.0, 12.0, 12)       # Max residual amplitude ceiling
SLOPE_GRID = np.linspace(0.1, 8.0, 80)        # Sigmoid transition slope s
SHIFT_GRID = np.linspace(-14.0, 2.0, 33)      # Sigmoid center h_0


def sigmoid_map(H_raw, Y_max, slope, shift):
    """Physically bounded activation map to prevent negative radial residuals."""
    x = np.clip(-slope * (H_raw - shift), -500, 500)
    return Y_max / (1.0 + np.exp(x))


def compute_rmse(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def signed_square(v):
    """SPARC convention: gas velocity contribution may be signed."""
    return np.sign(v) * (v ** 2)


def load_rotmod(radial_file):
    """
    Load a SPARC *_rotmod.dat file.

    Expected first columns:
        Rad, Vobs, errV, Vgas, Vdisk, Vbul

    Uses sep=r"\s+" instead of deprecated delim_whitespace.
    """
    rad_df = pd.read_csv(
        radial_file,
        sep=r"\s+",
        comment="#",
        names=["Rad", "Vobs", "errV", "Vgas", "Vdisk", "Vbul"],
        usecols=[0, 1, 2, 3, 4, 5],
        engine="python",
    )

    for col in ["Rad", "Vobs", "errV", "Vgas", "Vdisk", "Vbul"]:
        rad_df[col] = pd.to_numeric(rad_df[col], errors="coerce")

    rad_df = rad_df.dropna(subset=["Rad", "Vobs", "Vgas", "Vdisk", "Vbul"]).copy()

    if rad_df.empty:
        raise ValueError(f"No valid radial rows found in {radial_file}")

    return rad_df


def baryon_vsq(rad_df):
    Vgas = rad_df["Vgas"].values
    Vdisk = rad_df["Vdisk"].values
    Vbul = rad_df["Vbul"].values

    Vb_sq = signed_square(Vgas) + Vdisk**2 + Vbul**2
    Vb_sq = np.maximum(Vb_sq, EPS)
    return Vb_sq


def main():
    catalog_path = "sparc_sfh_depth_catalog.csv"
    data_dir = "Rotmod_LTG-main"
    output_json = "flux_unified_radial_wide_scan_summary.json"
    output_csv = "flux_unified_radial_wide_scan_per_galaxy.csv"

    targets = ["DDO154", "NGC2403", "NGC3198", "UGC00128", "UGC02885"]

    if not os.path.exists(catalog_path):
        print(f"Error: Could not find {catalog_path}")
        return

    df_cat = pd.read_csv(catalog_path)

    required_catalog_cols = ["Galaxy", "Rb", "Sigma_b", "D_hist_sfh"]
    missing = [c for c in required_catalog_cols if c not in df_cat.columns]
    if missing:
        raise ValueError(f"Catalog missing required columns: {missing}")

    D_star = float(np.nanmedian(df_cat["D_hist_sfh"]))
    target_data = {}

    for gal in targets:
        if "match_name" in df_cat.columns:
            row_df = df_cat[df_cat["match_name"].astype(str).str.upper() == gal.upper()]
            if row_df.empty:
                row_df = df_cat[df_cat["Galaxy"].astype(str).str.upper() == gal.upper()]
        else:
            row_df = df_cat[df_cat["Galaxy"].astype(str).str.upper() == gal.upper()]

        if row_df.empty:
            print(f"Warning: {gal} not found in catalog.")
            continue

        row = row_df.iloc[0]

        Rb = float(row["Rb"])
        Sigma_b = max(float(row["Sigma_b"]), EPS)
        D_hist = max(float(row["D_hist_sfh"]), EPS)

        L_chi = K_CLOSURE * Rb * (Sigma_b / SIGMA_REF) ** (-GAMMA)

        D_norm = np.log1p(D_hist / max(D_star, EPS))
        G_bound = (D_hist**Q_VAL) / (D_hist**Q_VAL + D_star**Q_VAL)
        H_raw = A_D_FIXED * D_norm + A_S_FIXED * G_bound

        radial_file = os.path.join(data_dir, f"{gal}_rotmod.dat")
        if not os.path.exists(radial_file):
            print(f"Warning: {radial_file} not found.")
            continue

        rad_df = load_rotmod(radial_file)
        Vb_sq = baryon_vsq(rad_df)

        target_data[gal] = {
            "L_chi": float(L_chi),
            "H_raw": float(H_raw),
            "D_norm": float(D_norm),
            "G_bound": float(G_bound),
            "Rad": rad_df["Rad"].values,
            "Vobs": rad_df["Vobs"].values,
            "Vb_sq": Vb_sq,
        }

    if not target_data:
        print("No target data loaded. Exiting.")
        return

    all_vobs = []
    all_vb = []
    for _, data in target_data.items():
        all_vobs.extend(data["Vobs"])
        all_vb.extend(np.sqrt(data["Vb_sq"]))

    baryon_rmse = compute_rmse(np.array(all_vobs), np.array(all_vb))

    print("\n" + "=" * 70)
    print("RUNNING WIDE-GRID BOUNDED MAP SCAN")
    print(
        f"Grid Size: {len(Y_MAX_GRID)} x {len(SLOPE_GRID)} x {len(SHIFT_GRID)} "
        f"= {len(Y_MAX_GRID) * len(SLOPE_GRID) * len(SHIFT_GRID)} permutations"
    )
    print("=" * 70)

    print("\nSanity check: per-galaxy baryon RMSE")
    for gal, data in target_data.items():
        rb = compute_rmse(data["Vobs"], np.sqrt(data["Vb_sq"]))
        print(f"  {gal:<10} {rb:.3f} km/s | H_raw={data['H_raw']:.3f}")
    print(f"  {'GLOBAL':<10} {baryon_rmse:.3f} km/s")

    best_global_rmse = float("inf")
    best_params = None
    best_per_galaxy = None

    for y_max in Y_MAX_GRID:
        for slope in SLOPE_GRID:
            for shift in SHIFT_GRID:
                all_obs = []
                all_pred = []
                per_gal = []

                for gal, data in target_data.items():
                    Y_amp = sigmoid_map(data["H_raw"], y_max, slope, shift)
                    G_radial = data["Rad"] / (data["Rad"] + data["L_chi"])

                    V_flux_sq = data["Vb_sq"] * (1.0 + Y_amp * G_radial)
                    V_flux_sq = np.maximum(V_flux_sq, EPS)
                    V_flux = np.sqrt(V_flux_sq)

                    all_obs.extend(data["Vobs"])
                    all_pred.extend(V_flux)

                    rb = compute_rmse(data["Vobs"], np.sqrt(data["Vb_sq"]))
                    rf = compute_rmse(data["Vobs"], V_flux)

                    per_gal.append({
                        "Galaxy": gal,
                        "RMSE_baryon": rb,
                        "RMSE_flux": rf,
                        "beats_baryon": bool(rf < rb),
                        "Y_amp": float(Y_amp),
                        "H_raw": float(data["H_raw"]),
                        "D_norm": float(data["D_norm"]),
                        "G_bound": float(data["G_bound"]),
                        "L_chi": float(data["L_chi"]),
                    })

                current_rmse = compute_rmse(np.array(all_obs), np.array(all_pred))

                is_edge = (
                    y_max == Y_MAX_GRID[-1]
                    or y_max == Y_MAX_GRID[0]
                    or slope == SLOPE_GRID[-1]
                    or slope == SLOPE_GRID[0]
                    or shift == SHIFT_GRID[-1]
                    or shift == SHIFT_GRID[0]
                )

                if current_rmse < best_global_rmse:
                    best_global_rmse = current_rmse
                    best_params = {
                        "Y_max": float(y_max),
                        "slope": float(slope),
                        "shift": float(shift),
                        "RMSE_flux": float(current_rmse),
                        "RMSE_baryon": float(baryon_rmse),
                        "improvement": float(baryon_rmse - current_rmse),
                        "is_edge_solution": bool(is_edge),
                    }
                    best_per_galaxy = per_gal

    print("\n" + "=" * 70)
    print("WIDE GRID SCAN RESULTS")
    print("=" * 70)
    print(f"Baryon-only RMSE : {baryon_rmse:.4f} km/s")
    print(f"Best Flux RMSE   : {best_params['RMSE_flux']:.4f} km/s")
    print(f"Improvement      : {best_params['improvement']:.4f} km/s")
    print("-" * 70)
    print("Best Bounded Map Parameters:")
    print(f"  Y_max (Ceiling) : {best_params['Y_max']:.2f}")
    print(f"  slope (s)       : {best_params['slope']:.2f}")
    print(f"  shift (h_0)     : {best_params['shift']:.2f}")
    print("-" * 70)

    if best_params["is_edge_solution"]:
        print("WARNING: Optimum is hitting a boundary of the search grid.")
        print("Regularization or wider bounds may still be required.")
    else:
        print("SUCCESS: Optimum found inside the parameter space. The map is stable.")

    print("-" * 70)
    for row in best_per_galaxy:
        verdict = "beats baryon" if row["beats_baryon"] else "worse"
        print(
            f"{row['Galaxy']:<10} "
            f"RMSE baryon={row['RMSE_baryon']:.3f} "
            f"RMSE flux={row['RMSE_flux']:.3f} "
            f"Y_amp={row['Y_amp']:.3f} "
            f"H={row['H_raw']:.3f} "
            f"{verdict}"
        )

    print("=" * 70 + "\n")

    with open(output_json, "w", encoding="utf-8") as f:
        json.dump({"best_params": best_params, "per_galaxy": best_per_galaxy}, f, indent=4)

    pd.DataFrame(best_per_galaxy).to_csv(output_csv, index=False)

    print(f"Wrote: {output_json}")
    print(f"Wrote: {output_csv}")


if __name__ == "__main__":
    main()
