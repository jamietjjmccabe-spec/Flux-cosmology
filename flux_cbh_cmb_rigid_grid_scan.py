#!/usr/bin/env python3
"""
flux_cbh_cmb_rigid_grid_scan.py

CMB-rigid Flux-CBH grid scan.

Purpose
-------
Test whether a compressed CBH-entropy dark-sector history can raise the
acoustic-angle inferred H0 while keeping early-universe physical densities fixed.

This is the referee-safe scan:

    omega_b = Omega_b h^2 fixed
    omega_m = Omega_m h^2 fixed
    omega_r = Omega_r h^2 fixed

Then scan only:

    H0
    p

where:

    F_metric(z) = F_flux(z)^p

and:

    E^2(z) =
        Omega_r (1+z)^4
      + Omega_m (1+z)^3
      + Omega_flux,0 F_metric(z)

with flat closure:

    Omega_flux,0 = 1 - Omega_m - Omega_r

Input
-----
    flux_cbh_schechter_entropy.csv

Required columns:
    z
    Omega_flux_4term_fit

Outputs
-------
    flux_cbh_cmb_rigid_grid_results.csv
    flux_cbh_cmb_rigid_best_models.csv
    flux_cbh_cmb_rigid_theta_heatmap.png
    flux_cbh_cmb_rigid_h0_vs_p.png
    flux_cbh_cmb_rigid_best_history.png
    flux_cbh_cmb_rigid_q_history.png

Important
---------
This is a lightweight acoustic-angle diagnostic, not a CLASS/CAMB replacement.
It is designed to test whether the Flux-CBH history moves the late-time
distance integral in the right direction while leaving r_s nearly intact.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d


# ============================================================
# 1. Constants and scan controls
# ============================================================

C_KM_S = 299792.458

# Planck-like acoustic angle target
THETA_STAR_TARGET = 0.0104132

# Loose diagnostic tolerance.
# Planck precision is far tighter, but this toy integral is not CLASS/CAMB.
THETA_TOL = 2.5e-5

Z_STAR = 1089.92

# Physical densities fixed: CMB-rigid test
OMEGA_B_H2 = 0.02237
OMEGA_M_H2 = 0.1430
OMEGA_R_H2 = 4.18e-5

# Grid ranges
H0_MIN = 65.0
H0_MAX = 80.0
N_H0 = 301

P_MIN = 0.05
P_MAX = 0.50
N_P = 181

# For full distance integral
N_Z_DIST = 9000
N_Z_SOUND = 12000

# Local comparison target, not used as a fit prior
LOCAL_H0_TARGET = 73.0


# ============================================================
# 2. Flux history loading
# ============================================================

def normalize_flux_history(z: np.ndarray, F_raw: np.ndarray, p: float) -> np.ndarray:
    """
    Build compressed metric response:

        F_metric(z) = [F_raw(z)/F_raw(0)]^p

    then renormalize so F_metric(0)=1.

    p=1   -> raw response
    p<1   -> compressed/screened response
    """
    z = np.asarray(z, dtype=float)
    F = np.asarray(F_raw, dtype=float)

    idx0 = int(np.argmin(np.abs(z - 0.0)))
    F0 = F[idx0]

    if not np.isfinite(F0) or F0 <= 0:
        raise ValueError(f"Invalid F(0) in Flux history: {F0}")

    F = F / F0
    F = np.clip(F, 0.0, None)
    F = F ** p

    F0p = F[idx0]
    if not np.isfinite(F0p) or F0p <= 0:
        raise ValueError(f"Invalid compressed F(0): {F0p}")

    return F / F0p


def load_flux_csv(path: str | Path = "flux_cbh_schechter_entropy.csv") -> pd.DataFrame:
    path = Path(path)

    if not path.exists():
        print("\nERROR: Required CSV not found.")
        print(f"Looking for: {path.resolve()}")
        print(f"Current working folder: {Path.cwd()}")
        print("\nFix: place flux_cbh_schechter_entropy.csv beside this script,")
        print("or edit CSV_PATH in main().")
        raise FileNotFoundError(path)

    df = pd.read_csv(path)

    required = {"z", "Omega_flux_4term_fit"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = df.sort_values("z").reset_index(drop=True)

    # Clean pathological values
    df["Omega_flux_4term_fit"] = pd.to_numeric(df["Omega_flux_4term_fit"], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["z", "Omega_flux_4term_fit"])

    if len(df) < 10:
        raise ValueError("Flux CSV has too few valid rows after cleaning.")

    return df


def make_flux_interp(df: pd.DataFrame, p: float):
    z = df["z"].values
    F_raw = df["Omega_flux_4term_fit"].values
    F_metric = normalize_flux_history(z, F_raw, p)

    interp = interp1d(
        z,
        F_metric,
        kind="linear",
        bounds_error=False,
        # z is ascending: fill below min-z with F[0], above max-z with F[-1].
        fill_value=(F_metric[0], F_metric[-1]),
    )

    return interp, F_metric


# ============================================================
# 3. Background cosmology
# ============================================================

def params_from_H0(H0: float) -> dict | None:
    h = H0 / 100.0

    Om_b = OMEGA_B_H2 / h**2
    Om_m = OMEGA_M_H2 / h**2
    Om_r = OMEGA_R_H2 / h**2

    Om_flux = 1.0 - Om_m - Om_r

    if Om_flux <= 0:
        return None

    return {
        "H0": H0,
        "h": h,
        "Omega_b": Om_b,
        "Omega_m": Om_m,
        "Omega_r": Om_r,
        "Omega_flux0": Om_flux,
    }


def E_lcdm(z: np.ndarray, par: dict) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    return np.sqrt(
        par["Omega_r"] * (1.0 + z)**4
        + par["Omega_m"] * (1.0 + z)**3
        + par["Omega_flux0"]
    )


def E_flux(z: np.ndarray, par: dict, flux_interp) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    F = flux_interp(z)
    F = np.clip(F, 0.0, None)

    return np.sqrt(
        par["Omega_r"] * (1.0 + z)**4
        + par["Omega_m"] * (1.0 + z)**3
        + par["Omega_flux0"] * F
    )


def omega_flux_z(z: np.ndarray, par: dict, flux_interp) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    F = np.clip(flux_interp(z), 0.0, None)
    E2 = (
        par["Omega_r"] * (1.0 + z)**4
        + par["Omega_m"] * (1.0 + z)**3
        + par["Omega_flux0"] * F
    )
    return par["Omega_flux0"] * F / E2


def omega_m_z(z: np.ndarray, par: dict, flux_interp=None) -> np.ndarray:
    z = np.asarray(z, dtype=float)

    if flux_interp is None:
        E2 = (
            par["Omega_r"] * (1.0 + z)**4
            + par["Omega_m"] * (1.0 + z)**3
            + par["Omega_flux0"]
        )
    else:
        F = np.clip(flux_interp(z), 0.0, None)
        E2 = (
            par["Omega_r"] * (1.0 + z)**4
            + par["Omega_m"] * (1.0 + z)**3
            + par["Omega_flux0"] * F
        )

    return par["Omega_m"] * (1.0 + z)**3 / E2


def q_flux_z(z: np.ndarray, par: dict, flux_interp) -> np.ndarray:
    """
    Approximate q(z).

    For a dynamic dark sector this exact q requires dlnH/dlna. Here we compute it
    numerically from E(z):

        q(z) = -1 + (1+z)/E(z) * dE/dz

    This is more honest than assuming w=-1 for the Flux term.
    """
    z = np.asarray(z, dtype=float)
    E = E_flux(z, par, flux_interp)
    dE_dz = np.gradient(E, z)
    return -1.0 + (1.0 + z) * dE_dz / E


def find_acceleration_onset_z(par: dict, flux_interp, zmax: float = 5.0) -> float:
    z_grid = np.linspace(0.0, zmax, 2000)
    q = q_flux_z(z_grid, par, flux_interp)

    # q < 0 today, q > 0 in matter era. Find first crossing from low z upward.
    sign = np.sign(q)
    crossings = np.where(np.diff(sign) != 0)[0]

    if len(crossings) == 0:
        return np.nan

    i = crossings[0]
    z0, z1 = z_grid[i], z_grid[i + 1]
    q0, q1 = q[i], q[i + 1]

    if q1 == q0:
        return z0

    return z0 - q0 * (z1 - z0) / (q1 - q0)


# ============================================================
# 4. CMB acoustic quantities
# ============================================================

def sound_speed_over_c(z: np.ndarray) -> np.ndarray:
    """
    Photon-baryon sound speed / c.

        c_s/c = 1 / sqrt[3(1 + R)]

    Approx:
        R(z) = 31500 * omega_b / (1+z)
    """
    z = np.asarray(z, dtype=float)
    R = 31500.0 * OMEGA_B_H2 / (1.0 + z)
    return 1.0 / np.sqrt(3.0 * (1.0 + R))


def compute_theta_star(H0: float, flux_interp) -> tuple[float, float, float]:
    par = params_from_H0(H0)
    if par is None:
        return np.nan, np.nan, np.nan

    # Distance to last scattering: 0 -> z*
    z_dist = np.linspace(0.0, Z_STAR, N_Z_DIST)
    E_dist = E_flux(z_dist, par, flux_interp)
    D_M = (C_KM_S / H0) * trapezoid(1.0 / E_dist, z_dist)

    # Sound horizon: z* -> high z, log grid
    z_sound = np.logspace(np.log10(Z_STAR), 6.0, N_Z_SOUND)
    E_sound = E_flux(z_sound, par, flux_interp)
    cs = sound_speed_over_c(z_sound)
    r_s = (C_KM_S / H0) * trapezoid(cs / E_sound, z_sound)

    theta = r_s / D_M
    return theta, r_s, D_M


def compute_theta_star_lcdm(H0: float) -> tuple[float, float, float]:
    par = params_from_H0(H0)
    if par is None:
        return np.nan, np.nan, np.nan

    z_dist = np.linspace(0.0, Z_STAR, N_Z_DIST)
    E_dist = E_lcdm(z_dist, par)
    D_M = (C_KM_S / H0) * trapezoid(1.0 / E_dist, z_dist)

    z_sound = np.logspace(np.log10(Z_STAR), 6.0, N_Z_SOUND)
    E_sound = E_lcdm(z_sound, par)
    cs = sound_speed_over_c(z_sound)
    r_s = (C_KM_S / H0) * trapezoid(cs / E_sound, z_sound)

    theta = r_s / D_M
    return theta, r_s, D_M


# ============================================================
# 5. Grid scan
# ============================================================

def run_grid_scan(df_flux: pd.DataFrame) -> pd.DataFrame:
    H0_grid = np.linspace(H0_MIN, H0_MAX, N_H0)
    p_grid = np.linspace(P_MIN, P_MAX, N_P)

    rows = []

    print("\nRunning CMB-rigid Flux-CBH grid scan")
    print("=" * 78)
    print(f"H0 grid: {H0_grid[0]:.2f} -> {H0_grid[-1]:.2f} with {len(H0_grid)} points")
    print(f"p grid:  {p_grid[0]:.3f} -> {p_grid[-1]:.3f} with {len(p_grid)} points")
    print(f"Total models: {len(H0_grid) * len(p_grid):,}")
    print()

    # Precompute flux interpolation for each p
    for j, p in enumerate(p_grid):
        if j % max(1, len(p_grid)//10) == 0:
            print(f"  p step {j+1}/{len(p_grid)}: p={p:.4f}")

        flux_interp, _ = make_flux_interp(df_flux, p)

        for H0 in H0_grid:
            par = params_from_H0(H0)
            if par is None:
                continue

            theta, r_s, D_M = compute_theta_star(H0, flux_interp)

            if not np.isfinite(theta):
                continue

            theta_err = theta - THETA_STAR_TARGET
            theta_abs_err = abs(theta_err)
            theta_ppm = theta_err / THETA_STAR_TARGET * 1e6

            q0 = q_flux_z(np.array([0.0, 0.001, 0.002]), par, flux_interp)[0]
            z_acc = find_acceleration_onset_z(par, flux_interp)

            rows.append({
                "p": p,
                "H0": H0,
                "theta_star": theta,
                "theta_err": theta_err,
                "theta_abs_err": theta_abs_err,
                "theta_err_ppm": theta_ppm,
                "r_s_Mpc": r_s,
                "D_M_Mpc": D_M,
                "Omega_m0": par["Omega_m"],
                "Omega_b0": par["Omega_b"],
                "Omega_r0": par["Omega_r"],
                "Omega_flux0": par["Omega_flux0"],
                "q0": q0,
                "z_acc": z_acc,
                "accepted_theta": theta_abs_err <= THETA_TOL,
                "H0_distance_to_73": abs(H0 - LOCAL_H0_TARGET),
            })

    out = pd.DataFrame(rows)

    if out.empty:
        raise RuntimeError("Grid scan produced no valid rows.")

    # Combined lightweight score:
    #   prioritize theta agreement, then closeness to local H0=73, then reasonable z_acc.
    # This is not a likelihood; it is a sorting diagnostic.
    out["score"] = (
        (out["theta_abs_err"] / THETA_TOL) ** 2
        + ((out["H0"] - LOCAL_H0_TARGET) / 2.0) ** 2
        + ((out["z_acc"] - 0.63) / 0.35) ** 2
    )

    return out.sort_values("score").reset_index(drop=True)


def lcdm_reference_scan() -> pd.DataFrame:
    H0_grid = np.linspace(H0_MIN, H0_MAX, N_H0)
    rows = []

    for H0 in H0_grid:
        par = params_from_H0(H0)
        if par is None:
            continue

        theta, r_s, D_M = compute_theta_star_lcdm(H0)
        rows.append({
            "H0": H0,
            "theta_star": theta,
            "theta_err": theta - THETA_STAR_TARGET,
            "theta_abs_err": abs(theta - THETA_STAR_TARGET),
            "r_s_Mpc": r_s,
            "D_M_Mpc": D_M,
            "Omega_m0": par["Omega_m"],
            "Omega_flux0": par["Omega_flux0"],
        })

    return pd.DataFrame(rows).sort_values("theta_abs_err").reset_index(drop=True)


# ============================================================
# 6. Plotting
# ============================================================

def make_heatmap(grid: pd.DataFrame, value_col: str, title: str, output: str):
    piv = grid.pivot_table(index="p", columns="H0", values=value_col, aggfunc="mean")

    plt.figure(figsize=(12, 7))
    plt.imshow(
        piv.values,
        origin="lower",
        aspect="auto",
        extent=[piv.columns.min(), piv.columns.max(), piv.index.min(), piv.index.max()],
    )
    plt.colorbar(label=value_col)
    plt.xlabel(r"$H_0$ [km/s/Mpc]")
    plt.ylabel(r"Flux compression power $p$")
    plt.title(title)

    # Target H0 and expected viable p window
    plt.axvline(73.0, color="white", linestyle=":", alpha=0.8, label="H0=73")
    plt.axhspan(0.15, 0.25, color="white", alpha=0.08, label="p≈0.15–0.25")

    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(output, dpi=180)
    plt.close()


def make_plots(df_flux: pd.DataFrame, grid: pd.DataFrame, best: pd.DataFrame, lcdm_ref: pd.DataFrame):
    plt.style.use("dark_background")

    # Heatmap of theta error ppm
    make_heatmap(
        grid,
        "theta_err_ppm",
        r"Flux-CBH CMB-Rigid Grid: $\theta_\star$ Error [ppm]",
        "flux_cbh_cmb_rigid_theta_heatmap.png",
    )

    # Best H0 for each p
    by_p = (
        grid.sort_values("theta_abs_err")
        .groupby("p")
        .first()
        .reset_index()
        .sort_values("p")
    )

    plt.figure(figsize=(10, 6))
    plt.plot(by_p["p"], by_p["H0"], color="lime", marker="o", markersize=3, linewidth=2,
             label=r"Best $\theta_\star$ H0 at each p")
    plt.axhline(73.0, color="white", linestyle=":", label="H0=73")
    plt.axhline(lcdm_ref.iloc[0]["H0"], color="cyan", linestyle="--",
                label=rf"ΛCDM ref H0={lcdm_ref.iloc[0]['H0']:.2f}")
    plt.axvspan(0.15, 0.25, color="white", alpha=0.08, label="p≈0.15–0.25")
    plt.xlabel(r"Flux compression power $p$")
    plt.ylabel(r"Best acoustic-angle $H_0$ [km/s/Mpc]")
    plt.title("CMB-Rigid Flux-CBH: H0 Shift vs Compression Power")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig("flux_cbh_cmb_rigid_h0_vs_p.png", dpi=180)
    plt.close()

    # Best history
    best_row = best.iloc[0]
    p_best = float(best_row["p"])
    H0_best = float(best_row["H0"])
    par_best = params_from_H0(H0_best)
    flux_interp, F_best = make_flux_interp(df_flux, p_best)

    z = df_flux["z"].values
    F_raw = normalize_flux_history(z, df_flux["Omega_flux_4term_fit"].values, p=1.0)

    plt.figure(figsize=(10, 6))
    plt.plot(z, F_raw, "--", color="gray", linewidth=2,
             label="Raw 4-term Flux history, F(0)=1")
    plt.plot(z, F_best, color="white", linewidth=2.8,
             label=rf"Best metric response $F^p$, p={p_best:.4f}")
    plt.gca().invert_xaxis()
    plt.xlim(8, 0)
    plt.xlabel("Redshift z")
    plt.ylabel("F(z), normalized to today")
    plt.title(rf"Best CMB-Rigid Flux-CBH History: H0={H0_best:.2f}, p={p_best:.4f}")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig("flux_cbh_cmb_rigid_best_history.png", dpi=180)
    plt.close()

    # q history for best model
    zq = np.linspace(0.0, 5.0, 1200)
    q = q_flux_z(zq, par_best, flux_interp)
    OmF = omega_flux_z(zq, par_best, flux_interp)
    OmM = omega_m_z(zq, par_best, flux_interp)

    plt.figure(figsize=(10, 6))
    plt.plot(zq, q, color="red", linewidth=2, label="q(z)")
    plt.plot(zq, -q, color="orange", linestyle=":", linewidth=2, label="-q(z)")
    plt.plot(zq, OmF, color="magenta", linewidth=2, label=r"$\Omega_{\rm flux}(z)$")
    plt.plot(zq, OmM, color="cyan", linewidth=2, label=r"$\Omega_m(z)$")
    plt.axhline(0.0, color="white", alpha=0.4)
    if np.isfinite(best_row["z_acc"]):
        plt.axvline(best_row["z_acc"], color="white", linestyle="--",
                    label=rf"z_acc={best_row['z_acc']:.3f}")
    plt.gca().invert_xaxis()
    plt.xlim(5, 0)
    plt.xlabel("Redshift z")
    plt.ylabel("Cosmological diagnostic")
    plt.title("Best Flux-CBH Acceleration History")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig("flux_cbh_cmb_rigid_q_history.png", dpi=180)
    plt.close()


# ============================================================
# 7. Main
# ============================================================

def main():
    CSV_PATH = "flux_cbh_schechter_entropy.csv"

    df_flux = load_flux_csv(CSV_PATH)

    grid = run_grid_scan(df_flux)
    grid.to_csv("flux_cbh_cmb_rigid_grid_results.csv", index=False)

    accepted = grid[grid["accepted_theta"]].copy()

    # Save top model list:
    # if accepted rows exist, use them first; otherwise use score-sorted grid.
    if not accepted.empty:
        best = accepted.sort_values(["H0_distance_to_73", "theta_abs_err"]).head(100).reset_index(drop=True)
    else:
        best = grid.head(100).reset_index(drop=True)

    best.to_csv("flux_cbh_cmb_rigid_best_models.csv", index=False)

    lcdm_ref = lcdm_reference_scan()

    make_plots(df_flux, grid, best, lcdm_ref)

    print("\nCMB-Rigid Flux-CBH Grid Scan Complete")
    print("=" * 82)
    print("Physical densities locked:")
    print(f"  omega_b = Omega_b h^2 = {OMEGA_B_H2:.6f}")
    print(f"  omega_m = Omega_m h^2 = {OMEGA_M_H2:.6f}")
    print(f"  omega_r = Omega_r h^2 = {OMEGA_R_H2:.6e}")
    print()
    print(f"Theta target: {THETA_STAR_TARGET:.8f}")
    print(f"Theta tolerance used for accepted list: {THETA_TOL:.2e}")
    print()
    print("LCDM reference, same acoustic toy integral:")
    print(lcdm_ref.head(5)[[
        "H0", "theta_star", "theta_abs_err", "r_s_Mpc", "D_M_Mpc", "Omega_m0", "Omega_flux0"
    ]].to_string(index=False, float_format=lambda x: f"{x:.8g}"))
    print()

    print(f"Accepted theta models: {len(accepted):,} / {len(grid):,}")
    print()
    print("Best Flux-CBH models:")
    print(best.head(12)[[
        "p", "H0", "theta_star", "theta_abs_err", "theta_err_ppm",
        "r_s_Mpc", "D_M_Mpc", "Omega_m0", "Omega_flux0", "q0", "z_acc", "score"
    ]].to_string(index=False, float_format=lambda x: f"{x:.8g}"))
    print()
    print("Wrote:")
    print("  flux_cbh_cmb_rigid_grid_results.csv")
    print("  flux_cbh_cmb_rigid_best_models.csv")
    print("  flux_cbh_cmb_rigid_theta_heatmap.png")
    print("  flux_cbh_cmb_rigid_h0_vs_p.png")
    print("  flux_cbh_cmb_rigid_best_history.png")
    print("  flux_cbh_cmb_rigid_q_history.png")


if __name__ == "__main__":
    main()
