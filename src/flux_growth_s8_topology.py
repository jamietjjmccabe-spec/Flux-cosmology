#!/usr/bin/env python3
"""
flux_growth_s8_topology.py

Topology-aware Flux growth diagnostic.

Purpose
-------
The earlier smooth growth check showed that Flux-DE phantom-lite background evolution
does not automatically lower S8. This script adds the missing Flux ingredient:
substrate topology.

Core idea
---------
Retained metric memory should not act as a homogeneous μ_eff > 1 correction.
It should activate mainly inside collapsed basins / halo-like topological nodes.

Escaped entropy should not act only as a background density.
It can also produce a direct topological friction term in open substrate channels
(voids / IGM / diffuse cosmic web).

Growth equation
---------------
D'' + [2 + dlnH/dN + Gamma_topo(z,k)] D'
    - 3/2 * mu_eff(z,k) * Omega_m(z) * D = 0

where:
    mu_eff(z,k) = 1 + A_mu * G_collapse(k) * G_ret(z)
    Gamma_topo(z,k) = Gamma0 * G_void(k) * G_leak(z) * Sdot_norm(z)

Outputs
-------
flux_growth_s8_topology_results.csv
flux_growth_s8_topology_summary.csv
flux_growth_s8_topology.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

EPS = 1e-12

# -----------------------------
# Base cosmology / normalization
# -----------------------------
OMEGA_M0 = 0.30
OMEGA_DE0 = 0.70
SIGMA8_LCDM_ANCHOR = 0.811
S8_TARGET_LOW = 0.77
S8_TARGET_HIGH = 0.80

# Growth integration range
Z_INIT = 100.0
N_INIT = -math.log(1.0 + Z_INIT)
N_END = 0.0


# -----------------------------
# Star formation / throughput
# -----------------------------
def madau_dickinson_sfr(z: np.ndarray) -> np.ndarray:
    """
    Cosmic star formation rate density proxy.
    Madau & Dickinson-like shape, used here only as a normalized structural entropy proxy.
    """
    z = np.asarray(z)
    return 0.015 * ((1.0 + z) ** 2.7) / (1.0 + ((1.0 + z) / 2.9) ** 5.6)


def leakage_gate(z: np.ndarray, z_leak: float = 12.0, m: float = 8.0) -> np.ndarray:
    """
    Smooth activation of cosmic leakage, anchored to structure maturity / reionization proxy.
    """
    z = np.asarray(z)
    return 1.0 / (1.0 + ((1.0 + z) / (1.0 + z_leak)) ** m)


def simple_C_sigma(z: np.ndarray, C_high: float = 0.92, zc: float = 2.0, sharpness: float = 2.0) -> np.ndarray:
    """
    Toy coherence profile.
    High coherence at early times, lower coherence at late times / diffuse substrate.
    """
    z = np.asarray(z)
    # At high z: close to C_high. At low z: reduced.
    return C_high * ((1.0 + z) ** sharpness) / (((1.0 + z) ** sharpness) + ((1.0 + zc) ** sharpness))


def escape_gate(z: np.ndarray, n: float = 1.0, z_leak: float = 12.0, m: float = 8.0) -> np.ndarray:
    C = simple_C_sigma(z)
    return np.clip((1.0 - C) ** n, 0.0, 1.0) * leakage_gate(z, z_leak=z_leak, m=m)


def sdot_norm(z: np.ndarray, z_leak: float = 12.0, m: float = 8.0, n: float = 1.0) -> np.ndarray:
    """
    Normalized structural entropy leakage source.
    """
    src = madau_dickinson_sfr(z) * escape_gate(z, n=n, z_leak=z_leak, m=m)
    mx = np.nanmax(src)
    if mx <= 0 or not np.isfinite(mx):
        return np.zeros_like(src)
    return src / mx


# -----------------------------
# Background H(z)
# -----------------------------
def E_lcdm_z(z: np.ndarray, omega_m0: float = OMEGA_M0) -> np.ndarray:
    z = np.asarray(z)
    omega_de0 = 1.0 - omega_m0
    return np.sqrt(omega_m0 * (1.0 + z) ** 3 + omega_de0)


def w_flux_z(z: np.ndarray) -> np.ndarray:
    """
    Approximate phantom-lite fossil hump from previous bridge result.
    A compact toy profile:
    - w0 about -1.014
    - extra dip near z~1-3
    - returns closer to -1 at high z.
    """
    z = np.asarray(z)
    base = -1.0 - 0.014 * np.exp(-(z / 0.8) ** 2)
    hump = -0.020 * np.exp(-0.5 * ((z - 2.0) / 0.9) ** 2)
    return base + hump


def rho_de_flux_relative(z_grid: np.ndarray) -> np.ndarray:
    """
    Build rho_DE(z)/rho_DE(0) from w(z):
        d ln rho / d ln a = -3(1+w)
    Integrate from z=0 upward:
        d ln rho / dz = 3(1+w)/(1+z)
    """
    z = np.asarray(z_grid)
    order = np.argsort(z)
    zs = z[order]
    w = w_flux_z(zs)
    integrand = 3.0 * (1.0 + w) / (1.0 + zs)

    ln_rho = np.zeros_like(zs)
    for i in range(1, len(zs)):
        dz = zs[i] - zs[i-1]
        ln_rho[i] = ln_rho[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * dz

    rho = np.exp(ln_rho)
    out = np.empty_like(rho)
    out[order] = rho
    return out


def E_flux_z(z: np.ndarray, omega_m0: float = OMEGA_M0) -> np.ndarray:
    z = np.asarray(z)
    rho_rel = rho_de_flux_relative(z)
    omega_de0 = 1.0 - omega_m0
    return np.sqrt(omega_m0 * (1.0 + z) ** 3 + omega_de0 * rho_rel)


def background_arrays(model: str, n_grid: int = 3000) -> Dict[str, np.ndarray]:
    """
    Build N, z, E, Omega_m, dlnH/dN arrays.
    """
    N = np.linspace(N_INIT, N_END, n_grid)
    a = np.exp(N)
    z = 1.0 / a - 1.0

    if model == "lcdm":
        E = E_lcdm_z(z)
    elif model == "flux":
        E = E_flux_z(z)
    else:
        raise ValueError(f"Unknown model: {model}")

    Om = OMEGA_M0 * (1.0 + z) ** 3 / (E ** 2)
    dlnH_dN = np.gradient(np.log(E), N)

    return {"N": N, "a": a, "z": z, "E": E, "Omega_m": Om, "dlnH_dN": dlnH_dN}


# -----------------------------
# Topological gates
# -----------------------------
def sigmoid_logk(k: float, k0: float, sharpness: float) -> float:
    """
    Logistic sigmoid in log10(k/k0).
    """
    x = sharpness * math.log10(max(k, EPS) / max(k0, EPS))
    return 1.0 / (1.0 + math.exp(-x))


def G_collapse_k(k: float, k_collapse: float = 1.0, sharpness: float = 5.0) -> float:
    """
    Collapse gate active on small scales / high k.
    """
    return sigmoid_logk(k, k_collapse, sharpness)


def G_void_k(k: float, k_void: float = 0.25, sharpness: float = 5.0) -> float:
    """
    Void/open-channel gate active on large linear scales / low k.
    """
    return 1.0 - sigmoid_logk(k, k_void, sharpness)


def G_ret_z(z: np.ndarray) -> np.ndarray:
    """
    Retention gate: more retention at later collapsed states but suppressed on diffuse early modes.
    """
    z = np.asarray(z)
    # Smooth late activation. This is deliberately conservative.
    return 1.0 / (1.0 + ((1.0 + z) / 3.0) ** 3.0)


def mu_eff_zk(z: np.ndarray, k: float, A_mu: float = 0.15, k_collapse: float = 1.0) -> np.ndarray:
    """
    Retained-memory gravitational reinforcement.
    Active mainly for collapsed/high-k modes.
    """
    return 1.0 + A_mu * G_collapse_k(k, k_collapse=k_collapse) * G_ret_z(z)


def gamma_topo_zk(
    z: np.ndarray,
    k: float,
    Gamma0: float = 0.0,
    k_void: float = 0.25,
    z_leak: float = 12.0,
    m: float = 8.0,
    n: float = 1.0,
) -> np.ndarray:
    """
    Topological friction from open leakage channels.
    Active mainly on large linear scales and during structure entropy leakage.
    """
    return Gamma0 * G_void_k(k, k_void=k_void) * sdot_norm(z, z_leak=z_leak, m=m, n=n)


# -----------------------------
# Growth solver
# -----------------------------
def interp(x: float, xp: np.ndarray, fp: np.ndarray) -> float:
    return float(np.interp(x, xp, fp))


def solve_growth(
    model: str,
    k: float,
    Gamma0: float = 0.0,
    A_mu: float = 0.0,
    k_void: float = 0.25,
    k_collapse: float = 1.0,
    z_leak: float = 12.0,
    m: float = 8.0,
    n: float = 1.0,
) -> Dict[str, np.ndarray | float]:
    bg = background_arrays("lcdm" if model == "lcdm" else "flux")
    N = bg["N"]
    z = bg["z"]
    Om = bg["Omega_m"]
    dlnH = bg["dlnH_dN"]

    gamma_topo_arr = gamma_topo_zk(z, k, Gamma0=Gamma0, k_void=k_void, z_leak=z_leak, m=m, n=n)
    mu_arr = mu_eff_zk(z, k, A_mu=A_mu, k_collapse=k_collapse)

    def ode(Nval, y):
        D, Dp = y
        om = interp(Nval, N, Om)
        drag = 2.0 + interp(Nval, N, dlnH) + interp(Nval, N, gamma_topo_arr)
        mu = interp(Nval, N, mu_arr)
        return [Dp, -drag * Dp + 1.5 * mu * om * D]

    # Matter-era initial conditions: D ~ a, D' ~ D.
    D_init = math.exp(N_INIT)
    y0 = [D_init, D_init]

    sol = solve_ivp(ode, (N_INIT, N_END), y0, t_eval=N, rtol=1e-8, atol=1e-10)

    D_raw = sol.y[0]
    Dp_raw = sol.y[1]

    # Do NOT normalize D(a=1)=1 for S8 comparison.
    # Instead calibrate LCDM to sigma8 anchor and compare absolute growth ratios.
    f = Dp_raw / np.maximum(D_raw, EPS)
    gamma_index = np.full_like(f, np.nan)
    valid = (f > 0) & (Om > 0.01) & (Om < 0.99)
    gamma_index[valid] = np.log(f[valid]) / np.log(Om[valid])

    return {
        "N": N,
        "z": z,
        "a": bg["a"],
        "D": D_raw,
        "f": f,
        "Omega_m": Om,
        "dlnH_dN": dlnH,
        "gamma": gamma_index,
        "Gamma_topo": gamma_topo_arr,
        "mu_eff": mu_arr,
        "G_void": G_void_k(k, k_void=k_void),
        "G_collapse": G_collapse_k(k, k_collapse=k_collapse),
    }


def value_at_z(arr_z: np.ndarray, arr_y: np.ndarray, ztarget: float) -> float:
    order = np.argsort(arr_z)
    return float(np.interp(ztarget, arr_z[order], arr_y[order]))


def evaluate_case(
    name: str,
    model: str,
    k: float,
    Gamma0: float,
    A_mu: float,
    D_lcdm_today: float,
    k_void: float = 0.25,
    k_collapse: float = 1.0,
) -> Dict[str, float]:
    data = solve_growth(model=model, k=k, Gamma0=Gamma0, A_mu=A_mu, k_void=k_void, k_collapse=k_collapse)

    D_today = data["D"][-1]
    sigma8 = SIGMA8_LCDM_ANCHOR * (D_today / D_lcdm_today)
    S8 = sigma8 * math.sqrt(OMEGA_M0 / 0.3)

    zarr = data["z"]
    fs8_z05 = value_at_z(zarr, data["f"] * sigma8 * data["D"] / max(D_today, EPS), 0.5)
    fs8_z10 = value_at_z(zarr, data["f"] * sigma8 * data["D"] / max(D_today, EPS), 1.0)
    gamma_z0 = float(data["gamma"][-1])
    gamma_z05 = value_at_z(zarr, data["gamma"], 0.5)

    return {
        "case": name,
        "model": model,
        "k": k,
        "Gamma0": Gamma0,
        "A_mu": A_mu,
        "G_void_k": float(data["G_void"]),
        "G_collapse_k": float(data["G_collapse"]),
        "D_today": float(D_today),
        "sigma8_today": float(sigma8),
        "S8": float(S8),
        "f_sigma8_z0p5": fs8_z05,
        "f_sigma8_z1p0": fs8_z10,
        "gamma_z0": gamma_z0,
        "gamma_z0p5": gamma_z05,
        "Gamma_topo_z0": float(data["Gamma_topo"][-1]),
        "Gamma_topo_z1": value_at_z(zarr, data["Gamma_topo"], 1.0),
        "mu_eff_z0": float(data["mu_eff"][-1]),
        "mu_eff_z1": value_at_z(zarr, data["mu_eff"], 1.0),
    }


def scan_smoothing_window(
    gamma_values,
    A_mu_linear: float = 0.0,
    k_linear: float = 0.1,
    k_halo: float = 10.0,
    A_mu_halo: float = 0.25,
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    # LCDM baseline at linear scale.
    lcdm_linear = solve_growth(model="lcdm", k=k_linear, Gamma0=0.0, A_mu=0.0)
    D_lcdm_today = float(lcdm_linear["D"][-1])

    rows = []
    for Gamma0 in gamma_values:
        rows.append(evaluate_case(
            name="B_flux_topology_linear",
            model="flux",
            k=k_linear,
            Gamma0=Gamma0,
            A_mu=A_mu_linear,
            D_lcdm_today=D_lcdm_today,
        ))
        rows.append(evaluate_case(
            name="C_flux_topology_halo",
            model="flux",
            k=k_halo,
            Gamma0=Gamma0,
            A_mu=A_mu_halo,
            D_lcdm_today=D_lcdm_today,
        ))

    # Add baseline cases for reference.
    baseline_rows = [
        evaluate_case("A_lcdm_linear", "lcdm", k_linear, 0.0, 0.0, D_lcdm_today),
        evaluate_case("A_lcdm_halo", "lcdm", k_halo, 0.0, 0.0, D_lcdm_today),
    ]
    df = pd.DataFrame(baseline_rows + rows)

    return df, {"lcdm_linear": lcdm_linear, "D_lcdm_today": D_lcdm_today}


def pick_best(df: pd.DataFrame) -> pd.Series:
    lin = df[df["case"] == "B_flux_topology_linear"].copy()
    lin["target_distance"] = np.minimum(
        np.abs(lin["S8"] - S8_TARGET_LOW),
        np.abs(lin["S8"] - S8_TARGET_HIGH),
    )
    # prefer inside range, then closest, then lower Gamma0
    lin["inside_target"] = ((lin["S8"] >= S8_TARGET_LOW) & (lin["S8"] <= S8_TARGET_HIGH)).astype(int)
    lin = lin.sort_values(["inside_target", "target_distance", "Gamma0"], ascending=[False, True, True])
    return lin.iloc[0]


def save_plot(best_gamma0: float, outdir: Path) -> None:
    k_linear = 0.1
    k_halo = 10.0

    lcdm = solve_growth("lcdm", k_linear, 0.0, 0.0)
    flux_lin = solve_growth("flux", k_linear, best_gamma0, 0.0)
    flux_halo = solve_growth("flux", k_halo, best_gamma0, 0.25)

    # Convert z ordering for plots
    def sorted_xy(data, ykey):
        z = data["z"]
        y = data[ykey]
        order = np.argsort(z)
        return z[order], y[order]

    plt.figure(figsize=(12, 9))

    plt.subplot(2, 2, 1)
    z, D_lcdm = sorted_xy(lcdm, "D")
    _, D_lin = sorted_xy(flux_lin, "D")
    _, D_halo = sorted_xy(flux_halo, "D")
    # normalize for shape visualization only
    plt.plot(z, D_lcdm / D_lcdm[-1], label="LCDM linear")
    plt.plot(z, D_lin / D_lin[-1], label="Flux linear/topology")
    plt.plot(z, D_halo / D_halo[-1], label="Flux halo/topology")
    plt.xlim(0, 5)
    plt.xlabel("redshift z")
    plt.ylabel("D(z) / D(0)")
    plt.title("Growth history shape")
    plt.legend()

    plt.subplot(2, 2, 2)
    z, gamma_topo = sorted_xy(flux_lin, "Gamma_topo")
    _, mu_lin = sorted_xy(flux_lin, "mu_eff")
    _, mu_halo = sorted_xy(flux_halo, "mu_eff")
    plt.plot(z, gamma_topo, label="Gamma_topo linear")
    plt.plot(z, mu_lin - 1.0, label="mu_eff-1 linear")
    plt.plot(z, mu_halo - 1.0, label="mu_eff-1 halo")
    plt.xlim(0, 5)
    plt.xlabel("redshift z")
    plt.ylabel("gate strength")
    plt.title("Topology gates")
    plt.legend()

    plt.subplot(2, 2, 3)
    z, f_lcdm = sorted_xy(lcdm, "f")
    _, f_lin = sorted_xy(flux_lin, "f")
    _, f_halo = sorted_xy(flux_halo, "f")
    plt.plot(z, f_lcdm, label="LCDM linear")
    plt.plot(z, f_lin, label="Flux linear/topology")
    plt.plot(z, f_halo, label="Flux halo/topology")
    plt.xlim(0, 5)
    plt.xlabel("redshift z")
    plt.ylabel("f = dlnD/dlna")
    plt.title("Growth rate")
    plt.legend()

    plt.subplot(2, 2, 4)
    z, gam_lcdm = sorted_xy(lcdm, "gamma")
    _, gam_lin = sorted_xy(flux_lin, "gamma")
    _, gam_halo = sorted_xy(flux_halo, "gamma")
    plt.plot(z, gam_lcdm, label="LCDM linear")
    plt.plot(z, gam_lin, label="Flux linear/topology")
    plt.plot(z, gam_halo, label="Flux halo/topology")
    plt.xlim(0, 5)
    plt.ylim(0.45, 0.75)
    plt.xlabel("redshift z")
    plt.ylabel("growth index gamma")
    plt.title("Growth index")
    plt.legend()

    plt.tight_layout()
    plt.savefig(outdir / "flux_growth_s8_topology.png", dpi=180)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--gamma-min", type=float, default=0.0)
    parser.add_argument("--gamma-max", type=float, default=0.35)
    parser.add_argument("--gamma-steps", type=int, default=36)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gamma_values = np.linspace(args.gamma_min, args.gamma_max, args.gamma_steps)

    df, aux = scan_smoothing_window(gamma_values)
    best = pick_best(df)
    best_gamma0 = float(best["Gamma0"])

    # Add deltas relative to LCDM linear.
    lcdm_s8 = float(df.loc[df["case"] == "A_lcdm_linear", "S8"].iloc[0])
    df["delta_S8_vs_LCDM_percent"] = 100.0 * (df["S8"] - lcdm_s8) / lcdm_s8

    verdict = "PASS_SMOOTHING_WINDOW" if (S8_TARGET_LOW <= best["S8"] <= S8_TARGET_HIGH) else "NO_TARGET_WINDOW"
    summary = pd.DataFrame([{
        "verdict": verdict,
        "best_Gamma0": best_gamma0,
        "best_linear_S8": float(best["S8"]),
        "best_linear_sigma8": float(best["sigma8_today"]),
        "best_linear_delta_S8_percent": 100.0 * (float(best["S8"]) - lcdm_s8) / lcdm_s8,
        "lcdm_S8": lcdm_s8,
        "target_low": S8_TARGET_LOW,
        "target_high": S8_TARGET_HIGH,
        "interpretation": (
            "Topology friction can suppress linear S8 while the high-k halo track keeps retained-memory reinforcement active."
            if verdict == "PASS_SMOOTHING_WINDOW"
            else "The scanned topology friction range did not land inside the target S8 window."
        ),
    }])

    df.to_csv(outdir / "flux_growth_s8_topology_results.csv", index=False)
    summary.to_csv(outdir / "flux_growth_s8_topology_summary.csv", index=False)
    save_plot(best_gamma0, outdir)

    print("\nFlux topology growth scan")
    print("=" * 72)
    print(summary.to_string(index=False))
    print("\nReference rows")
    cols = ["case", "k", "Gamma0", "A_mu", "G_void_k", "G_collapse_k", "S8", "delta_S8_vs_LCDM_percent", "f_sigma8_z0p5", "gamma_z0", "Gamma_topo_z1", "mu_eff_z1"]
    print(df[(df["case"].str.startswith("A_")) | ((df["Gamma0"] == best_gamma0) & (df["case"].str.startswith(("B_", "C_"))))][cols].to_string(index=False))


if __name__ == "__main__":
    main()
