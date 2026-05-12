#!/usr/bin/env python3
"""
scan8_fixed_dynamic_exponent.py

Fix for scan8.py.

Original scan8 likely found no models because it used S_CBH_density_kB_mpc3
directly while the successful CMB scans used Omega_flux_4term_fit, and because
the known-good constant-p model was not cleanly included.

This version uses a running exponent:

    p_eff(z) = p_high + (p0 - p_high) / [1 + exp((z - z_t)/w)]

Metric response:

    F_metric(z) = F_raw(z) ** p_eff(z)

Outputs:
    scan8_fixed_results.csv
    scan8_fixed_best_models.csv
    scan8_fixed_best_plot.png
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d

c = 299792.458
z_star = 1089.92
theta_star_planck = 0.0104132

H0_flux = 73.0
omega_m_flux = 0.144
omega_b_flux = 0.02237
omega_r_flux = 4.183e-5

h_flux = H0_flux / 100.0
Om_flux = omega_m_flux / h_flux**2
Or_flux = omega_r_flux / h_flux**2
Oflux_0 = 1.0 - Om_flux - Or_flux

H0_lcdm = 67.4
Om_lcdm = 0.315
Ol_lcdm = 1.0 - Om_lcdm

THETA_PPM_SOFT = 2500.0
THETA_PPM_STRICT = 800.0

z_sn_eval = np.linspace(0.01, 2.3, 90)
anchor_mask = (z_sn_eval >= 0.01) & (z_sn_eval <= 0.05)

z_dist_cmb = np.concatenate([
    np.linspace(0.0, 8.0, 1800),
    np.linspace(8.0, z_star, 2500)
])
z_sound = np.logspace(np.log10(z_star), 6.0, 4500)

p0_vals = np.linspace(0.17, 0.25, 17)
p_high_vals = np.linspace(0.02, 0.22, 21)
zt_vals = np.linspace(0.6, 2.8, 23)
w_vals = np.linspace(0.35, 1.40, 8)

print("Loading Schechter / Flux entropy data...")
df = pd.read_csv("flux_cbh_schechter_entropy.csv").sort_values("z")

if "Omega_flux_4term_fit" in df.columns:
    source_col = "Omega_flux_4term_fit"
elif "S_CBH_density_kB_mpc3" in df.columns:
    source_col = "S_CBH_density_kB_mpc3"
else:
    raise ValueError("CSV needs Omega_flux_4term_fit or S_CBH_density_kB_mpc3")

z_data = df["z"].values
F_raw_data = df[source_col].values.astype(float)

idx0 = np.argmin(np.abs(z_data - 0.0))
F0 = F_raw_data[idx0]
if F0 <= 0:
    raise ValueError("F(0) is non-positive.")

F_norm_data = np.clip(F_raw_data / F0, 1e-90, None)

F_interp = interp1d(
    z_data,
    F_norm_data,
    kind="linear",
    bounds_error=False,
    fill_value=(F_norm_data[0], F_norm_data[-1]),
)

def F_raw(z):
    return np.clip(F_interp(z), 1e-90, None)

def p_eff(z, p0, p_high, zt, w):
    z = np.asarray(z, dtype=float)
    return p_high + (p0 - p_high) / (1.0 + np.exp((z - zt) / w))

def F_metric(z, p0, p_high, zt, w):
    F = F_raw(z)
    pe = p_eff(z, p0, p_high, zt, w)
    out = F ** pe
    F00 = F_raw(0.0) ** p_eff(0.0, p0, p_high, zt, w)
    return out / F00

def E_flux(z, p0, p_high, zt, w):
    z = np.asarray(z, dtype=float)
    flux_term = Oflux_0 * F_metric(z, p0, p_high, zt, w)
    return np.sqrt(Or_flux * (1.0 + z)**4 + Om_flux * (1.0 + z)**3 + flux_term)

def E_lcdm(z):
    z = np.asarray(z, dtype=float)
    return np.sqrt(Om_lcdm * (1.0 + z)**3 + Ol_lcdm)

def sound_speed_over_c(z):
    z = np.asarray(z, dtype=float)
    R = 31500.0 * omega_b_flux / (1.0 + z)
    return 1.0 / np.sqrt(3.0 * (1.0 + R))

def evaluate_cmb(p0, p_high, zt, w):
    E_d = E_flux(z_dist_cmb, p0, p_high, zt, w)
    D_M = (c / H0_flux) * trapezoid(1.0 / E_d, z_dist_cmb)

    E_s = E_flux(z_sound, p0, p_high, zt, w)
    r_s = (c / H0_flux) * trapezoid(sound_speed_over_c(z_sound) / E_s, z_sound)

    theta = r_s / D_M
    err_ppm = (theta - theta_star_planck) / theta_star_planck * 1e6
    return theta, err_ppm, r_s, D_M

def luminosity_distance_grid(z_grid, E_func, H0):
    z_grid = np.asarray(z_grid)
    invE = 1.0 / E_func(z_grid)
    Dc = np.zeros_like(z_grid)
    dz = np.diff(z_grid)
    avg = 0.5 * (invE[1:] + invE[:-1])
    Dc[1:] = np.cumsum(dz * avg)
    return (1.0 + z_grid) * (c / H0) * Dc

z_sn_int = np.linspace(0.0, 2.3, 1200)
Dl_lcdm_int = luminosity_distance_grid(z_sn_int, E_lcdm, H0_lcdm)
mu_lcdm_int = 5.0 * np.log10(np.clip(Dl_lcdm_int, 1e-30, None)) + 25.0
mu_lcdm_interp = interp1d(z_sn_int, mu_lcdm_int, bounds_error=False, fill_value="extrapolate")
mu_lcdm = mu_lcdm_interp(z_sn_eval)

def evaluate_sn_shape(p0, p_high, zt, w):
    def E_model(z):
        return E_flux(z, p0, p_high, zt, w)

    Dl_flux_int = luminosity_distance_grid(z_sn_int, E_model, H0_flux)
    mu_flux_int = 5.0 * np.log10(np.clip(Dl_flux_int, 1e-30, None)) + 25.0
    mu_flux = interp1d(z_sn_int, mu_flux_int, bounds_error=False, fill_value="extrapolate")(z_sn_eval)

    delta = mu_flux - mu_lcdm
    offset = np.mean(delta[anchor_mask])
    shape = delta - offset

    return {
        "sn_shape_max_abs": float(np.max(np.abs(shape))),
        "sn_shape_end": float(shape[-1]),
        "sn_offset": float(offset),
    }

print("\nKnown-good constant p sanity checks:")
for ptest in [0.190, 0.214]:
    theta, err_ppm, rs, dm = evaluate_cmb(ptest, ptest, 1.0, 1.0)
    sn = evaluate_sn_shape(ptest, ptest, 1.0, 1.0)
    print(
        f"  p={ptest:.3f}: theta={theta:.8f}, err={err_ppm:+.1f} ppm, "
        f"r_s={rs:.3f}, D_M={dm:.3f}, SN_shape={sn['sn_shape_max_abs']:.4f}"
    )

print("\nStarting dynamic-exponent scan (p0, p_high, z_t, w)...")
print(f"Source column: {source_col}")
print(f"Models: {len(p0_vals) * len(p_high_vals) * len(zt_vals) * len(w_vals):,}")

results = []

for p0 in p0_vals:
    print(f"  p0={p0:.4f}")
    for p_high in p_high_vals:
        for zt in zt_vals:
            for w in w_vals:
                theta, err_ppm, r_s, D_M = evaluate_cmb(p0, p_high, zt, w)
                abs_err_ppm = abs(err_ppm)

                if abs_err_ppm <= THETA_PPM_SOFT:
                    sn = evaluate_sn_shape(p0, p_high, zt, w)

                    score = (
                        (abs_err_ppm / 800.0)**2
                        + (sn["sn_shape_max_abs"] / 0.06)**2
                        + ((p0 - 0.214) / 0.05)**2
                    )

                    results.append({
                        "p0": p0,
                        "p_high": p_high,
                        "zt": zt,
                        "w": w,
                        "theta_star": theta,
                        "theta_err_ppm": err_ppm,
                        "theta_abs_err_ppm": abs_err_ppm,
                        "r_s_Mpc": r_s,
                        "D_M_Mpc": D_M,
                        "sn_shape_max_abs": sn["sn_shape_max_abs"],
                        "sn_shape_end": sn["sn_shape_end"],
                        "sn_offset": sn["sn_offset"],
                        "strict_cmb_pass": abs_err_ppm <= THETA_PPM_STRICT,
                        "score": score,
                    })

if not results:
    print("\nNo models passed even the soft CMB filter.")
    print("Check the known-good sanity values above.")
    raise SystemExit

df_res = pd.DataFrame(results).sort_values("score").reset_index(drop=True)
df_res.to_csv("scan8_fixed_results.csv", index=False)

strict = df_res[df_res["strict_cmb_pass"]].copy()
if strict.empty:
    print("\nNo strict CMB pass, saving nearest soft candidates.")
    best = df_res.head(100).copy()
else:
    best = strict.sort_values(["sn_shape_max_abs", "theta_abs_err_ppm"]).head(100).copy()

best.to_csv("scan8_fixed_best_models.csv", index=False)

print("\n=== TOP 12 DYNAMIC-EXPONENT MODELS ===")
print(best.head(12).to_string(index=False, float_format=lambda x: f"{x:.6g}"))

row = best.iloc[0]
p0_b = row["p0"]
ph_b = row["p_high"]
zt_b = row["zt"]
w_b = row["w"]

z_plot = np.linspace(0.0, 8.0, 800)
F_raw_plot = F_raw(z_plot)
F_best_plot = F_metric(z_plot, p0_b, ph_b, zt_b, w_b)
p_plot = p_eff(z_plot, p0_b, ph_b, zt_b, w_b)

def E_best(z):
    return E_flux(z, p0_b, ph_b, zt_b, w_b)

Dl_flux_int = luminosity_distance_grid(z_sn_int, E_best, H0_flux)
mu_flux_int = 5.0 * np.log10(np.clip(Dl_flux_int, 1e-30, None)) + 25.0
mu_flux = interp1d(z_sn_int, mu_flux_int, bounds_error=False, fill_value="extrapolate")(z_sn_eval)
delta = mu_flux - mu_lcdm
offset = np.mean(delta[anchor_mask])
shape = delta - offset

plt.style.use("dark_background")
fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=False)

axes[0].plot(z_plot, F_raw_plot, "--", color="gray", label="Raw Flux history F(0)=1")
axes[0].plot(z_plot, F_best_plot, color="white", linewidth=2.5, label="Best dynamic metric response")
axes[0].invert_xaxis()
axes[0].set_xlim(8, 0)
axes[0].set_ylabel("F(z)")
axes[0].set_title(
    f"Best scan8 fixed: p0={p0_b:.3f}, p_high={ph_b:.3f}, zt={zt_b:.2f}, w={w_b:.2f}"
)
axes[0].grid(alpha=0.25)
axes[0].legend()

axes[1].plot(z_plot, p_plot, color="cyan", linewidth=2.5)
axes[1].invert_xaxis()
axes[1].set_xlim(8, 0)
axes[1].set_ylabel(r"$p_{\rm eff}(z)$")
axes[1].grid(alpha=0.25)

axes[2].axhline(0, color="white", linestyle="--", alpha=0.6)
axes[2].axhspan(-0.05, 0.05, color="green", alpha=0.15, label="±0.05 mag")
axes[2].axhspan(-0.10, 0.10, color="gray", alpha=0.15, label="±0.10 mag")
axes[2].plot(z_sn_eval, shape, color="magenta", linewidth=2.5, label="SN shape residual")
axes[2].set_xlim(0, 2.3)
axes[2].set_xlabel("Redshift z")
axes[2].set_ylabel(r"$\Delta\mu_{\rm shape}$")
axes[2].grid(alpha=0.25)
axes[2].legend()

plt.tight_layout()
plt.savefig("scan8_fixed_best_plot.png", dpi=220)
print("\nWrote:")
print("  scan8_fixed_results.csv")
print("  scan8_fixed_best_models.csv")
print("  scan8_fixed_best_plot.png")
