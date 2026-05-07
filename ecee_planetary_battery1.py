#!/usr/bin/env python3
"""
ecee_planetary_battery.py

Flux-cosmology exoplanet screener for the "planetary battery" upgrade.

ECEM = Earth Coherence Envelope Match
    Similarity score: does the planet/star system resemble Earth's measured solution?

ECEE = Earth Coherence Envelope Equivalence
    Functional score: can different inputs produce an Earth-equivalent surface envelope?

Planetary Battery = core/internal-maintenance layer
    Estimates whether internal flux, atmosphere retention, magnetic/dynamo persistence,
    rocky boundary validity, and stellar stress combine into a long-lived coherent shell.

The point of this script is to surface the hidden class:
    larger/farther-out rocky or dense super-Earths that do not look exactly like Earth,
    but may reproduce an Earth-equivalent envelope through stronger internal retention.

Runs standalone with a built-in candidate/control catalog, or accepts a flexible CSV.

Examples:
    python ecee_planetary_battery.py
    python ecee_planetary_battery.py --input exoplanets.csv --top 50

Common NASA Exoplanet Archive aliases are accepted:
    pl_name, hostname, st_spectype, st_mass, st_rad, st_teff, st_lum,
    st_age, pl_bmasse, pl_rade, pl_orbsmax, pl_orbper, pl_orbeccen,
    pl_insol, pl_eqt, pl_dens

Outputs:
    ecee_planetary_battery_ranked.csv
    ecee_planetary_battery_scores.png
    ecee_usable_battery_scores.png
    ecee_core_to_stellar.png
    ecee_battery_vs_ecee.png
    ecee_mass_radius_battery.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False

EPS = 1e-12

EARTH = {
    "name": "Earth",
    "star_class": "G",
    "star_mass_solar": 1.0,
    "star_radius_solar": 1.0,
    "star_teff_k": 5772.0,
    "star_luminosity_solar": 1.0,
    "star_age_gyr": 4.57,
    "planet_mass_earth": 1.0,
    "planet_radius_earth": 1.0,
    "semi_major_axis_au": 1.0,
    "orbital_period_days": 365.25,
    "eccentricity": 0.0167,
    "insolation_earth": 1.0,
    "equilibrium_temp_k": 255.0,
    "density_earth": 1.0,
    "gravity_earth": 1.0,
    "atmosphere_score": 1.0,
    "magnetosphere_score": 1.0,
    "stellar_activity_score": 1.0,
    "rotation_period_days": 1.0,
    "volatile_flag": "",
}

STAR_CLASS_COHERENCE = {
    "O": 0.10, "B": 0.18, "A": 0.38, "F": 0.65,
    "G": 1.00, "K": 1.05, "M": 0.70, "UNKNOWN": 0.72,
}

ECEM_WEIGHTS = {
    "star": 0.16,
    "orbit": 0.18,
    "planet": 0.18,
    "thermal": 0.14,
    "atmosphere": 0.10,
    "magnetosphere": 0.10,
    "tidal": 0.08,
    "sigma": 0.06,
}

ECEE_WEIGHTS = {
    "thermal_balance": 0.22,
    "retention_capacity": 0.24,
    "rocky_boundary": 0.18,
    "magnetosphere_dynamo": 0.14,
    "tidal_stability": 0.10,
    "stellar_coherence": 0.08,
    "flux_history": 0.04,
}

BATTERY_WEIGHTS = {
    "core_heat": 0.22,
    "core_to_stellar_balance": 0.18,
    "geological_activity": 0.18,
    "retention_capacity": 0.16,
    "magnetosphere_dynamo": 0.12,
    "rocky_boundary": 0.10,
    "stellar_stress_resistance": 0.04,
}

ALIASES = {
    "name": ["name", "pl_name", "planet_name"],
    "host_name": ["host_name", "hostname", "star_name"],
    "star_class": ["star_class", "spectral_type", "st_spectype", "sp_type"],
    "star_mass_solar": ["star_mass_solar", "st_mass", "sy_smass"],
    "star_radius_solar": ["star_radius_solar", "st_rad"],
    "star_teff_k": ["star_teff_k", "st_teff", "teff"],
    "star_luminosity_solar": ["star_luminosity_solar", "st_luminosity_solar", "st_lum", "luminosity"],
    "star_age_gyr": ["star_age_gyr", "st_age", "age_gyr"],
    "planet_mass_earth": ["planet_mass_earth", "pl_bmasse", "pl_masse", "mass_earth"],
    "planet_radius_earth": ["planet_radius_earth", "pl_rade", "radius_earth"],
    "semi_major_axis_au": ["semi_major_axis_au", "pl_orbsmax", "a_au"],
    "orbital_period_days": ["orbital_period_days", "pl_orbper", "period_days"],
    "eccentricity": ["eccentricity", "pl_orbeccen", "e"],
    "insolation_earth": ["insolation_earth", "pl_insol", "insolation"],
    "equilibrium_temp_k": ["equilibrium_temp_k", "pl_eqt", "teq_k"],
    "density_earth": ["density_earth", "pl_dens_earth", "density_earth_units"],
    "density_g_cm3": ["density_g_cm3", "pl_dens", "density"],
    "gravity_earth": ["gravity_earth", "gravity", "surface_gravity_earth"],
    "atmosphere_score": ["atmosphere_score", "atm_score"],
    "magnetosphere_score": ["magnetosphere_score", "mag_score"],
    "stellar_activity_score": ["stellar_activity_score", "activity_score"],
    "rotation_period_days": ["rotation_period_days", "rot_period_days", "pl_rotper"],
    "volatile_flag": ["volatile_flag", "sub_neptune_flag", "volatile_rich"],
}


def finite_or_nan(x: Any) -> float:
    try:
        y = float(x)
        return y if np.isfinite(y) else np.nan
    except Exception:
        return np.nan


def clamp01(x: float) -> float:
    if not np.isfinite(x):
        return 0.0
    return float(np.clip(x, 0.0, 1.0))


def safe_log(x: float) -> float:
    return float(np.log(max(float(x), EPS)))


def log_gauss_ratio(value: float, target: float = 1.0, width: float = 0.5) -> float:
    if not np.isfinite(value) or not np.isfinite(target) or value <= 0 or target <= 0:
        return 0.0
    return clamp01(np.exp(-((np.log(value / target) / width) ** 2)))


def linear_gauss(value: float, target: float, width: float) -> float:
    if not np.isfinite(value):
        return 0.0
    return clamp01(np.exp(-(((value - target) / width) ** 2)))


def weighted_product(scores: Dict[str, float], weights: Dict[str, float]) -> float:
    total_w = sum(weights.values())
    if total_w <= 0:
        return 0.0
    acc = 0.0
    for key, w in weights.items():
        s = max(clamp01(scores.get(key, 0.0)), 1e-6)
        acc += (w / total_w) * safe_log(s)
    return clamp01(np.exp(acc))


def weighted_sum(scores: Dict[str, float], weights: Dict[str, float]) -> float:
    total_w = sum(weights.values())
    if total_w <= 0:
        return 0.0
    return clamp01(sum(clamp01(scores.get(k, 0.0)) * w for k, w in weights.items()) / total_w)


def asymmetric_log_score(value: float, target: float = 1.0, width_low: float = 1.0, width_high: float = 2.0) -> float:
    """Log-space Gaussian with different tolerance above/below the target."""
    if not np.isfinite(value) or value <= 0 or target <= 0:
        return 0.0
    ratio = value / target
    width = width_high if ratio >= 1 else width_low
    return clamp01(np.exp(-((np.log(ratio) / width) ** 2)))


def get_alias(row: pd.Series, canonical: str, default: Any = np.nan) -> Any:
    for col in ALIASES.get(canonical, [canonical]):
        if col in row.index and pd.notna(row[col]):
            return row[col]
    return default


def infer_star_class(value: Any, teff: float = np.nan) -> str:
    if isinstance(value, str) and value.strip():
        s = value.strip().upper()
        for c in ["O", "B", "A", "F", "G", "K", "M"]:
            if s.startswith(c) or c in s[:3]:
                return c
    if np.isfinite(teff):
        if teff >= 30000: return "O"
        if teff >= 10000: return "B"
        if teff >= 7500: return "A"
        if teff >= 6000: return "F"
        if teff >= 5200: return "G"
        if teff >= 3700: return "K"
        return "M"
    return "UNKNOWN"


def mass_from_radius_rocky(radius: float) -> float:
    if not np.isfinite(radius) or radius <= 0:
        return np.nan
    if radius <= 1.5:
        return radius ** 3.7
    if radius <= 2.0:
        return 2.7 * radius ** 1.3
    return 5.5 * (radius / 2.0) ** 0.7


def infer_luminosity(mass: float, radius: float, teff: float) -> float:
    if np.isfinite(radius) and np.isfinite(teff) and radius > 0 and teff > 0:
        return (radius ** 2) * ((teff / EARTH["star_teff_k"]) ** 4)
    if np.isfinite(mass) and mass > 0:
        if mass < 0.43:
            return 0.23 * mass ** 2.3
        if mass < 2.0:
            return mass ** 4.0
        return 1.5 * mass ** 3.5
    return np.nan


def infer_semi_major_axis(period_days: float, star_mass: float) -> float:
    if not np.isfinite(period_days) or period_days <= 0:
        return np.nan
    m = star_mass if np.isfinite(star_mass) and star_mass > 0 else 1.0
    period_years = period_days / 365.25
    return (m * period_years ** 2) ** (1.0 / 3.0)


def infer_period_days(a_au: float, star_mass: float) -> float:
    if not np.isfinite(a_au) or a_au <= 0:
        return np.nan
    m = star_mass if np.isfinite(star_mass) and star_mass > 0 else 1.0
    return 365.25 * math.sqrt(a_au ** 3 / m)


def tidal_strength_earth(row: Dict[str, Any] | pd.Series) -> float:
    mstar = finite_or_nan(row.get("star_mass_solar", np.nan))
    a = finite_or_nan(row.get("semi_major_axis_au", np.nan))
    if not np.isfinite(mstar) or mstar <= 0:
        mstar = 1.0
    if not np.isfinite(a) or a <= 0:
        return np.nan
    return mstar / max(a ** 3, EPS)


def infer_atmosphere_score(row: Dict[str, Any]) -> float:
    m = finite_or_nan(row.get("planet_mass_earth")); r = finite_or_nan(row.get("planet_radius_earth"))
    rho = finite_or_nan(row.get("density_earth")); insol = finite_or_nan(row.get("insolation_earth"))
    if not np.isfinite(m): m = 1.0
    if not np.isfinite(r): r = 1.0
    if not np.isfinite(rho): rho = m / max(r ** 3, EPS)
    if not np.isfinite(insol): insol = 1.0

    retention = log_gauss_ratio(max(m / max(r, EPS), EPS), target=1.0, width=1.5)
    rocky = 1.0 if rho >= 0.65 and r <= 2.0 else max(0.2, 1.0 - 0.35 * max(r - 1.8, 0.0))
    flux_window = 1.0
    if insol > 1.35: flux_window *= math.exp(-1.1 * (insol - 1.35) ** 2)
    if insol < 0.20: flux_window *= math.exp(-1.4 * (0.20 - insol) ** 2)
    return clamp01(0.55 * retention + 0.25 * rocky + 0.20 * flux_window)


def infer_magnetosphere_score(row: Dict[str, Any]) -> float:
    m = finite_or_nan(row.get("planet_mass_earth")); r = finite_or_nan(row.get("planet_radius_earth"))
    star_class = str(row.get("star_class", "UNKNOWN")).upper()
    a = finite_or_nan(row.get("semi_major_axis_au"))
    if not np.isfinite(m): m = 1.0
    if not np.isfinite(r): r = 1.0
    dyn = log_gauss_ratio(max(m / max(r ** 2, EPS), EPS), 1.0, 1.4)
    if 1.5 <= m <= 6.5 and r <= 2.0:
        dyn = max(dyn, 0.74)
    tidal = tidal_strength_earth(row)
    tidal_penalty = 1.0 if not np.isfinite(tidal) else math.exp(-0.13 * max(np.log10(max(tidal, 1.0)), 0.0) ** 2)
    class_penalty = 0.82 if star_class == "M" else 1.0
    if np.isfinite(a) and a > 0.65 and star_class in ["G", "K"]:
        class_penalty *= 1.03
    return clamp01(dyn * tidal_penalty * class_penalty)


def normalize_catalog(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        out: Dict[str, Any] = {}
        out["name"] = str(get_alias(row, "name", "Unknown"))
        out["host_name"] = str(get_alias(row, "host_name", ""))
        out["star_teff_k"] = finite_or_nan(get_alias(row, "star_teff_k"))
        out["star_class"] = infer_star_class(get_alias(row, "star_class", ""), out["star_teff_k"])
        out["star_mass_solar"] = finite_or_nan(get_alias(row, "star_mass_solar"))
        out["star_radius_solar"] = finite_or_nan(get_alias(row, "star_radius_solar"))
        out["star_luminosity_solar"] = finite_or_nan(get_alias(row, "star_luminosity_solar"))
        out["star_age_gyr"] = finite_or_nan(get_alias(row, "star_age_gyr"))
        out["planet_mass_earth"] = finite_or_nan(get_alias(row, "planet_mass_earth"))
        out["planet_radius_earth"] = finite_or_nan(get_alias(row, "planet_radius_earth"))
        out["semi_major_axis_au"] = finite_or_nan(get_alias(row, "semi_major_axis_au"))
        out["orbital_period_days"] = finite_or_nan(get_alias(row, "orbital_period_days"))
        out["eccentricity"] = finite_or_nan(get_alias(row, "eccentricity"))
        out["insolation_earth"] = finite_or_nan(get_alias(row, "insolation_earth"))
        out["equilibrium_temp_k"] = finite_or_nan(get_alias(row, "equilibrium_temp_k"))
        out["density_earth"] = finite_or_nan(get_alias(row, "density_earth"))
        density_g = finite_or_nan(get_alias(row, "density_g_cm3"))
        if not np.isfinite(out["density_earth"]) and np.isfinite(density_g):
            out["density_earth"] = density_g / 5.514
        out["gravity_earth"] = finite_or_nan(get_alias(row, "gravity_earth"))
        out["atmosphere_score"] = finite_or_nan(get_alias(row, "atmosphere_score"))
        out["magnetosphere_score"] = finite_or_nan(get_alias(row, "magnetosphere_score"))
        out["stellar_activity_score"] = finite_or_nan(get_alias(row, "stellar_activity_score"))
        out["rotation_period_days"] = finite_or_nan(get_alias(row, "rotation_period_days"))
        out["volatile_flag"] = str(get_alias(row, "volatile_flag", ""))

        if not np.isfinite(out["planet_mass_earth"]) and np.isfinite(out["planet_radius_earth"]):
            out["planet_mass_earth"] = mass_from_radius_rocky(out["planet_radius_earth"])
        if not np.isfinite(out["density_earth"]) and np.isfinite(out["planet_mass_earth"]) and np.isfinite(out["planet_radius_earth"]):
            out["density_earth"] = out["planet_mass_earth"] / max(out["planet_radius_earth"] ** 3, EPS)
        if not np.isfinite(out["gravity_earth"]) and np.isfinite(out["planet_mass_earth"]) and np.isfinite(out["planet_radius_earth"]):
            out["gravity_earth"] = out["planet_mass_earth"] / max(out["planet_radius_earth"] ** 2, EPS)
        if not np.isfinite(out["star_luminosity_solar"]):
            out["star_luminosity_solar"] = infer_luminosity(out["star_mass_solar"], out["star_radius_solar"], out["star_teff_k"])
        if not np.isfinite(out["semi_major_axis_au"]):
            out["semi_major_axis_au"] = infer_semi_major_axis(out["orbital_period_days"], out["star_mass_solar"])
        if not np.isfinite(out["orbital_period_days"]):
            out["orbital_period_days"] = infer_period_days(out["semi_major_axis_au"], out["star_mass_solar"])
        if not np.isfinite(out["insolation_earth"]) and np.isfinite(out["star_luminosity_solar"]) and np.isfinite(out["semi_major_axis_au"]):
            out["insolation_earth"] = out["star_luminosity_solar"] / max(out["semi_major_axis_au"] ** 2, EPS)
        if not np.isfinite(out["equilibrium_temp_k"]) and np.isfinite(out["insolation_earth"]):
            out["equilibrium_temp_k"] = 255.0 * out["insolation_earth"] ** 0.25
        if not np.isfinite(out["eccentricity"]):
            out["eccentricity"] = 0.05
        if not np.isfinite(out["star_age_gyr"]):
            out["star_age_gyr"] = {"F": 2.5, "G": 4.5, "K": 5.5, "M": 5.0}.get(out["star_class"], 4.5)
        if not np.isfinite(out["stellar_activity_score"]):
            out["stellar_activity_score"] = {"F": 0.70, "G": 0.95, "K": 0.98, "M": 0.55}.get(out["star_class"], 0.70)
        if not np.isfinite(out["atmosphere_score"]):
            out["atmosphere_score"] = infer_atmosphere_score(out)
        if not np.isfinite(out["magnetosphere_score"]):
            out["magnetosphere_score"] = infer_magnetosphere_score(out)
        rows.append(out)
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# ECEM: Earth resemblance
# -----------------------------------------------------------------------------

def star_score(row: pd.Series) -> float:
    c = str(row.get("star_class", "UNKNOWN")).upper()
    return weighted_sum(
        {
            "class": STAR_CLASS_COHERENCE.get(c, STAR_CLASS_COHERENCE["UNKNOWN"]),
            "mass": log_gauss_ratio(row.get("star_mass_solar", np.nan), 1.0, 0.45),
            "lum": log_gauss_ratio(row.get("star_luminosity_solar", np.nan), 1.0, 0.75),
            "teff": linear_gauss(row.get("star_teff_k", np.nan), 5772.0, 1200.0),
        },
        {"class": 0.35, "mass": 0.25, "lum": 0.25, "teff": 0.15},
    )


def orbit_score(row: pd.Series) -> float:
    ecc = finite_or_nan(row.get("eccentricity", 0.05))
    ecc_score = math.exp(-((ecc / 0.25) ** 2)) if np.isfinite(ecc) else 0.8
    return weighted_sum(
        {
            "a": log_gauss_ratio(row.get("semi_major_axis_au", np.nan), 1.0, 0.8),
            "period": log_gauss_ratio(row.get("orbital_period_days", np.nan), 365.25, 1.0),
            "insol": log_gauss_ratio(row.get("insolation_earth", np.nan), 1.0, 0.6),
            "ecc": ecc_score,
        },
        {"a": 0.25, "period": 0.15, "insol": 0.45, "ecc": 0.15},
    )


def planet_similarity_score(row: pd.Series) -> float:
    return weighted_sum(
        {
            "mass": log_gauss_ratio(row.get("planet_mass_earth", np.nan), 1.0, 0.9),
            "radius": log_gauss_ratio(row.get("planet_radius_earth", np.nan), 1.0, 0.55),
            "density": log_gauss_ratio(row.get("density_earth", np.nan), 1.0, 0.7),
            "gravity": log_gauss_ratio(row.get("gravity_earth", np.nan), 1.0, 0.7),
        },
        {"mass": 0.25, "radius": 0.30, "density": 0.20, "gravity": 0.25},
    )


def thermal_similarity_score(row: pd.Series) -> float:
    return weighted_sum(
        {
            "insol": log_gauss_ratio(row.get("insolation_earth", np.nan), 1.0, 0.55),
            "teq": log_gauss_ratio(row.get("equilibrium_temp_k", np.nan), 255.0, 0.30),
        },
        {"insol": 0.65, "teq": 0.35},
    )


def tidal_similarity_score(row: pd.Series) -> float:
    t = tidal_strength_earth(row)
    if not np.isfinite(t):
        return 0.5
    return log_gauss_ratio(t, 1.0, 1.25)


def score_sigma_proxy(row: pd.Series) -> float:
    insol = finite_or_nan(row.get("insolation_earth")); grav = finite_or_nan(row.get("gravity_earth"))
    atm = finite_or_nan(row.get("atmosphere_score")); mag = finite_or_nan(row.get("magnetosphere_score"))
    tidal = tidal_strength_earth(row)
    if not np.isfinite(insol): insol = 1.0
    if not np.isfinite(grav): grav = 1.0
    if not np.isfinite(atm): atm = 0.5
    if not np.isfinite(mag): mag = 0.5
    if not np.isfinite(tidal): tidal = 1.0
    sigma_proxy = (grav * (0.5 + atm) * (0.5 + mag)) / max((insol ** 0.2) * (tidal ** 0.08), EPS)
    return log_gauss_ratio(sigma_proxy, 2.0, 1.2)


def venus_branch_penalty(row: pd.Series) -> float:
    insol = finite_or_nan(row.get("insolation_earth", 1.0)); atm = finite_or_nan(row.get("atmosphere_score", 0.5))
    mag = finite_or_nan(row.get("magnetosphere_score", 0.5)); rot = finite_or_nan(row.get("rotation_period_days", np.nan))
    penalty = 1.0
    if np.isfinite(insol) and insol > 1.35:
        penalty *= np.exp(-2.5 * (insol - 1.35) ** 2)
    if np.isfinite(insol) and insol > 1.1 and np.isfinite(atm) and atm < 0.45:
        penalty *= 0.55
    if np.isfinite(insol) and insol > 1.1 and np.isfinite(mag) and mag < 0.45:
        penalty *= 0.70
    if np.isfinite(rot) and rot > 100:
        penalty *= 0.75
    return clamp01(penalty)


def compute_ecem_flux(row: pd.Series) -> Dict[str, float]:
    component = {
        "star": star_score(row),
        "orbit": orbit_score(row),
        "planet": planet_similarity_score(row),
        "thermal": thermal_similarity_score(row),
        "atmosphere": clamp01(row.get("atmosphere_score", 0.5)),
        "magnetosphere": clamp01(row.get("magnetosphere_score", 0.5)),
        "tidal": tidal_similarity_score(row),
        "sigma": score_sigma_proxy(row),
    }
    ecem_base = weighted_sum(component, ECEM_WEIGHTS)
    flux_tension = flux_tension_proxy(row)
    score_flux_tension = log_gauss_ratio(flux_tension, 1.0, 1.0)
    history = flux_history_score(row)
    surface = venus_branch_penalty(row)
    flux_adjustment = weighted_product(
        {"flux": score_flux_tension, "history": history, "surface": surface},
        {"flux": 0.35, "history": 0.35, "surface": 0.30},
    )
    return {
        **{f"score_ecem_{k}": v for k, v in component.items()},
        "ECEM_base": ecem_base,
        "flux_tension_proxy": flux_tension,
        "score_flux_tension": score_flux_tension,
        "score_flux_history": history,
        "score_surface_state": surface,
        "flux_adjustment": flux_adjustment,
        "ECEM_flux": clamp01(ecem_base * flux_adjustment),
    }


# -----------------------------------------------------------------------------
# ECEE: functional envelope equivalence
# -----------------------------------------------------------------------------

def flux_tension_proxy(row: pd.Series) -> float:
    insol = finite_or_nan(row.get("insolation_earth", 1.0)); activity = finite_or_nan(row.get("stellar_activity_score", 0.7))
    tidal = tidal_strength_earth(row)
    age = finite_or_nan(row.get("star_age_gyr", 4.5))
    if not np.isfinite(insol): insol = 1.0
    if not np.isfinite(activity): activity = 0.7
    if not np.isfinite(tidal): tidal = 1.0
    if not np.isfinite(age): age = 4.5
    activity_stress = 1.0 / max(activity, 0.05)
    tidal_stress = max(tidal, 1.0) ** 0.12
    age_factor = 1.0 + 0.04 * max(age - 4.57, 0.0)
    return float(max(insol * activity_stress * tidal_stress * age_factor, EPS))


def flux_history_score(row: pd.Series) -> float:
    c = str(row.get("star_class", "UNKNOWN")).upper()
    age = finite_or_nan(row.get("star_age_gyr", 4.5))
    activity = finite_or_nan(row.get("stellar_activity_score", 0.7))
    class_life = {"F": 0.60, "G": 0.92, "K": 1.00, "M": 0.78}.get(c, 0.72)
    if not np.isfinite(age): age = 4.5
    if not np.isfinite(activity): activity = 0.7
    maturity = 1.0 - math.exp(-age / 1.5)
    old_penalty = math.exp(-0.025 * max(age - 8.0, 0.0) ** 2)
    return clamp01(class_life * activity * maturity * old_penalty)


def rocky_boundary_score(row: pd.Series) -> float:
    r = finite_or_nan(row.get("planet_radius_earth")); m = finite_or_nan(row.get("planet_mass_earth"))
    rho = finite_or_nan(row.get("density_earth")); volatile = str(row.get("volatile_flag", "")).lower()
    if not np.isfinite(r): r = 1.0
    if not np.isfinite(m): m = mass_from_radius_rocky(r)
    if not np.isfinite(rho): rho = m / max(r ** 3, EPS)
    score = 1.0
    if r > 2.0: score *= np.exp(-1.5 * (r - 2.0) ** 2)
    if r > 2.5: score *= 0.35
    if rho < 0.65: score *= np.exp(-2.0 * (0.65 - rho) ** 2)
    if rho < 0.45: score *= 0.35
    if m > 10.0: score *= 0.55
    if any(x in volatile for x in ["sub", "neptune", "volatile", "h/he", "gas"]):
        score *= 0.25
    return clamp01(score)


def retention_raw(row: pd.Series) -> float:
    m = finite_or_nan(row.get("planet_mass_earth")); g = finite_or_nan(row.get("gravity_earth"))
    insol = finite_or_nan(row.get("insolation_earth")); atm = finite_or_nan(row.get("atmosphere_score"))
    mag = finite_or_nan(row.get("magnetosphere_score")); tidal = tidal_strength_earth(row)
    if not np.isfinite(m): m = 1.0
    if not np.isfinite(g): g = 1.0
    if not np.isfinite(insol): insol = 1.0
    if not np.isfinite(atm): atm = 0.5
    if not np.isfinite(mag): mag = 0.5
    if not np.isfinite(tidal): tidal = 1.0
    return float(max((m * g * atm * mag) / max(insol * (tidal ** 0.15), EPS), EPS))


def retention_capacity_score(row: pd.Series) -> float:
    raw = retention_raw(row)
    earth_raw = retention_raw(pd.Series(EARTH))
    return asymmetric_log_score(raw, earth_raw, width_low=1.0, width_high=2.5)


def thermal_envelope_proxy(row: pd.Series) -> float:
    insol = finite_or_nan(row.get("insolation_earth")); raw = retention_raw(row)
    if not np.isfinite(insol): insol = 1.0
    greenhouse_retention = 1.0 + 0.35 * np.log1p(raw)
    return float(max(insol * greenhouse_retention, EPS))


def thermal_balance_score(row: pd.Series) -> float:
    envelope = thermal_envelope_proxy(row)
    earth_envelope = thermal_envelope_proxy(pd.Series(EARTH))
    score = log_gauss_ratio(envelope, earth_envelope, 0.75)
    # Surface runaway/freeze-out sanity gates.
    insol = finite_or_nan(row.get("insolation_earth")); surface = venus_branch_penalty(row)
    if np.isfinite(insol) and insol < 0.25:
        score *= math.exp(-1.0 * (0.25 - insol) ** 2)
    return clamp01(score * max(surface, 0.10))


def magnetosphere_dynamo_score(row: pd.Series) -> float:
    mag = clamp01(row.get("magnetosphere_score", 0.5))
    core = core_heat_relative(row)
    core_support = asymmetric_log_score(core, 1.0, width_low=1.0, width_high=1.7)
    tidal = tidal_strength_earth(row)
    tidal_penalty = 1.0 if not np.isfinite(tidal) else math.exp(-0.10 * max(np.log10(max(tidal, 1.0)), 0.0) ** 2)
    return clamp01(0.60 * mag + 0.40 * core_support * tidal_penalty)


def tidal_stability_score(row: pd.Series) -> float:
    tidal = tidal_strength_earth(row)
    ecc = finite_or_nan(row.get("eccentricity", 0.05))
    if not np.isfinite(tidal): tidal = 1.0
    if not np.isfinite(ecc): ecc = 0.05
    low_penalty = asymmetric_log_score(tidal, 1.0, width_low=2.5, width_high=1.2)
    ecc_penalty = math.exp(-((ecc / 0.28) ** 2))
    # Do not eliminate all compact M branches; just classify them as different.
    return clamp01(0.75 * low_penalty + 0.25 * ecc_penalty)


def stellar_coherence_score(row: pd.Series) -> float:
    c = str(row.get("star_class", "UNKNOWN")).upper()
    activity = finite_or_nan(row.get("stellar_activity_score", 0.7))
    if not np.isfinite(activity): activity = 0.7
    return clamp01(STAR_CLASS_COHERENCE.get(c, 0.72) * activity)


def compute_ecee(row: pd.Series) -> Dict[str, float]:
    scores = {
        "thermal_balance": thermal_balance_score(row),
        "retention_capacity": retention_capacity_score(row),
        "rocky_boundary": rocky_boundary_score(row),
        "magnetosphere_dynamo": magnetosphere_dynamo_score(row),
        "tidal_stability": tidal_stability_score(row),
        "stellar_coherence": stellar_coherence_score(row),
        "flux_history": flux_history_score(row),
    }
    return {
        **{f"score_ecee_{k}": v for k, v in scores.items()},
        "retention_raw": retention_raw(row),
        "thermal_envelope_proxy": thermal_envelope_proxy(row),
        "tidal_strength_earth": tidal_strength_earth(row),
        "ECEE_score": weighted_product(scores, ECEE_WEIGHTS),
    }


# -----------------------------------------------------------------------------
# Planetary battery: internal flux / core maintenance layer
# -----------------------------------------------------------------------------

def core_heat_relative(row: pd.Series) -> float:
    """
    Relative internal heat-flux capacity.

    First-order proxy:
        heat production/storage ~ mass and radiogenic inventory
        escape area ~ radius^2
        old systems cool slowly; apply a mild radiogenic decay with age
        dense rocky bodies get a small boost; volatile low-density worlds are penalized
    """
    m = finite_or_nan(row.get("planet_mass_earth")); r = finite_or_nan(row.get("planet_radius_earth"))
    rho = finite_or_nan(row.get("density_earth")); age = finite_or_nan(row.get("star_age_gyr"))
    if not np.isfinite(m): m = 1.0
    if not np.isfinite(r): r = 1.0
    if not np.isfinite(rho): rho = m / max(r ** 3, EPS)
    if not np.isfinite(age): age = 4.57

    area_flux = m / max(r ** 2, EPS)
    thermal_inertia = max(m, EPS) ** 0.20
    age_decay = math.exp(-(age - EARTH["star_age_gyr"]) / 12.0)
    density_factor = float(np.clip(rho ** 0.25, 0.55, 1.35))
    rocky_gate = rocky_boundary_score(row) ** 0.35
    return float(max(area_flux * thermal_inertia * age_decay * density_factor * rocky_gate, EPS))


def core_to_stellar_flux_ratio_relative(row: pd.Series) -> float:
    """Core/starlight ratio relative to Earth's core/starlight ratio."""
    core = core_heat_relative(row)
    insol = finite_or_nan(row.get("insolation_earth"))
    if not np.isfinite(insol): insol = 1.0
    return float(max(core / max(insol, EPS), EPS))


def battery_lifetime_proxy(row: pd.Series) -> float:
    """
    Proxy for how long the planetary battery can keep supplying maintenance energy.
    Larger rocky planets have more reservoir and cool more slowly, but volatile/gas worlds
    and very high-mass bodies are penalized.
    """
    m = finite_or_nan(row.get("planet_mass_earth")); r = finite_or_nan(row.get("planet_radius_earth"))
    rho = finite_or_nan(row.get("density_earth")); age = finite_or_nan(row.get("star_age_gyr"))
    if not np.isfinite(m): m = 1.0
    if not np.isfinite(r): r = 1.0
    if not np.isfinite(rho): rho = m / max(r ** 3, EPS)
    if not np.isfinite(age): age = 4.57
    reservoir = max(m, EPS) ** 0.55 / max(r ** 0.25, EPS)
    cooling = max(r, EPS) ** 0.35
    radiogenic_remaining = math.exp(-(age - EARTH["star_age_gyr"]) / 18.0)
    rocky = rocky_boundary_score(row)
    if m > 8.0:
        rocky *= 0.70
    return float(max(reservoir * cooling * radiogenic_remaining * (0.35 + 0.65 * rocky) * (rho ** 0.15), EPS))


def score_core_heat(row: pd.Series) -> float:
    # Extra internal heat is tolerated; too little drops quickly. Too much is mildly penalized.
    return asymmetric_log_score(core_heat_relative(row), 1.0, width_low=0.85, width_high=1.85)


def score_core_to_stellar_balance(row: pd.Series) -> float:
    # Higher core/star ratio can support farther-out branches, but extreme ratios may imply ice shell/ocean world branch.
    return asymmetric_log_score(core_to_stellar_flux_ratio_relative(row), 1.0, width_low=0.85, width_high=1.65)


def score_geological_activity(row: pd.Series) -> float:
    core = core_heat_relative(row)
    lifetime = battery_lifetime_proxy(row)
    rocky = rocky_boundary_score(row)
    gravity = finite_or_nan(row.get("gravity_earth"))
    if not np.isfinite(gravity): gravity = 1.0
    pressure_gate = 1.0
    if gravity > 2.5:
        pressure_gate *= math.exp(-0.25 * (gravity - 2.5) ** 2)
    if gravity > 3.5:
        pressure_gate *= 0.60
    core_window = asymmetric_log_score(core, 1.0, width_low=0.90, width_high=1.90)
    lifetime_score = asymmetric_log_score(lifetime, 1.0, width_low=1.20, width_high=2.20)
    return clamp01((0.45 * core_window + 0.35 * lifetime_score + 0.20 * rocky) * pressure_gate)


def score_stellar_stress_resistance(row: pd.Series) -> float:
    activity = finite_or_nan(row.get("stellar_activity_score", 0.7))
    tidal = tidal_strength_earth(row)
    mag = magnetosphere_dynamo_score(row)
    if not np.isfinite(activity): activity = 0.7
    if not np.isfinite(tidal): tidal = 1.0
    stress = (1.0 / max(activity, 0.05)) * max(tidal, 1.0) ** 0.10
    resistance = (0.45 + 0.55 * mag) / stress
    return asymmetric_log_score(resistance, 1.0, width_low=1.4, width_high=1.8)


def compute_planetary_battery(row: pd.Series) -> Dict[str, float]:
    scores = {
        "core_heat": score_core_heat(row),
        "core_to_stellar_balance": score_core_to_stellar_balance(row),
        "geological_activity": score_geological_activity(row),
        "retention_capacity": retention_capacity_score(row),
        "magnetosphere_dynamo": magnetosphere_dynamo_score(row),
        "rocky_boundary": rocky_boundary_score(row),
        "stellar_stress_resistance": score_stellar_stress_resistance(row),
    }
    battery = weighted_product(scores, BATTERY_WEIGHTS)

    # A large volatile/sub-Neptune world can carry a big internal reservoir while
    # still being the wrong *usable* boundary for an Earth-equivalent surface
    # coherence shell.  The rocky boundary gate converts raw battery capacity into
    # usable rocky-envelope battery capacity.
    rocky_boundary_gate = scores["rocky_boundary"]
    usable_battery_score = clamp01(battery * rocky_boundary_gate)

    ecee = finite_or_nan(row.get("ECEE_score", 0.0))
    ecem_flux = finite_or_nan(row.get("ECEM_flux", 0.0))
    # Battery-enhanced equivalence. Similarity remains visible but no longer dominates.
    # Use the usable battery so wrong-boundary controls cannot dominate the composite.
    battery_coherence_score = clamp01(0.25 * ecem_flux + 0.35 * ecee + 0.40 * usable_battery_score)
    return {
        **{f"score_battery_{k}": v for k, v in scores.items()},
        "core_heat_relative": core_heat_relative(row),
        "core_to_stellar_ratio_relative": core_to_stellar_flux_ratio_relative(row),
        "battery_lifetime_proxy": battery_lifetime_proxy(row),
        "rocky_boundary_gate": rocky_boundary_gate,
        "planetary_battery_score": battery,
        "usable_battery_score": usable_battery_score,
        "battery_coherence_score": battery_coherence_score,
    }


def classify(row: pd.Series) -> str:
    name = str(row.get("name", "")).lower()
    star_class = str(row.get("star_class", "UNKNOWN")).upper()
    r = finite_or_nan(row.get("planet_radius_earth")); m = finite_or_nan(row.get("planet_mass_earth"))
    insol = finite_or_nan(row.get("insolation_earth")); tidal = finite_or_nan(row.get("tidal_strength_earth"))
    rocky = finite_or_nan(row.get("score_ecee_rocky_boundary"))
    battery = finite_or_nan(row.get("planetary_battery_score"))
    usable_battery = finite_or_nan(row.get("usable_battery_score"))
    ecee = finite_or_nan(row.get("ECEE_score")); ecem = finite_or_nan(row.get("ECEM_flux"))
    core_star = finite_or_nan(row.get("core_to_stellar_ratio_relative")); surface = finite_or_nan(row.get("score_surface_state"))
    volatile = str(row.get("volatile_flag", "")).lower()

    if any(x in volatile for x in ["sub", "neptune", "volatile", "h/he", "gas"]) or (np.isfinite(r) and r > 2.4 and np.isfinite(rocky) and rocky < 0.45):
        return "volatile/sub-Neptune wrong-boundary control"
    if "venus" in name or (np.isfinite(insol) and insol > 1.25 and np.isfinite(surface) and surface < 0.25):
        return "Venus/runaway over-flux branch"
    if "mars" in name or (np.isfinite(insol) and insol < 0.50 and np.isfinite(m) and m < 0.5):
        return "Mars/freeze-out weak-battery branch"
    if np.isfinite(ecem) and ecem > 0.70 and np.isfinite(ecee) and ecee > 0.70:
        return "Earth-like and Earth-equivalent benchmark"
    if star_class in ["G", "K"] and np.isfinite(r) and 1.2 <= r <= 2.0 and np.isfinite(usable_battery) and usable_battery > 0.60 and np.isfinite(core_star) and core_star > 1.0:
        return "outer super-Earth usable rocky-battery candidate"
    if star_class == "K" and np.isfinite(usable_battery) and usable_battery > 0.58:
        return "K-dwarf long-coherence usable-battery branch"
    if star_class == "M" and np.isfinite(tidal) and tidal > 100 and np.isfinite(usable_battery) and usable_battery > 0.45:
        return "dense M-dwarf retained usable-battery branch"
    if star_class == "M" and np.isfinite(tidal) and tidal > 100:
        return "compact M-dwarf high-tidal branch"
    if np.isfinite(ecee) and ecee > 0.60 and np.isfinite(usable_battery) and usable_battery > 0.55:
        return "non-Earth coherent usable-envelope candidate"
    return "partial / insufficient-data coherence candidate"


# -----------------------------------------------------------------------------
# Built-in catalog
# -----------------------------------------------------------------------------

def built_in_catalog() -> pd.DataFrame:
    data = [
        EARTH,
        {"name":"Venus", "star_class":"G", "star_mass_solar":1.0, "star_radius_solar":1.0, "star_teff_k":5772, "star_luminosity_solar":1.0, "star_age_gyr":4.57, "planet_mass_earth":0.815, "planet_radius_earth":0.949, "semi_major_axis_au":0.723, "orbital_period_days":224.7, "eccentricity":0.007, "insolation_earth":1.91, "equilibrium_temp_k":260, "density_earth":0.95, "gravity_earth":0.904, "atmosphere_score":0.08, "magnetosphere_score":0.10, "stellar_activity_score":1.0, "rotation_period_days":243},
        {"name":"Mars", "star_class":"G", "star_mass_solar":1.0, "star_radius_solar":1.0, "star_teff_k":5772, "star_luminosity_solar":1.0, "star_age_gyr":4.57, "planet_mass_earth":0.107, "planet_radius_earth":0.532, "semi_major_axis_au":1.524, "orbital_period_days":687, "eccentricity":0.093, "insolation_earth":0.43, "equilibrium_temp_k":210, "density_earth":0.71, "gravity_earth":0.379, "atmosphere_score":0.18, "magnetosphere_score":0.12, "stellar_activity_score":1.0, "rotation_period_days":1.03},
        {"name":"Kepler-452 b", "star_class":"G", "star_mass_solar":1.04, "star_radius_solar":1.11, "star_teff_k":5757, "star_luminosity_solar":1.20, "star_age_gyr":6.0, "planet_mass_earth":5.0, "planet_radius_earth":1.63, "semi_major_axis_au":1.046, "orbital_period_days":384.8, "eccentricity":0.05, "insolation_earth":1.10, "equilibrium_temp_k":265, "density_earth":1.15, "gravity_earth":1.88, "atmosphere_score":0.72, "magnetosphere_score":0.75, "stellar_activity_score":0.92},
        {"name":"Kepler-442 b", "star_class":"K", "star_mass_solar":0.61, "star_radius_solar":0.60, "star_teff_k":4402, "star_luminosity_solar":0.12, "star_age_gyr":5.0, "planet_mass_earth":2.36, "planet_radius_earth":1.34, "semi_major_axis_au":0.409, "orbital_period_days":112.3, "eccentricity":0.04, "insolation_earth":0.70, "equilibrium_temp_k":233, "density_earth":0.98, "gravity_earth":1.31, "atmosphere_score":0.72, "magnetosphere_score":0.70, "stellar_activity_score":0.96},
        {"name":"Kepler-62 e", "star_class":"K", "star_mass_solar":0.69, "star_radius_solar":0.64, "star_teff_k":4925, "star_luminosity_solar":0.21, "star_age_gyr":7.0, "planet_mass_earth":4.5, "planet_radius_earth":1.61, "semi_major_axis_au":0.427, "orbital_period_days":122.4, "eccentricity":0.05, "insolation_earth":1.20, "equilibrium_temp_k":267, "density_earth":1.08, "gravity_earth":1.74, "atmosphere_score":0.68, "magnetosphere_score":0.73, "stellar_activity_score":0.95},
        {"name":"Kepler-62 f", "star_class":"K", "star_mass_solar":0.69, "star_radius_solar":0.64, "star_teff_k":4925, "star_luminosity_solar":0.21, "star_age_gyr":7.0, "planet_mass_earth":3.0, "planet_radius_earth":1.41, "semi_major_axis_au":0.718, "orbital_period_days":267.3, "eccentricity":0.05, "insolation_earth":0.41, "equilibrium_temp_k":204, "density_earth":1.07, "gravity_earth":1.51, "atmosphere_score":0.66, "magnetosphere_score":0.72, "stellar_activity_score":0.95},
        {"name":"Outer Super-Earth Test", "star_class":"K", "star_mass_solar":0.75, "star_radius_solar":0.75, "star_teff_k":5000, "star_luminosity_solar":0.35, "star_age_gyr":6.0, "planet_mass_earth":4.5, "planet_radius_earth":1.60, "semi_major_axis_au":0.85, "orbital_period_days":330, "eccentricity":0.05, "insolation_earth":0.48, "equilibrium_temp_k":212, "density_earth":1.10, "gravity_earth":1.76, "atmosphere_score":0.75, "magnetosphere_score":0.80, "stellar_activity_score":0.97},
        {"name":"Cold Dense Super-Earth Test", "star_class":"K", "star_mass_solar":0.72, "star_radius_solar":0.72, "star_teff_k":4850, "star_luminosity_solar":0.28, "star_age_gyr":6.5, "planet_mass_earth":5.2, "planet_radius_earth":1.72, "semi_major_axis_au":0.95, "orbital_period_days":398, "eccentricity":0.04, "insolation_earth":0.31, "equilibrium_temp_k":190, "density_earth":1.02, "gravity_earth":1.76, "atmosphere_score":0.78, "magnetosphere_score":0.82, "stellar_activity_score":0.96},
        {"name":"LHS 1140 b", "star_class":"M", "star_mass_solar":0.18, "star_radius_solar":0.21, "star_teff_k":3216, "star_luminosity_solar":0.0044, "star_age_gyr":5.0, "planet_mass_earth":5.6, "planet_radius_earth":1.73, "semi_major_axis_au":0.094, "orbital_period_days":24.7, "eccentricity":0.04, "insolation_earth":0.50, "equilibrium_temp_k":230, "density_earth":1.08, "gravity_earth":1.87, "atmosphere_score":0.68, "magnetosphere_score":0.55, "stellar_activity_score":0.62},
        {"name":"TOI-700 d", "star_class":"M", "star_mass_solar":0.42, "star_radius_solar":0.42, "star_teff_k":3480, "star_luminosity_solar":0.023, "star_age_gyr":2.0, "planet_mass_earth":1.7, "planet_radius_earth":1.19, "semi_major_axis_au":0.163, "orbital_period_days":37.4, "eccentricity":0.05, "insolation_earth":0.86, "equilibrium_temp_k":246, "density_earth":1.01, "gravity_earth":1.20, "atmosphere_score":0.72, "magnetosphere_score":0.50, "stellar_activity_score":0.68},
        {"name":"TRAPPIST-1 e", "star_class":"M", "star_mass_solar":0.089, "star_radius_solar":0.12, "star_teff_k":2559, "star_luminosity_solar":0.00055, "star_age_gyr":7.6, "planet_mass_earth":0.77, "planet_radius_earth":0.92, "semi_major_axis_au":0.029, "orbital_period_days":6.1, "eccentricity":0.006, "insolation_earth":0.66, "equilibrium_temp_k":230, "density_earth":0.99, "gravity_earth":0.91, "atmosphere_score":0.55, "magnetosphere_score":0.35, "stellar_activity_score":0.38},
        {"name":"Proxima Centauri b", "star_class":"M", "star_mass_solar":0.122, "star_radius_solar":0.154, "star_teff_k":3042, "star_luminosity_solar":0.00155, "star_age_gyr":4.8, "planet_mass_earth":1.27, "planet_radius_earth":1.08, "semi_major_axis_au":0.0485, "orbital_period_days":11.2, "eccentricity":0.11, "insolation_earth":0.65, "equilibrium_temp_k":234, "density_earth":1.00, "gravity_earth":1.09, "atmosphere_score":0.48, "magnetosphere_score":0.30, "stellar_activity_score":0.28},
        {"name":"K2-18 b", "star_class":"M", "star_mass_solar":0.50, "star_radius_solar":0.47, "star_teff_k":3457, "star_luminosity_solar":0.023, "star_age_gyr":2.4, "planet_mass_earth":8.6, "planet_radius_earth":2.61, "semi_major_axis_au":0.142, "orbital_period_days":33, "eccentricity":0.05, "insolation_earth":1.0, "equilibrium_temp_k":255, "density_earth":0.48, "gravity_earth":1.26, "atmosphere_score":0.40, "magnetosphere_score":0.45, "stellar_activity_score":0.55, "volatile_flag":"sub-Neptune"},
    ]
    return pd.DataFrame(data)


# -----------------------------------------------------------------------------
# Scoring/output
# -----------------------------------------------------------------------------

def score_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    norm = normalize_catalog(df)
    rows = []
    for _, row in norm.iterrows():
        base = row.to_dict()
        ecem = compute_ecem_flux(row)
        tmp1 = pd.Series({**base, **ecem})
        ecee = compute_ecee(tmp1)
        tmp2 = pd.Series({**base, **ecem, **ecee})
        battery = compute_planetary_battery(tmp2)
        out = {**base, **ecem, **ecee, **battery}
        out["branch_classification"] = classify(pd.Series(out))
        rows.append(out)
    outdf = pd.DataFrame(rows)
    outdf["rank_ECEM_flux"] = outdf["ECEM_flux"].rank(ascending=False, method="min").astype(int)
    outdf["rank_ECEE"] = outdf["ECEE_score"].rank(ascending=False, method="min").astype(int)
    outdf["rank_battery"] = outdf["planetary_battery_score"].rank(ascending=False, method="min").astype(int)
    outdf["rank_usable_battery"] = outdf["usable_battery_score"].rank(ascending=False, method="min").astype(int)
    outdf["rank_battery_coherence"] = outdf["battery_coherence_score"].rank(ascending=False, method="min").astype(int)
    return outdf.sort_values(["rank_battery_coherence", "rank_usable_battery", "rank_ECEE"]).reset_index(drop=True)


def save_plots(df: pd.DataFrame, outdir: Path, prefix: str, top: int) -> None:
    if not HAS_MPL:
        return
    outdir.mkdir(parents=True, exist_ok=True)
    topdf = df.nsmallest(top, "rank_battery_coherence").copy()

    plt.figure(figsize=(12, max(6, 0.38 * len(topdf))))
    y = np.arange(len(topdf))
    plt.barh(y - 0.25, topdf["ECEM_flux"], height=0.24, label="ECEM_flux")
    plt.barh(y, topdf["ECEE_score"], height=0.24, label="ECEE")
    plt.barh(y + 0.25, topdf["usable_battery_score"], height=0.24, label="Usable battery")
    plt.yticks(y, topdf["name"])
    plt.xlabel("score")
    plt.title("ECEM vs ECEE vs usable planetary battery")
    plt.gca().invert_yaxis()
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_scores.png", dpi=180)
    plt.close()

    plt.figure(figsize=(12, max(6, 0.38 * len(topdf))))
    y = np.arange(len(topdf))
    plt.barh(y - 0.18, topdf["planetary_battery_score"], height=0.34, label="Raw battery capacity")
    plt.barh(y + 0.18, topdf["usable_battery_score"], height=0.34, label="Usable rocky-boundary battery")
    plt.yticks(y, topdf["name"])
    plt.xlabel("score")
    plt.title("Raw battery vs usable rocky-boundary battery")
    plt.gca().invert_yaxis()
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "ecee_usable_battery_scores.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.scatter(df["core_to_stellar_ratio_relative"], df["usable_battery_score"], s=80)
    for _, r in topdf.iterrows():
        plt.annotate(str(r["name"]), (r["core_to_stellar_ratio_relative"], r["usable_battery_score"]), fontsize=8, xytext=(4,4), textcoords="offset points")
    plt.xscale("log")
    plt.xlabel("core-to-stellar flux ratio / Earth")
    plt.ylabel("usable_battery_score")
    plt.title("Internal flux can compensate for lower stellar flux")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / "ecee_core_to_stellar.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.scatter(df["ECEE_score"], df["usable_battery_score"], s=80)
    for _, r in topdf.iterrows():
        plt.annotate(str(r["name"]), (r["ECEE_score"], r["usable_battery_score"]), fontsize=8, xytext=(4,4), textcoords="offset points")
    plt.xlabel("ECEE_score")
    plt.ylabel("usable_battery_score")
    plt.title("Surface equivalence vs internal maintenance engine")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / "ecee_battery_vs_ecee.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(df["planet_radius_earth"], df["planet_mass_earth"], c=df["usable_battery_score"], s=90)
    plt.colorbar(sc, label="usable_battery_score")
    for _, r in topdf.iterrows():
        plt.annotate(str(r["name"]), (r["planet_radius_earth"], r["planet_mass_earth"]), fontsize=8, xytext=(4,4), textcoords="offset points")
    plt.xlabel("planet radius / Earth")
    plt.ylabel("planet mass / Earth")
    plt.title("Usable rocky/super-Earth batteries vs wrong-boundary controls")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / "ecee_mass_radius_battery.png", dpi=180)
    plt.close()


def print_summary(df: pd.DataFrame, top: int) -> None:
    cols = [
        "rank_battery_coherence", "rank_usable_battery", "rank_battery", "rank_ECEE", "rank_ECEM_flux",
        "name", "star_class", "battery_coherence_score", "usable_battery_score", "planetary_battery_score", "rocky_boundary_gate",
        "ECEE_score", "ECEM_flux", "core_heat_relative", "core_to_stellar_ratio_relative",
        "battery_lifetime_proxy", "retention_raw", "score_ecee_rocky_boundary",
        "score_surface_state", "branch_classification",
    ]
    view = df[cols].head(top).copy()
    with pd.option_context("display.max_rows", top, "display.max_columns", None, "display.width", 260):
        print("\nECEE × Planetary Battery ranking")
        print("=" * 96)
        print(view.to_string(index=False, float_format=lambda x: f"{x:.6f}"))


def main() -> None:
    parser = argparse.ArgumentParser(description="ECEM/ECEE planetary-battery ranker for flux-coherence exoplanet screening.")
    parser.add_argument("--input", type=str, default=None, help="Optional input CSV. If omitted, uses built-in candidate/control catalog.")
    parser.add_argument("--output", type=str, default="ecee_planetary_battery_ranked.csv", help="Output ranked CSV path.")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory for plots.")
    parser.add_argument("--top", type=int, default=25, help="Number of rows to print/annotate.")
    parser.add_argument("--prefix", type=str, default="ecee_planetary_battery", help="Prefix for the main score plot.")
    args = parser.parse_args()

    if args.input:
        inpath = Path(args.input)
        if not inpath.exists():
            raise FileNotFoundError(f"Input CSV not found: {inpath.resolve()}")
        df = pd.read_csv(inpath)
    else:
        df = built_in_catalog()

    scored = score_dataframe(df)
    outpath = Path(args.output)
    scored.to_csv(outpath, index=False)
    save_plots(scored, Path(args.outdir), args.prefix, args.top)
    print_summary(scored, args.top)
    print(f"\nWrote: {outpath.resolve()}")
    if HAS_MPL:
        print(f"Wrote plots to: {Path(args.outdir).resolve()}")


if __name__ == "__main__":
    main()
