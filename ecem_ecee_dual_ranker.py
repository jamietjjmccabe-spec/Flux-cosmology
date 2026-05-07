#!/usr/bin/env python3
"""
ecem_ecee_dual_ranker.py

Dual-ranker for Flux Cosmology exoplanet screening.

ECEM  = Earth Coherence Envelope Match
        Similarity metric: how closely the system resembles Earth's measured solution.

ECEE  = Earth Coherence Envelope Equivalence
        Functional metric: whether different physical inputs can reproduce an
        Earth-equivalent surface/coherence envelope.

This file is deliberately standalone:
- runs with a built-in demonstration catalog
- accepts user/catalog CSV files
- writes ranked CSV + plots

Example:
    python ecem_ecee_dual_ranker.py
    python ecem_ecee_dual_ranker.py --input exoplanets.csv --top 50

Expected CSV columns are flexible. Common aliases are supported, including
NASA Exoplanet Archive-like names:
    pl_name, hostname, st_spectype, st_mass, st_rad, st_teff, st_lum,
    pl_bmasse, pl_rade, pl_orbsmax, pl_orbper, pl_orbeccen,
    pl_insol, pl_eqt, pl_dens

Missing values are inferred conservatively where possible.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False


# -----------------------------------------------------------------------------
# Constants / reference values
# -----------------------------------------------------------------------------

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
}

STAR_CLASS_COHERENCE = {
    "O": 0.10,
    "B": 0.18,
    "A": 0.38,
    "F": 0.65,
    "G": 1.00,
    "K": 1.05,
    "M": 0.70,
    "UNKNOWN": 0.72,
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


# -----------------------------------------------------------------------------
# Utility math
# -----------------------------------------------------------------------------

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
    """Gaussian in log-ratio space. Good for scale quantities."""
    if not np.isfinite(value) or not np.isfinite(target) or value <= 0 or target <= 0:
        return 0.0
    return clamp01(np.exp(-((np.log(value / target) / width) ** 2)))


def linear_gauss(value: float, target: float, width: float) -> float:
    if not np.isfinite(value):
        return 0.0
    return clamp01(np.exp(-(((value - target) / width) ** 2)))


def weighted_product(scores: Dict[str, float], weights: Dict[str, float]) -> float:
    """Weighted product in log-space. Prevents one catastrophic failure being hidden."""
    total_w = sum(weights.values())
    if total_w <= 0:
        return 0.0
    acc = 0.0
    for key, w in weights.items():
        s = clamp01(scores.get(key, 0.0))
        acc += (w / total_w) * safe_log(max(s, 1e-6))
    return clamp01(np.exp(acc))


def weighted_sum(scores: Dict[str, float], weights: Dict[str, float]) -> float:
    total_w = sum(weights.values())
    if total_w <= 0:
        return 0.0
    return clamp01(sum(clamp01(scores.get(k, 0.0)) * w for k, w in weights.items()) / total_w)


# -----------------------------------------------------------------------------
# Input normalization and inference
# -----------------------------------------------------------------------------

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
        if teff >= 30000:
            return "O"
        if teff >= 10000:
            return "B"
        if teff >= 7500:
            return "A"
        if teff >= 6000:
            return "F"
        if teff >= 5200:
            return "G"
        if teff >= 3700:
            return "K"
        return "M"
    return "UNKNOWN"


def mass_from_radius_rocky(radius: float) -> float:
    """Conservative mass-radius estimate in Earth units."""
    if not np.isfinite(radius) or radius <= 0:
        return np.nan
    if radius <= 1.5:
        return radius ** 3.7
    if radius <= 2.0:
        return 2.7 * radius ** 1.3
    # larger planets become ambiguous; don't let mass explode blindly
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
        out["volatile_flag"] = get_alias(row, "volatile_flag", "")

        # Inferences.
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
            # Conservative class baseline. M dwarfs can be flare-active.
            out["stellar_activity_score"] = {"F": 0.70, "G": 0.95, "K": 0.98, "M": 0.55}.get(out["star_class"], 0.70)
        if not np.isfinite(out["atmosphere_score"]):
            out["atmosphere_score"] = infer_atmosphere_score(out)
        if not np.isfinite(out["magnetosphere_score"]):
            out["magnetosphere_score"] = infer_magnetosphere_score(out)

        rows.append(out)
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Proxy inference
# -----------------------------------------------------------------------------

def tidal_strength_earth(row: Dict[str, Any] | pd.Series) -> float:
    mstar = finite_or_nan(row.get("star_mass_solar", np.nan))
    a = finite_or_nan(row.get("semi_major_axis_au", np.nan))
    if not np.isfinite(mstar) or mstar <= 0:
        mstar = 1.0
    if not np.isfinite(a) or a <= 0:
        return np.nan
    return mstar / max(a ** 3, EPS)


def infer_atmosphere_score(row: Dict[str, Any]) -> float:
    m = finite_or_nan(row.get("planet_mass_earth"))
    r = finite_or_nan(row.get("planet_radius_earth"))
    rho = finite_or_nan(row.get("density_earth"))
    insol = finite_or_nan(row.get("insolation_earth"))
    if not np.isfinite(m):
        m = 1.0
    if not np.isfinite(r):
        r = 1.0
    if not np.isfinite(rho):
        rho = m / max(r ** 3, EPS)
    if not np.isfinite(insol):
        insol = 1.0

    retention = log_gauss_ratio(max(m / max(r, EPS), EPS), target=1.0, width=1.5)
    rocky = 1.0 if rho >= 0.65 and r <= 2.0 else max(0.2, 1.0 - 0.35 * max(r - 1.8, 0.0))
    runaway = 1.0
    if insol > 1.35:
        runaway *= math.exp(-1.1 * (insol - 1.35) ** 2)
    if insol < 0.25:
        runaway *= math.exp(-1.4 * (0.25 - insol) ** 2)
    return clamp01(0.55 * retention + 0.25 * rocky + 0.20 * runaway)


def infer_magnetosphere_score(row: Dict[str, Any]) -> float:
    m = finite_or_nan(row.get("planet_mass_earth"))
    r = finite_or_nan(row.get("planet_radius_earth"))
    star_class = str(row.get("star_class", "UNKNOWN")).upper()
    a = finite_or_nan(row.get("semi_major_axis_au"))
    if not np.isfinite(m):
        m = 1.0
    if not np.isfinite(r):
        r = 1.0
    dyn = log_gauss_ratio(max(m / max(r ** 2, EPS), EPS), 1.0, 1.4)
    if m >= 1.5 and m <= 6.0 and r <= 2.0:
        dyn = max(dyn, 0.72)
    tidal = tidal_strength_earth(row)
    tidal_penalty = 1.0 if not np.isfinite(tidal) else math.exp(-0.13 * max(np.log10(max(tidal, 1.0)), 0.0) ** 2)
    class_penalty = 0.82 if star_class == "M" else 1.0
    if np.isfinite(a) and a > 0.65 and star_class in ["G", "K"]:
        class_penalty *= 1.03
    return clamp01(dyn * tidal_penalty * class_penalty)


# -----------------------------------------------------------------------------
# ECEM scoring: Earth similarity
# -----------------------------------------------------------------------------

def star_score(row: pd.Series) -> float:
    c = str(row.get("star_class", "UNKNOWN")).upper()
    class_score = STAR_CLASS_COHERENCE.get(c, STAR_CLASS_COHERENCE["UNKNOWN"])
    mass_score = log_gauss_ratio(row.get("star_mass_solar", np.nan), 1.0, 0.45)
    lum_score = log_gauss_ratio(row.get("star_luminosity_solar", np.nan), 1.0, 0.75)
    teff_score = linear_gauss(row.get("star_teff_k", np.nan), 5772.0, 1200.0)
    return weighted_sum(
        {"class": class_score, "mass": mass_score, "lum": lum_score, "teff": teff_score},
        {"class": 0.35, "mass": 0.25, "lum": 0.25, "teff": 0.15},
    )


def orbit_score(row: pd.Series) -> float:
    a_score = log_gauss_ratio(row.get("semi_major_axis_au", np.nan), 1.0, 0.8)
    p_score = log_gauss_ratio(row.get("orbital_period_days", np.nan), 365.25, 1.0)
    insol_score = log_gauss_ratio(row.get("insolation_earth", np.nan), 1.0, 0.6)
    ecc = finite_or_nan(row.get("eccentricity", 0.05))
    ecc_score = math.exp(-((ecc / 0.25) ** 2)) if np.isfinite(ecc) else 0.8
    return weighted_sum(
        {"a": a_score, "period": p_score, "insol": insol_score, "ecc": ecc_score},
        {"a": 0.25, "period": 0.15, "insol": 0.45, "ecc": 0.15},
    )


def planet_similarity_score(row: pd.Series) -> float:
    mass_score = log_gauss_ratio(row.get("planet_mass_earth", np.nan), 1.0, 0.9)
    radius_score = log_gauss_ratio(row.get("planet_radius_earth", np.nan), 1.0, 0.55)
    density_score = log_gauss_ratio(row.get("density_earth", np.nan), 1.0, 0.8)
    gravity_score = log_gauss_ratio(row.get("gravity_earth", np.nan), 1.0, 0.75)
    return weighted_sum(
        {"mass": mass_score, "radius": radius_score, "density": density_score, "gravity": gravity_score},
        {"mass": 0.25, "radius": 0.35, "density": 0.20, "gravity": 0.20},
    )


def thermal_similarity_score(row: pd.Series) -> float:
    insol = finite_or_nan(row.get("insolation_earth", np.nan))
    teq = finite_or_nan(row.get("equilibrium_temp_k", np.nan))
    s1 = log_gauss_ratio(insol, 1.0, 0.55)
    s2 = linear_gauss(teq, 255.0, 45.0)
    return weighted_sum({"insol": s1, "teq": s2}, {"insol": 0.7, "teq": 0.3})


def tidal_score_ecem(row: pd.Series) -> float:
    t = tidal_strength_earth(row)
    if not np.isfinite(t) or t <= 0:
        return 0.5
    return log_gauss_ratio(t, 1.0, 1.0)


def sigma_proxy(row: pd.Series) -> float:
    """Toy residual flux/coherence scalar proxy centered on Earth."""
    insol = finite_or_nan(row.get("insolation_earth", np.nan))
    tidal = tidal_strength_earth(row)
    activity = finite_or_nan(row.get("stellar_activity_score", 0.75))
    retention = finite_or_nan(row.get("atmosphere_score", 0.6)) * finite_or_nan(row.get("magnetosphere_score", 0.6))
    if not np.isfinite(insol):
        insol = 1.0
    if not np.isfinite(tidal):
        tidal = 1.0
    if not np.isfinite(activity):
        activity = 0.75
    return float(max(insol, EPS) * max(activity, EPS) / max(tidal ** 0.08, EPS) * max(retention, EPS) ** 0.2)


def score_sigma(row: pd.Series) -> float:
    return log_gauss_ratio(sigma_proxy(row), 1.0, 0.85)


def venus_branch_penalty(row: pd.Series) -> float:
    insol = finite_or_nan(row.get("insolation_earth", np.nan))
    atm = finite_or_nan(row.get("atmosphere_score", 0.5))
    mag = finite_or_nan(row.get("magnetosphere_score", 0.5))
    rot = finite_or_nan(row.get("rotation_period_days", np.nan))
    penalty = 1.0
    if np.isfinite(insol) and insol > 1.35:
        penalty *= np.exp(-2.5 * (insol - 1.35) ** 2)
    if np.isfinite(insol) and insol > 1.1 and atm < 0.45:
        penalty *= 0.55
    if np.isfinite(insol) and insol > 1.1 and mag < 0.45:
        penalty *= 0.70
    if np.isfinite(rot) and rot > 100:
        penalty *= 0.75
    return clamp01(penalty)


def surface_state_score(row: pd.Series) -> float:
    insol = finite_or_nan(row.get("insolation_earth", np.nan))
    atm = finite_or_nan(row.get("atmosphere_score", 0.6))
    mag = finite_or_nan(row.get("magnetosphere_score", 0.6))
    if not np.isfinite(insol):
        insol = 1.0
    thermal_window = log_gauss_ratio(insol, 1.0, 0.80)
    if insol > 1.7:
        thermal_window *= 0.35
    if insol < 0.25:
        thermal_window *= 0.45
    return clamp01(thermal_window * venus_branch_penalty(row) * (0.55 + 0.25 * atm + 0.20 * mag))


def flux_tension_proxy(row: pd.Series) -> float:
    insol = finite_or_nan(row.get("insolation_earth", np.nan))
    activity = finite_or_nan(row.get("stellar_activity_score", np.nan))
    tidal = tidal_strength_earth(row)
    age = finite_or_nan(row.get("star_age_gyr", np.nan))
    c = str(row.get("star_class", "UNKNOWN")).upper()
    if not np.isfinite(insol):
        insol = 1.0
    if not np.isfinite(activity):
        activity = {"F": 1.2, "G": 1.0, "K": 0.85, "M": 1.4}.get(c, 1.0)
    if not np.isfinite(tidal):
        tidal = 1.0
    if not np.isfinite(age):
        age = 4.57
    age_factor = 1.0 + 0.05 * np.log1p(max(age, 0.0))
    return float(max(insol, EPS) * max(activity, EPS) * max(tidal, EPS) ** 0.12 * age_factor)


def flux_history_score(row: pd.Series) -> float:
    tension = flux_tension_proxy(row)
    age = finite_or_nan(row.get("star_age_gyr", 4.57))
    c = str(row.get("star_class", "UNKNOWN")).upper()
    base = log_gauss_ratio(tension, 1.0, 1.15)
    class_bonus = {"G": 1.0, "K": 1.05, "M": 0.86, "F": 0.82}.get(c, 0.85)
    age_score = log_gauss_ratio(max(age, 0.1), 4.57, 1.2)
    # K/M stars can remain coherent for longer; don't punish older age too hard.
    if c in ["K", "M"] and age > 4.57:
        age_score = max(age_score, 0.82)
    return clamp01(base * class_bonus * (0.7 + 0.3 * age_score))


def compute_ecem_flux(row: pd.Series) -> Dict[str, float]:
    scores = {
        "star": star_score(row),
        "orbit": orbit_score(row),
        "planet": planet_similarity_score(row),
        "thermal": thermal_similarity_score(row),
        "atmosphere": finite_or_nan(row.get("atmosphere_score", 0.5)),
        "magnetosphere": finite_or_nan(row.get("magnetosphere_score", 0.5)),
        "tidal": tidal_score_ecem(row),
        "sigma": score_sigma(row),
    }
    base = weighted_sum(scores, ECEM_WEIGHTS)
    surf = surface_state_score(row)
    history = flux_history_score(row)
    tension_score = log_gauss_ratio(flux_tension_proxy(row), 1.0, 1.0)
    adjustment = weighted_product(
        {"surface": surf, "history": history, "tension": tension_score},
        {"surface": 0.40, "history": 0.35, "tension": 0.25},
    )
    return {
        **{f"score_ecem_{k}": clamp01(v) for k, v in scores.items()},
        "ECEM_base": clamp01(base),
        "score_surface_state": surf,
        "score_flux_history": history,
        "score_flux_tension": tension_score,
        "flux_adjustment": adjustment,
        "ECEM_flux": clamp01(base * adjustment),
        "sigma_proxy": sigma_proxy(row),
        "flux_tension_proxy": flux_tension_proxy(row),
        "tidal_strength_earth": tidal_strength_earth(row),
    }


# -----------------------------------------------------------------------------
# ECEE scoring: functional equivalence
# -----------------------------------------------------------------------------

def retention_raw(row: pd.Series) -> float:
    m = finite_or_nan(row.get("planet_mass_earth", np.nan))
    g = finite_or_nan(row.get("gravity_earth", np.nan))
    atm = finite_or_nan(row.get("atmosphere_score", np.nan))
    mag = finite_or_nan(row.get("magnetosphere_score", np.nan))
    insol = finite_or_nan(row.get("insolation_earth", np.nan))
    tide = tidal_strength_earth(row)
    if not np.isfinite(m):
        m = 1.0
    if not np.isfinite(g):
        g = 1.0
    if not np.isfinite(atm):
        atm = 0.55
    if not np.isfinite(mag):
        mag = 0.55
    if not np.isfinite(insol):
        insol = 1.0
    if not np.isfinite(tide):
        tide = 1.0
    return float(max(m * g * atm * mag, EPS) / max(insol * (tide ** 0.15), EPS))


def retention_capacity_score(row: pd.Series) -> float:
    rr = retention_raw(row)
    earth_rr = 1.0
    if rr >= earth_rr:
        # Extra retention is allowed for outer super-Earths unless extreme.
        return log_gauss_ratio(rr, earth_rr, 2.5)
    # Too little retention fails quickly.
    return log_gauss_ratio(rr, earth_rr, 1.0)


def rocky_boundary_score(row: pd.Series) -> float:
    r = finite_or_nan(row.get("planet_radius_earth", np.nan))
    m = finite_or_nan(row.get("planet_mass_earth", np.nan))
    rho = finite_or_nan(row.get("density_earth", np.nan))
    volatile_flag = str(row.get("volatile_flag", "")).strip().lower()

    if not np.isfinite(r):
        r = 1.5
    if not np.isfinite(m):
        m = mass_from_radius_rocky(r)
    if not np.isfinite(rho):
        rho = m / max(r ** 3, EPS)

    score = 1.0
    # ECEE allows larger rocky worlds, but not volatile/sub-Neptune boundary layers.
    if r > 2.0:
        score *= np.exp(-1.5 * (r - 2.0) ** 2)
    if r > 2.5:
        score *= 0.35
    if rho < 0.65:
        score *= np.exp(-2.0 * (0.65 - rho) ** 2)
    if rho < 0.45:
        score *= 0.35
    if m > 10.0:
        score *= 0.55
    if any(x in volatile_flag for x in ["true", "yes", "sub", "volatile", "neptune"]):
        score *= 0.30
    # Sweet spot bonus for dense outer super-Earths.
    if 1.25 <= r <= 2.0 and 2.0 <= m <= 8.0 and rho >= 0.65:
        score *= 1.08
    return clamp01(score)


def thermal_balance_score_ecee(row: pd.Series) -> float:
    insol = finite_or_nan(row.get("insolation_earth", np.nan))
    if not np.isfinite(insol):
        insol = 1.0
    rr = retention_raw(row)
    # Greenhouse/retention compensation. Saturated to avoid making gas dwarfs win.
    greenhouse_retention = 1.0 + 0.35 * np.log1p(min(rr, 12.0))
    thermal_envelope = insol * greenhouse_retention
    # ECEE accepts cooler outer-HZ bodies if retention compensates.
    score = log_gauss_ratio(thermal_envelope, 1.0, 0.75)
    # Hard runaway/freeze guards.
    if insol > 1.8:
        score *= np.exp(-1.5 * (insol - 1.8) ** 2)
    if insol < 0.18:
        score *= np.exp(-2.0 * (0.18 - insol) ** 2)
    return clamp01(score)


def magnetosphere_dynamo_score(row: pd.Series) -> float:
    mag = finite_or_nan(row.get("magnetosphere_score", np.nan))
    if not np.isfinite(mag):
        mag = infer_magnetosphere_score(row.to_dict())
    m = finite_or_nan(row.get("planet_mass_earth", 1.0))
    r = finite_or_nan(row.get("planet_radius_earth", 1.0))
    g = finite_or_nan(row.get("gravity_earth", 1.0))
    if not np.isfinite(m):
        m = 1.0
    if not np.isfinite(r):
        r = 1.0
    if not np.isfinite(g):
        g = m / max(r ** 2, EPS)
    dyn_capacity = log_gauss_ratio(max(m * g / max(r, EPS), EPS), 1.0, 2.0)
    if 2 <= m <= 8 and r <= 2.0:
        dyn_capacity = max(dyn_capacity, 0.75)
    return clamp01(0.55 * mag + 0.45 * dyn_capacity)


def tidal_stability_score_ecee(row: pd.Series) -> float:
    tide = tidal_strength_earth(row)
    c = str(row.get("star_class", "UNKNOWN")).upper()
    if not np.isfinite(tide) or tide <= 0:
        return 0.55
    # ECEE does not require Earth-level tides, but extreme compact systems are distinct branches.
    if tide <= 5.0:
        return clamp01(log_gauss_ratio(tide, 1.0, 2.0) * 1.05)
    penalty = math.exp(-0.18 * (np.log10(tide / 5.0) ** 2))
    if c == "M":
        penalty *= 0.82
    return clamp01(penalty)


def stellar_coherence_score_ecee(row: pd.Series) -> float:
    c = str(row.get("star_class", "UNKNOWN")).upper()
    class_score = {"F": 0.72, "G": 1.00, "K": 1.08, "M": 0.72}.get(c, 0.72)
    activity = finite_or_nan(row.get("stellar_activity_score", np.nan))
    if not np.isfinite(activity):
        activity = {"F": 0.70, "G": 0.95, "K": 0.98, "M": 0.55}.get(c, 0.70)
    age = finite_or_nan(row.get("star_age_gyr", 4.57))
    age_score = log_gauss_ratio(max(age, 0.1), 4.57, 1.4)
    if c == "K" and age > 4:
        age_score = max(age_score, 0.86)
    if c == "M" and age > 2:
        age_score = max(age_score, 0.72)
    return clamp01(class_score * (0.65 + 0.35 * activity) * (0.75 + 0.25 * age_score))


def compute_ecee(row: pd.Series) -> Dict[str, float]:
    scores = {
        "thermal_balance": thermal_balance_score_ecee(row),
        "retention_capacity": retention_capacity_score(row),
        "rocky_boundary": rocky_boundary_score(row),
        "magnetosphere_dynamo": magnetosphere_dynamo_score(row),
        "tidal_stability": tidal_stability_score_ecee(row),
        "stellar_coherence": stellar_coherence_score_ecee(row),
        "flux_history": flux_history_score(row),
    }
    ecee = weighted_product(scores, ECEE_WEIGHTS)
    return {
        **{f"score_ecee_{k}": clamp01(v) for k, v in scores.items()},
        "retention_raw": retention_raw(row),
        "thermal_envelope_proxy": finite_or_nan(row.get("insolation_earth", 1.0)) * (1.0 + 0.35 * np.log1p(min(retention_raw(row), 12.0))),
        "ECEE_score": ecee,
    }


# -----------------------------------------------------------------------------
# Classification
# -----------------------------------------------------------------------------

def classify(row: pd.Series) -> str:
    name = str(row.get("name", "")).lower()
    c = str(row.get("star_class", "UNKNOWN")).upper()
    r = finite_or_nan(row.get("planet_radius_earth", np.nan))
    rho = finite_or_nan(row.get("density_earth", np.nan))
    m = finite_or_nan(row.get("planet_mass_earth", np.nan))
    insol = finite_or_nan(row.get("insolation_earth", np.nan))
    tide = finite_or_nan(row.get("tidal_strength_earth", tidal_strength_earth(row)))
    ecem = finite_or_nan(row.get("ECEM_flux", 0.0))
    ecee = finite_or_nan(row.get("ECEE_score", 0.0))
    rocky = finite_or_nan(row.get("score_ecee_rocky_boundary", 0.0))
    surface = finite_or_nan(row.get("score_surface_state", 0.0))

    if "earth" == name:
        return "reference Earth solution"
    if r > 2.2 and (not np.isfinite(rho) or rho < 0.65):
        return "sub-Neptune / volatile boundary control"
    if rocky < 0.35:
        return "wrong boundary layer / volatile-rich branch"
    if np.isfinite(insol) and insol > 1.25 and surface < 0.35:
        return "Venus/runaway failed coherence branch"
    if np.isfinite(insol) and insol < 0.45 and ecee < 0.45 and (m < 1.2 if np.isfinite(m) else True):
        return "Mars/freeze-out weak-envelope branch"
    if ecem >= 0.68 and ecee >= 0.68:
        return "high Earth-envelope analog"
    if c in ["G", "K"] and ecee >= 0.62 and ecem < 0.68:
        return "outer super-Earth envelope candidate"
    if c == "K" and ecee >= 0.55:
        return "K-dwarf long-coherence envelope candidate"
    if c == "M" and np.isfinite(tide) and tide > 50:
        if ecee >= 0.48 and m >= 2.0 and rocky >= 0.65:
            return "dense M-dwarf retained-envelope branch"
        return "compact M-dwarf high-tidal branch"
    if ecee >= 0.55 and ecem < 0.50:
        return "non-Earth equivalent coherence branch"
    if ecem >= 0.55:
        return "candidate Earth-like coherence analog"
    if ecee >= 0.45:
        return "partial equivalence candidate / needs data"
    return "low coherence-envelope match"


# -----------------------------------------------------------------------------
# Built-in demo catalog
# -----------------------------------------------------------------------------

def demo_catalog() -> pd.DataFrame:
    data = [
        EARTH,
        {"name":"Venus", "star_class":"G", "star_mass_solar":1.0, "star_radius_solar":1.0, "star_teff_k":5772, "star_luminosity_solar":1.0, "star_age_gyr":4.57, "planet_mass_earth":0.815, "planet_radius_earth":0.949, "semi_major_axis_au":0.723, "orbital_period_days":224.7, "eccentricity":0.0068, "insolation_earth":1.91, "equilibrium_temp_k":328, "density_earth":0.95, "gravity_earth":0.904, "atmosphere_score":0.08, "magnetosphere_score":0.10, "stellar_activity_score":1.0, "rotation_period_days":243},
        {"name":"Mars", "star_class":"G", "star_mass_solar":1.0, "star_radius_solar":1.0, "star_teff_k":5772, "star_luminosity_solar":1.0, "star_age_gyr":4.57, "planet_mass_earth":0.107, "planet_radius_earth":0.532, "semi_major_axis_au":1.524, "orbital_period_days":687, "eccentricity":0.093, "insolation_earth":0.43, "equilibrium_temp_k":210, "density_earth":0.71, "gravity_earth":0.379, "atmosphere_score":0.18, "magnetosphere_score":0.12, "stellar_activity_score":1.0, "rotation_period_days":1.03},
        {"name":"Kepler-452 b", "star_class":"G", "star_mass_solar":1.04, "star_radius_solar":1.11, "star_teff_k":5757, "star_luminosity_solar":1.20, "star_age_gyr":6.0, "planet_mass_earth":5.0, "planet_radius_earth":1.63, "semi_major_axis_au":1.046, "orbital_period_days":384.8, "eccentricity":0.05, "insolation_earth":1.10, "equilibrium_temp_k":265, "density_earth":1.15, "gravity_earth":1.88, "atmosphere_score":0.72, "magnetosphere_score":0.75, "stellar_activity_score":0.92},
        {"name":"Kepler-442 b", "star_class":"K", "star_mass_solar":0.61, "star_radius_solar":0.60, "star_teff_k":4402, "star_luminosity_solar":0.12, "star_age_gyr":5.0, "planet_mass_earth":2.36, "planet_radius_earth":1.34, "semi_major_axis_au":0.409, "orbital_period_days":112.3, "eccentricity":0.04, "insolation_earth":0.70, "equilibrium_temp_k":233, "density_earth":0.98, "gravity_earth":1.31, "atmosphere_score":0.72, "magnetosphere_score":0.70, "stellar_activity_score":0.96},
        {"name":"Kepler-62 e", "star_class":"K", "star_mass_solar":0.69, "star_radius_solar":0.64, "star_teff_k":4925, "star_luminosity_solar":0.21, "star_age_gyr":7.0, "planet_mass_earth":4.5, "planet_radius_earth":1.61, "semi_major_axis_au":0.427, "orbital_period_days":122.4, "eccentricity":0.05, "insolation_earth":1.20, "equilibrium_temp_k":267, "density_earth":1.08, "gravity_earth":1.74, "atmosphere_score":0.68, "magnetosphere_score":0.73, "stellar_activity_score":0.95},
        {"name":"Kepler-62 f", "star_class":"K", "star_mass_solar":0.69, "star_radius_solar":0.64, "star_teff_k":4925, "star_luminosity_solar":0.21, "star_age_gyr":7.0, "planet_mass_earth":3.0, "planet_radius_earth":1.41, "semi_major_axis_au":0.718, "orbital_period_days":267.3, "eccentricity":0.05, "insolation_earth":0.41, "equilibrium_temp_k":204, "density_earth":1.07, "gravity_earth":1.51, "atmosphere_score":0.66, "magnetosphere_score":0.72, "stellar_activity_score":0.95},
        {"name":"Outer Super-Earth Test", "star_class":"K", "star_mass_solar":0.75, "star_radius_solar":0.75, "star_teff_k":5000, "star_luminosity_solar":0.35, "star_age_gyr":6.0, "planet_mass_earth":4.5, "planet_radius_earth":1.60, "semi_major_axis_au":0.85, "orbital_period_days":330, "eccentricity":0.05, "insolation_earth":0.48, "equilibrium_temp_k":212, "density_earth":1.10, "gravity_earth":1.76, "atmosphere_score":0.75, "magnetosphere_score":0.80, "stellar_activity_score":0.97},
        {"name":"LHS 1140 b", "star_class":"M", "star_mass_solar":0.18, "star_radius_solar":0.21, "star_teff_k":3216, "star_luminosity_solar":0.0044, "star_age_gyr":5.0, "planet_mass_earth":5.6, "planet_radius_earth":1.73, "semi_major_axis_au":0.094, "orbital_period_days":24.7, "eccentricity":0.04, "insolation_earth":0.50, "equilibrium_temp_k":230, "density_earth":1.08, "gravity_earth":1.87, "atmosphere_score":0.68, "magnetosphere_score":0.55, "stellar_activity_score":0.62},
        {"name":"TOI-700 d", "star_class":"M", "star_mass_solar":0.42, "star_radius_solar":0.42, "star_teff_k":3480, "star_luminosity_solar":0.023, "star_age_gyr":2.0, "planet_mass_earth":1.7, "planet_radius_earth":1.19, "semi_major_axis_au":0.163, "orbital_period_days":37.4, "eccentricity":0.05, "insolation_earth":0.86, "equilibrium_temp_k":246, "density_earth":1.01, "gravity_earth":1.20, "atmosphere_score":0.72, "magnetosphere_score":0.50, "stellar_activity_score":0.68},
        {"name":"TRAPPIST-1 e", "star_class":"M", "star_mass_solar":0.089, "star_radius_solar":0.12, "star_teff_k":2559, "star_luminosity_solar":0.00055, "star_age_gyr":7.6, "planet_mass_earth":0.77, "planet_radius_earth":0.92, "semi_major_axis_au":0.029, "orbital_period_days":6.1, "eccentricity":0.006, "insolation_earth":0.66, "equilibrium_temp_k":230, "density_earth":0.99, "gravity_earth":0.91, "atmosphere_score":0.55, "magnetosphere_score":0.35, "stellar_activity_score":0.38},
        {"name":"Proxima Centauri b", "star_class":"M", "star_mass_solar":0.122, "star_radius_solar":0.154, "star_teff_k":3042, "star_luminosity_solar":0.00155, "star_age_gyr":4.8, "planet_mass_earth":1.27, "planet_radius_earth":1.08, "semi_major_axis_au":0.0485, "orbital_period_days":11.2, "eccentricity":0.11, "insolation_earth":0.65, "equilibrium_temp_k":234, "density_earth":1.00, "gravity_earth":1.09, "atmosphere_score":0.48, "magnetosphere_score":0.30, "stellar_activity_score":0.28},
        {"name":"K2-18 b", "star_class":"M", "star_mass_solar":0.50, "star_radius_solar":0.47, "star_teff_k":3457, "star_luminosity_solar":0.023, "star_age_gyr":2.4, "planet_mass_earth":8.6, "planet_radius_earth":2.61, "semi_major_axis_au":0.142, "orbital_period_days":33, "eccentricity":0.05, "insolation_earth":1.0, "equilibrium_temp_k":255, "density_earth":0.48, "gravity_earth":1.26, "atmosphere_score":0.40, "magnetosphere_score":0.45, "stellar_activity_score":0.55, "volatile_flag":"sub-Neptune"},
    ]
    return pd.DataFrame(data)


# -----------------------------------------------------------------------------
# Output / plotting
# -----------------------------------------------------------------------------

def score_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    norm = normalize_catalog(df)
    rows = []
    for _, row in norm.iterrows():
        d = row.to_dict()
        ecem = compute_ecem_flux(row)
        tmp = pd.Series({**d, **ecem})
        ecee = compute_ecee(tmp)
        out = {**d, **ecem, **ecee}
        out_series = pd.Series(out)
        out["dual_coherence_score"] = clamp01(0.45 * out["ECEM_flux"] + 0.55 * out["ECEE_score"])
        out_series = pd.Series(out)
        out["branch_classification"] = classify(out_series)
        rows.append(out)
    outdf = pd.DataFrame(rows)
    outdf["rank_ECEM_flux"] = outdf["ECEM_flux"].rank(ascending=False, method="min").astype(int)
    outdf["rank_ECEE"] = outdf["ECEE_score"].rank(ascending=False, method="min").astype(int)
    outdf["rank_dual"] = outdf["dual_coherence_score"].rank(ascending=False, method="min").astype(int)
    return outdf.sort_values(["rank_dual", "rank_ECEE", "rank_ECEM_flux"]).reset_index(drop=True)


def save_plots(df: pd.DataFrame, outdir: Path, prefix: str, top: int) -> None:
    if not HAS_MPL:
        return
    outdir.mkdir(parents=True, exist_ok=True)
    topdf = df.nsmallest(top, "rank_dual").copy()

    # 1. Dual ranking bar plot.
    plt.figure(figsize=(12, max(6, 0.35 * len(topdf))))
    y = np.arange(len(topdf))
    plt.barh(y - 0.2, topdf["ECEM_flux"], height=0.38, label="ECEM_flux similarity")
    plt.barh(y + 0.2, topdf["ECEE_score"], height=0.38, label="ECEE equivalence")
    plt.yticks(y, topdf["name"])
    plt.xlabel("score")
    plt.title("ECEM similarity vs ECEE functional equivalence")
    plt.gca().invert_yaxis()
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_dual_scores.png", dpi=180)
    plt.close()

    # 2. ECEM vs ECEE scatter.
    plt.figure(figsize=(8, 6))
    plt.scatter(df["ECEM_flux"], df["ECEE_score"], s=70)
    for _, r in topdf.iterrows():
        plt.annotate(str(r["name"]), (r["ECEM_flux"], r["ECEE_score"]), fontsize=8, xytext=(4, 4), textcoords="offset points")
    plt.xlabel("ECEM_flux: Earth resemblance")
    plt.ylabel("ECEE_score: envelope equivalence")
    plt.title("Similarity vs functional equivalence")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_ecem_vs_ecee.png", dpi=180)
    plt.close()

    # 3. Retention vs insolation.
    plt.figure(figsize=(8, 6))
    plt.scatter(df["insolation_earth"], df["retention_raw"], s=70)
    for _, r in topdf.iterrows():
        plt.annotate(str(r["name"]), (r["insolation_earth"], r["retention_raw"]), fontsize=8, xytext=(4, 4), textcoords="offset points")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("insolation / Earth")
    plt.ylabel("retention raw proxy")
    plt.title("ECEE: lower flux can be compensated by higher retention")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_retention_vs_insolation.png", dpi=180)
    plt.close()

    # 4. Rocky boundary / radius-mass.
    plt.figure(figsize=(8, 6))
    sc = plt.scatter(df["planet_radius_earth"], df["planet_mass_earth"], c=df["ECEE_score"], s=80)
    plt.colorbar(sc, label="ECEE_score")
    for _, r in topdf.iterrows():
        plt.annotate(str(r["name"]), (r["planet_radius_earth"], r["planet_mass_earth"]), fontsize=8, xytext=(4, 4), textcoords="offset points")
    plt.xlabel("planet radius / Earth")
    plt.ylabel("planet mass / Earth")
    plt.title("Rocky/super-Earth envelope candidates vs volatile controls")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / f"{prefix}_mass_radius_ecee.png", dpi=180)
    plt.close()


def print_summary(df: pd.DataFrame, top: int) -> None:
    cols = [
        "rank_dual", "rank_ECEM_flux", "rank_ECEE", "name", "star_class", 
        "dual_coherence_score", "ECEM_flux", "ECEE_score", "ECEM_base",
        "flux_tension_proxy", "tidal_strength_earth", "retention_raw", 
        "score_ecee_rocky_boundary", "score_surface_state", "branch_classification",
    ]
    view = df[cols].head(top).copy()
    with pd.option_context("display.max_rows", top, "display.max_columns", None, "display.width", 220):
        print("\nECEM × ECEE dual ranking")
        print("=" * 88)
        print(view.to_string(index=False, float_format=lambda x: f"{x:.6f}"))


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Dual ECEM/ECEE exoplanet coherence-envelope ranker.")
    parser.add_argument("--input", "-i", type=str, default=None, help="Optional input CSV. If omitted, uses built-in demo catalog.")
    parser.add_argument("--output", "-o", type=str, default="ecem_ecee_dual_ranked.csv", help="Output ranked CSV path.")
    parser.add_argument("--outdir", type=str, default=".", help="Directory for plots.")
    parser.add_argument("--prefix", type=str, default="ecem_ecee", help="Output plot filename prefix.")
    parser.add_argument("--top", type=int, default=25, help="Number of rows to print/label.")
    parser.add_argument("--no-plots", action="store_true", help="Disable plot generation.")
    args = parser.parse_args()

    if args.input:
        input_path = Path(args.input)
        if not input_path.exists():
            raise FileNotFoundError(
                f"Input CSV not found: {input_path}\n"
                f"Run without --input to use the built-in demo catalog, or provide a full path."
            )
        df = pd.read_csv(input_path)
    else:
        df = demo_catalog()

    ranked = score_dataframe(df)
    output_path = Path(args.output)
    ranked.to_csv(output_path, index=False)
    print_summary(ranked, min(args.top, len(ranked)))
    print(f"\nSaved ranked CSV: {output_path.resolve()}")

    if not args.no_plots:
        save_plots(ranked, Path(args.outdir), args.prefix, min(args.top, len(ranked)))
        if HAS_MPL:
            print(f"Saved plots to: {Path(args.outdir).resolve()}")
        else:
            print("matplotlib unavailable; skipped plots.")


if __name__ == "__main__":
    main()
