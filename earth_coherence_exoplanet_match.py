#!/usr/bin/env python3
"""
earth_coherence_exoplanet_match.py

Flux Cosmology / Earth Coherence Envelope exoplanet ranking tool.

Purpose
-------
This script treats Earth not merely as a planet in a habitable zone,
but as a reference coherence-envelope state: a coupled star--planet--orbit--
atmosphere--magnetosphere--tidal stability system.

It computes an Earth Coherence Envelope Match (ECEM) score for exoplanets.

Core idea
---------
A planet is Earth-like in this model if its boundary conditions allow a
similar local actualization/coherence pocket:

    star stability + orbital stability + rocky mass/radius + Earth-like flux
    + atmosphere retention + magnetic/tidal stability + sigma/coherence balance

The score is intentionally modular. Each observable sector returns a value in
[0, 1]. Higher = closer to Earth envelope.

Inputs
------
1. Built-in demonstration catalog if no CSV is provided.
2. Optional CSV file with some or all of these columns:

    name
    planet_mass_earth
    planet_radius_earth
    orbital_period_days
    semi_major_axis_au
    eccentricity
    stellar_mass_solar
    stellar_radius_solar
    stellar_teff_k
    stellar_luminosity_solar
    stellar_age_gyr
    equilibrium_temp_k
    insolation_earth
    rotation_period_days
    has_large_moon
    atmosphere_score
    magnetosphere_score

Missing values are handled with neutral penalties rather than hard failure.

Outputs
-------
- Ranked CSV: earth_coherence_ranked.csv
- Plots:
    earth_coherence_scores.png
    earth_coherence_mass_radius.png
    earth_coherence_flux_orbit.png

Usage
-----
    python earth_coherence_exoplanet_match.py
    python earth_coherence_exoplanet_match.py --input exoplanets.csv
    python earth_coherence_exoplanet_match.py --top 25 --no-plots

Notes
-----
This is not a claim that the built-in sample data is a complete or current
exoplanet catalog. It is a scoring harness. For serious use, feed it a current
NASA Exoplanet Archive / TEPCat / confirmed-planet CSV export.
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, asdict
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Earth reference values
# ============================================================

EARTH = {
    "planet_mass_earth": 1.0,
    "planet_radius_earth": 1.0,
    "density_earth": 1.0,
    "surface_gravity_earth": 1.0,
    "semi_major_axis_au": 1.0,
    "eccentricity": 0.0167,
    "stellar_mass_solar": 1.0,
    "stellar_radius_solar": 1.0,
    "stellar_teff_k": 5772.0,
    "stellar_luminosity_solar": 1.0,
    "stellar_age_gyr": 4.57,
    "insolation_earth": 1.0,
    "equilibrium_temp_k": 255.0,
    "rotation_period_days": 1.0,
    "has_large_moon": 1.0,
    "atmosphere_score": 1.0,
    "magnetosphere_score": 1.0,
}


# ============================================================
# Model configuration
# ============================================================

@dataclass
class ScoreWeights:
    """Sector weights for the final ECEM score."""

    star: float = 0.16
    orbit: float = 0.20
    planet: float = 0.22
    thermal: float = 0.14
    atmosphere: float = 0.10
    magnetosphere: float = 0.08
    tidal: float = 0.05
    sigma: float = 0.05

    def normalized(self) -> "ScoreWeights":
        total = sum(asdict(self).values())
        if total <= 0:
            raise ValueError("Score weights must sum to a positive number.")
        vals = {k: v / total for k, v in asdict(self).items()}
        return ScoreWeights(**vals)


@dataclass
class EnvelopeWidths:
    """Soft tolerance widths around Earth values.

    Widths are deliberately broad. This keeps the script useful for sparse
    exoplanet data and prevents overfitting to Earth as a cosmetic twin.
    """

    mass_log_width: float = 0.55          # log-space width for mass
    radius_log_width: float = 0.35        # log-space width for radius
    density_log_width: float = 0.45       # log-space width for density
    gravity_log_width: float = 0.40       # log-space width for surface g
    insolation_log_width: float = 0.45    # log-space width around 1 S_earth
    sma_log_width: float = 0.70           # orbit distance width
    eccentricity_width: float = 0.18
    stellar_mass_width: float = 0.35
    stellar_teff_width: float = 900.0
    stellar_age_width: float = 4.0
    eq_temp_width: float = 45.0
    rotation_log_width: float = 1.20
    sigma_width: float = 0.35


# ============================================================
# Utility functions
# ============================================================

def safe_float(value, default: float = np.nan) -> float:
    try:
        if value is None:
            return default
        if isinstance(value, str) and value.strip() == "":
            return default
        out = float(value)
        if math.isfinite(out):
            return out
        return default
    except Exception:
        return default


def bounded01(x: float) -> float:
    if not math.isfinite(x):
        return 0.0
    return max(0.0, min(1.0, x))


def gaussian_delta(delta: float, width: float) -> float:
    """Gaussian closeness score for a linear difference."""
    if not math.isfinite(delta) or width <= 0:
        return 0.5
    return math.exp(-0.5 * (delta / width) ** 2)


def gaussian_log_ratio(value: float, reference: float, width: float) -> float:
    """Gaussian closeness score in log-space.

    Useful for mass, radius, flux, semi-major axis, and rotation period where
    ratios matter more than absolute offsets.
    """
    value = safe_float(value)
    reference = safe_float(reference)
    if value <= 0 or reference <= 0 or width <= 0:
        return 0.5
    return math.exp(-0.5 * (math.log(value / reference) / width) ** 2)


def weighted_mean(scores: Iterable[Tuple[float, float]]) -> float:
    """Weighted mean that ignores non-finite scores or non-positive weights."""
    numerator = 0.0
    denominator = 0.0
    for score, weight in scores:
        if math.isfinite(score) and weight > 0:
            numerator += bounded01(score) * weight
            denominator += weight
    if denominator <= 0:
        return 0.5
    return bounded01(numerator / denominator)


def infer_luminosity_from_star(row: pd.Series) -> float:
    """Infer stellar luminosity L/Lsun from radius and temperature if needed.

    L/Lsun ≈ (R/Rsun)^2 (T/5772 K)^4
    """
    lum = safe_float(row.get("stellar_luminosity_solar"))
    if math.isfinite(lum) and lum > 0:
        return lum

    radius = safe_float(row.get("stellar_radius_solar"))
    teff = safe_float(row.get("stellar_teff_k"))
    if radius > 0 and teff > 0:
        return radius ** 2 * (teff / EARTH["stellar_teff_k"]) ** 4

    return np.nan


def infer_insolation(row: pd.Series) -> float:
    """Infer incident stellar flux in Earth units.

    S/Searth ≈ L/Lsun / a_AU^2
    """
    insolation = safe_float(row.get("insolation_earth"))
    if math.isfinite(insolation) and insolation > 0:
        return insolation

    lum = infer_luminosity_from_star(row)
    sma = safe_float(row.get("semi_major_axis_au"))
    if lum > 0 and sma > 0:
        return lum / (sma ** 2)

    return np.nan


def infer_semi_major_axis(row: pd.Series) -> float:
    """Infer semi-major axis from period and stellar mass if absent.

    Kepler approximation in solar/Earth units:
        a_AU ≈ (Mstar * P_year^2)^(1/3)
    """
    sma = safe_float(row.get("semi_major_axis_au"))
    if math.isfinite(sma) and sma > 0:
        return sma

    period_days = safe_float(row.get("orbital_period_days"))
    mstar = safe_float(row.get("stellar_mass_solar"))
    if period_days > 0 and mstar > 0:
        p_year = period_days / 365.25
        return (mstar * p_year ** 2) ** (1.0 / 3.0)

    return np.nan


def infer_density_earth_units(row: pd.Series) -> float:
    """Density relative to Earth from mass/radius."""
    m = safe_float(row.get("planet_mass_earth"))
    r = safe_float(row.get("planet_radius_earth"))
    if m > 0 and r > 0:
        return m / (r ** 3)
    return np.nan


def infer_surface_gravity_earth_units(row: pd.Series) -> float:
    """Surface gravity relative to Earth from mass/radius."""
    m = safe_float(row.get("planet_mass_earth"))
    r = safe_float(row.get("planet_radius_earth"))
    if m > 0 and r > 0:
        return m / (r ** 2)
    return np.nan


def infer_equilibrium_temp(row: pd.Series) -> float:
    """Very rough Earth-albedo equilibrium temperature estimate.

    T_eq ≈ 255 K * (S/Searth)^0.25
    """
    teq = safe_float(row.get("equilibrium_temp_k"))
    if math.isfinite(teq) and teq > 0:
        return teq

    ins = infer_insolation(row)
    if ins > 0:
        return EARTH["equilibrium_temp_k"] * ins ** 0.25

    return np.nan


def neutral_if_missing(value: float, missing_score: float = 0.55) -> float:
    """Use a cautious neutral score where observations are unavailable."""
    if not math.isfinite(value):
        return missing_score
    return bounded01(value)


# ============================================================
# Coherence-envelope sector scores
# ============================================================

def star_score(row: pd.Series, widths: EnvelopeWidths) -> float:
    mstar = safe_float(row.get("stellar_mass_solar"))
    teff = safe_float(row.get("stellar_teff_k"))
    age = safe_float(row.get("stellar_age_gyr"))

    m_score = neutral_if_missing(
        gaussian_delta(mstar - EARTH["stellar_mass_solar"], widths.stellar_mass_width)
        if math.isfinite(mstar) else np.nan
    )
    t_score = neutral_if_missing(
        gaussian_delta(teff - EARTH["stellar_teff_k"], widths.stellar_teff_width)
        if math.isfinite(teff) else np.nan
    )
    age_score = neutral_if_missing(
        gaussian_delta(age - EARTH["stellar_age_gyr"], widths.stellar_age_width)
        if math.isfinite(age) else np.nan
    )

    # Penalize very low-mass flare-prone M-dwarf environments in this toy layer.
    # This is a heuristic, not a hard astrophysical verdict.
    dwarf_penalty = 1.0
    if math.isfinite(mstar) and mstar < 0.55:
        dwarf_penalty = 0.70
    if math.isfinite(teff) and teff < 3900:
        dwarf_penalty *= 0.85

    return bounded01(weighted_mean([
        (m_score, 0.45),
        (t_score, 0.40),
        (age_score, 0.15),
    ]) * dwarf_penalty)


def orbit_score(row: pd.Series, widths: EnvelopeWidths) -> float:
    sma = infer_semi_major_axis(row)
    ecc = safe_float(row.get("eccentricity"))
    ins = infer_insolation(row)

    sma_score = neutral_if_missing(
        gaussian_log_ratio(sma, EARTH["semi_major_axis_au"], widths.sma_log_width)
        if math.isfinite(sma) else np.nan
    )
    ecc_score = neutral_if_missing(
        gaussian_delta(ecc - EARTH["eccentricity"], widths.eccentricity_width)
        if math.isfinite(ecc) else np.nan
    )
    flux_score = neutral_if_missing(
        gaussian_log_ratio(ins, EARTH["insolation_earth"], widths.insolation_log_width)
        if math.isfinite(ins) else np.nan
    )

    return weighted_mean([
        (sma_score, 0.20),
        (ecc_score, 0.30),
        (flux_score, 0.50),
    ])


def planet_score(row: pd.Series, widths: EnvelopeWidths) -> float:
    mass = safe_float(row.get("planet_mass_earth"))
    radius = safe_float(row.get("planet_radius_earth"))
    density = infer_density_earth_units(row)
    gravity = infer_surface_gravity_earth_units(row)

    mass_score = neutral_if_missing(
        gaussian_log_ratio(mass, EARTH["planet_mass_earth"], widths.mass_log_width)
        if math.isfinite(mass) else np.nan
    )
    radius_score = neutral_if_missing(
        gaussian_log_ratio(radius, EARTH["planet_radius_earth"], widths.radius_log_width)
        if math.isfinite(radius) else np.nan
    )
    density_score = neutral_if_missing(
        gaussian_log_ratio(density, EARTH["density_earth"], widths.density_log_width)
        if math.isfinite(density) else np.nan
    )
    gravity_score = neutral_if_missing(
        gaussian_log_ratio(gravity, EARTH["surface_gravity_earth"], widths.gravity_log_width)
        if math.isfinite(gravity) else np.nan
    )

    # Soft rockiness gate. Very large radii probably mean volatile-rich/sub-Neptune.
    rockiness_gate = 1.0
    if math.isfinite(radius) and radius > 1.8:
        rockiness_gate *= math.exp(-0.5 * ((radius - 1.8) / 0.45) ** 2)
    if math.isfinite(mass) and mass > 10.0:
        rockiness_gate *= math.exp(-0.5 * ((math.log(mass / 10.0)) / 0.6) ** 2)

    return bounded01(weighted_mean([
        (mass_score, 0.30),
        (radius_score, 0.30),
        (density_score, 0.25),
        (gravity_score, 0.15),
    ]) * rockiness_gate)


def thermal_score(row: pd.Series, widths: EnvelopeWidths) -> float:
    ins = infer_insolation(row)
    teq = infer_equilibrium_temp(row)

    ins_score = neutral_if_missing(
        gaussian_log_ratio(ins, EARTH["insolation_earth"], widths.insolation_log_width)
        if math.isfinite(ins) else np.nan
    )
    teq_score = neutral_if_missing(
        gaussian_delta(teq - EARTH["equilibrium_temp_k"], widths.eq_temp_width)
        if math.isfinite(teq) else np.nan
    )

    return weighted_mean([
        (ins_score, 0.65),
        (teq_score, 0.35),
    ])


def atmosphere_score(row: pd.Series) -> float:
    """Atmosphere retention/probability score.

    If a user/catalog provides atmosphere_score, use it directly.
    Otherwise infer weakly from gravity, radius, and insolation.
    """
    provided = safe_float(row.get("atmosphere_score"))
    if math.isfinite(provided):
        return bounded01(provided)

    g = infer_surface_gravity_earth_units(row)
    radius = safe_float(row.get("planet_radius_earth"))
    ins = infer_insolation(row)

    g_score = neutral_if_missing(gaussian_log_ratio(g, 1.0, 0.65) if math.isfinite(g) else np.nan)
    r_score = neutral_if_missing(gaussian_log_ratio(radius, 1.0, 0.55) if math.isfinite(radius) else np.nan)
    ins_score = neutral_if_missing(gaussian_log_ratio(ins, 1.0, 0.75) if math.isfinite(ins) else np.nan)

    return weighted_mean([
        (g_score, 0.45),
        (r_score, 0.25),
        (ins_score, 0.30),
    ])


def magnetosphere_score(row: pd.Series, widths: EnvelopeWidths) -> float:
    """Magnetosphere/dynamo proxy.

    If provided, use magnetosphere_score directly.
    Otherwise infer from mass/radius/density/rotation.
    """
    provided = safe_float(row.get("magnetosphere_score"))
    if math.isfinite(provided):
        return bounded01(provided)

    density = infer_density_earth_units(row)
    mass = safe_float(row.get("planet_mass_earth"))
    rotation = safe_float(row.get("rotation_period_days"))

    density_score = neutral_if_missing(
        gaussian_log_ratio(density, 1.0, widths.density_log_width)
        if math.isfinite(density) else np.nan
    )
    mass_score = neutral_if_missing(
        gaussian_log_ratio(mass, 1.0, widths.mass_log_width)
        if math.isfinite(mass) else np.nan
    )
    rotation_score = neutral_if_missing(
        gaussian_log_ratio(rotation, 1.0, widths.rotation_log_width)
        if math.isfinite(rotation) else np.nan,
        missing_score=0.50,
    )

    return weighted_mean([
        (density_score, 0.35),
        (mass_score, 0.30),
        (rotation_score, 0.35),
    ])


def tidal_score(row: pd.Series, widths: EnvelopeWidths) -> float:
    """Tidal/rotational stability score.

    A large stabilizing moon scores high if known. Most exoplanets have unknown
    moon status, so missing moon data is treated as neutral.
    """
    moon = safe_float(row.get("has_large_moon"))
    rotation = safe_float(row.get("rotation_period_days"))
    ecc = safe_float(row.get("eccentricity"))

    moon_score = 0.55 if not math.isfinite(moon) else bounded01(moon)
    rot_score = neutral_if_missing(
        gaussian_log_ratio(rotation, 1.0, widths.rotation_log_width)
        if math.isfinite(rotation) else np.nan,
        missing_score=0.50,
    )
    ecc_score = neutral_if_missing(
        gaussian_delta(ecc - EARTH["eccentricity"], widths.eccentricity_width)
        if math.isfinite(ecc) else np.nan
    )

    return weighted_mean([
        (moon_score, 0.45),
        (rot_score, 0.20),
        (ecc_score, 0.35),
    ])


def sigma_coherence_score(row: pd.Series, widths: EnvelopeWidths) -> Tuple[float, float]:
    """Speculative Flux Cosmology local sigma score.

    This is a dimensionless proxy for how close the local environment is to
    Earth's coherence/actualization envelope.

    We define sigma_proxy as a combination of:
      - orbital flux closeness
      - stellar stability closeness
      - planetary gravity/density closeness
      - eccentricity/orderliness
      - magnetic/atmospheric retention proxies

    Earth should evaluate close to sigma_proxy = 1.
    """
    ins = infer_insolation(row)
    mstar = safe_float(row.get("stellar_mass_solar"))
    density = infer_density_earth_units(row)
    gravity = infer_surface_gravity_earth_units(row)
    ecc = safe_float(row.get("eccentricity"))

    flux_term = neutral_if_missing(
        gaussian_log_ratio(ins, 1.0, widths.insolation_log_width)
        if math.isfinite(ins) else np.nan
    )
    star_term = neutral_if_missing(
        gaussian_delta(mstar - 1.0, widths.stellar_mass_width)
        if math.isfinite(mstar) else np.nan
    )
    density_term = neutral_if_missing(
        gaussian_log_ratio(density, 1.0, widths.density_log_width)
        if math.isfinite(density) else np.nan
    )
    gravity_term = neutral_if_missing(
        gaussian_log_ratio(gravity, 1.0, widths.gravity_log_width)
        if math.isfinite(gravity) else np.nan
    )
    order_term = neutral_if_missing(
        gaussian_delta(ecc - EARTH["eccentricity"], widths.eccentricity_width)
        if math.isfinite(ecc) else np.nan
    )
    atm_term = atmosphere_score(row)
    mag_term = magnetosphere_score(row, widths)

    sigma_proxy = weighted_mean([
        (flux_term, 0.24),
        (star_term, 0.16),
        (density_term, 0.14),
        (gravity_term, 0.12),
        (order_term, 0.12),
        (atm_term, 0.10),
        (mag_term, 0.12),
    ])

    sigma_score = gaussian_delta(sigma_proxy - 1.0, widths.sigma_width)
    return bounded01(sigma_score), bounded01(sigma_proxy)


# ============================================================
# Main scoring routine
# ============================================================

def score_catalog(df: pd.DataFrame, weights: ScoreWeights, widths: EnvelopeWidths) -> pd.DataFrame:
    weights = weights.normalized()
    rows: List[Dict[str, float]] = []

    for _, row in df.iterrows():
        s_star = star_score(row, widths)
        s_orbit = orbit_score(row, widths)
        s_planet = planet_score(row, widths)
        s_thermal = thermal_score(row, widths)
        s_atm = atmosphere_score(row)
        s_mag = magnetosphere_score(row, widths)
        s_tidal = tidal_score(row, widths)
        s_sigma, sigma_proxy = sigma_coherence_score(row, widths)

        final = weighted_mean([
            (s_star, weights.star),
            (s_orbit, weights.orbit),
            (s_planet, weights.planet),
            (s_thermal, weights.thermal),
            (s_atm, weights.atmosphere),
            (s_mag, weights.magnetosphere),
            (s_tidal, weights.tidal),
            (s_sigma, weights.sigma),
        ])

        out = dict(row)
        out.update({
            "semi_major_axis_au_inferred": infer_semi_major_axis(row),
            "insolation_earth_inferred": infer_insolation(row),
            "density_earth_inferred": infer_density_earth_units(row),
            "surface_gravity_earth_inferred": infer_surface_gravity_earth_units(row),
            "equilibrium_temp_k_inferred": infer_equilibrium_temp(row),
            "score_star": s_star,
            "score_orbit": s_orbit,
            "score_planet": s_planet,
            "score_thermal": s_thermal,
            "score_atmosphere": s_atm,
            "score_magnetosphere": s_mag,
            "score_tidal": s_tidal,
            "sigma_proxy": sigma_proxy,
            "score_sigma": s_sigma,
            "ECEM_score": final,
        })
        rows.append(out)

    result = pd.DataFrame(rows)
    result = result.sort_values("ECEM_score", ascending=False).reset_index(drop=True)
    result.insert(0, "rank", np.arange(1, len(result) + 1))
    return result


# ============================================================
# Demo catalog
# ============================================================

def demo_catalog() -> pd.DataFrame:
    """Small demonstration catalog with well-known approximate examples.

    The values are intentionally approximate placeholders. Replace with a
    current catalog export for research use.
    """
    return pd.DataFrame([
        {
            "name": "Earth",
            "planet_mass_earth": 1.0,
            "planet_radius_earth": 1.0,
            "orbital_period_days": 365.25,
            "semi_major_axis_au": 1.0,
            "eccentricity": 0.0167,
            "stellar_mass_solar": 1.0,
            "stellar_radius_solar": 1.0,
            "stellar_teff_k": 5772,
            "stellar_luminosity_solar": 1.0,
            "stellar_age_gyr": 4.57,
            "equilibrium_temp_k": 255,
            "insolation_earth": 1.0,
            "rotation_period_days": 1.0,
            "has_large_moon": 1.0,
            "atmosphere_score": 1.0,
            "magnetosphere_score": 1.0,
        },
        {
            "name": "Kepler-452 b",
            "planet_mass_earth": 5.0,
            "planet_radius_earth": 1.63,
            "orbital_period_days": 384.8,
            "semi_major_axis_au": 1.05,
            "eccentricity": np.nan,
            "stellar_mass_solar": 1.04,
            "stellar_radius_solar": 1.11,
            "stellar_teff_k": 5757,
            "stellar_luminosity_solar": 1.20,
            "stellar_age_gyr": 6.0,
            "insolation_earth": 1.1,
        },
        {
            "name": "Kepler-186 f",
            "planet_mass_earth": 1.4,
            "planet_radius_earth": 1.17,
            "orbital_period_days": 129.9,
            "semi_major_axis_au": 0.43,
            "eccentricity": np.nan,
            "stellar_mass_solar": 0.54,
            "stellar_radius_solar": 0.52,
            "stellar_teff_k": 3755,
            "stellar_luminosity_solar": 0.041,
            "insolation_earth": 0.32,
        },
        {
            "name": "TRAPPIST-1 e",
            "planet_mass_earth": 0.69,
            "planet_radius_earth": 0.92,
            "orbital_period_days": 6.10,
            "semi_major_axis_au": 0.029,
            "eccentricity": 0.005,
            "stellar_mass_solar": 0.089,
            "stellar_radius_solar": 0.12,
            "stellar_teff_k": 2566,
            "stellar_luminosity_solar": 0.00055,
            "insolation_earth": 0.66,
            "rotation_period_days": 6.10,
        },
        {
            "name": "Proxima Centauri b",
            "planet_mass_earth": 1.27,
            "planet_radius_earth": np.nan,
            "orbital_period_days": 11.19,
            "semi_major_axis_au": 0.0485,
            "eccentricity": 0.02,
            "stellar_mass_solar": 0.122,
            "stellar_radius_solar": 0.154,
            "stellar_teff_k": 3042,
            "stellar_luminosity_solar": 0.00155,
            "insolation_earth": 0.65,
            "rotation_period_days": 11.19,
        },
        {
            "name": "TOI-700 d",
            "planet_mass_earth": 1.7,
            "planet_radius_earth": 1.14,
            "orbital_period_days": 37.4,
            "semi_major_axis_au": 0.163,
            "eccentricity": 0.03,
            "stellar_mass_solar": 0.416,
            "stellar_radius_solar": 0.42,
            "stellar_teff_k": 3480,
            "stellar_luminosity_solar": 0.023,
            "insolation_earth": 0.86,
            "rotation_period_days": 37.4,
        },
        {
            "name": "Mars",
            "planet_mass_earth": 0.107,
            "planet_radius_earth": 0.532,
            "orbital_period_days": 687,
            "semi_major_axis_au": 1.524,
            "eccentricity": 0.0934,
            "stellar_mass_solar": 1.0,
            "stellar_radius_solar": 1.0,
            "stellar_teff_k": 5772,
            "stellar_luminosity_solar": 1.0,
            "stellar_age_gyr": 4.57,
            "insolation_earth": 0.43,
            "rotation_period_days": 1.03,
            "has_large_moon": 0.0,
            "atmosphere_score": 0.25,
            "magnetosphere_score": 0.05,
        },
        {
            "name": "Venus",
            "planet_mass_earth": 0.815,
            "planet_radius_earth": 0.949,
            "orbital_period_days": 224.7,
            "semi_major_axis_au": 0.723,
            "eccentricity": 0.0068,
            "stellar_mass_solar": 1.0,
            "stellar_radius_solar": 1.0,
            "stellar_teff_k": 5772,
            "stellar_luminosity_solar": 1.0,
            "stellar_age_gyr": 4.57,
            "insolation_earth": 1.91,
            "rotation_period_days": 243.0,
            "has_large_moon": 0.0,
            "atmosphere_score": 0.45,
            "magnetosphere_score": 0.05,
        },
        {
            "name": "Sub-Neptune control",
            "planet_mass_earth": 8.0,
            "planet_radius_earth": 2.5,
            "orbital_period_days": 80,
            "semi_major_axis_au": 0.38,
            "eccentricity": 0.05,
            "stellar_mass_solar": 0.8,
            "stellar_radius_solar": 0.75,
            "stellar_teff_k": 5000,
            "stellar_luminosity_solar": 0.35,
            "insolation_earth": 2.4,
        },
    ])


# ============================================================
# Plotting
# ============================================================

def plot_scores(ranked: pd.DataFrame, outdir: str, top: int = 20) -> None:
    os.makedirs(outdir, exist_ok=True)
    top_df = ranked.head(top).iloc[::-1]

    plt.figure(figsize=(10, max(5, 0.38 * len(top_df))))
    plt.barh(top_df["name"].astype(str), top_df["ECEM_score"])
    plt.xlabel("Earth Coherence Envelope Match score")
    plt.ylabel("Planet")
    plt.title("Top Earth coherence-envelope matches")
    plt.xlim(0, 1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "earth_coherence_scores.png"), dpi=180)
    plt.close()

    plt.figure(figsize=(8, 6))
    x = ranked["planet_radius_earth"].astype(float)
    y = ranked["planet_mass_earth"].astype(float)
    c = ranked["ECEM_score"].astype(float)
    plt.scatter(x, y, s=80 + 220 * c, c=c)
    for _, r in ranked.head(top).iterrows():
        rx = safe_float(r.get("planet_radius_earth"))
        ry = safe_float(r.get("planet_mass_earth"))
        if math.isfinite(rx) and math.isfinite(ry):
            plt.text(rx, ry, str(r["name"]), fontsize=8)
    plt.scatter([1.0], [1.0], marker="*", s=220, label="Earth reference")
    plt.xlabel("Radius [Earth radii]")
    plt.ylabel("Mass [Earth masses]")
    plt.title("Mass-radius coherence space")
    plt.colorbar(label="ECEM score")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "earth_coherence_mass_radius.png"), dpi=180)
    plt.close()

    plt.figure(figsize=(8, 6))
    x = ranked["semi_major_axis_au_inferred"].astype(float)
    y = ranked["insolation_earth_inferred"].astype(float)
    c = ranked["ECEM_score"].astype(float)
    plt.scatter(x, y, s=80 + 220 * c, c=c)
    for _, r in ranked.head(top).iterrows():
        rx = safe_float(r.get("semi_major_axis_au_inferred"))
        ry = safe_float(r.get("insolation_earth_inferred"))
        if math.isfinite(rx) and math.isfinite(ry):
            plt.text(rx, ry, str(r["name"]), fontsize=8)
    plt.scatter([1.0], [1.0], marker="*", s=220, label="Earth reference")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Semi-major axis [AU]")
    plt.ylabel("Insolation [Earth = 1]")
    plt.title("Flux-orbit coherence space")
    plt.colorbar(label="ECEM score")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "earth_coherence_flux_orbit.png"), dpi=180)
    plt.close()


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank exoplanets by Earth Coherence Envelope Match score."
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Optional input CSV. If omitted, uses a built-in demo catalog.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=".",
        help="Output directory for CSV and plots.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of top-ranked planets to print/plot.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Disable plot generation.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.input:
        df = pd.read_csv(args.input)
    else:
        df = demo_catalog()

    if "name" not in df.columns:
        df["name"] = [f"planet_{i}" for i in range(len(df))]

    weights = ScoreWeights()
    widths = EnvelopeWidths()
    ranked = score_catalog(df, weights, widths)

    out_csv = os.path.join(args.outdir, "earth_coherence_ranked.csv")
    ranked.to_csv(out_csv, index=False)

    columns = [
        "rank",
        "name",
        "ECEM_score",
        "score_star",
        "score_orbit",
        "score_planet",
        "score_thermal",
        "score_atmosphere",
        "score_magnetosphere",
        "score_tidal",
        "score_sigma",
        "sigma_proxy",
        "insolation_earth_inferred",
        "density_earth_inferred",
        "surface_gravity_earth_inferred",
    ]
    visible_columns = [c for c in columns if c in ranked.columns]

    print("\nEarth Coherence Envelope Match ranking")
    print("=" * 48)
    print(ranked[visible_columns].head(args.top).to_string(index=False))
    print(f"\nSaved ranked catalog: {out_csv}")

    if not args.no_plots:
        plot_scores(ranked, args.outdir, top=args.top)
        print("Saved plots:")
        print(f"  {os.path.join(args.outdir, 'earth_coherence_scores.png')}")
        print(f"  {os.path.join(args.outdir, 'earth_coherence_mass_radius.png')}")
        print(f"  {os.path.join(args.outdir, 'earth_coherence_flux_orbit.png')}")


if __name__ == "__main__":
    main()
