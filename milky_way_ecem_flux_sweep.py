#!/usr/bin/env python3
"""
milky_way_ecem_flux_sweep.py

Flux Cosmology / ECEM extension for Milky Way stellar-exoplanet systems.

Purpose
-------
This script extends the Earth Coherence Envelope Match (ECEM) idea from a
static Earth-envelope similarity score into a flux-dynamics sweep:

    exoplanet catalog -> ECEM_base -> flux tension/history adjustment -> ECEM_flux

The goal is not to prove habitability. The goal is to rank systems by how
closely their star--planet--orbit--thermal--atmosphere--magnetosphere--tidal
boundary conditions resemble Earth's coherence envelope under a local
flux-tension interpretation.

Core outputs
------------
- ECEM_base: Earth-envelope similarity from ordinary system parameters.
- flux_tension_proxy: local stress term from insolation, stellar activity,
  tidal forcing, and age/history.
- flux_history_score: whether the system has likely occupied a stable,
  Earth-like coherence branch long enough.
- surface_state_score: rejects Venus-like failed branches.
- ECEM_flux: base ECEM multiplied by the flux/history adjustment.
- classification: broad interpretive class for repo triage.

Inputs
------
Use the built-in approximate Milky Way candidate catalog:

    python milky_way_ecem_flux_sweep.py

Or pass your own CSV:

    python milky_way_ecem_flux_sweep.py --input my_exoplanets.csv --top 50

Recognized direct columns:

    name
    star_class
    star_mass_solar
    star_radius_solar
    star_teff_k
    star_luminosity_solar
    star_age_gyr
    stellar_activity_score
    planet_mass_earth
    planet_radius_earth
    semi_major_axis_au
    orbital_period_days
    eccentricity
    insolation_earth
    equilibrium_temp_k
    rotation_period_days
    atmosphere_score
    magnetosphere_score
    has_large_moon
    surface_state_score

Common NASA Exoplanet Archive aliases are also auto-mapped, e.g. pl_name,
pl_bmasse, pl_rade, pl_orbper, pl_orbsmax, pl_orbeccen, st_mass, st_rad,
st_teff, st_lum. For NASA `st_lum`, the value is treated as log10(L/Lsun)
when mapped.

Outputs
-------
Default output names:

    milky_way_ecem_flux_ranked.csv
    milky_way_ecem_flux_scores.png
    milky_way_ecem_tidal_vs_flux.png
    milky_way_ecem_starclass.png
    milky_way_ecem_mass_radius.png
    milky_way_ecem_flux_orbit.png

Notes
-----
The built-in catalog values are approximate placeholders for model testing.
For research use, feed a current NASA Exoplanet Archive / Planetary Systems
Composite Parameters CSV export and inspect missing fields.
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, asdict
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Earth reference envelope
# ============================================================

EARTH = {
    "planet_mass_earth": 1.0,
    "planet_radius_earth": 1.0,
    "density_earth": 1.0,
    "surface_gravity_earth": 1.0,
    "semi_major_axis_au": 1.0,
    "orbital_period_days": 365.25,
    "eccentricity": 0.0167,
    "star_mass_solar": 1.0,
    "star_radius_solar": 1.0,
    "star_teff_k": 5772.0,
    "star_luminosity_solar": 1.0,
    "star_age_gyr": 4.57,
    "insolation_earth": 1.0,
    "equilibrium_temp_k": 255.0,
    "rotation_period_days": 1.0,
    "atmosphere_score": 1.0,
    "magnetosphere_score": 1.0,
    "surface_state_score": 1.0,
}

STAR_CLASS_COHERENCE = {
    "O": 0.10,
    "B": 0.20,
    "A": 0.45,
    "F": 0.65,
    "G": 1.00,
    "K": 1.05,
    "M": 0.70,
    "UNKNOWN": 0.75,
}


# ============================================================
# Configuration
# ============================================================

@dataclass
class BaseWeights:
    star: float = 0.15
    orbit: float = 0.18
    planet: float = 0.20
    thermal: float = 0.13
    atmosphere: float = 0.10
    magnetosphere: float = 0.08
    tidal: float = 0.06
    surface: float = 0.06
    sigma: float = 0.04

    def normalized(self) -> "BaseWeights":
        total = sum(asdict(self).values())
        vals = {k: v / total for k, v in asdict(self).items()}
        return BaseWeights(**vals)


@dataclass
class FluxWeights:
    stellar_class: float = 0.16
    stellar_activity: float = 0.15
    flux_tension: float = 0.22
    tidal_damage: float = 0.17
    surface_state: float = 0.18
    history: float = 0.12

    def normalized(self) -> "FluxWeights":
        total = sum(asdict(self).values())
        vals = {k: v / total for k, v in asdict(self).items()}
        return FluxWeights(**vals)


@dataclass
class Widths:
    mass_log_width: float = 0.55
    radius_log_width: float = 0.35
    density_log_width: float = 0.45
    gravity_log_width: float = 0.42
    insolation_log_width: float = 0.45
    semi_major_axis_log_width: float = 0.70
    eccentricity_width: float = 0.18
    star_mass_width: float = 0.35
    star_teff_width: float = 900.0
    star_age_width: float = 4.0
    eq_temp_width: float = 45.0
    rotation_log_width: float = 1.20
    sigma_width: float = 0.35
    flux_tension_log_width: float = 0.90
    tidal_log_width: float = 2.40


# ============================================================
# Utilities
# ============================================================

def safe_float(value, default: float = np.nan) -> float:
    try:
        if value is None:
            return default
        if isinstance(value, str) and value.strip() == "":
            return default
        x = float(value)
        return x if math.isfinite(x) else default
    except Exception:
        return default


def bounded01(x: float) -> float:
    if not math.isfinite(x):
        return 0.0
    return max(0.0, min(1.0, x))


def gaussian_delta(delta: float, width: float, missing: float = 0.55) -> float:
    if not math.isfinite(delta) or width <= 0:
        return missing
    return bounded01(math.exp(-0.5 * (delta / width) ** 2))


def gaussian_log_ratio(value: float, reference: float, width: float, missing: float = 0.55) -> float:
    value = safe_float(value)
    reference = safe_float(reference)
    if value <= 0 or reference <= 0 or width <= 0:
        return missing
    return bounded01(math.exp(-0.5 * (math.log(value / reference) / width) ** 2))


def weighted_mean(scores: Iterable[Tuple[float, float]], missing: float = 0.55) -> float:
    num = 0.0
    den = 0.0
    for score, weight in scores:
        if math.isfinite(score) and weight > 0:
            num += bounded01(score) * weight
            den += weight
    return bounded01(num / den) if den > 0 else missing


def first_valid(row: pd.Series, keys: List[str], default=np.nan):
    for key in keys:
        if key in row.index:
            val = row.get(key)
            if isinstance(val, str):
                if val.strip():
                    return val
            elif pd.notna(val):
                return val
    return default


# ============================================================
# Input column normalization
# ============================================================

def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Map common archive columns onto the script's internal schema."""
    out = df.copy()

    mappings = {
        "name": ["name", "pl_name", "planet_name", "Planet", "planet"],
        "planet_mass_earth": ["planet_mass_earth", "pl_bmasse", "pl_masse", "mass_earth", "M_earth"],
        "planet_radius_earth": ["planet_radius_earth", "pl_rade", "radius_earth", "R_earth"],
        "orbital_period_days": ["orbital_period_days", "pl_orbper", "period_days", "P_days"],
        "semi_major_axis_au": ["semi_major_axis_au", "pl_orbsmax", "a_au", "sma_au"],
        "eccentricity": ["eccentricity", "pl_orbeccen", "ecc"],
        "star_mass_solar": ["star_mass_solar", "stellar_mass_solar", "st_mass", "Mstar"],
        "star_radius_solar": ["star_radius_solar", "stellar_radius_solar", "st_rad", "Rstar"],
        "star_teff_k": ["star_teff_k", "stellar_teff_k", "st_teff", "Teff"],
        "star_age_gyr": ["star_age_gyr", "stellar_age_gyr", "st_age", "age_gyr"],
        "insolation_earth": ["insolation_earth", "pl_insol", "insolation", "S_earth"],
        "equilibrium_temp_k": ["equilibrium_temp_k", "pl_eqt", "teq_k", "T_eq"],
        "rotation_period_days": ["rotation_period_days", "rot_period_days", "P_rot"],
        "atmosphere_score": ["atmosphere_score", "atm_score"],
        "magnetosphere_score": ["magnetosphere_score", "mag_score", "dynamo_score"],
        "surface_state_score": ["surface_state_score", "surface_score"],
        "stellar_activity_score": ["stellar_activity_score", "activity_score"],
        "has_large_moon": ["has_large_moon", "moon_score", "large_moon"],
        "star_class": ["star_class", "spectral_type", "st_spectype"],
    }

    for dest, aliases in mappings.items():
        if dest not in out.columns:
            for alias in aliases:
                if alias in out.columns:
                    out[dest] = out[alias]
                    break

    # NASA Exoplanet Archive st_lum is usually log10(L/Lsun).
    if "star_luminosity_solar" not in out.columns:
        if "stellar_luminosity_solar" in out.columns:
            out["star_luminosity_solar"] = out["stellar_luminosity_solar"]
        elif "st_lum" in out.columns:
            vals = pd.to_numeric(out["st_lum"], errors="coerce")
            out["star_luminosity_solar"] = np.where(np.isfinite(vals), 10.0 ** vals, np.nan)
        elif "luminosity_solar" in out.columns:
            out["star_luminosity_solar"] = out["luminosity_solar"]

    if "name" not in out.columns:
        out["name"] = [f"planet_{i}" for i in range(len(out))]

    return out


# ============================================================
# Inference helpers
# ============================================================

def infer_star_class(row: pd.Series) -> str:
    raw = first_valid(row, ["star_class"], default="")
    if isinstance(raw, str) and raw.strip():
        s = raw.strip().upper()
        for c in ["O", "B", "A", "F", "G", "K", "M"]:
            if s.startswith(c):
                return c

    teff = safe_float(row.get("star_teff_k"))
    if not math.isfinite(teff):
        return "UNKNOWN"
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


def infer_luminosity(row: pd.Series) -> float:
    lum = safe_float(row.get("star_luminosity_solar"))
    if lum > 0:
        return lum
    r = safe_float(row.get("star_radius_solar"))
    t = safe_float(row.get("star_teff_k"))
    if r > 0 and t > 0:
        return r ** 2 * (t / EARTH["star_teff_k"]) ** 4
    return np.nan


def infer_sma(row: pd.Series) -> float:
    sma = safe_float(row.get("semi_major_axis_au"))
    if sma > 0:
        return sma
    period = safe_float(row.get("orbital_period_days"))
    mstar = safe_float(row.get("star_mass_solar"))
    if period > 0 and mstar > 0:
        return (mstar * (period / 365.25) ** 2) ** (1.0 / 3.0)
    return np.nan


def infer_insolation(row: pd.Series) -> float:
    ins = safe_float(row.get("insolation_earth"))
    if ins > 0:
        return ins
    lum = infer_luminosity(row)
    sma = infer_sma(row)
    if lum > 0 and sma > 0:
        return lum / (sma ** 2)
    return np.nan


def infer_eq_temp(row: pd.Series) -> float:
    teq = safe_float(row.get("equilibrium_temp_k"))
    if teq > 0:
        return teq
    ins = infer_insolation(row)
    if ins > 0:
        return EARTH["equilibrium_temp_k"] * ins ** 0.25
    return np.nan


def infer_density(row: pd.Series) -> float:
    m = safe_float(row.get("planet_mass_earth"))
    r = safe_float(row.get("planet_radius_earth"))
    if m > 0 and r > 0:
        return m / (r ** 3)
    return np.nan


def infer_gravity(row: pd.Series) -> float:
    m = safe_float(row.get("planet_mass_earth"))
    r = safe_float(row.get("planet_radius_earth"))
    if m > 0 and r > 0:
        return m / (r ** 2)
    return np.nan


def infer_stellar_activity_raw(row: pd.Series) -> float:
    """Raw activity/stress proxy, Earth near 1. Higher = harsher."""
    provided = safe_float(row.get("stellar_activity_raw"))
    if provided > 0:
        return provided

    cls = infer_star_class(row)
    mstar = safe_float(row.get("star_mass_solar"))
    age = safe_float(row.get("star_age_gyr"))

    base = {
        "F": 1.35,
        "G": 1.00,
        "K": 0.85,
        "M": 1.45,
        "A": 2.30,
        "B": 4.00,
        "O": 6.00,
        "UNKNOWN": 1.20,
    }.get(cls, 1.20)

    if cls == "M" and math.isfinite(mstar):
        if mstar < 0.20:
            base *= 1.35
        elif mstar > 0.45:
            base *= 0.90

    # Very young stars are usually more active. Old stars are calmer but can
    # represent longer flux-history exhaustion depending on luminosity evolution.
    if math.isfinite(age):
        if age < 1.0:
            base *= 1.45
        elif age > 8.0:
            base *= 1.08

    return max(base, 0.05)


def stellar_activity_score(row: pd.Series) -> float:
    provided = safe_float(row.get("stellar_activity_score"))
    if math.isfinite(provided):
        return bounded01(provided)
    raw = infer_stellar_activity_raw(row)
    return gaussian_log_ratio(raw, 1.0, 0.70, missing=0.55)


def tidal_strength(row: pd.Series) -> float:
    """Earth-normalized tidal forcing proxy: Mstar/a^3."""
    mstar = safe_float(row.get("star_mass_solar"))
    sma = infer_sma(row)
    if mstar > 0 and sma > 0:
        return mstar / (sma ** 3)
    return np.nan


# ============================================================
# Base ECEM sectors
# ============================================================

def star_score(row: pd.Series, widths: Widths) -> float:
    cls = infer_star_class(row)
    mstar = safe_float(row.get("star_mass_solar"))
    teff = safe_float(row.get("star_teff_k"))
    age = safe_float(row.get("star_age_gyr"))
    class_mult = STAR_CLASS_COHERENCE.get(cls, STAR_CLASS_COHERENCE["UNKNOWN"])

    m_score = gaussian_delta(mstar - 1.0, widths.star_mass_width) if math.isfinite(mstar) else 0.55
    t_score = gaussian_delta(teff - EARTH["star_teff_k"], widths.star_teff_width) if math.isfinite(teff) else 0.55
    age_score = gaussian_delta(age - EARTH["star_age_gyr"], widths.star_age_width) if math.isfinite(age) else 0.55

    return bounded01(weighted_mean([(m_score, 0.42), (t_score, 0.40), (age_score, 0.18)]) * class_mult)


def orbit_score(row: pd.Series, widths: Widths) -> float:
    sma = infer_sma(row)
    ecc = safe_float(row.get("eccentricity"))
    ins = infer_insolation(row)
    sma_score = gaussian_log_ratio(sma, 1.0, widths.semi_major_axis_log_width) if math.isfinite(sma) else 0.55
    ecc_score = gaussian_delta(ecc - EARTH["eccentricity"], widths.eccentricity_width) if math.isfinite(ecc) else 0.55
    ins_score = gaussian_log_ratio(ins, 1.0, widths.insolation_log_width) if math.isfinite(ins) else 0.55
    return weighted_mean([(sma_score, 0.20), (ecc_score, 0.25), (ins_score, 0.55)])


def planet_score(row: pd.Series, widths: Widths) -> float:
    m = safe_float(row.get("planet_mass_earth"))
    r = safe_float(row.get("planet_radius_earth"))
    dens = infer_density(row)
    g = infer_gravity(row)

    m_score = gaussian_log_ratio(m, 1.0, widths.mass_log_width) if math.isfinite(m) else 0.55
    r_score = gaussian_log_ratio(r, 1.0, widths.radius_log_width) if math.isfinite(r) else 0.55
    d_score = gaussian_log_ratio(dens, 1.0, widths.density_log_width) if math.isfinite(dens) else 0.55
    g_score = gaussian_log_ratio(g, 1.0, widths.gravity_log_width) if math.isfinite(g) else 0.55

    rock_gate = 1.0
    if math.isfinite(r) and r > 1.75:
        rock_gate *= math.exp(-0.5 * ((r - 1.75) / 0.42) ** 2)
    if math.isfinite(m) and m > 10.0:
        rock_gate *= math.exp(-0.5 * (math.log(m / 10.0) / 0.55) ** 2)
    if math.isfinite(dens) and dens < 0.35:
        rock_gate *= 0.55

    return bounded01(weighted_mean([(m_score, 0.30), (r_score, 0.32), (d_score, 0.24), (g_score, 0.14)]) * rock_gate)


def thermal_score(row: pd.Series, widths: Widths) -> float:
    ins = infer_insolation(row)
    teq = infer_eq_temp(row)
    ins_score = gaussian_log_ratio(ins, 1.0, widths.insolation_log_width) if math.isfinite(ins) else 0.55
    teq_score = gaussian_delta(teq - 255.0, widths.eq_temp_width) if math.isfinite(teq) else 0.55
    return weighted_mean([(ins_score, 0.65), (teq_score, 0.35)])


def atmosphere_score(row: pd.Series, widths: Widths) -> float:
    provided = safe_float(row.get("atmosphere_score"))
    if math.isfinite(provided):
        return bounded01(provided)
    g = infer_gravity(row)
    r = safe_float(row.get("planet_radius_earth"))
    ins = infer_insolation(row)
    cls = infer_star_class(row)

    g_score = gaussian_log_ratio(g, 1.0, 0.65) if math.isfinite(g) else 0.55
    r_score = gaussian_log_ratio(r, 1.0, 0.55) if math.isfinite(r) else 0.55
    ins_score = gaussian_log_ratio(ins, 1.0, 0.80) if math.isfinite(ins) else 0.55
    flare_penalty = 0.78 if cls == "M" and stellar_activity_score(row) < 0.55 else 1.0

    return bounded01(weighted_mean([(g_score, 0.45), (r_score, 0.20), (ins_score, 0.35)]) * flare_penalty)


def magnetosphere_score(row: pd.Series, widths: Widths) -> float:
    provided = safe_float(row.get("magnetosphere_score"))
    if math.isfinite(provided):
        return bounded01(provided)
    dens = infer_density(row)
    m = safe_float(row.get("planet_mass_earth"))
    rot = safe_float(row.get("rotation_period_days"))
    cls = infer_star_class(row)

    d_score = gaussian_log_ratio(dens, 1.0, widths.density_log_width) if math.isfinite(dens) else 0.55
    m_score = gaussian_log_ratio(m, 1.0, widths.mass_log_width) if math.isfinite(m) else 0.55
    rot_score = gaussian_log_ratio(rot, 1.0, widths.rotation_log_width, missing=0.50) if math.isfinite(rot) else 0.50

    # Synchronous/slow rotation can reduce dynamo confidence in this proxy.
    lock_penalty = 1.0
    if math.isfinite(rot) and rot > 20.0:
        lock_penalty *= 0.82
    if cls == "M" and math.isfinite(rot) and rot > 5.0:
        lock_penalty *= 0.90

    return bounded01(weighted_mean([(d_score, 0.38), (m_score, 0.27), (rot_score, 0.35)]) * lock_penalty)


def tidal_score(row: pd.Series, widths: Widths) -> float:
    moon = safe_float(row.get("has_large_moon"))
    ecc = safe_float(row.get("eccentricity"))
    rot = safe_float(row.get("rotation_period_days"))
    tau = tidal_strength(row)

    moon_score = 0.55 if not math.isfinite(moon) else bounded01(moon)
    ecc_score = gaussian_delta(ecc - EARTH["eccentricity"], widths.eccentricity_width) if math.isfinite(ecc) else 0.55
    rot_score = gaussian_log_ratio(rot, 1.0, widths.rotation_log_width, missing=0.50) if math.isfinite(rot) else 0.50
    tau_score = gaussian_log_ratio(tau, 1.0, widths.tidal_log_width, missing=0.55) if math.isfinite(tau) else 0.55

    return weighted_mean([(moon_score, 0.25), (ecc_score, 0.25), (rot_score, 0.20), (tau_score, 0.30)])


def venus_branch_penalty(row: pd.Series) -> float:
    ins = infer_insolation(row)
    atm = atmosphere_score(row, Widths())
    mag = magnetosphere_score(row, Widths())
    rot = safe_float(row.get("rotation_period_days"))
    teq = infer_eq_temp(row)

    penalty = 1.0
    if math.isfinite(ins) and ins > 1.35:
        penalty *= math.exp(-2.5 * (ins - 1.35) ** 2)
    if math.isfinite(teq) and teq > 285:
        penalty *= math.exp(-0.5 * ((teq - 285.0) / 55.0) ** 2)
    if math.isfinite(ins) and ins > 1.10 and atm < 0.50:
        penalty *= 0.55
    if math.isfinite(ins) and ins > 1.10 and mag < 0.45:
        penalty *= 0.70
    if math.isfinite(rot) and rot > 100:
        penalty *= 0.75
    return bounded01(penalty)


def surface_state_score(row: pd.Series, widths: Widths) -> float:
    provided = safe_float(row.get("surface_state_score"))
    if math.isfinite(provided):
        return bounded01(provided)
    therm = thermal_score(row, widths)
    atm = atmosphere_score(row, widths)
    mag = magnetosphere_score(row, widths)
    penalty = venus_branch_penalty(row)
    return bounded01(weighted_mean([(therm, 0.45), (atm, 0.35), (mag, 0.20)]) * penalty)


def sigma_score(row: pd.Series, widths: Widths) -> Tuple[float, float]:
    flux = thermal_score(row, widths)
    star = star_score(row, widths)
    planet = planet_score(row, widths)
    atm = atmosphere_score(row, widths)
    mag = magnetosphere_score(row, widths)
    tidal = tidal_score(row, widths)
    surface = surface_state_score(row, widths)
    sigma_proxy = weighted_mean([
        (flux, 0.20), (star, 0.16), (planet, 0.17),
        (atm, 0.12), (mag, 0.12), (tidal, 0.10), (surface, 0.13),
    ])
    return gaussian_delta(sigma_proxy - 1.0, widths.sigma_width), sigma_proxy


def ecem_base_score(row: pd.Series, base_weights: BaseWeights, widths: Widths) -> Dict[str, float]:
    w = base_weights.normalized()
    s_star = star_score(row, widths)
    s_orbit = orbit_score(row, widths)
    s_planet = planet_score(row, widths)
    s_thermal = thermal_score(row, widths)
    s_atm = atmosphere_score(row, widths)
    s_mag = magnetosphere_score(row, widths)
    s_tidal = tidal_score(row, widths)
    s_surface = surface_state_score(row, widths)
    s_sigma, sigma_proxy = sigma_score(row, widths)

    base = weighted_mean([
        (s_star, w.star), (s_orbit, w.orbit), (s_planet, w.planet),
        (s_thermal, w.thermal), (s_atm, w.atmosphere), (s_mag, w.magnetosphere),
        (s_tidal, w.tidal), (s_surface, w.surface), (s_sigma, w.sigma),
    ])
    return {
        "score_star": s_star,
        "score_orbit": s_orbit,
        "score_planet": s_planet,
        "score_thermal": s_thermal,
        "score_atmosphere": s_atm,
        "score_magnetosphere": s_mag,
        "score_tidal": s_tidal,
        "score_surface_state": s_surface,
        "score_sigma": s_sigma,
        "sigma_proxy": sigma_proxy,
        "ECEM_base": base,
    }


# ============================================================
# Flux dynamics extension
# ============================================================

def flux_tension_proxy(row: pd.Series) -> float:
    """Local flux-tension proxy; Earth is near 1.

    Higher values mean harsher or more compressed coherence processing.
    This combines insolation, activity, weak tidal contribution, and age drift.
    """
    ins = infer_insolation(row)
    activity_raw = infer_stellar_activity_raw(row)
    tau = tidal_strength(row)
    age = safe_float(row.get("star_age_gyr"))

    ins_term = ins if ins > 0 else 1.0
    tau_term = tau ** 0.12 if tau > 0 else 1.0
    age_term = 1.0
    if math.isfinite(age):
        # Young stars: unsettled. Very old systems: mild exhaustion/evolution drift.
        if age < 1.0:
            age_term = 1.25
        elif age > 8.0:
            age_term = 1.12

    return max(ins_term * activity_raw * tau_term * age_term, 1e-6)


def flux_tension_score(row: pd.Series, widths: Widths) -> float:
    phi = flux_tension_proxy(row)
    return gaussian_log_ratio(phi, 1.0, widths.flux_tension_log_width, missing=0.55)


def tidal_damage_score(row: pd.Series, widths: Widths) -> float:
    """Tidal stress survival score.

    A compact M-dwarf planet can remain interesting, but it is not treated as
    Earth-equivalent. This score tracks envelope damage risk, not just orbital
    stability.
    """
    tau = tidal_strength(row)
    cls = infer_star_class(row)
    mag = magnetosphere_score(row, widths)
    atm = atmosphere_score(row, widths)

    if not math.isfinite(tau):
        base = 0.55
    else:
        # Broad tolerance allows compact settled branches but penalizes extremes.
        base = gaussian_log_ratio(tau, 1.0, widths.tidal_log_width, missing=0.55)
        if tau > 100:
            base *= 0.78
        if tau > 1000:
            base *= 0.70

    if cls == "M":
        base *= 0.90
    if mag < 0.45:
        base *= 0.88
    if atm < 0.45:
        base *= 0.88
    return bounded01(base)


def flux_history_score(row: pd.Series, widths: Widths) -> float:
    age = safe_float(row.get("star_age_gyr"))
    cls = infer_star_class(row)
    activity = stellar_activity_score(row)
    surface = surface_state_score(row, widths)
    star_mult = STAR_CLASS_COHERENCE.get(cls, 0.75)

    if not math.isfinite(age):
        age_score = 0.58
    elif age < 0.8:
        age_score = 0.35
    elif age < 1.5:
        age_score = 0.62
    elif age <= 7.5:
        age_score = 0.95
    elif age <= 10.0:
        age_score = 0.78
    else:
        age_score = 0.62

    # K dwarfs can remain high at older ages; hot stars decay quickly in viability.
    if cls == "K" and math.isfinite(age) and age > 6.0:
        age_score = max(age_score, 0.88)
    if cls in ["A", "B", "O"]:
        age_score *= 0.45
    if cls == "F" and math.isfinite(age) and age > 4.5:
        age_score *= 0.82

    return bounded01(weighted_mean([
        (age_score, 0.38),
        (activity, 0.24),
        (bounded01(star_mult), 0.18),
        (surface, 0.20),
    ]))


def flux_adjustment(row: pd.Series, flux_weights: FluxWeights, widths: Widths) -> Dict[str, float]:
    w = flux_weights.normalized()
    cls = infer_star_class(row)
    class_score = bounded01(STAR_CLASS_COHERENCE.get(cls, 0.75))
    activity = stellar_activity_score(row)
    tension = flux_tension_score(row, widths)
    tidal_damage = tidal_damage_score(row, widths)
    surface = surface_state_score(row, widths)
    history = flux_history_score(row, widths)

    adj = weighted_mean([
        (class_score, w.stellar_class),
        (activity, w.stellar_activity),
        (tension, w.flux_tension),
        (tidal_damage, w.tidal_damage),
        (surface, w.surface_state),
        (history, w.history),
    ])
    return {
        "star_class_inferred": cls,
        "star_class_coherence": class_score,
        "stellar_activity_raw": infer_stellar_activity_raw(row),
        "score_stellar_activity": activity,
        "tidal_strength_earth": tidal_strength(row),
        "flux_tension_proxy": flux_tension_proxy(row),
        "score_flux_tension": tension,
        "score_tidal_damage": tidal_damage,
        "score_flux_history": history,
        "flux_adjustment": adj,
    }


def classify(row_out: Dict[str, float]) -> str:
    base = row_out.get("ECEM_base", 0.0)
    flux = row_out.get("ECEM_flux", 0.0)
    r = row_out.get("planet_radius_earth", np.nan)
    ins = row_out.get("insolation_earth_inferred", np.nan)
    tau = row_out.get("tidal_strength_earth", np.nan)
    cls = row_out.get("star_class_inferred", "UNKNOWN")
    surface = row_out.get("score_surface_state", 0.0)

    if math.isfinite(r) and r > 1.8:
        return "sub-Neptune / volatile-rich control"
    if surface < 0.35 and math.isfinite(ins) and ins > 1.1:
        return "Venus-like failed coherence branch"
    if surface < 0.35 and math.isfinite(ins) and ins < 0.75:
        return "Mars/freeze-out weak-envelope branch"
    if cls == "M" and math.isfinite(tau) and tau > 100:
        return "compact M-dwarf high-tidal branch"
    if flux >= 0.72:
        return "high Earth-envelope flux match"
    if base >= 0.60 and flux >= 0.50:
        return "candidate coherence-envelope analog"
    if base >= 0.50:
        return "partial ECEM analog / needs data"
    return "low ECEM / control or weak analog"


# ============================================================
# Catalog scoring
# ============================================================

def score_catalog(df: pd.DataFrame, base_weights: BaseWeights, flux_weights: FluxWeights, widths: Widths) -> pd.DataFrame:
    df = normalize_columns(df)
    outputs: List[Dict[str, float]] = []

    for _, row in df.iterrows():
        base = ecem_base_score(row, base_weights, widths)
        flux = flux_adjustment(row, flux_weights, widths)

        out = dict(row)
        out.update({
            "semi_major_axis_au_inferred": infer_sma(row),
            "insolation_earth_inferred": infer_insolation(row),
            "equilibrium_temp_k_inferred": infer_eq_temp(row),
            "density_earth_inferred": infer_density(row),
            "surface_gravity_earth_inferred": infer_gravity(row),
        })
        out.update(base)
        out.update(flux)

        # The strict version requested: flux score is base ECEM multiplied by
        # a flux/history adjustment. This keeps Earth at 1 and demotes systems
        # that are superficially habitable but coherence-divergent.
        out["ECEM_flux"] = bounded01(out["ECEM_base"] * out["flux_adjustment"])
        out["classification"] = classify(out)
        outputs.append(out)

    ranked = pd.DataFrame(outputs).sort_values("ECEM_flux", ascending=False).reset_index(drop=True)
    ranked.insert(0, "rank_flux", np.arange(1, len(ranked) + 1))
    ranked.insert(1, "rank_base", ranked["ECEM_base"].rank(ascending=False, method="min").astype(int))
    return ranked


# ============================================================
# Built-in approximate Milky Way test catalog
# ============================================================

def demo_catalog() -> pd.DataFrame:
    """Approximate candidate/control catalog for model testing.

    Values are deliberately approximate and incomplete where appropriate.
    Replace with a current archive export for serious scoring.
    """
    return pd.DataFrame([
        # Solar-system controls
        dict(name="Earth", star_class="G", star_mass_solar=1.0, star_radius_solar=1.0, star_teff_k=5772, star_luminosity_solar=1.0, star_age_gyr=4.57, planet_mass_earth=1.0, planet_radius_earth=1.0, semi_major_axis_au=1.0, orbital_period_days=365.25, eccentricity=0.0167, insolation_earth=1.0, equilibrium_temp_k=255, rotation_period_days=1.0, atmosphere_score=1.0, magnetosphere_score=1.0, surface_state_score=1.0, has_large_moon=1.0),
        dict(name="Venus", star_class="G", star_mass_solar=1.0, star_radius_solar=1.0, star_teff_k=5772, star_luminosity_solar=1.0, star_age_gyr=4.57, planet_mass_earth=0.815, planet_radius_earth=0.949, semi_major_axis_au=0.723, orbital_period_days=224.7, eccentricity=0.0068, insolation_earth=1.91, rotation_period_days=243, atmosphere_score=0.35, magnetosphere_score=0.05, surface_state_score=0.10, has_large_moon=0.0),
        dict(name="Mars", star_class="G", star_mass_solar=1.0, star_radius_solar=1.0, star_teff_k=5772, star_luminosity_solar=1.0, star_age_gyr=4.57, planet_mass_earth=0.107, planet_radius_earth=0.532, semi_major_axis_au=1.524, orbital_period_days=687, eccentricity=0.0934, insolation_earth=0.43, rotation_period_days=1.03, atmosphere_score=0.25, magnetosphere_score=0.05, surface_state_score=0.20, has_large_moon=0.0),

        # G/K benchmark candidates
        dict(name="Kepler-452 b", star_class="G", star_mass_solar=1.04, star_radius_solar=1.11, star_teff_k=5757, star_luminosity_solar=1.20, star_age_gyr=6.0, planet_mass_earth=5.0, planet_radius_earth=1.63, semi_major_axis_au=1.05, orbital_period_days=384.8, eccentricity=np.nan, insolation_earth=1.10),
        dict(name="Kepler-442 b", star_class="K", star_mass_solar=0.61, star_radius_solar=0.60, star_teff_k=4400, star_luminosity_solar=0.12, star_age_gyr=3.0, planet_mass_earth=2.3, planet_radius_earth=1.34, semi_major_axis_au=0.409, orbital_period_days=112.3, eccentricity=0.04, insolation_earth=0.70),
        dict(name="Kepler-62 e", star_class="K", star_mass_solar=0.69, star_radius_solar=0.64, star_teff_k=4925, star_luminosity_solar=0.21, star_age_gyr=7.0, planet_mass_earth=4.5, planet_radius_earth=1.61, semi_major_axis_au=0.427, orbital_period_days=122.4, eccentricity=np.nan, insolation_earth=1.20),
        dict(name="Kepler-62 f", star_class="K", star_mass_solar=0.69, star_radius_solar=0.64, star_teff_k=4925, star_luminosity_solar=0.21, star_age_gyr=7.0, planet_mass_earth=2.8, planet_radius_earth=1.41, semi_major_axis_au=0.718, orbital_period_days=267.3, eccentricity=np.nan, insolation_earth=0.41),

        # M-dwarf / compact HZ candidates
        dict(name="TOI-700 d", star_class="M", star_mass_solar=0.416, star_radius_solar=0.42, star_teff_k=3480, star_luminosity_solar=0.023, star_age_gyr=1.5, planet_mass_earth=1.7, planet_radius_earth=1.14, semi_major_axis_au=0.163, orbital_period_days=37.4, eccentricity=0.03, insolation_earth=0.86, rotation_period_days=37.4),
        dict(name="TOI-700 e", star_class="M", star_mass_solar=0.416, star_radius_solar=0.42, star_teff_k=3480, star_luminosity_solar=0.023, star_age_gyr=1.5, planet_mass_earth=0.95, planet_radius_earth=0.95, semi_major_axis_au=0.134, orbital_period_days=27.8, eccentricity=0.05, insolation_earth=1.27, rotation_period_days=27.8),
        dict(name="Kepler-186 f", star_class="M", star_mass_solar=0.54, star_radius_solar=0.52, star_teff_k=3755, star_luminosity_solar=0.041, star_age_gyr=4.0, planet_mass_earth=1.4, planet_radius_earth=1.17, semi_major_axis_au=0.43, orbital_period_days=129.9, eccentricity=np.nan, insolation_earth=0.32),
        dict(name="Proxima Centauri b", star_class="M", star_mass_solar=0.122, star_radius_solar=0.154, star_teff_k=3042, star_luminosity_solar=0.00155, star_age_gyr=4.8, planet_mass_earth=1.27, planet_radius_earth=np.nan, semi_major_axis_au=0.0485, orbital_period_days=11.19, eccentricity=0.02, insolation_earth=0.65, rotation_period_days=11.19, stellar_activity_raw=2.3),
        dict(name="TRAPPIST-1 e", star_class="M", star_mass_solar=0.089, star_radius_solar=0.12, star_teff_k=2566, star_luminosity_solar=0.00055, star_age_gyr=7.6, planet_mass_earth=0.69, planet_radius_earth=0.92, semi_major_axis_au=0.029, orbital_period_days=6.10, eccentricity=0.005, insolation_earth=0.66, rotation_period_days=6.10),
        dict(name="TRAPPIST-1 f", star_class="M", star_mass_solar=0.089, star_radius_solar=0.12, star_teff_k=2566, star_luminosity_solar=0.00055, star_age_gyr=7.6, planet_mass_earth=1.04, planet_radius_earth=1.05, semi_major_axis_au=0.0385, orbital_period_days=9.21, eccentricity=0.01, insolation_earth=0.38, rotation_period_days=9.21),
        dict(name="TRAPPIST-1 g", star_class="M", star_mass_solar=0.089, star_radius_solar=0.12, star_teff_k=2566, star_luminosity_solar=0.00055, star_age_gyr=7.6, planet_mass_earth=1.32, planet_radius_earth=1.13, semi_major_axis_au=0.0469, orbital_period_days=12.35, eccentricity=0.003, insolation_earth=0.26, rotation_period_days=12.35),
        dict(name="LHS 1140 b", star_class="M", star_mass_solar=0.18, star_radius_solar=0.21, star_teff_k=3131, star_luminosity_solar=0.0044, star_age_gyr=5.0, planet_mass_earth=5.6, planet_radius_earth=1.73, semi_major_axis_au=0.094, orbital_period_days=24.7, eccentricity=0.04, insolation_earth=0.46, rotation_period_days=24.7),
        dict(name="Ross 128 b", star_class="M", star_mass_solar=0.17, star_radius_solar=0.20, star_teff_k=3192, star_luminosity_solar=0.0036, star_age_gyr=5.0, planet_mass_earth=1.4, planet_radius_earth=np.nan, semi_major_axis_au=0.049, orbital_period_days=9.87, eccentricity=0.04, insolation_earth=1.38, rotation_period_days=9.87),
        dict(name="Wolf 1061 c", star_class="M", star_mass_solar=0.29, star_radius_solar=0.31, star_teff_k=3342, star_luminosity_solar=0.0079, star_age_gyr=5.0, planet_mass_earth=4.3, planet_radius_earth=1.6, semi_major_axis_au=0.089, orbital_period_days=17.9, eccentricity=0.11, insolation_earth=1.0, rotation_period_days=17.9),
        dict(name="GJ 1061 d", star_class="M", star_mass_solar=0.12, star_radius_solar=0.16, star_teff_k=2953, star_luminosity_solar=0.0017, star_age_gyr=7.0, planet_mass_earth=1.64, planet_radius_earth=np.nan, semi_major_axis_au=0.054, orbital_period_days=13.0, eccentricity=0.05, insolation_earth=0.58, rotation_period_days=13.0),
        dict(name="Teegarden b", star_class="M", star_mass_solar=0.089, star_radius_solar=0.107, star_teff_k=2904, star_luminosity_solar=0.00073, star_age_gyr=8.0, planet_mass_earth=1.05, planet_radius_earth=np.nan, semi_major_axis_au=0.0252, orbital_period_days=4.91, eccentricity=0.0, insolation_earth=1.15, rotation_period_days=4.91),
        dict(name="Teegarden c", star_class="M", star_mass_solar=0.089, star_radius_solar=0.107, star_teff_k=2904, star_luminosity_solar=0.00073, star_age_gyr=8.0, planet_mass_earth=1.11, planet_radius_earth=np.nan, semi_major_axis_au=0.0443, orbital_period_days=11.4, eccentricity=0.0, insolation_earth=0.37, rotation_period_days=11.4),

        # False positive / boundary-layer control
        dict(name="K2-18 b", star_class="M", star_mass_solar=0.36, star_radius_solar=0.41, star_teff_k=3457, star_luminosity_solar=0.023, star_age_gyr=2.4, planet_mass_earth=8.6, planet_radius_earth=2.6, semi_major_axis_au=0.143, orbital_period_days=33.0, eccentricity=0.20, insolation_earth=1.0, rotation_period_days=33.0, atmosphere_score=0.45, magnetosphere_score=0.35, surface_state_score=0.20),
    ])


# ============================================================
# Plotting
# ============================================================

def plot_outputs(ranked: pd.DataFrame, outdir: str, top: int) -> None:
    os.makedirs(outdir, exist_ok=True)
    top_df = ranked.head(top).iloc[::-1]

    plt.figure(figsize=(11, max(6, 0.36 * len(top_df))))
    plt.barh(top_df["name"].astype(str), top_df["ECEM_flux"], label="ECEM_flux")
    plt.barh(top_df["name"].astype(str), top_df["ECEM_base"], alpha=0.35, label="ECEM_base")
    plt.xlabel("Score")
    plt.ylabel("Planet/system")
    plt.title("Milky Way ECEM × flux-history ranking")
    plt.xlim(0, 1.05)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "milky_way_ecem_flux_scores.png"), dpi=180)
    plt.close()

    plt.figure(figsize=(9, 6))
    x = pd.to_numeric(ranked["tidal_strength_earth"], errors="coerce")
    y = pd.to_numeric(ranked["flux_tension_proxy"], errors="coerce")
    c = pd.to_numeric(ranked["ECEM_flux"], errors="coerce")
    plt.scatter(x, y, s=70 + 260 * c.fillna(0), c=c)
    for _, r in ranked.head(top).iterrows():
        rx = safe_float(r.get("tidal_strength_earth"))
        ry = safe_float(r.get("flux_tension_proxy"))
        if rx > 0 and ry > 0:
            plt.text(rx, ry, str(r["name"]), fontsize=8)
    plt.scatter([1.0], [1.0], marker="*", s=220, label="Earth reference")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Tidal strength proxy [Earth = 1]")
    plt.ylabel("Flux tension proxy [Earth ≈ 1]")
    plt.title("Flux-tension vs tidal-compression coherence space")
    plt.colorbar(label="ECEM_flux")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "milky_way_ecem_tidal_vs_flux.png"), dpi=180)
    plt.close()

    plt.figure(figsize=(8, 5))
    class_order = ["F", "G", "K", "M", "UNKNOWN"]
    data = []
    labels = []
    for cls in class_order:
        vals = ranked.loc[ranked["star_class_inferred"] == cls, "ECEM_flux"].astype(float).values
        if len(vals):
            data.append(vals)
            labels.append(cls)
    if data:
        plt.boxplot(data, tick_labels=labels, showmeans=True)
    plt.ylabel("ECEM_flux")
    plt.xlabel("Stellar class")
    plt.title("Flux-adjusted ECEM by stellar class")
    plt.ylim(0, 1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "milky_way_ecem_starclass.png"), dpi=180)
    plt.close()

    plt.figure(figsize=(8, 6))
    x = pd.to_numeric(ranked["planet_radius_earth"], errors="coerce")
    y = pd.to_numeric(ranked["planet_mass_earth"], errors="coerce")
    c = pd.to_numeric(ranked["ECEM_flux"], errors="coerce")
    plt.scatter(x, y, s=70 + 260 * c.fillna(0), c=c)
    for _, r in ranked.head(top).iterrows():
        rx = safe_float(r.get("planet_radius_earth"))
        ry = safe_float(r.get("planet_mass_earth"))
        if math.isfinite(rx) and math.isfinite(ry):
            plt.text(rx, ry, str(r["name"]), fontsize=8)
    plt.scatter([1.0], [1.0], marker="*", s=220, label="Earth reference")
    plt.xlabel("Radius [Earth radii]")
    plt.ylabel("Mass [Earth masses]")
    plt.title("Mass-radius flux-coherence space")
    plt.colorbar(label="ECEM_flux")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "milky_way_ecem_mass_radius.png"), dpi=180)
    plt.close()

    plt.figure(figsize=(8, 6))
    x = pd.to_numeric(ranked["semi_major_axis_au_inferred"], errors="coerce")
    y = pd.to_numeric(ranked["insolation_earth_inferred"], errors="coerce")
    c = pd.to_numeric(ranked["ECEM_flux"], errors="coerce")
    plt.scatter(x, y, s=70 + 260 * c.fillna(0), c=c)
    for _, r in ranked.head(top).iterrows():
        rx = safe_float(r.get("semi_major_axis_au_inferred"))
        ry = safe_float(r.get("insolation_earth_inferred"))
        if rx > 0 and ry > 0:
            plt.text(rx, ry, str(r["name"]), fontsize=8)
    plt.scatter([1.0], [1.0], marker="*", s=220, label="Earth reference")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Semi-major axis [AU]")
    plt.ylabel("Insolation [Earth = 1]")
    plt.title("Flux-orbit ECEM space")
    plt.colorbar(label="ECEM_flux")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "milky_way_ecem_flux_orbit.png"), dpi=180)
    plt.close()


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Milky Way ECEM × flux-dynamics sweep.")
    parser.add_argument("--input", type=str, default=None, help="Optional input exoplanet CSV.")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory.")
    parser.add_argument("--top", type=int, default=25, help="Number of top systems to print/plot.")
    parser.add_argument("--no-plots", action="store_true", help="Disable plot generation.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.input:
        df = pd.read_csv(args.input)
    else:
        df = demo_catalog()

    ranked = score_catalog(df, BaseWeights(), FluxWeights(), Widths())
    out_csv = os.path.join(args.outdir, "milky_way_ecem_flux_ranked.csv")
    ranked.to_csv(out_csv, index=False)

    display_cols = [
        "rank_flux", "rank_base", "name", "star_class_inferred",
        "ECEM_flux", "ECEM_base", "flux_adjustment", "flux_tension_proxy",
        "tidal_strength_earth", "score_surface_state", "score_flux_history",
        "classification",
    ]
    display_cols = [c for c in display_cols if c in ranked.columns]

    print("\nMilky Way ECEM × flux-dynamics ranking")
    print("=" * 72)
    with pd.option_context("display.max_colwidth", 42, "display.width", 180):
        print(ranked[display_cols].head(args.top).to_string(index=False))
    print(f"\nSaved ranked catalog: {out_csv}")

    if not args.no_plots:
        plot_outputs(ranked, args.outdir, args.top)
        print("Saved plots:")
        for fname in [
            "milky_way_ecem_flux_scores.png",
            "milky_way_ecem_tidal_vs_flux.png",
            "milky_way_ecem_starclass.png",
            "milky_way_ecem_mass_radius.png",
            "milky_way_ecem_flux_orbit.png",
        ]:
            print(f"  {os.path.join(args.outdir, fname)}")


if __name__ == "__main__":
    main()
