"""
flux_mip_v3q_model.py

Standalone implementation of the Flux/MIP v3q galactic archive model.
Units:
    G is in kpc (km/s)^2 / Msun.
    Radii are kpc.
    Masses are Msun.
    Velocities are km/s.
"""
from __future__ import annotations
import numpy as np

G_KPC_KMS2_PER_MSUN = 4.30091e-6
SIGMA_REF = 4.08e9

# Frozen v3o / v3q constants
Q_STRENGTH = 7.0
P_RADIUS = 4.6
S_TEXTURE = 1.36
ETA0 = 8.0
D_CRIT = 0.68
M_GATE = 8.0
ALPHA_HSB = 7.9
CHI_REM = 3.0
ZETA_REM = 0.5


def spatial_archive_depth(m_bar: float, r_disk: float, sigma_ref: float = SIGMA_REF) -> float:
    """Global baryonic archive depth D_spatial."""
    sigma_bar = m_bar / (2.0 * np.pi * r_disk**2)
    return float((sigma_bar / sigma_ref) ** 0.08)


def eta_v3o(D: np.ndarray | float) -> np.ndarray | float:
    """Saturating archival efficiency, the v3o causal sluice."""
    D = np.asarray(D)
    return (ETA0 * D ** (-Q_STRENGTH)) / (1.0 + (D / D_CRIT) ** M_GATE)


def memory_radius(D: np.ndarray | float, r_disk: np.ndarray | float) -> np.ndarray | float:
    """Historical wake radial scale L_chi."""
    return ALPHA_HSB * np.asarray(r_disk) * np.asarray(D) ** (-P_RADIUS)


def remnant_proxy(f_gas: np.ndarray | float, T: np.ndarray | float) -> np.ndarray | float:
    """Raw compact-remnant archive proxy."""
    f = np.clip(np.asarray(f_gas), 0.0, 1.0)
    t_factor = np.clip((10.0 - np.asarray(T)) / 10.0, 0.0, 1.0)
    return (1.0 - f) * t_factor


def retention_gate(vmax: np.ndarray | float, v_ref: float = 200.0, n: float = 4.0) -> np.ndarray | float:
    """Potential-depth retention gate for compact archive nodes."""
    ratio = np.maximum(np.asarray(vmax), 1e-12) / v_ref
    return ratio**n / (1.0 + ratio**n)


def eta_effective(D, f_gas, T, vmax=None, use_retention=False, v_ref=200.0, n=4.0):
    """v3q effective archival efficiency."""
    R_rem = remnant_proxy(f_gas, T)
    if use_retention:
        if vmax is None:
            raise ValueError("vmax is required when use_retention=True")
        R_eff = R_rem * retention_gate(vmax, v_ref=v_ref, n=n)
    else:
        R_eff = R_rem
    return eta_v3o(D) * (1.0 + CHI_REM * R_eff ** ZETA_REM)


def historical_velocity_squared(r, m_bar, r_disk, D, eta_eff):
    """Return V_hist^2 contribution from the historical wake."""
    r = np.asarray(r)
    L = memory_radius(D, r_disk)
    x = r / L
    M_hist = eta_eff * m_bar * (x ** S_TEXTURE) / ((1.0 + x) ** (S_TEXTURE - 1.0))
    return G_KPC_KMS2_PER_MSUN * M_hist / np.maximum(r, 1e-30)


def baryonic_velocity_squared(vgas, vdisk, vbul=0.0, ups_disk=0.5, ups_bulge=0.7):
    """SPARC baryonic velocity convention. Vgas may be signed."""
    vgas = np.asarray(vgas)
    return np.abs(vgas) * vgas + ups_disk * np.asarray(vdisk)**2 + ups_bulge * np.asarray(vbul)**2


def flux_velocity(r, vgas, vdisk, vbul, m_bar, r_disk, f_gas, T, vmax=None, use_retention=False):
    """Compute the v3q predicted circular velocity."""
    D = spatial_archive_depth(float(m_bar), float(r_disk))
    eta_eff = eta_effective(D, f_gas=f_gas, T=T, vmax=vmax, use_retention=use_retention)
    vbar2 = baryonic_velocity_squared(vgas, vdisk, vbul)
    vh2 = historical_velocity_squared(r, m_bar, r_disk, D, eta_eff)
    return np.sqrt(np.maximum(vbar2 + vh2, 0.0))