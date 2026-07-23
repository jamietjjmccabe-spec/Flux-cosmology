"""Reference utilities for the 23 July 2026 MIP laboratory framework.

This module is intentionally narrow. It does not implement cosmology, variable
constants, metric-memory gravity, or a universal Planck-scale update clock.
It provides the four frozen objects used by the current candidate model:

1. environmental distinguishability m(t), represented by trace distance;
2. positive record-production source s(t) = max(0, dm/dt);
3. finite-memory registration load q(t);
4. a bounded stochastic actualisation rate gamma_MIP(q).

The stochastic trajectory helper is a numerical research scaffold, not evidence
that one particular quantum-state unravelling is ontically real.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class RegistrationParameters:
    """Parameters of the finite-memory, bounded-rate registration model."""

    tau_r: float
    gamma_0: float = 0.0
    gamma_max: float = 1.0
    q_star: float = 1.0
    power: float = 2.0

    def __post_init__(self) -> None:
        if self.tau_r <= 0:
            raise ValueError("tau_r must be positive")
        if self.gamma_0 < 0 or self.gamma_max < 0:
            raise ValueError("gamma values must be non-negative")
        if self.q_star <= 0:
            raise ValueError("q_star must be positive")
        if self.power <= 0:
            raise ValueError("power must be positive")


def _as_square_complex(matrix: Array, name: str) -> Array:
    value = np.asarray(matrix, dtype=np.complex128)
    if value.ndim != 2 or value.shape[0] != value.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} contains non-finite values")
    return value


def _is_hermitian(matrix: Array, atol: float = 1e-10) -> bool:
    return bool(np.allclose(matrix, matrix.conj().T, atol=atol, rtol=0.0))


def trace_distance(rho_a: Array, rho_b: Array, atol: float = 1e-10) -> float:
    """Return D_tr(rho_a, rho_b) = 1/2 ||rho_a-rho_b||_1.

    Inputs are checked for square shape and Hermiticity. The function does not
    silently renormalise states, because missing trace is itself an audit error.
    """

    a = _as_square_complex(rho_a, "rho_a")
    b = _as_square_complex(rho_b, "rho_b")
    if a.shape != b.shape:
        raise ValueError("rho_a and rho_b must have the same shape")
    if not _is_hermitian(a, atol) or not _is_hermitian(b, atol):
        raise ValueError("density matrices must be Hermitian")
    if not np.isclose(np.trace(a).real, 1.0, atol=atol, rtol=0.0):
        raise ValueError("rho_a must have unit trace")
    if not np.isclose(np.trace(b).real, 1.0, atol=atol, rtol=0.0):
        raise ValueError("rho_b must have unit trace")

    delta = 0.5 * ((a - b) + (a - b).conj().T)
    eigenvalues = np.linalg.eigvalsh(delta)
    distance = 0.5 * float(np.sum(np.abs(eigenvalues)))
    return float(np.clip(distance, 0.0, 1.0))


def positive_record_source(times: Array, distinguishability: Array) -> Array:
    """Calculate s(t)=max(0, dm/dt) on an arbitrary increasing time grid."""

    t = np.asarray(times, dtype=float)
    m = np.asarray(distinguishability, dtype=float)
    if t.ndim != 1 or m.ndim != 1 or t.shape != m.shape:
        raise ValueError("times and distinguishability must be equal-length vectors")
    if len(t) < 2:
        raise ValueError("at least two samples are required")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(m)):
        raise ValueError("inputs must be finite")
    if np.any(np.diff(t) <= 0):
        raise ValueError("times must be strictly increasing")
    if np.any((m < -1e-12) | (m > 1.0 + 1e-12)):
        raise ValueError("trace-distance samples must lie in [0, 1]")

    derivative = np.gradient(m, t, edge_order=1)
    return np.maximum(derivative, 0.0)


def integrate_registration_load(
    times: Array,
    source: Array,
    tau_r: float,
    q0: float = 0.0,
) -> Array:
    """Integrate tau_R dq/dt = -q + s(t) with an exact frozen-source step.

    Over each interval the source is approximated by its endpoint average, while
    the linear relaxation is integrated exactly. This preserves q >= 0 for a
    non-negative source and initial condition.
    """

    t = np.asarray(times, dtype=float)
    s = np.asarray(source, dtype=float)
    if t.ndim != 1 or s.ndim != 1 or t.shape != s.shape:
        raise ValueError("times and source must be equal-length vectors")
    if len(t) < 2:
        raise ValueError("at least two samples are required")
    if tau_r <= 0:
        raise ValueError("tau_r must be positive")
    if q0 < 0:
        raise ValueError("q0 must be non-negative")
    if np.any(np.diff(t) <= 0):
        raise ValueError("times must be strictly increasing")
    if np.any(~np.isfinite(s)) or np.any(s < 0):
        raise ValueError("source must be finite and non-negative")

    q = np.empty_like(s)
    q[0] = q0
    for index, dt in enumerate(np.diff(t), start=1):
        decay = np.exp(-dt / tau_r)
        source_mid = 0.5 * (s[index - 1] + s[index])
        q[index] = q[index - 1] * decay + source_mid * (1.0 - decay)
    return q


def bounded_actualisation_rate(q: Array | float, params: RegistrationParameters) -> Array:
    """Return gamma_0 + gamma_max q^p/(q_star^p+q^p)."""

    load = np.asarray(q, dtype=float)
    if np.any(~np.isfinite(load)) or np.any(load < 0):
        raise ValueError("q must be finite and non-negative")
    numerator = np.power(load, params.power)
    denominator = np.power(params.q_star, params.power) + numerator
    return params.gamma_0 + params.gamma_max * numerator / denominator


def gaussian_spatial_kernel(distances: Array, xi: float) -> Array:
    """Return a normalised isotropic Gaussian kernel sampled at given radii.

    This helper supplies only weights. Experimental use must include the actual
    apparatus geometry, electromagnetic, phononic, thermal, and common-mode
    nuisance propagation.
    """

    r = np.asarray(distances, dtype=float)
    if xi <= 0:
        raise ValueError("xi must be positive")
    if np.any(~np.isfinite(r)) or np.any(r < 0):
        raise ValueError("distances must be finite and non-negative")
    weights = np.exp(-0.5 * (r / xi) ** 2)
    total = float(np.sum(weights))
    if total <= 0:
        raise ValueError("kernel underflowed to zero")
    return weights / total


def stochastic_actualisation_step(
    psi: Array,
    hamiltonian: Array,
    observable: Array,
    gamma: float,
    dt: float,
    d_wiener: Optional[float] = None,
    rng: Optional[np.random.Generator] = None,
) -> Array:
    """Take one normalised Euler-Maruyama step of a diffusive trajectory.

    The convention is

        d|psi> = [-i H dt - gamma/2 (A-<A>)^2 dt
                  + sqrt(gamma) (A-<A>) dW] |psi>.

    This is a compact scaffold for tests and small prototypes. Production audits
    should use a completely-positive finite-step trajectory method and explicit
    timestep convergence checks.
    """

    state = np.asarray(psi, dtype=np.complex128)
    if state.ndim != 1:
        raise ValueError("psi must be a state vector")
    h = _as_square_complex(hamiltonian, "hamiltonian")
    a = _as_square_complex(observable, "observable")
    if h.shape != a.shape or h.shape[0] != state.shape[0]:
        raise ValueError("state and operators have incompatible dimensions")
    if not _is_hermitian(h) or not _is_hermitian(a):
        raise ValueError("hamiltonian and observable must be Hermitian")
    if gamma < 0 or dt <= 0:
        raise ValueError("gamma must be non-negative and dt must be positive")

    norm = np.linalg.norm(state)
    if not np.isfinite(norm) or norm == 0:
        raise ValueError("psi must have finite non-zero norm")
    state = state / norm

    if d_wiener is None:
        generator = rng if rng is not None else np.random.default_rng()
        d_wiener = float(generator.normal(0.0, np.sqrt(dt)))
    if not np.isfinite(d_wiener):
        raise ValueError("d_wiener must be finite")

    expectation = float(np.vdot(state, a @ state).real)
    centered = a - expectation * np.eye(a.shape[0], dtype=np.complex128)
    drift = -1j * (h @ state) - 0.5 * gamma * (centered @ centered @ state)
    diffusion = np.sqrt(gamma) * (centered @ state)
    updated = state + drift * dt + diffusion * d_wiener
    updated_norm = np.linalg.norm(updated)
    if not np.isfinite(updated_norm) or updated_norm == 0:
        raise FloatingPointError("trajectory step produced an invalid state")
    return updated / updated_norm
