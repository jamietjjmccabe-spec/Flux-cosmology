# MIP Laboratory Experimental Programme — 23 July 2026

## Primary hypothesis

At matched deposited energy and matched conventional noise, processes that generate more distinguishable and durable environmental records produce a different residual stochastic actualisation rate.

## Frozen observables

- \(m(t)=D_{tr}[\rho_E^{(0)}(t),\rho_E^{(1)}(t)]\)
- \(s(t)=\max(0,dm/dt)\)
- \(q(t)\), governed by \(\tau_R\dot q=-q+s\)
- bounded \(\gamma_{MIP}(q)\)
- dephasing/locking observables and witness response
- standard nuisance parameters for bath memory, coupling direction, temperature, leakage, heating, and hardware relaxation

## Stage A — reversible/irreversible discrimination

Compare a coherent reversible transfer with an irreversible dump while matching energy, pulse envelope, and hardware occupancy as closely as possible.

**Pass condition:** a residual follows environmental distinguishability rather than deposited energy.

**Failure condition:** standard cQED/open-system parameters explain the difference, or no stable residual remains.

## Stage B — temporal memory

Sweep pulse spacing, pulse count, pulse shape, and record strength.

**Pass condition:** one common \(\tau_R\) and bounded \(\gamma_{MIP}(q)\) explain all temporal sweeps.

**Failure condition:** each configuration requires a different recovery time or ordinary saturation/non-Markovian bath models fit equally well.

## Stage C — spatial memory

Use a separated witness and vary distance, orientation, shielding, substrate, and common-mode geometry.

**Pass condition:** one common spatial kernel \(K_\xi(r)\) survives electromagnetic, thermal, phononic, material, and readout cross-talk controls.

**Failure condition:** the effect follows known propagation channels or the inferred \(\xi\) changes with apparatus details.

## Stage D — nuisance-space orthogonality

Perform controlled sweeps of:

- bath coupling operator and direction;
- coupling strength and bath memory time;
- temperature and occupation;
- detuning;
- resonator and qubit relaxation/dephasing times;
- reversible/irreversible pathway;
- phase-ladder orientation and geometry.

Fit the full standard model and the MIP-augmented model globally.

**Pass condition:** a single MIP parameter vector remains outside the standard nuisance tangent space across all datasets.

**Failure condition:** the residual rotates into the nuisance span or the MIP parameters drift between configurations.

## Four-state phase-ladder control

Use four blinded apparatus states chosen so that the candidate dephasing signal is even under the relevant geometric reversal while ordinary phase recalibration is odd.

The analysis must pre-register:

- state labels and blinding map;
- parity combinations;
- nuisance model and priors;
- exclusion thresholds;
- stopping rule;
- global-fit statistic;
- criteria for unblinding.

## Required simulation audits before hardware claims

1. total-excitation rather than rectangular product-Fock truncation;
2. mixed-state SME for genuinely unmonitored channels;
3. complete co-rotation of coupling, dissipator tensor, and output noises;
4. equality of trajectory distributions under symmetry-equivalent rotations;
5. timestep, cutoff, and Monte Carlo convergence;
6. record-basis invariance of unconditional dynamics;
7. synthetic-data recovery against a full structured-bath null model.

## Decision standard

No single nonlinear curve, locking event, or visibility loss is sufficient. The bridge is supported only by a common global residual that tracks record distinguishability and survives every standard-physics control above.
