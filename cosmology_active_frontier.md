# Cosmology Active Frontier

## Status

**Active frontier:** cosmological flux sector  
**Date of transition:** 2026-03-25  
**Reason:** the scalar galactic bridge chapter reached a documented architectural ceiling, while the cosmological branch remains structurally live and has already produced interpretable benchmark behavior.

This file defines the current active research frontier for Flux-Cosmology after closure of the scalar galactic chapter.

---

## Why cosmology is now the live branch

The galactic scalar chapter established a bounded success regime and a hard limit:

- meaningful predictive success in the MW-band
- non-generalization across the full 33-galaxy validation set
- ceiling reached for the one-sector scalar galactic bridge with universal parameters

That result does **not** close the broader framework. The cosmological branch is structurally independent enough that it remains a legitimate live research path.

The cosmological sector already has:

- a functioning background solver
- a functioning linear growth solver
- GR-limit recovery behavior
- a weak-deviation branch and a visible modified-growth branch
- a benchmark comparison structure (`with xi`, `without xi`, `GR limit`)

This makes cosmology the correct next active frontier.

---

## Current cosmology baseline

The present cosmology work is built around a flux/scalar background with linear growth evolution on top of it.

### Core variables

- \(N = \ln a\)
- \(a\) = scale factor
- \(z = a^{-1} - 1\)
- scalar / flux field \(\sigma(N)\)
- Hubble rate \(H(N)\)
- effective matter fraction \(\Omega_m(N)\)
- effective flux/scalar fraction \(\Omega_\sigma(N)\)
- effective gravity modifier \(\mu_{\rm eff}(N)\)

### Current potential ansatz

The active toy potential is of the form

\[
V(\sigma) = V_0\left(1 - e^{-\sigma/\mu_V}\right)^2
\]

with benchmark-level parameter choices of the type:

- `Omega_m0 ≈ 0.30`
- `Omega_r0 ≈ 9e-5`
- `V0 ≈ 0.70`
- `mu_V ≈ 2.0`
- `sigma_init ≈ 5.0`
- `sigma_prime_init ≈ 0.0`

These are working benchmark values, not final fundamental constants.

---

## Current benchmark branches

The cosmology codebase currently distinguishes three important comparison branches.

### Branch A — Flux + xi

- `beta_mu > 0`
- `use_xi_factor = True`

This branch gives a near-GR growth response because the coherence factor suppresses \(\mu_{\rm eff} - 1\) almost completely.

Interpretation:
- useful as a GR-recovery / safe baseline
- physically conservative
- currently too suppressed to generate a strong linear-growth signal

### Branch B — Flux, no xi

- `beta_mu > 0`
- `use_xi_factor = False`

This branch produces a visible modified-growth signal.

Interpretation:
- useful as the first clearly testable modified-growth branch
- demonstrates that the flux sector can affect linear structure growth when the coherence suppression is removed
- currently functions as the main “visible deviation” branch

### Branch C — GR limit

- `beta_mu = 0`

Interpretation:
- control branch
- validates solver behavior
- allows direct comparison between flux and GR-limit evolution

---

## What has already been established

The following points are considered established at the benchmark level.

### 1. The background solver runs and remains stable

The scalar/flux background evolution can be integrated successfully over the chosen range in \(N\), and the reconstructed quantities \(H\), \(\Omega_m\), and \(\Omega_\sigma\) behave smoothly in the tested benchmark runs.

### 2. The growth solver runs and produces interpretable output

The linear growth equation produces sensible values for:

- \(D(a)\)
- \(f = d\ln D / d\ln a\)
- \(f\sigma_8(z)\)

The previous plotting and derivative issues were resolved to the point that the current growth behavior is interpretable.

### 3. The apparent “hump” in \(f\sigma_8\) is not a numerical failure

The diagnostic work showed that the low-redshift shape in the benchmark run was compatible with the product structure

\[
f\sigma_8(z) = f(z)\,D(z)\,\sigma_{8,0}
\]

rather than being an endpoint-artifact or broken derivative extraction.

### 4. The xi-weighted branch self-silences

The coherence factor currently used in the growth-sector modification suppresses \(\mu_{\rm eff} - 1\) to the level of roughly a few \(\times 10^{-5}\) in benchmark runs.

Interpretation:
- mathematically bounded and stable
- phenomenologically too weak to produce a strong linear-growth deviation

### 5. Removing xi produces a real modified-growth signal

The no-xi branch gives a clear separation from the GR-limit branch in

- \(f(z)\)
- \(f\sigma_8(z)\)
- \(\mu_{\rm eff}(z) - 1\)

Interpretation:
- the growth-sector machinery is capable of producing a visible cosmological signature
- the current question is not “can the framework deviate?” but “what is the physically correct bounded form?”

---

## The central open problem

The central cosmology problem is now clear:

**What is the correct shared cosmological coherence / flux variable and bounded coupling law that preserves stability while avoiding complete self-silencing?**

This problem mirrors the broader theory issue identified earlier: different sectors still use related but not fully unified proxies for the same underlying quantity.

At cosmology level, the active practical version of that problem is:

- the present `xi` form is too suppressive
- removing `xi` entirely gives visible deviations but may be too unconstrained
- therefore the next progress depends on finding a better bounded coupling form

---

## What remains live in the cosmology branch

The following are still active research targets.

### 1. Coherence-factor redesign

The current `xi` factor has not been promoted to canon. It remains a benchmark suppressor only.

The main live question is whether a softer bounded factor can preserve:

- GR recovery in frozen-field regimes
- numerical stability
- visible but not excessive growth deviation

Candidate families include softened bounded forms such as:

- rational forms
- tanh-like forms
- bounded exponential saturation forms

These are still open and not yet benchmark-frozen.

### 2. Better normalization / interpretation of \(\Omega_m\) and \(\Omega_\sigma\)

The bookkeeping of the growth source term versus the Hubble-normalized energy content has already been identified as an area where structural clarity matters.

This is still a live interpretation / implementation checkpoint rather than a settled canonical point.

### 3. Cosmological-to-observational comparison

The code is now ready for more direct comparison against observationally relevant growth behavior. This includes:

- broader \(f\sigma_8(z)\) comparison
- low-redshift growth trends
- benchmark comparison to GR-limit behavior over the same solver stack

### 4. Lensing / refractive extension

A future extension may build a flux-lensing or refractive-actualization module using the same field logic, but this is not yet the active benchmark target. It remains a downstream branch, contingent on a cleaner cosmological field/coupling definition.

---

## What is currently *not* active

The following are not active frontier tasks right now.

### 1. Reopening the scalar galactic benchmark as if it were unresolved

That chapter is closed in its present scalar form. Any future galactic reopening must begin as a new architecture chapter, not as a continuation of the old scalar benchmark search.

### 2. Additional ad hoc scalar tweaks without benchmark discipline

Cosmology work must avoid the same uncontrolled branch proliferation that was necessary in the galactic chapter. New changes must be explicitly tagged as:

- numerical
- structural
- physical / interpretive
- benchmarking only

### 3. Grand-claim unification without a weak-field benchmark

The framework should not jump directly to claims about replacing dark matter in lensing or strong-field quasar optics until the weak-field cosmological / refractive benchmark structure is explicitly written and tested.

---

## Immediate next tests

The next active tests should proceed in this order.

### Test 1 — Freeze the current cosmology benchmark outputs

Create or preserve a benchmark artifact that records:

- Flux + xi
- Flux, no xi
- GR limit
- present-day values of \(f(z=0)\), \(f\sigma_8(z=0)\), and \(\mu_{\rm eff}(z=0)-1\)
- the benchmark comparison plots

This prevents future ambiguity about what the current cosmology branch has already established.

### Test 2 — Rescan the coherence-factor family

Hold the background and solver structure fixed, and test alternative bounded forms replacing the current `xi` suppression.

This is the highest-value next research step in cosmology.

### Test 3 — Quantify the visible-deviation window

Determine whether there is a bounded-coupling regime that gives:

- visible separation from GR
- no obvious instability
- no grossly excessive late-time enhancement

### Test 4 — Prepare a weak-field flux-lensing benchmark

Only after the coupling law is cleaner should the refractive / lensing branch be formalized into a dedicated module.

---

## Falsifiability targets

The cosmology branch remains valuable only if it stays tied to explicit failure conditions.

Examples of meaningful falsifiability targets include:

- the bounded coupling law always self-silences or always blows up, with no viable intermediate regime
- the flux branch cannot produce a consistent low-redshift growth signature distinguishable from GR while remaining stable
- the same field logic cannot be extended even approximately into a weak-field refractive-lensing benchmark

If any of these become unavoidable, the cosmology branch must be revised structurally, not cosmetically.

---

## Current interpretation discipline

At this stage, the correct interpretation is:

- the cosmological flux sector is **still live**
- it has **not** yet reached the architectural ceiling that closed the scalar galactic chapter
- the present benchmark already shows both a conservative GR-like branch and a visible modified-growth branch
- the next progress depends on identifying the correct bounded coupling, not on inventing a new galaxy fit family

---

## Working summary

The scalar galactic chapter is closed.  
The cosmological branch is now the active frontier.

The current best reading is:

- Flux-Cosmology has a functioning cosmological benchmark stack
- the growth sector can produce both GR-like and visibly modified branches
- the principal open problem is the correct bounded coherence / coupling law
- the next frontier is therefore **cosmological coupling refinement and observational comparison**, not further scalar galactic tuning

