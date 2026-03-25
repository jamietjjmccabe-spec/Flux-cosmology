# SPARC Subsample Protocol

## Purpose

This protocol defines the first frozen-parameter SPARC subsample test for the flux-halo galactic bridge. Its purpose is to test whether **Closure D**, with parameters fixed in advance, can reproduce the qualitative and semi-quantitative structure of real rotation curves across a small but deliberately diverse galaxy set.

This is **not** a refit stage. It is a controlled benchmark stage.

The goal is to answer one question clearly:

**Does the frozen Closure D law survive first contact with real galaxies without per-galaxy retuning?**

---

## Frozen model under test

### Closure D benchmark law

The galactic flux sector is defined by a baryon-seeded, nonlocal sourced potential with coherence scale

\[
L_\chi = 18\,R_b\left(\frac{\Sigma_b}{\Sigma_{\rm ref}}\right)^{-0.08}
\]

where:

- \(R_b\) is the baryonic structural scale used in the benchmark
- \(\Sigma_b\) is the baryonic surface density proxy
- \(\Sigma_{\rm ref}\) is the benchmark reference surface density

### Frozen benchmark parameters

Use the currently frozen benchmark values exactly:

- `k = 18`
- `gamma = 0.08`
- `eta = 12`
- `n = 0.7`

No per-galaxy tuning is allowed in this subsample stage.

### Frozen interpretation

The current working interpretation is:

- the coherence scale is **primarily geometric**
- the surface-density term is a **weak but real correction**
- the flux sector acts through **nonlocal sourced potential accumulation**, not through a local multiplier on baryonic acceleration

---

## Scope of this subsample

The subsample is intentionally small. It is not intended to prove the full theory. It is intended to stress-test the frozen benchmark across a representative range of galaxy structure.

Target size:

- **6 to 10 SPARC galaxies**

Selection must span the following structural regimes as cleanly as possible:

1. **Diffuse dwarf / low-surface-brightness system**
2. **MW-like spiral**
3. **High-surface-brightness compact system**
4. **Bulgeless or weak-bulge late-type disk**
5. **Gas-rich low-surface-density disk**
6. **At least one massive, extended high-mass system**

The point is not sample size. The point is structural diversity.

---

## Galaxy selection rules

A galaxy may be included only if all of the following are true:

### Data quality

- SPARC rotation curve is sufficiently sampled over inner and outer radii
- Baryonic decomposition is available and usable
- Inclination and distance uncertainties are not extreme relative to the rest of the candidate list

### Structural usefulness

Each selected galaxy must help stress at least one axis relevant to Closure D:

- baryonic size
- baryonic surface density
- bulge dominance vs disk dominance
- gas dominance vs stellar dominance
- compactness vs diffuseness

### Exclusion rule

Do **not** choose galaxies simply because they look easy to fit.

This subsample is a stress test, not a showcase.

---

## Inputs to be used

For each selected SPARC galaxy, record and freeze the following inputs before running the model:

- galaxy name
- distance used
- inclination used
- stellar disk mass or mass proxy
- bulge mass or mass proxy, if present
- gas mass or gas contribution
- structural radius proxy used for \(R_b\)
- surface density proxy used for \(\Sigma_b\)
- observed rotation curve data
- baryonic contribution curve or reconstructable baryonic components

These must be written to a run table before any model comparison is interpreted.

---

## Structural definitions to keep fixed

Unless a benchmark revision is explicitly approved, use the same structural definitions across the full subsample.

### Baryonic scale

Use the currently benchmarked baryonic scale rule:

- `R_b = max(R_disk, a_bulge)`

If a galaxy lacks a bulge, this reduces naturally to the disk scale.

### Surface density proxy

Use the same surface-density construction across the sample. Do not change definitions galaxy by galaxy.

If a later revision is needed, it must be applied to the **entire subsample** and logged as a benchmark revision.

---

## Metrics to compute

For every galaxy, compute the following model diagnostics.

### 1. Outer-slope metric

Primary shape metric:

- logarithmic outer slope of the model rotation curve
- compare against observed outer behavior

This remains the first-pass discriminator for whether the model produces halo-like persistence.

### 2. Residual structure

Compute residuals between model and observed rotation curve:

- inner region
- transition region
- outer region

Look not only at RMS mismatch, but at whether the residual pattern is structured or random.

### 3. Global fit score

For each galaxy, store a scalar summary such as:

- RMS residual
- mean absolute residual
- optionally weighted chi-like score if uncertainty handling is stable

Use the same score definition for all galaxies.

### 4. Flux-scale sanity checks

Track whether the inferred coherence scale and flux sector remain physically sensible:

- \(L_\chi\) relative to baryonic scale
- \(L_\chi\) relative to halo-scale proxy if available
- total flux mass implied by the sourced potential
- whether the outer tail is smooth and halo-like

### 5. Failure notes

Each galaxy run must also include a short plain-language note:

- pass
- tension
- fail
- likely reason

This note is mandatory.

---

## Pass / tension / fail criteria

### Pass

A galaxy is classified as **pass** if all of the following hold:

- outer curve shape is qualitatively correct
- no obvious pathological rise or collapse in the outer rotation curve
- residuals are moderate and not strongly structured
- no per-galaxy tuning was introduced
- coherence scale and implied flux behavior remain physically interpretable

### Tension

A galaxy is classified as **tension** if:

- the model gets the broad shape right but misses amplitude or transition details
- or residuals are systematic but not catastrophic
- or a scoring metric is borderline while the physical curve still looks plausible

Tension does **not** justify parameter retuning at this stage.

### Fail

A galaxy is classified as **fail** if any of the following occur:

- outer curve shape is clearly wrong
- model becomes strongly rising or strongly Keplerian where the data are not
- residual structure is large and systematic across most radii
- model requires hidden per-galaxy adjustments to appear viable
- coherence scale or sourced potential becomes physically unreasonable

---

## Rules of discipline

These rules are mandatory for the subsample phase.

### Rule 1 — no parameter drift

Do not alter `k`, `gamma`, `eta`, or `n` during the subsample test.

### Rule 2 — no galaxy-by-galaxy rescue

If one galaxy fails, record the failure. Do not patch the model for that galaxy alone.

### Rule 3 — separate scoring artifacts from physics failures

If a run scores poorly, inspect whether the problem is:

- real curve mismatch
- a pathological inner component
- a penalty-function artifact
- a bad structural proxy

Do not let a scoring artifact rewrite the benchmark without diagnosis.

### Rule 4 — preserve the audit trail

Every galaxy run must record:

- inputs
- outputs
- classification
- whether the issue was physical, numerical, or scoring-related

### Rule 5 — lensing and merger constraints remain hard constraints

A rotation-curve pass does **not** promote the theory by itself.

Any successful SPARC subsample result must still be followed by:

- lensing mass-profile checks
- merger geometry / Bullet-style separation logic

---

## Output artifact requirements

The subsample stage should produce a single summary artifact with:

1. selected galaxy list
2. frozen parameter statement
3. per-galaxy plots
4. per-galaxy metrics table
5. pass / tension / fail classification
6. summary of common residual patterns
7. recommendation: promote, revise, or reject

---

## Promotion criteria after subsample

Closure D may be promoted from benchmark to broader observational testing only if the subsample shows all of the following:

- no hidden per-galaxy tuning
- majority of galaxies classified as pass or mild tension
- no consistent failure mode across one structural class
- no obvious contradiction between the galactic operator and the benchmark physical interpretation

If the subsample instead shows a class-specific failure pattern, that pattern becomes the next research target.

---

## Current status

At the moment of writing:

- Closure D is frozen as the galactic benchmark
- the fourth-archetype corrected rerun supports `gamma = 0.08`
- the next stage is the frozen-parameter SPARC subsample
- no benchmark revision is authorized until that subsample is run and logged

