# Flux Cosmology / Minimum Interface Point Research Archive

This repository is an open research archive for Flux cosmology, registration-depth phenomenology, Minimum Interface Point (MIP) concepts, numerical toy models, observational screening work, and falsifiable laboratory proposals.

## Scientific status

The repository contains work at several levels of maturity. Every result should be read according to its declared class:

- derived mathematics;
- phenomenological closure;
- toy simulation;
- heuristic analogy;
- observational fit or data-screening result;
- prospective experimental protocol.

No toy output is evidence for the complete theory.

## Active July 2026 module: Flux-MIP V2.1

The current laboratory branch is a frozen curvature-sourced Lindblad benchmark:

- regularized inverse-square finite-source precursor;
- local curvature response `S_reg = -ell^2 nabla^2 D_reg`;
- bounded susceptibility `F_reg = tanh(S_reg/S_*)`;
- linked phase, dephasing, diffusion, and heating predictions;
- conservative 7.4112 kg tungsten-plate apparatus;
- blinded four-state Phase Ladder null test.

**Benchmark signature:** `f819067c8feca48e6d0a51b790e01f5a1be599f4afa22c6accc46f3967827670`  
**Optimizer signature:** `ed3bde6ab4adcef685d3cee33942b0cc5467118dd82f86ad4d2863e2603ef6b0`

Start here:

- [`papers/flux_mip_v2_1/Flux_MIP_V2_1_Final_Paper.md`](papers/flux_mip_v2_1/Flux_MIP_V2_1_Final_Paper.md)
- [`docs/2026-07/PROJECT_STATUS_2026-07-12.md`](docs/2026-07/PROJECT_STATUS_2026-07-12.md)
- [`docs/2026-07/FLUX_MIP_V2_1_THEORY.md`](docs/2026-07/FLUX_MIP_V2_1_THEORY.md)
- [`docs/2026-07/PHASE_LADDER_PROTOCOL.md`](docs/2026-07/PHASE_LADDER_PROTOCOL.md)
- [`UPDATE_MANIFEST_2026-07-12.md`](UPDATE_MANIFEST_2026-07-12.md)

## Empirical status

Flux-MIP V2.1 has not been detected. The available Overstreet Figure 2 tables contain processed phase points and cannot independently constrain residual visibility. The dedicated apparatus is prospective.

## Earlier June 2026 update

The prior development index is preserved in [`UPDATE_MANIFEST_2026-06-30.md`](UPDATE_MANIFEST_2026-06-30.md), including MIP terminology, registration depth, metric-memory hypotheses, active cosmology boundaries, ECEE screening outputs, and toy experiments.

## Repository map

- `papers/` - manuscript sources and paper packages.
- `docs/` - theory records, status boundaries, and protocols.
- `src/` - active numerical implementations.
- `experiments/` - geometry studies, null-test tools, and toy experiments.
- `outputs/` - compact reproducible research outputs.
- `scans/` and `figures/` - observational scan products and figures from earlier modules.

## Build the V2.1 paper

```bash
cd papers/flux_mip_v2_1
latexmk -pdf -interaction=nonstopmode -halt-on-error Flux_MIP_V2_1_Final_Paper.tex
```

## Python requirements for the V2.1 pipeline

```bash
pip install numpy scipy pandas matplotlib
```

## Interpretation rule

The laboratory paper is a test specification. A positive signal would require independent replication and would not by itself establish the wider cosmological interpretation. A null result constrains the frozen V2.1 channel and may not be evaded by changing the preregistered kernel or background rule after unblinding.
