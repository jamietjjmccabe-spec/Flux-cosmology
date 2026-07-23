# Framework Update Manifest — 23 July 2026

## Purpose

This update changes the repository from a broad cross-scale Flux/MIP claim set into a status-controlled research archive with one active laboratory frontier.

## Added

- `docs/2026-07/CURRENT_FRAMEWORK_2026-07-23.md`
- `docs/2026-07/BRANCH_STATUS_LEDGER_2026-07-23.md`
- `docs/2026-07/EXPERIMENTAL_PROGRAM_2026-07-23.md`
- `LEGACY_MODEL_STATUS.md`
- `src/mip_record_conditioned_actualisation.py`
- `src/__init__.py`
- `tests/test_mip_record_conditioned_actualisation.py`

## Updated

- `README.md` now points to the active laboratory framework and clearly separates active, speculative, quarantined, retired, and heuristic branches.

## Scientific changes

- MIP remains “Minimum Interface Points,” but no Planck length or universal update interval is assumed without derivation.
- Fixed MIP frame-rate models are retired.
- The active source variable is environmental trace-distance production.
- Registration uses finite temporal recovery and an optional spatial kernel.
- The nonlinear rate is bounded.
- Individual stochastic trajectories are required to discuss outcomes, while unravelling objectivity remains unresolved.
- Structured-bath degeneracy and local identifiability are both recorded.
- Gravity, dark matter, cosmology, variable constants, and cross-domain analogies are removed from the active quantum-to-classical bridge.
- ECEM/ECEE is retained as a separate heuristic project.

## Validation

The reference module was checked with:

```bash
python -m unittest discover -s tests -v
```

Six tests passed, covering trace distance, positive record production, temporal recovery, bounded rates, spatial-kernel normalisation, and trajectory-state normalisation.

## Merge policy

This update should be reviewed as a framework/status correction. It does not claim that the active MIP model has been experimentally confirmed.
