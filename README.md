# Flux/MIP Research Archive

**Current framework date:** 23 July 2026  
**Author:** Jamie McCabe  
**Status:** speculative laboratory quantum-to-classical research programme with archived cosmology, gravity, and astrophysical branches

## Start here

The active framework is no longer a single cross-scale claim linking quantum actualisation directly to dark matter, gravity, variable constants, and cosmology.

Read these documents first:

1. [`docs/2026-07/CURRENT_FRAMEWORK_2026-07-23.md`](docs/2026-07/CURRENT_FRAMEWORK_2026-07-23.md)
2. [`docs/2026-07/BRANCH_STATUS_LEDGER_2026-07-23.md`](docs/2026-07/BRANCH_STATUS_LEDGER_2026-07-23.md)
3. [`docs/2026-07/EXPERIMENTAL_PROGRAM_2026-07-23.md`](docs/2026-07/EXPERIMENTAL_PROGRAM_2026-07-23.md)
4. [`LEGACY_MODEL_STATUS.md`](LEGACY_MODEL_STATUS.md)
5. [`UPDATE_MANIFEST_2026-07-23.md`](UPDATE_MANIFEST_2026-07-23.md)

## Active research question

Can environmental record distinguishability generate a finite-memory registration load that modulates a bounded stochastic actualisation channel in a way that cannot be reproduced by standard open-system physics?

The frozen candidate chain is:

```text
unitary system–environment dynamics
        ↓
trace-distance record production
        ↓
finite temporal/spatial registration load
        ↓
bounded stochastic trajectory channel
        ↓
controlled test against the full nuisance space
```

The reference implementation is in:

- [`src/mip_record_conditioned_actualisation.py`](src/mip_record_conditioned_actualisation.py)
- [`tests/test_registration_model.py`](tests/test_registration_model.py)

Run its checks with:

```bash
python -m unittest discover -s tests -v
```

## Meaning of MIP

MIP means **Minimum Interface Points**. In the current framework this is a hypothetical interface interpretation, not a derived Planck-scale lattice and not a universal “frame rate of reality.”

## Current scientific boundary

The active model is **not experimentally confirmed**. Numerical trajectory audits show that monitored open-system models can lock branches while preserving Born-compatible ensemble statistics, and that record geometry can alter locking kinetics. They do not identify one unravelling as objectively real.

A standard structured bath can mimic important nonlinear signatures. A later tangent-space result found local identifiability in selected operating regions. Both results are part of the current status.

## Archived branches

This repository still contains the JWST Metric Archive paper package, Closure D galaxy phenomenology, cosmology scripts, variable-constant toy models, compact-object and jet simulations, visual maps, and ECEM/ECEE outputs.

They are retained for provenance and possible separate research use. They are not evidence for the active MIP bridge unless a newer status document explicitly promotes them.

### JWST Metric Archive package

The original manuscript package remains available, including its LaTeX/PDF paper, figures, scans, and scripts. It is now classified as historical phenomenology rather than a validated quantum-to-gravity derivation.

### ECEM/ECEE planetary battery

The planetary-battery work remains a separate heuristic exoplanet-screening project. Its scores are not habitability probabilities, measured planetary properties, or tests of MIP physics.

## Previous development record

The June update remains indexed in [`UPDATE_MANIFEST_2026-06-30.md`](UPDATE_MANIFEST_2026-06-30.md). Statements in that update that conflict with the 23 July framework are superseded by the July status ledger.

## Interpretation standard

Every result should be labelled as one of:

- standard-physics validation;
- numerical implementation result;
- local identifiability result;
- standard-model degeneracy;
- speculative physical postulate;
- quarantined phenomenology;
- retired branch;
- separate heuristic project.

No simulation output should be described as confirmation of objective collapse, dark matter, dark energy, or a theory of everything without an independent experimental discriminator.
