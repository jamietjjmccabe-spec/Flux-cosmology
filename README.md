# Flux/MIP Metric Archive JWST Paper Package

This repository package contains a working manuscript and supporting scan files for:

**A Scale-Gated Metric Archive Model for JWST Early Galaxy Formation**

## Contents

- `flux_metric_archive_jwst.tex` - full LaTeX manuscript source.
- `flux_metric_archive_jwst.pdf` - compiled manuscript PDF.
- `figures/flux_env_filter.png` - environmental gate containment plot.
- `figures/id12385_zplane_stability.png` - ID 12385 redshift-plane stability plot.
- `scans/flux_containment_audit.csv` - CMB/RSD/galactic containment scale values.
- `scans/phase10_2_classification_summary.csv` - final Phase-10.2 classification counts.
- `scans/flux_uncover_final_stability_top15.csv` - top candidates by registration score from the local audit output.
- `scans/id12385_redshift_plane_scan.csv` - ID 12385 magnification and registration values by redshift plane.
- `scripts/flux_containment_audit.py` - reproduces containment table.
- `scripts/flux_phase10_sample_magnification_map.py` - samples FITS magnification maps at source RA/DEC.
- `scripts/flux_phase10_merge_stability.py` - merges redshift-plane sampled catalogues into final stability classifications.

## June 2026 development update

The 28-day development record ending 30 June 2026 is indexed in [`UPDATE_MANIFEST_2026-06-30.md`](UPDATE_MANIFEST_2026-06-30.md). It includes:

- the active Minimum Interface Point terminology and registration-depth bridge;
- the metric-memory dark-matter differential;
- the frozen pre-CLASS integration boundaries;
- the full project progress record for 2-30 June 2026;
- the ECEM/ECEE planetary-battery candidate tables, compact SVG figures, and regeneration script;
- the three-body flux-fragment toy experiment.

The repository's existing `cmbr.py` remains the active version. It contains the later FLRW-aware coherence-field implementation and was not replaced by an older local file with the same name.

## Notes

The numerical candidate results in this package come from local working audits provided during development. Before formal submission, re-run the scripts against the newest official source catalogues and lensing maps, then update the CSVs and manuscript tables.

## Build

From this directory:

```bash
latexmk -pdf flux_metric_archive_jwst.tex
```

## Required Python packages for scripts

```bash
pip install numpy pandas astropy scipy matplotlib
```

## Conservative interpretation

The paper deliberately avoids claiming proof of Flux/MIP cosmology. The central claim is that the scale-gated Metric Archive model provides a selective, falsifiable filter that identifies a small registered minority of compact, de-lensed high-redshift candidates while rejecting the majority.

The June 2026 ECEM/ECEE scores are likewise heuristic screening quantities, not habitability probabilities or measurements of atmospheres, magnetic fields, composition, or internal heat.
