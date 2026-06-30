# 28-day update manifest - 30 June 2026

This update records and adds the work products completed during the 28-day interval ending 30 June 2026 that were not represented in the May project index.

## Theory records

- `docs/2026-06/PROJECT_PROGRESS_2026-06-02_TO_2026-06-30.md`
- `docs/2026-06/MINIMUM_INTERFACE_POINTS_AND_REGISTRATION_DEPTH.md`
- `docs/2026-06/DARK_MATTER_AS_METRIC_MEMORY.md`
- `docs/2026-06/ACTIVE_BASELINE_AND_FREEZE_BOUNDARIES.md`

## Code and experiments

- `src/ecee_planetary_battery.py` - existing full ECEM/ECEE source, byte-verified against the local project copy.
- `src/2.py` - existing cross-domain universal-flux conceptual prototype; retained as heuristic material, not part of the active astrophysical baseline.
- `experiments/three_body_flux_fragments.py` - descriptively named three-body toy experiment.

Other local cosmology helpers were also byte-verified against their repository copies and were not duplicated.

## ECEE outputs

- `outputs/ecee/nasa_ecee_usable_battery_candidates_compact.csv` - 57-candidate compact table.
- `outputs/ecee/nasa_ecee_usable_battery_top50_compact.csv` - top-50 compact table.
- `outputs/ecee/make_candidate_figures.py` - dependency-light SVG generator.
- `outputs/ecee/figures/` - three compact SVG figures.
- `outputs/ecee/README.md` - interpretation and provenance notes.

The generated 6,224-row, approximately 6.8 MB intermediate score table is not duplicated; the source screener and compact research outputs preserve the reproducible result.

## Versioning decision

The repository's existing `cmbr.py` was retained. It contains a later FLRW-aware coherence-field implementation and should not be replaced by the older local plateau-potential script with the same filename.
