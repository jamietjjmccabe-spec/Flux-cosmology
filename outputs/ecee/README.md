# June 2026 ECEE planetary-battery outputs

This directory contains the reproducible, compact output set from the June 2026 NASA Exoplanet Archive screening experiment.

## Files

- `nasa_ecee_usable_battery_candidates_compact.csv` - 57-object compact candidate table with principal observables, scores, classifications, and follow-up notes.
- `nasa_ecee_usable_battery_top50_compact.csv` - top-50 compact ranking table.
- `make_candidate_figures.py` - regenerates the three SVG figures.
- `figures/nasa_ecee_mass_radius_candidates.svg` - candidate mass-radius search space.
- `figures/nasa_ecee_core_to_stellar_candidates.svg` - insolation versus core-to-stellar proxy.
- `figures/nasa_ecee_usable_battery_top25.svg` - top-25 usable-score ranking.

The full local scored table contains 6,224 rows and 75 scoring/data columns. It is intentionally not duplicated here because it is a generated intermediate of roughly 6.8 MB; the compact candidate tables and scoring script preserve the research result and can be regenerated from the source catalogue.

## Interpretation warning

ECEM, ECEE, and planetary-battery values are heuristic screening scores. Many input quantities are inferred or imputed. A high score does not establish habitability, a magnetic field, an atmosphere, a rocky composition, or an internal heat budget. The output is a target-prioritization experiment only.

Candidates with inferred radius or mass, radial-velocity \(M\sin i\), weak stellar-age constraints, or uncertain volatile content require special caution. The `search_note` field records several of these limitations.
