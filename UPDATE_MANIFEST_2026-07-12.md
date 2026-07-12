# July 2026 Flux-MIP V2.1 update manifest

This update records the development from the zero-amplitude/re-registration motivation to a fixed open-system benchmark and a prospective dedicated atom-interferometer test.

## Final paper

- `papers/flux_mip_v2_1/Flux_MIP_V2_1_Final_Paper.tex`
- `papers/flux_mip_v2_1/Flux_MIP_V2_1_Final_Paper.md`
- `papers/flux_mip_v2_1/README.md`

The PDF and DOCX are generated release artifacts. The LaTeX source is the canonical paper source.

## Frozen theory records

- `docs/2026-07/FLUX_MIP_V2_1_THEORY.md`
- `docs/2026-07/PHASE_LADDER_PROTOCOL.md`
- `docs/2026-07/PROJECT_STATUS_2026-07-12.md`

## Executable code

- `src/flux_lindblad_mip_v2_1.py`
- `experiments/v2_1_geometry_optimizer.py`
- `experiments/flux_mip_bound.py`
- `experiments/dataverse_fetch_manifest.py`

## Compact outputs

- `outputs/flux_mip_v2_1/flux_mip_v2_1_background_rule.json`
- `outputs/flux_mip_v2_1/v2_1_plate_validated_designs.csv`
- `outputs/flux_mip_v2_1/v2_1_plate_optimizer_coarse_top.csv`
- `outputs/flux_mip_v2_1/v2_1_plate_optimizer_coarse_summary.json`
- `outputs/flux_mip_v2_1/overstreet_2022_rx0_v2_1_results.json`
- `outputs/flux_mip_v2_1/Fig2_manifest_summary.json`
- `outputs/flux_mip_v2_1/Fig2_tab_manifest_summary.json`

Large intermediate NPZ arrays are intentionally not added. They can be regenerated from the committed scripts.

## Frozen signatures

- Benchmark: `f819067c8feca48e6d0a51b790e01f5a1be599f4afa22c6accc46f3967827670`
- Optimizer: `ed3bde6ab4adcef685d3cee33942b0cc5467118dd82f86ad4d2863e2603ef6b0`
