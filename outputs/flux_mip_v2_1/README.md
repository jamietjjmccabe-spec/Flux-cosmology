# Flux-MIP V2.1 compact outputs

These files preserve the frozen benchmark declaration, validated plate designs, the nominal Overstreet forward result, and the audit of the uploaded Figure 2 tables.

## Interpretation

- `flux_mip_v2_1_background_rule.json` freezes the curvature high-pass rule.
- `v2_1_plate_validated_designs.csv` contains the mathematical maximum and the conservative 5 mm-clearance design.
- `overstreet_2022_rx0_v2_1_results.json` is a forward-model result, not an empirical fit.
- `Fig2_manifest_summary.json` and `Fig2_tab_manifest_summary.json` record that the available Figure 2 tables contain processed phase data and cannot support an independent residual-contrast bound.

Large NPZ arrays and exhaustive candidate tables are intentionally omitted. They can be regenerated from the committed evaluator and optimizer.
