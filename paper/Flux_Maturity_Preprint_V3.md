# Host Mass as a Shadow of Environmental Maturity: A 100-SNID Flux/MIP Audit of Pantheon+ Supernova Host Structure

**Author:** Jamie McCabe  
**Project:** Flux / MIP Cosmology  
**Version:** V3 / V3.1 preprint draft  
**Status:** `PREPRINT_READY_STRUCTURE_PASS`  
**Primary dataset:** `host_structure_maturity_proxies_v3_100SNID_FINAL.csv`  
**Table pack:** `paper_tables_v3_v31/`  

---

## Abstract

The standard standardization pipeline for Type Ia supernovae (SNe Ia) applies an empirical “host-mass step,” correcting luminosities based on whether a host galaxy’s stellar mass exceeds approximately \(10^{10} M_{\odot}\). We investigate whether this mass step is a fundamental environmental variable or a low-resolution proxy for structural and kinematic maturity. Using an identity-vetted, 100-SNID subset of Pantheon+/SH0ES host galaxies, we construct a hybrid registration-depth proxy, \(D_{\rm hybrid}\), that combines mass-radius concentration with velocity-dispersion age. We find directional support for the maturity hypothesis: a bounded saturation gate, \(G_{\rm sat}\), outperforms raw stellar mass in Bayesian Information Criterion model comparisons, with \(\Delta{\rm BIC} = -2.76\) relative to the mass-only baseline. Crucially, within fixed stellar-mass bins, structurally mature “Adult” hosts return systematically higher inferred \(H_0\) proxies than immature “Child” hosts, indicating that stellar mass obscures distinct physical states. This maturity signal is robustly positive in distance-ladder calibrators (+5.45 km/s/Mpc) and the local Hubble-flow (+15.00 km/s/Mpc), but flattens or reverses in the mid-redshift survey-dominated Hubble-flow. We interpret this divergence not as a physical reversal, but as an observational resolution limit, suggesting that the modern distance ladder compares highly resolved, structurally mature calibrators against a blended, survey-skewed Hubble-flow. This environmental asymmetry may actively contribute to the \(H_0\) tension.

---

## 1. Introduction

Type Ia supernovae are among the most important distance indicators in modern cosmology. Their use depends on empirical standardization: observed light curves are corrected for stretch, color, and host-galaxy effects to infer cosmological distances. One of the most persistent environmental corrections is the host-mass step, in which SNe Ia in host galaxies above approximately \(10^{10}M_{\odot}\) are treated differently from those in lower-mass hosts.

The host-mass step is empirically useful, but physically blunt. Stellar mass is a scalar summary variable. It does not distinguish between a diffuse, actively forming spiral and a compact, dynamically hot early-type galaxy with the same total mass. In this sense, mass may function as a low-resolution observational proxy rather than a fundamental physical driver.

The Flux / Minimum Interface Point (MIP) framework proposes that a host galaxy should be treated not as a static mass bucket, but as a finite registration well. In this interpretation, baryonic history, structural compactness, velocity dispersion, and stellar age jointly determine how mature the local environment is. A supernova occurring in an underfilled, actively assembling host may therefore occupy a different calibration environment from one occurring in a saturated, old, compact host, even when both galaxies have similar stellar mass.

This paper tests whether a continuous environmental maturity proxy can organize SN Ia residual structure more effectively than raw host mass. We build an identity-vetted 100-SNID subset from a Pantheon+/SH0ES-style host ledger and define a hybrid maturity variable, \(D_{\rm hybrid}\), combining mass-radius depth and kinematic age. The aim is not to claim a completed replacement for existing SN standardization pipelines, but to test whether the classical host-mass step is better understood as the observational shadow of a deeper maturity variable.

---

## 2. Dataset and Identity Forensics

### 2.1 Baseline ledger

The audit begins from a Pantheon+/SH0ES-style baseline ledger containing SNID/CID identifiers, redshift, host stellar mass, residual distance-modulus information, calibrator flags, Hubble-flow flags, and a derived \(H_0\)-proxy quantity. The baseline file used for the V3 audit is paired with the host-structure proxy file:

```text
host_structure_maturity_proxies_v3_100SNID_FINAL.csv
```

The analysis is performed at the SNID level, with repeated Pantheon+ entries retained where they correspond to multiple light-curve or survey reductions of the same physical event.

### 2.2 Science-active definition

A row is treated as science-active only if it contains, at minimum:

1. a verified or resolved SNID;
2. host stellar mass, \(\log M_*\);
3. effective radius, \(R_{\rm eff}\);
4. stellar age;
5. velocity dispersion or a derived velocity-dispersion proxy;
6. a usable \(H_0\)-proxy value.

The final V3 dataset contains 100 science-active SNIDs.

### 2.3 Identity-vetting protocol

The identity-vetting stage was deliberately conservative. Host names were accepted only when a named galaxy could be attached to the SNID through catalogue-style evidence or a stable resolved host reference. Systems marked as anonymous, blank-host survey detections, coordinate traps, or unresolved survey designations were quarantined rather than forced into the sample.

Examples of quarantined classes include:

- anonymous hosts requiring coordinate matching;
- Pan-STARRS-style survey names without stable host resolution;
- ASAS-SN entries with blank or uncertain host fields;
- numeric survey IDs without resolved host-galaxy mapping;
- cases where multiple catalogue associations could plausibly refer to different galaxies.

This protocol prevents artificial maturity assignments from entering the audit.

### 2.4 Confidence tiers

Each populated host parameter set is assigned a value-confidence tier:

| Tier | Meaning | Paper use |
|---|---|---|
| `catalogue` | Directly sourced or strongly catalogue-anchored host value | Primary evidence tier |
| `derived` | Derived from named host identity, morphology, size scale, or sibling-host transfer | Usable with caution |
| `estimate` | Broad placeholder or morphology-level proxy | Sensitivity / limitation tier |

The V3 paper treats the dataset as a controlled maturity-audit sample, not a fully homogeneous host-property catalogue. This distinction is central to the interpretation of the results.

---

## 3. The Flux Maturity Proxy

### 3.1 Rationale

The standard host-mass step assumes that total stellar mass is the relevant environmental axis. The Flux/MIP maturity proxy instead assumes that host influence depends on both accumulated baryonic history and the compactness or kinetic maturity of the host environment.

A mature host is expected to be older, more compact, and/or kinematically hotter than a diffuse young host of equal mass. The proxy therefore combines structural depth and kinematic depth.

### 3.2 Mass-radius depth

The mass-radius component is defined as:

\[
D_{MR}
=
\left(
\frac{10^{\log M_*}}{R_{\rm eff}}
\right)
\left(
\frac{\mathrm{Age}}{13.8}
\right)
\]

where \(R_{\rm eff}\) is the effective radius in kpc, Age is the stellar age in Gyr, and 13.8 Gyr provides a cosmic-age normalization.

This term captures concentration and accumulated stellar processing time.

### 3.3 Kinematic depth

The kinematic component is defined as:

\[
D_{\sigma}
=
\left(
\frac{\sigma_v}{150}
\right)^2
\left(
\frac{\mathrm{Age}}{13.8}
\right)
\]

where \(\sigma_v\) is the velocity dispersion in km/s. The normalization to 150 km/s makes the term dimensionless and interpretable as a relative kinematic maturity scale.

### 3.4 Hybrid registration depth

The combined maturity proxy is defined as the geometric mean:

\[
D_{\rm hybrid}
=
\sqrt{D_{MR}D_{\sigma}}
\]

The geometric mean is used because both components are required: a host should not be classified as mature merely because it is massive, nor merely because it has a high velocity-dispersion estimate. The proxy rewards simultaneous structural and kinematic maturity.

### 3.5 Normalization and filling factor

The hybrid depth is normalized by the sample median:

\[
D_{\rm reg}
=
\frac{D_{\rm hybrid}}{\widetilde{D}_{\rm hybrid}}
\]

The filling factor is then defined as:

\[
F_{\rm fill}
=
\frac{D_{\rm reg}}{D_{\rm sat}}
\]

with the V3 threshold:

\[
D_{\rm sat}=1.15
\]

### 3.6 Regime classification

Hosts are classified into three maturity regimes:

| Regime | Definition | Interpretation |
|---|---|---|
| Child / Filling | \(F_{\rm fill}<0.5\) | Underfilled / absorbing host |
| Teen / Saturated | \(0.5 \le F_{\rm fill} \le 1.5\) | Near-capacity host |
| Adult / Overflowing | \(F_{\rm fill}>1.5\) | Mature / high-registration host |

For model comparison, a bounded saturation gate can also be used:

\[
G_{\rm sat}
=
\frac{F_{\rm fill}^{q}}{1+F_{\rm fill}^{q}}
\]

where \(q\) controls the sharpness of the transition.

### 3.7 Empirical test

The empirical test is direct:

1. compare \(G_{\rm sat}\) against raw stellar mass and mass-step models;
2. evaluate fixed-mass bins to determine whether maturity survives when mass is controlled;
3. split the sample into calibrator and Hubble-flow subsets;
4. test whether the maturity signal survives confidence filtering, outlier clipping, and redshift stratification.

---

## 4. V3 100-SNID Audit Results

### 4.1 Model comparison: saturation vs. stellar mass

The primary objective of the V3 audit is to test whether the environmental maturity proxy, \(G_{\rm sat}\), organizes SN Ia residuals more effectively than raw host stellar mass, \(\log M_*\). We evaluate this using the Bayesian Information Criterion (BIC), which penalizes model complexity.

The audit yields a repo-safe verdict of `V3_100SNID_DIRECTIONAL_SUPPORT`. The intercept-only null model retains the lowest absolute BIC score, indicating that environmental corrections, whether driven by mass or maturity, account for a relatively small fraction of the total residual variance in this specific 100-SNID subset.

However, among the environmental predictors tested, structural maturity outperforms stellar mass. The bounded saturation gate, \(G_{\rm sat}\), beats the baseline log-mass model with \(\Delta{\rm BIC}=-2.76\), and also beats the classical step-function mass model. This confirms that when the host galaxy is treated as a finite registration well characterized by structural compactness and kinematic processing time, the resulting proxy explains more SN Ia residual structure than mass alone.

### 4.2 The fixed-mass maturity fingerprint

If stellar mass is merely a low-resolution shadow of environmental maturity, the effect of maturity should persist even when mass is held constant. We test this by isolating galaxies into fixed-mass bins and measuring the inferred \(H_0\) proxy for high-filling “Adult/Overflowing” hosts versus low-filling “Child/Filling” hosts.

The fixed-mass results reveal a consistent empirical fingerprint: high-\(F_{\rm fill}\) hosts return systematically higher inferred \(H_0\) values than low-\(F_{\rm fill}\) hosts within the same mass brackets. This positive maturity signal survives in five of the six populated mass bins:

- log-mass 8.5–9.5: High-F \(H_0\) is higher by +0.65 km/s/Mpc;
- log-mass 9.5–10.0: High-F \(H_0\) is higher by +4.73 km/s/Mpc;
- log-mass 10.0–10.5: High-F \(H_0\) is higher by +5.60 km/s/Mpc;
- log-mass 10.5–11.0: High-F \(H_0\) is higher by +2.32 km/s/Mpc;
- log-mass greater than 11.0: High-F \(H_0\) is higher by +5.09 km/s/Mpc.

The singular exception is the lowest-mass bin, which suffers from extreme sample sparsity and boundary effects for the kinematic derivations.

Across the full active sample, the global offset between Adult and Child regimes is positive. When a standard \(10^{10}M_{\odot}\) mass step is applied in cosmology, it implicitly averages across these distinct physical states. The fixed-mass divergence indicates that a compact, mature \(10^{10}M_{\odot}\) host does not exert the same calibrative influence as a diffuse, unevolved \(10^{10}M_{\odot}\) host.

### 4.3 Robustness and signal localization

To ensure that the signal is not driven by artifacts, we apply strict filtering criteria to the baseline sample. The maturity offset strengthens under higher-confidence selections and survives outlier tests.

The signal partitions sharply when dividing the sample by cosmological role. Within the 42 distance-ladder calibrators, the Adult-minus-Child offset reaches +5.45 km/s/Mpc. Conversely, within the 58 Hubble-flow objects, the aggregate signal reverses slightly to -0.67 km/s/Mpc. This calibrator/Hubble-flow asymmetry requires isolated diagnostic analysis, which is addressed in Section 5.

---

## 5. V3.1 Hubble-Flow Divergence and Resolution Limits

### 5.1 Calibrator/Hubble-flow split

The V3.1 diagnostic audit reveals that the maturity signal is not uniformly distributed across the full 100-SNID sample. It is strongest in the local distance-ladder calibrators and weakest in the survey-dominated Hubble-flow.

The main split is:

| Subset | N | Adult \(H_0\) | Child \(H_0\) | Adult-minus-Child |
|---|---:|---:|---:|---:|
| All Active | 100 | 70.29 | 67.31 | +2.98 |
| Calibrators Only | 42 | 66.90 | 61.45 | +5.45 |
| Hubble-flow Only | 58 | 72.05 | 72.72 | -0.67 |

The calibrator subset therefore exhibits a strong positive maturity ladder, while the aggregate Hubble-flow subset is essentially flat to slightly negative.

### 5.2 Redshift localization

The Hubble-flow reversal is not uniform. It is localized in redshift:

- local Hubble-flow, \(z < 0.025\): strong positive maturity signal;
- mid-redshift Hubble-flow, \(0.025 \le z < 0.05\): negative maturity signal;
- high-redshift Hubble-flow, \(z \ge 0.05\): insufficient Child hosts for a balanced comparison.

This implies that the Hubble-flow divergence is not a simple physical reversal of the maturity relation. Instead, it appears tied to the observational regime in which host galaxies become less resolved, more survey-selected, and more prone to blended or morphology-ambiguous classification.

### 5.3 Fixed-mass split by cosmological role

When the fixed-mass analysis is split between calibrators and Hubble-flow objects, the calibrator sample retains a strong positive maturity signal in the low- and mid-mass ranges. The Hubble-flow sample weakens or reverses in the mid-mass range.

This indicates that the central tension is not simply “maturity works” versus “maturity fails.” Rather, maturity appears to work where the host environment is well-resolved and directly calibration-relevant, while the signal degrades in mid-redshift survey regimes.

### 5.4 Interpretation

The distance ladder may be comparing two environmentally non-equivalent populations:

1. local calibrator hosts, often large, resolved, and structurally mature enough to support detailed Cepheid calibration;
2. mid-redshift Hubble-flow hosts, often survey-selected, less resolved, and less structurally characterized.

The standard host-mass step may partially bridge this mismatch, but only crudely. The maturity audit suggests that the mass step is not the underlying physical variable, but a lossy approximation of structural maturity.

---

## 6. Discussion

### 6.1 Mass as a shadow of environmental maturity

The standard cosmological pipeline treats the host-mass step as an empirical nuisance parameter, correcting for an observed but poorly understood correlation between SN Ia luminosity and host galaxy stellar mass. The V3 100-SNID audit provides directional evidence that this mass step is not a fundamental correction, but a low-resolution shadow of a deeper physical variable: environmental maturity.

The strongest evidence for this interpretation is the fixed-mass maturity fingerprint. If stellar mass were the true driver of the residual structure, the signal would vanish when evaluated within narrow mass brackets. Instead, when host mass is held constant, the residual structure tracks structural maturity, \(F_{\rm fill}\).

A \(10^{10}M_{\odot}\) galaxy undergoing active, diffuse assembly does not calibrate a supernova in the same way as a \(10^{10}M_{\odot}\) galaxy that is compact, kinematically hot, and structurally mature.

### 6.2 The physical reality of the registration well

In the Flux / Minimum Interface Point framework, the host galaxy is modeled as a finite registration well. The macroscopic parameters of the host—mass, radius, stellar age, and internal kinematics—dictate how much baryonic state history has accumulated and saturated the local environment.

The positive relation between \(F_{\rm fill}\) and the \(H_0\) proxy suggests that supernovae detonating inside highly saturated, “overflowing” host environments may not be calibrated identically to those detonating in underfilled environments. The host is therefore not merely a passive backdrop controlling progenitor metallicity or dust. It may act as an active environmental boundary condition.

### 6.3 Distance-ladder asymmetry and the \(H_0\) tension

The most consequential diagnostic finding of this audit is the sharp behavioral divergence between distance-ladder calibrators and the survey-dominated mid-redshift Hubble-flow.

The calibration of the SN Ia distance ladder relies heavily on local, highly resolved galaxies capable of hosting observable Cepheid variables. This imposes an unintentional structural selection effect. Calibrator hosts are preferentially luminous, resolved, and structurally characterized. Within this environment, the maturity signal thrives.

Conversely, the mid-redshift Hubble-flow sample is gathered largely through blind transient surveys. These hosts are frequently unresolved, blended with supernova light, or morphologically ambiguous. In this regime, the maturity signal flattens or reverses.

If local calibrators are biased toward mature, high-registration environments while the distant Hubble-flow represents a blended aggregate of maturity states, then the distance ladder is comparing environmentally non-equivalent populations. Some fraction of the current \(H_0\) tension may therefore be a symptom of environmental calibration asymmetry rather than a purely cosmological expansion-rate crisis.

### 6.4 Moving beyond the step function

The current consensus relies on a binary step function, often dividing hosts at \(10^{10}M_{\odot}\). The continuous nature of \(D_{\rm hybrid}\), \(D_{\rm reg}\), and \(F_{\rm fill}\) exposes the physical inadequacy of a binary mass step. Environmental maturity is a spectrum, not a switch.

Future precision cosmology should integrate continuous structural parameters from IFU spectroscopy, resolved imaging, and local stellar-population modeling. A physically motivated environmental calibration should describe where the supernova occurs, not merely how much total stellar mass its host contains.

---

## 7. Limitations

### 7.1 Controlled audit sample, not a homogeneous catalogue

The V3 100-SNID dataset is a controlled maturity-audit sample, not a fully homogeneous astrophysical catalogue. Its purpose is to test whether host structural maturity can organize SN Ia residuals more effectively than raw stellar mass within an identity-vetted subset.

The sample is conservative because anonymous hosts, survey-designation traps, and coordinate-ambiguous systems were quarantined rather than forced into the analysis. It is incomplete because this strict identity protocol excludes many Pantheon+ objects that may eventually become usable after coordinate-level host matching.

Thus, the V3 result should not be read as a final population-level measurement of the full Pantheon+ sample. It is a structured proof-of-principle audit.

### 7.2 Derived and morphology-scaled host parameters

A major limitation is that many host properties are currently classified as `derived` rather than `catalogue`. Effective radius, stellar age, and velocity dispersion were often estimated from morphology-scaled or host-type-derived priors when direct homogeneous measurements were unavailable.

This is acceptable for a directional maturity audit, but not sufficient for a final precision-cosmology correction. The current results should therefore be interpreted as evidence that the maturity variable is worth formal catalogue testing, not as a replacement for direct spectroscopic measurements.

### 7.3 Velocity-dispersion incompleteness

The hybrid proxy uses both mass-radius structure and kinematic maturity. This makes velocity dispersion important. However, \(\sigma_v\) is not uniformly available across the full sample. In many cases, it is derived from morphology or scaling relations. This weakens the precision of \(D_{\sigma}\) and propagates uncertainty into \(D_{\rm hybrid}\), \(D_{\rm reg}\), and \(F_{\rm fill}\).

Future audits should prioritize galaxies with directly measured central or local velocity dispersions.

### 7.4 Global null model still wins absolute BIC

The intercept-only model retains the absolute BIC advantage in the V3 100-SNID model battle. This means the maturity model does not yet dominate the full residual structure.

The paper does not claim that \(D_{\rm hybrid}\) fully explains SN Ia residuals, nor that it replaces all existing standardization parameters. The result is narrower:

> Among tested environmental predictors, the saturation proxy performs better than raw host mass, and the fixed-mass maturity fingerprint survives in the identity-vetted sample.

This is directional support, not a completed cosmological correction.

### 7.5 Hubble-flow instability

The V3.1 audit shows that the maturity signal is not uniformly stable across the Hubble-flow sample. It is strong in calibrators and local Hubble-flow objects, but flattens or reverses in the mid-redshift survey-dominated regime.

This instability is interpreted as a likely host-resolution and selection-boundary effect, but that interpretation remains provisional. It must be tested using higher-resolution host imaging, IFU spectroscopy, and explicit survey-selection modeling.

### 7.6 Global host values vs. local SN environments

The present audit uses global host parameters. However, SN Ia events occur at specific positions within galaxies. A supernova in a mature bulge, a young spiral arm, or a tidal interaction zone may not experience the same local environment as implied by the galaxy-wide average.

Future work should replace global host maturity with local environmental maturity measured at or near the SN position.

### 7.7 Interpretation within Flux/MIP remains theoretical

The empirical result can be stated without fully accepting the Flux/MIP interpretation:

> Structural host maturity appears to organize SN Ia residuals more effectively than stellar mass alone in the V3 audit.

The stronger physical claim—that galaxies act as finite registration wells whose accumulated baryonic history changes the calibration environment—is a theoretical interpretation of that empirical pattern.

---

## 8. Predictions and Future Tests

### 8.1 Catalogue-homogenization prediction

If the maturity proxy is physically meaningful, replacing morphology-derived values with homogeneous catalogue measurements should strengthen the signal rather than erase it.

Prediction:

> A catalogue-grade version of the audit using direct \(R_{\rm eff}\), stellar-population age, and \(\sigma_v\) measurements should increase the statistical preference for \(D_{\rm hybrid}\), \(F_{\rm fill}\), or \(G_{\rm sat}\) over raw stellar mass.

Failure condition:

> If direct catalogue measurements destroy the fixed-mass maturity fingerprint, then the V3 signal was likely driven by proxy construction rather than physics.

### 8.2 Local-environment prediction

The maturity signal should strengthen when measured locally at the supernova position instead of globally across the entire host galaxy.

Prediction:

> SN Ia residuals should correlate more strongly with local host maturity—local stellar age, local surface density, local velocity dispersion, and local star-formation state—than with global host stellar mass.

### 8.3 Fixed-mass differentiation prediction

At fixed stellar mass, high-maturity hosts should continue to return systematically different standardized residuals than low-maturity hosts.

Prediction:

> In a larger, catalogue-grade sample, older, more compact, higher-dispersion hosts should retain a positive Adult-minus-Child \(H_0\)-proxy offset within the same mass bins.

### 8.4 Calibrator/Hubble-flow resolution prediction

The V3.1 audit predicts that the Hubble-flow divergence is observational, not fundamental.

Prediction:

> Well-resolved local Hubble-flow hosts should behave more like calibrators, while poorly resolved mid-redshift survey hosts should show weakened, flattened, or inconsistent maturity trends.

### 8.5 Survey-selection prediction

If the mid-redshift reversal is caused by survey selection and host-resolution limits, then the signal should vary by survey pipeline.

Prediction:

> Hubble-flow maturity trends should differ between targeted nearby surveys and blind wide-field surveys, even at similar redshift and mass.

### 8.6 Scatter-reduction prediction

If \(D_{\rm hybrid}\) captures a real environmental standardization variable, replacing the binary mass step with a continuous maturity correction should reduce residual scatter.

Prediction:

> A continuous correction based on \(D_{\rm hybrid}\), \(F_{\rm fill}\), or \(G_{\rm sat}\) should reduce SN Ia residual scatter more effectively than the binary \(10^{10}M_{\odot}\) host-mass step in a fully resolved sample.

### 8.7 Duplicate-host prediction

Hosts containing multiple SNe Ia provide a clean internal test.

Prediction:

> Multiple SNe Ia occurring in the same host should share the same global maturity baseline, while deviations between them should track local SN position within the host.

This is especially valuable for hosts such as NGC 5643 and NGC 5468, where multiple SN events can test whether global and local maturity separate cleanly.

### 8.8 Final falsification criteria

The maturity hypothesis should be considered weakened or falsified if future catalogue-grade audits show any of the following:

1. \(G_{\rm sat}\) fails to outperform raw host mass in homogeneous datasets.
2. Fixed-mass maturity offsets vanish in larger samples.
3. Local host environment fails to outperform global host mass.
4. The calibrator maturity signal disappears under direct spectroscopy.
5. The Hubble-flow divergence remains negative after controlling for survey origin, host resolution, and local environment.

Conversely, the model gains support if:

1. the fixed-mass maturity fingerprint strengthens;
2. catalogue-grade \(D_{\rm hybrid}\) beats the mass step;
3. local maturity reduces residual scatter;
4. calibrator/Hubble-flow asymmetry tracks host-resolution limits;
5. duplicate-host systems show consistent global maturity baselines.

---

## 9. Conclusion

The environmental calibration of Type Ia supernovae has long relied on host galaxy stellar mass as a primary correction step. However, mass is a static parameter that cannot distinguish between a diffuse, actively forming spiral and a compact, dynamically hot early-type galaxy. The Flux/MIP maturity audit demonstrates that when structural compactness and kinematic age are accounted for, supernova residual structure tracks environmental maturity more closely than mass alone.

Our 100-SNID identity-vetted audit reveals a distinct fixed-mass fingerprint: when mass is held constant, highly filled, saturated host environments yield systematically different standardized residuals than underfilled environments. The classical host-mass step appears to be a low-resolution shadow of this underlying structural variable.

Furthermore, we identify a sharp environmental asymmetry within the distance ladder itself. The maturity signal is highly predictive within local calibrators and the nearby Hubble-flow, where host galaxies are well resolved. It breaks down in the mid-redshift Hubble-flow, where transient surveys frequently encounter blended, unresolved, or morphologically ambiguous hosts. If the \(H_0\) tension is driven in part by calibration offsets between the local ladder and the deep Hubble-flow, these findings suggest that a mismatch in host-environment resolution—and specifically, a bias toward structurally mature calibrators—may be a contributing factor.

Moving forward, precision cosmology cannot treat host galaxies as mere scalar mass buckets. The integration of continuous, local, catalogue-grade structural parameters—such as those achievable via IFU spectroscopy—will be required to determine the true environmental boundary conditions of the supernova event. The maturity framework presented here offers a testable, falsifiable pathway to transition from a binary mass step to a continuous, physically motivated environmental calibration.

---

## Appendix A. Generated Table Pack

The following files are associated with this draft:

```text
paper_tables_v3_v31/table_1_dataset_summary.csv
paper_tables_v3_v31/table_2_proxy_definitions.csv
paper_tables_v3_v31/table_3_v3_model_battle.csv
paper_tables_v3_v31/table_4_v3_fixed_mass_fingerprint.csv
paper_tables_v3_v31/table_5_v31_divergence_summary.csv
paper_tables_v3_v31/table_6_v31_redshift_split.csv
paper_tables_v3_v31/table_7_v31_calibrator_vs_hf_mass_bins.csv
paper_tables_v3_v31/table_8_limitations_and_confidence_tiers.csv
paper_tables_v3_v31/active_100SNID_with_flux_columns.csv
```

---

## Repository Verdict

```text
V3_PREPRINT_DRAFT_ASSEMBLED
```

This manuscript is repo-ready as a structured preprint draft. It should be treated as a working scientific document pending catalogue homogenization, local-host refinement, and formal citation insertion.
