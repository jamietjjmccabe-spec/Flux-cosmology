# Curvature-Sourced Lindblad Dephasing from a Baryonic Registration Field

## Flux-MIP V2.1 and a Four-State Atom-Interferometer Null Test

**Author:** Jamie Thomas McCabe  
**Status:** Final prospective paper / frozen phenomenological benchmark  
**Date:** 12 July 2026  
**Benchmark signature:** `f819067c8feca48e6d0a51b790e01f5a1be599f4afa22c6accc46f3967827670`  
**Optimizer signature:** `ed3bde6ab4adcef685d3cee33942b0cc5467118dd82f86ad4d2863e2603ef6b0`

> Flux-MIP V2.1 is a falsifiable phenomenological extension of open quantum dynamics. It is not empirical evidence, a completed cosmological theory, or a derivation from quantum gravity.

## Abstract

We formulate and preregister a phenomenological modification of open quantum dynamics in which a baryon-sourced registration field contributes a Hermitian phase term and a Lindblad dephasing channel. A regularized inverse-square precursor is spatially high-pass filtered through its Laplacian, producing a bounded response that rejects constant and locally affine backgrounds while remaining sensitive to structured matter. For two interferometer arms, the theory predicts an odd-parity phase shift proportional to `I_phi = integral Delta F dt` and an even-parity visibility loss `V/V_standard = exp(-gamma_MIP X)`, where `X = 1/2 integral (Delta F)^2 dt`. The same operator fixes momentum-diffusion and heating predictions. A conservative dedicated geometry uses a 7.4112 kg, 8 mm-thick tungsten plate between arms separated by 40 mm, with 27 mm and 5 mm clearances. The frozen model gives `X = 0.085181848 s`, `|I_phi| = 0.412750911 s`, and `J_grad = 1.2141e4 s m^-2`. A four-state mirror-modulation protocol separates even dephasing from odd phase shifts and common-mode technical loss. A one-sided 1% residual-contrast sensitivity would imply `gamma_MIP < 0.118 s^-1`. The framework is mathematically admissible and falsifiable but is not empirically supported or relativistically complete.

## 1. Claim boundary

The additional channel is

\[
\frac{d\rho}{d\tau}=-\frac{i}{\hbar}[H_0+\hbar\beta_{\rm MIP}F_{\rm reg}(\hat x),\rho]
-\frac{\gamma_{\rm MIP}}{2}[F_{\rm reg}(\hat x),[F_{\rm reg}(\hat x),\rho]].
\]

The additional-channel null is

\[
\gamma_{\rm MIP}=0,\qquad \beta_{\rm MIP}=0.
\]

Ordinary environmental decoherence, gravity-gradient wave-packet mismatch, pulse inefficiency, atom loss, and detection effects remain in the standard model and nuisance structure.

The paper does not claim that the kernel is derived from a fundamental action, that the model is covariant, that the laboratory channel proves metric-memory cosmology, or that processed Overstreet figure data constrain `gamma_MIP`.

## 2. Frozen field

\[
G_\ell(r)=\frac{\ell^2}{r^2+\ell^2},
\]

\[
D_{\rm reg}(\mathbf x)=\frac{\kappa_{\rm reg}}{\rho_*V_*}
\int \rho_b(\mathbf r')G_\ell(|\mathbf x-\mathbf r'|)d^3r',
\qquad
V_*=\frac{4\pi}{3}\ell^3.
\]

The local susceptibility is

\[
S_{\rm reg}=-\ell^2\nabla^2D_{\rm reg},
\qquad
F_{\rm reg}=\tanh(S_{\rm reg}/S_*).
\]

For the point kernel,

\[
-\ell^2\nabla^2G_\ell(r)=
\frac{2\ell^4(3\ell^2-r^2)}{(r^2+\ell^2)^3}.
\]

It changes sign at `r = sqrt(3) ell` and decays as `r^-4`. Constant and uniform-gradient precursor backgrounds are rejected by construction. Nearby apparatus is not ignored and must be modeled.

### Frozen values

| Parameter | Value |
|---|---:|
| `ell` | `0.010 m` |
| kernel power | `2` |
| `rho_*` | `1000 kg m^-3` |
| `V_*` | `4.188790204786391e-6 m^3` |
| `S_*` | `10` |
| `kappa_reg` | `1` |
| ordinary-matter factor | `1` |
| adaptive threshold | prohibited |
| temporal high-pass | prohibited |

## 3. Observables

For paths `x1(t)` and `x2(t)`,

\[
\Delta F(t)=F(x_1(t))-F(x_2(t)),
\]

\[
I_\phi=\int \Delta F(t)dt,
\qquad
X=\frac12\int[\Delta F(t)]^2dt.
\]

The predicted phase and visibility are

\[
\delta\phi_{\rm MIP}=\beta_{\rm MIP}I_\phi,
\qquad
\frac{V}{V_{\rm standard}}=e^{-\gamma_{\rm MIP}X}.
\]

The linked momentum diffusion is

\[
J_{ij}=\frac12\int(\partial_iF_1\partial_jF_1+\partial_iF_2\partial_jF_2)dt,
\]

\[
\Delta\langle p_ip_j\rangle=\gamma_{\rm MIP}\hbar^2J_{ij},
\qquad
\Delta E_k=\frac{\gamma_{\rm MIP}\hbar^2}{2m}\operatorname{Tr}J.
\]

## 4. Pre-data failure checks

A 1 cm Gaussian V1 field was too short-ranged to sample the Overstreet geometry. A direct inverse-square V2 gate was vulnerable to absolute-background saturation. V2.1 fixed the curvature high-pass rule before prospective contrast analysis.

The nominal Overstreet V2.1 exposure was `X = 2.33103e-5 s`. The available uploaded `Fig2.csv` and `Fig2.tab.tsv` files contained processed phase points rather than shot-level populations or visibility estimates. No retrospective bound was invented.

## 5. Dedicated conservative geometry

| Quantity | Value |
|---|---:|
| arm separation | `40 mm` |
| interaction duration | `1.0 s` |
| modeled path length | `200 mm` |
| tungsten dimensions | `160 x 300 x 8 mm` |
| plate mass | `7.4112 kg` |
| clearances | `27 mm` and `5 mm` |
| `X` | `0.085181848 s` |
| `|I_phi|` | `0.412750911 s` |
| `J_grad` | `12140.990296 s m^-2` |
| `Delta E` for Rb-87 and `gamma=1 s^-1` | `4.678016e-40 J` |

The mathematical maximum used a 27 mm-thick plate with 12 mm and 1 mm clearances and gave `X = 0.699776 s`; it is not the primary engineering design.

## 6. Four-state Phase Ladder

- **L0 - retracted null:** `X <= 1e-4 X_signal`.
- **L1 - symmetric control:** plate centered, 16 mm/16 mm clearances, `Delta F ~= 0`.
- **L2 - asymmetric signal:** 27 mm/5 mm clearances.
- **L3 - mirrored signal:** 5 mm/27 mm clearances.

The frozen parity relations are

\[
X_{L2}=X_{L3},
\qquad
I_{\phi,L2}=-I_{\phi,L3}.
\]

Use balanced Latin ordering, randomized supercycles, identical motion and settling profiles, blinding, and both signs of `k_eff`.

## 7. Estimators

\[
y_2=-\ln(V_{L2}/V_{L1}),
\qquad
y_3=-\ln(V_{L3}/V_{L1}),
\]

\[
\widehat\gamma_{\rm MIP}=\frac{y_2+y_3}{2X_{\rm signal}}.
\]

The odd visibility null is

\[
y_{\rm odd}=\frac12\ln(V_{L3}/V_{L2}),
\]

which must be consistent with zero.

At block level,

\[
-\ln V_{bsk}=a_b+\gamma_{\rm MIP}X_s+\boldsymbol\eta^T\mathbf z_{bsk}+\epsilon_{bsk}.
\]

## 8. Projected sensitivity

\[
\gamma_{\rm MIP}^{95}=\frac{-\ln(1-\epsilon)}{0.085181848\,\mathrm{s}}.
\]

| One-sided residual-loss limit | Projected limit |
|---:|---:|
| `0.1%` | `0.01175 s^-1` |
| `0.5%` | `0.05885 s^-1` |
| `1.0%` | `0.1180 s^-1` |

These are projections, not measurements.

## 9. Required controls

The experiment must independently model Newtonian phase and phase-space closure and monitor geometry, cloud size and centroid, atom number, pulse efficiency, magnetics, electrostatics, source motion, temperature, vibration, optics, and detection gain. The cloud qualification is

\[
3\sigma_y+r_{\rm centroid,error}\le4\,\mathrm{mm}.
\]

All shields and supports must be included in both the registration and conventional forward models.

## 10. Falsification criteria

A candidate anomaly requires a positive even-parity slope versus frozen `X`, an odd visibility null, invariance under `k_eff` reversal, no cloud clipping or overlap deficit, survival of preregistered covariates, and scaling with a second frozen `X` configuration.

A null result is reported as `0 <= gamma_MIP < gamma_95`. The kernel, scale, gate, background rule, and geometry exposure cannot be modified after unblinding.

## 11. Relation to cosmology

The wider Flux programme interprets accumulated registration as possible metric memory. The laboratory benchmark does not derive or validate that cosmological sector. A complete cosmology still requires a covariant action, energy accounting, perturbation equations, and agreement with Solar-System, CMB, lensing, and structure-growth constraints.

## 12. Limitations

The benchmark lacks a fundamental action, relativistic completion, derivation of the Markov approximation, a first-principles relation between `beta_MIP` and `gamma_MIP`, and a demonstrated separation from every possible unobserved environmental channel. The parity signature is discriminating but not unique.

## Conclusion

Flux-MIP V2.1 is a fixed, executable open-system hypothesis with linked phase, visibility, diffusion, and heating predictions. The dedicated four-state experiment provides a direct route to either a source-correlated anomaly or an upper bound on the frozen channel. The present paper is a test specification, not evidence for the effect.

## References

1. G. Lindblad, "On the generators of quantum dynamical semigroups," *Communications in Mathematical Physics* **48**, 119-130 (1976). DOI: 10.1007/BF01608499.
2. V. Gorini, A. Kossakowski, and E. C. G. Sudarshan, "Completely positive dynamical semigroups of N-level systems," *Journal of Mathematical Physics* **17**, 821-825 (1976). DOI: 10.1063/1.522979.
3. C. Overstreet et al., "Observation of a gravitational Aharonov-Bohm effect," *Science* **375**, 226-229 (2022). DOI: 10.1126/science.abl7152.
4. C. D. Panda et al., "Measuring gravitational attraction with a lattice atom interferometer," *Nature* **631**, 515-520 (2024). DOI: 10.1038/s41586-024-07561-3.
5. D. M. Harber et al., "Measurement of the Casimir-Polder force through center-of-mass oscillations of a Bose-Einstein condensate," *Physical Review A* **72**, 033610 (2005). DOI: 10.1103/PhysRevA.72.033610.
