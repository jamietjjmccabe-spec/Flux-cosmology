# Flux-MIP V2.1 Phase Ladder protocol

**Frozen benchmark:** `f819067c8feca48e6d0a51b790e01f5a1be599f4afa22c6accc46f3967827670`  
**Frozen optimizer:** `ed3bde6ab4adcef685d3cee33942b0cc5467118dd82f86ad4d2863e2603ef6b0`

## Conservative geometry

- arm separation: `40 mm`
- interrogation duration: `1.0 s`
- tungsten plate: `160 x 300 x 8 mm`
- plate mass: `7.4112 kg`
- asymmetric clearances: `27 mm / 5 mm`
- `X_signal = 0.085181848 s`
- `|I_phase| = 0.412750911 s`

## Four states

1. `L0`: retracted null, `X <= 1e-4 X_signal`.
2. `L1`: centered common-mode control, `16 mm / 16 mm` clearances.
3. `L2`: asymmetric signal, `27 mm / 5 mm`.
4. `L3`: mirrored signal, `5 mm / 27 mm`.

Required parity:

```math
X_{L2}=X_{L3},\qquad I_{\phi,L2}=-I_{\phi,L3}.
```

## Acquisition

Use balanced Latin ordering, randomized supercycles, identical movement and settling profiles, both signs of `k_eff`, and one frozen contrast estimator.

## Primary estimator

```math
y_2=-\ln(V_{L2}/V_{L1}),\qquad y_3=-\ln(V_{L3}/V_{L1}),
```

```math
\widehat\gamma_{\rm MIP}=\frac{y_2+y_3}{2X_{\rm signal}}.
```

The odd visibility null is

```math
y_{\rm odd}=\frac12\ln(V_{L3}/V_{L2}).
```

## Qualification

```math
3\sigma_y+r_{\rm centroid,error}\le4\,\mathrm{mm}.
```

Monitor atom number, centroid, width, pulse efficiency, magnetics, electrostatics, vibration, source position, temperature, optics, and detector gain.

## Decision rule

A candidate requires an even-parity slope versus frozen `X`, odd visibility null, `k_eff` invariance, no clipping or overlap deficit, survival of preregistered covariates, and confirmation with a second preregistered `X`.
