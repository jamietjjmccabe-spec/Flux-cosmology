# Dark Matter as the Metric Shadow of Dead History

## A Flux/MIP Archive Model for Galactic Rotation Curves

**Author:** Jamie McCabe, with AI-assisted formalization  
**Version:** v3q paper bundle, 2026-05-12  

---

## Abstract

The standard interpretation of galactic rotation curves introduces a non-luminous collisionless mass component, while modified-gravity alternatives usually encode the observed excess acceleration through a universal acceleration scale. This paper develops a different interpretation within the Flux/MIP framework: galactic dark-sector behaviour is treated as the metric response to historical registration structure. In this view, baryonic matter is not only a source of mass-energy but also a persistent archive of prior actualization events. We define a global spatial archive depth, introduce a saturating archival-efficiency function, and extend the model with a compact-remnant archive amplifier motivated by white dwarfs, neutron stars, and stellar/rogue black holes as high historical-depth nodes. The resulting v3q model is tested against the official SPARC175 mass-model catalogue. In the present implementation, Newtonian baryons alone give a galaxy-mean MAE of approximately **39.98 km/s**, the v3o spatial archive model gives **27.96 km/s**, and the v3q compact-remnant archive model gives **18.71 km/s**. The model does not yet solve SPARC at MOND-level precision, but it substantially improves over the baryonic baseline and isolates a structured residual pathway associated with old, gas-poor, remnant-rich systems. The central thesis is: **dark matter is the metric shadow of dead history.**

---

## 1. Introduction

The galactic rotation-curve problem is usually stated as a missing-mass problem. Observed circular velocities remain too high at large radii compared with the Newtonian prediction from visible baryons. The standard cosmological answer is cold dark matter: an extended collisionless halo surrounding galaxies. However, the tight empirical coupling between baryons and acceleration suggests that the dark-sector signal is not arbitrary. It tracks baryonic structure with remarkable fidelity.

The Flux/MIP programme reframes this coupling. Instead of asking where unseen mass is located, we ask how the metric stores, exports, and responds to historical registration. In this framework, matter is interpreted as bound history: a stable residue of measurement-induced actualization events. Mass is therefore not merely substance in space, but persistent causal record.

The present paper develops the v3q version of this idea for galaxy rotation curves. The model starts with a global archive-depth scalar derived from baryonic surface compression, then adds a saturating efficiency law and a compact-remnant wake amplifier. The SPARC175 dataset provides a first falsification benchmark.

---

## 2. Flux/MIP Ontology

The Flux/MIP framework treats physical reality as a system of actualization events. MIP has two associated meanings:

1. **Measurement-Induced Projection:** the transition from unresolved quantum potential into classical registration.
2. **Minimum Intersection Point:** the local condition under which a possibility becomes a persistent event in the metric.

Matter is the stable residue of repeated actualization. The nucleus of an atom is interpreted as a hardened archive of prior registrations, while electron shells form the active interface through which the object exchanges information and entropy with its environment.

The galactic model therefore distinguishes between:

- **spatial compression**, the present baryonic density state;
- **temporal persistence**, the accumulated duration of stable registration;
- **compact-remnant texture**, the hidden population of high-depth archive nodes.

---

## 3. Mass as Bound History

In standard dynamics, inertial mass is a scalar resistance to acceleration. In the Flux/MIP interpretation, inertial mass is the work required to re-register an already persistent historical record. This gives a compact statement:

> Mass is persistent bound history.

Gravity is then the metric response to gradients and concentrations of this bound historical structure.

This interpretation does not require changing local inertial mass directly. The model instead modifies the effective nonlocal wake produced by the galaxy-level archive. This avoids immediately violating equivalence-principle constraints, because the correction is treated as a metric response rather than a composition-dependent inertial-mass change.

---

## 4. Black Holes as Horizon-Locked Historical Archives

Black holes are treated as the limiting case of historical registration. A collapse event occurs when accumulated temporal depth and coherence exceed the local entropy-export capacity:

```text
D_t C_sigma / S_leak > Theta_BH
```

where:

- `D_t` is temporal depth,
- `C_sigma` is coherence/compression,
- `S_leak` is entropy export capacity,
- `Theta_BH` is the horizon-lock threshold.

The event horizon is interpreted as a causal vault: a topological boundary that quarantines saturated historical registration from the surrounding metric. In this language, black holes are not singularities of matter but singularities of accumulated history.

This motivates v3q. Mature galaxies contain compact remnants: white dwarfs, neutron stars, stellar black holes, and potentially rogue black holes. They need not dominate baryonic mass, but they may dominate historical-depth texture.

---

## 5. The v3q Master Equations

### 5.1 Global Spatial Archive Depth

The galaxy is assigned a global baryonic archive depth:

```text
D_spatial = (Sigma_bar / Sigma_ref)^0.08
```

with:

```text
Sigma_bar = M_bar / (2 pi R_disk^2)
M_bar = M_star + M_gas
M_star = Upsilon_star L_3.6
M_gas = 1.33 M_HI
Sigma_ref = 4.08e9 Msun kpc^-2
Upsilon_star = 0.5
```

### 5.2 Saturating Archival Efficiency

The v3o causal-sluice efficiency is:

```text
eta(D) = 8.0 D^-7 / [1 + (D / 0.68)^8]
```

This term prevents runaway archive leakage at high depth.

### 5.3 Memory Radius

The radial expression scale of the wake is:

```text
L_chi = 7.9 R_disk D^-4.6
```

Low-depth systems generate broader wakes; high-depth systems compress their archive more tightly.

### 5.4 Universal Wake Texture

The enclosed historical wake mass is:

```text
M_hist(<r) = eta_eff M_bar x^s / (1+x)^(s-1)
```

with:

```text
x = r / L_chi
s = 1.36
```

### 5.5 Compact-Remnant Archive Proxy

The raw remnant proxy is:

```text
R_rem = (1 - f_gas) ((10 - T) / 10)
f_gas = M_gas / M_bar
```

where `T` is Hubble type. The term targets old, gas-poor systems.

### 5.6 Retention Gate

A physically conservative variant requires a sufficiently deep potential to retain and organize compact archive nodes:

```text
G_ret = (V_max / V_ref)^n / [1 + (V_max / V_ref)^n]
```

A useful physical gate explored here is:

```text
V_ref = 200 km/s
n = 4
```

The retained remnant archive is:

```text
R_eff = R_rem G_ret
```

### 5.7 Effective Archival Efficiency

The v3q remnant amplifier is:

```text
eta_eff = eta(D) [1 + chi_rem R_eff^zeta]
```

with the empirical best ungated values:

```text
chi_rem = 3.0
zeta = 0.5
```

The square-root dependence implies sublinear amplification:

```text
eta_eff = eta(D) [1 + 3 sqrt(R_eff)]
```

Compact remnants therefore act as a coherence amplifier, not as a linear mass reservoir.

### 5.8 Total Rotation Curve

The historical acceleration is:

```text
g_hist(r) = G M_hist(<r) / r^2
```

The baryonic velocity contribution is computed from SPARC mass-model components:

```text
V_bar^2 = |V_gas| V_gas + Upsilon_disk V_disk^2 + Upsilon_bulge V_bul^2
```

with:

```text
Upsilon_disk = 0.5
Upsilon_bulge = 0.7
```

The final prediction is:

```text
V_Flux^2(r) = V_bar^2(r) + G eta_eff M_bar / r * x^1.36 / (1+x)^0.36
```


