# Flux Cosmology & Quantum Coherence Model

**Author:** Jamie McCabe  
**License:** Code – MIT, Paper – CC BY 4.0

This repository contains the research code and main paper for a **closed-flux cosmology** and **flux–coherence gravity** framework. The model treats cosmic acceleration, structure formation, and compact-object feedback as emergent phenomena of a flux–scalar field coupled to gravity and entropy.

---

## Contents

### Core Cosmology

- `flux_universe.py`  
  Baseline closed-flux cosmology background integrator: scale factor, flux / entropy variables, and effective equation of state.

- `M_I_P_Universe model using closed quantum flux cosmology.py`  
  “Measurement–Induced Projection (MIP)” version of the universe model, including explicit hardened / flux state variables.

- `cmbr.py`, `cmbr4.py`, `cmbr_universe_backend.py`  
  Toy CMB / background comparison modules (sound horizon, distance measures, growth index style outputs in the flux model).

### Local / Solar-System Tests

- `Lunar Orbit Using Quantum Flux closed model.py`  
  Lunar orbital evolution in the flux framework (tests against observed recession and orbital period drift).

- `Solarfluxmap.py`  
  Solar-neighbourhood flux mapping / visualisation tools.

### Compact Objects & Jets

- `Blackhole.py`, `M_I_P_Blackhole.py`  
  Flux-modified black-hole potential and horizon-scale toy models.

- `M_I_P_Eccretion.py`, `flux_jet_full.py`, `jet_3d.py`, `mip_jet_2d.py`, `jetsim.py`  
  Accretion + jet toy models, including flux-driven jet power and feedback.

### Structure Formation & Galaxies

- `flux_spiral_galaxy.py`, `spiral_galaxy_3d.py`, `spiral_galaxy_flux_coupled.py`  
  Flux-driven spiral galaxy / disc formation models.

- `structure_seeds.py`, `supercluster_flux.py`, `universe3d.py`, `supernova.py`, `quazar.py`, `redgiant.py`  
  Stellar evolution, supernova, quasar, and supercluster-scale flux experiments.

### Mapping & Visualisation

- `map2.py` – `map29.py`, `maps22.py`, `milkyway1.py`, `milky_way.py`, `milky_way2.py`, `milky way formation.png`  
  2D / 3D flux-field, galaxy and large-scale-structure maps used to generate the figures in the paper.

### Paper & Notes

- `The flux--scalar cosmological model.tex`  
- `The flux--scalar cosmological model.pdf`  
- `Information.txt`, `Gravity magnitude from ∇V.txt`, `start.txt`  

The TeX source and compiled PDF of the main paper, plus supporting notes.

---

## Quick Start

### 1. Requirements

You will need:

- Python 3.10+ (3.8+ should work)
- `numpy`
- `scipy`
- `matplotlib`

Install dependencies with:

```bash
pip install numpy scipy matplotlib
