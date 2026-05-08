# Fibonacci / 3D Universe Flux Model Copy

This zip is a runnable subset of the Flux Cosmology repository focused on the 3D universe visualisation layer.

## Main entry point

```bash
python universe3d.py
```

## Included support modules

- `flux_field3d.py` — builds the 3D flux / sigma field and gradients.
- `structure_seeds.py` — generates black-hole / structure seeds.
- `flux_background_api.py` — bridges the background cosmology into the 3D visual layer.
- `cmbr4.py`, `cmbr.py`, `cmbr_universe_backend.py` — background CMB / scalar-field model support.
- `spiral_galaxy_3d.py`, `spiral_galaxy_flux_coupled.py`, `jet_3d.py` — related 3D visual/toy models.

## Dependencies

```bash
pip install numpy scipy matplotlib vispy pyqt5
```

## Note

This is packaged as the Fibonacci/3D universe model copy, but the current local files use the flux/sigma/phi notation. If the Fibonacci-gated version is in another local branch, send that branch/file and it can be folded into this package.
