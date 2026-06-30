# Active baseline and freeze boundaries - June 2026

## Frozen pre-integration baseline

The following elements are held fixed until a full Boltzmann-code implementation is available:

- Newtonian-gauge perturbation variables \(\{\delta V,\delta\rho_{\rm reg},q_{\rm reg}\}\);
- diagnostics \(\{\epsilon_Q,\epsilon_{\rm birth},R_{\rm force}\}\);
- explicit conservation through the build/binding reservoir;
- no artificial numerical floors that hide a singularity or instability;
- the target cosmology uses \(\Omega_{\rm cdm}=0\).

Background-only agreement, a fitted CPL fluid, or an acoustic-scale integral is not sufficient to claim CMB viability.

## Required integration sequence

1. Implement the exact background equations in an unmodified writable CLASS source tree.
2. Verify background closure and energy transfer at \(\Omega_{\rm cdm}=0\).
3. Implement the exact Newtonian-gauge perturbation equations.
4. Run the full physical-wavenumber stability scan.
5. Compute TT, TE, EE, lensing, matter power, and growth observables.
6. Compare against a matched \(\Lambda\)CDM baseline with the same nuisance and primordial assumptions.
7. Reject the model if stability or observational agreement requires hidden floors, scale-by-scale retuning, or unphysical initial conditions.

## Scalar-mediated local-force path

The active local path treats the scalar as a mediator:

\[
\mathbf a_\phi=-q_b\nabla\phi.
\]

It is subject to the simultaneous requirements of galactic amplitude, local suppression, and healthy field dynamics. No parameter point is accepted until all three are demonstrated in one model.

## Antimatter constraint

Matter and antimatter contribute positively to the gravitational energy-history ledger. Any odd matter-antimatter quantity is tracked separately as a current/asymmetry variable and is not used as negative gravitational mass.

## Interpretation rule

Repository documents and plots must distinguish:

- **derived equations**;
- **phenomenological closures**;
- **toy simulations**;
- **heuristic analogies**;
- **observational fits**.

No toy output should be described as confirmation of the full theory.
