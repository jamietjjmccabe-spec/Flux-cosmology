# Flux-MIP V2.1 frozen theory specification

## Status

Phenomenological open-system benchmark. Not empirically established and not derived from a fundamental action.

## Source field

```math
G_\ell(r)=\frac{\ell^2}{r^2+\ell^2},
```

```math
D_{\rm reg}(\mathbf x)=\frac{\kappa_{\rm reg}}{\rho_*V_*}\int \rho_b(\mathbf r')G_\ell(|\mathbf x-\mathbf r'|)d^3r',
```

```math
S_{\rm reg}=-\ell^2\nabla^2D_{\rm reg},\qquad F_{\rm reg}=\tanh(S_{\rm reg}/S_*).
```

The curvature rule rejects constant and affine precursor backgrounds. It does not permit post-hoc adaptive saturation.

## Open-system dynamics

```math
\dot\rho=-\frac{i}{\hbar}[H_0+\hbar\beta_{\rm MIP}F_{\rm reg}(\hat x),\rho]
-\frac{\gamma_{\rm MIP}}{2}[F_{\rm reg}(\hat x),[F_{\rm reg}(\hat x),\rho]].
```

For two narrow paths,

```math
I_\phi=\int\Delta Fdt,\qquad X=\frac12\int(\Delta F)^2dt,
```

```math
\delta\phi_{\rm MIP}=\beta_{\rm MIP}I_\phi,\qquad V/V_{\rm standard}=e^{-\gamma_{\rm MIP}X}.
```

## Frozen values

- `ell = 0.010 m`
- `p = 2`
- `rho_star = 1000 kg m^-3`
- `V_star = 4.188790204786391e-6 m^3`
- `S_star = 10`
- `kappa_reg = 1`
- ordinary matter factor `= 1`

Benchmark signature: `f819067c8feca48e6d0a51b790e01f5a1be599f4afa22c6accc46f3967827670`

## Claim boundary

A positive result would establish an unexplained source-correlated open-system anomaly under the declared controls. It would not automatically establish MIP ontology, metric-memory dark matter, or a complete theory of quantum gravity.
