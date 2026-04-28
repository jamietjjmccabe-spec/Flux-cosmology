# Flux Weak-Field Limit

## Purpose

This note defines a first falsifiable weak-field limit for Flux Cosmology.

The goal is not to claim a completed replacement for General Relativity or dark matter. The goal is narrower and testable:

> Define a screened weak-field correction whose strength is controlled by the local available coherence fraction, then test whether a single parameter region can remain suppressed in the Solar System while becoming significant in galactic outskirts.

In this limit, standard Newtonian gravity is recovered in highly actualized, high-coherence environments, while an additional flux-coherence gradient may appear in diffuse, low-density, low-acceleration regions.

---

# 1. Weak-Field Potential

In the weak-field, slow-motion regime, Flux Cosmology is approximated by a Newtonian potential plus a screened flux-coherence correction:

\[
\Phi_{\rm total}(r)=\Phi_N(r)+\Phi_\sigma(r)
\]

with

\[
\Phi_N(r)=-\frac{GM}{r}
\]

and

\[
\Phi_\sigma(r)=-\frac{GM}{r}\alpha C_\sigma(r)e^{-r/\lambda_\sigma}.
\]

Therefore,

\[
\Phi_{\rm total}(r)
=
-\frac{GM}{r}
\left[
1+\alpha C_\sigma(r)e^{-r/\lambda_\sigma}
\right].
\]

Here:

- \(G\) is Newton's gravitational constant,
- \(M\) is the baryonic source mass,
- \(r\) is radial distance from the source,
- \(\alpha\) is the flux-coupling strength,
- \(\lambda_\sigma\) is the flux-coherence envelope length,
- \(C_\sigma(r)\) is the local coherence-availability factor.

The physical radial acceleration is

\[
g(r)=-\frac{d\Phi_{\rm total}}{dr}.
\]

This can be split into the Newtonian part and the flux correction:

\[
g(r)=g_N(r)+g_\sigma(r)
\]

where

\[
g_N(r)=\frac{GM}{r^2}.
\]

---

# 2. Flux-Native Definition of \(C_\sigma\)

The key Flux-native definition is that \(C_\sigma\) is not an arbitrary interpolation function.

It directly represents the local fraction of unrealized coherence, or available flux:

\[
C_\sigma(r)
\equiv
\frac{S_{\max}(r)-S(r)}{S_{\max}(r)}.
\]

Equivalently:

\[
C_\sigma(r)
=
\frac{\text{available flux}}{\text{entropy capacity}}.
\]

Where:

- \(S_{\max}(r)\) is the local entropy capacity,
- \(S(r)\) is the actualized entropy already occupying that capacity,
- \(S_{\max}(r)-S(r)\) is the remaining unspent coherence potential.

This makes \(C_\sigma\) the weak-field version of the same closed-system bookkeeping already used in the MIP/cosmology scripts.

In the global cosmology language:

\[
\Phi_{\rm available}=S_{\max}-S
\]

and

\[
\text{flux fraction}
=
\frac{\Phi_{\rm available}}{S_{\max}}.
\]

The weak-field coherence factor is the local radial version of that same quantity.

---

# 3. Physical Interpretation

## 3.1 Highly Actualized Regions

When a local region is highly actualized and coherence-saturated,

\[
S(r)\rightarrow S_{\max}(r),
\]

then

\[
C_\sigma(r)\rightarrow 0.
\]

The flux correction vanishes:

\[
\Phi_\sigma(r)\rightarrow 0.
\]

Therefore,

\[
\Phi_{\rm total}(r)\rightarrow \Phi_N(r)
\]

and standard Newtonian/GR gravity is recovered.

This is the Solar System regime.

In Flux language:

> The local physical state has already spent most of its available coherence. There is little unactualized flux left to contribute an additional weak-field gradient.

---

## 3.2 Diffuse Coherence Envelopes

When a region is diffuse, low-density, or not fully actualized,

\[
S(r)<S_{\max}(r),
\]

then

\[
C_\sigma(r)>0.
\]

The flux correction becomes nonzero:

\[
\Phi_\sigma(r)\neq 0.
\]

This produces an additional weak-field acceleration:

\[
g_\sigma(r)>0
\]

over some range of radii.

This is the galactic-envelope regime.

In Flux language:

> The region still contains available coherence capacity. That unspent capacity behaves as a weak-field flux reservoir, producing an additional effective gravitational gradient.

---

# 4. Practical Screening Approximation

The primary definition remains:

\[
C_\sigma(r)=\frac{S_{\max}(r)-S(r)}{S_{\max}(r)}.
\]

However, for numerical testing, this can be approximated by an acceleration-screening function:

\[
C_\sigma(r)
\approx
\frac{1}
{
1+\left(
\frac{a_N(r)}{a_\sigma}
\right)^n
}.
\]

The Newtonian acceleration is

\[
a_N(r)=\frac{GM}{r^2}.
\]

The flux acceleration scale is taken to be tied to the cosmic background:

\[
a_\sigma=\eta cH_0.
\]

Where:

- \(c\) is the speed of light,
- \(H_0\) is the present Hubble expansion rate,
- \(\eta\) is a dimensionless scaling parameter,
- \(n\) controls the sharpness of the screening transition.

This gives:

\[
C_\sigma(r)
=
\frac{1}
{
1+
\left(
\frac{GM}{a_\sigma r^2}
\right)^n
}.
\]

---

# 5. Limiting Behaviour

## 5.1 High-Acceleration Limit

In high-acceleration environments,

\[
a_N\gg a_\sigma.
\]

Then,

\[
\left(\frac{a_N}{a_\sigma}\right)^n\gg 1
\]

and therefore

\[
C_\sigma(r)\approx
\left(\frac{a_\sigma}{a_N}\right)^n
\ll 1.
\]

So the flux correction is heavily screened.

This is required for:

- Mercury,
- Earth,
- the Moon,
- inner planetary orbits,
- laboratory-scale gravity tests.

The model must satisfy:

\[
g_\sigma \ll g_N.
\]

---

## 5.2 Low-Acceleration Limit

In low-acceleration environments,

\[
a_N\ll a_\sigma.
\]

Then,

\[
\left(\frac{a_N}{a_\sigma}\right)^n\ll 1
\]

and therefore

\[
C_\sigma(r)\approx 1.
\]

So the flux correction becomes available:

\[
\Phi_\sigma(r)
\approx
-\frac{GM}{r}\alpha e^{-r/\lambda_\sigma}.
\]

This is the intended galactic-envelope regime.

The model becomes interesting only if it can produce:

\[
g_\sigma \sim g_N
\]

in galactic outskirts, while still satisfying Solar System suppression.

---

# 6. Exact Weak-Field Acceleration

The flux potential is

\[
\Phi_\sigma(r)
=
-\frac{GM}{r}\alpha C_\sigma(r)e^{-r/\lambda_\sigma}.
\]

The physical acceleration contribution is

\[
g_\sigma(r)
=
-\frac{d\Phi_\sigma}{dr}.
\]

Carrying out the derivative gives:

\[
g_\sigma(r)
=
\alpha GM e^{-r/\lambda_\sigma}
\left[
\frac{C_\sigma(r)}{r^2}
-
\frac{1}{r}\frac{dC_\sigma}{dr}
+
\frac{C_\sigma(r)}{\lambda_\sigma r}
\right].
\]

Therefore, the total acceleration is

\[
g_{\rm total}(r)
=
\frac{GM}{r^2}
+
\alpha GM e^{-r/\lambda_\sigma}
\left[
\frac{C_\sigma(r)}{r^2}
-
\frac{1}{r}\frac{dC_\sigma}{dr}
+
\frac{C_\sigma(r)}{\lambda_\sigma r}
\right].
\]

This expression should be used for numerical testing.

The derivative term is essential.

If \(C_\sigma(r)\) is treated as constant, the approximation becomes:

\[
g_\sigma(r)
\approx
g_N(r)\alpha C_\sigma(r)e^{-r/\lambda_\sigma}
\left[
1+\frac{r}{\lambda_\sigma}
\right].
\]

But this approximation is incomplete because \(C_\sigma\) is precisely the spatially varying coherence-availability field.

---

# 7. Derivative of the Screening Function

For the acceleration-screened approximation,

\[
C_\sigma(r)
=
\frac{1}
{
1+
\left(
\frac{a_N(r)}{a_\sigma}
\right)^n
}
\]

with

\[
a_N(r)=\frac{GM}{r^2},
\]

the radial derivative is

\[
\frac{dC_\sigma}{dr}
=
\frac{2n}{r}C_\sigma(r)\left[1-C_\sigma(r)\right].
\]

This allows the exact weak-field acceleration to be computed without numerical differentiation.

Substituting into the flux acceleration gives:

\[
g_\sigma(r)
=
\alpha GM e^{-r/\lambda_\sigma}
\left[
\frac{C_\sigma}{r^2}
-
\frac{2n}{r^2}C_\sigma(1-C_\sigma)
+
\frac{C_\sigma}{\lambda_\sigma r}
\right].
\]

Or equivalently:

\[
g_\sigma(r)
=
\alpha GM e^{-r/\lambda_\sigma}
\left[
\frac{C_\sigma}{r^2}
\left(
1-2n(1-C_\sigma)
\right)
+
\frac{C_\sigma}{\lambda_\sigma r}
\right].
\]

This is the preferred numerical form for the first parameter scan.

---

# 8. Why This Is Not Just MOND

MOND introduces a phenomenological acceleration scale and modifies dynamics below that scale.

Flux Cosmology instead interprets the transition as a coherence-availability effect:

\[
C_\sigma
=
\frac{\text{available flux}}{\text{entropy capacity}}.
\]

The acceleration-screened form is not the fundamental definition. It is a practical approximation to a deeper closed-flux bookkeeping rule.

The conceptual chain is:

\[
S_{\max}-S
\rightarrow
\text{available flux}
\rightarrow
C_\sigma
\rightarrow
\Phi_\sigma
\rightarrow
g_\sigma.
\]

Thus, the weak-field correction is not inserted by hand as a galaxy-fitting switch.

It inherits its meaning from the same entropy-capacity framework used in:

- MIP actualization,
- horizon saturation,
- cosmic expansion,
- flux exhaustion,
- jet and black-hole toy models,
- background scalar-field evolution.

---

# 9. Relation to the Background Scalar Field

The flux acceleration scale is tied to the cosmological background:

\[
a_\sigma=\eta cH_0.
\]

This is motivated by the idea that the weak-field coherence threshold should be set by the large-scale flux background, not by a purely local arbitrary constant.

In the scalar-field cosmology model, the background field \(\sigma\) evolves according to a potential:

\[
V(\sigma)
=
V_0
\left(
1-e^{-\sigma/\mu}
\right)^2.
\]

The scalar behaves as a residual flux reservoir.

At late times, when the field remains high on the plateau, the flux sector behaves similarly to a dark-energy-like component.

The weak-field parameter \(a_\sigma\) should eventually be derived more directly from the background scalar quantities:

\[
\sigma(a), \quad V(\sigma), \quad \Omega_\sigma(a), \quad w_\sigma(a).
\]

For the first numerical scan, however, the phenomenological link

\[
a_\sigma=\eta cH_0
\]

is sufficient.

---

# 10. Falsifiability Conditions

The weak-field model is useful only if a single reasonable parameter region satisfies both local and galactic constraints.

## 10.1 Solar System Constraint

In the Solar System:

\[
g_\sigma \ll g_N.
\]

A first crude requirement is:

\[
\left|\frac{g_\sigma}{g_N}\right|_{\rm Moon}
\ll 1
\]

and

\[
\left|\frac{g_\sigma}{g_N}\right|_{\rm Mercury}
\ll 1.
\]

Initial sandbox thresholds may use:

\[
\left|\frac{g_\sigma}{g_N}\right|<10^{-12}
\]

or tighter.

A successful model should not noticeably disturb:

- lunar laser ranging,
- Mercury perihelion constraints,
- planetary ephemerides,
- Earth-Moon tidal recession modelling.

The existing lunar orbit script should therefore be treated as a constraint sandbox, not as an independent prediction of lunar recession.

---

## 10.2 Galactic Constraint

In galactic outskirts:

\[
g_\sigma \sim g_N
\]

or at least large enough to increase circular velocity.

The circular velocity is

\[
v_c(r)=\sqrt{r g_{\rm total}(r)}.
\]

The Newtonian-only value is

\[
v_N(r)=\sqrt{r g_N(r)}.
\]

A useful diagnostic is the boost factor:

\[
B_v(r)=\frac{v_c(r)}{v_N(r)}.
\]

Since

\[
v_c^2=r(g_N+g_\sigma),
\]

the boost satisfies:

\[
B_v(r)
=
\sqrt{
1+\frac{g_\sigma}{g_N}
}.
\]

A 30 percent velocity boost requires:

\[
B_v\approx 1.3
\]

which implies:

\[
\frac{g_\sigma}{g_N}
\approx 0.69.
\]

So a rough target for galaxy outskirts is:

\[
0.3 \lesssim \frac{g_\sigma}{g_N} \lesssim 1.
\]

---

# 11. First Parameter Set to Explore

The initial free parameters are:

\[
\alpha,\eta,\lambda_\sigma,n.
\]

Where:

- \(\alpha\) controls coupling strength,
- \(\eta\) sets the acceleration threshold through \(a_\sigma=\eta cH_0\),
- \(\lambda_\sigma\) controls the coherence-envelope range,
- \(n\) controls screening sharpness.

A first coarse scan can explore:

\[
\alpha \in [10^{-3},10]
\]

\[
\eta \in [10^{-4},1]
\]

\[
\lambda_\sigma \in [10^{18},10^{23}] \ \text{m}
\]

\[
n \in [1,4].
\]

The model is viable only if overlap exists between:

1. Solar System suppression,
2. galactic enhancement,
3. non-pathological acceleration profiles.

---

# 12. Environment-Dependent Coherence Length

A fixed \(\lambda_\sigma\) may be too rigid.

The next extension is to make the coherence-envelope length environmentally dependent:

\[
\lambda_\sigma \rightarrow \lambda_\sigma(\Sigma_b)
\]

where \(\Sigma_b\) is baryonic surface density.

A possible first closure is:

\[
\lambda_\sigma
=
\lambda_0
\left(
\frac{\Sigma_b}{\Sigma_{\rm ref}}
\right)^{-q}.
\]

Where:

- \(\lambda_0\) is a reference coherence length,
- \(\Sigma_{\rm ref}\) is a reference baryonic surface density,
- \(q\) is a small positive scaling exponent.

This encodes the idea that dense actualized regions shrink the coherence envelope, while diffuse regions allow it to extend.

In Flux language:

> Dense baryonic systems spend coherence capacity quickly; diffuse systems leave larger unresolved coherence envelopes.

---

# 13. Numerical Implementation

A minimal implementation should include:

```python
import numpy as np

G = 6.67430e-11
c = 299792458.0
H0 = 70_000 / 3.085677581e22
a_cosmic = c * H0

def C_sigma(r, M, a_sigma, n):
    """
    Coherence availability:
    C_sigma = available flux / entropy capacity.

    Practical acceleration-screened approximation:
    C = 1 / [1 + (a_N / a_sigma)^n]
    """
    r = np.asarray(r, dtype=float)
    a_N = G * M / r**2
    return 1.0 / (1.0 + (a_N / a_sigma)**n)


def dC_sigma_dr(r, M, a_sigma, n):
    """
    Analytic derivative of the acceleration-screened C_sigma.
    """
    r = np.asarray(r, dtype=float)
    C = C_sigma(r, M, a_sigma, n)
    return (2.0 * n / r) * C * (1.0 - C)


def g_flux_exact(r, M, alpha, a_sigma, lambda_sigma, n):
    """
    Exact weak-field radial acceleration including dC/dr.

    Phi_sigma = -GM/r * alpha * C_sigma(r) * exp(-r/lambda_sigma)

    g_sigma = -dPhi_sigma/dr
    """
    r = np.asarray(r, dtype=float)

    C = C_sigma(r, M, a_sigma, n)
    dCdr = dC_sigma_dr(r, M, a_sigma, n)
    E = np.exp(-r / lambda_sigma)

    g_N = G * M / r**2

    g_sigma = alpha * G * M * E * (
        C / r**2
        - dCdr / r
        + C / (lambda_sigma * r)
    )

    g_total = g_N + g_sigma

    return g_total, g_N, g_sigma, C14. Solar System Diagnostic

For a Solar System body at radius r, compute:

R
solar
	​

(r)=
	​

g
N
	​

(r)
g
σ
	​

(r)
	​

	​

.

The model passes the crude screen if:

R
solar
	​

(r)<ϵ
solar
	​

.

Example values:

M_sun = 1.98847e30

AU = 1.495978707e11
r_mercury = 0.387 * AU
r_earth = 1.0 * AU
r_moon = 3.844e8

A useful first test:

epsilon_solar = 1e-12

The scan should reject any parameter set where the Solar System correction is too large.

15. Galactic Diagnostic

For a crude galaxy test, use a baryonic mass scale:

M_gal = 6e10 * M_sun

and galactic radii:

kpc = 3.085677581e19
r_inner = 5 * kpc
r_outer = 30 * kpc

Compute:

g_total, g_N, g_sigma, C = g_flux_exact(
    r_outer,
    M_gal,
    alpha,
    a_sigma,
    lambda_sigma,
    n
)

boost = np.sqrt(g_total / g_N)

A crude useful condition is:

boost > 1.2

or more aggressively:

boost > 1.3

This corresponds to meaningful flattening support in the outskirts.

16. Interpretation of a Successful Overlap

If a single parameter set gives:

	​

g
N
	​

g
σ
	​

	​

	​

Solar
	​

≪1

and

g
N
	​

g
σ
	​

	​

galaxy
	​

∼0.3−1,

then the weak-field Flux model has a viable numerical window.

That does not prove the theory.

It would show only that the framework is not immediately ruled out by the simplest local-versus-galactic comparison.

The correct claim would be:

A screened Flux weak-field limit has a non-empty parameter region where Solar System corrections are suppressed while galactic-envelope corrections become dynamically relevant.

17. Interpretation of a Failed Overlap

If no overlap exists, then at least one assumption must be revised.

Possible failure modes:

Fixed λ
σ
	​

 is too rigid.
The acceleration-screened approximation is too simple.
n makes the transition too sharp.
a
σ
	​

=ηcH
0
	​

 is too restrictive.
The model requires explicit surface-density dependence.
C
σ
	​

 must be computed from actual S
max
	​

−S, not approximated through a
N
	​

.
The weak-field correction is too constrained to explain galaxy-scale missing gravity.

A failed scan would still be useful because it would identify which closure assumption breaks first.

18. Next Development Steps
Step 1 — Add the weak-field note

Save this file as:

flux_weak_field_note.md
Step 2 — Add a parameter scan

Create:

parameter_collision.py

The script should scan over:

α,η,λ
σ
	​

,n.

It should reject parameter sets that fail Solar System suppression and keep those that produce galactic velocity boosts.

Step 3 — Mark lunar script as sandbox

Rename or document the lunar script as:

constraint_sandbox_lunar.py

It should be explicitly described as a constraint calibration tool, because it calibrates tidal torque to match the observed lunar recession rate.

It does not independently predict the 3.8 cm/year recession.

Step 4 — Add environmental λ
σ
	​


Introduce:

λ
σ
	​

(Σ
b
	​

)

for galaxy testing.

Step 5 — Replace acceleration-screened C
σ
	​


Eventually replace:

C
σ
	​

(r)≈
1+(a
N
	​

/a
σ
	​

)
n
1
	​


with an explicit local entropy-capacity computation:

C
σ
	​

(r)=
S
max
	​

(r)
S
max
	​

(r)−S(r)
	​

.

This will make the weak-field law fully Flux-native instead of partly phenomenological.

19. Minimal Claim

The minimal defensible claim is:

Flux Cosmology admits a weak-field screened-coherence limit in which the correction to Newtonian gravity is controlled by the available coherence fraction C
σ
	​

=(S
max
	​

−S)/S
max
	​

. In this limit, highly actualized regions suppress the correction, while diffuse coherence envelopes may produce an additional effective gravitational gradient. The model is numerically falsifiable by searching for overlap between Solar System suppression and galactic rotation support.

20. Stronger Future Claim

The stronger future claim, only if the scans support it, would be:

The same entropy-capacity bookkeeping used in MIP and background Flux Cosmology may also generate a viable screened weak-field correction, allowing galaxy-scale missing-gravity behaviour without requiring the correction to appear in high-actualization environments such as the Solar System.

This claim should not be made until the parameter scan is complete.

21. Summary

The weak-field Flux limit is:

Φ
total
	​

(r)=−
r
GM
	​

[1+αC
σ
	​

(r)e
−r/λ
σ
	​

].

The Flux-native coherence factor is:

C
σ
	​

(r)=
S
max
	​

(r)
S
max
	​

(r)−S(r)
	​

.

The practical approximation is:

C
σ
	​

(r)≈
1+(
a
σ
	​

a
N
	​

(r)
	​

)
n
1
	​

.

The exact weak-field acceleration is:

g
total
	​

(r)=
r
2
GM
	​

+αGMe
−r/λ
σ
	​

[
r
2
C
σ
	​

(r)
	​

−
r
1
	​

dr
dC
σ
	​

	​

+
λ
σ
	​

r
C
σ
	​

(r)
	​

].

The acceleration-screened derivative is:

dr
dC
σ
	​

	​

=
r
2n
	​

C
σ
	​

(1−C
σ
	​

).

The core falsifiability test is:

g
σ
	​

≪g
N
	​


in the Solar System, while

g
σ
	​

∼g
N
	​


in galactic outskirts.

If both can be satisfied by the same reasonable parameter region, the weak-field Flux model survives its first numerical collision test.

If not, the closure must be revised.