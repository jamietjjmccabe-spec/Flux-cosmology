import numpy as np
from cmbr import (N_vals, a_vals, w_sigma_vals, Omega_sigma_0)  # import from your model

# We only care about a range where dark energy is relevant, say a in [0.2, 1]
mask = (a_vals >= 0.2) & (a_vals <= 1.0)
a_fit = a_vals[mask]
w_fit = w_sigma_vals[mask]

# CPL form: w(a) = w0 + wa*(1 - a)
# Put into linear regression form: w = (w0 + wa) + (-wa)*a
# Unknowns: A = w0 + wa, B = -wa
A = np.vstack([np.ones_like(a_fit), a_fit]).T   # shape (n,2)

# Solve least squares A*[A0, B0]^T ≈ w_fit
coeff, *_ = np.linalg.lstsq(A, w_fit, rcond=None)
A0, B0 = coeff

wa = -B0
w0 = A0 + wa

print("Fitted CPL parameters:")
print("w0 =", w0)
print("wa =", wa)
print("Omega_fld =", Omega_sigma_0)
