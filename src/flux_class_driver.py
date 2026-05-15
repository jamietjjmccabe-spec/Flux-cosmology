from classy import Class
from fit_w_cpl import w0, wa, Omega_sigma_0  # from step 2

# 1. Define cosmological parameters
params = {
    # Base cosmology (tweak to taste / Planck)
    'h': 0.67,
    'omega_b': 0.0224,      # Ω_b h^2
    'omega_cdm': 0.12,      # Ω_c h^2
    'A_s': 2.1e-9,
    'n_s': 0.965,
    'tau_reio': 0.054,

    # Turn off vanilla Λ, we replace it with "fld"
    'Omega_Lambda': 0.0,

    # Your flux field as a CPL dark energy fluid
    'Omega_fld': Omega_sigma_0,  # from your model
    'w0_fld': w0,
    'wa_fld': wa,
    'use_ppf': 'yes',            # stable crossing near w=-1

    # Outputs
    'output': 'tCl,lCl,pCl',
    'l_max_scalars': 2500
}

# 2. Run CLASS
cosmo = Class()
cosmo.set(params)
cosmo.compute()

# 3. Get C_ell^TT spectrum
cls = cosmo.lensed_cl(2500)  # dict with 'tt', 'ee', 'te', ...
ell = cls['ell']
clTT = cls['tt']

# Optional: save to file for plotting elsewhere
import numpy as np
np.savetxt("flux_TT_cls.txt", np.column_stack([ell, clTT]))

# You can also print the acoustic scale, etc.
theta_s = cosmo.theta_s()  # sound horizon angle
print("theta_s =", theta_s)

cosmo.struct_cleanup()
cosmo.empty()

