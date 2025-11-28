# flux_background_api.py
import vispy
vispy.set_log_level('debug')
vispy.use('pyqt5')
import numpy as np
from cmbr4 import (
    N_vals, a_vals, H_vals, sigma_vals,
    Omega_sigma_vals, rho_m_vals, rho_tot_vals
)

# Reconstruct an effective cosmic time t(N) from background
t_vals = np.zeros_like(N_vals)
for i in range(1, len(N_vals)):
    dN = N_vals[i] - N_vals[i-1]
    Hmid = 0.5 * (H_vals[i] + H_vals[i-1])
    t_vals[i] = t_vals[i-1] + dN / max(Hmid, 1e-12)

# Effective gravity from flux fraction
g1 = -1e-6
G_eff_vals = 1.0 * (1.0 + g1 * Omega_sigma_vals)

def background():
    return {
        "t": t_vals,
        "N": N_vals,
        "a": a_vals,
        "H": H_vals,
        "sigma": sigma_vals,
        "Omega_sigma": Omega_sigma_vals,
        "G_eff": G_eff_vals,
        "Omega_m": rho_m_vals / rho_tot_vals,
    }
