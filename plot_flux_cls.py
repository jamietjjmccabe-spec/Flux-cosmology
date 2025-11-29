# plot_flux_cls.py
import numpy as np
import matplotlib.pyplot as plt

ell, clTT = np.loadtxt("flux_TT_cls.txt", unpack=True)

plt.figure()
plt.semilogx(ell[2:], ell[2:]*(ell[2:]+1)*clTT[2:]/(2.0*np.pi))
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell^{TT}/2\pi$')
plt.title("CMB TT spectrum with flux-driven dark energy")
plt.show()
