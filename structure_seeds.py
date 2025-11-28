# structure_seeds.py
import numpy as np
import vispy
vispy.set_log_level('debug')
vispy.use('pyqt5')
import vispy
vispy.set_log_level('debug')
vispy.use('pyqt5')


def generate_bh_seeds(n_bh, box_size, M_min=1.0, M_max=50.0, rng=None):
    rng = np.random.default_rng(rng)
    positions = rng.uniform(-box_size, box_size, size=(n_bh, 3))
    masses = rng.uniform(M_min, M_max, size=n_bh)
    return positions, masses
