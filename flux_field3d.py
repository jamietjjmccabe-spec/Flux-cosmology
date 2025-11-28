import numpy as np

def f_noise(rng, shape, freq):
    """Simple periodic noise – lightweight substitute for Perlin."""
    return rng.normal(0, 1, shape) * freq

def build_flux_field(box_size, grid_n=64, amp=1.5, rng=None):
    """
    Multi-octave fractal flux field producing sharp filaments and voids.
    """
    rng = np.random.default_rng(rng)

    xs = np.linspace(-box_size, box_size, grid_n)
    ys = np.linspace(-box_size, box_size, grid_n)
    zs = np.linspace(-box_size, box_size, grid_n)
    dx = xs[1] - xs[0]

    # --- Multi-octave 3D field ---
    field = np.zeros((grid_n, grid_n, grid_n))

    # Octaves: low → medium → high frequency structure
    octaves = [
        (1.0, 0.5),   # (frequency_scale, amplitude_scale)
        (2.0, 0.35),
        (4.0, 0.25),
        (8.0, 0.12),
    ]

    for freq, weight in octaves:
        noise = rng.normal(0, 1, field.shape)
        # Smooth the noise (simple blur)
        kernel = np.array([1, 4, 6, 4, 1]) / 16

        def smooth(a, axis):
            return np.apply_along_axis(
                lambda m: np.convolve(m, kernel, mode='same'),
                axis,
                a
            )

        for axis in range(3):
            noise = smooth(noise, axis)

        field += weight * noise

    # Normalise
    field -= field.mean()
    field /= (field.std() + 1e-9)
    field *= amp

    # Non-linear shaping → sharp ridges + deep valleys
    field = np.tanh(1.2 * field)

    # --- Gradient field ---
    dsdx, dsdy, dsdz = np.gradient(field, dx, dx, dx, edge_order=2)

    return {
        "xs": xs,
        "ys": ys,
        "zs": zs,
        "sigma": field,
        "grad_x": dsdx,
        "grad_y": dsdy,
        "grad_z": dsdz,
        "dx": dx,
        "box": box_size,
        "N": grid_n,
    }
