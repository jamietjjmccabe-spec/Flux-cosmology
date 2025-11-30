import numpy as np
import matplotlib.pyplot as plt

try:
    import imageio.v2 as imageio   # for GIF export (optional)
    HAS_IMAGEIO = True
except ImportError:
    HAS_IMAGEIO = False

# ==========================================================
#  FLUX DRIVER HOOKS
# ==========================================================
#
# In your main project you can replace `flux_driver(t)` with
# an interpolation from your cosmology / jet models, e.g.:
#
# from M_I_P_Universe_model_using_closed_quantum_flux_cosmology import run_unified
# data = run_unified()
# t_series = data["t"]
# phi_series = data["flux_frac"]
# flux = np.interp(t, t_series, phi_series)
#
# For now we just use a simple toy function 0..1..0.

def flux_driver(t, t_max):
    """Toy flux driver in [0,1] – replace with real flux_frac(t)."""
    x = t / max(t_max, 1e-9)
    return np.clip(4 * x * (1 - x), 0.0, 1.0)


# ==========================================================
#  2D RED SPIDER SLICE (r,z) WITH VELOCITY
# ==========================================================

def red_spider_slice(r, z, t, t_max, flux):
    """
    Build a 2D emissivity + LOS velocity slice for given time t and flux driver.

    r, z : 1D arrays
    t    : current time
    t_max: total time (for scaling)
    flux : 0..1 driver from your flux cosmology
    """
    R, Z = np.meshgrid(r, z)
    NR, NZ = len(r), len(z)

    # --- global scaling from time / flux ---
    # lobe extent grows with time, turbulence with flux
    v_lobe = 1.4 + 0.6 * flux
    t_eff  = 1.8 + 0.8 * (t / max(t_max, 1e-9))
    z_extent = v_lobe * t_eff

    waist_radius = 0.35
    curvature    = 1.7
    shell_thick  = 0.03

    # base bipolar shell (smooth)
    r_base = waist_radius + (np.abs(Z) / z_extent) ** curvature

    # --- precession (S-shape) controlled by flux ---
    A_p     = 0.25 + 0.2 * flux
    omega_p = 0.45
    x_shift = A_p * np.sin(omega_p * Z) + 0.1 * (Z / z_extent)
    R_pre   = R - x_shift

    # --- fractal / turbulent perturbations ---
    def fractal_noise(Z, octaves=5):
        noise = np.zeros_like(Z)
        amp   = 0.16 + 0.06 * flux
        freq  = 1.3
        for _ in range(octaves):
            phase1 = 2*np.pi*np.random.rand()
            phase2 = 2*np.pi*np.random.rand()
            noise += amp * np.sin(freq * np.abs(Z)       + phase1)
            noise += amp * 0.7 * np.cos(freq * np.abs(Z)**1.15 + phase2)
            amp  *= 0.55
            freq *= 1.8
        return noise

    fractal   = fractal_noise(Z)
    wrinkles  = 0.04 * np.sin(18 * Z + 2*np.pi*np.random.rand())
    turb      = fractal + wrinkles

    # strongest near waist, fade outward
    turb *= (1.2 * np.exp(-np.abs(Z)/3.0) * (1.0 + np.exp(-(Z**2)/0.7)))

    r_shell = r_base + turb

    shell = np.abs(np.abs(R_pre) - r_shell) < shell_thick
    shell &= (np.abs(Z) < z_extent)

    # --- inner hot bubble ---
    inner = (R_pre**2 + (0.55*Z)**2) < (0.55**2)

    # --- dust lane ---
    dust_depth = 0.75
    dust_sigma = 0.25
    dust_rad   = 0.9

    dust_mask = np.exp(-(Z**2) / (2*dust_sigma**2)) * (np.abs(R_pre) < dust_rad)
    dust_atten = 1.0 - dust_depth * dust_mask

    # --- jet knots (bullets) along axis ---
    base_knots = np.array([0.8, 1.6, 2.4, 3.2])
    # spread in time a bit
    shift = 0.3 * (t / max(t_max, 1e-9) - 0.5)
    knots_z = np.concatenate([base_knots + shift, -(base_knots + shift)])

    knot_rad_sig = 0.09
    knot_z_sig   = 0.13

    emiss_knots = np.zeros_like(R)
    for z0 in knots_z:
        gauss = (
            np.exp(-(R_pre**2) / (2*knot_rad_sig**2)) *
            np.exp(-((Z - z0)**2) / (2*knot_z_sig**2))
        )
        emiss_knots += gauss

    emiss_knots /= emiss_knots.max() + 1e-12
    emiss_knots *= 0.7

    # --- base emission ---
    emiss = np.zeros_like(R)
    emiss[shell] = 1.0
    emiss[inner] = 0.4 + 0.3 * flux
    emiss += emiss_knots
    emiss *= dust_atten
    emiss = np.clip(emiss, 0.0, None)
    emiss /= emiss.max() + 1e-12

    # --- LOS velocity (toy) ---
    v0 = 1.0 + 0.4 * flux
    v_z = v0 * np.tanh(Z / 1.0)
    v_z += 0.25 * v0 * np.sin(omega_p * Z)  # precession contribution
    v_norm = v_z / (np.max(np.abs(v_z)) + 1e-12)

    return emiss, v_norm, R_pre, Z


# ==========================================================
#  3D VOLUME FROM 2D SLICE (REVOLUTION)
# ==========================================================

def revolve_to_3d(emiss_2d, vel_2d, r, z, n_phi=128):
    """
    Revolve (r,z) slice around z-axis to build 3D voxel volumes.

    Returns:
      density  : (nz, n_phi, nr) emissivity
      velocity : same shape, LOS component (we keep v_z everywhere)
    """
    nz, nr = emiss_2d.shape
    phi = np.linspace(0, 2*np.pi, n_phi, endpoint=False)

    density = np.repeat(emiss_2d[:, None, :], n_phi, axis=1)
    vel     = np.repeat(vel_2d[:, None, :], n_phi, axis=1)

    # If later you want full 3D vector velocity:
    # you can compute vx, vy from v_z + some toroidal component.

    return density, vel, phi


# ==========================================================
#  SIMPLE RAYMARCH RENDERER
# ==========================================================

def raymarch_render(density, vel, cmap="inferno", v_cmap="seismic"):
    """
    Very simple front-to-back raymarch along 'phi' axis:
    produce 2D image in (r,z) plane from 3D volume.

    density, vel : (nz, n_phi, nr)
    """
    nz, n_phi, nr = density.shape

    # alpha compositing along phi
    alpha = density / (density.max() + 1e-12)
    # small optical depth
    alpha *= 0.12

    # emissive color scalar
    col_scalar = density.copy()

    # integrate along phi
    img = np.zeros((nz, nr))
    opacity = np.zeros((nz, nr))

    for i in range(n_phi):
        a = alpha[:, i, :]
        c = col_scalar[:, i, :]
        img   = img + (1 - opacity) * a * c
        opacity = opacity + (1 - opacity) * a

    img /= img.max() + 1e-12
    return img


# ==========================================================
#  ANIMATION / EXPORT
# ==========================================================

def generate_animation(
    n_frames=40,
    nr=400,
    nz=400,
    n_phi=128,
    save_prefix="red_spider",
    make_gif=True
):
    r = np.linspace(-2.0, 2.0, nr)
    z = np.linspace(-4.0, 4.0, nz)
    t_max = 1.0

    frame_paths = []

    for k in range(n_frames):
        t = t_max * k / max(n_frames - 1, 1)
        flux = flux_driver(t, t_max)

        emiss2d, v2d, R_pre, Z = red_spider_slice(r, z, t, t_max, flux)
        density3d, vel3d, phi = revolve_to_3d(emiss2d, v2d, r, z, n_phi=n_phi)

        # simple ray-marched image from volume
        img = raymarch_render(density3d, vel3d)

        # plot and save
        plt.figure(figsize=(5, 8))
        extent = [r.min(), r.max(), z.min(), z.max()]
        plt.imshow(
            img,
            origin="lower",
            extent=extent,
            cmap="inferno",
            aspect="auto"
        )
        plt.title(f"Red Spider Nebula – t={t:.2f}, flux={flux:.2f}")
        plt.xlabel("x (symmetric)")
        plt.ylabel("z")
        fname = f"{save_prefix}_frame_{k:03d}.png"
        plt.tight_layout()
        plt.savefig(fname, dpi=150)
        plt.close()
        frame_paths.append(fname)
        print(f"Saved {fname}")

    if make_gif and HAS_IMAGEIO:
        imgs = [imageio.imread(p) for p in frame_paths]
        gif_name = f"{save_prefix}.gif"
        imageio.mimsave(gif_name, imgs, fps=10)
        print(f"Saved GIF {gif_name}")
    elif make_gif:
        print("imageio not installed – skipping GIF export.")

    return frame_paths


# ==========================================================
#  VOXEL EXPORT FOR 3D APPS
# ==========================================================

def export_volume_npz(density, vel, filename="red_spider_volume.npz"):
    """
    Save density + velocity volumes to NPZ.

    Later you can:
      arr = np.load(filename)
      density = arr["density"]
      vel     = arr["vel"]
    """
    np.savez_compressed(filename, density=density, vel=vel)
    print(f"Saved volume to {filename}")


# ==========================================================
#  DEMO MAIN
# ==========================================================

if __name__ == "__main__":
    # 1) single snapshot + volume export
    nr = 400
    nz = 400
    n_phi = 128

    r = np.linspace(-2.0, 2.0, nr)
    z = np.linspace(-4.0, 4.0, nz)
    t_max = 1.0
    t0 = 0.6
    flux0 = flux_driver(t0, t_max)

    emiss2d, v2d, R_pre, Z = red_spider_slice(r, z, t0, t_max, flux0)
    density3d, vel3d, phi = revolve_to_3d(emiss2d, v2d, r, z, n_phi=n_phi)

    # quick 2D slice plot
    plt.figure(figsize=(6, 10))
    plt.imshow(
        emiss2d,
        origin="lower",
        extent=[r.min(), r.max(), z.min(), z.max()],
        cmap="inferno",
        aspect="auto"
    )
    plt.colorbar(label="Emission")
    plt.title("Red Spider Nebula – 2D slice (snapshot)")
    plt.xlabel("r (symmetric)")
    plt.ylabel("z")
    plt.tight_layout()
    plt.show()

    # 3D ray-march preview
    img = raymarch_render(density3d, vel3d)
    plt.figure(figsize=(6, 10))
    plt.imshow(
        img,
        origin="lower",
        extent=[r.min(), r.max(), z.min(), z.max()],
        cmap="inferno",
        aspect="auto"
    )
    plt.title("Red Spider Nebula – Raymarched 3D Volume")
    plt.xlabel("x (symmetric)")
    plt.ylabel("z")
    plt.tight_layout()
    plt.show()

    # export voxel grid
    export_volume_npz(density3d, vel3d, filename="red_spider_volume.npz")

    # 2) animation (comment out if you only want snapshot)
    # generate_animation(
    #     n_frames=40,
    #     nr=300,
    #     nz=300,
    #     n_phi=96,
    #     save_prefix="red_spider_anim",
    #     make_gif=True
    # )
