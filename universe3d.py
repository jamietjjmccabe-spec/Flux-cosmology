import os
os.environ['VISPY_GL_BACKEND'] = 'gl2'

import vispy
vispy.use('pyqt5')

import numpy as np
from vispy import app, scene

from flux_background_api import background
from structure_seeds import generate_bh_seeds
from flux_field3d import build_flux_field


# ============================================
# LOAD BACKGROUND COSMOLOGY
# ============================================
bg = background()

t_arr   = bg["t"]
a_arr   = bg["a"]
H_arr   = bg["H"]
G_arr   = bg["G_eff"]
Omega_s = bg["Omega_sigma"]

N_steps = len(t_arr)
current = int(N_steps * 0.8)     # start late so a ~ 1


# ============================================
# SIMULATION PARAMETERS
# ============================================
N_PARTICLES = 200000
N_BH        = 40
BOX         = 1.0
SOFT        = 0.05


# ============================================
# BUILD FLUX FIELD
# ============================================
flux_field = build_flux_field(BOX, grid_n=64, amp=1.5, rng=123)
GRID_N     = flux_field["N"]
sigma_grid = flux_field["sigma"]# reusable colour buffer
colors = np.zeros((N_PARTICLES, 4), dtype=np.float32)
grad_x     = flux_field["grad_x"]
grad_y     = flux_field["grad_y"]
grad_z     = flux_field["grad_z"]
colors = np.zeros((N_PARTICLES, 4), dtype=np.float32)


def sample_sigma(xc, yc, zc):
    """Sample σ field at comoving coords."""
    N = GRID_N
    L = 2.0 * BOX

    ix = ((xc + BOX) / L * (N - 1)).astype(np.int32)
    iy = ((yc + BOX) / L * (N - 1)).astype(np.int32)
    iz = ((zc + BOX) / L * (N - 1)).astype(np.int32)

    ix = np.clip(ix, 0, N - 1)
    iy = np.clip(iy, 0, N - 1)
    iz = np.clip(iz, 0, N - 1)

    return sigma_grid[ix, iy, iz]


# ============================================
# INITIAL PARTICLES
# ============================================
rng = np.random.default_rng(77)

x_c = rng.uniform(-BOX, BOX, size=N_PARTICLES)
y_c = rng.uniform(-BOX, BOX, size=N_PARTICLES)
z_c = rng.uniform(-BOX, BOX, size=N_PARTICLES)

vx = rng.normal(0, 0.02, size=N_PARTICLES)
vy = rng.normal(0, 0.02, size=N_PARTICLES)
vz = rng.normal(0, 0.02, size=N_PARTICLES)


# ============================================
# BLACK HOLE SEEDS
# ============================================
bh_pos, bh_mass = generate_bh_seeds(
    N_BH, BOX, M_min=5.0, M_max=80.0, rng=42
)


# ============================================
# FORCE FUNCTIONS
# ============================================
def bh_forces(xc, yc, zc, a_now, G_eff):
    """BH gravity in comoving coords."""
    x = a_now * xc
    y = a_now * yc
    z = a_now * zc

    ax = np.zeros_like(x)
    ay = np.zeros_like(y)
    az = np.zeros_like(z)

    for (bx_c, by_c, bz_c), M in zip(bh_pos, bh_mass):
        bx = a_now * bx_c
        by = a_now * by_c
        bz = a_now * bz_c

        dx = x - bx
        dy = y - by
        dz = z - bz

        r2 = dx*dx + dy*dy + dz*dz + SOFT**2
        inv_r3 = 1.0 / (r2 * np.sqrt(r2))

        coeff = -G_eff * M * inv_r3
        ax += coeff * dx
        ay += coeff * dy
        az += coeff * dz

    # convert back to comoving accelerations
    return ax / a_now, ay / a_now, az / a_now


def flux_forces_sigma(xc, yc, zc, strength=3.0):
    """Flux-gradient forces from ∇σ."""
    N = GRID_N
    L = 2.0 * BOX

    ix = ((xc + BOX) / L * (N - 1)).astype(np.int32)
    iy = ((yc + BOX) / L * (N - 1)).astype(np.int32)
    iz = ((zc + BOX) / L * (N - 1)).astype(np.int32)

    ix = np.clip(ix, 0, N - 1)
    iy = np.clip(iy, 0, N - 1)
    iz = np.clip(iz, 0, N - 1)

    gx = grad_x[ix, iy, iz]
    gy = grad_y[ix, iy, iz]
    gz = grad_z[ix, iy, iz]

    return -strength * gx, -strength * gy, -strength * gz


# ============================================
# VISUAL SETUP
# ============================================
canvas = scene.SceneCanvas(keys='interactive',
                           bgcolor='black',
                           size=(1600, 900))
view   = canvas.central_widget.add_view()
view.camera = scene.cameras.TurntableCamera(fov=55)
view.camera.distance = 2.0

a0     = a_arr[current]
pos0   = np.vstack((a0*x_c, a0*y_c, a0*z_c)).T.astype(np.float32)

scatter = scene.visuals.Markers(parent=view.scene)
scatter.set_data(pos0, size=2.5, face_color=[0.3, 0.7, 1.0, 0.9])

bh_pos_phys0 = a0 * bh_pos

# bright white core
bh_core = scene.visuals.Markers(parent=view.scene)
bh_core.set_data(bh_pos_phys0.astype(np.float32),
                 size=12,
                 face_color=[1.0, 1.0, 1.0, 1.0])

# big soft golden halo
bh_halo = scene.visuals.Markers(parent=view.scene)
bh_halo.set_data(bh_pos_phys0.astype(np.float32),
                 size=28,
                 face_color=[1.0, 0.85, 0.2, 0.18])  # RGBA (low alpha)


text = scene.Text("", color='white',
                  parent=view.scene, pos=(20, 40))


# ============================================
# UPDATE LOOP
# ============================================
def update(ev):
    global current, x_c, y_c, z_c, vx, vy, vz, colors

    if current >= N_steps - 1:
        return

    t0 = t_arr[current]
    t1 = t_arr[current+1]
    dt = t1 - t0
    a_now = a_arr[current]
    G_eff = G_arr[current]

    # Total acceleration = BH gravity + flux-gradient
    ax_g, ay_g, az_g = bh_forces(x_c, y_c, z_c, a_now, G_eff)
    ax_s, ay_s, az_s = flux_forces_sigma(x_c, y_c, z_c, strength=3.0)

    ax = ax_g + ax_s
    ay = ay_g + ay_s
    az = az_g + az_s

    vx[:] += ax * dt
    vy[:] += ay * dt
    vz[:] += az * dt

    x_c[:] += vx * dt
    y_c[:] += vy * dt
    z_c[:] += vz * dt

    # periodic box
    x_c[:] = (x_c + BOX) % (2*BOX) - BOX
    y_c[:] = (y_c + BOX) % (2*BOX) - BOX
    z_c[:] = (z_c + BOX) % (2*BOX) - BOX

    # ---------- σ-based colouring ----------
    sigma_now = sample_sigma(x_c, y_c, z_c)
    sig_min = sigma_now.min()
    sig_max = sigma_now.max()
    sigma_norm = (sigma_now - sig_min) / (sig_max - sig_min + 1e-9)

    # low σ (wells / filaments) → bright, high σ (voids) → dark
    lum = 1.0 - sigma_norm

    colors[:,0] = lum**2         # red
    colors[:,1] = lum**0.5       # green
    colors[:,2] = 1.0 - lum      # blue
    colors[:, 3] = 1.0      # alpha
    # ---------------------------------------

    pos = np.vstack((a_now*x_c, a_now*y_c, a_now*z_c)).T.astype(np.float32)
    scatter.set_data(pos, size=2.8, face_color=colors)

    # update BH markers (core + halo)
    bh_pos_phys = a_now * bh_pos

    bh_core.set_data(bh_pos_phys.astype(np.float32),
                     size=12,
                     face_color=[1.0, 1.0, 1.0, 1.0])

    bh_halo.set_data(bh_pos_phys.astype(np.float32),
                     size=28,
                     face_color=[1.0, 0.85, 0.2, 0.18])

    # status text
    text.text = f"step {current}/{N_steps}  a={a_now:.4e}  G_eff={G_eff:.6f}"

    current += 1



timer = app.Timer(interval=0.016, connect=update, start=True)
canvas.show()
app.run()
