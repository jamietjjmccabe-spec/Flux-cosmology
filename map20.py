import numpy as np
import matplotlib.pyplot as plt

# Three-body ejection toy model (in A-layer), interpreted as 3 Φ-fragments

G = 1.0  # gravitational constant in toy units

# Masses: two heavy fragments + one lighter fragment
m = np.array([1.0, 1.0, 0.3])

# Initial positions (x,y)
r0 = np.array([
    [-1.0, 0.0],   # fragment A
    [ 1.0, 0.0],   # fragment B
    [ 0.0, 0.5],   # fragment C (lighter, slightly above)
])

# Initial velocities (x,y)
# Give A and B small opposite y velocities, C small x velocity to stir chaos
v0 = np.array([
    [0.2,  0.1],   # A
    [-0.2, 0.1],   # B
    [0.0, -0.4],   # C
])

def accelerations(r):
    """Compute gravitational accelerations on each body."""
    N = len(m)
    a = np.zeros_like(r)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            diff = r[j] - r[i]
            dist2 = np.dot(diff, diff) + 1e-6  # softening
            dist3 = dist2 * np.sqrt(dist2)
            a[i] += G * m[j] * diff / dist3
    return a

# Time integration (leapfrog)
dt = 0.01
steps = 8000

r = r0.copy()
v = v0.copy()

traj = np.zeros((steps, 3, 2))

# Initial half-step velocity
a = accelerations(r)
v += 0.5 * dt * a

for n in range(steps):
    # Drift
    r += dt * v
    # Kick
    a = accelerations(r)
    v += dt * a
    traj[n] = r

# Extract trajectories
xA, yA = traj[:, 0, 0], traj[:, 0, 1]
xB, yB = traj[:, 1, 0], traj[:, 1, 1]
xC, yC = traj[:, 2, 0], traj[:, 2, 1]

plt.figure(figsize=(6,6))
plt.plot(xA, yA, label="Fragment A (m=1.0)")
plt.plot(xB, yB, label="Fragment B (m=1.0)")
plt.plot(xC, yC, label="Fragment C (m=0.3)")
plt.scatter([xA[0], xB[0], xC[0]], [yA[0], yB[0], yC[0]], marker='o')  # starting points
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.title("Three-body interaction: lighter fragment ejection")
plt.legend()
plt.tight_layout()
plt.show()
