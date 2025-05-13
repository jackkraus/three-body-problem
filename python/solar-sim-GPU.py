import cupy as cp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
G = 6.673e-11  # Gravitational constant
mass_sun   = 1.98855e30
mass_earth = 5.972e24
mass_moon  = 7.347e22
num_bodies = 3

# Initialize positions, velocities, and masses (as GPU arrays)
def init():
    r = cp.zeros((num_bodies, 2))  # Positions [x,y]
    v = cp.zeros((num_bodies, 2))  # Velocities [vx, vy]
    m = cp.zeros(num_bodies)       # Masses

    # Earth
    r[0] = cp.array([1.4709e11, 0.0])
    v[0] = cp.array([0.0, 30290.0])
    m[0] = mass_earth

    # Moon
    r[1] = cp.array([1.4709e11, 382500000.0])
    v[1] = cp.array([1022.0, 30290.0])
    m[1] = mass_moon

    # Sun
    r[2] = cp.array([0.0, 0.0])
    v[2] = cp.array([0.0, 0.0])
    m[2] = mass_sun

    return r, v, m

# Compute acceleration due to gravity on the GPU
def acceleration(r, m):
    a = cp.zeros_like(r)
    for i in range(num_bodies):
        for j in range(num_bodies):
            if i != j:
                dist = r[j] - r[i]
                dist_mag = cp.linalg.norm(dist)
                a[i] += G * m[j] * dist / (dist_mag**3 + 1e-8)
    return a

# Runge-Kutta 4th order integration (GPU version)
def rk4_step(r, v, m, h):
    a = acceleration(r, m)
    k1_r = v
    k1_v = a

    k2_r = v + 0.5 * h * k1_v
    k2_v = acceleration(r + 0.5 * h * k1_r, m)

    k3_r = v + 0.5 * h * k2_v
    k3_v = acceleration(r + 0.5 * h * k2_r, m)

    k4_r = v + h * k3_v
    k4_v = acceleration(r + h * k3_r, m)

    r_next = r + (h / 6) * (k1_r + 2*k2_r + 2*k3_r + k4_r)
    v_next = v + (h / 6) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
    return r_next, v_next

# Main simulation loop running on the GPU
def simulate():
    Nt = 1000  # Number of time steps (reduce if needed)
    ht = 130   # Time step size

    r, v, m = init()
    positions = cp.zeros((Nt, num_bodies, 2))
    for n in range(Nt):
        positions[n] = r
        r, v = rk4_step(r, v, m, ht)
    # Transfer positions back to host (CPU) memory for plotting
    return cp.asnumpy(positions)

# Animate the trajectories
def animate_trajectories(positions):
    # Radii of the bodies for plotting
    radii = [6.371e6, 1.737e6, 6.9634e8]  # Earth, Moon, Sun
    scale_factor = 1e-7
    marker_sizes = [r * scale_factor for r in radii]

    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_xlim(-2e11, 2e11)
    ax.set_ylim(-2e11, 2e11)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Three-Body Problem GPU Animation")
    ax.grid()
    ax.axis("equal")

    lines = [ax.plot([], [], label=lbl)[0] for lbl in ["Earth", "Moon", "Sun"]]
    scatters = [ax.scatter([], [], s=marker_sizes[i], zorder=5) for i in range(num_bodies)]

    def init_anim():
        for line, scatter in zip(lines, scatters):
            line.set_data([], [])
            scatter.set_offsets([])
        return lines + scatters

    def update(frame):
        for i, (line, scatter) in enumerate(zip(lines, scatters)):
            line.set_data(positions[:frame, i, 0], positions[:frame, i, 1])
            scatter.set_offsets(positions[frame, i])
        return lines + scatters

    ani = FuncAnimation(fig, update, frames=len(positions), init_func=init_anim, blit=True, interval=20)
    plt.legend()
    plt.show()
    # Optionally save animation:
    # ani.save("three_body_simulation_gpu.mp4", fps=30, writer="ffmpeg")

if __name__ == "__main__":
    pos = simulate()
    animate_trajectories(pos)