"""
vicsek_animation.py
===================
Real-time animation of the Vicsek flocking model.

Physics
-------
N self-propelled particles move at constant speed v0 in a 2D periodic box of
side L.  At each time step every particle aligns its direction with the average
direction of all neighbours within radius R, then a uniform random noise
η ∈ [-π, π] is added to the angle.

    θᵢ(t+1) = ⟨θⱼ⟩_{|rᵢ-rⱼ|≤R}  +  η · ξ,    ξ ~ U(-π, π)

Low η  →  ordered flocking.
High η →  disordered, random motion.

Controls
--------
Adjust the parameters in the PARAMETERS block at the bottom of this file.

Usage
-----
    python vicsek_animation.py
"""

import math
import random

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


# ---------------------------------------------------------------------------
# Core simulation functions
# ---------------------------------------------------------------------------

def initialise(N, L):
    """
    Randomly place N particles in a square box of side L.

    Each particle is assigned a random position (x, y) ∈ [0, L)²  and a
    unit velocity vector (vx, vy) pointing in a random direction.

    Parameters
    ----------
    N : int    Number of particles.
    L : float  Box side length.

    Returns
    -------
    x, y   : ndarray, shape (N,)  Positions.
    vx, vy : ndarray, shape (N,)  Unit velocity components.
    """
    x  = np.array([random.random() * L for _ in range(N)])
    y  = np.array([random.random() * L for _ in range(N)])
    th = np.array([random.random() * 2 * np.pi for _ in range(N)])
    vx = np.cos(th)
    vy = np.sin(th)
    return x, y, vx, vy


def move(x, y, vx, vy, v0, L):
    """
    Advance all particle positions by one time step (Δt = 1).

    Applies periodic boundary conditions so particles re-enter on the
    opposite side when they leave the box.

    Parameters
    ----------
    x, y   : ndarray  Positions (modified in place).
    vx, vy : ndarray  Unit velocity components.
    v0     : float    Speed (distance moved per time step).
    L      : float    Box side length.
    """
    x += vx * v0
    y += vy * v0
    # Periodic boundary conditions
    x %= L
    y %= L


def align(x, y, vx, vy, eta, R):
    """
    Update all velocity directions by neighbour-alignment + noise.

    For each particle k, compute the average direction of all particles
    (including k itself) within radius R, then perturb by noise η.

    Parameters
    ----------
    x, y   : ndarray  Positions.
    vx, vy : ndarray  Unit velocity components (updated in place).
    eta    : float    Noise amplitude ∈ [0, 1].  Effective noise = η·π.
    R      : float    Interaction radius.
    """
    N       = len(x)
    new_vx  = np.zeros(N)
    new_vy  = np.zeros(N)

    for k in range(N):
        # Find neighbours within radius R
        dx   = x - x[k]
        dy   = y - y[k]
        mask = dx ** 2 + dy ** 2 <= R ** 2

        # Average direction of neighbours
        mean_vx = vx[mask].sum()
        mean_vy = vy[mask].sum()
        norm    = math.sqrt(mean_vx ** 2 + mean_vy ** 2)

        # Add angular noise
        theta      = math.atan2(mean_vy / norm, mean_vx / norm)
        theta     += eta * np.random.uniform(-np.pi, np.pi)
        new_vx[k]  = math.cos(theta)
        new_vy[k]  = math.sin(theta)

    vx[:] = new_vx
    vy[:] = new_vy


# ---------------------------------------------------------------------------
# Animation
# ---------------------------------------------------------------------------

def run_animation(N, L, eta, R, v0, n_frames):
    """
    Open an interactive window showing the flocking animation.

    Parameters
    ----------
    N        : int    Number of particles.
    L        : float  Box side length.
    eta      : float  Noise amplitude ∈ [0, 1].
    R        : float  Interaction radius.
    v0       : float  Particle speed.
    n_frames : int    Number of animation frames.
    """
    x, y, vx, vy = initialise(N, L)

    fig, axis = plt.subplots(figsize=(6, 6))
    axis.set_xlim(0, L)
    axis.set_ylim(0, L)
    axis.set_aspect('equal')
    axis.set_title(f'Vicsek model   N={N},  η={eta},  R={R:.1f},  v₀={v0:.2f}')
    axis.set_xlabel('x')
    axis.set_ylabel('y')

    colors   = np.random.uniform(-np.pi, np.pi, size=N)
    quiver   = axis.quiver(x, y, vx, vy, colors, cmap='hsv', clim=(-np.pi, np.pi))

    def _step(frame):
        move(x, y, vx, vy, v0, L)
        align(x, y, vx, vy, eta, R)
        quiver.set_offsets(np.column_stack((x, y)))
        quiver.set_UVC(vx, vy)
        return quiver,

    anim = animation.FuncAnimation(
        fig, _step, frames=np.arange(n_frames), interval=1, blit=True
    )
    plt.tight_layout()
    plt.show()
    return anim   # keep reference so animation is not garbage-collected


# ---------------------------------------------------------------------------
# Entry point — adjust parameters here
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    # ---- Simulation parameters ----
    N        = 300     # number of particles
    L        = 100.0   # box side length
    eta      = 0.2     # noise amplitude  (0 = no noise, 1 = maximum noise)
    R        = 0.10*L  # interaction radius
    v0       = 0.01*L  # particle speed
    n_frames = 200     # number of animation frames

    run_animation(N, L, eta, R, v0, n_frames)

    