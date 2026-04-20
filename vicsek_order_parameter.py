"""
vicsek_order_parameter.py
=========================
Order parameter dynamics and noise-driven phase transition in the Vicsek model.

Physics
-------
The order parameter φ measures the degree of collective alignment:

    φ(t) = (1/N) |Σᵢ vᵢ(t)|

φ ≈ 1  →  all particles move in the same direction (ordered / flocking phase).
φ ≈ 0  →  random, disordered motion.

This script computes two things:

1. φ(t) as a function of time for a single noise value η  (live animated plot).
2. Time-averaged ⟨φ⟩ as a function of η, sweeping from η=0 to η=1, revealing
   the flocking phase transition.

Results are saved to results/.

Usage
-----
    python vicsek_order_parameter.py
"""

import math
import os
import random

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


# ---------------------------------------------------------------------------
# Core simulation functions  (shared with vicsek_animation.py)
# ---------------------------------------------------------------------------

def initialise(N, L):
    """
    Randomly place N particles in a square box of side L.

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
    Advance all positions by one time step, with periodic boundary conditions.

    Parameters
    ----------
    x, y   : ndarray  Positions (modified in place).
    vx, vy : ndarray  Unit velocity components.
    v0     : float    Speed.
    L      : float    Box side length.
    """
    x += vx * v0
    y += vy * v0
    x %= L
    y %= L


def align(x, y, vx, vy, eta, R):
    """
    Update velocity directions by neighbour-alignment + noise.

    Parameters
    ----------
    x, y   : ndarray  Positions.
    vx, vy : ndarray  Unit velocity components (updated in place).
    eta    : float    Noise amplitude ∈ [0, 1].
    R      : float    Interaction radius.
    """
    N      = len(x)
    new_vx = np.zeros(N)
    new_vy = np.zeros(N)

    for k in range(N):
        dx   = x - x[k]
        dy   = y - y[k]
        mask = dx ** 2 + dy ** 2 <= R ** 2

        mean_vx = vx[mask].sum()
        mean_vy = vy[mask].sum()
        norm    = math.sqrt(mean_vx ** 2 + mean_vy ** 2)

        theta      = math.atan2(mean_vy / norm, mean_vx / norm)
        theta     += eta * np.random.uniform(-np.pi, np.pi)
        new_vx[k]  = math.cos(theta)
        new_vy[k]  = math.sin(theta)

    vx[:] = new_vx
    vy[:] = new_vy


def order_parameter(vx, vy):
    """
    Compute the scalar order parameter φ = (1/N)|Σ vᵢ|.

    Parameters
    ----------
    vx, vy : ndarray  Unit velocity components of all particles.

    Returns
    -------
    float  φ ∈ [0, 1].
    """
    N = len(vx)
    return math.sqrt(vx.sum() ** 2 + vy.sum() ** 2) / N


# ---------------------------------------------------------------------------
# Part 1 — Animated order parameter vs time for a single η
# ---------------------------------------------------------------------------

def animate_order_parameter(N, L, eta, R, v0,
                             t_total, t_avg_start,
                             results_dir):
    """
    Show a live plot of φ(t) and print ⟨φ⟩ averaged over [t_avg_start, t_total].
    Saves the final plot to results/.

    Parameters
    ----------
    N           : int    Number of particles.
    L           : float  Box side length.
    eta         : float  Noise amplitude.
    R           : float  Interaction radius.
    v0          : float  Particle speed.
    t_total     : int    Total number of time steps.
    t_avg_start : int    Time step from which averaging begins.
    results_dir : str    Directory where the figure is saved.
    """
    x, y, vx, vy = initialise(N, L)

    # Accumulators for time-averaged φ
    phi_sum   = 0.0
    avg_count = 0

    time_data = []
    phi_data  = []

    plt.rc('font', family='serif', size=11)
    fig, axis = plt.subplots(figsize=(8, 4))
    axis.set_xlim(0, t_total)
    axis.set_ylim(-0.05, 1.15)
    axis.set_xlabel('Time  (steps)', fontsize=13)
    axis.set_ylabel('Order parameter  φ(t)', fontsize=13)
    fig.suptitle(rf'Vicsek model — order parameter   η = {eta},  N = {N}', fontsize=13)
    line, = axis.plot([], [], 'steelblue', linewidth=1.2)
    # Horizontal line marking the averaging window start
    axis.axvline(t_avg_start, color='gray', linewidth=0.8,
                 linestyle='--', label=f'Averaging starts at t={t_avg_start}')
    axis.legend(fontsize=9)

    def _step(t):
        nonlocal phi_sum, avg_count

        move(x, y, vx, vy, v0, L)
        align(x, y, vx, vy, eta, R)

        phi = order_parameter(vx, vy)
        time_data.append(t)
        phi_data.append(phi)
        line.set_xdata(time_data)
        line.set_ydata(phi_data)

        # Accumulate average after transient
        if t >= t_avg_start:
            phi_sum   += phi
            avg_count += 1

        return line,

    anim = animation.FuncAnimation(
        fig, _step, frames=range(t_total), interval=1, blit=True, repeat=False
    )
    plt.tight_layout()
    plt.show()

    # Save figure after animation completes
    fname = os.path.join(results_dir, f'order_parameter_eta{eta:.2f}.png')
    fig.savefig(fname, dpi=150)
    print(f'  Saved: results/order_parameter_eta{eta:.2f}.png')

    avg_phi = phi_sum / avg_count if avg_count > 0 else float('nan')
    print(f'  ⟨φ⟩ (averaged over t = {t_avg_start} – {t_total}) = {avg_phi:.4f}')

    return anim   # keep reference


# ---------------------------------------------------------------------------
# Part 2 — Noise sweep: ⟨φ⟩ vs η
# ---------------------------------------------------------------------------

def noise_sweep(N, L, eta_values, R, v0,
                t_total, t_avg_start,
                results_dir):
    """
    Run the simulation for each η in eta_values, compute ⟨φ⟩, and plot
    the phase transition curve.

    Parameters
    ----------
    N           : int       Number of particles.
    L           : float     Box side length.
    eta_values  : array     Noise amplitudes to sweep over.
    R           : float     Interaction radius.
    v0          : float     Particle speed.
    t_total     : int       Time steps per run.
    t_avg_start : int       Step from which averaging begins.
    results_dir : str       Directory where the figure is saved.

    Returns
    -------
    avg_phi : list  Time-averaged order parameter for each η.
    """
    avg_phi = []

    print(f'\nNoise sweep:  {len(eta_values)} values,  N={N},  '
          f't_total={t_total},  t_avg from {t_avg_start}')
    print(f'  {"η":>6}   {"⟨φ⟩":>8}')
    print('  ' + '-' * 18)

    for eta in eta_values:
        x, y, vx, vy = initialise(N, L)
        phi_sum   = 0.0
        avg_count = 0

        for t in range(t_total):
            move(x, y, vx, vy, v0, L)
            align(x, y, vx, vy, eta, R)

            if t >= t_avg_start:
                phi_sum   += order_parameter(vx, vy)
                avg_count += 1

        mean_phi = phi_sum / avg_count if avg_count > 0 else float('nan')
        avg_phi.append(mean_phi)
        print(f'  {eta:>6.3f}   {mean_phi:>8.4f}')

    # ---- Plot ⟨φ⟩ vs η ----
    plt.rc('font', family='serif', size=11)
    fig, axis = plt.subplots(figsize=(6, 4))
    axis.plot(eta_values, avg_phi, 'o-', color='steelblue',
              markersize=5, linewidth=1.4)
    axis.set_xlabel(r'Noise amplitude  $\eta$', fontsize=13)
    axis.set_ylabel(r'Time-averaged order parameter  $\langle\varphi\rangle$', fontsize=13)
    axis.set_title(f'Vicsek model — flocking phase transition\n'
                   f'N = {N},  R = {R:.1f},  v₀ = {v0:.3f}', fontsize=11)
    axis.set_xlim(0, max(eta_values) * 1.05)
    axis.set_ylim(-0.05, 1.10)
    axis.axhline(0, color='k', linewidth=0.5, linestyle=':')
    plt.tight_layout()

    fname = os.path.join(results_dir, 'noise_vs_order_parameter.png')
    fig.savefig(fname, dpi=150)
    print(f'\n  Saved: results/noise_vs_order_parameter.png')
    plt.show()
    plt.close(fig)

    return avg_phi


# ---------------------------------------------------------------------------
# Entry point — adjust parameters here
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    # ---- Simulation parameters ----
    N           = 200      # number of particles
    L           = 100.0    # box side length
    eta         = 0.0      # noise for the animated single run
    R           = 0.10*L   # interaction radius
    v0          = 0.02*L   # particle speed
    t_total     = 1200     # total simulation steps
    t_avg_start = 400      # discard transient; average from this step onward

    # ---- Noise sweep parameters ----
    eta_sweep   = np.linspace(0.0, 1.0, 20)   # η values for the phase diagram

    # ---- Output directory ----
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(results_dir, exist_ok=True)

    # ---- Part 1: animated φ(t) for chosen η ----
    print(f'--- Part 1: animated order parameter   η = {eta} ---')
    anim = animate_order_parameter(
        N, L, eta, R, v0, t_total, t_avg_start, results_dir
    )

    # ---- Part 2: noise sweep ⟨φ⟩ vs η ----
    print('\n--- Part 2: noise sweep ---')
    avg_phi = noise_sweep(
        N, L, eta_sweep, R, v0, t_total, t_avg_start, results_dir
    )

    print(f'\nAll figures saved in: {results_dir}')