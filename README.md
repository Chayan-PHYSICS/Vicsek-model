# Vicsek-Model

Simulation and analysis of the **Vicsek flocking model** — a minimal model of collective motion in active matter.

Includes a real-time animation of emergent flocking and a quantitative study of the flocking phase transition via the order parameter.

---

## Physics

The Vicsek model describes N self-propelled particles moving at constant speed v₀ in a 2D periodic box of side L. At each discrete time step every particle:

1. **Moves** in its current direction by distance v₀.
2. **Aligns** with the average direction of all neighbours within radius R.
3. **Perturbs** its angle by a uniform random noise ξ ∈ [−π, π] scaled by η.

```
θᵢ(t+1) = ⟨θⱼ⟩_{|rᵢ−rⱼ|≤R}  +  η · ξ
xᵢ(t+1) = xᵢ(t) + v₀ cos θᵢ(t+1)
yᵢ(t+1) = yᵢ(t) + v₀ sin θᵢ(t+1)
```

The **order parameter** φ measures global alignment:

```
φ(t) = (1/N) |Σᵢ vᵢ(t)|     ∈ [0, 1]
```

- φ ≈ 1 → all particles move together (ordered / flocking phase)
- φ ≈ 0 → disordered, random motion

Increasing the noise η drives a phase transition from the ordered to the disordered phase, analogous to a ferromagnetic transition.

---

## Repository structure

```
Vicsek-Model/
├── vicsek_animation.py         # Real-time flocking animation
├── vicsek_order_parameter.py   # Order parameter dynamics + noise sweep
├── requirements.txt
├── results/                    # Auto-created on first run
│   ├── order_parameter_eta0.00.png
│   └── noise_vs_order_parameter.png
└── README.md
```

---

## Installation

```bash
git clone https://github.com/Chayan-PHYSICS/Vicsek-Model.git
cd Vicsek-Model
pip install -r requirements.txt
```

---

## Usage

### 1 — Flocking animation

```bash
python vicsek_animation.py
```

Opens a real-time animated quiver plot of all N particles coloured by their direction. Adjust parameters at the bottom of the file:

| Parameter | Default | Description                   |
|-----------|---------|-------------------------------|
| `N`       | 300     | Number of particles            |
| `L`       | 100     | Box side length                |
| `eta`     | 0.2     | Noise amplitude ∈ [0, 1]      |
| `R`       | 0.1·L   | Interaction radius             |
| `v0`      | 0.01·L  | Particle speed                 |
| `n_frames`| 200     | Animation frames               |

### 2 — Order parameter and phase transition

```bash
python vicsek_order_parameter.py
```

Runs in two stages:

**Part 1** — animated live plot of φ(t) for a single chosen η. The averaging window [t_avg_start, t_total] is marked on the plot. Saves `results/order_parameter_eta{η}.png`.

**Part 2** — sweeps η from 0 to 1 (configurable), computes ⟨φ⟩ for each value, and plots the phase transition curve. Saves `results/noise_vs_order_parameter.png`.

Example terminal output:

```
--- Part 2: noise sweep ---
Noise sweep:  20 values,  N=200,  t_total=1200,  t_avg from 400
      η      ⟨φ⟩
  ------------------
   0.000    0.9981
   0.053    0.9864
   ...
   0.895    0.1203
   1.000    0.0891

  Saved: results/noise_vs_order_parameter.png
```

---

## Key results

| η (noise) | Phase       | ⟨φ⟩ |
|-----------|-------------|------|
| Low       | Ordered     | ≈ 1  |
| Critical  | Transition  | ~0.5 |
| High      | Disordered  | ≈ 0  |

The transition is sharp at intermediate densities (ρ = N/L²), consistent with the original Vicsek et al. (1995) result.

---

## Implementation notes

- **Periodic boundary conditions** are implemented as `x %= L`, ensuring particles wrap correctly in both directions.
- The order parameter is averaged over the **steady-state window** `[t_avg_start, t_total]` to exclude the initial transient.
- Vectorised NumPy operations are used for neighbour detection (`dx² + dy² ≤ R²`), replacing explicit Python loops in the original.

---

## Reference

Vicsek, T., Czirók, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995).
*Novel type of phase transition in a system of self-driven particles.*
Physical Review Letters, 75(6), 1226.

---

## Related work

See also the [Numerov-QHO](https://github.com/Chayan-PHYSICS/Numerov-QHO) repository for a numerical quantum mechanics project, and the [TEBD-spin-chain](https://github.com/Chayan-PHYSICS/TEBD-spin-chain) repository for many-body quantum dynamics.
