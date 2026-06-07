# Configuration reference

All settings live in a single `config.yaml` file that you copy into your
simulation directory and edit.

---

## `temperature`

```yaml
temperature: 298.15
```

Simulation temperature in Kelvin. Used in:

- The PLUMED `TEMP=` parameter
- The MDP `ref_t` and `gen_temp` fields
- The reweighting factor $k_BT$
- The Jacobian correction $2RT\ln(r)$

---

## `forcefield`

```yaml
forcefield: martini   # 'martini', 'opls', or 'gaff2' — REQUIRED, no default
```

Required — there is no default, so you must set it to `martini`, `opls`, or
`gaff2` (any other value is rejected). Selects the built-in production MDP
template written by `metadpmf run`:

| Value | Template | Timestep | Electrostatics | Constraints |
|---|---|---|---|---|
| `martini` | `md_martini.mdp` | 20 fs | reaction-field | none |
| `opls` | `md_opls.mdp` | 2 fs | PME | h-bonds |
| `gaff2` | `md_gaff2.mdp` | 2 fs | PME | h-bonds |

The `mdp.dt` default follows this choice automatically (0.020 for Martini,
0.002 for OPLS/GAFF2). Override it under `mdp:` if you need a different timestep.
Ignored if you supply your own MDP via `paths.mdp`.

---

## `molecule`

```yaml
molecule:
  mol1_atoms: [1, 2, 3]
  mol2_atoms: [4, 5, 6]
```

1-based atom indices (GROMACS/PLUMED convention) of the beads/atoms that
define the centre of each molecule. These are used in:

- `WHOLEMOLECULES ENTITY0=... ENTITY1=...` — ensures each molecule is kept
  whole across periodic boundaries
- `CENTER ATOMS=...` — defines the collective variable as the distance
  between the geometric centres

For a Martini benzene dimer, all three ring beads of each molecule are
listed.

---

## `plumed`

```yaml
plumed:
  pace:         500
  height:       1.0
  sigma:        0.05
  biasfactor:   5
  grid_min:     0.1
  grid_max:     3.0
  grid_bin:     300
  wall_at:      2.0
  wall_kappa:   200.0
  print_stride: 5000
```

| Key | Description |
|---|---|
| `pace` | Deposit a Gaussian every N MD steps |
| `height` | Initial Gaussian height (kJ/mol) |
| `sigma` | Gaussian width (nm) |
| `biasfactor` | Well-tempered scaling factor ($\gamma$); higher = more exploration |
| `grid_min/max` | CV grid bounds (nm); should span the full distance range of interest |
| `grid_bin` | Number of grid bins |
| `wall_at` | Upper wall position (nm); keeps sampling efficient |
| `wall_kappa` | Upper wall force constant (kJ/mol/nm²) |
| `print_stride` | Write COLVAR every N steps |

### How the grid, wall, and FES ranges relate

Three distance ranges appear in the config, and they must be nested correctly:

- **METAD grid** (`grid_min`, `grid_max`, `grid_bin`) — PLUMED stores the
  accumulated bias on this grid instead of re-summing every hill each step. The
  CV must *never* leave the grid: if the distance exceeds `grid_max`, PLUMED
  errors out. The grid is therefore the **widest** of the ranges.
- **Upper wall** (`wall_at`, `wall_kappa`) — a harmonic restraint that pushes
  the molecules back once they pass `wall_at`. It is needed because at large
  separation the PMF is flat: nothing pulls the molecules back, and
  metadynamics would otherwise keep filling the basin and drive them apart
  indefinitely (and off the grid). The wall caps the explorable range.
- **FES histogram** (`fes.min`, `fes.max`) — the range you actually analyse. It
  must stay *inside* the wall, because near and beyond the wall the sampling is
  distorted by the wall potential and the PMF there is contaminated.

The required ordering is:

```
grid_min ≤ fes.min  <  fes.max  ≤  wall_at  <  grid_max
   0.1      0.2          1.7        2.0        3.0
```

`metadpmf` validates `fes.max ≤ wall_at` and that `pmf.shift_range` lies within
`[fes.min, fes.max]`.

**Choosing `wall_at` for your box.** Under the minimum-image convention a
centre-of-mass distance is only meaningful up to **half the shortest box edge**
(`L/2`); beyond that a molecule starts seeing periodic images of its partner.
Place the wall *below* `L/2` with margin, because both the soft-wall
fluctuations (the CV strays ~0.2–0.3 nm above `wall_at` before being pushed
back) and the grid headroom (`grid_max > wall_at`) must still fit under `L/2`.
For a cubic box of edge 5.15 nm (`L/2 ≈ 2.58`), a wall around 2.2 nm with
`grid_max ≈ 2.5` is a safe choice. Below `wall_at` the wall potential is
exactly zero, so the FES is energetically clean up to the wall — but expect
poorer statistics (larger error bars) in the ~0.1 nm just below it as sampling
density drops.

---

## `mdp`

```yaml
mdp:
  nsteps: 50000000
  # dt:   0.020   # default follows forcefield (0.020 martini / 0.002 opls,gaff2)
```

Applied to the built-in MDP template selected by `forcefield`. `dt` defaults
to the forcefield-appropriate value but can be overridden here. If you supply
your own MDP via `paths.mdp`, these keys are ignored.

`nsteps` is 50M by default: 1 µs for Martini (20 fs) or 100 ns for OPLS (2 fs).

---

## `fes`

```yaml
fes:
  min:       0.2
  max:       1.7
  bins:      51
  block_max: 1000
```

Parameters for the block FES analysis (`metadpmf fes`).

| Key | Description |
|---|---|
| `min` / `max` | Histogram range (nm); should match the sampled CV range |
| `bins` | Number of histogram bins |
| `block_max` | Scan block sizes `range(1, block_max, 10)`; inspect `errors.block` to verify convergence |

---

## `pmf`

```yaml
pmf:
  shift_range: [1.5, 1.7]
  xrange: null
  yrange: null
```

| Key | Description |
|---|---|
| `shift_range` | Distance range (nm) used to define the PMF zero. Should be a plateau region where the molecules no longer interact. |
| `xrange` | x-axis limits for the plot in nm, e.g. `[0.3, 2.0]`; `null` = auto |
| `yrange` | y-axis limits in kJ/mol, e.g. `[-10, 5]`; `null` = auto |

---

## `cv2`

```yaml
cv2:
  enabled: false
  type:    costheta
  bias:    false
  sigma:   0.1
  fold:    none
  label:   "cos θ"
  min:    -1.0
  max:     1.0
  bins:    51
  # mol1_plane: [1, 2, 3]
  # mol2_plane: [4, 5, 6]
  plumed_analysis: plumed_analysis_2d.dat
```

A second collective variable, used to build a 2D free-energy surface. The
metadynamics bias **always** acts on the centre-of-mass distance; `cv2` adds an
orientational coordinate that is either projected onto or also biased.

### Three modes

`enabled` and `bias` together select the mode:

| Mode | `enabled` | `bias` | Bias acts on | Output |
|---|---|---|---|---|
| **A** | `false` | — | distance | 1D PMF |
| **B** — project | `true` | `false` | distance | 2D FES (cv2 unbiased) |
| **C** — bias both | `true` | `true` | distance + cv2 | 2D FES |

In Mode B the production run is identical to Mode A; cv2 is computed only during
reweighting. In Mode C cv2 is biased by METAD, so it appears in both the
production and analysis PLUMED inputs.

### Keys

| Key | Description |
|---|---|
| `enabled` | Turn the second CV on |
| `type` | `costheta` = built-in generator (cosine of the angle between the two ring-plane normals, built from `mol1_plane`/`mol2_plane`); `custom` = supply your own PLUMED via `plumed_analysis` |
| `bias` | `false` = project (Mode B); `true` = also bias cv2 (Mode C) |
| `sigma` | metadynamics Gaussian width on cv2 (only used when `bias: true`) |
| `fold` | `none`, `abs` (→ \|cosθ\|), or `sq` (→ cos²θ); see below |
| `label` | y-axis label in the 2D plot |
| `min` / `max` | cv2 histogram range; also the metad grid bounds on cv2 when biased |
| `bins` | number of cv2 histogram (and grid) bins |
| `mol1_plane` / `mol2_plane` | 3 atoms per molecule defining each ring plane (for `type: costheta`). Default to `molecule.mol*_atoms` when those have exactly 3 atoms |
| `plumed_analysis` | for `type: custom` only: your PLUMED file that defines cv2 and prints `dist, cv2, metad.bias` |

### cosθ sign degeneracy and `fold`

For two stacked rings, `cosθ = +1` and `cosθ = -1` both mean "planes parallel":
the sign just reflects which face of each ring the normal points from. For a
**planar** molecule the two faces are equivalent (the molecular plane is a mirror,
true even for in-plane-asymmetric molecules like phenol), so the two halves of
the cv2 axis are physically redundant.

- `fold: none` (default) keeps the full `[-1, 1]` range. For projection (Mode B)
  this is useful: a converged FES should be symmetric about 0, a free
  convergence check.
- `fold: abs` or `sq` collapses the redundant half. This is mainly worth it when
  **biasing** a planar molecule (Mode C), where it avoids filling two equivalent
  regions. Set `min: 0.0` to match the folded range. Do **not** fold a molecule
  whose two faces differ (non-planar, or different groups above/below the ring) —
  there `+1` and `-1` are genuinely distinct states.

---

## `paths`

```yaml
paths:
  gro:  start.gro
  top:  system.top
  mdp:  null
  traj: metad.xtc
```

All paths are relative to the config file location.

| Key | Description |
|---|---|
| `gro` | Equilibrated starting structure |
| `top` | GROMACS topology |
| `mdp` | Custom MDP file; `null` = use built-in template selected by `forcefield` |
| `traj` | Trajectory filename produced by mdrun (default `metad.xtc`) |

---

## `backend` / `backends`

```yaml
backend: local

backends:
  local:
    gmx:    gmx
    plumed: plumed

  slurm:
    gmx:    gmx_mpi
    plumed: plumed
    mdrun:  "srun gmx_mpi mdrun"
    grompp: "srun gmx_mpi grompp"
    header: slurm_header.sh
```

Set `backend: slurm` on HPC clusters. The `header` key must point to a file
containing your `#SBATCH` directives and `module load` commands.
