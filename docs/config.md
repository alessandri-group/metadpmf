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
forcefield: martini   # 'martini' or 'opls' — REQUIRED, no default
```

Required — there is no default, so you must set it to `martini` or `opls`
(any other value is rejected). Selects the built-in production MDP template
written by `metadpmf run`:

| Value | Template | Timestep | Electrostatics | Constraints |
|---|---|---|---|---|
| `martini` | `md_martini.mdp` | 20 fs | reaction-field | none |
| `opls` | `md_opls.mdp` | 2 fs | PME | h-bonds |

The `mdp.dt` default follows this choice automatically (0.020 for Martini,
0.002 for OPLS). Override it under `mdp:` if you need a different timestep.
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

---

## `mdp`

```yaml
mdp:
  nsteps: 50000000
  # dt:   0.020   # default follows forcefield (0.020 martini / 0.002 opls)
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
