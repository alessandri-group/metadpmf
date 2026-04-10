# Configuration reference

All settings live in a single `config.yaml` file that you copy into your
simulation directory and edit.

---

## `temperature`

```yaml
temperature: 300.0
```

Simulation temperature in Kelvin. Used in:

- The PLUMED `TEMP=` parameter
- The MDP `ref_t` and `gen_temp` fields
- The reweighting factor $k_BT$
- The Jacobian correction $2RT\ln(r)$

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
  dt:     0.020
```

Applied to the built-in Martini CG MDP template. If you supply your own MDP
via `paths.mdp`, these keys are ignored.

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
| `mdp` | Custom MDP file; `null` = use built-in Martini CG template |
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
