# Tutorial

This tutorial walks through a complete metadynamics PMF calculation,
from a prepared starting structure to a corrected PMF plot.

---

## 0. Prepare your simulation directory

Create a working directory for your system and copy in the config:

```bash
mkdir MOL_in_SOLVENT && cd MOL_in_SOLVENT
cp /path/to/metadpmf/config.yaml .
```

You need two input files in this directory:

```
MOL_in_SOLVENT/
├── config.yaml
├── start.gro   # equilibrated starting structure of the dimer in solvent
└── system.top            # GROMACS topology (+ any included .itp files)
```

The starting structure should already be NPT-equilibrated. The two molecules
of interest must be present in the box with their atoms numbered consecutively
so you can identify them in the `molecule` section of `config.yaml`.

---

## 1. Edit config.yaml

Set the atom indices of the beads/atoms that define each molecule's centre,
and the temperature:

```yaml
temperature: 298.15

molecule:
  mol1_atoms: [1, 2, 3]
  mol2_atoms: [4, 5, 6]
```

Tune the PLUMED metadynamics parameters for your system, in particular
`plumed.grid_min/max` (must span the full distance range you want to explore)
and `plumed.wall_at` (upper wall to keep sampling efficient).

Choose your backend. For a workstation:

```yaml
backend: local
backends:
  local:
    gmx:    gmx
    plumed: plumed
```

For a SLURM cluster:

```yaml
backend: slurm
backends:
  slurm:
    gmx:    gmx_mpi
    plumed: plumed
    mdrun:  "srun gmx_mpi mdrun"
    grompp: "srun gmx_mpi grompp"
    header: slurm_header.sh
```

`slurm_header.sh` contains your `#SBATCH` directives and `module load` lines.
This is cluster-specific; consult your HPC documentation.

---

## 2. Generate run files

```bash
metadpmf run
```

This writes:

- `plumed.dat` — PLUMED input for the production metadynamics run, generated
  from your config (atom indices, METAD parameters, upper wall)
- `metad.mdp` — GROMACS MDP with your temperature, timestep, and nsteps
  (or skipped if you supply your own MDP via `paths.mdp`)
- `run_metad.sh` — shell script: grompp followed by mdrun with `--plumed`

---

## 3. Run the metadynamics simulation

```bash
bash run_metad.sh    # local
# or: sbatch run_metad.sh   # slurm
```

This produces `metad.xtc` (compressed trajectory), `HILLS` (deposited
Gaussians), and `COLVAR` (CV and bias values during the run).

Monitor `HILLS` to confirm that Gaussians are being deposited and the bias
is building up. The simulation is converged when the FES no longer changes
significantly with additional sampling.

---

## 4. Reweight the trajectory

```bash
metadpmf reweight
```

This writes:

- `analysis/plumed_analysis.dat` — PLUMED input for `plumed driver`:
  reads `HILLS` to reconstruct the bias on each frame, `PACE=1e9` so no new
  hills are added during the analysis
- `reweight.sh` — copies `HILLS` into `analysis/` and runs `plumed driver`

```bash
bash reweight.sh     # local
# or: sbatch reweight.sh   # slurm
```

Produces `analysis/COLVAR` with three columns: `time  distance  bias`.

---

## 5. Block FES analysis

```bash
metadpmf fes
```

Reads `analysis/COLVAR`, computes per-frame reweighting weights
$w_i = \exp\!\bigl[(V_i - V_\text{max})/k_BT\bigr]$,
then runs block analysis for block sizes 1 to `fes.block_max` (step 10).

Output in `analysis/`:

| File | Contents |
|---|---|
| `dist.weight` | Per-frame (distance, reweighting weight) |
| `errors.block` | (block size, mean FES error) — inspect for convergence |
| `fes.dat` | FES from the largest block size |
| `blocks/fes.{N}.dat` | FES for every block size N |

Inspect `errors.block` to verify the error has plateaued at large block
sizes (sign of a well-converged simulation). If you want to use a specific
block size directly:

```bash
metadpmf fes --block-size 200
```

---

## 6. Apply Jacobian correction and plot

```bash
metadpmf pmf
```

Reads `analysis/fes.dat` and applies a Jacobian (translational entropy)
correction. When the PMF is computed along a radial distance $r$ in 3D, the
spherical volume element scales as $r^2$, meaning there is more accessible
phase-space volume at larger separations even for a flat interaction. The
correction removes this purely geometric bias:

$$\Delta G(r) = \text{FES}(r) + 2RT\ln(r)$$

The corrected profile is then shifted to zero by subtracting the mean value
over the reference plateau region defined by `pmf.shift_range` (default
1.5–1.7 nm). This range should correspond to distances where the two
molecules no longer interact (flat, noise-dominated region of the FES).

Output in `analysis/`:

| File | Contents |
|---|---|
| `fes_corr_and_shifted.dat` | Jacobian-corrected and shifted PMF (nm, kJ/mol, error) |
| `pmf.pdf` | Plot of the shifted PMF with statistical error bars |

!!! tip "Adjusting the zero reference"
    If the default `shift_range` does not fall in a flat region of your FES,
    adjust it in `config.yaml`:
    ```yaml
    pmf:
      shift_range: [1.8, 2.2]
    ```

## Directory layout after a complete run

```
MOL_in_SOLVENT/
├── config.yaml
├── start.gro
├── system.top
├── plumed.dat
├── metad.mdp
├── run_metad.sh
├── metad.tpr / metad.xtc / metad.log
├── HILLS
├── COLVAR
├── reweight.sh
└── analysis/
    ├── plumed_analysis.dat
    ├── HILLS
    ├── COLVAR
    ├── dist.weight
    ├── errors.block
    ├── fes.dat
    ├── fes_corr_and_shifted.dat
    ├── pmf.pdf
    └── blocks/
        ├── fes.1.dat
        ├── fes.11.dat
        └── ...
```
