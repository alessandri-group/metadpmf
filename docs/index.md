<p align="center">
  <img src="assets/logo.png" width="280" alt="metadpmf logo"/>
</p>

# metadpmf

**metadpmf** is a Python package for running and analysing well-tempered
metadynamics simulations in GROMACS + PLUMED to compute potentials of mean
force (PMFs) along the centre-of-mass distance between two molecules.

It supports both coarse-grained Martini and atomistic OPLS-AA systems, offers
optional projection onto a second collective variable for a 2D free-energy
surface, and generates ready-to-run shell scripts rather than running
simulations directly, so you remain in control of submission.

---

## Requirements

- Python ≥ 3.9
- GROMACS compiled with PLUMED patch (`gmx_mpi` with PLUMED support)
- PLUMED ≥ 2.4
- `pip`

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/alessandri-group/metadpmf.git
cd metadpmf
pip install -e .
```

Verify the installation:

```bash
metadpmf --help
```

For the docs site:

```bash
pip install mkdocs mkdocs-material
mkdocs serve        # preview at http://127.0.0.1:8000
mkdocs build        # build static site in site/
```

---

## What you need before starting

| File | Description |
|---|---|
| `start.gro` | NPT-equilibrated simulation box containing both molecules in solvent |
| `system.top` | GROMACS topology (with `#include` paths to your `.itp` files) |

Copy `config.yaml` from the package root into your simulation directory,
edit it for your system, and follow the [tutorial](tutorial.md).
