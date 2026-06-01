<p align="center">
  <img src="docs/assets/logo.png" width="200" alt="metadpmf logo"/>
</p>

<h1 align="center">metadpmf</h1>

<p align="center">
  Metadynamics PMF workflow for GROMACS + PLUMED
</p>

---

**metadpmf** is a Python CLI package for running and analysing well-tempered metadynamics simulations in GROMACS + PLUMED to compute potentials of mean force (PMFs) along the centre-of-mass distance between two molecules. It supports both coarse-grained Martini and atomistic OPLS-AA systems, optional projection onto a second collective variable for a 2D FES, and generates ready-to-run shell scripts rather than executing simulations directly.

## Features

- Single `config.yaml` per simulation
- Built-in production `mdp` templates for Martini CG and OPLS-AA (selected via `forcefield`)
- 1D PMF, or 2D FES by projecting onto a second (unbiased) collective variable

## Installation

```bash
git clone https://github.com/alessandri-group/metadpmf.git
cd metadpmf
pip install -e .
```

## Quick start

```bash
# 1. Generate PLUMED input and run script
metadpmf run && bash run_metad.sh

# 2. Reweight trajectory
metadpmf reweight && bash reweight.sh

# 3. Block FES analysis
metadpmf fes

# 4. Jacobian correction + plot
metadpmf pmf
```

## Documentation

Build and serve locally:

```bash
pip install mkdocs mkdocs-material
mkdocs serve   # → http://127.0.0.1:8000
```

## Requirements

- Python ≥ 3.9
- GROMACS compiled with PLUMED patch
- PLUMED ≥ 2.4

## Credits

- Logo: [Gemini](https://gemini.google.com)
- See commit history for contributors, coding assisted by [Claude](https://claude.ai)

## License

[MIT](LICENSE)
