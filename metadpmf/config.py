"""Load, validate, and render configuration for metadpmf."""

from pathlib import Path

import yaml

TEMPLATES_DIR = Path(__file__).parent / "templates"

R = 0.008314  # kJ/mol/K


def load(config_path) -> dict:
    config_path = Path(config_path).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    cfg = _apply_defaults(cfg)
    _validate(cfg)
    cfg["_dir"] = config_path.parent
    return cfg


def _apply_defaults(cfg: dict) -> dict:
    cfg.setdefault("temperature", 298.15)

    cfg.setdefault("molecule", {})
    cfg["molecule"].setdefault("mol1_atoms", [])
    cfg["molecule"].setdefault("mol2_atoms", [])

    cfg.setdefault("plumed", {})
    p = cfg["plumed"]
    p.setdefault("pace", 500)
    p.setdefault("height", 1.0)
    p.setdefault("sigma", 0.05)
    p.setdefault("biasfactor", 5)
    p.setdefault("grid_min", 0.1)
    p.setdefault("grid_max", 3.0)
    p.setdefault("grid_bin", 300)
    p.setdefault("wall_at", 2.0)
    p.setdefault("wall_kappa", 200.0)
    p.setdefault("print_stride", 5000)

    cfg.setdefault("mdp", {})
    cfg["mdp"].setdefault("nsteps", 50000000)
    cfg["mdp"].setdefault("dt", 0.020)

    cfg.setdefault("fes", {})
    cfg["fes"].setdefault("min", 0.2)
    cfg["fes"].setdefault("max", 1.7)
    cfg["fes"].setdefault("bins", 51)
    cfg["fes"].setdefault("block_max", 1000)

    cfg.setdefault("pmf", {})
    cfg["pmf"].setdefault("shift_range", [1.5, 1.7])
    cfg["pmf"].setdefault("xrange", None)
    cfg["pmf"].setdefault("yrange", None)

    cfg.setdefault("paths", {})
    cfg["paths"].setdefault("gro", "start.gro")
    cfg["paths"].setdefault("top", "system.top")
    cfg["paths"].setdefault("mdp", None)
    cfg["paths"].setdefault("traj", "metad.xtc")

    cfg.setdefault("backend", "local")
    cfg.setdefault("backends", {})
    cfg["backends"].setdefault("local", {})
    cfg["backends"]["local"].setdefault("gmx", "gmx")
    cfg["backends"]["local"].setdefault("plumed", "plumed")
    cfg["backends"]["local"].setdefault("gmxrc", None)

    return cfg


def _validate(cfg: dict) -> None:
    errors = []

    mol = cfg.get("molecule", {})
    if not mol.get("mol1_atoms"):
        errors.append("molecule.mol1_atoms is required and must not be empty")
    if not mol.get("mol2_atoms"):
        errors.append("molecule.mol2_atoms is required and must not be empty")

    backend = cfg.get("backend", "local")
    if backend not in ("local", "slurm"):
        errors.append("backend must be 'local' or 'slurm'")
    if backend == "slurm":
        header = cfg.get("backends", {}).get("slurm", {}).get("header")
        if not header:
            errors.append("backends.slurm.header is required when backend=slurm")

    if errors:
        raise ValueError("Config errors:\n" + "\n".join(f"  - {e}" for e in errors))


def backend_cmds(cfg: dict) -> dict:
    b = cfg["backend"]
    if b == "local":
        be = cfg["backends"]["local"]
        return {
            "gmx":    be.get("gmx", "gmx_mpi"),
            "plumed": be.get("plumed", "plumed"),
            "gmxrc":  be.get("gmxrc"),
        }
    else:
        be = cfg["backends"]["slurm"]
        return {
            "gmx":    be.get("gmx", "gmx_mpi"),
            "plumed": be.get("plumed", "plumed"),
            "gmxrc":  be.get("gmxrc"),
        }


def kbt(cfg: dict) -> float:
    return R * cfg["temperature"]


def _atoms_str(atoms: list) -> str:
    """Convert a list of atom indices to a comma-separated PLUMED atom string."""
    return ",".join(str(a) for a in atoms)


def render_plumed(cfg: dict) -> str:
    """Generate plumed.dat content for the production metadynamics run."""
    mol = cfg["molecule"]
    p = cfg["plumed"]
    a1 = _atoms_str(mol["mol1_atoms"])
    a2 = _atoms_str(mol["mol2_atoms"])
    T = cfg["temperature"]

    lines = [
        "# treat each molecule as whole",
        f"WHOLEMOLECULES ENTITY0={a1} ENTITY1={a2}",
        "",
        "# define atoms for distance",
        f"c1: CENTER ATOMS={a1}",
        f"c2: CENTER ATOMS={a2}",
        "dist: DISTANCE ATOMS=c1,c2",
        "",
        "# upper wall to limit sampling range",
        f"uwall: UPPER_WALLS ARG=dist AT={p['wall_at']} KAPPA={p['wall_kappa']}",
        "",
        "# metadynamics",
        (
            f"metad: METAD ARG=dist"
            f" PACE={p['pace']}"
            f" HEIGHT={p['height']}"
            f" SIGMA={p['sigma']}"
            f" FILE=HILLS"
            f" BIASFACTOR={p['biasfactor']}"
            f" GRID_MIN={p['grid_min']}"
            f" GRID_MAX={p['grid_max']}"
            f" TEMP={T}"
            f" GRID_BIN={p['grid_bin']}"
        ),
        "",
        "# output",
        f"PRINT STRIDE={p['print_stride']} ARG=dist,metad.bias,uwall.bias FILE=COLVAR",
    ]
    return "\n".join(lines) + "\n"


def render_plumed_analysis(cfg: dict) -> str:
    """Generate plumed_analysis.dat content for the plumed driver reweighting step."""
    mol = cfg["molecule"]
    p = cfg["plumed"]
    a1 = _atoms_str(mol["mol1_atoms"])
    a2 = _atoms_str(mol["mol2_atoms"])
    T = cfg["temperature"]

    lines = [
        "# plumed input for reweighting analysis",
        "RESTART",
        "",
        "# treat each molecule as whole",
        f"WHOLEMOLECULES ENTITY0={a1} ENTITY1={a2}",
        "",
        "# define atoms for distance",
        f"c1: CENTER ATOMS={a1}",
        f"c2: CENTER ATOMS={a2}",
        "dist: DISTANCE ATOMS=c1,c2",
        "",
        "# upper wall (same geometry as production run)",
        f"uwall: UPPER_WALLS ARG=dist AT={p['wall_at']} KAPPA={p['wall_kappa']}",
        "",
        "# read deposited Gaussians from HILLS (PACE >> nframes so no new hills added)",
        (
            f"metad: METAD ARG=dist"
            f" PACE=1000000000"
            f" HEIGHT={p['height']}"
            f" SIGMA={p['sigma']}"
            f" FILE=HILLS"
            f" BIASFACTOR={p['biasfactor']}"
            f" GRID_MIN={p['grid_min']}"
            f" GRID_MAX={p['grid_max']}"
            f" TEMP={T}"
            f" GRID_BIN={p['grid_bin']}"
        ),
        "",
        "# print distance and bias for every frame",
        "PRINT STRIDE=1 ARG=dist,metad.bias FILE=COLVAR",
    ]
    return "\n".join(lines) + "\n"


def render_mdp(cfg: dict) -> str:
    """Return MDP file content, either from config-specified file or built-in template."""
    mdp_path = cfg["paths"].get("mdp")
    if mdp_path is not None:
        return Path(cfg["_dir"] / mdp_path).read_text()

    template = (TEMPLATES_DIR / "md.mdp").read_text()
    replacements = {
        "NSTEPS": str(cfg["mdp"]["nsteps"]),
        "TEMP":   str(cfg["temperature"]),
        "DT":     str(cfg["mdp"]["dt"]),
    }
    for token, value in replacements.items():
        template = template.replace(token, value)
    return template
