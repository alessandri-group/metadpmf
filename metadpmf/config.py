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

    # forcefield is required (no default) — validated in _validate.
    # dt default is forcefield-aware: Martini CG runs at 20 fs, atomistic at 2 fs.
    # nsteps stays 50M for both → 1 us (Martini) vs 100 ns (atomistic), sensible per method.
    _dt_default = {"martini": 0.020, "opls": 0.002, "gaff2": 0.002}.get(cfg.get("forcefield"))
    cfg.setdefault("mdp", {})
    cfg["mdp"].setdefault("nsteps", 50000000)
    if _dt_default is not None:
        cfg["mdp"].setdefault("dt", _dt_default)

    cfg.setdefault("fes", {})
    cfg["fes"].setdefault("min", 0.2)
    cfg["fes"].setdefault("max", 1.7)
    cfg["fes"].setdefault("bins", 51)
    cfg["fes"].setdefault("block_max", 1000)

    cfg.setdefault("pmf", {})
    cfg["pmf"].setdefault("shift_range", [1.5, 1.7])
    cfg["pmf"].setdefault("xrange", None)
    cfg["pmf"].setdefault("yrange", None)
    cfg["pmf"].setdefault("levels_2d", [-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10])

    cfg.setdefault("cv2", {})
    cv2 = cfg["cv2"]
    cv2.setdefault("enabled", False)
    cv2.setdefault("type", "costheta")   # 'costheta' (built-in generator) or 'custom'
    cv2.setdefault("bias", False)        # False = project only; True = also bias cv2
    cv2.setdefault("sigma", 0.1)         # metad Gaussian width on cv2 (only if bias)
    cv2.setdefault("fold", "none")       # none | abs | sq — see NOTES.md
    cv2.setdefault("label", "cos θ")
    cv2.setdefault("min", -1.0)
    cv2.setdefault("max", 1.0)
    cv2.setdefault("bins", 51)
    cv2.setdefault("plumed_analysis", "plumed_analysis_2d.dat")
    # plane-defining atoms for the cosθ generator default to the molecule's
    # centre atoms when those number exactly 3 (the common 3-bead ring case).
    _m1 = cfg["molecule"].get("mol1_atoms", [])
    _m2 = cfg["molecule"].get("mol2_atoms", [])
    cv2.setdefault("mol1_plane", list(_m1) if len(_m1) == 3 else None)
    cv2.setdefault("mol2_plane", list(_m2) if len(_m2) == 3 else None)

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

    ff = cfg.get("forcefield")
    if ff is None:
        errors.append("forcefield is required: set it to 'martini', 'opls', or 'gaff2'")
    elif ff not in ("martini", "opls", "gaff2"):
        errors.append(f"forcefield '{ff}' is not known: must be 'martini', 'opls', or 'gaff2'")

    backend = cfg.get("backend", "local")
    if backend not in ("local", "slurm"):
        errors.append("backend must be 'local' or 'slurm'")
    if backend == "slurm":
        header = cfg.get("backends", {}).get("slurm", {}).get("header")
        if not header:
            errors.append("backends.slurm.header is required when backend=slurm")

    fes      = cfg.get("fes", {})
    pmf      = cfg.get("pmf", {})
    p        = cfg.get("plumed", {})
    fes_min  = fes.get("min", 0.2)
    fes_max  = fes.get("max", 1.7)
    wall_at  = p.get("wall_at", 2.0)
    sr       = pmf.get("shift_range", [1.5, 1.7])

    if fes_max > wall_at:
        errors.append(
            f"fes.max ({fes_max}) > plumed.wall_at ({wall_at}): "
            "the FES histogram extends beyond the upper wall where sampling is suppressed"
        )
    if sr[0] < fes_min or sr[1] > fes_max:
        errors.append(
            f"pmf.shift_range {sr} must lie within [fes.min, fes.max] = [{fes_min}, {fes_max}]"
        )
    if sr[0] >= sr[1]:
        errors.append(
            f"pmf.shift_range {sr}: first value must be less than second"
        )

    cv2 = cfg.get("cv2", {})
    if cv2.get("enabled"):
        ctype = cv2.get("type", "costheta")
        if ctype not in ("costheta", "custom"):
            errors.append(f"cv2.type '{ctype}' is not known: must be 'costheta' or 'custom'")
        if cv2.get("fold") not in ("none", "abs", "sq"):
            errors.append(f"cv2.fold '{cv2.get('fold')}' is not known: must be 'none', 'abs', or 'sq'")
        if ctype == "costheta":
            for key in ("mol1_plane", "mol2_plane"):
                atoms = cv2.get(key)
                if not atoms or len(atoms) != 3:
                    errors.append(
                        f"cv2.{key} must list exactly 3 atoms for cv2.type=costheta (got {atoms})"
                    )
        if cv2.get("bias"):
            if cv2.get("min") >= cv2.get("max"):
                errors.append(f"cv2.min ({cv2.get('min')}) must be < cv2.max ({cv2.get('max')})")
            if cv2.get("sigma", 0) <= 0:
                errors.append("cv2.sigma must be > 0 when cv2.bias is true")

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
        gmx = be.get("gmx", "gmx_mpi")
        return {
            "gmx":    gmx,
            "plumed": be.get("plumed", "plumed"),
            "gmxrc":  be.get("gmxrc"),
            "mdrun":  be.get("mdrun", f"{gmx} mdrun"),
            "grompp": be.get("grompp", f"{gmx} grompp"),
        }


def kbt(cfg: dict) -> float:
    return R * cfg["temperature"]


def _atoms_str(atoms: list) -> str:
    """Convert a list of atom indices to a comma-separated PLUMED atom string."""
    return ",".join(str(a) for a in atoms)


def cv2_is_costheta(cfg: dict) -> bool:
    """True when the second CV is enabled and uses the built-in cosθ generator."""
    cv2 = cfg["cv2"]
    return bool(cv2.get("enabled")) and cv2.get("type", "costheta") == "costheta"


def cv2_is_biased(cfg: dict) -> bool:
    """True for Mode C: cosθ is generated *and* biased by METAD."""
    return cv2_is_costheta(cfg) and bool(cfg["cv2"].get("bias"))


def _render_cv2_block(cfg: dict):
    """Return (plumed_lines, cv_name) defining cosθ between the two ring-plane
    normals. The normal of each ring is the cross product of two edge vectors
    built from the three plane atoms. Only valid for cv2.type == 'costheta'."""
    cv2 = cfg["cv2"]
    a = cv2["mol1_plane"]
    b = cv2["mol2_plane"]
    lines = [
        "# --- CV2: cosθ between ring-plane normals ---",
        "# edge vectors within each ring",
        f"d12: DISTANCE ATOMS={a[0]},{a[1]} COMPONENTS",
        f"d23: DISTANCE ATOMS={a[1]},{a[2]} COMPONENTS",
        f"d45: DISTANCE ATOMS={b[0]},{b[1]} COMPONENTS",
        f"d56: DISTANCE ATOMS={b[1]},{b[2]} COMPONENTS",
        "",
        "# ring-plane normals via cross products  n = e1 x e2",
        "n1x: CUSTOM ARG=d12.y,d12.z,d23.y,d23.z VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO",
        "n1y: CUSTOM ARG=d12.z,d12.x,d23.z,d23.x VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO",
        "n1z: CUSTOM ARG=d12.x,d12.y,d23.x,d23.y VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO",
        "n2x: CUSTOM ARG=d45.y,d45.z,d56.y,d56.z VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO",
        "n2y: CUSTOM ARG=d45.z,d45.x,d56.z,d56.x VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO",
        "n2z: CUSTOM ARG=d45.x,d45.y,d56.x,d56.y VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO",
        "",
        "# cosθ = (n1·n2)/(|n1||n2|), range [-1, 1]",
        (
            "costheta: CUSTOM ARG=n1x,n1y,n1z,n2x,n2y,n2z VAR=a,b,c,d,f,g "
            "FUNC=(a*d+b*f+c*g)/sqrt((a^2+b^2+c^2)*(d^2+f^2+g^2)) PERIODIC=NO"
        ),
    ]
    fold = cv2.get("fold", "none")
    if fold == "abs":
        lines.append("cv2f: CUSTOM ARG=costheta VAR=a FUNC=abs(a) PERIODIC=NO")
        return lines, "cv2f"
    if fold == "sq":
        lines.append("cv2f: CUSTOM ARG=costheta VAR=a FUNC=a*a PERIODIC=NO")
        return lines, "cv2f"
    return lines, "costheta"


def _metad_line(cfg: dict, pace, cv_name=None) -> str:
    """Build a METAD line, optionally biasing a second CV (comma-joined args)."""
    p = cfg["plumed"]
    cv2 = cfg["cv2"]
    T = cfg["temperature"]
    if cv_name is not None:
        arg   = f"dist,{cv_name}"
        sigma = f"{p['sigma']},{cv2['sigma']}"
        gmin  = f"{p['grid_min']},{cv2['min']}"
        gmax  = f"{p['grid_max']},{cv2['max']}"
        gbin  = f"{p['grid_bin']},{cv2['bins']}"
    else:
        arg, sigma = "dist", str(p["sigma"])
        gmin, gmax, gbin = p["grid_min"], p["grid_max"], p["grid_bin"]
    return (
        f"metad: METAD ARG={arg}"
        f" PACE={pace}"
        f" HEIGHT={p['height']}"
        f" SIGMA={sigma}"
        f" FILE=HILLS"
        f" BIASFACTOR={p['biasfactor']}"
        f" GRID_MIN={gmin}"
        f" GRID_MAX={gmax}"
        f" TEMP={T}"
        f" GRID_BIN={gbin}"
    )


def render_plumed(cfg: dict) -> str:
    """Generate plumed.dat content for the production metadynamics run.

    Mode C (cv2 enabled + biased) injects the cosθ block and biases it; Modes A
    and B leave the production run on the distance CV only (cosθ, if any, is
    computed later during reweighting)."""
    mol = cfg["molecule"]
    p = cfg["plumed"]
    a1 = _atoms_str(mol["mol1_atoms"])
    a2 = _atoms_str(mol["mol2_atoms"])

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
    ]

    if cv2_is_biased(cfg):
        block, cv_name = _render_cv2_block(cfg)
        lines += ["", *block]
        metad = _metad_line(cfg, p["pace"], cv_name)
        print_arg = f"dist,{cv_name},metad.bias,uwall.bias"
    else:
        metad = _metad_line(cfg, p["pace"])
        print_arg = "dist,metad.bias,uwall.bias"

    lines += [
        "",
        "# metadynamics",
        metad,
        "",
        "# output",
        f"PRINT STRIDE={p['print_stride']} ARG={print_arg} FILE=COLVAR",
    ]
    return "\n".join(lines) + "\n"


def render_plumed_analysis(cfg: dict) -> str:
    """Generate plumed_analysis.dat content for the plumed driver reweighting step.

    For cv2.type == 'costheta' the cosθ block is generated and printed (Modes B
    and C); in Mode C the METAD also biases cosθ to reconstruct the 2D bias."""
    mol = cfg["molecule"]
    p = cfg["plumed"]
    a1 = _atoms_str(mol["mol1_atoms"])
    a2 = _atoms_str(mol["mol2_atoms"])

    two_d  = cv2_is_costheta(cfg)
    biased = cv2_is_biased(cfg)

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
    ]

    cv_name = None
    if two_d:
        block, cv_name = _render_cv2_block(cfg)
        lines += ["", *block]

    metad = _metad_line(cfg, 1000000000, cv_name if biased else None)
    print_arg = f"dist,{cv_name},metad.bias" if two_d else "dist,metad.bias"

    lines += [
        "",
        "# read deposited Gaussians from HILLS (PACE >> nframes so no new hills added)",
        metad,
        "",
        "# print distance (and cv2) and bias for every frame",
        f"PRINT STRIDE=1 ARG={print_arg} FILE=COLVAR",
    ]
    return "\n".join(lines) + "\n"


def render_mdp(cfg: dict) -> str:
    """Return MDP file content, either from config-specified file or built-in template."""
    mdp_path = cfg["paths"].get("mdp")
    if mdp_path is not None:
        return Path(cfg["_dir"] / mdp_path).read_text()

    template = (TEMPLATES_DIR / f"md_{cfg['forcefield']}.mdp").read_text()
    replacements = {
        "NSTEPS": str(cfg["mdp"]["nsteps"]),
        "TEMP":   str(cfg["temperature"]),
        "DT":     str(cfg["mdp"]["dt"]),
    }
    for token, value in replacements.items():
        template = template.replace(token, value)
    return template
