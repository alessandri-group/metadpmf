"""Step 1: write plumed.dat and run_metad.sh.

Generates the PLUMED input file for the production metadynamics run and a
shell script that calls grompp followed by mdrun with --plumed. The script
sources GMXRC if configured (for clusters with a PLUMED-patched GROMACS).
"""

from pathlib import Path

from metadpmf.config import backend_cmds, render_mdp, render_plumed


def run(cfg: dict) -> None:
    base = Path(cfg["_dir"])
    cmds = backend_cmds(cfg)

    # write plumed.dat
    plumed_dat = base / "plumed.dat"
    plumed_dat.write_text(render_plumed(cfg))
    print(f"Written: plumed.dat")

    # write metad.mdp (only if not supplied by user)
    if cfg["paths"].get("mdp") is None:
        mdp_out = base / "metad.mdp"
        mdp_out.write_text(render_mdp(cfg))
        print(f"Written: metad.mdp")
        mdp_name = "metad.mdp"
    else:
        mdp_name = cfg["paths"]["mdp"]

    gro  = cfg["paths"]["gro"]
    top  = cfg["paths"]["top"]
    gmx  = cmds["gmx"]

    if cfg["backend"] == "slurm":
        _write_slurm(cfg, base, gmx, gro, top, mdp_name, cmds)
    else:
        _write_local(base, gmx, gro, top, mdp_name, cmds)


def _source_line(gmxrc):
    return f"source {gmxrc}" if gmxrc else "# (GROMACS already in PATH)"


def _write_local(base, gmx, gro, top, mdp_name, cmds):
    gmxrc = cmds.get("gmxrc")
    lines = [
        "#!/bin/bash",
        "# metadpmf: production metadynamics run",
        'cd "$(dirname "$0")"',
        "",
        _source_line(gmxrc),
        "",
        f'{gmx} grompp -p {top} -f {mdp_name} -c {gro} -maxwarn 5 -o metad -po metad',
        f'{gmx} mdrun -v -deffnm metad -plumed plumed.dat',
    ]
    script = base / "run_metad.sh"
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o755)
    print(f"Written: run_metad.sh")
    print(f"Run with: bash run_metad.sh")


def _write_slurm(cfg, base, gmx, gro, top, mdp_name, cmds):
    gmxrc   = cmds.get("gmxrc")
    mdrun   = cmds.get("mdrun", f"{gmx} mdrun")
    grompp  = cmds.get("grompp", f"{gmx} grompp")
    header_text = (base / cfg["backends"]["slurm"]["header"]).read_text().rstrip()

    lines = [
        "#!/bin/bash",
        "#SBATCH --job-name=metad",
        "#SBATCH --output=metad.%J.out",
        "",
        header_text,
        "",
        'cd "$(dirname "$0")"',
        "",
        _source_line(gmxrc),
        "",
        f'{grompp} -p {top} -f {mdp_name} -c {gro} -maxwarn 5 -o metad -po metad',
        f'{mdrun} -v -deffnm metad -plumed plumed.dat',
    ]
    script = base / "run_metad.sh"
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o755)
    print(f"Written: run_metad.sh")
    print(f"Submit with: sbatch run_metad.sh")
