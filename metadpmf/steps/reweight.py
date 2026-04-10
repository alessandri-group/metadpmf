"""Step 2: write analysis/plumed_analysis.dat and reweight.sh.

Creates the analysis/ subdirectory, writes the PLUMED input for reweighting
(PACE=1e9 so no new hills are deposited, STRIDE=1 so every frame is printed),
and generates a shell script that copies HILLS into analysis/ then runs
plumed driver to produce analysis/COLVAR with columns: time, distance, bias.
"""

from pathlib import Path

from metadpmf.config import backend_cmds, render_plumed_analysis


def run(cfg: dict) -> None:
    base = Path(cfg["_dir"])
    analysis_dir = base / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    # write plumed_analysis.dat into analysis/
    plumed_dat = analysis_dir / "plumed_analysis.dat"
    plumed_dat.write_text(render_plumed_analysis(cfg))
    print(f"Written: analysis/plumed_analysis.dat")

    cmds = backend_cmds(cfg)
    traj = cfg["paths"].get("traj", "metad.xtc")

    if cfg["backend"] == "slurm":
        _write_slurm(cfg, base, traj, cmds)
    else:
        _write_local(base, traj, cmds)


def _source_line(gmxrc):
    return f"source {gmxrc}" if gmxrc else "# (GROMACS/PLUMED already in PATH)"


def _write_local(base, traj, cmds):
    gmxrc  = cmds.get("gmxrc")
    plumed = cmds["plumed"]

    lines = [
        "#!/bin/bash",
        "# metadpmf: reweight trajectory with plumed driver",
        'cd "$(dirname "$0")"',
        "",
        _source_line(gmxrc),
        "",
        "# copy HILLS into analysis/ so plumed driver can read it",
        "cp HILLS analysis/",
        "",
        "cd analysis",
        f"{plumed} driver --plumed plumed_analysis.dat --mf_xtc ../{traj}",
    ]
    script = base / "reweight.sh"
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o755)
    print(f"Written: reweight.sh")
    print(f"Run with: bash reweight.sh")


def _write_slurm(cfg, base, traj, cmds):
    gmxrc       = cmds.get("gmxrc")
    plumed      = cmds["plumed"]
    header_text = (base / cfg["backends"]["slurm"]["header"]).read_text().rstrip()

    lines = [
        "#!/bin/bash",
        "#SBATCH --job-name=reweight",
        "#SBATCH --output=reweight.%J.out",
        "",
        header_text,
        "",
        'cd "$(dirname "$0")"',
        "",
        _source_line(gmxrc),
        "",
        "cp HILLS analysis/",
        "",
        "cd analysis",
        f"{plumed} driver --plumed plumed_analysis.dat --mf_xtc ../{traj}",
    ]
    script = base / "reweight.sh"
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o755)
    print(f"Written: reweight.sh")
    print(f"Submit with: sbatch reweight.sh")
