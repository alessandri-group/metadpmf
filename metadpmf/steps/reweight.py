"""Step 2: write analysis/plumed_analysis.dat and reweight.sh.

Creates the analysis/ subdirectory, writes the PLUMED input for reweighting
(PACE=1e9 so no new hills are deposited, STRIDE=1 so every frame is printed),
and generates a shell script that copies HILLS into analysis/ then runs
plumed driver to produce analysis/COLVAR.

1D mode (cv2.enabled: false):
  COLVAR columns: time, distance, bias
  plumed_analysis.dat is generated automatically from config.

2D mode (cv2.enabled: true):
  COLVAR columns: time, distance, cv2, bias
  cv2.type == 'costheta': plumed_analysis.dat is generated from config (the
    cosθ block is built from cv2.mol1_plane / cv2.mol2_plane).
  cv2.type == 'custom':   plumed_analysis.dat is copied from the user-provided
    file (cv2.plumed_analysis); the user defines CV2 and prints
    dist, cv2, metad.bias to COLVAR.
"""

import shutil
from pathlib import Path

from metadpmf.config import backend_cmds, render_plumed_analysis


def run(cfg: dict) -> None:
    base     = Path(cfg["_dir"])
    cv2      = cfg["cv2"]
    two_d    = cv2["enabled"]
    custom   = two_d and cv2.get("type", "costheta") == "custom"
    dir_name = "analysis_2d" if two_d else "analysis"
    anal_dir = base / dir_name
    anal_dir.mkdir(exist_ok=True)

    plumed_dat = anal_dir / "plumed_analysis.dat"

    if custom:
        # type: custom — copy user-provided PLUMED input into analysis_2d/
        src = base / cv2["plumed_analysis"]
        if not src.exists():
            raise FileNotFoundError(
                f"cv2.plumed_analysis file not found: {src}\n"
                "Create a PLUMED input file that defines CV2 and prints "
                "dist, cv2, metad.bias to COLVAR, then set cv2.plumed_analysis "
                "to its path in config.yaml (or use cv2.type: costheta)."
            )
        shutil.copy(src, plumed_dat)
        print(f"Copied: {cv2['plumed_analysis']} → {dir_name}/plumed_analysis.dat")
        print("Note: COLVAR will have columns: time, distance, cv2, bias")
    else:
        # 1D, or 2D with the built-in cosθ generator: write from config
        plumed_dat.write_text(render_plumed_analysis(cfg))
        print(f"Written: {dir_name}/plumed_analysis.dat")
        if two_d:
            print("Note: COLVAR will have columns: time, distance, cv2, bias")

    cmds = backend_cmds(cfg)
    traj = cfg["paths"].get("traj", "metad.xtc")

    if cfg["backend"] == "slurm":
        _write_slurm(cfg, base, traj, cmds, dir_name)
    else:
        _write_local(base, traj, cmds, dir_name)


def _source_line(gmxrc):
    return f"source {gmxrc}" if gmxrc else "# (GROMACS/PLUMED already in PATH)"


def _write_local(base, traj, cmds, dir_name):
    gmxrc  = cmds.get("gmxrc")
    plumed = cmds["plumed"]

    lines = [
        "#!/bin/bash",
        "# metadpmf: reweight trajectory with plumed driver",
        'cd "$(dirname "$0")"',
        "",
        _source_line(gmxrc),
        "",
        f"# copy HILLS into {dir_name}/ so plumed driver can read it",
        f"cp HILLS {dir_name}/",
        "",
        f"cd {dir_name}",
        f"{plumed} driver --plumed plumed_analysis.dat --mf_xtc ../{traj}",
    ]
    script = base / "reweight.sh"
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o755)
    print(f"Written: reweight.sh")
    print(f"Run with: bash reweight.sh")


def _write_slurm(cfg, base, traj, cmds, dir_name):
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
        f"cp HILLS {dir_name}/",
        "",
        f"cd {dir_name}",
        f"{plumed} driver --plumed plumed_analysis.dat --mf_xtc ../{traj}",
    ]
    script = base / "reweight.sh"
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o755)
    print(f"Written: reweight.sh")
    print(f"Submit with: sbatch reweight.sh")
