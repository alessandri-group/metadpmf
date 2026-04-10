"""Step 4: apply Jacobian correction, shift PMF to zero, and plot.

The raw FES from the block analysis is a 1D free-energy profile along the
centre-of-mass distance r. Because the spherical volume element scales as
r^2, there is more accessible phase-space volume at larger separations even
for a flat interaction. The Jacobian (translational entropy) correction
removes this purely geometric bias:

    DeltaG(r) = FES(r) + 2 R T ln(r)

The corrected profile is then shifted to zero by subtracting the mean value
of DeltaG(r) over a configurable reference range (pmf.shift_range, default
1.5–1.7 nm). This range should correspond to a plateau region where the two
molecules no longer interact; the mean rather than a single point is used to
reduce sensitivity to noise at large separations.

Reads:  analysis/fes.dat
Writes: analysis/fes_corr_and_shifted.dat
        analysis/pmf.pdf
"""

import math
from pathlib import Path

import matplotlib.pyplot as plt

from metadpmf.config import R


def run(cfg: dict, temperature: float) -> None:
    base     = Path(cfg["_dir"])
    anal_dir = base / "analysis"
    fes_path = anal_dir / "fes.dat"

    if not fes_path.exists():
        raise FileNotFoundError(
            "analysis/fes.dat not found. Run 'metadpmf fes' first."
        )

    distances, fes, errors = _parse_fes(fes_path)

    # Jacobian correction
    corrected = [f + 2 * R * temperature * math.log(r)
                 for r, f in zip(distances, fes)]

    # shift to zero using mean corrected value in shift_range
    shift_range = cfg["pmf"].get("shift_range", [1.5, 1.7])
    r_lo, r_hi  = shift_range
    ref_vals    = [c for r, c in zip(distances, corrected) if r_lo <= r <= r_hi]

    if not ref_vals:
        raise ValueError(
            f"No FES bins found in shift_range [{r_lo}, {r_hi}] nm. "
            "Adjust pmf.shift_range in config.yaml."
        )

    shift   = sum(ref_vals) / len(ref_vals)
    shifted = [c - shift for c in corrected]

    _write(anal_dir / "fes_corr_and_shifted.dat", distances, shifted, errors)
    print("Written: analysis/fes_corr_and_shifted.dat")

    _plot(
        distances, shifted, errors, temperature,
        anal_dir / "pmf.pdf",
        xrange=cfg["pmf"].get("xrange"),
        yrange=cfg["pmf"].get("yrange"),
    )
    print("Written: analysis/pmf.pdf")

    min_val = min(shifted)
    r_min   = distances[shifted.index(min_val)]
    print(f"\nDelta G (well depth) = {min_val:.2f} kJ/mol  at r = {r_min:.3f} nm")
    print(f"Zero reference: mean corrected FES in [{r_lo:.2f}, {r_hi:.2f}] nm")


def _parse_fes(path: Path):
    distances, fes, errors = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            distances.append(float(parts[0]))
            fes.append(float(parts[1]))
            errors.append(float(parts[2]))
    return distances, fes, errors


def _write(path: Path, distances, values, errors):
    with open(path, "w") as f:
        for r, v, e in zip(distances, values, errors):
            f.write(f"{r:.9f}  {v:.9f}  {e:.9f}\n")


def _plot(distances, shifted, errors, temperature, out_path, xrange=None, yrange=None):
    fig, ax = plt.subplots(figsize=(6, 5))

    ax.errorbar(
        distances, shifted, yerr=errors,
        fmt="-o", markersize=3, linewidth=1.5, capsize=2,
        color="steelblue", ecolor="lightsteelblue", elinewidth=1,
    )
    ax.axhline(0, color="black", linewidth=0.8, linestyle=":")

    if xrange is not None:
        ax.set_xlim(xrange)
    if yrange is not None:
        ax.set_ylim(yrange)

    ax.set_xlabel("r (nm)", fontsize=14)
    ax.set_ylabel(r"$\Delta G$ (kJ mol$^{-1}$)", fontsize=14)
    ax.tick_params(labelsize=12, direction="in", top=True, right=True,
                   which="both", length=4)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
