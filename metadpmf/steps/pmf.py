"""Step 4: apply Jacobian correction, shift PMF to zero, and plot.

The raw FES from the block analysis is a free-energy profile along the
centre-of-mass distance r. Because the spherical volume element scales as
r^2, there is more accessible phase-space volume at larger separations even
for a flat interaction. The Jacobian (translational entropy) correction
removes this purely geometric bias:

    DeltaG(r) = FES(r) + 2 R T ln(r)

The corrected profile is then shifted to zero by subtracting the mean value
of DeltaG over a configurable reference range (pmf.shift_range, default
1.5–1.7 nm). This range should correspond to a plateau region where the two
molecules no longer interact; the mean rather than a single point is used to
reduce sensitivity to noise at large separations.

1D mode (cv2.enabled: false):
  Reads:  analysis/fes.dat            (columns: dist, fes, error)
  Writes: analysis/fes_corr_and_shifted.dat
          analysis/pmf.pdf

2D mode (cv2.enabled: true):
  Reads:  analysis/fes.dat            (columns: cv1, cv2, fes, error)
  The Jacobian correction is applied on the CV1 (distance) axis only.
  Writes: analysis/fes_2d_corr.dat
          analysis/pmf_2d.pdf
"""

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from metadpmf.config import R


def run(cfg: dict, temperature: float, marginal: bool = False) -> None:
    base      = Path(cfg["_dir"])
    two_d     = cfg["cv2"]["enabled"]
    dir_name  = "analysis_2d" if two_d else "analysis"
    anal_dir  = base / dir_name

    if marginal and not two_d:
        raise ValueError(
            "--1d only applies to a 2D run (cv2.enabled: true). For a 1D run "
            "just use 'metadpmf pmf' without --1d."
        )

    if marginal:
        fes_path = anal_dir / "fes_1d.dat"
        if not fes_path.exists():
            raise FileNotFoundError(
                f"{dir_name}/fes_1d.dat not found. Run 'metadpmf fes --1d' first."
            )
        _run_1d(cfg, temperature, anal_dir, fes_path,
                corr_name="fes_1d_corr_and_shifted.dat", pdf_name="pmf_1d.pdf")
        return

    fes_path = anal_dir / "fes.dat"
    if not fes_path.exists():
        raise FileNotFoundError(
            f"{dir_name}/fes.dat not found. Run 'metadpmf fes' first."
        )

    if two_d:
        _run_2d(cfg, temperature, anal_dir, fes_path)
    else:
        _run_1d(cfg, temperature, anal_dir, fes_path)


# ---------------------------------------------------------------------------
# 1D PMF
# ---------------------------------------------------------------------------

def _run_1d(cfg, temperature, anal_dir, fes_path,
            corr_name="fes_corr_and_shifted.dat", pdf_name="pmf.pdf"):
    dname = anal_dir.name
    distances, fes, errors = _parse_fes_1d(fes_path)

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

    _write_1d(anal_dir / corr_name, distances, shifted, errors)
    print(f"Written: {dname}/{corr_name}")

    _plot_1d(
        distances, shifted, errors,
        anal_dir / pdf_name,
        xrange=cfg["pmf"].get("xrange"),
        yrange=cfg["pmf"].get("yrange"),
    )
    print(f"Written: {dname}/{pdf_name}")

    min_val = min(shifted)
    r_min   = distances[shifted.index(min_val)]
    print(f"\nDelta G (well depth) = {min_val:.2f} kJ/mol  at r = {r_min:.3f} nm")
    print(f"Zero reference: mean corrected FES in [{r_lo:.2f}, {r_hi:.2f}] nm")


# ---------------------------------------------------------------------------
# 2D PMF
# ---------------------------------------------------------------------------

def _run_2d(cfg, temperature, anal_dir, fes_path):
    nbin1   = cfg["fes"]["bins"]
    nbin2   = cfg["cv2"]["bins"]
    cv2_cfg = cfg["cv2"]

    FES, ERR, cv1_edges, cv2_edges = _parse_fes_2d(fes_path, nbin1, nbin2)

    KBT = R * temperature

    # Jacobian correction on CV1 (distance) axis only
    for i1 in range(nbin1):
        r = cv1_edges[i1]
        if r > 0:
            FES[i1, :] += 2.0 * KBT * math.log(r)

    # shift: mean of corrected FES over pmf.shift_range × all CV2 bins
    shift_range = cfg["pmf"].get("shift_range", [1.5, 1.7])
    r_lo, r_hi  = shift_range
    ref_vals    = []
    for i1 in range(nbin1):
        if r_lo <= cv1_edges[i1] <= r_hi:
            for i2 in range(nbin2):
                if math.isfinite(FES[i1, i2]):
                    ref_vals.append(FES[i1, i2])

    if not ref_vals:
        raise ValueError(
            f"No finite FES bins found in shift_range [{r_lo}, {r_hi}] nm. "
            "Adjust pmf.shift_range in config.yaml."
        )

    shift = sum(ref_vals) / len(ref_vals)
    FES   = FES - shift
    print(f"2D FES shifted by {shift:.4f} kJ/mol")

    _write_2d(anal_dir / "fes_2d_corr.dat", FES, cv1_edges, cv2_edges, nbin1, nbin2)
    print("Written: analysis/fes_2d_corr.dat")

    _plot_2d(FES, cv1_edges, cv2_edges, cv2_cfg, anal_dir / "pmf_2d.pdf", cfg["pmf"])
    print("Written: analysis/pmf_2d.pdf")


# ---------------------------------------------------------------------------
# Internal helpers — 1D
# ---------------------------------------------------------------------------

def _parse_fes_1d(path: Path):
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


def _write_1d(path: Path, distances, values, errors):
    with open(path, "w") as f:
        for r, v, e in zip(distances, values, errors):
            f.write(f"{r:.9f}  {v:.9f}  {e:.9f}\n")


def _plot_1d(distances, shifted, errors, out_path, xrange=None, yrange=None):
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


# ---------------------------------------------------------------------------
# Internal helpers — 2D
# ---------------------------------------------------------------------------

def _parse_fes_2d(path: Path, nbin1: int, nbin2: int):
    """Read 2D fes.dat into arrays of shape (nbin1, nbin2).

    File format (from _write_fes_2d): cv1 fast-varying, cv2 slow-varying,
    blank lines between CV2 groups. Values of 'Infinity' become np.inf.
    """
    FES       = np.full((nbin1, nbin2), np.inf)
    ERR       = np.full((nbin1, nbin2), np.inf)
    cv1_edges = np.zeros(nbin1)
    cv2_edges = np.zeros(nbin2)

    count = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            i1 = count % nbin1
            i2 = count // nbin1
            if i2 >= nbin2:
                break
            cv1_edges[i1] = float(parts[0])
            cv2_edges[i2] = float(parts[1])
            if parts[2] != "Infinity":
                FES[i1, i2] = float(parts[2])
                if len(parts) >= 4 and parts[3] != "Infinity":
                    ERR[i1, i2] = float(parts[3])
            count += 1

    return FES, ERR, cv1_edges, cv2_edges


def _write_2d(path: Path, FES, cv1_edges, cv2_edges, nbin1: int, nbin2: int):
    with open(path, "w") as f:
        count = 0
        for i2 in range(nbin2):
            for i1 in range(nbin1):
                v = FES[i1, i2]
                if math.isfinite(v):
                    f.write(f"{cv1_edges[i1]:.9f}  {cv2_edges[i2]:.9f}  {v:.9f}\n")
                else:
                    f.write(f"{cv1_edges[i1]:.9f}  {cv2_edges[i2]:.9f}  Infinity\n")
                count += 1
            if i2 < nbin2 - 1:
                f.write("\n")


def _plot_2d(FES, cv1_edges, cv2_edges, cv2_cfg, out_path, pmf_cfg):
    """Contour plot of the 2D corrected PMF."""
    from matplotlib.colors import TwoSlopeNorm

    # Unsampled bins (FES = inf) are pushed above the top level so that, with
    # extend='both', they saturate to the over-colour rather than leaving the
    # facecolor (grey) showing through.
    levels = pmf_cfg.get("levels_2d", [-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10])
    FES_plot = FES.copy()
    FES_plot[~np.isfinite(FES_plot)] = max(levels) + 1.0

    # meshgrid: X varies along CV1, Y along CV2
    X, Y = np.meshgrid(cv1_edges, cv2_edges)  # shape: (nbin2, nbin1)

    cmap   = plt.get_cmap("coolwarm").copy()
    # Out-of-range values saturate to the colormap ends instead of going
    # unfilled (which showed the grey facecolor). set_over/under make this
    # explicit so the over-region (incl. unsampled bins) is solid dark red.
    cmap.set_over(cmap(1.0))
    cmap.set_under(cmap(0.0))
    # TwoSlopeNorm pins vcenter=0 to the colormap midpoint so the red/blue
    # transition is always at zero regardless of asymmetric level choices
    norm   = TwoSlopeNorm(vmin=min(levels), vcenter=0, vmax=max(levels))

    fig, ax = plt.subplots(figsize=(9.6, 7))
    ax.set_facecolor(cmap(norm(0.0)))

    cf = ax.contourf(X, Y, FES_plot.T, levels=levels, cmap=cmap, norm=norm,
                     extend="both")
    cb = fig.colorbar(cf, extend="both")
    cb.set_label(r"$\Delta G$ (kJ mol$^{-1}$)", fontsize=13)

    ax.set_xlabel("r (nm)", fontsize=14)
    ax.set_ylabel(cv2_cfg.get("label", "S"), fontsize=14)
    ax.tick_params(labelsize=12, direction="in", top=True, right=True,
                   which="both", length=4)

    # y-axis always spans the CV2 grid — use cv2.min/max, not pmf.yrange
    # (pmf.yrange controls the ΔG axis in the 1D plot only)
    xrange = pmf_cfg.get("xrange")

    fig.tight_layout()

    if xrange is not None:
        ax.set_xlim(xrange)
    ax.set_ylim([cv2_cfg["min"], cv2_cfg["max"]])

    fig.savefig(out_path, dpi=150)
    plt.close(fig)
