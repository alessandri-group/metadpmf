"""Step 3: block free-energy surface analysis.

Reads analysis/COLVAR (columns: time, distance, bias), computes per-frame
reweighting factors, then performs block analysis to estimate statistical
errors on the FES at a range of block sizes.

For each block size the FES and its standard error are computed via the
block-average method (Flyvbjerg & Petersen, J. Chem. Phys. 91, 461, 1989).
The average error across all FES bins is stored in errors.block so the user
can inspect convergence and select an appropriate block size.

Output written to analysis/:
  dist.weight          — (distance, reweighting_weight) pairs
  errors.block         — (block_size, mean_error) convergence data
  fes.dat              — FES from largest finite block size (use for pmf step)
  blocks/fes.{N}.dat   — FES for every block size N (distance, FES, error)
"""

import math
from pathlib import Path

import numpy as np

from metadpmf.config import kbt as get_kbt


def run(cfg: dict, block_size: int = None) -> None:
    base     = Path(cfg["_dir"])
    anal_dir = base / "analysis"
    colvar   = anal_dir / "COLVAR"

    if not colvar.exists():
        raise FileNotFoundError(
            f"analysis/COLVAR not found. Run 'metadpmf reweight' first."
        )

    KBT = get_kbt(cfg)
    fes_cfg = cfg["fes"]
    gmin    = fes_cfg["min"]
    gmax    = fes_cfg["max"]
    nbin    = fes_cfg["bins"]
    bmax_bs = fes_cfg["block_max"]

    # --- load COLVAR ---
    print("Reading analysis/COLVAR ...")
    distances, weights = _load_colvar(colvar, KBT)
    print(f"  {len(distances)} frames loaded")

    # write dist.weight
    dw_path = anal_dir / "dist.weight"
    np.savetxt(dw_path, np.column_stack([distances, weights]), fmt="%.9f")
    print(f"Written: analysis/dist.weight")

    # --- block FES analysis ---
    blocks_dir = anal_dir / "blocks"
    blocks_dir.mkdir(exist_ok=True)

    if block_size is not None:
        # single block size: write fes.dat directly
        fes_data = _block_fes(distances, weights, gmin, gmax, nbin, KBT, block_size)
        _write_fes(anal_dir / "fes.dat", fes_data)
        print(f"Written: analysis/fes.dat  (block size {block_size})")
        return

    # scan block sizes: range(1, block_max, 10)
    block_sizes = list(range(1, bmax_bs, 10))
    errors = []
    last_finite_fes = None

    print(f"Scanning {len(block_sizes)} block sizes (1 to {bmax_bs - 1}, step 10) ...")

    for bs in block_sizes:
        fes_data = _block_fes(distances, weights, gmin, gmax, nbin, KBT, bs)
        # save to blocks/
        _write_fes(blocks_dir / f"fes.{bs}.dat", fes_data)
        # compute mean error over finite bins
        finite = [(r, f, e) for r, f, e in fes_data if math.isfinite(f)]
        if finite:
            mean_err = sum(e for _, _, e in finite) / len(finite)
            errors.append((bs, mean_err))
            last_finite_fes = [(r, f, e) for r, f, e in finite]

    np.savetxt(anal_dir / "errors.block", errors, fmt="%.9f")
    print(f"Written: analysis/errors.block")
    print(f"Written: analysis/blocks/  ({len(block_sizes)} files)")

    if last_finite_fes:
        _write_fes(anal_dir / "fes.dat", last_finite_fes)
        print(f"Written: analysis/fes.dat  (block size {block_sizes[-1]})")
    else:
        print("Warning: no finite FES bins found; analysis/fes.dat not written.")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_colvar(path: Path, kbt: float):
    """Read COLVAR and return (distances, reweighting_weights) arrays.

    COLVAR columns: time  distance  bias
    Weight = exp((bias - max_bias) / kbt)  — standard reweighting formula.
    """
    data = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("!"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            data.append((float(parts[1]), float(parts[2])))

    if not data:
        raise ValueError("COLVAR file appears empty or has no valid data rows.")

    distances = np.array([d for d, _ in data])
    biases    = np.array([b for _, b in data])
    bmax      = biases.max()
    weights   = np.exp((biases - bmax) / kbt)
    return distances, weights


def _block_fes(distances, weights, gmin, gmax, nbin, kbt, block_size):
    """Compute FES and block-error estimate for a given block size.

    Returns list of (r, fes, error) tuples, one per grid bin.
    Bins with no data have fes=inf, error=inf.
    """
    dx = (gmax - gmin) / float(nbin - 1)

    # bin each frame
    keys = []
    ws   = []
    for r, w in zip(distances, weights):
        k = int(round((r - gmin) / dx))
        if 0 <= k < nbin:
            keys.append(k)
            ws.append(w)

    ndata  = len(keys)
    nblock = ndata // block_size if block_size > 0 else 0

    histo_ave  = {}
    histo_ave2 = {}

    for iblock in range(nblock):
        i0 = iblock * block_size
        i1 = i0 + block_size
        histo = {}
        for i in range(i0, i1):
            k = keys[i]
            histo[k] = histo.get(k, 0.0) + ws[i]
        for k in histo:
            histo[k] /= float(block_size)
        for k, v in histo.items():
            if k in histo_ave:
                histo_ave[k]  += v
                histo_ave2[k] += v * v
            else:
                histo_ave[k]  = v
                histo_ave2[k] = v * v

    result = []
    nb = float(nblock) if nblock > 1 else None

    for i in range(nbin):
        r = gmin + i * dx
        if nb is not None and i in histo_ave:
            aveh = histo_ave[i] / nb
            s2h  = (histo_ave2[i] / nb - aveh * aveh) * nb / (nb - 1.0)
            errh = math.sqrt(max(s2h, 0.0) / nb)
            fes  = -kbt * math.log(aveh)
            errf = kbt / aveh * errh
            result.append((r, fes, errf))
        else:
            result.append((r, math.inf, math.inf))

    return result


def _write_fes(path: Path, fes_data):
    with open(path, "w") as f:
        for r, fes, err in fes_data:
            if math.isfinite(fes):
                f.write(f"{r:.9f}  {fes:.9f}  {err:.9f}\n")
