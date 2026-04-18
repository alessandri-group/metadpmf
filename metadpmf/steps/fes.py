"""Step 3: block free-energy surface analysis.

Reads analysis/COLVAR (produced by plumed driver), computes per-frame
reweighting factors, then performs block analysis to estimate statistical
errors on the FES at a range of block sizes.

For each block size the FES and its standard error are computed via the
block-average method (Flyvbjerg & Petersen, J. Chem. Phys. 91, 461, 1989).
The average error across all FES bins is stored in errors.block so the user
can inspect convergence and select an appropriate block size.

1D mode (cv2.enabled: false):
  COLVAR columns: time  distance  bias
  Output written to analysis/:
    dist.weight          — (distance, reweighting_weight) pairs
    errors.block         — (block_size, mean_error) convergence data
    fes.dat              — FES from largest finite block size; columns: dist, fes, error
    blocks/fes.{N}.dat   — FES for every block size N

2D mode (cv2.enabled: true):
  COLVAR columns: time  distance  cv2  bias
  Output written to analysis/:
    dist.cv2.weight      — (distance, cv2, reweighting_weight) triples
    errors.block         — (block_size, mean_error) convergence data
    fes.dat              — 2D FES from largest finite block size; columns: cv1, cv2, fes, error
    blocks/fes.{N}.dat   — 2D FES for every block size N
"""

import math
from pathlib import Path

import numpy as np

from metadpmf.config import kbt as get_kbt


def run(cfg: dict, block_size: int = None) -> None:
    base      = Path(cfg["_dir"])
    two_d     = cfg["cv2"]["enabled"]
    dir_name  = "analysis_2d" if two_d else "analysis"
    anal_dir  = base / dir_name
    colvar    = anal_dir / "COLVAR"

    if not colvar.exists():
        raise FileNotFoundError(
            f"{dir_name}/COLVAR not found. Run 'metadpmf reweight' first."
        )

    KBT     = get_kbt(cfg)
    fes_cfg = cfg["fes"]
    gmin    = fes_cfg["min"]
    gmax    = fes_cfg["max"]
    nbin    = fes_cfg["bins"]
    bmax_bs = fes_cfg["block_max"]

    if two_d:
        _run_2d(cfg, anal_dir, colvar, KBT, gmin, gmax, nbin, bmax_bs, block_size)
    else:
        _run_1d(anal_dir, colvar, KBT, gmin, gmax, nbin, bmax_bs, block_size)


# ---------------------------------------------------------------------------
# 1D workflow
# ---------------------------------------------------------------------------

def _run_1d(anal_dir, colvar, KBT, gmin, gmax, nbin, bmax_bs, block_size):
    print("Reading analysis/COLVAR ...")
    distances, weights = _load_colvar_1d(colvar, KBT)
    print(f"  {len(distances)} frames loaded")

    dw_path = anal_dir / "dist.weight"
    np.savetxt(dw_path, np.column_stack([distances, weights]), fmt="%.9f")
    print("Written: analysis/dist.weight")

    blocks_dir = anal_dir / "blocks"
    blocks_dir.mkdir(exist_ok=True)

    if block_size is not None:
        fes_data = _block_fes_1d(distances, weights, gmin, gmax, nbin, KBT, block_size)
        _write_fes_1d(anal_dir / "fes.dat", fes_data)
        print(f"Written: analysis/fes.dat  (block size {block_size})")
        return

    block_sizes = list(range(1, bmax_bs, 10))
    errors = []
    last_finite_fes = None

    print(f"Scanning {len(block_sizes)} block sizes (1 to {bmax_bs - 1}, step 10) ...")

    for bs in block_sizes:
        fes_data = _block_fes_1d(distances, weights, gmin, gmax, nbin, KBT, bs)
        _write_fes_1d(blocks_dir / f"fes.{bs}.dat", fes_data)
        finite = [(r, f, e) for r, f, e in fes_data if math.isfinite(f)]
        if finite:
            mean_err = sum(e for _, _, e in finite) / len(finite)
            errors.append((bs, mean_err))
            last_finite_fes = finite

    np.savetxt(anal_dir / "errors.block", errors, fmt="%.9f")
    print("Written: analysis/errors.block")
    print(f"Written: analysis/blocks/  ({len(block_sizes)} files)")

    if last_finite_fes:
        _write_fes_1d(anal_dir / "fes.dat", last_finite_fes)
        print(f"Written: analysis/fes.dat  (block size {block_sizes[-1]})")
    else:
        print("Warning: no finite FES bins found; analysis/fes.dat not written.")


# ---------------------------------------------------------------------------
# 2D workflow
# ---------------------------------------------------------------------------

def _run_2d(cfg, anal_dir, colvar, KBT, gmin1, gmax1, nbin1, bmax_bs, block_size):
    cv2_cfg = cfg["cv2"]
    gmin2   = cv2_cfg["min"]
    gmax2   = cv2_cfg["max"]
    nbin2   = cv2_cfg["bins"]

    print("Reading analysis/COLVAR (2D mode) ...")
    distances, cv2_vals, weights = _load_colvar_2d(colvar, KBT)
    print(f"  {len(distances)} frames loaded")

    dcw_path = anal_dir / "dist.cv2.weight"
    np.savetxt(dcw_path, np.column_stack([distances, cv2_vals, weights]), fmt="%.9f")
    print("Written: analysis/dist.cv2.weight")

    blocks_dir = anal_dir / "blocks"
    blocks_dir.mkdir(exist_ok=True)

    if block_size is not None:
        fes_data = _block_fes_2d(distances, cv2_vals, weights,
                                  gmin1, gmax1, nbin1, gmin2, gmax2, nbin2, KBT, block_size)
        _write_fes_2d(anal_dir / "fes.dat", fes_data, nbin1)
        print(f"Written: analysis/fes.dat  (block size {block_size})")
        return

    block_sizes = list(range(1, bmax_bs, 10))
    errors = []
    last_fes = None

    print(f"Scanning {len(block_sizes)} block sizes (1 to {bmax_bs - 1}, step 10) ...")

    for bs in block_sizes:
        fes_data = _block_fes_2d(distances, cv2_vals, weights,
                                  gmin1, gmax1, nbin1, gmin2, gmax2, nbin2, KBT, bs)
        _write_fes_2d(blocks_dir / f"fes.{bs}.dat", fes_data, nbin1)
        finite = [(cv1, cv2, f, e) for cv1, cv2, f, e in fes_data if math.isfinite(f)]
        if finite:
            mean_err = sum(e for _, _, _, e in finite) / len(finite)
            errors.append((bs, mean_err))
            last_fes = fes_data

    np.savetxt(anal_dir / "errors.block", errors, fmt="%.9f")
    print("Written: analysis/errors.block")
    print(f"Written: analysis/blocks/  ({len(block_sizes)} files)")

    if last_fes is not None:
        _write_fes_2d(anal_dir / "fes.dat", last_fes, nbin1)
        print(f"Written: analysis/fes.dat  (block size {block_sizes[-1]})")
    else:
        print("Warning: no finite FES bins found; analysis/fes.dat not written.")


# ---------------------------------------------------------------------------
# Internal helpers — 1D
# ---------------------------------------------------------------------------

def _load_colvar_1d(path: Path, kbt: float):
    """Read 1D COLVAR and return (distances, weights).

    COLVAR columns: time  distance  bias
    Weight = exp((bias - max_bias) / kbt)
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


def _block_fes_1d(distances, weights, gmin, gmax, nbin, kbt, block_size):
    """Compute 1D FES and block-error estimate for a given block size.

    Returns list of (r, fes, error) tuples, one per grid bin.
    """
    dx = (gmax - gmin) / float(nbin - 1)

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


def _write_fes_1d(path: Path, fes_data):
    with open(path, "w") as f:
        for r, fes, err in fes_data:
            if math.isfinite(fes):
                f.write(f"{r:.9f}  {fes:.9f}  {err:.9f}\n")


# ---------------------------------------------------------------------------
# Internal helpers — 2D
# ---------------------------------------------------------------------------

def _load_colvar_2d(path: Path, kbt: float):
    """Read 2D COLVAR and return (distances, cv2_vals, weights).

    COLVAR columns: time  distance  cv2  bias
    Weight = exp((bias - max_bias) / kbt)
    """
    data = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("!"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            data.append((float(parts[1]), float(parts[2]), float(parts[3])))

    if not data:
        raise ValueError("COLVAR file appears empty or has no valid data rows.")

    distances = np.array([d for d, _, _ in data])
    cv2_vals  = np.array([s for _, s, _ in data])
    biases    = np.array([b for _, _, b in data])
    bmax      = biases.max()
    weights   = np.exp((biases - bmax) / kbt)
    return distances, cv2_vals, weights


def _block_fes_2d(distances, cv2_vals, weights,
                  gmin1, gmax1, nbin1,
                  gmin2, gmax2, nbin2,
                  kbt, block_size):
    """Compute 2D FES and block-error estimate for a given block size.

    Iteration order: CV2 slow (outer), CV1 fast (inner) — matches do_block_fes.py
    and the reading convention in the raw plot script.

    Returns list of (cv1, cv2, fes, error) with length nbin1 * nbin2.
    """
    dx1 = (gmax1 - gmin1) / float(nbin1 - 1)
    dx2 = (gmax2 - gmin2) / float(nbin2 - 1)

    keys = []
    ws   = []
    for r, s, w in zip(distances, cv2_vals, weights):
        k1 = int(round((r - gmin1) / dx1))
        k2 = int(round((s - gmin2) / dx2))
        if 0 <= k1 < nbin1 and 0 <= k2 < nbin2:
            keys.append((k1, k2))
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

    nb = float(nblock) if nblock > 1 else None

    result = []
    for i2 in range(nbin2):
        for i1 in range(nbin1):
            key  = (i1, i2)
            cv1v = gmin1 + i1 * dx1
            cv2v = gmin2 + i2 * dx2
            if nb is not None and key in histo_ave:
                aveh = histo_ave[key] / nb
                s2h  = (histo_ave2[key] / nb - aveh * aveh) * nb / (nb - 1.0)
                errh = math.sqrt(max(s2h, 0.0) / nb)
                fes  = -kbt * math.log(aveh)
                errf = kbt / aveh * errh
                result.append((cv1v, cv2v, fes, errf))
            else:
                result.append((cv1v, cv2v, math.inf, math.inf))

    return result


def _write_fes_2d(path: Path, fes_data, nbin1: int):
    """Write 2D fes.dat with blank lines separating CV2 groups (gnuplot convention)."""
    with open(path, "w") as f:
        for i, (cv1, cv2, fes, err) in enumerate(fes_data):
            if math.isfinite(fes):
                f.write(f"{cv1:.9f}  {cv2:.9f}  {fes:.9f}  {err:.9f}\n")
            else:
                f.write(f"{cv1:.9f}  {cv2:.9f}  Infinity  Infinity\n")
            # blank line after each complete CV1 strip (CV2 group boundary)
            if (i + 1) % nbin1 == 0 and i + 1 < len(fes_data):
                f.write("\n")
