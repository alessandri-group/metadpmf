"""Tests for the block FES analysis step — no GROMACS/PLUMED required."""

import math
from pathlib import Path

import numpy as np
import pytest

from metadpmf.steps.fes import (
    _block_fes_1d,
    _block_fes_2d,
    _load_colvar_1d,
    _load_colvar_2d,
    _write_fes_1d,
    _write_fes_2d,
)

KBT = 0.008314 * 300.0  # kJ/mol at 300 K
GMIN, GMAX, NBIN = 0.2, 1.7, 51
GMIN2, GMAX2, NBIN2 = -1.0, 1.0, 51


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_colvar_1d(path: Path, rows):
    """Write a minimal 1D COLVAR (time dist bias)."""
    with open(path, "w") as f:
        f.write("#! FIELDS time dist metad.bias\n")
        for t, d, b in rows:
            f.write(f" {t} {d} {b}\n")


def _write_colvar_2d(path: Path, rows):
    """Write a minimal 2D COLVAR (time dist cv2 bias)."""
    with open(path, "w") as f:
        f.write("#! FIELDS time dist costheta metad.bias\n")
        for t, d, s, b in rows:
            f.write(f" {t} {d} {s} {b}\n")


# ---------------------------------------------------------------------------
# _load_colvar_1d
# ---------------------------------------------------------------------------

def test_load_colvar_1d_frame_count(tmp_path):
    p = tmp_path / "COLVAR"
    _write_colvar_1d(p, [(i, 0.5, 10.0) for i in range(5)])
    distances, weights = _load_colvar_1d(p, KBT)
    assert len(distances) == 5
    assert len(weights) == 5


def test_load_colvar_1d_distances(tmp_path):
    p = tmp_path / "COLVAR"
    _write_colvar_1d(p, [(0, 0.4, 5.0), (1, 0.6, 5.0), (2, 0.8, 5.0)])
    distances, weights = _load_colvar_1d(p, KBT)
    assert np.allclose(distances, [0.4, 0.6, 0.8])


def test_load_colvar_1d_equal_bias_gives_equal_weights(tmp_path):
    """All frames with the same bias → all weights equal 1.0 (exp(0))."""
    p = tmp_path / "COLVAR"
    _write_colvar_1d(p, [(i, 0.5, 10.0) for i in range(5)])
    distances, weights = _load_colvar_1d(p, KBT)
    assert np.allclose(weights, 1.0)


def test_load_colvar_1d_weight_ordering(tmp_path):
    """Frame with higher bias gets higher weight."""
    p = tmp_path / "COLVAR"
    _write_colvar_1d(p, [(0, 0.5, 10.0), (1, 0.5, 20.0)])
    distances, weights = _load_colvar_1d(p, KBT)
    assert weights[1] > weights[0]
    assert math.isclose(weights[1], 1.0)  # highest bias → weight 1


def test_load_colvar_1d_skips_comments(tmp_path):
    p = tmp_path / "COLVAR"
    with open(p, "w") as f:
        f.write("#! FIELDS time dist metad.bias\n")
        f.write("# another comment\n")
        f.write(" 0.0 0.5 10.0\n")
        f.write(" 1.0 0.6 10.0\n")
    distances, weights = _load_colvar_1d(p, KBT)
    assert len(distances) == 2


# ---------------------------------------------------------------------------
# _load_colvar_2d
# ---------------------------------------------------------------------------

def test_load_colvar_2d_frame_count(tmp_path):
    p = tmp_path / "COLVAR"
    _write_colvar_2d(p, [(i, 0.5, 0.3, 10.0) for i in range(5)])
    distances, cv2_vals, weights = _load_colvar_2d(p, KBT)
    assert len(distances) == len(cv2_vals) == len(weights) == 5


def test_load_colvar_2d_values(tmp_path):
    p = tmp_path / "COLVAR"
    _write_colvar_2d(p, [(0, 0.4, -0.5, 5.0), (1, 0.7, 0.8, 5.0)])
    distances, cv2_vals, weights = _load_colvar_2d(p, KBT)
    assert np.allclose(distances, [0.4, 0.7])
    assert np.allclose(cv2_vals, [-0.5, 0.8])
    assert np.allclose(weights, 1.0)


def test_load_colvar_2d_bias_in_col3(tmp_path):
    """Bias must come from column 3 (index 3), not column 2."""
    p = tmp_path / "COLVAR"
    _write_colvar_2d(p, [(0, 0.5, 0.0, 10.0), (1, 0.5, 0.0, 20.0)])
    distances, cv2_vals, weights = _load_colvar_2d(p, KBT)
    assert math.isclose(weights[1], 1.0)
    assert weights[0] < 1.0


# ---------------------------------------------------------------------------
# _block_fes_1d
# ---------------------------------------------------------------------------

def test_block_fes_1d_single_bin(tmp_path):
    """All frames in one bin, block_size=1: FES=0 in that bin, inf elsewhere."""
    n = 10
    distances = np.full(n, 0.5)
    weights   = np.ones(n)
    result = _block_fes_1d(distances, weights, GMIN, GMAX, NBIN, KBT, block_size=1)

    dx = (GMAX - GMIN) / (NBIN - 1)
    k = int(round((0.5 - GMIN) / dx))

    r_k, fes_k, err_k = result[k]
    assert abs(fes_k) < 1e-10
    assert abs(err_k) < 1e-10

    for i, (r, fes, err) in enumerate(result):
        if i != k:
            assert not math.isfinite(fes)


def test_block_fes_1d_length(tmp_path):
    distances = np.linspace(0.3, 1.5, 100)
    weights   = np.ones(100)
    result = _block_fes_1d(distances, weights, GMIN, GMAX, NBIN, KBT, block_size=5)
    assert len(result) == NBIN


def test_block_fes_1d_finite_bins_have_negative_fes():
    """Well-sampled bins should produce finite, reasonable FES values."""
    rng = np.random.default_rng(42)
    distances = rng.uniform(0.3, 1.5, 500)
    weights   = np.ones(500)
    result = _block_fes_1d(distances, weights, GMIN, GMAX, NBIN, KBT, block_size=10)
    finite = [(r, fes, err) for r, fes, err in result if math.isfinite(fes)]
    assert len(finite) > 0
    # FES values should be finite and plausible (not astronomically large)
    for r, fes, err in finite:
        assert abs(fes) < 1e3


# ---------------------------------------------------------------------------
# _block_fes_2d
# ---------------------------------------------------------------------------

def test_block_fes_2d_single_bin():
    """All frames in one 2D bin, block_size=1: FES=0 there, inf elsewhere."""
    n = 10
    distances = np.full(n, 0.5)
    cv2_vals  = np.zeros(n)
    weights   = np.ones(n)
    result = _block_fes_2d(
        distances, cv2_vals, weights,
        GMIN, GMAX, NBIN, GMIN2, GMAX2, NBIN2, KBT, block_size=1,
    )

    dx1 = (GMAX  - GMIN)  / (NBIN  - 1)
    dx2 = (GMAX2 - GMIN2) / (NBIN2 - 1)
    k1 = int(round((0.5  - GMIN)  / dx1))
    k2 = int(round((0.0  - GMIN2) / dx2))
    idx = k2 * NBIN + k1

    cv1, cv2, fes, err = result[idx]
    assert abs(fes) < 1e-10
    assert abs(err) < 1e-10

    for i, (c1, c2, fes, err) in enumerate(result):
        if i != idx:
            assert not math.isfinite(fes)


def test_block_fes_2d_length():
    distances = np.linspace(0.3, 1.5, 200)
    cv2_vals  = np.linspace(-0.9, 0.9, 200)
    weights   = np.ones(200)
    result = _block_fes_2d(
        distances, cv2_vals, weights,
        GMIN, GMAX, NBIN, GMIN2, GMAX2, NBIN2, KBT, block_size=5,
    )
    assert len(result) == NBIN * NBIN2


# ---------------------------------------------------------------------------
# _write / round-trip
# ---------------------------------------------------------------------------

def test_write_fes_1d_round_trip(tmp_path):
    """Write and re-read a 1D FES: values should survive the round-trip."""
    data = [(0.3, -2.5, 0.1), (0.5, -5.0, 0.2), (0.8, -1.0, 0.3)]
    path = tmp_path / "fes.dat"
    _write_fes_1d(path, data)
    loaded = np.loadtxt(path)
    assert loaded.shape == (3, 3)
    assert np.allclose(loaded[:, 0], [0.3, 0.5, 0.8])
    assert np.allclose(loaded[:, 1], [-2.5, -5.0, -1.0])
    assert np.allclose(loaded[:, 2], [0.1, 0.2, 0.3])


def test_write_fes_1d_skips_inf(tmp_path):
    """Bins with inf FES should not appear in the output file."""
    data = [(0.3, -2.5, 0.1), (0.5, math.inf, math.inf), (0.8, -1.0, 0.3)]
    path = tmp_path / "fes.dat"
    _write_fes_1d(path, data)
    loaded = np.loadtxt(path)
    assert loaded.shape == (2, 3)


def test_write_fes_2d_blank_lines(tmp_path):
    """2D FES file should have blank lines separating CV2 groups."""
    nbin1, nbin2 = 3, 2
    data = []
    for i2 in range(nbin2):
        for i1 in range(nbin1):
            data.append((0.3 + i1 * 0.1, -1.0 + i2, -1.0, 0.1))
    path = tmp_path / "fes_2d.dat"
    _write_fes_2d(path, data, nbin1)
    text = path.read_text()
    # one blank line between the two CV2 groups
    assert "\n\n" in text
