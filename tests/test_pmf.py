"""Tests for the PMF correction step — no GROMACS/PLUMED required.

Regression tests use real FES data from a 2 µs benzene-in-water simulation
(martini300_metad_2us) as ground truth. Both 1D and 2D corrections are
verified to floating-point precision (~1e-9 kJ/mol).
"""

from pathlib import Path

import numpy as np
import pytest

from metadpmf.config import _apply_defaults
from metadpmf.steps.pmf import _parse_fes_1d, _parse_fes_2d, _run_1d, _run_2d

DATA  = Path(__file__).parent / "data"
T     = 300.0  # temperature used for the reference simulation


def _cfg_1d(tmp_path):
    cfg = _apply_defaults({"temperature": T})
    cfg["_dir"] = tmp_path
    return cfg


def _cfg_2d(tmp_path):
    cfg = _apply_defaults({"temperature": T, "cv2": {"enabled": True}})
    cfg["_dir"] = tmp_path
    return cfg


# ---------------------------------------------------------------------------
# _parse_fes_1d
# ---------------------------------------------------------------------------

def test_parse_fes_1d_length():
    distances, fes, errors = _parse_fes_1d(DATA / "fes_1d.dat")
    assert len(distances) == len(fes) == len(errors)
    assert len(distances) > 0


def test_parse_fes_1d_cv1_range():
    distances, fes, errors = _parse_fes_1d(DATA / "fes_1d.dat")
    assert min(distances) >= 0.2
    assert max(distances) <= 1.7


def test_parse_fes_1d_finite_values():
    distances, fes, errors = _parse_fes_1d(DATA / "fes_1d.dat")
    assert all(np.isfinite(fes))
    assert all(np.isfinite(errors))


# ---------------------------------------------------------------------------
# _parse_fes_2d
# ---------------------------------------------------------------------------

def test_parse_fes_2d_shape():
    FES, ERR, cv1_edges, cv2_edges = _parse_fes_2d(DATA / "fes_2d.dat", 51, 51)
    assert FES.shape == (51, 51)
    assert ERR.shape == (51, 51)
    assert len(cv1_edges) == 51
    assert len(cv2_edges) == 51


def test_parse_fes_2d_cv1_range():
    FES, ERR, cv1_edges, cv2_edges = _parse_fes_2d(DATA / "fes_2d.dat", 51, 51)
    assert abs(cv1_edges[0]  - 0.2)  < 1e-6
    assert abs(cv1_edges[-1] - 1.7)  < 1e-6


def test_parse_fes_2d_cv2_range():
    FES, ERR, cv1_edges, cv2_edges = _parse_fes_2d(DATA / "fes_2d.dat", 51, 51)
    assert abs(cv2_edges[0]  - (-1.0)) < 1e-6
    assert abs(cv2_edges[-1] - 1.0)   < 1e-6


def test_parse_fes_2d_has_finite_bins():
    FES, ERR, cv1_edges, cv2_edges = _parse_fes_2d(DATA / "fes_2d.dat", 51, 51)
    assert np.isfinite(FES).any()


# ---------------------------------------------------------------------------
# 1D PMF — regression against ground truth
# ---------------------------------------------------------------------------

def test_pmf_1d_regression(tmp_path):
    """Jacobian correction + shift must reproduce stored output to 1e-6 kJ/mol."""
    cfg      = _cfg_1d(tmp_path)
    anal_dir = tmp_path / "analysis"
    anal_dir.mkdir()
    (anal_dir / "fes.dat").write_text((DATA / "fes_1d.dat").read_text())

    _run_1d(cfg, T, anal_dir, anal_dir / "fes.dat")

    result       = np.loadtxt(anal_dir / "fes_corr_and_shifted.dat")
    ground_truth = np.loadtxt(DATA / "fes_1d_corr.dat")

    assert result.shape == ground_truth.shape
    assert np.allclose(result, ground_truth, atol=1e-6)


def test_pmf_1d_well_depth(tmp_path):
    """Well depth should be ~-3.13 kJ/mol at r ≈ 0.5 nm."""
    cfg      = _cfg_1d(tmp_path)
    anal_dir = tmp_path / "analysis"
    anal_dir.mkdir()
    (anal_dir / "fes.dat").write_text((DATA / "fes_1d.dat").read_text())

    _run_1d(cfg, T, anal_dir, anal_dir / "fes.dat")

    data       = np.loadtxt(anal_dir / "fes_corr_and_shifted.dat")
    well_depth = data[:, 1].min()
    r_min      = data[data[:, 1].argmin(), 0]

    assert abs(well_depth - (-3.1339)) < 0.01
    assert abs(r_min - 0.5) < 0.05


def test_pmf_1d_zero_at_shift_range(tmp_path):
    """Mean of corrected FES over shift_range [1.5, 1.7] should be ~0."""
    cfg      = _cfg_1d(tmp_path)
    anal_dir = tmp_path / "analysis"
    anal_dir.mkdir()
    (anal_dir / "fes.dat").write_text((DATA / "fes_1d.dat").read_text())

    _run_1d(cfg, T, anal_dir, anal_dir / "fes.dat")

    data = np.loadtxt(anal_dir / "fes_corr_and_shifted.dat")
    r, pmf = data[:, 0], data[:, 1]
    plateau = pmf[(r >= 1.5) & (r <= 1.7)]
    assert abs(plateau.mean()) < 0.1


def test_pmf_1d_creates_plot(tmp_path):
    cfg      = _cfg_1d(tmp_path)
    anal_dir = tmp_path / "analysis"
    anal_dir.mkdir()
    (anal_dir / "fes.dat").write_text((DATA / "fes_1d.dat").read_text())

    _run_1d(cfg, T, anal_dir, anal_dir / "fes.dat")

    assert (anal_dir / "pmf.pdf").exists()


# ---------------------------------------------------------------------------
# 2D PMF — regression against ground truth
# ---------------------------------------------------------------------------

def _read_corr_2d(path: Path) -> dict:
    """Read a 2D corrected FES file into a {(cv1, cv2): value} dict."""
    result = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3 or parts[2] == "Infinity":
                continue
            result[(float(parts[0]), float(parts[1]))] = float(parts[2])
    return result


def test_pmf_2d_regression(tmp_path):
    """2D Jacobian correction + shift must reproduce stored output to 1e-6 kJ/mol."""
    cfg      = _cfg_2d(tmp_path)
    anal_dir = tmp_path / "analysis_2d"
    anal_dir.mkdir()
    (anal_dir / "fes.dat").write_text((DATA / "fes_2d.dat").read_text())

    _run_2d(cfg, T, anal_dir, anal_dir / "fes.dat")

    computed = _read_corr_2d(anal_dir / "fes_2d_corr.dat")
    expected = _read_corr_2d(DATA / "fes_2d_corr.dat")

    assert len(computed) == len(expected)
    for key, exp_val in expected.items():
        assert key in computed, f"missing bin {key}"
        assert abs(computed[key] - exp_val) < 1e-6, (
            f"bin {key}: got {computed[key]:.6f}, expected {exp_val:.6f}"
        )


def test_pmf_2d_creates_plot(tmp_path):
    cfg      = _cfg_2d(tmp_path)
    anal_dir = tmp_path / "analysis_2d"
    anal_dir.mkdir()
    (anal_dir / "fes.dat").write_text((DATA / "fes_2d.dat").read_text())

    _run_2d(cfg, T, anal_dir, anal_dir / "fes.dat")

    assert (anal_dir / "pmf_2d.pdf").exists()
