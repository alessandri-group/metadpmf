"""Tests for metadpmf.config — no GROMACS/PLUMED required."""

from metadpmf.config import _apply_defaults, kbt, R, render_plumed, render_plumed_analysis


def _cfg(**overrides):
    """Return a config with defaults applied, with optional overrides merged in."""
    base = {"molecule": {"mol1_atoms": [1, 2, 3], "mol2_atoms": [4, 5, 6]}}
    base.update(overrides)
    return _apply_defaults(base)


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

def test_temperature_default():
    cfg = _apply_defaults({})
    assert cfg["temperature"] == 298.15


def test_fes_defaults():
    cfg = _apply_defaults({})
    assert cfg["fes"]["min"] == 0.2
    assert cfg["fes"]["max"] == 1.7
    assert cfg["fes"]["bins"] == 51
    assert cfg["fes"]["block_max"] == 1000


def test_pmf_defaults():
    cfg = _apply_defaults({})
    assert cfg["pmf"]["shift_range"] == [1.5, 1.7]
    assert cfg["pmf"]["xrange"] is None
    assert cfg["pmf"]["yrange"] is None


def test_cv2_defaults():
    cfg = _apply_defaults({})
    assert cfg["cv2"]["enabled"] is False
    assert cfg["cv2"]["min"] == -1.0
    assert cfg["cv2"]["max"] == 1.0
    assert cfg["cv2"]["bins"] == 51
    assert cfg["cv2"]["label"] == "cos θ"
    assert cfg["cv2"]["plumed_analysis"] == "plumed_analysis_2d.dat"


def test_backend_defaults():
    cfg = _apply_defaults({})
    assert cfg["backend"] == "local"
    assert cfg["backends"]["local"]["gmx"] == "gmx"
    assert cfg["backends"]["local"]["plumed"] == "plumed"


# ---------------------------------------------------------------------------
# kbt
# ---------------------------------------------------------------------------

def test_kbt_at_300():
    cfg = _apply_defaults({"temperature": 300.0})
    assert abs(kbt(cfg) - R * 300.0) < 1e-12


def test_kbt_at_default():
    cfg = _apply_defaults({})
    assert abs(kbt(cfg) - R * 298.15) < 1e-12


# ---------------------------------------------------------------------------
# render_plumed
# ---------------------------------------------------------------------------

def test_render_plumed_contains_wholemolecules():
    text = render_plumed(_cfg())
    assert "WHOLEMOLECULES" in text
    assert "ENTITY0=1,2,3" in text
    assert "ENTITY1=4,5,6" in text


def test_render_plumed_contains_metad():
    text = render_plumed(_cfg())
    assert "METAD ARG=dist" in text


def test_render_plumed_contains_print():
    text = render_plumed(_cfg())
    assert "PRINT" in text
    assert "dist,metad.bias" in text


def test_render_plumed_analysis_restart():
    text = render_plumed_analysis(_cfg())
    assert "RESTART" in text
    assert "PACE=1000000000" in text
    assert "STRIDE=1" in text
