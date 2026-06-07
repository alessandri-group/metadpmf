"""Tests for metadpmf.config — no GROMACS/PLUMED required."""

import pytest

from metadpmf.config import _apply_defaults, _validate, kbt, R, render_plumed, render_plumed_analysis, render_mdp


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
    assert cfg["cv2"]["type"] == "costheta"
    assert cfg["cv2"]["bias"] is False
    assert cfg["cv2"]["fold"] == "none"
    assert cfg["cv2"]["sigma"] == 0.1
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


# ---------------------------------------------------------------------------
# Validation — CV range consistency
# ---------------------------------------------------------------------------

def _valid_cfg(**overrides):
    """Minimal valid config (molecule + forcefield required for _validate to pass)."""
    base = {
        "molecule": {"mol1_atoms": [1, 2], "mol2_atoms": [3, 4]},
        "forcefield": "martini",
    }
    base.update(overrides)
    return _apply_defaults(base)


def test_validate_ok():
    """Default config should pass validation without errors."""
    _validate(_valid_cfg())


def test_validate_fes_max_exceeds_wall():
    cfg = _valid_cfg()
    cfg["fes"]["max"] = 2.5          # > wall_at default of 2.0
    with pytest.raises(ValueError, match="fes.max"):
        _validate(cfg)


def test_validate_shift_range_above_fes_max():
    cfg = _valid_cfg()
    cfg["pmf"]["shift_range"] = [1.8, 2.0]  # 2.0 > fes.max default of 1.7
    with pytest.raises(ValueError, match="shift_range"):
        _validate(cfg)


def test_validate_shift_range_below_fes_min():
    cfg = _valid_cfg()
    cfg["pmf"]["shift_range"] = [0.1, 0.15]  # 0.1 < fes.min default of 0.2
    with pytest.raises(ValueError, match="shift_range"):
        _validate(cfg)


def test_validate_shift_range_inverted():
    cfg = _valid_cfg()
    cfg["pmf"]["shift_range"] = [1.6, 1.5]
    with pytest.raises(ValueError, match="shift_range"):
        _validate(cfg)


def test_validate_bad_forcefield():
    cfg = _valid_cfg()
    cfg["forcefield"] = "amber"
    with pytest.raises(ValueError, match="not known"):
        _validate(cfg)


def test_validate_missing_forcefield():
    cfg = _valid_cfg()
    del cfg["forcefield"]
    with pytest.raises(ValueError, match="forcefield is required"):
        _validate(cfg)


# ---------------------------------------------------------------------------
# Force field selection — dt default and MDP template
# ---------------------------------------------------------------------------

def test_forcefield_has_no_default():
    """forcefield is required; _apply_defaults must not invent one."""
    assert "forcefield" not in _apply_defaults({})


def test_dt_default_martini():
    assert _apply_defaults({"forcefield": "martini"})["mdp"]["dt"] == 0.020


def test_dt_default_opls():
    assert _apply_defaults({"forcefield": "opls"})["mdp"]["dt"] == 0.002


def test_dt_default_gaff2():
    assert _apply_defaults({"forcefield": "gaff2"})["mdp"]["dt"] == 0.002


def test_dt_unset_when_forcefield_missing():
    """No forcefield → no dt default (validation will reject the config)."""
    assert "dt" not in _apply_defaults({})["mdp"]


def test_dt_explicit_override_wins():
    cfg = _apply_defaults({"forcefield": "opls", "mdp": {"dt": 0.001}})
    assert cfg["mdp"]["dt"] == 0.001


def test_render_mdp_martini_reaction_field():
    text = render_mdp(_cfg(forcefield="martini"))
    assert "reaction-field" in text
    assert "dt                       = 0.02" in text


def test_render_mdp_opls_pme():
    text = render_mdp(_cfg(forcefield="opls"))
    assert "coulombtype              = PME" in text
    assert "constraints              = h-bonds" in text
    assert "dt                       = 0.002" in text


def test_render_mdp_gaff2_pme():
    text = render_mdp(_cfg(forcefield="gaff2"))
    assert "coulombtype              = PME" in text
    assert "dt                       = 0.002" in text


def test_render_mdp_substitutes_temperature():
    text = render_mdp(_cfg(forcefield="martini", temperature=310.0))
    assert "310.0" in text
    assert "TEMP" not in text


# ---------------------------------------------------------------------------
# Second CV (cosθ) — three modes
# ---------------------------------------------------------------------------

def _cv2_cfg(**cv2):
    base = {"enabled": True}
    base.update(cv2)
    return _cfg(forcefield="martini", cv2=base)


def test_cv2_plane_defaults():
    """Plane atoms default to the molecule's centre atoms iff there are 3."""
    assert _cfg()["cv2"]["mol1_plane"] == [1, 2, 3]   # 3-bead molecule
    four = _apply_defaults({"molecule": {"mol1_atoms": [1, 2, 3, 4], "mol2_atoms": [5, 6, 7, 8]}})
    assert four["cv2"]["mol1_plane"] is None           # not 3 → no default


def test_mode_b_production_stays_1d():
    """Projection mode must not touch the production run (bias on distance only)."""
    text = render_plumed(_cv2_cfg(bias=False))
    assert "METAD ARG=dist " in text
    assert "costheta" not in text


def test_mode_b_analysis_projects_costheta_unbiased():
    text = render_plumed_analysis(_cv2_cfg(bias=False))
    assert "costheta: CUSTOM" in text
    assert "d12: DISTANCE ATOMS=1,2 COMPONENTS" in text
    assert "METAD ARG=dist " in text                   # 1D bias
    assert "ARG=dist,costheta,metad.bias" in text       # 4-column COLVAR


def test_mode_c_production_biases_two_cvs():
    """Order matters: distance first, cv2 second, in every 2-valued field."""
    text = render_plumed(_cv2_cfg(bias=True))
    assert "costheta: CUSTOM" in text
    assert "METAD ARG=dist,costheta" in text
    assert "SIGMA=0.05,0.1" in text
    assert "GRID_MIN=0.1,-1.0" in text
    assert "GRID_MAX=3.0,1.0" in text
    assert "GRID_BIN=300,51" in text


def test_mode_c_analysis_biases_two_cvs():
    text = render_plumed_analysis(_cv2_cfg(bias=True))
    assert "METAD ARG=dist,costheta" in text
    assert "ARG=dist,costheta,metad.bias" in text


def test_cv2_fold_changes_biased_variable():
    abs_text = render_plumed(_cv2_cfg(bias=True, fold="abs", min=0.0))
    assert "cv2f: CUSTOM ARG=costheta VAR=a FUNC=abs(a)" in abs_text
    assert "METAD ARG=dist,cv2f" in abs_text
    sq_text = render_plumed(_cv2_cfg(bias=True, fold="sq", min=0.0))
    assert "cv2f: CUSTOM ARG=costheta VAR=a FUNC=a*a" in sq_text


def test_cv2_custom_type_skips_generator():
    """type: custom must not inject the cosθ block (user supplies the file)."""
    text = render_plumed_analysis(_cv2_cfg(type="custom"))
    assert "costheta" not in text
    assert "ARG=dist,metad.bias" in text


def test_validate_cv2_bad_type():
    with pytest.raises(ValueError, match="cv2.type"):
        _validate(_cv2_cfg(type="foo"))


def test_validate_cv2_plane_missing():
    """The main footgun: type=costheta with a molecule that isn't 3 atoms."""
    cfg = _apply_defaults({
        "molecule": {"mol1_atoms": [1, 2, 3, 4], "mol2_atoms": [5, 6, 7, 8]},
        "forcefield": "martini",
        "cv2": {"enabled": True},
    })
    with pytest.raises(ValueError, match="mol1_plane"):
        _validate(cfg)
