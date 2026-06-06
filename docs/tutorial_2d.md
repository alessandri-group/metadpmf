# Tutorial: 2D free-energy surface

This tutorial extends the [standard 1D workflow](tutorial.md) to a **2D
free-energy surface**, adding an orientational second collective variable
(CV2) on top of the centre-of-mass distance.

CV2 is the cosine of the angle between the two ring-plane normals — it
distinguishes stacked/parallel arrangements (`cosθ ≈ ±1`) from T-shaped ones
(`cosθ ≈ 0`). `metadpmf` builds the (somewhat convoluted) PLUMED for it
automatically from the ring atoms; you only choose the atoms and the mode.

---

## Three modes

The metadynamics bias always acts on the distance. CV2 is added in one of two
ways, selected by `enabled` + `bias`:

| Mode | `cv2.enabled` | `cv2.bias` | Bias acts on | Use when |
|---|---|---|---|---|
| **A** | `false` | — | distance | 1D PMF (the [standard tutorial](tutorial.md)) |
| **B** — project | `true` | `false` | distance | CV2 explores freely on its own; you just want to *see* it |
| **C** — bias both | `true` | `true` | distance + CV2 | CV2 does not sample well unbiased (common for atomistic systems) |

The analysis (reweighting, block FES, Jacobian, plot) is identical for B and C;
only the PLUMED inputs differ. This tutorial covers **B and C** — the only
config difference between them is the `bias` flag.

---

## 1. Enable CV2 in `config.yaml`

```yaml
cv2:
  enabled: true
  type:    costheta      # built-in cosθ generator
  bias:    false         # false → Mode B (project); true → Mode C (bias both)
  sigma:   0.1           # metad width on cosθ (only used in Mode C)
  label:   "cos θ"
  min:    -1.0
  max:     1.0
  bins:    51
  # mol1_plane: [1, 2, 3]   # 3 atoms defining each ring plane;
  # mol2_plane: [4, 5, 6]   # default to molecule.mol*_atoms when those have 3
```

For a 3-bead Martini ring the plane atoms are just the molecule's centre atoms,
so `mol1_plane`/`mol2_plane` can be omitted. For atomistic rings, give three
non-collinear atoms of each ring explicitly.

The cosθ normal of each ring is the cross product of two edge vectors built from
those three atoms, so atom **order matters** only in that it must be consistent
between runs (it is, since it comes from the same config).

---

## 2. Run the production simulation

```bash
metadpmf run && bash run_metad.sh
```

- **Mode B:** `plumed.dat` is identical to the 1D case — the bias is on distance
  only. cosθ is reconstructed later from the trajectory.
- **Mode C:** `plumed.dat` includes the generated cosθ block and a 2D
  `METAD ARG=dist,costheta` with per-CV `SIGMA` and grid. cosθ is biased.

---

## 3. Reweight the trajectory

```bash
metadpmf reweight && bash reweight.sh
```

Because `cv2.enabled: true`, metadpmf creates `analysis_2d/` and **generates**
`analysis_2d/plumed_analysis.dat` from your config (no hand-written PLUMED). The
METAD in it mirrors the production run — 1D in Mode B, 2D in Mode C — so the
bias is reconstructed correctly in both cases.

Produces `analysis_2d/COLVAR` with **four** columns: `time  dist  cv2  bias`.

---

## 4. Block FES analysis

```bash
metadpmf fes
```

Reads the four-column `analysis_2d/COLVAR`, computes reweighting weights
`exp((bias − bias_max)/k_BT)`, and runs the block analysis over the 2D grid.

Output in `analysis_2d/`:

| File | Contents |
|---|---|
| `dist.cv2.weight` | Per-frame (distance, CV2, reweighting weight) |
| `errors.block` | (block size, mean 2D FES error) |
| `fes.dat` | 2D FES from the largest block size (4 columns: cv1, cv2, fes, error) |
| `blocks/fes.{N}.dat` | 2D FES for every block size N |

`fes.dat` follows the gnuplot convention: rows grouped by CV2 value, blank lines
between groups, CV1 varying within each group.

---

## 5. Apply Jacobian correction and plot

```bash
metadpmf pmf
```

The Jacobian correction is applied on the **CV1 (distance) axis only**, because
the geometric sampling bias arises from the spherical volume element in $r$:

$$\Delta G(r, \text{CV2}) = \text{FES}(r, \text{CV2}) + 2RT\ln(r)$$

The corrected surface is shifted to zero by subtracting the mean over
`pmf.shift_range` (default 1.5–1.7 nm) averaged across **all CV2 bins**.

Output in `analysis_2d/`:

| File | Contents |
|---|---|
| `fes_2d_corr.dat` | Jacobian-corrected and shifted 2D FES (cv1, cv2, ΔG) |
| `pmf_2d.pdf` | Filled contour plot with coolwarm colormap |

Contour levels and axis limits are set in `config.yaml`:

```yaml
pmf:
  levels_2d: [-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10]
  xrange: [0.2, 1.7]
  yrange: [-1.0, 1.0]
```

---

## Also getting the 1D PMF from a 2D-biased run

A 2D-biased run (Mode C) still contains the full 1D PMF: reweighting removes the
*entire* 2D bias, so the weighted frames can be histogrammed along distance
alone (marginalising over cosθ). Add `--1d` to `fes` and `pmf`:

```bash
metadpmf fes --1d      # → analysis_2d/fes_1d.dat (+ dist_1d.weight, blocks_1d/)
metadpmf pmf --1d      # → analysis_2d/pmf_1d.pdf, fes_1d_corr_and_shifted.dat
```

This reuses the existing `analysis_2d/COLVAR` (no re-reweighting needed) and
writes alongside the 2D outputs, so you get both the 2D FES and the 1D PMF from
a single simulation. The flag is a no-op error on a 1D run — there, plain
`metadpmf fes` / `pmf` already give the 1D PMF.

---

## cosθ sign degeneracy (`fold`)

`cosθ = +1` and `cosθ = -1` both mean the ring planes are parallel — the sign
only reflects which face each normal points from. For a **planar** molecule the
two faces are equivalent (this holds even for in-plane-asymmetric molecules such
as phenol), so the two halves of the cosθ axis are physically redundant.

- Keep `fold: none` (default) for projection (Mode B): a converged FES should be
  symmetric about 0, which is a useful sanity check.
- Set `fold: abs` (→ |cosθ|) or `fold: sq` (→ cos²θ) when **biasing** a planar
  molecule (Mode C) to avoid filling two equivalent regions; also set `min: 0.0`.
- Do **not** fold a molecule whose two faces genuinely differ (non-planar, or
  different substituents above/below the ring): there `+1` and `-1` are distinct
  states.

---

## Advanced: a custom CV2 (`type: custom`)

If you want a different second CV (not cosθ), set `type: custom` and write the
analysis PLUMED yourself:

```yaml
cv2:
  enabled: true
  type:    custom
  plumed_analysis: plumed_analysis_2d.dat
```

`metadpmf reweight` then copies your file into `analysis_2d/plumed_analysis.dat`
instead of generating it. The file must define your CV2 and print the four
columns in order:

```
PRINT STRIDE=1 ARG=dist,<your_cv2>,metad.bias FILE=COLVAR
```

with `STRIDE=1` and a METAD block matching your production run (same
`BIASFACTOR`, `TEMP`, `SIGMA`, `FILE=HILLS`). This is the escape hatch; for ring
orientation the built-in `type: costheta` is the easy path.

---

## Directory layout after a complete 2D run

```
MOL_in_SOLVENT/
├── config.yaml
├── start.gro
├── system.top
├── plumed.dat
├── metad.mdp
├── run_metad.sh
├── metad.tpr / metad.xtc / metad.log
├── HILLS
├── COLVAR
├── reweight.sh
└── analysis_2d/
    ├── plumed_analysis.dat   ← generated (type: costheta) or copied (type: custom)
    ├── HILLS
    ├── COLVAR
    ├── dist.cv2.weight
    ├── errors.block
    ├── fes.dat
    ├── fes_2d_corr.dat
    ├── pmf_2d.pdf
    ├── fes_1d.dat            ← only if you ran `fes --1d`
    ├── pmf_1d.pdf            ← only if you ran `pmf --1d`
    └── blocks/
        ├── fes.1.dat
        ├── fes.11.dat
        └── ...
```
