# Tutorial: 2D free-energy surface

This tutorial extends the [standard 1D workflow](tutorial.md) to produce a
**2D free-energy surface** by projecting the trajectory onto a second,
unbiased collective variable (CV2).

The metadynamics bias still acts on the centre-of-mass distance $r$ (CV1)
only. CV2 is a structural order parameter that is tracked during the
simulation but not biased — for example, the cosine of the angle between the
normal vectors of two ring molecules, which captures their relative
orientation. Because the reweighting weights depend only on the bias, the 2D
FES is obtained by binning frames simultaneously in both CV1 and CV2.

---

## Prerequisites

Complete steps 0–3 of the [standard tutorial](tutorial.md): run the
metadynamics simulation and produce `metad.xtc` and `HILLS`. The production
run is identical; the only differences come in the analysis steps.

---

## 1. Enable CV2 in config.yaml

Add the `cv2` block to your `config.yaml`:

```yaml
cv2:
  enabled: true
  label:   "cos θ"        # y-axis label in the 2D PMF plot
  min:    -1.0            # CV2 grid lower bound
  max:     1.0            # CV2 grid upper bound
  bins:    51             # number of CV2 histogram bins
  plumed_analysis: plumed_analysis_2d.dat
```

`label`, `min`, `max`, and `bins` describe your CV2. For the ring-orientation
example below, $\cos\theta$ ranges from −1 (perpendicular planes) to +1
(parallel/anti-parallel planes).

---

## 2. Write plumed_analysis_2d.dat

Unlike the 1D case — where `metadpmf reweight` generates the PLUMED input
automatically — the 2D case requires you to provide the PLUMED file yourself,
because the definition of CV2 is system-specific.

Place `plumed_analysis_2d.dat` in your simulation directory (next to
`config.yaml`). Below is a template for two 3-bead ring molecules where CV2
is $\cos\theta$, the cosine of the angle between the two plane normals. It
uses only standard PLUMED actions (`DISTANCE COMPONENTS` + `CUSTOM`) and does
not require any optional modules.

The plane normal for each molecule is computed as the cross product of two
edge vectors: for molecule 1, $\mathbf{n}_1 = \overrightarrow{12} \times
\overrightarrow{23}$, and similarly for molecule 2.

```plumed
RESTART

WHOLEMOLECULES ENTITY0=1,2,3 ENTITY1=4,5,6

c1:   CENTER ATOMS=1,2,3
c2:   CENTER ATOMS=4,5,6
dist: DISTANCE ATOMS=c1,c2

# edge vectors within each ring
d12: DISTANCE ATOMS=1,2 COMPONENTS
d23: DISTANCE ATOMS=2,3 COMPONENTS
d45: DISTANCE ATOMS=4,5 COMPONENTS
d56: DISTANCE ATOMS=5,6 COMPONENTS

# normal to ring 1: n1 = d12 x d23
n1x: CUSTOM ARG=d12.y,d12.z,d23.y,d23.z VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO
n1y: CUSTOM ARG=d12.z,d12.x,d23.z,d23.x VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO
n1z: CUSTOM ARG=d12.x,d12.y,d23.x,d23.y VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO

# normal to ring 2: n2 = d45 x d56
n2x: CUSTOM ARG=d45.y,d45.z,d56.y,d56.z VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO
n2y: CUSTOM ARG=d45.z,d45.x,d56.z,d56.x VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO
n2z: CUSTOM ARG=d45.x,d45.y,d56.x,d56.y VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO

# cos(theta) = (n1 . n2) / (|n1| |n2|),  range [-1, 1]
costheta: CUSTOM ARG=n1x,n1y,n1z,n2x,n2y,n2z VAR=a,b,c,d,f,g FUNC=(a*d+b*f+c*g)/sqrt((a^2+b^2+c^2)*(d^2+f^2+g^2)) PERIODIC=NO

uwall: UPPER_WALLS ARG=dist AT=2.0 KAPPA=200.0

metad: METAD ARG=dist PACE=1000000000 HEIGHT=1.0 SIGMA=0.05 FILE=HILLS BIASFACTOR=5 GRID_MIN=0.1 GRID_MAX=3.0 TEMP=298.15 GRID_BIN=300

PRINT STRIDE=1 ARG=dist,costheta,metad.bias FILE=COLVAR
```

Key requirements for the `PRINT` line:

- Column order must be `dist, <cv2>, metad.bias` (bias last).
- `STRIDE=1` so every frame is written.
- `FILE=COLVAR` (the name `metadpmf fes` expects).

Adjust atom indices, METAD parameters, and `UPPER_WALLS` to match your
production `plumed.dat` exactly.

!!! tip "Verify against your plumed.dat"
    The `METAD` block must use the same `BIASFACTOR`, `TEMP`, `SIGMA`, and
    `FILE=HILLS` as the production run, otherwise the reconstructed bias will
    be incorrect.

---

## 3. Reweight the trajectory

```bash
metadpmf reweight
```

Because `cv2.enabled: true`, metadpmf creates `analysis_2d/` and copies your
`plumed_analysis_2d.dat` into `analysis_2d/plumed_analysis.dat`.

```bash
bash reweight.sh
```

Produces `analysis_2d/COLVAR` with **four** columns: `time  dist  cv2  bias`.

---

## 4. Block FES analysis

```bash
metadpmf fes
```

Reads the four-column `analysis_2d/COLVAR`, computes reweighting weights, and
runs the block analysis over the 2D grid simultaneously in CV1 and CV2.

Output in `analysis_2d/`:

| File | Contents |
|---|---|
| `dist.cv2.weight` | Per-frame (distance, CV2, reweighting weight) |
| `errors.block` | (block size, mean 2D FES error) |
| `fes.dat` | 2D FES from the largest block size (4 columns: cv1, cv2, fes, error) |
| `blocks/fes.{N}.dat` | 2D FES for every block size N |

`fes.dat` follows the gnuplot convention: rows are grouped by CV2 value,
with blank lines between groups and CV1 varying within each group.

---

## 5. Apply Jacobian correction and plot

```bash
metadpmf pmf
```

The Jacobian correction is applied on the **CV1 (distance) axis only**,
because the geometric sampling bias arises from the spherical volume element
in $r$:

$$\Delta G(r, \text{CV2}) = \text{FES}(r, \text{CV2}) + 2RT\ln(r)$$

The corrected surface is then shifted to zero by subtracting the mean value
over `pmf.shift_range` (default 1.5–1.7 nm) averaged across **all CV2 bins**.

Output in `analysis_2d/`:

| File | Contents |
|---|---|
| `fes_2d_corr.dat` | Jacobian-corrected and shifted 2D FES (cv1, cv2, ΔG) |
| `pmf_2d.pdf` | Filled contour plot with coolwarm colormap |

The contour levels can be adjusted in `config.yaml`:

```yaml
pmf:
  levels_2d: [-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10]
  xrange: [0.2, 1.7]   # optional axis limits
  yrange: [-1.0, 1.0]
```

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
├── plumed_analysis_2d.dat    ← written by you
├── reweight.sh
└── analysis_2d/
    ├── plumed_analysis.dat   ← copied from plumed_analysis_2d.dat
    ├── HILLS
    ├── COLVAR
    ├── dist.cv2.weight
    ├── errors.block
    ├── fes.dat
    ├── fes_2d_corr.dat
    ├── pmf_2d.pdf
    └── blocks/
        ├── fes.1.dat
        ├── fes.11.dat
        └── ...
```
