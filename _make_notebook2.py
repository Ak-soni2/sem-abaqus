"""Generate SEM_TO_ABAQUS_MULTI_ABRASIVE.ipynb.

A second, separate notebook. ``SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb`` and
everything it embeds is left exactly as it was verified: this file writes a new
notebook that carries the same ``semgrit`` package plus ``semgrit_multi``, and
adds the two things the single-abrasive notebook cannot do.

1. The chip thickness is **computed for every element** by sweeping the grit
   trajectories, instead of being handed to the subroutine as four constants
   that describe one wedge. Any number of abrasives, and a later grit correctly
   sees the groove an earlier one left.
2. A **second, purely local criterion** in ``vumat_grind2.for``:
   ``W_p L_c >= PSI Kc^2/E``, which needs no chip thickness, no coordinates and
   no kinematics at all, and triggers on history so a point turns brittle as
   the cut deepens under it.

Both are gated. ``verify_envelope.py`` proves the swept field reproduces the
closed-form single-grit wedge it generalises, and ``verify_vumat_grind2.py``
proves the new subroutine is still bit-identical to ``vumat_grind.for`` wherever
the new criterion is off.

    python _make_notebook2.py
"""
import base64
import glob
import gzip
import io
import json
import os
import tarfile

# Everything the notebook needs to run and to check itself. semgrit is included
# unmodified -- semgrit_multi builds on it rather than replacing it.
FILES = (sorted(glob.glob('semgrit/*.py'))
         + sorted(glob.glob('semgrit_multi/*.py'))
         + ['verify_rigid_deck.py',
            'verify_rigid_deck2.py',
            'verify_pipeline_A.py',
            'verify_hybrid_deck.py',
            'verify_vumat_grind.py',
            'verify_vumat_grind2.py',
            'verify_envelope.py',
            '_derive_grind2.py',
            'vumat_grind.for',
            'vumat_grind2.for',
            'vumat_jh2.for',
            '_hybrid_test/driver.f',
            '_hybrid_test/vaba_param.inc'])


def payload() -> str:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode='w') as tf:
        for f in FILES:
            ti = tarfile.TarInfo(f.replace(os.sep, '/'))
            data = open(f, 'rb').read()
            ti.size, ti.mtime, ti.mode = len(data), 0, 0o644
            tf.addfile(ti, io.BytesIO(data))
    return base64.b64encode(gzip.compress(buf.getvalue(), 9)).decode()


def _lines(text):
    ls = text.strip('\n').split('\n')
    return [l + '\n' for l in ls[:-1]] + [ls[-1]]


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": _lines(text)}


def code(text):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": _lines(text)}


CELLS = []

CELLS.append(md(r"""
# SEM image -> multi-abrasive grinding wheel, ductile and brittle

The companion notebook, `SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb`, builds a
single-abrasive deck whose chip thickness is described by four constants in the
material card: `h(u) = H0 + HG u - u^2/2R`. That is exact for one grit and
cannot be anything else, because one wedge is one grit.

This notebook removes that limit two different ways, and you can use either or
both.

## 1 - the chip thickness, computed instead of prescribed

The wheel is one discrete rigid body driven by a prescribed velocity boundary
condition. It cannot deflect, slow down or be pushed back, so **every grit's
position at every instant of the step is known in closed form before Abaqus
starts**. Sweeping that motion gives the undeformed chip thickness for every
element of the workpiece:

```
for each grit, in the order it crosses the block:
    d_i(u,z) = how deep its surface reaches at that station
    h_i      = max(0, d_i - D)        <- what is left for it to remove
    material between D and d_i is removed by grit i, and its chip
    thickness is h_i
    D <- max(D, d_i)
```

The running `D` is what makes several grits correct rather than merely possible:
the second grit over a station only removes what the first one left. For one
grit `D` starts at zero and the whole thing collapses to the closed-form wedge
-- **which is checked**: `verify_envelope.py` compares the swept field against
that wedge station by station and requires agreement to the sweep's own stated
depth resolution. It currently agrees to **0.2 nm over 132 stations**.

The field is written into the deck as `*Initial Conditions, type=FIELD`, and
`vumat_grind.for` reads field variable 1 when `PROPS(56) = 1`. **So the verified
subroutine does not change at all.**

## 2 - a criterion that needs no geometry whatsoever

`vumat_grind2.for` adds a second switch that a material point can evaluate
entirely from its own history:

$$W_p \, L_c \;\ge\; \Psi \frac{K_c^2}{E} \qquad\Rightarrow\qquad \text{brittle}$$

`W_p` is the accumulated plastic work per unit volume and `L_c` the element's
own characteristic length, which Abaqus hands the VUMAT for free. No
coordinates, no kinematics, no field variable, no grit count -- and it works
unchanged for a second pass over the same groove, or for a wheel whose motion is
not prescribed at all.

It comes from the same Griffith balance `dc` does, but **the exponents do not
match and this notebook does not pretend they do**: the pointwise balance gives
`dc = Psi (H/E)^1 (Kc/H)^2`, where the two published geometric forms use `+0.5`
and `-1`. So `PSI` is defaulted to the value that makes the local criterion trip
at exactly the `dc` the deck already chose,

$$\Psi = \frac{d_c E H}{K_c^2} \qquad\Longrightarrow\qquad W_p L_c \ge H d_c$$

which reads as plainly as it should: brittle once the plastic work per unit area
exceeds the cost of plastically removing a layer of thickness `dc` at flow
stress `H`.

Unlike the geometric switch this one **triggers on history**, so a point starts
ductile and turns brittle as the cut deepens under it. That is the physical
transition rather than a line drawn in advance.

> It is regularised by `L_c`, so it is mesh-dependent by construction, as every
> energy-based failure criterion is. Halving the element halves the work density
> needed to trigger. **`PSI` is therefore calibrated for a mesh** -- quote the
> element size with it.

## What runs what

| | |
|---|---|
| `SWMODE = 0` | geometric only. `vumat_grind2.for` is then bit-identical to `vumat_grind.for`, which is checked on every history. |
| `SWMODE = 1` | energy only. The chip thickness is never consulted; you need no field at all. |
| `SWMODE = 2` | both: brittle if either says so. |

**Nothing in the single-abrasive notebook, in `semgrit/`, or in
`vumat_grind.for` is modified by any of this.** Run the cells top to bottom.
"""))

CELLS.append(code('''
#@title 1 - Setup: unpack the pipeline (run once) { display-mode: "form" }
# semgrit, semgrit_multi, both VUMATs and all four gates are embedded below, so
# this notebook is self-contained.
import base64, gzip, io, os, subprocess, sys, tarfile

PAYLOAD = (
__PAYLOAD__
)

WORK = "/content" if os.path.isdir("/content") else os.getcwd()
os.chdir(WORK)
with tarfile.open(fileobj=io.BytesIO(gzip.decompress(base64.b64decode(PAYLOAD)))) as tf:
    tf.extractall(WORK)
if WORK not in sys.path:
    sys.path.insert(0, WORK)

missing = []
for mod, pip in [("numpy", "numpy"), ("scipy", "scipy"),
                 ("skimage", "scikit-image"),
                 ("cv2", "opencv-python-headless"), ("shapely", "shapely"),
                 ("PIL", "pillow"), ("mapbox_earcut", "mapbox-earcut"),
                 ("matplotlib", "matplotlib")]:
    try:
        __import__(mod)
    except ImportError:
        missing.append(pip)
if missing:
    print("installing:", " ".join(missing))
    subprocess.run([sys.executable, "-m", "pip", "-q", "install", *missing],
                   check=True)

import semgrit.build_deck as _bd
import semgrit_multi as _sm            # noqa: F401
print("pipeline ready in", WORK)
print("modules   : %d in semgrit/, %d in semgrit_multi/"
      % (len([f for f in os.listdir("semgrit") if f.endswith(".py")]),
         len([f for f in os.listdir("semgrit_multi") if f.endswith(".py")])))
print("subroutines: vumat_grind.for, vumat_grind2.for, vumat_jh2.for")

# The energy-criterion gate compiles Fortran. Colab has gcc but not always
# gfortran; without it the gates still run and say which part they skipped.
HAVE_FC = False
try:
    from verify_vumat_grind import find_gfortran
    print("gfortran  :", find_gfortran())
    HAVE_FC = True
except SystemExit:
    print("gfortran  : not found. The deck gates still run; the parts that")
    print("            compile the subroutine will say they were skipped.")
    print("            To enable them:  !apt-get -qq install gfortran")


def need(names, where):
    """Stop with the cell to run, instead of a NameError on an unfamiliar name."""
    absent = [n for n in names.split() if n not in globals()]
    if absent:
        raise SystemExit("run %s first - this cell needs %s"
                         % (where, ", ".join(absent)))
'''))

CELLS.append(code('''
#@title 2 - Where are your SEM images? { display-mode: "form" }
SOURCE = "upload"  #@param ["upload", "google drive", "already on disk"]
IMAGE_PATH = "/content/drive/MyDrive/sem/*.tif"  #@param {type:"string"}

import glob, os
IMAGES = []
if SOURCE == "upload":
    from google.colab import files
    for name in files.upload():
        IMAGES.append(os.path.abspath(name))
else:
    if SOURCE == "google drive":
        from google.colab import drive
        drive.mount("/content/drive", force_remount=False)
    pat = IMAGE_PATH
    if os.path.isdir(pat):
        pat = os.path.join(pat, "*.tif")
    IMAGES = sorted(glob.glob(pat))

if not IMAGES:
    raise SystemExit("no images found - check SOURCE / IMAGE_PATH")
print("%d image(s):" % len(IMAGES))
for p in IMAGES:
    print("  ", p, "  %.1f MB" % (os.path.getsize(p) / 1e6))
'''))

CELLS.append(code('''
#@title 3 - Measure the grains { display-mode: "form" }
#@markdown Identical to the companion notebook: the pixel size comes from the
#@markdown Zeiss metadata, the segmentation settings are the tuned ones, and the
#@markdown library is cached so re-running a wheel change costs nothing.
import os
from semgrit.quick import SIMPLE_MEASURE, library_summary, measure_images

need("IMAGES", "cell 2")
OUT_MEAS = os.path.join(WORK, "1_measurements")
MEASURED = measure_images(IMAGES, OUT_MEAS, pixel_size_um=0.0, **SIMPLE_MEASURE)
SOLIDS, ALL_GRAINS = MEASURED["solids"], MEASURED["grains"]
print()
LIB = library_summary(SOLIDS)
'''))

CELLS.append(code('''
#@title 4 - Settings { display-mode: "form" }
MA_NAME = "multi_abrasive"       #@param {type:"string"}

#@markdown ### The wheel and how many abrasives
MA_DIAMETER_MM = 50.0            #@param {type:"number"}
MA_SLICE_MM = 2.0                #@param {type:"number"}
MA_GRIT_MODE = "count"           #@param ["count", "single", "areal_density", "concentration"]
MA_GRIT_COUNT = 12               #@param {type:"integer"}
MA_AREAL_DENSITY = 5000.0        #@param {type:"number"}
MA_CONCENTRATION = 100.0         #@param {type:"number"}
MA_ARC_WINDOW_MM = 0.1           #@param {type:"number"}
#@markdown &nbsp;&nbsp;Dress only this much of the arc, centred. `0` = the whole
#@markdown slice. This is the knob that decides **how many grits actually cross
#@markdown the block**: spread 12 grits over 2 mm and they are 167 um apart, so
#@markdown only one reaches a 48 um block. Put them in a 0.1 mm window and they
#@markdown are 8 um apart, and ten of them cross.
MA_FACE_WINDOW_MM = 0.0          #@param {type:"number"}
MA_SEED = 20260731               #@param {type:"integer"}

#@markdown ### The workpiece
MA_WP_LENGTH_MM = 0.048          #@param {type:"number"}
MA_WP_WIDTH_MM = 0.015           #@param {type:"number"}
MA_WP_DEPTH_MM = 0.006           #@param {type:"number"}
MA_ELEMENT_UM = 0.30             #@param {type:"number"}
MA_ELEMENT_AXIAL_UM = 0.0        #@param {type:"number"}
MA_ELEMENT_DEPTH_UM = 0.0        #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = use the base size. Refining only the **depth** is
#@markdown what lets a nanometre-scale cut be resolved at all, and it is nearly
#@markdown free: the stable increment follows the smallest dimension, so keep
#@markdown the axial size coarse to pay for it.
MA_SURFACE_LAYER_UM = 0.0        #@param {type:"number"}
MA_DEPTH_GROWTH = 1.3            #@param {type:"number"}
#@markdown &nbsp;&nbsp;A **graded** depth mesh: `MA_ELEMENT_DEPTH_UM` elements for
#@markdown the first `MA_SURFACE_LAYER_UM`, then growing by `MA_DEPTH_GROWTH` into
#@markdown the body. `0` = ungraded. This is what makes a 0.03 um surface layer
#@markdown affordable -- the increment was already going to follow the smallest
#@markdown element, and grading stops the coarse body below costing 150 layers.
MA_PROTRUSION_STD = 0.12         #@param {type:"number"}
#@markdown &nbsp;&nbsp;Spread of grit protrusion. **Small = well dressed**, so many
#@markdown grits stand at nearly one height and several cut at a shallow infeed.
#@markdown At the default 0.12 only the tallest grit reaches the work at a
#@markdown sub-micron depth of cut and the multi-abrasive case collapses to the
#@markdown single-abrasive one.
MA_STANDOFF_UM = 0.0             #@param {type:"number"}
MA_DEPTH_OF_CUT_UM = 0.0         #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = automatic (85% of the grain protrusion). **This is
#@markdown the knob that decides whether anything is ductile at all.** For a
#@markdown rock with dc of a few nanometres, a micron-scale infeed is entirely
#@markdown brittle -- which is the real physics, not a modelling artefact.
MA_SURFACE_SPEED_M_S = 30.0      #@param {type:"number"}

#@markdown ### The workpiece material
MA_MATERIAL = "sandstone"  #@param ["sandstone", "silicon_carbide"]
#@markdown &nbsp;&nbsp;Picks the **whole** card together: the 17 JH-2 constants, the
#@markdown density, the ductile Johnson-Cook constants, and the hardness and
#@markdown toughness `dc` comes from. `silicon_carbide` is the **SiC-N** card
#@markdown (rho 3163, G 183 GPa, HEL 14.457 GPa, K1 204.785 GPa, A 0.96,
#@markdown B 0.35, N 0.65, T 0.37 GPa, D1 = D2 = 0.48), supplied labelled
#@markdown "monocrystalline silicon" but actually the published silicon
#@markdown **carbide** numbers. It runs about 7x slower than sandstone on the
#@markdown same mesh: its wave speed is 1.19e7 mm/s against 1.76e6, so the
#@markdown stable increment is 7x smaller.
MA_OVERRIDE_MATERIAL = False  #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;Off: every Johnson-Cook / SGE / damage / hardness number
#@markdown below is taken from the material and the fields are ignored. On:
#@markdown hand-enter them. Cell 5 prints which source it used.

#@markdown ### The transition
MA_DC_NM = 0.0                   #@param {type:"number"}
MA_DC_FORM = "Bifano: lambda_c (E/H) (Kc/H)^2"  #@param ["Bifano: lambda_c (E/H) (Kc/H)^2", "lambda_c (H/E)^0.5 (Kc/H)^2"]
#@markdown &nbsp;&nbsp;Bifano is the default: it is the form with a published
#@markdown calibrated `lambda_c` and it is what the `RUN_ME*` packages use. The
#@markdown other form gives 5.3 nm on sandstone and 0.7 nm on SiC, below what any
#@markdown mesh resolves, so everything comes out brittle.
MA_LAMBDA_C = 0.15               #@param {type:"number"}
MA_HARDNESS_MPA = 1000.0         #@param {type:"number"}
MA_KIC_MPA_SQRT_M = 0.30         #@param {type:"number"}

#@markdown ### Which criterion
MA_SWMODE = "geometric: h vs dc"  #@param ["geometric: h vs dc", "energy: W_p L_c vs PSI Kc^2/E", "both"]
MA_PSI = 0.0                     #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = derive it from `dc` so the two criteria agree.
#@markdown Remember it is calibrated **for a mesh**.

#@markdown ### Ductile branch - Johnson-Cook (PLACEHOLDERS, calibrate them)
MA_JC_A_MPA = 90.0               #@param {type:"number"}
MA_JC_B_MPA = 50.0               #@param {type:"number"}
MA_JC_N = 0.50                   #@param {type:"number"}
MA_JC_C = 0.020                  #@param {type:"number"}
MA_JC_M = 1.0                    #@param {type:"number"}
MA_E_MPA = 6500.0                #@param {type:"number"}
MA_NU = 0.21                     #@param {type:"number"}
MA_DENSITY_KG_M3 = 2350.0        #@param {type:"number"}
MA_CP_J_KGK = 800.0              #@param {type:"number"}
MA_TMELT_K = 1473.15             #@param {type:"number"}
MA_BURGERS_NM = 0.50             #@param {type:"number"}
MA_TAYLOR_M = 3.0                #@param {type:"number"}
MA_ALPHA = 0.30                  #@param {type:"number"}
MA_LAMBDA_SGE = 1.0              #@param {type:"number"}
MA_R_PRIME = 2.0                 #@param {type:"number"}
MA_D1 = 0.0                      #@param {type:"number"}
MA_D2 = 0.15                     #@param {type:"number"}
MA_D3 = -1.5                     #@param {type:"number"}
MA_D4 = 0.0                      #@param {type:"number"}
MA_D5 = 0.0                      #@param {type:"number"}
MA_DCRIT = 1.0                   #@param {type:"number"}

#@markdown ### To reproduce a clearly visible transition
#@markdown &nbsp;&nbsp;The defaults above are a general-purpose deck and give only
#@markdown a few percent ductile. These are the settings the `RUN_ME*` packages
#@markdown ship, measured by sweeping `ae` until both regions were tens of
#@markdown elements wide:
#@markdown
#@markdown | | sandstone | silicon carbide |
#@markdown |---|---|---|
#@markdown | `MA_DC_FORM` | Bifano | Bifano |
#@markdown | `dc` | 87.75 nm | 52.92 nm |
#@markdown | `MA_DEPTH_OF_CUT_UM` | **0.40** | **0.36** |
#@markdown | `MA_ELEMENT_UM` | 0.30 | 0.30 |
#@markdown | `MA_ELEMENT_AXIAL_UM` | 1.5 | 1.5 |
#@markdown | `MA_ELEMENT_DEPTH_UM` | 0.03 | 0.03 |
#@markdown | `MA_SURFACE_LAYER_UM` | 0.45 | 0.45 |
#@markdown | `MA_DEPTH_GROWTH` | 1.45 | 1.45 |
#@markdown | `MA_PROTRUSION_STD` | 0.015 | 0.015 |
#@markdown | `MA_WP_WIDTH_MM` | 0.009 | 0.009 |
#@markdown | `MA_GRIT_COUNT` | 12 | 12 |
#@markdown | `MA_ARC_WINDOW_MM` | 0.1 | 0.1 |
#@markdown | result | 99 ductile / 136 brittle | 56 ductile / 103 brittle |
#@markdown
#@markdown &nbsp;&nbsp;`ae` is **not** a fixed multiple of `dc` (4.6x here, 6.8x
#@markdown there): the mesh is absolute and `dc` is not, so `ae` also has to be
#@markdown deep enough that enough elements are cut at all. Scaling `ae` down with
#@markdown `dc` for SiC leaves 24 cut elements, all ductile -- a deck that runs and
#@markdown shows one regime. Use the depth-of-cut curve below to retune for a new
#@markdown material rather than assuming a ratio.

#@markdown ### How finely to sweep
MA_DEPTH_RESOLUTION_NM = 0.20    #@param {type:"number"}
#@markdown &nbsp;&nbsp;The blur the sweep allows in the recorded depth. The
#@markdown sample count follows from it, so this is a resolution you state
#@markdown rather than a sample count you guess. Keep it well under `dc`.
MA_FACET_SUBDIVISION = 2         #@param {type:"integer"}

#@markdown ### The depth-of-cut sweep curve
MA_SWEEP_DEPTHS = "0.10,0.20,0.30,0.40,0.60,0.80"  #@param {type:"string"}
#@markdown &nbsp;&nbsp;Depths of cut, in microns, to plot the ductile share
#@markdown against. Each one is a real sweep and costs a few seconds. Leave it
#@markdown empty to skip the curve.
print("settings captured")
'''))

CELLS.append(code('''
#@title 5 - Sweep it and look at it. Nothing is written { display-mode: "form" }
#@markdown Runs the same sweep the build will, so the split shown here is the
#@markdown split the deck gets. Change anything above and re-run.
import matplotlib.pyplot as plt

from semgrit import materials
from semgrit.hybrid import HybridParams, kic_from_mpa_sqrt_m
import dataclasses
from semgrit_multi.build import MultiParams, plan_multi, summary_text
from semgrit_multi.envelope import EnvelopeParams
from semgrit_multi.plot import field_slice_figure, preview_figure

need("SOLIDS", "cell 3")

MA_FORM = 2 if MA_DC_FORM.startswith("Bifano") else 1
print(materials.summary_text(MA_MATERIAL))
print()
if MA_OVERRIDE_MATERIAL:
    print("MA_OVERRIDE_MATERIAL is ON: ductile constants from the cell-4")
    print("fields, not from the material above. The JH-2 card still comes")
    print("from the material.")
    MA_HP = HybridParams(
        enabled=True, h_source=1,
        a_mpa=MA_JC_A_MPA, b_mpa=MA_JC_B_MPA, n=MA_JC_N, c=MA_JC_C, m=MA_JC_M,
        youngs_mpa=MA_E_MPA, poisson=MA_NU, density_kg_m3=MA_DENSITY_KG_M3,
        specific_heat_j_kgk=MA_CP_J_KGK, tmelt_k=MA_TMELT_K,
        burgers_mm=MA_BURGERS_NM * 1e-6, taylor_factor=MA_TAYLOR_M,
        alpha=MA_ALPHA, sge_exponent=MA_LAMBDA_SGE, r_prime=MA_R_PRIME,
        d1=MA_D1, d2=MA_D2, d3=MA_D3, d4=MA_D4, d5=MA_D5, dcrit=MA_DCRIT,
        dc_mm=MA_DC_NM * 1e-6, lambda_c=MA_LAMBDA_C,
        hardness_mpa=MA_HARDNESS_MPA,
        kic=kic_from_mpa_sqrt_m(MA_KIC_MPA_SQRT_M), dc_form=MA_FORM)
    _ref = materials.hybrid_params(MA_MATERIAL)
    _bad = [k for k in ("youngs_mpa", "poisson", "density_kg_m3", "a_mpa")
            if abs(getattr(MA_HP, k) - getattr(_ref, k)) > 1e-9 * max(
                1.0, abs(getattr(_ref, k)))]
    if _bad:
        print()
        print("!! WARNING: %s disagree with %s, so the brittle branch and the"
              % (", ".join(_bad), MA_MATERIAL))
        print("!! ductile branch are describing DIFFERENT materials. That is")
        print("!! legitimate only if you meant it.")
else:
    print("ductile constants taken from the material card above. Set")
    print("MA_OVERRIDE_MATERIAL = True in cell 4 to hand-enter them instead.")
    MA_HP = materials.hybrid_params(
        MA_MATERIAL, h_source=1, dc_form=MA_FORM,
        dc_mm=MA_DC_NM * 1e-6, lambda_c=MA_LAMBDA_C)

MA_PARAMS = MultiParams(
    material=MA_MATERIAL,
    name=MA_NAME, diameter_mm=MA_DIAMETER_MM, arc_length_mm=MA_SLICE_MM,
    grit_mode=MA_GRIT_MODE, grit_count=MA_GRIT_COUNT,
    areal_density_per_mm2=MA_AREAL_DENSITY, concentration=MA_CONCENTRATION,
    grit_arc_window_mm=MA_ARC_WINDOW_MM, grit_width_window_mm=MA_FACE_WINDOW_MM,
    seed=MA_SEED,
    wp_length_mm=MA_WP_LENGTH_MM, wp_width_mm=MA_WP_WIDTH_MM,
    wp_depth_mm=MA_WP_DEPTH_MM, element_um=MA_ELEMENT_UM,
    element_axial_um=MA_ELEMENT_AXIAL_UM,
    element_depth_um=MA_ELEMENT_DEPTH_UM,
    surface_layer_um=MA_SURFACE_LAYER_UM, depth_growth=MA_DEPTH_GROWTH,
    protrusion_std=MA_PROTRUSION_STD,
    standoff_um=MA_STANDOFF_UM, depth_of_cut_um=MA_DEPTH_OF_CUT_UM,
    surface_speed_m_s=MA_SURFACE_SPEED_M_S, hybrid=MA_HP,
    envelope=EnvelopeParams(depth_resolution_mm=MA_DEPTH_RESOLUTION_NM * 1e-6,
                            facet_subdivision=MA_FACET_SUBDIVISION))

MA_PLAN = plan_multi(MA_PARAMS, SOLIDS, log=print)
MA_ENV, MA_DC = MA_PLAN["envelope"], MA_PLAN["dc_mm"]
print()
print("=" * 78)
print(summary_text({"split": MA_PLAN["split"], "envelope": MA_ENV.stats,
                    "dc_nm": MA_DC * 1e6}, MA_WP_LENGTH_MM))
print("=" * 78)
# # transition visuals + measured trajectory (multi abrasive)
# Four pictures, in the order the questions get asked.
from semgrit_multi.plot import trajectory_figure

# 1. Where did the abrasives go, and where does that cross dc? The depth of cut
#    the transition happens at is read straight off the vertical axis.
fig = trajectory_figure(MA_PLAN["plan"]["_place"], MA_PLAN["motion"],
                        MA_PLAN["plan"]["_wp"], MA_DC,
                        step_time_s=MA_PLAN["step_time_s"],
                        paths=MA_PATHS if "MA_PATHS" in dir() else None)
plt.show()
# 2. The groove, the chip thickness against dc, and the map of the ground face.
fig = preview_figure(MA_ENV, MA_DC, MA_PLAN["plan"]["_wp"], title=MA_NAME)
plt.show()
# 3. The field through the depth, with dc as a contour.
fig2 = field_slice_figure(MA_ENV, MA_DC, MA_PLAN["plan"]["_wp"])
plt.show()

# 4. How the split moves with the depth of cut. Each point is a real sweep, not
#    an extrapolation, so it costs a few seconds each -- turn it off while you
#    are iterating on something else.
if MA_SWEEP_DEPTHS:
    from semgrit_multi.plot import dc_sweep_figure
    _aes, _fracs = [], []
    for _ae in [float(x) for x in MA_SWEEP_DEPTHS.split(",") if x.strip()]:
        _p = dataclasses.replace(MA_PARAMS, depth_of_cut_um=_ae)
        try:
            _r = plan_multi(_p, SOLIDS, log=lambda *a: None)
        except Exception as _exc:
            print("  ae = %.3f um: %s" % (_ae, str(_exc)[:60]))
            continue
        _aes.append(_ae)
        _fracs.append(_r["split"]["ductile_fraction_of_cut"])
        print("  ae = %5.3f um -> %6s cut, %5.1f%% ductile"
              % (_ae, format(_r["split"]["n_cut"], ","),
                 100 * _r["split"]["ductile_fraction_of_cut"]))
    if _aes:
        fig3 = dc_sweep_figure(_aes, _fracs, dc_nm=MA_DC * 1e6,
                               chosen_um=MA_PLAN["plan"]["depth_of_cut_um"])
        plt.show()

print("Nothing written. When the split is what you want, run the next cell.")
'''))

CELLS.append(code('''
#@title 5b - Optional: replay a MEASURED trajectory { display-mode: "form" }
#@markdown After a real single-grit experiment you have the groove, and the
#@markdown groove *is* the measurement. Give it here and the simulated abrasive
#@markdown is put exactly where the real one went.
#@markdown
#@markdown **What this changes and what it does not.** It replaces the chip
#@markdown thickness the switch reads, so the ductile/brittle split follows your
#@markdown measurement. It does **not** make Abaqus drive the wheel along that
#@markdown path -- the deck still turns the wheel with its velocity boundary
#@markdown condition. Driving the wheel along a measured path replaces the
#@markdown rotation-plus-infeed BC with a prescribed displacement and changes
#@markdown the mechanics of the run, so it is deliberately a separate decision;
#@markdown `semgrit_multi.trajectory.deck_amplitudes` writes those tables if you
#@markdown want them.
MA_USE_MEASURED_PATH = False   #@param {type:"boolean"}
MA_PATH_SOURCE = "csv"         #@param ["csv", "image"]

#@markdown ### From a table of coordinates
MA_PATH_CSV = "/content/groove.csv"   #@param {type:"string"}
MA_PATH_COLUMNS = "u,depth"    #@param ["u,depth", "u,z,depth", "t,u,z,depth", "auto"]
MA_PATH_SCALE = "um"           #@param ["mm", "um", "nm"]
#@markdown &nbsp;&nbsp;The units your table is in. Getting this wrong is the
#@markdown single most likely mistake: a profile in microns read as millimetres
#@markdown lands a thousand times too deep, and the loader will refuse it.
MA_PATH_DEPTH_SIGN = "groove is positive"  #@param ["groove is positive", "groove is negative"]
#@markdown &nbsp;&nbsp;Most profilometers report a groove as negative.

#@markdown ### Or traced from a scaled image
MA_PATH_IMAGE = "/content/groove.png"  #@param {type:"string"}
MA_PATH_MM_PER_PX_X = 0.001    #@param {type:"number"}
MA_PATH_MM_PER_PX_Y = 0.0001   #@param {type:"number"}
MA_PATH_DARK_IS_MATERIAL = True  #@param {type:"boolean"}
MA_PATH_SURFACE_ROW = -1       #@param {type:"integer"}
#@markdown &nbsp;&nbsp;Pixel row of the uncut surface; `-1` takes the shallowest
#@markdown traced row, which assumes the trace starts outside the groove.

MA_PATH_GRIT = 0               #@param {type:"integer"}
#@markdown &nbsp;&nbsp;Which placed grit follows the path. The others keep their
#@markdown ideal arcs, so one measured scratch can be mixed with modelled ones.

import matplotlib.pyplot as plt
from semgrit_multi.plot import trajectory_check_figure
from semgrit_multi.trajectory import from_csv, from_points, from_profile_image

MA_PATHS = None
if MA_USE_MEASURED_PATH:
    need("MA_PLAN", "cell 5")
    _wp = MA_PLAN["plan"]["_wp"]
    _scale = {"mm": 1.0, "um": 1e-3, "nm": 1e-6}[MA_PATH_SCALE]
    _sign = -1.0 if MA_PATH_DEPTH_SIGN.endswith("negative") else 1.0
    _overlay = None
    if MA_PATH_SOURCE == "csv":
        MA_TRAJ = from_csv(MA_PATH_CSV, columns=MA_PATH_COLUMNS,
                           scale_mm=_scale, depth_sign=_sign)
    else:
        MA_TRAJ, _overlay = from_profile_image(
            MA_PATH_IMAGE, mm_per_px_x=MA_PATH_MM_PER_PX_X,
            mm_per_px_y=MA_PATH_MM_PER_PX_Y,
            dark_is_material=MA_PATH_DARK_IS_MATERIAL,
            surface_row=(None if MA_PATH_SURFACE_ROW < 0
                         else MA_PATH_SURFACE_ROW))
    print("trajectory read from %s" % MA_TRAJ.source)
    for _k, _v in MA_TRAJ.summary().items():
        print("  %-12s %s" % (_k, _v))
    for _n in MA_TRAJ.notes:
        print("  note: %s" % _n)

    # LOOK AT THIS before believing it. A units mistake is obvious here and
    # invisible three cells later.
    fig = trajectory_check_figure(MA_TRAJ, _wp, MA_DC, overlay=_overlay)
    plt.show()

    MA_TRAJ = MA_TRAJ.clipped_to_block(_wp).retimed(
        0.0, MA_PLAN["step_time_s"])
    MA_PATHS = {int(MA_PATH_GRIT): MA_TRAJ.samples}

    # Re-plan with the measured path in place, so cell 5's pictures and the
    # build below both describe the same thing.
    MA_PLAN = plan_multi(MA_PARAMS, SOLIDS, paths=MA_PATHS, log=print)
    MA_ENV, MA_DC = MA_PLAN["envelope"], MA_PLAN["dc_mm"]
    print()
    print(summary_text({"split": MA_PLAN["split"],
                        "envelope": MA_ENV.stats, "dc_nm": MA_DC * 1e6},
                       MA_WP_LENGTH_MM))
    print()
    print("Re-run cell 5 to redraw with the measured path, then build.")
else:
    print("using the ideal wheel kinematics; tick MA_USE_MEASURED_PATH to "
          "replay a measurement instead")
'''))

CELLS.append(code('''
#@title 6 - Build it, inject the field, verify it { display-mode: "form" }
#@markdown Writes the deck, sweeps the field into it, and runs three gates: the
#@markdown two deck verifiers, `verify_hybrid_deck.py`, and `verify_envelope.py`.
#@markdown The subroutine that goes with your `SWMODE` is copied alongside.
MA_BUILD = True                  #@param {type:"boolean"}
MA_DOWNLOAD = True               #@param {type:"boolean"}

import os, shutil, subprocess, sys, time
from semgrit_multi.build import build_multi
from semgrit.quick import bundle, verify_decks

_SW = {"geometric: h vs dc": 0, "energy: W_p L_c vs PSI Kc^2/E": 1,
       "both": 2}[MA_SWMODE]

if MA_BUILD:
    need("MA_PARAMS SOLIDS", "cell 5")
    MA_OUT = os.path.join(WORK, "2_multi_abrasive")
    _t0 = time.time()
    MA_INFO = build_multi(MA_PARAMS, SOLIDS, MA_OUT,
                          paths=(MA_PATHS if "MA_PATHS" in dir() else None))
    print("     built in %.0f s" % (time.time() - _t0))

    # Which subroutine to run, and the two extra constants it needs.
    if _SW == 0:
        _for = "vumat_grind.for"
        print()
        print("SWMODE 0: the geometric criterion only, so vumat_grind.for is")
        print("enough. vumat_grind2.for with SWMODE=0 is bit-identical to it.")
    else:
        _for = "vumat_grind2.for"
        print()
        print("SWMODE %d needs vumat_grind2.for and TWO EXTRA CONSTANTS on the"
              % _SW)
        print("*User Material line. Append them and change constants=56 to 58:")
        print()
        print("   ..., %g, %g" % (_SW, MA_PSI))
        print()
        print("   PROPS(57) = SWMODE, PROPS(58) = PSI (0 = derive it from dc).")
        print("Also raise *Depvar from 20 to 22 so SDV21 and SDV22 are kept.")
    for _f in (_for, "vumat_jh2.for"):
        if os.path.exists(os.path.join(WORK, _f)):
            shutil.copy(os.path.join(WORK, _f), os.path.join(MA_OUT, _f))

    # The two deck verifiers check GEOMETRY, and they are run on the
    # un-injected deck: it is the same model, and verify_envelope.py proves
    # below that injection adds the field block and changes nothing else. That
    # keeps two Abaqus-validated verifiers untouched instead of teaching them a
    # keyword they have never needed.
    print()
    print("geometry, on the un-injected deck (the field block is additive,")
    print("which verify_envelope.py proves line by line):")
    if not verify_decks(WORK, [MA_INFO["plain_path"] or MA_INFO["path"]]):
        raise SystemExit("the deck did not verify - read the FAIL lines above")

    # And the field itself: that it reproduces the closed-form wedge it
    # generalises, that the deck carries it exactly, and that the compiled
    # subroutine reads it and branches on it as the sweep predicted.
    print()
    print("#" * 78)
    print("# verify_envelope.py - the field, the injection, and the subroutine")
    print("#" * 78)
    _r = subprocess.run([sys.executable,
                         os.path.join(WORK, "verify_envelope.py"),
                         "--library",
                         os.path.join(OUT_MEAS, "grain_library.pkl")],
                        capture_output=True, text=True, cwd=WORK)
    print(_r.stdout[-6000:])
    if _r.stderr.strip():
        print(_r.stderr[-1500:])
    if _r.returncode != 0:
        raise SystemExit("verify_envelope.py failed - do not run this deck")

    print()
    print("to run it:")
    print("  abaqus job=%s input=%s user=%s double=both cpus=%d"
          % (MA_NAME, os.path.basename(MA_INFO["path"]), _for,
             MA_PARAMS.cores))
    if MA_DOWNLOAD:
        _zip = bundle(WORK, (MA_OUT,), MA_NAME)
        try:
            from google.colab import files
            files.download(_zip)
        except Exception as exc:
            print("(not on Colab - copy the zip yourself)", exc)
else:
    print("not built - tick MA_BUILD when cell 5 shows the split you want")
'''))

CELLS.append(code('''
#@title 7 - Optional: check the subroutines themselves { display-mode: "form" }
#@markdown Compiles both VUMATs and exercises them on a single material point:
#@markdown 112 checks on `vumat_grind.for` and 81 on `vumat_grind2.for`,
#@markdown including bit-identity with `vumat_jh2.for` and with each other.
#@markdown Needs gfortran; on Colab run `!apt-get -qq install gfortran` first.
RUN_SUBROUTINE_GATES = True      #@param {type:"boolean"}

import os, subprocess, sys

if RUN_SUBROUTINE_GATES:
    if not HAVE_FC:
        print("gfortran is not available, so these gates cannot compile")
        print("anything. Run:  !apt-get -qq install gfortran   and re-run")
        print("cell 1, then this cell.")
    else:
        for _g in ("verify_vumat_grind.py", "verify_vumat_grind2.py"):
            print("#" * 78)
            print("# " + _g)
            print("#" * 78)
            _r = subprocess.run([sys.executable, os.path.join(WORK, _g)],
                                capture_output=True, text=True, cwd=WORK)
            print(_r.stdout[-4000:])
            if _r.returncode != 0:
                print(_r.stderr[-1500:])
                raise SystemExit("%s FAILED" % _g)
            print()
else:
    print("subroutine gates skipped")
'''))

CELLS.append(md(r"""
---
## Running it, and reading it

```
abaqus verify -user_explicit
abaqus job=grind input=multi_abrasive_field.inp user=vumat_grind.for double=both cpus=8 interactive
```

`double=both` is not optional: the chip thickness is compared against a
threshold of nanometres on a wheel of 25 mm, a ratio of 1e-7, and single
precision does not have the digits.

For `SWMODE` other than 0, run `vumat_grind2.for` instead and append the two
constants the build cell prints.

### What to plot

| SDV | |
|---|---|
| **13** | the branch: 1 ductile, 2 brittle. Compare it with the map cell 5 drew. |
| **14** | `h` at that point, straight from the field |
| **15** | `dc` |
| **19** | the SGE amplification, 1 = no size effect |
| **21, 22** | plastic work and the energy ratio (`vumat_grind2.for` only). SDV22 reaching 1 is what flips a point. |
| 1, 2, 12 | damage, equivalent plastic strain, STATUS |

The interesting comparison is **SDV13 against the map from cell 5**. Cell 5
predicts, from kinematics alone, which elements should be ductile; SDV13 is what
the subroutine actually ran. They should be the same picture, and if they are not
the field did not reach the material points.

### The experiment this deck is for

Sweep `MA_DEPTH_OF_CUT_UM`. Everything else -- wheel, grains, seating, mesh --
is identical between builds, so any difference in force or chip morphology is
the constitutive law and nothing else. Three runs are worth having:

* deep enough that the whole cut is brittle,
* shallow enough that a real ductile zone appears,
* and `SWMODE = 1`, which reaches the same transition from history rather than
  from geometry. If the two disagree, that disagreement is a result.

### Honest limits

* **The Johnson-Cook constants are placeholders.** `A` is tied to the JH-2
  card's own quasi-static compressive strength so the two branches meet at the
  transition; `B, n, C, m` and `D1..D5` are order-of-magnitude values for a
  quartz-bonded rock. Calibrate them before quoting a force.
* **`lambda_c` belongs to whichever `dc` form you chose.** The two differ by
  `(E/H)^1.5`, about 17x on this rock -- and the energy criterion is a third
  member of the family again, with exponent `+1`.
* **The sweep assumes the wheel's motion is prescribed.** True here, since the
  wheel is a rigid body on a velocity boundary condition. It stops being true if
  the wheel is ever made deformable or its grits are allowed to wear during the
  run -- and that is exactly the case `SWMODE = 1` was added for, because the
  energy criterion needs no motion at all.
* **The sweep ignores elastic deflection of the workpiece.** Deliberately: the
  quantity `dc` is calibrated against is the *undeformed* chip thickness, a
  kinematic quantity by definition.
* **`PSI` is mesh-dependent.** Halving the element halves the work density that
  trips the criterion. Quote the element size with it.
* **The mesh must be able to hold the cut.** Cell 5 prints how many elements sit
  through the deepest cut; below about three, the chip is not resolved and the
  force reads low whichever law is running.
"""))

nb = {
    "nbformat": 4, "nbformat_minor": 0,
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True,
                  "name": "SEM_TO_ABAQUS_MULTI_ABRASIVE.ipynb"},
        "kernelspec": {"name": "python3", "display_name": "Python 3"},
        "language_info": {"name": "python"},
    },
    "cells": CELLS,
}

b64 = payload()
CHUNK = 120
chunks = ['    "%s"\n' % b64[i:i + CHUNK] for i in range(0, len(b64), CHUNK)]
chunks[-1] = chunks[-1].rstrip('\n') + '\n'
spliced = 0
for cell in nb["cells"]:
    if cell["cell_type"] != "code" or "__PAYLOAD__\n" not in cell["source"]:
        continue
    i = cell["source"].index("__PAYLOAD__\n")
    cell["source"][i:i + 1] = chunks
    spliced += 1
assert spliced == 1, 'payload marker found %d times' % spliced

OUT = 'SEM_TO_ABAQUS_MULTI_ABRASIVE.ipynb'
with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
    json.dump(nb, fh, indent=1)

import ast
for n, cell in enumerate(json.load(open(OUT, encoding='utf-8'))["cells"], 1):
    if cell["cell_type"] == "code":
        ast.parse(''.join(cell["source"]), filename='cell %d' % n)

check = json.load(open(OUT, encoding='utf-8'))
print('%s  %.0f KB  %d cells (%d code, %d markdown)'
      % (OUT, os.path.getsize(OUT) / 1024, len(check['cells']),
         sum(c['cell_type'] == 'code' for c in check['cells']),
         sum(c['cell_type'] == 'markdown' for c in check['cells'])))
print('the companion notebook is untouched: %s'
      % ('present' if os.path.exists('SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb')
         else 'MISSING'))
