"""Generate SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb, the self-contained Colab pipeline.

The whole semgrit package and both verifiers are embedded as one base64 blob, so the
notebook is the only file a future user needs: no pip install from a repo, no upload
of a zip, no path fixing. Regenerate it with this script whenever the package changes.
"""
import base64
import glob
import gzip
import io
import json
import os
import tarfile

# The subroutines and their gates ship with the notebook too: a deck that
# needs vumat_grind.for is useless without it, and a gate you have to go and
# find is a gate nobody runs.
FILES = (sorted(glob.glob('semgrit/*.py'))
         + sorted(glob.glob('semgrit_multi/*.py'))) + [
    'verify_rigid_deck.py',
    'verify_rigid_deck2.py',
    'verify_pipeline_A.py',
    'verify_hybrid_deck.py',
    'verify_vumat_grind.py',
    'vumat_grind.for',
    'vumat_jh2.for',
    '_hybrid_test/driver.f',
    '_hybrid_test/vaba_param.inc',
]


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
    """Split into nbformat `source` lines.

    Every line except the last must keep its trailing newline: nbformat concatenates
    the list verbatim, so stripping them collapses the whole cell onto one line.
    """
    ls = text.strip('\n').split('\n')
    return [l + '\n' for l in ls[:-1]] + [ls[-1]]


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": _lines(text)}


def code(text):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": _lines(text)}


CELLS = []

CELLS.append(md(r"""
# SEM image → Abaqus grinding wheel

Turns a scanning-electron micrograph of abrasive grit into a **verified Abaqus/Explicit
input deck** of a grinding wheel built from those measured grains.

```
SEM .tif ──▶ calibrate ──▶ segment ──▶ measure ──▶ 3D grain library
                                                        │
                                          ┌─────────────┴─────────────┐
                                          ▼                           ▼
                              wheel .inp + CAE loader          STEP / STL for CAD
                                          │
                                          ▼
                          84 checks (98 when run-ready)
```

## Which cells do I run?

**In a hurry — four cells.**

| | |
|---|---|
| **1** | Setup — run once |
| **2** | Point at your SEM images |
| **3** | Seven choices, then **look** at the model. Nothing is written. Change and re-run freely. |
| **4** | Build, verify, download |

Cell 3 is where you decide. It draws the wheel, the block and every grain, so if the
slice is too long or the workpiece is the wrong size you see it *before* a 25 MB file
is written. The grains are measured once and cached, so re-running cell 3 after a
change takes about a second.

**Full control — the `A` cells.** Same code underneath, every knob exposed.

| | | | |
|---|---|---|---|
| **A1** calibration & segmentation | **A2** measure | **A3** check the measurements | **A4** wheel & grits |
| **A5** workpiece, mesh, outputs | **A6** run-ready analysis | **A7** preview | **A8** abrasive heights & standoff |
| **A9** grinding theory | **A10–A12** 3-D views | **A13** build | **A14** APS (optional) |
| **A15** verify the deck | **A16** download | | |

Skip cells 3 and 4 if you are using the `A` path — set `RUN_SIMPLE` to false in cell 3.

**What comes out**

| file | what it is |
|---|---|
| `<name>.inp` | the Abaqus deck — geometry only, or fully run-ready |
| `<name>_import_into_cae.py` | run this in CAE (**File → Run Script**) to load the deck |
| `<name>_report.json` | every number the build decided, machine-readable |
| `<name>_placements.csv` | where each grit ended up |
| `<name>_postprocess_odb.py` | run after the job: forces, energy balance, material removed |
| `<name>.step` / `.stl` | optional CAD for SOLIDWORKS |
| `<name>_cad.glb` / `<name>_view.glb` | the model as glTF — what the in-notebook CAD viewer shows, and it opens in Blender, Windows 3D Viewer and PowerPoint |
| `*_grains.csv` | 25 measured descriptors per grain |

**Seeing it before you build it.** Cell **A12** is a CAD viewer running in the notebook —
shaded with edges, section planes on any axis, a parts tree, standard views plus one
that looks straight at the dressed face, and click-a-grain to read its protrusion,
size and volume. It draws the deck's own triangles, so it is not a preview of the
model, it *is* the model. No account and no API key; nothing is uploaded.

**The model it builds.** The whole wheel — bond rim **and** every grit — is one
discrete rigid body driven by a single reference node on the axis. The workpiece is
the only deformable part. So you rotate the wheel with **one** boundary condition, and
the bond contributes nothing to the stable time increment.

**Units** are mm, tonne, s, MPa, N throughout; the wheel axis is **Z**.

**Two output modes.** Leave `RUN_READY` off and you get geometry only, to finish in CAE.
Turn it on and the deck carries its own step, boundary conditions, contact, JH-2 material,
section controls, restart and output — **submit it straight from the terminal, no CAE at
all**:

```
abaqus job=grind input=<name>.inp user=vumat_jh2.for double=both cpus=8 interactive
```

> Run the cells top to bottom. Each settings cell is a form: change the boxes, don't
> edit code.
"""))

CELLS.append(code('''
#@title 🔧 1 · Setup — unpack the pipeline (run once) { display-mode: "form" }
# The entire semgrit package and both verifiers are embedded below, so this notebook
# is self-contained: nothing is downloaded and no repository has to still exist.
import base64, gzip, io, os, subprocess, sys, tarfile, textwrap

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
# mapbox_earcut is NOT preinstalled on Colab. Without it every grain fails to
# triangulate and the library comes out empty, which looks like a segmentation
# problem but is not one -- so it is listed here with the rest.
for mod, pip in [("numpy", "numpy"), ("scipy", "scipy"), ("skimage", "scikit-image"),
                 ("cv2", "opencv-python-headless"), ("shapely", "shapely"),
                 ("PIL", "pillow"), ("mapbox_earcut", "mapbox-earcut"),
                 ("matplotlib", "matplotlib"), ("plotly", "plotly"),
                 ("requests", "requests")]:
    try:
        __import__(mod)
    except ImportError:
        missing.append(pip)
if missing:
    print("installing:", " ".join(missing))
    subprocess.run([sys.executable, "-m", "pip", "-q", "install", *missing], check=True)

import semgrit.build_deck as _bd
import mapbox_earcut, shapely, skimage, cv2, PIL          # noqa: F401
print("pipeline ready in", WORK)
print("versions   : Pillow %s, skimage %s, cv2 %s, shapely %s, earcut ok"
      % (PIL.__version__, skimage.__version__, cv2.__version__, shapely.__version__))
print("modules   :", len([f for f in os.listdir("semgrit") if f.endswith(".py")]),
      "in semgrit/, 4 verifiers, and vumat_grind.for + vumat_jh2.for")


def need(names, where):
    """Stop with the cell to run, instead of a NameError on an unfamiliar name."""
    missing = [n for n in names.split() if n not in globals()]
    if missing:
        raise SystemExit("run %s first - this cell needs %s"
                         % (where, ", ".join(missing)))
'''))

CELLS.append(md(r"""
---
## 1 · Load your SEM image

Zeiss SmartSEM `.tif` files carry the exact pixel size in TIFF tag 34118, and the
pipeline reads it. **That is the calibration that matters** — the burnt-in scale bar is
only used as a cross-check, and the run stops if the two disagree by more than 5 %.

For a non-Zeiss image with no usable metadata, set `PIXEL_SIZE_UM` in the next cell.
"""))

CELLS.append(code('''
#@title 📷 2 · Where are your SEM images? { display-mode: "form" }
SOURCE = "upload"  #@param ["upload", "google drive", "already on disk"]
#@markdown Used for **google drive** / **already on disk** — a folder or a glob:
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
#@title ⚡ 3 · SIMPLE — set it up and look at it { display-mode: "form" }
#@markdown Seven choices. This cell **writes nothing** — it measures your grains, works
#@markdown out the model and shows it to you. Change anything and re-run; when the model
#@markdown looks right, the next cell builds and downloads it.
#@markdown
#@markdown Everything not asked for here is the configuration the two Abaqus-validated
#@markdown decks were built with. Skip both cells if you want the Advanced path below.
RUN_SIMPLE = True                #@param {type:"boolean"}

#@markdown ### 1 · the wheel
S_DIAMETER_MM = 50.0             #@param {type:"number"}
#@markdown ### 2 · how much of it to model
S_SLICE_MM = 2.0                 #@param {type:"number"}
#@markdown &nbsp;&nbsp;Arc length of the slice. It must be longer than the workpiece.
#@markdown ### 3 · how many abrasives
S_GRITS = "concentration"        #@param ["concentration", "a fixed number", "grains per mm2", "single grain"]
S_GRIT_VALUE = 100.0             #@param {type:"number"}
#@markdown &nbsp;&nbsp;C-number for *concentration*, a count for *a fixed number*,
#@markdown grains/mm² for *grains per mm2*; ignored for *single grain*.
#@markdown ### 4 · the workpiece
S_WORKPIECE = "small  48 x 15 x 6 um"  #@param ["small  48 x 15 x 6 um", "medium  100 x 40 x 20 um", "large  200 x 200 x 200 um", "custom"]
S_CUSTOM_MM = "0.048 x 0.015 x 0.006"  #@param {type:"string"}
#@markdown &nbsp;&nbsp;`length x width x depth` in mm, used only when **custom**.
#@markdown ### 5 · where it sits on the wheel
S_POSITION = "centred"           #@param ["centred", "first grit at entry", "under the tallest grit", "custom angle"]
#@markdown ### 6 · the gap between wheel and work
S_STANDOFF_UM = 0.0              #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = the tallest grain under the block just touches it. The depth
#@markdown of cut is chosen automatically to close this gap and then cut 85% of the way
#@markdown through the grain protrusion.
#@markdown ### 7 · what you want out
S_OUTPUT = "run-ready .inp + CAE deck"  #@param ["run-ready .inp + CAE deck", "run-ready .inp only", "CAE deck only", "run-ready .inp + CAE deck + CAD"]
S_NAME = "wheel"                 #@param {type:"string"}
S_SHOW_CAD = True                #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;Show the 3-D viewer as well as the drawings. Untick it if you are
#@markdown iterating quickly and only want the numbers.

import os, time
import matplotlib.pyplot as plt
from IPython.display import HTML, display
from semgrit.quick import (SIMPLE_MEASURE, WORKPIECE_SIZES, library_summary,
                           measure_images, simple_params)
from semgrit.build_deck import plan_deck
from semgrit.preview import preview, summary_text
from semgrit.cadviewer import build as build_cad_view

if RUN_SIMPLE:
    WORK = globals().get("WORK", "/content/semgrit_work")
    OUT_MEAS = os.path.join(WORK, "1_measurements")
    OUT_DECK = os.path.join(WORK, "2_abaqus")
    _t0 = time.time()

    print("=" * 78)
    print("1/3  the grains")
    print("=" * 78)
    # Pixel size 0 = read it from the image metadata. Simple mode does not offer an
    # override, because an override is exactly the thing that silently rescales every
    # grain and therefore the whole wheel. Measuring is cached, so re-running after a
    # wheel change costs nothing.
    MEASURED = measure_images(IMAGES, OUT_MEAS, pixel_size_um=0.0, **SIMPLE_MEASURE)
    SOLIDS, ALL_GRAINS = MEASURED["solids"], MEASURED["grains"]
    if not MEASURED["cached"]:
        print()
        library_summary(SOLIDS)

    if S_WORKPIECE == "custom":
        _wp = tuple(float(x) for x in S_CUSTOM_MM.lower().replace(",", "x").split("x"))
        if len(_wp) != 3:
            raise SystemExit("S_CUSTOM_MM must be 'length x width x depth' in mm, "
                             "got %r" % S_CUSTOM_MM)
    else:
        _wp = WORKPIECE_SIZES[S_WORKPIECE]

    PARAMS = simple_params(
        diameter_mm=S_DIAMETER_MM, slice_mm=S_SLICE_MM, grit_kind=S_GRITS,
        grit_value=S_GRIT_VALUE, workpiece_mm=_wp, wp_position=S_POSITION,
        standoff_um=S_STANDOFF_UM,
        run_ready=S_OUTPUT != "CAE deck only",
        cae_deck="CAE deck" in S_OUTPUT,
        cad="CAD" in S_OUTPUT, name=S_NAME)

    print()
    print("=" * 78)
    print("2/3  what this will be  (nothing written yet)")
    print("=" * 78)
    PLAN = plan_deck(PARAMS, SOLIDS)
    print(summary_text(PLAN))
    print()
    _c = PLAN.get("cost") or {}
    print("COST       about %.0f MB of .inp, and roughly %.1f h to solve on 8 cores"
          % (PLAN["estimated_mb"], (_c.get("est_hours") or {}).get("8", 0.0)))
    if PLAN["estimated_mb"] > 400:
        print("           that is a very large deck - consider a shorter slice or "
              "fewer grits")

    print()
    print("=" * 78)
    print("3/3  look at it")
    print("=" * 78)
    fig = preview(PLAN)
    plt.show()
    if S_SHOW_CAD:
        _glb = os.path.join(WORK, S_NAME + "_cad.glb")
        try:
            _html, _meta, _ci = build_cad_view(PLAN, _glb, mode="whole wheel",
                                               max_grits=0, height=720)
            print("CAD viewer: %d of %d grains, %s triangles.  Press Contact to dive "
                  "to the grains." % (_meta["grits_drawn"], _meta["grits_total"],
                                      format(_ci["triangles"], ",")))
            display(HTML(_html))
        except ValueError as exc:
            print("viewer not shown: %s" % exc)

    print()
    print("-" * 78)
    print("Nothing has been written. If the wheel or the block is the wrong size,")
    print("change a value above and re-run this cell - the grains are cached, so it")
    print("comes back in a second. When it looks right, run the next cell to build")
    print("and download.   (%.0f s so far)" % (time.time() - _t0))
    print("-" * 78)
else:
    print("simple mode skipped - use the Advanced cells below")
'''))

CELLS.append(code('''
#@title ⚡ 4 · SIMPLE — build it, verify it, download it { display-mode: "form" }
#@markdown Only run this once the cell above shows the model you want. This is the step
#@markdown that writes the `.inp`, so it is also the slow one.
BUILD_AND_DOWNLOAD = True        #@param {type:"boolean"}
AUTO_DOWNLOAD = True             #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;Untick to leave the zip in the runtime instead of downloading it.

import os, time
from semgrit.quick import bundle, verify_decks
from semgrit.build_deck import build_deck

if BUILD_AND_DOWNLOAD and RUN_SIMPLE:
    if "PARAMS" not in globals() or "SOLIDS" not in globals():
        raise SystemExit("run the SIMPLE setup cell above first - it is what decides "
                         "what to build")
    _t0 = time.time()
    print("=" * 78)
    print("building %s ... (a big deck takes a minute or two)" % PARAMS.name)
    print("=" * 78)
    INFO = build_deck(PARAMS, SOLIDS, OUT_DECK)
    print("wrote %s  (%.1f MB, %.0f s)"
          % (os.path.basename(INFO["path"]), INFO["size_bytes"] / 1e6,
             time.time() - _t0))

    print()
    _decks = [INFO["path"]] + ([INFO["cae_deck"]] if INFO.get("cae_deck") else [])
    if not verify_decks(WORK, _decks):
        raise SystemExit("the deck did not verify - read the FAIL lines above")

    print()
    _zip = bundle(WORK, (OUT_DECK, OUT_MEAS), S_NAME)
    print()
    print("to run it:           abaqus job=%s input=%s user=vumat_jh2.for "
          "double=both cpus=8" % (S_NAME, os.path.basename(INFO["path"])))
    if INFO.get("postprocess_script"):
        print("to read the result:  abaqus python %s %s.odb"
              % (os.path.basename(INFO["postprocess_script"]), S_NAME))
    if AUTO_DOWNLOAD:
        try:
            from google.colab import files
            files.download(_zip)
        except Exception as exc:
            print("(not on Colab - copy the zip yourself)", exc)
elif not RUN_SIMPLE:
    print("simple mode is off - use the Advanced cells below")
else:
    print("not built - tick BUILD_AND_DOWNLOAD when you are happy with the preview")
'''))

CELLS.append(code('''
#@title 🔬 A1 · Calibration, segmentation and grain-solid settings { display-mode: "form" }

#@markdown ### Calibration
PIXEL_SIZE_UM = 0.0  #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = read it from the SEM metadata (**recommended**). Any other
#@markdown value overrides the metadata — only do this for non-Zeiss images.

#@markdown ### Segmentation
THRESHOLD = "multiotsu"  #@param ["multiotsu", "otsu"]
MIN_GRAIN_UM = 0.9        #@param {type:"number"}
H_MAXIMA_UM = 0.12        #@param {type:"number"}
#@markdown &nbsp;&nbsp;Low on purpose: it over-segments, then boundaries without real
#@markdown image evidence are merged back. Raise it if grains are being split.
GRADIENT_WEIGHT = 1.0     #@param {type:"number"}
MIN_EDGE_STRENGTH = 1.5   #@param {type:"number"}
MIN_AREA_UM2 = 0.7        #@param {type:"number"}
INCLUDE_BORDER_GRAINS = False  #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;Grains cut by the frame edge have truncated outlines; including
#@markdown them biases the size statistics low.

#@markdown ### Outline → solid
SIMPLIFY_UM = 0.10        #@param {type:"number"}
MAX_VERTICES = 64         #@param {type:"integer"}
THICKNESS_RATIO = 0.70    #@param {type:"number"}
#@markdown &nbsp;&nbsp;Grain height as a fraction of its minimum Feret width. An SEM
#@markdown gives no depth, so height is modelled, not measured.
THICKNESS_STD = 0.12      #@param {type:"number"}
BASE_SCALE = 0.70         #@param {type:"number"}
MID_HEIGHT = 0.42         #@param {type:"number"}
TOP_SCALE = 0.30          #@param {type:"number"}
#@markdown &nbsp;&nbsp;The lofted profile: the outline is scaled to these fractions at
#@markdown the base, the waist and the tip.

#@markdown ### Cutting edge
EDGE_RADIUS_UM = 0.0      #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` leaves knife edges, which are stress singularities in FEA.
#@markdown A good starting point is ~10 % of the measured d50 (printed below after the
#@markdown first run, so you can come back and set it).
ARC_SEGMENTS = 3          #@param {type:"integer"}

MEASURE_SEED = 20260728   #@param {type:"integer"}
print("settings captured - run the next cell to measure")
'''))

CELLS.append(code('''
#@title ▶ A2 · Measure the grains and build the 3-D grain library { display-mode: "form" }
# The body of this lives in semgrit.quick so that Simple mode runs the same code.
import os, numpy as np
from semgrit.quick import measure_images, library_summary
from semgrit.segment import SegmentationParams
from semgrit.grain3d import HeightModel, LoftProfile

OUT_MEAS = os.path.join(WORK, "1_measurements")

MEASURED = measure_images(
    IMAGES, OUT_MEAS, pixel_size_um=PIXEL_SIZE_UM,
    seg_params=SegmentationParams(
        min_grain_um=MIN_GRAIN_UM, h_maxima_um=H_MAXIMA_UM,
        gradient_weight=GRADIENT_WEIGHT, min_edge_strength=MIN_EDGE_STRENGTH,
        min_area_um2=MIN_AREA_UM2, threshold_method=THRESHOLD),
    height_model=HeightModel(mean_ratio=THICKNESS_RATIO, std_ratio=THICKNESS_STD,
                             seed=MEASURE_SEED),
    profile=LoftProfile(base_scale=BASE_SCALE, top_scale=TOP_SCALE,
                        mid_height_fraction=MID_HEIGHT,
                        edge_radius_um=EDGE_RADIUS_UM, arc_segments=ARC_SEGMENTS),
    simplify_um=SIMPLIFY_UM, max_vertices=MAX_VERTICES,
    interior_only=not INCLUDE_BORDER_GRAINS)
SOLIDS, ALL_GRAINS = MEASURED["solids"], MEASURED["grains"]

print()
LIB = library_summary(SOLIDS)
print()
print("  -> a sensible EDGE_RADIUS_UM is ~10%% of the d50 width = %.3f um"
      % (0.10 * LIB["width_um"][1]))
print("     (set it in A1 and re-run if you want blunted cutting edges)")
'''))

CELLS.append(code('''
#@title ✅ A3 · Verify the measurements (optional but recommended) { display-mode: "form" }
#@markdown Checks the half of the pipeline the deck verifiers cannot see: that your image
#@markdown was **calibrated** and **measured** correctly. The pixel size is re-read
#@markdown straight from the raw TIFF bytes, the scale bar is re-measured and multiplied
#@markdown out to confirm it equals its printed label, and every grain descriptor is
#@markdown recomputed from the label mask in plain numpy. Run it once for a new kind of
#@markdown image; skip it on repeat runs.
import subprocess, sys, os

r = subprocess.run([sys.executable, os.path.join(WORK, "verify_pipeline_A.py"),
                    "--quick", *IMAGES], capture_output=True, text=True, cwd=WORK)
print(r.stdout)
if r.stderr.strip():
    print(r.stderr[-2000:])
print("=" * 78)
print("MEASUREMENTS VERIFIED" if r.returncode == 0 else
      "MEASUREMENT CHECKS FAILED - the deck would be built from bad numbers")
print("=" * 78)
'''))

CELLS.append(md(r"""
---
## 3 · Design the wheel

**Wheel extent** — give it whichever way you think in:

| `SECTOR_MODE` | you set | typical use |
|---|---|---|
| `arc` | arc length in mm | you care about how much surface engages |
| `angle` | degrees (30, 90, 180 …) | you want a named sector |
| `full` | nothing — 360° | the complete wheel |

**Will the arc look curved?** The bow across a chord is `sagitta = L²/8R`. A 2 mm arc on
a Ø50 wheel bows 20 µm; against a 12 µm rim that reads clearly as an arc. Make the arc
short *and* the rim deep and it renders as a rectangle — the verifier warns you when
`sagitta < rim depth`.

**Grit population** — four ways:

| `GRIT_MODE` | you set | notes |
|---|---|---|
| `concentration` | C-number (C100 = 25 vol %) | the real abrasive spec |
| `areal_density` | grains / mm² | direct control |
| `count` | exactly N grains | easiest to reason about cost |
| `single` | one grain | single-grit scratch test |

At true C100 with fine grit the implied density is tens of thousands per mm², which no
mesh can hold over a large sector. Grains are rejected where they would overlap and the
achieved density is reported — read it, don't assume you got what you asked for.
"""))

CELLS.append(code('''
#@title ⚙️ A4 · Wheel and grit settings { display-mode: "form" }

#@markdown ### Wheel body
DIAMETER_MM = 50.0        #@param {type:"number"}
SECTOR_MODE = "arc"       #@param ["arc", "angle", "full"]
ARC_LENGTH_MM = 2.0       #@param {type:"number"}
SECTOR_DEG = 30.0         #@param {type:"number"}
RIM_DEPTH_MM = 0.012      #@param {type:"number"}
WHEEL_WIDTH_MM = 0.030    #@param {type:"number"}
BOND_DENSITY_KG_M3 = 2700.0  #@param {type:"number"}

#@markdown ### Rigid-shell mesh (appearance and contact only — never the time increment)
SHELL_CIRCUMFERENTIAL_DIVISIONS = 200  #@param {type:"integer"}
SHELL_AXIAL_DIVISIONS = 6              #@param {type:"integer"}
SHELL_RADIAL_DIVISIONS = 1             #@param {type:"integer"}

#@markdown ### Grits
GRIT_MODE = "concentration"  #@param ["concentration", "areal_density", "count", "single"]
CONCENTRATION = 100.0        #@param {type:"number"}
AREAL_DENSITY_PER_MM2 = 5000.0  #@param {type:"number"}
GRIT_COUNT = 500             #@param {type:"integer"}
#@markdown &nbsp;&nbsp;For **single**: `-1` picks the largest grain in the library.
SINGLE_GRAIN_INDEX = -1      #@param {type:"integer"}
SINGLE_GRIT_OFFSET_MM = 0.015  #@param {type:"number"}
#@markdown &nbsp;&nbsp;How far along the block the lone grit starts. Positive puts it at
#@markdown the trailing end, so a wheel turning toward **decreasing θ** (`VR3 < 0`)
#@markdown drags it across the whole workpiece.

#@markdown ### Seating
PROTRUSION_MEAN = 0.55   #@param {type:"number"}
PROTRUSION_STD = 0.12    #@param {type:"number"}
PROTRUSION_MIN = 0.25    #@param {type:"number"}
PROTRUSION_MAX = 0.85    #@param {type:"number"}
MAX_TILT_DEG = 35.0      #@param {type:"number"}
SPACING_FACTOR = 1.05    #@param {type:"number"}
GRIT_ARC_WINDOW_MM = 0.0   #@param {type:"number"}
#@markdown &nbsp;&nbsp;Dress only this much of the arc, centred. `0` = the whole sector.
#@markdown A 13 mm arc at 5000/mm² is 65,000 grains and hundreds of MB, and only the arc
#@markdown the block sweeps can ever touch it — so dress a window and leave the rest bare.
GRIT_FACE_WINDOW_MM = 0.0  #@param {type:"number"}
#@markdown &nbsp;&nbsp;Dress only this much of the wheel's face, centred. `0` = the full
#@markdown width. Lets the slice be thick enough to look like a real chunk of wheel while
#@markdown the grains stay in the band the workpiece actually runs in.
INSET_GRIT_BAND = True   #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;Keeps whole grains inside the bond. Turn it off only if you want
#@markdown grits sliced by the sector cut faces.
WHEEL_SEED = 20260731    #@param {type:"integer"}
print("wheel settings captured")
'''))

CELLS.append(code('''
#@title 🧱 A5 · Workpiece, kinematics and output { display-mode: "form" }

#@markdown ### Workpiece — the only deformable part
INCLUDE_WORKPIECE = True   #@param {type:"boolean"}
WP_LENGTH_MM = 0.048       #@param {type:"number"}
WP_WIDTH_MM = 0.015        #@param {type:"number"}
WP_DEPTH_MM = 0.006        #@param {type:"number"}

#@markdown #### Mesh size — element type is fixed at **C3D8R**, only the size is yours
WP_ELEMENT_SIZE_MM = 0.0003  #@param {type:"number"}
#@markdown &nbsp;&nbsp;The base size, used for any direction left at `0` below. Cost
#@markdown scales as **1/h⁴** if you change all three together — halving it multiplies
#@markdown the run by ~16. Aim for 5–10 elements through the deepest cut a grit takes.
WP_ELEM_CUTTING_MM = 0.0    #@param {type:"number"}
WP_ELEM_AXIAL_MM = 0.0      #@param {type:"number"}
WP_ELEM_DEPTH_MM = 0.0      #@param {type:"number"}
#@markdown &nbsp;&nbsp;**Graded depth mesh** — fine where the chip forms, coarse in the
#@markdown body. The chip is removed *into the depth*, so this direction is what resolves
#@markdown chip thickness; and it is free in time, because `dt` follows the smallest
#@markdown element and the surface layer only needs to match the cutting size.
WP_SURFACE_LAYER_MM = 0.0   #@param {type:"number"}
#@markdown &nbsp;&nbsp;Depth of the finely meshed zone at the ground face. `0` = uniform.
#@markdown Make it 2-3x your depth of cut; `WP_ELEM_DEPTH_MM` then sets its layer size.
WP_DEPTH_GROWTH = 1.3       #@param {type:"number"}
WP_MAX_DEPTH_ELEM_MM = 0.0  #@param {type:"number"}
#@markdown &nbsp;&nbsp;Cap on layer thickness so the deep elements do not become slivers.
#@markdown &nbsp;&nbsp;Per-direction overrides (`0` = use the base size). The three
#@markdown directions do **not** cost the same. The stable time increment follows the
#@markdown *smallest* element dimension, so coarsening **axial** alone drops the element
#@markdown count without lengthening the run — the cheapest saving available. Coarsening
#@markdown **cutting** or **depth** blurs the chip and the damage zone, so do that last.
#@markdown The block keeps the dimensions you asked for, so a size that does not divide
#@markdown them exactly is rounded; the achieved sizes are printed after the build.
WP_MATERIAL = "STONE"      #@param {type:"string"}
WP_DENSITY_KG_M3 = 2650.0  #@param {type:"number"}
WP_YOUNGS_MPA = 50000.0    #@param {type:"number"}
WP_POISSON = 0.25          #@param {type:"number"}

#@markdown #### Where the block sits on the wheel
#@markdown The wheel turns so its surface travels toward **decreasing theta**, so
#@markdown grains arrive from the high-theta end. That end is the *entry*.
WP_POSITION = "centred"    #@param ["centred", "first grit at entry", "under the tallest grit", "custom angle"]
#@markdown &nbsp;&nbsp;**centred** — mid-arc, grain either side however the wheel turns.
#@markdown **first grit at entry** — the block's entry edge sits at the leading grain,
#@markdown so the pass starts with the first abrasive right at the edge and every grain
#@markdown downstream then sweeps across it. **under the tallest grit** — centred on the
#@markdown most protruding grain the block can reach, the one that takes the deepest
#@markdown cut. **custom angle** — you name it, below.
WP_POSITION_DEG = 0.0      #@param {type:"number"}
#@markdown &nbsp;&nbsp;Only used by **custom angle**. Measured from the global +X axis;
#@markdown the preview prints the angular span the grits occupy so you can aim at them.

#@markdown #### Standoff — the gap between wheel and workpiece
CLEARANCE_UM = 0.0         #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = the tallest grit that can reach the block is exactly
#@markdown tangent to it: contact at one point, zero initial overclosure. A positive
#@markdown value parks the block that many microns clear, so the infeed has to close
#@markdown the gap before anything cuts.
#@markdown
#@markdown **Cell **A8** reports how tall the abrasive actually stands** — minimum, maximum
#@markdown and mean protrusion above the bond — and the depth-of-cut window each
#@markdown standoff gives you. Run A8, read the numbers, then come back and set this.
#@markdown A standoff wider than the depth of cut means the wheel turns for the whole
#@markdown step and never touches the work; the build refuses rather than let that
#@markdown happen.

#@markdown ### Kinematics (used for the run-time estimate, not written into the deck)
SURFACE_SPEED_M_S = 30.0   #@param {type:"number"}
TRAVEL_MM = 0.0            #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = the block length plus the run-in below.
TRAVEL_MARGIN_MM = 0.006   #@param {type:"number"}
CORES = 8                  #@param {type:"integer"}

#@markdown ### Which files do you want?
MODEL_NAME = "wheel"       #@param {type:"string"}
#@markdown &nbsp;&nbsp;**Abaqus decks** — you can have both from one run. They are written
#@markdown from the same placed grits, so they are the same wheel.
WRITE_RUN_READY_INP = True   #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;`<name>.inp` — submit from the terminal, no CAE. Configured in the
#@markdown next cell.
WRITE_CAE_INP = True         #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;`<name>_cae.inp` + `<name>_import_into_cae.py` — geometry only, to
#@markdown assemble and set up by hand in CAE.

#@markdown &nbsp;&nbsp;**CAD of the assembled wheel** (SOLIDWORKS). STEP is a faceted B-rep
#@markdown and is far heavier per body than the FE mesh — cap it on a wheel with
#@markdown thousands of grits.
WRITE_WHEEL_STEP = False   #@param {type:"boolean"}
WRITE_WHEEL_STL = False    #@param {type:"boolean"}
STEP_MAX_GRAINS = 0        #@param {type:"integer"}
STL_MAX_GRAINS = 0         #@param {type:"integer"}

#@markdown &nbsp;&nbsp;**CAD of the individual grits**, laid out on a grid rather than at
#@markdown their wheel positions — this is what you open to inspect or measure one grain.
WRITE_GRAINS_STEP = False  #@param {type:"boolean"}
WRITE_GRAIN_STLS = False   #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;`WRITE_GRAIN_STLS` writes one `.stl` per measured grain into
#@markdown `grits_stl/` — handy, but it is one file per grain.
GRAINS_STEP_MAX = 200      #@param {type:"integer"}
print("workpiece and output settings captured")
'''))

CELLS.append(code('''
#@title 🚀 A6 · Run-ready analysis — submit from the terminal, no CAE { display-mode: "form" }
#@markdown These apply when **`WRITE_RUN_READY_INP`** is ticked in the previous cell.
RUN_READY = True  #@param {type:"boolean"}
#@markdown Leave on. Turning it off here also disables the run-ready deck:
#@markdown ```
#@markdown abaqus job=grind input=<name>.inp user=vumat_jh2.for double=both cpus=8 interactive
#@markdown ```

#@markdown ### Cutting
DEPTH_OF_CUT_UM = 0.0     #@param {type:"number"}
#@markdown &nbsp;&nbsp;**The one number that decides whether anything is ground at all.**
#@markdown Leave it at **`0` for automatic** — 85 % of whatever bond clearance this wheel
#@markdown turns out to have, which is always valid. Set a number to override.
#@markdown A wheel given only a rotation spins on its own axis: one grit grazes at t=0 and
#@markdown every grit behind it stays a micron below the surface for ever. The build refuses
#@markdown a depth greater than the bond-rim clearance, so the rim cannot hit the work.
STEP_TIME_S = 0.0         #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = travel / surface speed.
MASS_SCALING = 10.0       #@param {type:"number"}
BULK_VISCOSITY_LINEAR = 0.06     #@param {type:"number"}
BULK_VISCOSITY_QUADRATIC = 1.2   #@param {type:"number"}
#@markdown &nbsp;&nbsp;Multiplies density, so it lengthens `dt` by its square root and speeds
#@markdown the run up by the same factor — at the cost of distorting inertia. 1 disables it.
NLGEOM = True             #@param {type:"boolean"}

#@markdown ### Workpiece material
MATERIAL_MODEL = "jh2"    #@param ["jh2", "elastic"]
JH2_DENSITY_KG_M3 = 2350.0  #@param {type:"number"}
JH2_CONSTANTS = "3735.6, 2686, 1982, 1374, 8, 0.71, 0.30, 0.022, 0.55, 0.40, 1.0, 0.002, 1.20, 9000, 22000, 0.25, 912"  #@param {type:"string"}
#@markdown &nbsp;&nbsp;17 values in the order the VUMAT reads them:
#@markdown `K1 G HEL PHEL T A B C N M beta D1 D2 K2 K3 SFMAX SIGHEL`
N_DEPVAR = 12             #@param {type:"integer"}
ELEMENT_DELETION = True   #@param {type:"boolean"}
HOURGLASS = "ENHANCED"    #@param ["ENHANCED", "RELAX STIFFNESS", "STIFFNESS", "VISCOUS"]

#@markdown ### Contact and how the block is held
CONTACT_SCOPE = "engaging"  #@param ["engaging", "all exterior", "none"]
#@markdown &nbsp;&nbsp;`engaging` pairs only the grits that can reach the block — far cheaper
#@markdown than tracking half a million facets.
FRICTION = 0.2            #@param {type:"number"}
FIX_BACK_FACE = True      #@param {type:"boolean"}
FIX_ENDS = False          #@param {type:"boolean"}
FIX_SIDES = False         #@param {type:"boolean"}

#@markdown ### Output
FIELD_FRAMES = 60         #@param {type:"integer"}
RESTART_INTERVALS = 10    #@param {type:"integer"}
#@markdown &nbsp;&nbsp;Must be > 1 to be recoverable: with 1 the only restart state is written
#@markdown at the *end* of the step, so an interrupted run cannot be resumed at all.
ELEMENT_OUTPUT = "S, PEEQ, SDV, STATUS"  #@param {type:"string"}
NODE_OUTPUT = "U, V"      #@param {type:"string"}
HISTORY_PRESELECT = True  #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;Whole-model energies — `ALLKE` is how you confirm the wheel is
#@markdown actually turning, so leave this on.
ROTATION_REVERSED = False      #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;Turn the wheel the other way. The surface then travels toward
#@markdown *increasing* theta, so grains arrive from the block's **low**-theta end — and
#@markdown the `first grit at entry` placement follows, because the entry edge is
#@markdown whichever end the grains reach first. The deck header states the sense it
#@markdown actually applies, and the verifier checks the sentence against the sign.
HISTORY_REFERENCE_NODE = True  #@param {type:"boolean"}
HISTORY_INTERVALS = 200        #@param {type:"integer"}
#@markdown &nbsp;&nbsp;Reaction force and moment at the wheel's reference node, sampled
#@markdown this many times. **This is the grinding force** — `PRESELECT` does not include
#@markdown it, so with this off the job finishes and the `.odb` holds no force to plot.
#@markdown The post-processing script written with the deck reads exactly these.
print("analysis settings captured")
'''))

CELLS.append(code('''
#@title 👁 A7 · PREVIEW — see it before you build it { display-mode: "form" }
#@markdown Draws the whole assembly from the **same placement code the writer uses**, so
#@markdown what you see is what the deck will contain — without writing a single file.
#@markdown Change anything in the cells above and re-run this until it looks right.
SHOW_PREVIEW = True  #@param {type:"boolean"}

import math, os
import matplotlib.pyplot as plt
from semgrit.analysis import AnalysisParams
from semgrit.build_deck import DeckParams, plan_deck
from semgrit.preview import preview, summary_text


def make_params():
    """One place both the preview and the build read their settings from."""
    an = AnalysisParams(
        enabled=bool(RUN_READY and WRITE_RUN_READY_INP),
        step_time_s=STEP_TIME_S, nlgeom=NLGEOM,
        mass_scaling_factor=MASS_SCALING,
        bulk_viscosity=(BULK_VISCOSITY_LINEAR, BULK_VISCOSITY_QUADRATIC),
        depth_of_cut_um=DEPTH_OF_CUT_UM,
        material_model=MATERIAL_MODEL,
        jh2_constants=[float(x) for x in JH2_CONSTANTS.split(",")],
        jh2_density_kg_m3=JH2_DENSITY_KG_M3, n_depvar=N_DEPVAR,
        element_deletion=ELEMENT_DELETION, hourglass=HOURGLASS,
        contact_scope=CONTACT_SCOPE, friction=FRICTION,
        fix_back_face=FIX_BACK_FACE, fix_ends=FIX_ENDS, fix_sides=FIX_SIDES,
        field_frames=FIELD_FRAMES, restart_intervals=RESTART_INTERVALS,
        element_output=ELEMENT_OUTPUT, node_output=NODE_OUTPUT,
        history_preselect=HISTORY_PRESELECT,
        rotation_reversed=ROTATION_REVERSED,
        history_reference_node=HISTORY_REFERENCE_NODE,
        history_intervals=HISTORY_INTERVALS)
    return DeckParams(
        diameter_mm=DIAMETER_MM, sector_mode=SECTOR_MODE, sector_deg=SECTOR_DEG,
        arc_length_mm=ARC_LENGTH_MM, rim_depth_mm=RIM_DEPTH_MM, width_mm=WHEEL_WIDTH_MM,
        shell_circumferential_divisions=SHELL_CIRCUMFERENTIAL_DIVISIONS,
        shell_axial_divisions=SHELL_AXIAL_DIVISIONS,
        shell_radial_divisions=SHELL_RADIAL_DIVISIONS,
        bond_density_kg_m3=BOND_DENSITY_KG_M3,
        grit_mode=GRIT_MODE, concentration=CONCENTRATION,
        areal_density_per_mm2=AREAL_DENSITY_PER_MM2, grit_count=GRIT_COUNT,
        single_grain_index=SINGLE_GRAIN_INDEX,
        single_grit_offset_mm=SINGLE_GRIT_OFFSET_MM,
        grit_arc_window_mm=GRIT_ARC_WINDOW_MM, grit_width_window_mm=GRIT_FACE_WINDOW_MM,
        inset_grit_band=INSET_GRIT_BAND,
        protrusion_mean=PROTRUSION_MEAN, protrusion_std=PROTRUSION_STD,
        protrusion_min=PROTRUSION_MIN, protrusion_max=PROTRUSION_MAX,
        max_tilt_deg=MAX_TILT_DEG, spacing_factor=SPACING_FACTOR, seed=WHEEL_SEED,
        include_workpiece=INCLUDE_WORKPIECE, wp_length_mm=WP_LENGTH_MM,
        wp_width_mm=WP_WIDTH_MM, wp_depth_mm=WP_DEPTH_MM,
        wp_element_size_mm=WP_ELEMENT_SIZE_MM,
        wp_element_size_length_mm=WP_ELEM_CUTTING_MM,
        wp_element_size_width_mm=WP_ELEM_AXIAL_MM,
        wp_element_size_depth_mm=WP_ELEM_DEPTH_MM,
        wp_surface_layer_mm=WP_SURFACE_LAYER_MM, wp_depth_growth=WP_DEPTH_GROWTH,
        wp_max_depth_element_mm=WP_MAX_DEPTH_ELEM_MM,
        wp_material=WP_MATERIAL, wp_density_kg_m3=WP_DENSITY_KG_M3,
        wp_youngs_modulus_mpa=WP_YOUNGS_MPA, wp_poisson_ratio=WP_POISSON,
        clearance_um=CLEARANCE_UM, wp_position=WP_POSITION,
        wp_position_deg=WP_POSITION_DEG,
        surface_speed_mm_s=SURFACE_SPEED_M_S * 1000.0, travel_mm=TRAVEL_MM,
        travel_margin_mm=TRAVEL_MARGIN_MM, cores=CORES,
        analysis=an, also_write_cae_deck=WRITE_CAE_INP,
        name=MODEL_NAME, write_step=WRITE_WHEEL_STEP, write_stl=WRITE_WHEEL_STL,
        step_max_grains=STEP_MAX_GRAINS, stl_max_grains=STL_MAX_GRAINS,
        write_grain_stls=WRITE_GRAIN_STLS, write_grains_step=WRITE_GRAINS_STEP,
        grains_step_max=GRAINS_STEP_MAX)


need("SOLIDS", "A2 (measure the grains), or the SIMPLE cells")
PARAMS = make_params()

# This cell is the reset point for anything edited in the CAD viewer. Re-running the
# preview means you are driving from the widgets again, so viewer edits are dropped here
# rather than surviving invisibly into the build -- and it says when it drops some.
if globals().get("EDITED_PARAMS") is not None:
    print("note: the CAD viewer's edits (%s) are dropped -- this preview and the build"
          % ", ".join(globals().get("EDITED_CHANGED") or ["none"]))
    print("      now follow the widgets above. Re-run A12b to apply them again.")
    print()
EDITED_PARAMS = None
EDITED_BASE = None
EDITED_CHANGED = []
EDITED_SETTINGS = {}

if SHOW_PREVIEW:
    PLAN = plan_deck(PARAMS, SOLIDS)
    print(summary_text(PLAN))
    print()
    fig = preview(PLAN)
    plt.show()
    print("Happy with it? Run the next cell to build. Otherwise change a setting above")
    print("and re-run this cell - nothing has been written yet.")
else:
    print("preview skipped")
'''))

CELLS.append(code('''
#@title 📏 A8 · Abrasive heights, and what standoff to use { display-mode: "form" }
#@markdown How tall the grains actually stand, and the depth-of-cut window that
#@markdown follows. Read this, then set `WP_POSITION` and `CLEARANCE_UM` in A5.
#@markdown
#@markdown The standoff is measured from the **tallest grain under the block**, so a
#@markdown standoff of 0 means that grain touches the workpiece with zero overclosure.
#@markdown Every micron of standoff you add is a micron the infeed has to give back
#@markdown before anything cuts — the table below does that arithmetic for you.
STANDOFF_TABLE = True      #@param {type:"boolean"}

import numpy as _np

if "PLAN" not in globals():
    PLAN = plan_deck(PARAMS, SOLIDS)

_pa = PLAN["protrusion_um"]
_pu = PLAN["protrusion_under_block_um"]
_gh = PLAN["grain_height_um"]
print("ABRASIVE HEIGHT  (protrusion above the bond, microns)")
print("  %-26s %8s %8s %8s %8s %6s" % ("", "min", "median", "mean", "max", "n"))
for _lab, _d in (("every grain on the wheel", _pa),
                 ("grains under the block", _pu)):
    if _d["n"]:
        print("  %-26s %8.3f %8.3f %8.3f %8.3f %6d"
              % (_lab, _d["min"], _d["median"], _d["mean"], _d["max"], _d["n"]))
if _gh["n"]:
    print("  %-26s %8.3f %8.3f %8.3f %8.3f %6d"
          % ("grain height, as measured", _gh["min"], _gh["median"], _gh["mean"],
             _gh["max"], _gh["n"]))
_p = _np.asarray(PLAN["_place"]["protrusion_um"], dtype=float)
if _p.size:
    print("  percentiles  " + "  ".join(
        "%d%%=%.2f" % (q, _np.percentile(_p, q)) for q in (10, 25, 50, 75, 90)))

print()
print("WHERE THE BLOCK SITS")
print("  position          : %s" % PLAN["wp_position"])
print("  block spans theta : %.4f deg (entry) to %.4f deg, over %d grain(s)"
      % (PLAN["wp_entry_theta_deg"], PLAN["wp_exit_theta_deg"],
         PLAN["n_grits_under_block"]))
print("  grit spans theta  : %.4f to %.4f deg  (%.4f to %.4f within the block's "
      "width)" % (PLAN["grit_theta_range_deg"] + PLAN["grit_theta_reachable_deg"]))
print("  the surface travels toward DECREASING theta, so grains arrive from the")
print("  high-theta end - that end is the entry.")
if PLAN["wp_relocated"]:
    print("  NOTE the footprint you asked for held no grit, so the block was moved to")
    print("       the tallest grain it can reach.")

_s0 = PLAN["standoff_um"]
_f0 = PLAN["first_contact_um"]
_c0 = PLAN["depth_ceiling_um"]
if STANDOFF_TABLE and _f0 is not None:
    # A standoff only lifts the ground face; it shifts both ends of the window by
    # exactly the same amount, so the table is exact without rebuilding anything.
    print()
    print("DEPTH-OF-CUT WINDOW vs STANDOFF   (microns)")
    print("  %10s  %14s  %14s  %s" % ("standoff", "first contact", "bond hits",
                                      "auto ae"))
    _cand = sorted({0.0, round(_s0, 3), round(0.25 * _pa["max"], 3),
                    round(0.50 * _pa["max"], 3), round(_pa["max"], 3)})
    for _s in _cand:
        _lo, _hi = _f0 - _s0 + _s, _c0 - _s0 + _s
        print("  %10.3f  %14.3f  %14.3f  %.3f%s"
              % (_s, _lo, _hi, _s + 0.85 * (_c0 - _s0),
                 "   <- current" if abs(_s - _s0) < 1e-9 else ""))
    print("  Pick DEPTH_OF_CUT_UM strictly between the two middle columns.")
    print("  DEPTH_OF_CUT_UM = 0 asks for the automatic value in the last column.")
'''))

CELLS.append(code('''
#@title 📐 A9 · Grinding theory — is this a real grinding regime? { display-mode: "form" }
#@markdown Two columns. **Measured** is counted off the geometry the deck contains:
#@markdown grain density from the grains that were placed, active grains from the ones
#@markdown that reach the work at this infeed, mesh resolution from the elements that
#@markdown were written. **Classical** is the textbook expressions.
#@markdown
#@markdown Those formulas assume a *traverse* grind at a work speed. This deck is a
#@markdown plunge — fixed block, radial infeed — so give the traverse case you want to
#@markdown compare against. Leave it at `0` and the classical column reports only what
#@markdown needs no work speed, rather than quietly using zero.
SHOW_THEORY = True         #@param {type:"boolean"}
WORK_SPEED_MM_S = 0.0      #@param {type:"number"}
#@markdown &nbsp;&nbsp;Table speed of the equivalent traverse grind, mm/s. 0 = skip those rows.
CHIP_SHAPE_FACTOR = 10.0   #@param {type:"number"}
#@markdown &nbsp;&nbsp;Chip width-to-thickness ratio `r` in Malkin's `h_max`. Not
#@markdown measurable from this model and the literature spans about 5 to 20, so it is
#@markdown yours to state.

from semgrit.grinding_theory import format_report, report as theory_report

if SHOW_THEORY:
    if "PLAN" not in globals():
        PLAN = plan_deck(PARAMS, SOLIDS)
    THEORY = theory_report(PLAN, work_speed_mm_s=WORK_SPEED_MM_S,
                           shape_factor=CHIP_SHAPE_FACTOR)
    print(format_report(THEORY))
else:
    print("theory report skipped")
'''))

CELLS.append(code('''
#@title 🧊 A10 · Quick 3-D scatter view (Plotly) { display-mode: "form" }
#@markdown Drag to rotate, scroll to zoom. This draws the **same triangles the deck
#@markdown contains** — the rim shell, every measured grain, and the workpiece block —
#@markdown so what you orbit here is literally what Abaqus will read.
SHOW_3D = True             #@param {type:"boolean"}
VIEW_MODE = "contact"      #@param ["contact", "wheel"]
#@markdown &nbsp;&nbsp;`contact` clips to a window around the workpiece — the only zoom
#@markdown at which 3 µm grains are visible on a 50 mm wheel. `wheel` shows the whole
#@markdown sector for proportion, with the grits necessarily sub-pixel.
MAX_GRITS_DRAWN = 400      #@param {type:"integer"}
#@markdown &nbsp;&nbsp;A browser starts to struggle past ~100k triangles and a grain is
#@markdown ~116 of them. Grains nearest the block are drawn first, and the number
#@markdown actually drawn is reported.
VIEW_WINDOW_UM = 0         #@param {type:"number"}
#@markdown &nbsp;&nbsp;Size of the `contact` window. `0` = 1.8x the workpiece.
SHOW_BOND_IN_3D = True     #@param {type:"boolean"}
SHOW_WORKPIECE_IN_3D = True  #@param {type:"boolean"}

from semgrit.viewer import view3d

if SHOW_3D:
    if "PLAN" not in globals():
        PLAN = plan_deck(PARAMS, SOLIDS)
    FIG3D, _drew = view3d(PLAN, mode=VIEW_MODE, max_grits=MAX_GRITS_DRAWN,
                          window_um=VIEW_WINDOW_UM, show_bond=SHOW_BOND_IN_3D,
                          show_workpiece=SHOW_WORKPIECE_IN_3D)
    _tri = _drew.get("grit_triangles", 0) + _drew.get("bond_triangles", 0)
    print("drawing %s of %s grits (%s in this view), %s triangles"
          % (format(_drew.get("grits_drawn", 0), ","),
             format(_drew.get("grits_total", 0), ","),
             format(_drew.get("grits_in_view", 0), ","), format(_tri, ",")))
    if _drew.get("grits_in_view", 0) > _drew.get("grits_drawn", 0):
        print("  capped by MAX_GRITS_DRAWN - raise it to see the rest")
    if _tri > 150000:
        print("  that is a lot of triangles; if it is sluggish, lower MAX_GRITS_DRAWN")
    FIG3D.show()
else:
    print("3D view skipped")
'''))

CELLS.append(code('''
#@title ✨ A11 · glTF view (also opens in Blender and PowerPoint) { display-mode: "form" }
#@markdown Google's `<model-viewer>` renders a real **glTF** file with physically-based
#@markdown lighting, soft shadows and orbit controls. **No API key, no account, nothing
#@markdown uploaded** — the model is written here and rendered in your browser.
#@markdown
#@markdown The `.glb` it writes is a genuine CAD interchange file: it also opens in
#@markdown **Blender**, **Windows 3D Viewer** and **PowerPoint** — useful for a slide.
SHOW_GLTF = True          #@param {type:"boolean"}
GLTF_MODE = "contact"     #@param ["contact", "wheel"]
GLTF_MAX_GRITS = 400      #@param {type:"integer"}
GLTF_MAX_INLINE_MB = 12.0 #@param {type:"number"}
#@markdown &nbsp;&nbsp;The file is embedded in the output as a data URI, which inflates
#@markdown it by a third. Past this cap it refuses rather than bloating the notebook —
#@markdown lower `GLTF_MAX_GRITS`, or just download the `.glb` and open it in Blender.

import os
from IPython.display import HTML, display
from semgrit.glb import model_viewer_html, parts_from_plan, write_glb

if SHOW_GLTF:
    if "PLAN" not in globals():
        PLAN = plan_deck(PARAMS, SOLIDS)
    GLB_PATH = os.path.join(WORK, MODEL_NAME + "_view.glb")
    _i = write_glb(GLB_PATH, parts_from_plan(PLAN, mode=GLTF_MODE,
                                             max_grits=GLTF_MAX_GRITS))
    print("wrote %s  (%.2f MB, %d parts)" % (os.path.basename(GLB_PATH),
                                             _i["bytes"] / 1e6, _i["parts"]))
    try:
        display(HTML(model_viewer_html(GLB_PATH, max_inline_mb=GLTF_MAX_INLINE_MB)))
    except ValueError as exc:
        print("not embedded: %s" % exc)
else:
    print("glTF view skipped")
'''))

CELLS.append(code('''
#@title 🛠 A12 · CAD viewer — shaded with edges, section planes, click-to-inspect { display-mode: "form" }
#@markdown A full CAD viewer, built on **three.js**, running in this output cell. It is
#@markdown the same geometry the `.inp` contains — not a re-mesh, not an approximation —
#@markdown so what you inspect here is what Abaqus will solve.
#@markdown
#@markdown | | |
#@markdown |---|---|
#@markdown | **Shaded with edges** | the SolidWorks look: feature edges over a lit surface |
#@markdown | **Wheel / Contact** | jump between the whole 50 mm wheel and the grains on the work |
#@markdown | **Face / Axial** | look straight at the dressed surface, or down the wheel axis |
#@markdown | **Section plane** | cut on any axis and drag the slider through the model |
#@markdown | **Click a grain** | its id, protrusion, height, width, volume and position |
#@markdown | **Shift-click twice** | distance **and** ΔX ΔY ΔZ, plus radial / along-arc / across-face |
#@markdown | **Parts tree** | show or hide the bond, the grits, the workpiece |
#@markdown | **Boundary conditions** | ENCASTRE pins on the held faces, the infeed arrow, the rotation arc, the reference node, the contact surfaces — every symbol standing for a keyword the deck really writes |
#@markdown | **Drag block** (`G`) | drag the workpiece along the arc, shift-drag for standoff, arrow keys nudge by 0.1 µm, `Esc` cancels |
#@markdown | **Depth-of-cut band** | the valid window shaded green between *nothing touches* and *bond hits the work* — the two ways a run has already been wasted |
#@markdown | **Save PNG** | a figure for the report |
#@markdown
#@markdown The boundary conditions are read out of the deck, not decorated on: a
#@markdown geometry-only deck shows none, the rotation arc carries the sign of `VR3`, and
#@markdown the red grains are the deck's own `ES_GRITS_ENGAGE` set. Held-face symbols are
#@markdown sampled for legibility but the panel always states the true node count.
#@markdown
#@markdown No account, no API key, nothing uploaded. three.js loads from a CDN and the
#@markdown model is embedded in the page.
SHOW_CAD_VIEWER = True      #@param {type:"boolean"}
CAD_MODE = "whole wheel"    #@param ["whole wheel", "wheel", "contact"]
CAD_MAX_GRITS = 0           #@param {type:"integer"}
CAD_HEIGHT = 720            #@param {type:"integer"}
CAD_MAX_INLINE_MB = 24.0    #@param {type:"number"}
#@markdown &nbsp;&nbsp;**whole wheel** opens on the complete wheel — the only view where the
#@markdown curvature of a 2 mm slice is visible at all — with your slice on it and an
#@markdown orange marker at the contact; press **Contact** to dive to the grains. Both the
#@markdown ghost wheel and the marker are pointers, labelled as such in the parts tree, and
#@markdown they fade out as you zoom in. **wheel** is the modelled sector alone;
#@markdown **contact** is just the patch under the block, and is the fastest.
#@markdown
#@markdown &nbsp;&nbsp;`CAD_MAX_GRITS = 0` draws **every** grain. If that would exceed
#@markdown `CAD_MAX_INLINE_MB`, grains far from the contact are drawn as boxes rather than
#@markdown dropped, and the cell says how many — you always see the whole wheel.

import os
from IPython.display import HTML, display
from semgrit.cadviewer import build as build_cad_view

if SHOW_CAD_VIEWER:
    if "PLAN" not in globals():
        PLAN = plan_deck(PARAMS, SOLIDS)
    # If this is Colab, let the viewer's Apply button commit straight into the kernel.
    # Everywhere else the viewer falls back to exporting the settings, which is why the
    # feature is not built on this channel existing.
    try:
        from google.colab import output as _colab_out
        from semgrit.build_deck import DeckError
        from semgrit.editable import apply as _edit_apply
        from semgrit.editable import commit_reply as _cad_reply

        def _cad_commit(settings):
            # The return value must be commit_reply(...), not a plain dict: Colab runs it
            # through IPython's display formatter, and a dict formats to text/plain only,
            # which the viewer cannot read. That is what once made every Apply -- the
            # successful ones too -- come back as "Python refused it, no reason given".
            global PARAMS, PLAN
            try:
                got = _edit_apply(settings, PARAMS, SOLIDS)
            except DeckError as exc:
                return _cad_reply(False, error=str(exc))   # already user-facing prose
            except Exception as exc:
                # Anything else is a bug, not a rejected setting. Name it as one.
                return _cad_reply(False, error="%s: %s" % (type(exc).__name__, exc))
            try:
                PARAMS, PLAN = got["params"], got["plan"]
                with open(os.path.join(WORK, "viewer_settings.json"), "w") as fh:
                    import json as _json
                    _json.dump({"settings": got["settings"]}, fh, indent=1)
            except Exception as exc:
                # The edit is already live in PARAMS; only the record of it failed.
                return _cad_reply(True, message="applied, but writing "
                                  "viewer_settings.json failed: %s" % exc)
            return _cad_reply(
                True, message="%s changed (%s). Re-run this cell to redraw, then the "
                              "build cell to write the deck."
                              % (", ".join(got["changed"]) or "nothing", got["tier"]))

        _colab_out.register_callback("cad.commit", _cad_commit)
        print("Apply is live: edits commit straight into this kernel.")
    except Exception:
        print("Apply will export settings (no Colab kernel channel here).")

    CAD_GLB = os.path.join(WORK, MODEL_NAME + "_cad.glb")
    try:
        CAD_HTML, CAD_META, _ci = build_cad_view(
            PLAN, CAD_GLB, mode=CAD_MODE, max_grits=CAD_MAX_GRITS,
            height=CAD_HEIGHT, max_inline_mb=CAD_MAX_INLINE_MB)
        _nf = len(CAD_META["grains_far"])
        print("%s  |  %d parts, %s triangles, %d of %d grains drawn%s"
              % (os.path.basename(CAD_GLB), _ci["parts"],
                 format(_ci["triangles"], ","), CAD_META["grits_drawn"],
                 CAD_META["grits_total"],
                 " (%d of them as boxes)" % _nf if _nf else ", all in full detail"))
        for _n in CAD_META["notes"]:
            print("   note: %s" % _n)
        display(HTML(CAD_HTML))
    except ValueError as exc:
        print("viewer not shown: %s" % exc)
else:
    print("CAD viewer skipped")
'''))

CELLS.append(code('''
#@title ✏️ A12b · Rebuild from the viewer's edits { display-mode: "form" }
#@markdown Paste what the CAD viewer's **Copy JSON** gave you, or leave this blank and it
#@markdown reads `viewer_settings.json` from the working folder (what **Download** saves,
#@markdown and what a live **Apply** writes).
#@markdown
#@markdown Every edit goes through one Python function — `semgrit.editable.apply` — so a
#@markdown number typed in the browser reaches the deck by exactly the same path whether
#@markdown it arrived through the kernel, a file or your clipboard.
APPLY_VIEWER_EDITS = True   #@param {type:"boolean"}
PASTED_SETTINGS = ""        #@param {type:"string"}

import os
from semgrit.editable import apply as edit_apply
from semgrit.editable import load as edit_load
from semgrit.editable import param_block
from semgrit.preview import summary_text

_edits_src = ""
if APPLY_VIEWER_EDITS:
    need("PARAMS SOLIDS", "the SIMPLE cells or A7 (preview)")
    _src = PASTED_SETTINGS.strip() or os.path.join(WORK, "viewer_settings.json")
    if PASTED_SETTINGS.strip() or os.path.exists(_src):
        _edits_src = _src

if _edits_src:
    _base_for_block = PARAMS
    _got = edit_apply(edit_load(_edits_src), PARAMS, SOLIDS)
    PARAMS, PLAN = _got["params"], _got["plan"]
    # A13 rebuilds PARAMS from the widgets, which would throw this away. These four names
    # are how the build cell knows an edit is in force, and what it was based on.
    # EDITED_BASE stays the *widget* baseline across repeated applies, so re-running this
    # cell does not make the build cell think the widgets have drifted.
    EDITED_PARAMS = PARAMS
    if globals().get("EDITED_BASE") is None:
        EDITED_BASE = _base_for_block
    EDITED_CHANGED = list(_got["changed"])
    EDITED_SETTINGS = {_k: _got["settings"][_k] for _k in _got["changed"]}
    print("applied: %s" % (", ".join(_got["changed"]) or "nothing changed"))
    print("tier   : %s" % _got["tier"])
    print()
    print(summary_text(PLAN))
    print()
    print("These are now the settings the build cell will use. The widgets above still")
    print("show their old values -- they cannot be written back to -- so the numbers")
    print("printed here are the authoritative ones.")
    _blk = param_block(_got["settings"], _base_for_block)
    if _blk:
        print()
        print("To bring the form widgets back in step, paste these into the cells above")
        print("(names and units are the widgets' own -- note WHEEL_WIDTH_MM, and")
        print("SURFACE_SPEED_M_S which is in m/s):")
        print()
        for _l in _blk.split(chr(10)):
            print("    " + _l)
elif APPLY_VIEWER_EDITS:
    # Nothing to apply is the ordinary case -- most runs never touch the viewer. It is
    # not an error, and it must not stop the notebook: A13 below still has to build.
    print("no edits found, so the build below uses the widget values as they stand.")
    print("To edit from the viewer: Copy JSON there and paste it above, or Download and")
    print("upload viewer_settings.json to " + WORK)
else:
    print("viewer edits not applied")
'''))

CELLS.append(code('''
#@title ▶ A13 · Build the Abaqus deck { display-mode: "form" }
# Uses exactly the parameters the preview just drew.
import os, math
from semgrit.build_deck import build_deck

need("SOLIDS", "A2 (measure the grains), or the SIMPLE cells")
OUT_DECK = os.path.join(WORK, "2_abaqus")

# Where the settings come from. If A12b applied an edit, that edit is what gets built --
# rebuilding from the widgets here is exactly how the CAD viewer's edits used to be
# silently thrown away between "applied" and "written".
_widgets_now = make_params()
if globals().get("EDITED_SETTINGS"):
    from semgrit.editable import params_from_settings as _pfs
    if globals().get("EDITED_BASE") is not None and _widgets_now == EDITED_BASE:
        # Nothing moved underneath it, so build the very object A12b previewed.
        PARAMS = EDITED_PARAMS
        print("settings: the CAD viewer's edits from A12b (%s)"
              % (", ".join(EDITED_CHANGED) or "nothing changed"))
    else:
        # A widget changed after A12b ran. Re-apply the edited fields on top of the
        # current widgets so both survive -- the edit wins where they disagree.
        PARAMS = _pfs(EDITED_SETTINGS, _widgets_now)
        # Keep PLAN describing what is about to be built, so re-opening the viewer or the
        # standoff table after this shows the deck and not the state before the drift.
        from semgrit.build_deck import plan_deck as _plan_deck
        PLAN = _plan_deck(PARAMS, SOLIDS)
        EDITED_PARAMS, EDITED_BASE = PARAMS, _widgets_now
        print("settings: the widgets above, with the CAD viewer's edits from A12b applied")
        print("          on top (%s). The widgets changed after A12b ran, so the summary"
              % ", ".join(sorted(EDITED_SETTINGS)))
        print("          A12b printed is out of date -- what follows is the deck.")
else:
    PARAMS = _widgets_now
    print("settings: the form widgets above")

# Say it before it starts. A cell that sits silent for two minutes reads as hung.
import time as _time
_t0 = _time.time()
print("writing the deck - this is the slow step; a large one takes a minute or two ...")
INFO = build_deck(PARAMS, SOLIDS, OUT_DECK)
print("   done in %.0f s" % (_time.time() - _t0))
print()
R = INFO["outer_radius_mm"]

print("WHEEL   D%g, %s, arc %.4f mm, rim %.4f mm, width %g mm"
      % (2 * R, "FULL" if INFO["full_wheel"] else "%.4f deg" % INFO["resolved_sector_deg"],
         INFO["arc_length_mm"], INFO["rim_depth_mm"], PARAMS.width_mm))
if not INFO["full_wheel"]:
    print("        sagitta %.2f um = %.0f%% of the rim depth -> %s"
          % (INFO["sagitta_um"], 100 * INFO["sagitta_um"] / 1000 / INFO["rim_depth_mm"],
             "reads as an arc" if INFO["sagitta_um"] / 1000 > INFO["rim_depth_mm"]
             else "will look flat; lengthen the arc or thin the rim"))
print("        ONE rigid body: %s shell quads + %s grit facets, ref node %d"
      % (format(INFO["n_bond_shell_quads"], ","), format(INFO["n_grit_facets"], ","),
         INFO["wheel_ref_node"]))
print("GRITS   %s placed" % format(INFO["n_grits"], ","), end="")
if INFO.get("requested_grains"):
    print(" of %s requested" % format(INFO["requested_grains"], ","), end="")
if INFO.get("achieved_areal_density_per_mm2"):
    print("  (%.0f/mm2 achieved)" % INFO["achieved_areal_density_per_mm2"], end="")
print()
if INFO["has_workpiece"]:
    print("        %d can reach the block; tallest reaching protrusion %.4f um"
          % (INFO["n_grits_engaging"], INFO["max_engaging_protrusion_um"]))
    c = INFO["cost"]
    print("WP      %g x %g x %g mm -> %s C3D8R, %d x %d x %d (only deformable part)"
          % (PARAMS.wp_length_mm, PARAMS.wp_width_mm, PARAMS.wp_depth_mm,
             format(INFO["n_workpiece_elements"], ","), *c["element_divisions"]))
    print("        element %.4f cutting x %.4f axial x %.4f depth um; %.4f um sets dt"
          % (c["element_size_cutting_mm"] * 1000, c["element_size_axial_mm"] * 1000,
             c["element_size_depth_mm"] * 1000,
             c["governing_element_size_mm"] * 1000))
    print("        ground face r = %.6f mm, tangent to placement %s, penetration 0"
          % (INFO["workpiece_ground_radius_mm"], INFO["governing_grit_placement_id"]))
    if INFO["workpiece_relocated_to_tallest_grit"]:
        print("        NOTE: moved to theta = %.4f deg - the nominal angle had no grit"
              % INFO["theta_workpiece_deg"])
    print("RUN     dt = %.3e s, omega = %.1f rad/s (%.0f rpm), travel %.4f mm"
          % (c["stable_dt_s"], c["omega_rad_s"], c["rpm"], c["travel_mm"]))
    print("        step %.4e s = %s increments, %.2e element-increments"
          % (c["step_time_s"], format(int(c["increments"]), ","),
             c["element_increments"]))
    print("        estimate " + ", ".join("%s core %.1f h" % (k, v)
          for k, v in sorted(c["est_hours"].items(), key=lambda kv: int(kv[0]))))
else:
    print("        wheel-only deck; position your own ground face at r = %.6f mm"
          % INFO["tallest_tip_whole_arc_mm"])
for m in INFO["warnings"] + INFO["notes"]:
    print("  note: %s" % m)
if INFO.get("run_ready"):
    m = INFO["motion"]
    # Every number here comes off the deck, not off a widget: with DEPTH_OF_CUT_UM = 0
    # the infeed is chosen automatically, and after an edit in A12b the widget is stale.
    print("CUT     ae = %.3f um of infeed at a %.3f um standoff (tallest engaging grain "
          "%.3f um)" % (m["depth_of_cut_mm"] * 1000.0, INFO["clearance_um"],
                        INFO["max_engaging_protrusion_um"]))
    print("        %.1f rad/s = %.0f rpm; V1 = %.3f, V2 = %.3f mm/s inward; VR3 = %.1f;"
          " sweep %.4f mm" % (m["omega_rad_s"], m["rpm"], m["v1"], m["v2"], m["vr3"],
                              m["sweep_mm"]))
print()
print("FILES")
_wrote = [(INFO["path"], "run-ready Abaqus deck - submit from the terminal"
           if INFO.get("run_ready") else "Abaqus deck - geometry only, finish in CAE")]
if INFO.get("cae_deck"):
    _wrote.append((INFO["cae_deck"], "geometry-only twin, same wheel, for CAE"))
_wrote.append((os.path.join(OUT_DECK, MODEL_NAME + "_import_into_cae.py"),
               "run this in CAE: File > Run Script"))
for _k, _d in (("step", "assembled wheel, STEP for SOLIDWORKS"),
               ("stl", "assembled wheel, STL"),
               ("grains_step", "the grits themselves, laid out, STEP")):
    if INFO["cad"].get(_k):
        _wrote.append((os.path.join(OUT_DECK, MODEL_NAME +
                                    ("_grains.step" if _k == "grains_step"
                                     else "." + _k)), _d))
if INFO["cad"].get("grain_stls"):
    _wrote.append((INFO["cad"]["grain_stls"]["dir"],
                   "%d per-grain STL files" % INFO["cad"]["grain_stls"]["count"]))
_wrote.append((os.path.join(OUT_DECK, MODEL_NAME + "_placements.csv"),
               "where every grit ended up"))
_wrote.append((os.path.join(OUT_DECK, MODEL_NAME + "_report.json"),
               "every number this build decided"))
for _p, _d in _wrote:
    _sz = (os.path.getsize(_p) / 1e6) if os.path.isfile(_p) else 0.0
    print("  %-42s %8s  %s" % (os.path.basename(_p),
                               ("%.2f MB" % _sz) if _sz else "dir", _d))
if INFO.get("run_ready"):
    print()
    print("SUBMIT IT WITH:")
    print("  abaqus job=%s input=%s user=<your_vumat>.for double=both cpus=%d interactive"
          % (PARAMS.name, os.path.basename(INFO["path"]), PARAMS.cores))
    print()
    print("  The VUMAT must drive the deletion flag SDV%d to 0 once D reaches 1."
          % N_DEPVAR)
    print("  One that only ever writes 1 deletes nothing, and the result looks ductile.")
'''))

CELLS.append(code('''
#@title 🅰 A14 · Autodesk APS viewer (optional, uses your APS credits) { display-mode: "form" }
#@markdown **This one is billed.** The glTF cell above is free and needs no account;
#@markdown use this only if you specifically want Autodesk's renderer.
#@markdown The built-in 3-D viewer above draws the deck's own triangles and is verified
#@markdown vertex-for-vertex against the `.inp`. **This is an alternative**, not a
#@markdown replacement: it hands the geometry to Autodesk's renderer for nicer shading,
#@markdown section planes and a model tree — at the cost of a cloud round trip.
#@markdown
#@markdown Each view **uploads your model to Autodesk and runs a billed Model Derivative
#@markdown translation**, taking minutes. Use it for a final look, not for iterating.
#@markdown
#@markdown Credentials are a **Client ID and Client Secret** from an app you create at
#@markdown [aps.autodesk.com](https://aps.autodesk.com) (Create App → Custom Integration
#@markdown → enable Model Derivative + Data Management). An Autodesk account email and
#@markdown password will **not** work. They are typed into a password box and never saved.

APS_STEP = "probe only"  #@param ["probe only", "view the STL", "view the STEP"]
#@markdown &nbsp;&nbsp;**Run `probe only` first.** It just checks whether Colab's
#@markdown sandboxed output iframe will load the Autodesk viewer library at all. If that
#@markdown is blocked, nothing else here can work and you have spent nothing.
APS_MAX_UPLOAD_MB = 100.0   #@param {type:"number"}
APS_MAX_BODIES = 2000       #@param {type:"integer"}
#@markdown &nbsp;&nbsp;Hard caps. A multi-body STEP of a dressed wheel can carry
#@markdown thousands of solids and is slow and expensive to translate; the upload is
#@markdown refused above these rather than silently billed. STL is the cheap choice.
APS_BUCKET = ""             #@param {type:"string"}
APS_REGION = "US"           #@param ["US", "EMEA"]

from IPython.display import HTML, display
from semgrit import aps

if APS_STEP == "probe only":
    print("If the line below says LOADED, the sandbox permits the viewer and it is")
    print("worth creating an APS app. If it says BLOCKED, stop here.")
    display(HTML(aps.probe_html()))
else:
    import getpass, os
    _want = ".stl" if APS_STEP == "view the STL" else ".step"
    _f = os.path.join(OUT_DECK, MODEL_NAME + _want)
    if not os.path.exists(_f):
        raise SystemExit("%s was not written. Tick WRITE_WHEEL_%s in the outputs cell "
                         "and rebuild." % (_f, _want[1:].upper()))
    _cfg = aps.APSConfig(
        client_id=getpass.getpass("APS Client ID: ").strip(),
        client_secret=getpass.getpass("APS Client Secret: ").strip(),
        bucket_key=APS_BUCKET, region=APS_REGION,
        max_upload_mb=APS_MAX_UPLOAD_MB, max_bodies=APS_MAX_BODIES)
    APS_RESULT = aps.publish(_f, _cfg)
    display(HTML(aps.viewer_html(APS_RESULT["urn"], APS_RESULT["token"])))
'''))

CELLS.append(code('''
#@title ✅ A15 · Verify the deck — two independent verifiers { display-mode: "form" }
#@markdown Verifier A re-parses the file and re-derives every geometric claim from the
#@markdown node coordinates. Verifier B is a separate implementation that cross-checks
#@markdown the header and report against the mesh and integrates the mass and inertia
#@markdown numerically. Both must pass.
import os
from semgrit.quick import verify_decks

need("INFO", "A13 (build the deck)")
decks = [INFO["path"]]
if INFO.get("cae_deck"):
    decks.append(INFO["cae_deck"])
ok = verify_decks(WORK, decks)
'''))

CELLS.append(code('''
#@title 💾 A16 · Download everything { display-mode: "form" }
import os
from semgrit.quick import bundle as make_bundle

need("OUT_DECK OUT_MEAS", "A13 (build the deck)")
zip_path = make_bundle(WORK, (OUT_DECK, OUT_MEAS), MODEL_NAME)
try:
    from google.colab import files
    files.download(zip_path)
except Exception as exc:
    print("(not on Colab - copy the zip yourself)", exc)
'''))

CELLS.append(md(r"""
---
## 11 · Using the deck

### Run-ready: straight from the terminal

With `RUN_READY` on there is nothing to do in CAE. Put the `.inp` and your VUMAT in the
same folder and submit:

```
abaqus verify -user_explicit
abaqus job=grind input=<name>.inp user=vumat_jh2.for double=both cpus=8 interactive
```

The build cell prints this command with your own names filled in. Two things that will
cost you a run if you get them wrong:

* **The VUMAT filename must have no spaces or brackets.** `vumat (2).for` makes Abaqus
  read `(2).for` as a separate argument and abort.
* **The VUMAT must drive the deletion flag** — `stateNew(km,12) = 0` once `D >= 1`. The
  deck arms `*Depvar, delete=12` and `ELEMENT DELETION=YES` for you, but a VUMAT that
  only ever writes `1` deletes nothing, and the result looks ductile.

If it is interrupted, `RESTART_INTERVALS > 1` lets you resume:

```
abaqus job=grind2 oldjob=grind input=restart.inp user=vumat_jh2.for double=both cpus=8 interactive
```

where `restart.inp` is just `*Heading` followed by `*Restart, read, step=1` — no step
block, which tells Abaqus to finish the interrupted one. Use the **same `cpus`**.

### Geometry only: load it
**File → Run Script…** → `<name>_import_into_cae.py`

Do **not** use File → Import → Part: that reads the `*Part` blocks and skips the
`*Assembly`, so every grain arrives unplaced and the wheel looks bare. File → Import →
Model does not accept `.inp` at all. The script calls `mdb.ModelFromInputFile`, which
reads both.

### What you get
| name | what |
|---|---|
| `WHEEL-1` | one discrete rigid body — bond shell + every grit |
| `A_WHEEL_REF` | its reference node, on the axis at the origin |
| `A_GRITS_SURF` | grit facets only |
| `A_WHEEL_SURF` | grits + the bond's outer face |
| `A_GRITS_ENGAGE_SURF` | just the grits that can reach the block (cheaper contact) |
| `WP-1`, `A_WP_GROUND_SURF` | the workpiece and its ground face |
| `A_WP_BACK_FACE`, `A_WP_SIDE_A/B`, `A_WP_END_A/B` | node sets for fixing it |

### Driving the wheel
One velocity BC on `A_WHEEL_REF` does everything:

```
VR3 = -omega          rad/s   (negative = surface travels toward decreasing theta)
V1 = V2 = V3 = VR1 = VR2 = 0
```

`omega` and the equivalent rpm are printed by the build cell. To cut at depth `ae`,
add the radial infeed — the report gives `theta_workpiece_deg`, and radially inward is
`(-cos θ, -sin θ)`:

```
V1 = -cos(theta) * ae / t_step
V2 = -sin(theta) * ae / t_step
```

**Keep `ae` below the printed bond clearance**, or the bond rim itself hits the
workpiece.

### Step and contact
Use **Dynamic, Explicit**. General contact over the whole model is simplest;
`A_GRITS_ENGAGE_SURF` against `A_WP_GROUND_SURF` is much cheaper on a wheel with
thousands of grits. Element deletion must be **on** in Section Controls if you want
chips to separate.

### Brittle fracture with a JH-2 VUMAT
Two things silently prevent visible brittle fracture, both learned the hard way here:

1. **The VUMAT must drive the deletion flag.** It needs `stateNew(km,12) = 0` once
   `D >= 1`, together with `*Depvar, delete=12` and 12 state variables. A VUMAT that
   only ever writes `stateNew(km,12) = 1` can never delete an element: damaged
   material stays in the mesh carrying residual fractured strength and the result
   looks smeared and ductile instead of cracking.
2. **Bulking pressure must be added in compression only.** Applying the accumulated
   `Δp` in tension too can flip a stretched element to an apparently *compressive*
   pressure, which restores its fractured shear strength — so the crack never opens.

Also make sure element deletion is enabled in Section Controls; the flag alone does
nothing if the section has deletion switched off.

### Mesh size versus the chip
The workpiece element size is the single biggest lever on both cost and whether you see
fracture at all. If a grit takes a 2 µm cut and your elements are 1.5 µm, there is one
element through the chip and no fracture pattern can form. Aim for 5–10 elements through
the deepest cut. Shrinking all three directions together costs **1/h⁴** — 1/h³ more
elements and 1/h more increments.

The element type is fixed at C3D8R, but the size is yours per direction, and the three
are not equally expensive. The stable increment follows the **smallest** element
dimension, so coarsening the **axial** direction alone removes elements for free:

| mesh (cutting × axial × depth, µm) | elements | `dt` | est. 8-core |
|---|---|---|---|
| 0.30 × 0.30 × 0.30 | 160,000 | 6.30e-11 | 1.21 h |
| 0.30 × **1.50** × 0.30 | 32,000 | 6.30e-11 | 0.24 h |
| 0.30 × **1.50** × **0.60** | 16,000 | 6.30e-11 | 0.12 h |

Same cutting-direction resolution, same time increment, **5–10× less work**. Coarsen
cutting or depth only after that, since those blur the chip and the damage zone.

The block keeps the dimensions you ask for, so a size that does not divide them exactly
is rounded to a whole element count — the achieved sizes are printed after the build, and
`dt` is computed from the achieved minimum, not from what you typed.
"""))

# ==========================================================================
# Section B: the single-abrasive ductile/brittle model.
#
# Its own section, after A16, so nothing above it moves. It reuses the grain
# library the A cells (or the SIMPLE cells) already measured and builds its own
# deck, because the whole point is one grit.
# ==========================================================================

CELLS.append(md(r"""
---
# B - Single abrasive: ductile below the critical depth, brittle above it

Everything above builds a wheel and grinds it with **Johnson-Holmquist II**, which is a
brittle law: it damages, bulks and chips whatever the depth of cut. Real grinding does
not work that way. Below a critical depth of cut `dc` the material comes off by plastic
flow -- the ductile regime -- and only above it does it fracture.

This section builds a **single-abrasive** deck that carries both laws in one subroutine,
`vumat_grind.for`, and chooses between them **per material point**:

```
h <  dc  ->  Johnson-Cook + strain-gradient enhancement   (ductile)
h >= dc  ->  Johnson-Holmquist II                         (brittle)
```

### The critical depth of cut

Two published forms, both offered, because they are not interchangeable:

| | |
|---|---|
| `dc = lambda_c (H/E)^0.5 (Kc/H)^2` | the form on this project's slide |
| `dc = lambda_c (E/H) (Kc/H)^2` | Bifano, Dow & Scattergood (1991), whose calibrated `lambda_c` is **0.15** |

They differ by `(E/H)^1.5` -- about **17x** on this sandstone -- so `lambda_c` belongs to
one form and must not be carried over to the other. Say which one you used.

### How the subroutine knows h

A VUMAT is called at one material point and sees no kinematics, so `h` has to be handed
to it. With **one** grit the trajectory is exact, and `h` is a function of the point's
station `u` along the scratch and nothing else:

$$h(u) = H_0 + H_G\,u - \frac{u^2}{2R_{tip}}, \qquad H_G = -\frac{v_r}{\omega R_{tip}}$$

The linear term is the wedge every textbook draws for a grit trajectory -- rubbing, then
ploughing, then shearing -- produced here by the radial infeed rather than by a table
feed. The quadratic term is the sagitta of the grit's circular path: 15 nm across a
48 um block on a D50 mm wheel, which would be ignorable except that `dc` is of that same
order. The classical traverse form `h(theta) = L_g (v_w/v_s) sin(theta)` is the same
straight line over a block far shorter than the contact arc.

`H0`, `HG` and `RTIP` are computed in Python and written into the material card, so the
Fortran carries no process knowledge and stays verifiable on a single material point.
`H0` is pinned to the deck's own tangency -- the grit vertex the block was seated on has
`h` equal to minus the standoff, exactly -- so a sub-micron disagreement about which tip
is tallest cannot leak into a quantity being compared against a few nanometres.

### The strain-gradient term, and why it belongs here

$$\sigma_e = \sigma_{JC}\sqrt{1 + \left(\frac{r' \eta b (M\alpha G)^2}{\sigma_{JC}^2}\right)^{\Lambda}},
\qquad \eta = \frac{4\varepsilon^p_{eq}}{h}$$

As the cut gets thinner the strain gradient rises, geometrically necessary dislocations
accumulate, and the flow stress goes up. That size effect is *why* thin cuts are ductile
at all, so a Johnson-Cook branch without it would understate the ductile regime. At
`Lambda = 1` and `r' = 2` this is simultaneously eq. 7 of the blanking paper, eq. 25 of
the peening paper and eq. 8 of the micro-milling paper: the same Taylor/GND hardening
with a different characteristic length.

### What to plot afterwards

| SDV | |
|---|---|
| **13** | the branch: 1 ductile, 2 brittle -- this is the picture of the transition |
| **14** | `h` at that point |
| **15** | `dc` |
| **19** | the SGE amplification factor, 1 = no size effect |
| 1, 2 | damage and equivalent plastic strain, both branches |
| 12 | STATUS, the deletion flag |
"""))

CELLS.append(code('''
#@title B1 - Single-abrasive settings { display-mode: "form" }
RUN_SINGLE_ABRASIVE = True   #@param {type:"boolean"}
SA_NAME = "single_abrasive_hybrid"  #@param {type:"string"}

#@markdown ### The workpiece material
SA_MATERIAL = "sandstone"  #@param ["sandstone", "silicon_carbide"]
#@markdown &nbsp;&nbsp;Picks the **whole** card together: the 17 JH-2 constants for the
#@markdown brittle branch, the density, the ductile Johnson-Cook constants, and the
#@markdown hardness and toughness `dc` is computed from. Choosing three of those four
#@markdown by hand and forgetting the fourth is how a deck ends up silently mixing two
#@markdown materials, so they move as one.
#@markdown
#@markdown &nbsp;&nbsp;`silicon_carbide` is the **SiC-N** card (rho 3163, G 183 GPa,
#@markdown HEL 14.457 GPa, K1 204.785 GPa, A 0.96, B 0.35, N 0.65, T 0.37 GPa,
#@markdown D1 = D2 = 0.48). It was supplied labelled "monocrystalline silicon", but
#@markdown those are the published silicon **carbide** numbers -- silicon is
#@markdown rho 2329, E ~ 170 GPa, H ~ 11 GPa. Say silicon carbide in the paper.
SA_OVERRIDE_MATERIAL = False  #@param {type:"boolean"}
#@markdown &nbsp;&nbsp;Leave **off** and every Johnson-Cook / SGE / damage number below
#@markdown is taken from the material above and the fields are ignored. Turn it **on**
#@markdown to hand-enter them instead -- which is what you want once you have your own
#@markdown calibration. The cell prints which source it used either way.

#@markdown ### The abrasive and the block
SA_DIAMETER_MM = 50.0      #@param {type:"number"}
SA_SLICE_MM = 2.0          #@param {type:"number"}
SA_GRAIN_INDEX = -1        #@param {type:"integer"}
#@markdown &nbsp;&nbsp;`-1` picks the largest grain in the measured library.
SA_GRIT_OFFSET_MM = 0.015  #@param {type:"number"}
#@markdown &nbsp;&nbsp;Where the grit starts along the block. Positive puts it at the
#@markdown entry end, so the default rotation drags it across the whole workpiece.
SA_WP_LENGTH_MM = 0.048    #@param {type:"number"}
SA_WP_WIDTH_MM = 0.015     #@param {type:"number"}
SA_WP_DEPTH_MM = 0.006     #@param {type:"number"}
SA_ELEMENT_UM = 0.30       #@param {type:"number"}
SA_ELEMENT_AXIAL_UM = 0.0  #@param {type:"number"}
SA_ELEMENT_DEPTH_UM = 0.0  #@param {type:"number"}
SA_SURFACE_LAYER_UM = 0.0  #@param {type:"number"}
SA_DEPTH_GROWTH = 1.3      #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = use the base size. A **graded** depth mesh --
#@markdown `SA_ELEMENT_DEPTH_UM` elements for the first `SA_SURFACE_LAYER_UM`, then
#@markdown growing by `SA_DEPTH_GROWTH` -- is what lets a nanometre-scale cut be
#@markdown resolved affordably. The stable increment follows the smallest dimension
#@markdown either way, so keep the axial size coarse to pay for it.
SA_STANDOFF_UM = 0.0       #@param {type:"number"}
SA_DEPTH_OF_CUT_UM = 0.0   #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = automatic: close the standoff, then cut 85% of the grain
#@markdown protrusion. The depth of cut is what decides how much of the scratch is
#@markdown brittle, so this is the knob to sweep.

#@markdown ### The transition
SA_DC_NM = 0.0             #@param {type:"number"}
#@markdown &nbsp;&nbsp;`0` = compute it from the hardness and toughness below.
SA_DC_FORM = "Bifano: lambda_c (E/H) (Kc/H)^2"  #@param ["Bifano: lambda_c (E/H) (Kc/H)^2", "lambda_c (H/E)^0.5 (Kc/H)^2"]
#@markdown &nbsp;&nbsp;Bifano is the default because it is the form with a published
#@markdown calibrated `lambda_c`, and it is what the `RUN_ME*` packages use. The other
#@markdown form gives 5.3 nm on sandstone and 0.7 nm on SiC -- below anything a mesh can
#@markdown resolve, so every element comes out brittle and the switch shows nothing.
SA_LAMBDA_C = 0.15         #@param {type:"number"}
SA_HARDNESS_MPA = 1000.0   #@param {type:"number"}
SA_KIC_MPA_SQRT_M = 0.30   #@param {type:"number"}
#@markdown &nbsp;&nbsp;Toughness in the usual **MPa*sqrt(m)**; it is converted to the
#@markdown deck's MPa*sqrt(mm) for you. That conversion is a factor of 31.6, and getting
#@markdown it wrong scales `dc` by 1000. These two are used only when
#@markdown `SA_OVERRIDE_MATERIAL` is on; otherwise the material's own values are used
#@markdown (sandstone 1000 MPa / 0.30, SiC 25000 MPa / 3.5).
SA_SWITCH = "on dc"        #@param ["on dc", "force ductile everywhere", "force brittle everywhere"]
#@markdown &nbsp;&nbsp;The two overrides run the same deck as pure JC+SGE or pure JH-2,
#@markdown so you can see how much of a result the switch itself caused.

#@markdown ### To reproduce a clearly visible transition
#@markdown &nbsp;&nbsp;The defaults above are a general-purpose deck. These are the
#@markdown settings `RUN_ME/1_single_abrasive` and `RUN_ME_SIC/1_single_abrasive`
#@markdown ship, found by sweeping the depth of cut until both regimes were tens of
#@markdown elements wide:
#@markdown
#@markdown | | sandstone | silicon carbide |
#@markdown |---|---|---|
#@markdown | `SA_DC_FORM` | Bifano | Bifano |
#@markdown | `dc` | 87.75 nm | 52.92 nm |
#@markdown | `SA_DEPTH_OF_CUT_UM` | **0.40** | **0.36** |
#@markdown | `SA_ELEMENT_UM` | 0.30 | 0.30 |
#@markdown | `SA_ELEMENT_AXIAL_UM` | 1.5 | 1.5 |
#@markdown | `SA_ELEMENT_DEPTH_UM` | 0.03 | 0.03 |
#@markdown | `SA_SURFACE_LAYER_UM` | 0.45 | 0.45 |
#@markdown | `SA_DEPTH_GROWTH` | 1.45 | 1.45 |
#@markdown | `SA_WP_WIDTH_MM` | 0.009 | 0.009 |
#@markdown | transition lands at | u = +0.0042 mm | u = +0.0081 mm |
#@markdown | wall clock, 8 cores | ~0.22 h | ~1.54 h |
#@markdown
#@markdown &nbsp;&nbsp;SiC is 7x slower on the same mesh because its wave speed is
#@markdown 1.23e7 mm/s against sandstone's 1.76e6, so the stable increment is 7x
#@markdown smaller. Refining the mesh is not what costs the time and coarsening it
#@markdown will not buy it back.

#@markdown ### Ductile branch - Johnson-Cook
#@markdown &nbsp;&nbsp;`A` defaults to this material's own JH-2 quasi-static uniaxial
#@markdown compressive strength, so the two laws meet at the transition instead of
#@markdown stepping across it. The rest are PLACEHOLDERS: calibrate them.
SA_JC_A_MPA = 90.0         #@param {type:"number"}
SA_JC_B_MPA = 50.0         #@param {type:"number"}
SA_JC_N = 0.50             #@param {type:"number"}
SA_JC_C = 0.020            #@param {type:"number"}
SA_JC_M = 1.0              #@param {type:"number"}
SA_E_MPA = 6500.0          #@param {type:"number"}
SA_NU = 0.21               #@param {type:"number"}
#@markdown &nbsp;&nbsp;6500 MPa and 0.21 are exactly the JH-2 card's `K1` and `G`, so
#@markdown both branches share one elasticity and the stable increment is unambiguous.
SA_DENSITY_KG_M3 = 2350.0  #@param {type:"number"}
SA_CP_J_KGK = 800.0        #@param {type:"number"}
SA_TMELT_K = 1473.15       #@param {type:"number"}

#@markdown ### Strain-gradient enhancement
SA_BURGERS_NM = 0.50       #@param {type:"number"}
SA_TAYLOR_M = 3.0          #@param {type:"number"}
SA_ALPHA = 0.30            #@param {type:"number"}
SA_LAMBDA_SGE = 1.0        #@param {type:"number"}
SA_R_PRIME = 2.0           #@param {type:"number"}

#@markdown ### Ductile damage - Johnson-Cook
SA_D1 = 0.0                #@param {type:"number"}
SA_D2 = 0.15               #@param {type:"number"}
SA_D3 = -1.5               #@param {type:"number"}
SA_D4 = 0.0                #@param {type:"number"}
SA_D5 = 0.0                #@param {type:"number"}
SA_DCRIT = 1.0             #@param {type:"number"}
print("single-abrasive settings captured")
'''))

CELLS.append(code('''
#@title B2 - What the switch will do. Nothing is written { display-mode: "form" }
#@markdown Computes `dc`, the chip-thickness field and where along the scratch the
#@markdown transition lands, from the same placement code the writer uses. Change
#@markdown anything in B1 and re-run this until the split looks like the experiment you
#@markdown are modelling.
import dataclasses, math, os

import matplotlib.pyplot as plt

from semgrit import materials
from semgrit.analysis import wheel_motion
from semgrit.build_deck import hybrid_single_grit, plan_deck
from semgrit.hybrid import (HybridParams, kic_from_mpa_sqrt_m, plan_hybrid)
from semgrit.hybrid import summary_text as hybrid_summary
from semgrit_multi.plot import trajectory_figure

need("SOLIDS", "A2 (measure the grains), or the SIMPLE cells")

SA_H_SOURCE = {"on dc": 0, "force ductile everywhere": 2,
               "force brittle everywhere": 3}[SA_SWITCH]
SA_FORM = 2 if SA_DC_FORM.startswith("Bifano") else 1

print(materials.summary_text(SA_MATERIAL))
print()
if SA_OVERRIDE_MATERIAL:
    print("SA_OVERRIDE_MATERIAL is ON: the ductile constants come from the B1")
    print("fields, NOT from the material above. The JH-2 card still does.")
    SA_HP = HybridParams(
        enabled=True,
        a_mpa=SA_JC_A_MPA, b_mpa=SA_JC_B_MPA, n=SA_JC_N, c=SA_JC_C, m=SA_JC_M,
        youngs_mpa=SA_E_MPA, poisson=SA_NU,
        density_kg_m3=SA_DENSITY_KG_M3, specific_heat_j_kgk=SA_CP_J_KGK,
        tmelt_k=SA_TMELT_K,
        burgers_mm=SA_BURGERS_NM * 1e-6, taylor_factor=SA_TAYLOR_M,
        alpha=SA_ALPHA, sge_exponent=SA_LAMBDA_SGE, r_prime=SA_R_PRIME,
        d1=SA_D1, d2=SA_D2, d3=SA_D3, d4=SA_D4, d5=SA_D5, dcrit=SA_DCRIT,
        dc_mm=SA_DC_NM * 1e-6, lambda_c=SA_LAMBDA_C,
        hardness_mpa=SA_HARDNESS_MPA,
        kic=kic_from_mpa_sqrt_m(SA_KIC_MPA_SQRT_M),
        dc_form=SA_FORM, h_source=SA_H_SOURCE)
else:
    print("ductile constants taken from the material card above. Set")
    print("SA_OVERRIDE_MATERIAL = True in B1 to hand-enter them instead.")
    SA_HP = materials.hybrid_params(
        SA_MATERIAL, dc_form=SA_FORM, h_source=SA_H_SOURCE,
        dc_mm=SA_DC_NM * 1e-6, lambda_c=SA_LAMBDA_C)

SA_PARAMS = hybrid_single_grit(
    hybrid=SA_HP, name=SA_NAME,
    diameter_mm=SA_DIAMETER_MM, arc_length_mm=SA_SLICE_MM,
    single_grain_index=SA_GRAIN_INDEX, single_grit_offset_mm=SA_GRIT_OFFSET_MM,
    wp_length_mm=SA_WP_LENGTH_MM, wp_width_mm=SA_WP_WIDTH_MM,
    wp_depth_mm=SA_WP_DEPTH_MM, wp_element_size_mm=SA_ELEMENT_UM / 1000.0,
    wp_element_size_width_mm=SA_ELEMENT_AXIAL_UM / 1000.0,
    wp_element_size_depth_mm=SA_ELEMENT_DEPTH_UM / 1000.0,
    wp_surface_layer_mm=SA_SURFACE_LAYER_UM / 1000.0,
    wp_depth_growth=SA_DEPTH_GROWTH,
    clearance_um=SA_STANDOFF_UM)
SA_PARAMS.analysis.depth_of_cut_um = SA_DEPTH_OF_CUT_UM
# Moves the JH-2 card, the density and the *Material name together. Without it
# the brittle branch would stay on whatever material the preset shipped with.
materials.apply(SA_PARAMS, SA_MATERIAL,
                check_hybrid=not SA_OVERRIDE_MATERIAL)

if RUN_SINGLE_ABRASIVE:
    SA_PLAN = plan_deck(SA_PARAMS, SOLIDS)
    SA_FIELD, SA_DC = plan_hybrid(SA_PLAN, SA_HP)
    print("=" * 78)
    print("ONE ABRASIVE on a D%g wheel, %s C3D8R elements in the block"
          % (SA_DIAMETER_MM, format(SA_PLAN["n_workpiece_elements"], ",")))
    print("depth of cut %.4f um, standoff %.4f um, grain protrusion %.4f um"
          % (SA_PLAN["depth_of_cut_um"], SA_STANDOFF_UM,
             SA_PLAN["protrusion_um"]["max"]))
    print("=" * 78)
    print(hybrid_summary(SA_FIELD, SA_DC, SA_PLAN["_wp"], SA_HP))

    # # transition visuals (single abrasive)
    # The picture the numbers above describe: where the grit goes, how deep it
    # is at each station, and where that crosses dc. The depth of cut the
    # transition happens at is read straight off the vertical axis.
    _st = float((SA_PLAN.get("cost") or {}).get("step_time_s") or 0.0)
    _an = dataclasses.replace(SA_PARAMS.analysis,
                              depth_of_cut_um=float(SA_PLAN["depth_of_cut_um"]))
    _mot = wheel_motion(_an, SA_PLAN["_place"]["theta_c"],
                        SA_PARAMS.surface_speed_mm_s,
                        SA_PARAMS.outer_radius_mm, _st)
    fig = trajectory_figure(SA_PLAN["_place"], _mot, SA_PLAN["_wp"], SA_DC,
                            step_time_s=_st,
                            rotation_reversed=bool(
                                SA_PARAMS.analysis.rotation_reversed))
    plt.show()

    print()
    print("-" * 78)
    print("Nothing written. Sweep SA_DEPTH_OF_CUT_UM to move the transition;")
    print("when the split is what you want, run B3.")
    print("-" * 78)
else:
    print("single-abrasive section skipped")
'''))

CELLS.append(code('''
#@title B3 - Build it, verify it, download it { display-mode: "form" }
#@markdown Writes the deck, copies `vumat_grind.for` next to it, and runs three gates:
#@markdown the two deck verifiers plus `verify_hybrid_deck.py`, which is the only one
#@markdown that checks the deck and the subroutine agree about which elements are
#@markdown ductile. On Colab run `!apt-get -qq install gfortran` first if you want that
#@markdown last check to compile the subroutine rather than skip that part of itself.
SA_BUILD = True            #@param {type:"boolean"}
SA_DOWNLOAD = True         #@param {type:"boolean"}

import os, shutil, subprocess, sys, time

from semgrit.build_deck import build_deck
from semgrit.quick import bundle, verify_decks

if SA_BUILD and RUN_SINGLE_ABRASIVE:
    need("SA_PARAMS SOLIDS", "B2")
    SA_OUT = os.path.join(WORK, "3_single_abrasive")
    _t0 = time.time()
    print("writing %s ..." % SA_NAME)
    SA_INFO = build_deck(SA_PARAMS, SOLIDS, SA_OUT)
    _hy = SA_INFO["hybrid"]
    print("wrote %s  (%.1f MB, %.0f s)"
          % (os.path.basename(SA_INFO["path"]), SA_INFO["size_bytes"] / 1e6,
             time.time() - _t0))
    print("dc = %.4f nm; transition at u = %s"
          % (_hy["dc_nm"], _hy["chip_field"]["transition_u_mm"]))
    for _m in SA_INFO["warnings"] + SA_INFO["notes"]:
        print("  note: %s" % _m)

    # The subroutine has to travel with the deck, or the deck cannot be run.
    for _f in ("vumat_grind.for", "vumat_jh2.for"):
        if os.path.exists(os.path.join(WORK, _f)):
            shutil.copy(os.path.join(WORK, _f), os.path.join(SA_OUT, _f))

    print()
    if not verify_decks(WORK, [SA_INFO["path"]]):
        raise SystemExit("the deck did not verify - read the FAIL lines above")

    print()
    print("#" * 78)
    print("# verify_hybrid_deck.py - does the subroutine agree with the card?")
    print("#" * 78)
    _r = subprocess.run([sys.executable,
                         os.path.join(WORK, "verify_hybrid_deck.py"),
                         SA_INFO["path"]],
                        capture_output=True, text=True, cwd=WORK)
    print(_r.stdout[-6000:])
    if _r.stderr.strip():
        print(_r.stderr[-2000:])
    if _r.returncode != 0:
        raise SystemExit("the hybrid gate failed - do not run this deck")

    print()
    print("to run it:")
    print("  abaqus job=%s input=%s user=vumat_grind.for double=both cpus=%d"
          % (SA_NAME, os.path.basename(SA_INFO["path"]), SA_PARAMS.cores))
    if SA_INFO.get("postprocess_script"):
        print("to read the result:")
        print("  abaqus python %s %s.odb"
              % (os.path.basename(SA_INFO["postprocess_script"]), SA_NAME))
    if SA_DOWNLOAD:
        _zip = bundle(WORK, (SA_OUT,), SA_NAME)
        try:
            from google.colab import files
            files.download(_zip)
        except Exception as exc:
            print("(not on Colab - copy the zip yourself)", exc)
elif not RUN_SINGLE_ABRASIVE:
    print("single-abrasive section is off")
else:
    print("not built - tick SA_BUILD when B2 shows the split you want")
'''))

CELLS.append(md(r"""
---
## B4 - Running the single-abrasive deck, and reading it

```
abaqus verify -user_explicit
abaqus job=grind input=single_abrasive_hybrid.inp user=vumat_grind.for double=both cpus=8 interactive
```

`double=both` is not optional here: the chip thickness is compared against a threshold of
a few nanometres on a wheel of 25 mm, a ratio of 1e-7, and single precision does not have
the digits.

### The first three things to look at

1. **SDV13** on the ground face. It should be `1` (ductile) over the run-in and `2`
   (brittle) after the transition station the build printed. If it is uniform, the depth
   of cut never crossed `dc` -- that is a real answer, not a bug, but check that it is
   the answer you meant.
2. **SDV19**, the SGE amplification. It is `1` where the size effect does nothing and
   rises where the cut is thin. If it is `1` everywhere then `h` is far above the Burgers
   vector everywhere and the gradient term is inert.
3. **RF at `A_WHEEL_REF`**, which is the grinding force. The post-processing script
   written with the deck reads exactly this.

### Sweeping the transition

The cleanest experiment this deck supports is a depth-of-cut sweep: build three decks with
`SA_DEPTH_OF_CUT_UM` well below, near, and well above the value that puts the transition
mid-block, then compare the force traces and the chip morphology. Because the geometry is
identical between them -- same wheel, same grain, same seating, same mesh -- any
difference is the constitutive law and nothing else.

Two more runs worth having, from the same deck:

* `SA_SWITCH = "force brittle everywhere"` reproduces plain JH-2, and
  `verify_vumat_grind.py` proves that path is bit-identical to `vumat_jh2.for`;
* `SA_SWITCH = "force ductile everywhere"` is plain JC+SGE.

Those two bracket the hybrid result, so together they say how much of it the switch
caused.

### Honest limits

* The **Johnson-Cook constants are placeholders.** `A` is tied to the JH-2 card's own
  quasi-static compressive strength so the two branches meet, but `B, n, C, m` and
  `D1..D5` are order-of-magnitude values for a quartz-bonded rock. Calibrate them against
  nanoindentation or single-scratch data before quoting a force.
* **`lambda_c` is a calibration**, and it belongs to whichever `dc` form you chose. The
  two forms differ by `(E/H)^1.5`.
* **`h` is prescribed, not measured.** It comes from the grit trajectory this deck writes,
  which is exact for one grit and a constant radial infeed. It is not valid for many
  grits interacting, nor for a wheel whose grits wear during the run.
* The switch is **latched at the first increment** and does not migrate: a material point
  keeps the law its station implies for the whole run.
* Because `h` is compared with `dc` once per point, the transition is a sharp line
  between neighbouring elements. That is a bimaterial interface rather than a
  discontinuity inside an element, but it does mean the mesh has to be fine enough that
  the line lands where you want it.
"""))

nb = {
    "nbformat": 4, "nbformat_minor": 0,
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True,
                  "name": "SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb"},
        "kernelspec": {"name": "python3", "display_name": "Python 3"},
        "language_info": {"name": "python"},
    },
    "cells": CELLS,
}

# Splice the blob in as real Python string literals, one per line, so the cell is
# valid source. Substituting into the serialised JSON instead drops the quotes.
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

OUT = 'SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb'
with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
    json.dump(nb, fh, indent=1)

# Every code cell must be valid Python once its lines are concatenated.
import ast
for n, cell in enumerate(json.load(open(OUT, encoding='utf-8'))["cells"], 1):
    if cell["cell_type"] == "code":
        body = ''.join(cell["source"])
        # Colab form directives are plain comments, so this is real Python.
        ast.parse(body, filename='cell %d' % n)

check = json.load(open(OUT, encoding='utf-8'))
print('%s  %.0f KB  %d cells (%d code, %d markdown)'
      % (OUT, os.path.getsize(OUT) / 1024, len(check['cells']),
         sum(c['cell_type'] == 'code' for c in check['cells']),
         sum(c['cell_type'] == 'markdown' for c in check['cells'])))
