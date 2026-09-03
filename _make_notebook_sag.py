"""Generate SEM_TO_ABAQUS_SAG.ipynb -- shape adaptive grinding.

A third notebook. The other two model a RIGID wheel: the tool is three orders
stiffer than the cut, so its deformation is not the physics and it is emitted
as a discrete rigid surface. Shape adaptive grinding inverts that. The tool is
a compliant polyurethane layer, its deformation IS the process, and everything
downstream follows from it -- line contact spreads into an area, the load is
shared by every grain that area covers, the force per grain collapses by orders
of magnitude, and a brittle material can then be removed ductilely.

So this notebook cannot reuse the rigid writer, and it does not. It carries
``semgrit.sag`` (the contact chain), ``sagdeck`` (the two-scale planner),
``sagwrite`` / ``sagemit`` (the deformable-tool decks) and ``meshview`` (the
mesh in the CAD viewer), all of which are new.

It also inverts the transition criterion. ``vumat_grind.for`` compares a
PRESCRIBED chip thickness h(u) against dc, which needs a known trajectory; with
a compliant tool and hundreds of thousands of grains there is none, and the
load per grain is not knowable in advance -- it is what the contact solution
produces. ``vumat_grind2.for``'s local criterion needs no geometry at all:

    W_p * L_c >= PSI * Kc^2 / E

accumulated plastic work against a fracture energy, triggering on HISTORY. With
PSI left at 0 the subroutine derives dc*E*H/Kc^2, so the threshold is exactly
W_p*L_c >= H*dc and a MEASURED dc carries straight through.

    python _make_notebook_sag.py
"""
import base64
import glob
import gzip
import io
import json
import os
import tarfile

FILES = (sorted(glob.glob('semgrit/*.py'))
         + sorted(glob.glob('semgrit_multi/*.py'))
         + ['verify_rigid_deck.py',
            'verify_rigid_deck2.py',
            'verify_hybrid_deck.py',
            'verify_vumat_grind.py',
            'verify_vumat_grind2.py',
            'verify_sag_deck.py',
            '_make_sag_paper.py',
            '_make_sag_packages.py',
            '_derive_grind2.py',
            'vumat_grind.for',
            'vumat_grind2.for',
            'vumat_jh2.for',
            '_hybrid_test/driver.f',
            '_hybrid_test/vaba_param.inc'])


assert 'semgrit/sagfig.py' in [f.replace(os.sep, '/') for f in FILES], \
    'semgrit/sagfig.py is missing from the payload'
assert 'semgrit/meshview.py' in [f.replace(os.sep, '/') for f in FILES], \
    'semgrit/meshview.py is missing from the payload'
assert 'semgrit/sagemit.py' in [f.replace(os.sep, '/') for f in FILES], \
    'semgrit/sagemit.py is missing from the payload'


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


import _sagcells as _extra

CELLS = []

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
# SEM image → shape adaptive grinding, ductile and brittle

**Shape adaptive grinding (SAG)** replaces the rigid wheel with a *compliant*
one: a stiff hub, a polyurethane layer a few millimetres thick, and an abrasive
pad on the outside. Press it against the work and the layer squashes, so line
contact spreads into an **area**.

That single fact is the whole process. The contact load is shared by every
grain the patch covers — hundreds of thousands of them — so the force on each
one collapses to $10^{-5}$ N and the depth each takes collapses with it. A
material that fractures under a conventional wheel can then be removed by
**plastic flow**, which is how a brittle cermet reaches a 21 nm finish.

### What this notebook does

You give it SEM micrographs of your abrasive. It measures every grain,
reconstructs each as a 3-D solid, solves the compliant contact, and writes two
Abaqus decks:

| deck | question it answers | resolves $d_c$? |
|---|---|---|
| **MACRO** | the *contact* — patch size, pressure, engaged grains, load per grain | no |
| **MICRO** | the *transition* — SDV13, ductile against brittle | **yes**, at $d_c/5$ |

They are coupled by one number: the per-grain load MACRO computes is what MICRO
applies. Both decks print it, so the pair cannot be quoted out of step.

### How the transition is decided

The other two notebooks in this project compare a *prescribed* chip thickness
$h(u)$ against $d_c$. That needs a known trajectory. Here there isn't one — with
a compliant tool the load per grain is the *answer*, not an input. So SAG uses
the **local energy criterion**:

$$W_p \cdot L_c \;\ge\; \Psi\,\frac{K_c^2}{E}$$

accumulated plastic work per unit volume, times the element's own length,
against a fracture energy. It needs no geometry, and it triggers on **history**
— a point starts ductile and turns brittle as work accumulates under repeated
grain passes, which is what a polishing pad physically does.

With $\Psi = 0$ the subroutine derives $\Psi = d_c E H/K_c^2$, making the
threshold exactly $W_p L_c \ge H d_c$. So a **measured** $d_c$ carries straight
through with no new calibration.

> **One property to know before quoting a result.** The criterion is
> regularised by $L_c$, so it is mesh-dependent *by construction*: halving the
> element halves the work density needed to trigger. That is correct for a
> fracture-energy criterion, and it means $\Psi$ is calibrated **for a mesh**.
> Every deck states its element size. Cell 10 measures the sensitivity.

### Reference

Ghosh, Sidpara & Bandyopadhyay (2021), *Brittle-ductile transition in compliant
finishing of HVOF sprayed hard WC-Co coating*, Int. J. Refractory Metals and
Hard Materials **99**, 105610. The contact chain in cell 4 is that paper's
eqs. 1–16, and cell 11 rebuilds its experiment.
"""))

# ---------------------------------------------------------------------------
CELLS.append(code('''
#@title 1 - Setup: unpack the pipeline (run once) { display-mode: "form" }
# semgrit (including the four new SAG modules), semgrit_multi, both VUMATs and
# every gate are embedded below, so this notebook is self-contained.
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

from semgrit import sag as _sag
from semgrit import sagdeck as _sd
from semgrit import sagemit as _se       # noqa: F401
from semgrit import meshview as _mv      # noqa: F401
print("pipeline ready in", WORK)
print("SAG modules : sag (contact), sagdeck (planner), sagwrite + sagemit")
print("              (deformable-tool decks), meshview (mesh in the viewer)")
print("subroutine  : vumat_grind2.for -- 58 constants, energy criterion")


def need(names, where):
    """Stop with the cell to run, instead of a NameError on a stray name."""
    absent = [n for n in names.split() if n not in globals()]
    if absent:
        raise SystemExit("run %s first - this cell needs %s"
                         % (where, ", ".join(absent)))
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 2 · Your abrasive pad, under the microscope

SAG pads are characterised by two numbers the contact model needs: the **grain
size** $d_g$ and the **areal density** $C_0$ of grains on the pad. Both come
from SEM micrographs of the pad itself.

Upload your own images, or leave the default to use the B4C micrographs
embedded in this notebook.
"""))

CELLS.append(code('''
#@title 2 - Where are your SEM images? { display-mode: "form" }
SOURCE = "bundled"  #@param ["bundled", "upload", "google drive", "already on disk"]
IMAGE_PATH = "/content/drive/MyDrive/sem/*.tif"  #@param {type:"string"}
PIXEL_SIZE_UM = 0.0  #@param {type:"number"}
#@markdown `PIXEL_SIZE_UM = 0` reads the scale from the SEM databar.
import glob, os

if SOURCE == "upload":
    from google.colab import files
    up = files.upload()
    IMAGES = sorted(os.path.join(os.getcwd(), n) for n in up)
elif SOURCE == "google drive":
    from google.colab import drive
    drive.mount("/content/drive")
    IMAGES = sorted(glob.glob(IMAGE_PATH))
elif SOURCE == "bundled":
    IMAGES = sorted(glob.glob(os.path.join(WORK, "B4C_1*.tif")))
    if not IMAGES:
        raise SystemExit("no bundled images found; choose 'upload' instead")
else:
    IMAGES = sorted(glob.glob(IMAGE_PATH))

if not IMAGES:
    raise SystemExit("no images matched %r" % IMAGE_PATH)
print("%d image(s):" % len(IMAGES))
for p in IMAGES:
    print("   ", os.path.basename(p))
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 3 · Measure the grains

Every grain is segmented, measured (25 shape descriptors), and reconstructed as
a watertight 3-D solid whose maximum projected cross-section **is** the measured
outline. The figures below show every stage, so nothing is taken on trust.
"""))

CELLS.append(code('''
#@title 3 - Measure every grain, and show the work { display-mode: "form" }
SHOW_STAGES = True   #@param {type:"boolean"}
need("IMAGES", "cell 2")
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams["figure.dpi"] = 110
_show = plt.show

from semgrit import figures as figs
from semgrit.quick import measure_images

MEAS = measure_images(IMAGES, os.path.join(WORK, "_sag_meas"),
                      pixel_size_um=(PIXEL_SIZE_UM or None),
                      keep_stages=SHOW_STAGES, log=print)
SOLIDS = MEAS["solids"]
GRAINS = MEAS["grains"]
print("")
print("%d grain solids from %d image(s)" % (len(SOLIDS), len(IMAGES)))
hs = [s.height_um for s in SOLIDS]
print("heights %.2f to %.2f um (mean %.2f)"
      % (min(hs), max(hs), sum(hs) / len(hs)))

if SHOW_STAGES and MEAS.get("per_image"):
    rec = MEAS["per_image"][0]
    for fn in (figs.calibration, figs.segmentation_stages,
               figs.segmentation_overlay, figs.outline_fidelity,
               figs.solid_verification):
        try:
            fn(rec)
            _show()
        except Exception as exc:
            print("(%s skipped: %s)" % (fn.__name__, exc))
    figs.measurement_distributions(GRAINS)
    _show()
    figs.grain_gallery(SOLIDS)
    _show()
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 4 · The compliant contact

Now the SAG-specific physics, following the reference paper's eqs. 1–16.

The tool is pressed in by the **wheel compression** $T$, and Hertz gives the
load:

$$F_N = 1.44\,E_{eq}\,R^{1/2}\,T^{3/2}
\qquad
E_{eq} = \left(\frac{1-\nu_w^2}{E_w} + \frac{1-\nu_t^2}{E_t}\right)^{-1}$$

The patch area and length are empirical fits to measured finishing spots:

$$A_s = 138.22\,T^{0.151}N^{0.009}
\qquad
L_s = 17.69\,T^{0.232}N^{0.012}$$

The load is then divided among the grains the patch covers, and each grain's
indentation follows from the Brinell relation:

$$N_{abr} = C_a A_s
\qquad
F_n = \frac{F_N}{N_{abr}}
\qquad
d = \frac{d_g}{2} - \tfrac{1}{2}\sqrt{d_g^2 - d_i^2}$$

**Set your process here.** Everything downstream — patch size, per-grain load,
mesh, deck size, runtime — follows from these numbers.
"""))

CELLS.append(code('''
#@title 4 - Your SAG process { display-mode: "form" }
#@markdown ### The tool
WHEEL_DIAMETER_MM = 125.0   #@param {type:"number"}
WHEEL_WIDTH_MM = 10.0       #@param {type:"number"}
LAYER_THICKNESS_MM = 5.0    #@param {type:"number"}
#@markdown Polyurethane, neo-Hookean. `E = 6*C10`, so C10 = 0.16606 is ~1.0 MPa.
PU_C10_MPA = 0.16606        #@param {type:"number"}
PU_DENSITY_KG_M3 = 1100.0   #@param {type:"number"}
PU_PRONY_G = 0.11           #@param {type:"number"}
PU_PRONY_TAU_S = 0.01       #@param {type:"number"}

#@markdown ### The process
COMPRESSION_MM = 0.4        #@param {type:"number"}
SPEED_RPM = 1050.0          #@param {type:"number"}
FRICTION = 0.2              #@param {type:"number"}
GRAIN_UM = 6.0              #@param [6.0, 15.0, 30.0] {type:"raw", allow-input: true}
#@markdown Pad density in grains/mm2. 0 uses the measured value for 6/15/30 um.
PAD_DENSITY_PER_MM2 = 0.0   #@param {type:"number"}

#@markdown ### The workpiece
MATERIAL = "wc_co"          #@param ["wc_co", "silicon_carbide", "sandstone"]
CARBIDE_UM = 1.36           #@param {type:"number"}
BHN_KGF_MM2 = 581.0         #@param {type:"number"}

#@markdown ### Resolution and cost
ELEMENTS_PER_DC = 5.0       #@param {type:"number"}
MICRO_GRAINS = 1            #@param {type:"integer"}
MACRO_SECTOR_MODE = "contact"  #@param ["contact", "cap"]
MACRO_GRAIN_CAP = 400000    #@param {type:"integer"}
CORES = 8                   #@param {type:"integer"}

need("SOLIDS", "cell 3")
from semgrit.sagdeck import Polyurethane, SAGParams, plan

PU = Polyurethane(c10_mpa=PU_C10_MPA, density_kg_m3=PU_DENSITY_KG_M3,
                  prony_g=PU_PRONY_G, prony_tau_s=PU_PRONY_TAU_S,
                  thickness_mm=LAYER_THICKNESS_MM)
P = SAGParams(
    diameter_mm=WHEEL_DIAMETER_MM, width_mm=WHEEL_WIDTH_MM,
    polyurethane=PU, use_shore_modulus=False,
    compression_mm=COMPRESSION_MM, speed_rpm=SPEED_RPM, friction=FRICTION,
    grain_um=float(GRAIN_UM),
    pad_areal_per_mm2=PAD_DENSITY_PER_MM2,
    material=MATERIAL, carbide_um=CARBIDE_UM, bhn_kgf_mm2=BHN_KGF_MM2,
    elements_per_dc=ELEMENTS_PER_DC, micro_grains=MICRO_GRAINS,
    macro_sector_mode=MACRO_SECTOR_MODE, macro_grain_cap=MACRO_GRAIN_CAP,
    cores=CORES, name="sag_%gum" % float(GRAIN_UM))
PLAN = plan(P)
C = PLAN["contact"]

print(chr(10).join(_sd.macro_header(PLAN)))
print("")
print(chr(10).join(_sd.micro_header(PLAN)))
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 5 · The contact, in pictures

Four things worth seeing rather than reading:

1. **Why SAG works at all** — the per-grain load against wheel compression, for
   all three pads. The collapse is the process.
2. **The patch**, with its Hertzian pressure distribution.
3. **$d_c$ three ways** — the two published geometric forms and the energy
   criterion differ by orders of magnitude on the same material, which is why
   the deck records which one it used.
4. **The regime map** — where this operating point sits relative to $d_c$.
"""))

CELLS.append(code('''
#@title 5 - The contact, drawn { display-mode: "form" }
need("PLAN", "cell 4")
import numpy as np
from semgrit import sagfig

for fn in (sagfig.load_collapse, sagfig.contact_patch,
           sagfig.dc_comparison, sagfig.regime_map):
    fn(PLAN)
    _show()
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 6 · Write the decks

Two decks, both `*Dynamic, Explicit` with **general contact**.

General contact is required here, not merely convenient, for three independent
reasons: the VUMAT **deletes elements**, and deletion exposes interior faces
that a pre-declared contact pair would never see (a chip would separate and
then pass through the tool); **which grains touch is the answer**, so it cannot
be declared in advance; and a compliant layer at high compression can fold onto
**itself**.

The MACRO deck runs three steps, and the first two are timed by the layer's own
physics rather than chosen:

| step | what it does | why that duration |
|---|---|---|
| **PRESS** | push in by $T$ | slow enough that $v/c = 0.005$ in the layer — a fast ramp loads the patch *inertially* and its pressure is not the steady Hertzian one |
| **HOLD** | dwell | $3\tau$, so the polyurethane relaxes to its **long-term** modulus, which is the state a load-cell reading and the Hertz comparison both correspond to |
| **GRIND** | rotate | the process |
"""))

CELLS.append(code('''
#@title 6 - Write MACRO and MICRO { display-mode: "form" }
WRITE_MACRO = False   #@param {type:"boolean"}
#@markdown MACRO carries the full pad, so it is ~150 MB. MICRO is the deck that
#@markdown answers the transition; leave MACRO off unless you want the contact.
OUTDIR = "RUN_SAG_NB"  #@param {type:"string"}
need("PLAN SOLIDS", "cells 3 and 4")
import os
from semgrit import sagemit

os.makedirs(OUTDIR, exist_ok=True)
MICRO = sagemit.write_micro(os.path.join(OUTDIR, "micro.inp"), PLAN, SOLIDS)
print("MICRO  %s" % MICRO["path"])
print("  %s elements, %.1f nm depth element, %.2f MB"
      % (format(MICRO["elements"], ","), MICRO["element_depth_mm"] * 1e6,
         MICRO["bytes"] / 1e6))
print("  %d passes over one track, driven by %.4e N per grain"
      % (MICRO["n_passes"], MICRO["load_per_grain_n"]))
print("  energy threshold W_p*L_c >= %.4f MPa*mm = %.1f J/m2"
      % (MICRO["energy_threshold_mpa_mm"],
         MICRO["energy_threshold_mpa_mm"] * 1000.0))
print("  dc = %.1f nm (%s)"
      % (MICRO["dc_nm"], "MEASURED" if MICRO["dc_measured"] else "computed"))

MACRO = None
if WRITE_MACRO:
    MACRO = sagemit.write_macro(os.path.join(OUTDIR, "macro.inp"), PLAN,
                                SOLIDS)
    print("")
    print("MACRO  %s" % MACRO["path"])
    print("  %s elements (%s PU, %s work), %s grains, %.1f MB"
          % (format(MACRO["elements"], ","),
             format(MACRO["pu_elements"], ","),
             format(MACRO["work_elements"], ","),
             format(MACRO["grains"], ","), MACRO["bytes"] / 1e6))
    print("  sector %.3f deg, press %.1f mm/s (v/c = %.4f)"
          % (MACRO["sector_deg"], MACRO["press_velocity_mm_s"],
             PLAN["timing"]["press_mach"]))
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 7 · Look at it — CAD, mesh, and the numbers behind both

Everything above is arithmetic. This section is where you check it by eye, and
it is the same viewer the main notebook uses — not a reduced one.

| cell | what it shows |
|---|---|
| **A1** | a *viewable* placed model of the pad |
| **A2** | the **CAD viewer** — section planes, click-to-inspect, boundary conditions, explode, colour-by-property, 12 shortcuts |
| **A3** | the **mesh viewer** — element edges, quality per part, inverted elements refused |
| **A4** | abrasive heights against the depth this process actually cuts |
| **A5** | is this a real finishing regime? measured against textbook |
| **A6** | the pad's grain distribution, as a 3-D scatter |
| **A7** | download the lot |

> **A1 needs saying plainly.** The CAD viewer draws a *placed* model — bond,
> grains, workpiece, boundary conditions. The SAG planner does not produce one:
> its "bond" is a hyperelastic ring and its grain count runs to hundreds of
> thousands. So A1 builds a rigid-wheel plan of the **same tool geometry** —
> your diameter, the pad's own measured density, the SAG depth of cut — purely
> so there is something to inspect. It is a **visualisation of the pad**, not
> the deck that gets solved. The solved decks come from cell 6.
"""))

CELLS.append(code(_extra.CAD))
CELLS.append(code(_extra.CADVIEW))
CELLS.append(code(_extra.MESHVIEW))
CELLS.append(code(_extra.HEIGHTS))
CELLS.append(code(_extra.THEORY))
CELLS.append(code(_extra.PLOTLY))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 8 · A compact mesh preview

The same viewer, fed two different things.

**The CAD** is the geometry the deck describes. **The mesh** is where the
arguments are: whether $d_c$ is actually resolved, whether the compliant layer
has enough elements through its thickness to *bend* rather than merely shear,
whether anything is inverted. Element edges are drawn, and inverted elements
are refused rather than displayed — a viewer is the last place a human looks
before submitting a multi-day job, so it is the right place to stop a mesh that
cannot run.
"""))

CELLS.append(code('''
#@title 8 - Compact mesh preview { display-mode: "form" }
SHOW = "mesh"  #@param ["mesh", "cad"]
DRAW_EDGES = True  #@param {type:"boolean"}
need("PLAN", "cell 4")
from IPython.display import HTML, display

if SHOW == "mesh":
    from semgrit import meshview as mv
    from semgrit.sagwrite import build_block, build_compliant_ring
    p = PLAN["params"]
    r_out = 0.5 * p.diameter_mm
    r_in = r_out - p.polyurethane.thickness_mm
    sect = min(PLAN["macro"]["sector_deg"], 30.0)
    hub = build_compliant_ring(inner_r_mm=r_in - 2.5, outer_r_mm=r_in,
                               width_mm=p.width_mm, sector_deg=sect,
                               n_circ=24, n_rad=2, n_axial=6)
    pu = build_compliant_ring(inner_r_mm=r_in, outer_r_mm=r_out,
                              width_mm=p.width_mm, sector_deg=sect,
                              n_circ=24, n_rad=6, n_axial=6)
    mic = PLAN["micro"]
    wp = build_block(length_mm=mic["side_mm"], width_mm=mic["side_mm"],
                     depth_mm=mic["depth_mm"],
                     el_length_mm=mic["element_inplane_mm"],
                     el_width_mm=mic["element_inplane_mm"],
                     fine_depth_mm=mic["element_mm"],
                     band_mm=mic["depth_mm"] * 0.5, growth=1.3,
                     x0_mm=-0.5 * mic["side_mm"],
                     y0_mm=-0.5 * mic["side_mm"])
    html, meta, info = mv.build(
        [dict(name="hub", nodes=hub[0], conn=hub[1], color=mv.C_HUB),
         dict(name="polyurethane", nodes=pu[0], conn=pu[1],
              color=mv.C_COMPLIANT),
         dict(name="workpiece (MICRO)", nodes=wp[0], conn=wp[1],
              color=mv.C_WORK)],
        os.path.join(WORK, "_sagmesh.glb"), height=680, edges=DRAW_EDGES)
    for k, v in meta["stats"].items():
        print("%-20s %8s elements, aspect max %6.1f:1, inverted %d"
              % (k, format(v["elements"], ","), v["aspect_max"],
                 v["inverted"]))
    for n in meta["notes"]:
        print("note:", n)
    display(HTML(html))
else:
    print("The CAD view needs a placed rigid-wheel plan; SAG's tool is")
    print("deformable, so the mesh view above IS the model. Use the")
    print("grinding-wheel notebook for the rigid-wheel CAD.")
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 9 · Verify the deck

`verify_sag_deck.py` shares **no code** with the writer. It re-parses the
`.inp` text with its own keyword-grammar reader, re-measures the node
coordinates, recomputes every hex Jacobian, and re-interprets all 58 material
constants — so a bug in the writer cannot also be baked into its own verifier.

Among the things it checks: the energy threshold recomputed from the card must
equal $H d_c$; **Bifano's $d_c$ computed from that same card** must differ, to
catch a deck that quietly fell back on the 17×-too-large value; the press must
be a *velocity* whose product with the step time equals the compression; and
the passes must **alternate direction**, because a one-way slide leaves every
point with a single pass and could never accumulate to the threshold.
"""))

CELLS.append(code('''
#@title 9 - Verify, independently { display-mode: "form" }
need("MICRO", "cell 6")
import subprocess, sys

args = [sys.executable, "verify_sag_deck.py", MICRO["path"], "--no-converge"]
if MACRO:
    args.insert(3, MACRO["path"])
r = subprocess.run(args, capture_output=True, text=True)
print(r.stdout[-9000:])
if r.stderr.strip():
    print("stderr:", r.stderr[-2000:])
print("exit code", r.returncode,
      "-- 0 means every check passed" if r.returncode == 0 else "-- SEE ABOVE")
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 10 · Mesh convergence — read this before quoting a number

The energy criterion is regularised by the element length, so it is
**mesh-dependent by construction**. Halving the element halves the work
*density* needed to trigger.

That is not a defect; it is what an energy-based failure criterion does. The
quantity the criterion actually tests, $W_p \cdot L_c$, is mesh-*independent* —
and the cell below verifies that to $10^{-16}$ while the density it corresponds
to changes fourfold.

The consequence for a paper: **$\Psi$ is calibrated for a mesh**, and any
transition depth quoted from this model has to be quoted with the element size
that produced it.
"""))

CELLS.append(code('''
#@title 10 - How much does the mesh move the answer? { display-mode: "form" }
need("PLAN", "cell 4")
import subprocess, sys
r = subprocess.run([sys.executable, "-c",
                    "import verify_sag_deck as v; v.converge()"],
                   capture_output=True, text=True)
print(r.stdout)
if r.stderr.strip():
    print("stderr:", r.stderr[-1500:])
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 11 · Rebuild the reference paper

Everything above is your process. This cell rebuilds the *paper's* experiment —
all three pads at its best operating point — so the model can be tested against
a published result.

**One parameter is calibrated, and it is worth knowing which.** The paper gives
eq. 4 for the backing pad's modulus from its shore hardness, but never prints
the shore hardness. Two independent routes exist: a hand-built CAE deck for this
process carries C10 = 0.0575 MPa ($E$ = 0.345 MPa), and inverting the contact
chain for the modulus that reproduces the paper's *stated* per-grain forces
gives 0.43 MPa. Those agree to 25 % — a real corroboration.

Pinning it tighter uses the paper's headline result (6 µm pad, pure ductile,
60–100 nm chips) together with its 30 µm force ceiling, which leaves
**C10 = 0.16606 MPa**. Only that value satisfies both constraints.

### What this can and cannot test

| testable against the paper | |
|---|---|
| contact mechanics — groove width, per-grain force, $k$ ratio | **yes**, and they land in its bands |
| **transition ordering** — 30 µm brittle → 6 µm ductile | **yes. This is the test.** |
| force magnitudes | **no** — the WC-Co Johnson-Cook constants are placeholders except $A$ |
| surface roughness $S_a$ | **no** — needs ~20 000 grain crossings against the 11–24 simulated |

**SDV13, the branch map, is the result.** Everything else is diagnostic.
"""))

CELLS.append(code('''
#@title 11 - Build the paper's three decks { display-mode: "form" }
BUILD_PAPER = False  #@param {type:"boolean"}
PADS = "all"  #@param ["all", "6 um only", "30 um only"]
#@markdown Also write run.bat / run.sh / postprocessor / EXPECTED.md per folder.
MAKE_PACKAGES = True  #@param {type:"boolean"}
import subprocess, sys

if not BUILD_PAPER:
    print("Set BUILD_PAPER to see the calibration and build the decks.")
    print("Showing the calibration only:")
    r = subprocess.run([sys.executable, "_make_sag_paper.py", "--compare"],
                       capture_output=True, text=True)
    print(r.stdout)
else:
    args = [sys.executable, "_make_sag_paper.py"]
    if PADS == "all":
        args.append("--all")
    elif PADS == "30 um only":
        args += ["--all"]
    r = subprocess.run(args, capture_output=True, text=True)
    print(r.stdout[-6000:])
    if r.stderr.strip():
        print("stderr:", r.stderr[-1500:])
    if MAKE_PACKAGES and r.returncode == 0:
        q = subprocess.run([sys.executable, "_make_sag_packages.py"],
                           capture_output=True, text=True)
        print(q.stdout[-3000:])
'''))

# ---------------------------------------------------------------------------
CELLS.append(md(r"""
## 12 · Running the deck, and reading the result

```
abaqus job=micro input=micro.inp user=vumat_grind2.for double=both cpus=8 interactive
```

Three things about that command line are not optional.

**`double=both`.** $h$ and $d_c$ are compared at 80 nm against a millimetre
geometry — a ratio of $10^{-6}$. Single precision has ~7 decimal digits and
does not have them. The failure is **silent**: the branch flag comes out wrong
and the job does not crash.

**`vumat_grind2.for`, not `vumat_grind.for`.** This deck carries 58 constants
and the energy criterion; the other subroutine reads 56 and would misinterpret
the card.

**A datacheck first.** `cpus=1 datacheck` takes seconds and reads every keyword
and the material card. The one real submission this project ever made died
exactly there, on a `*User Material` card written four values to a line instead
of eight.

### What to plot

**SDV13 is the result**: 1 = ductile, 2 = brittle.

Plot it **after every pass**, not only at the end. The criterion accumulates,
so *when* a point flips is the physics — and it is what distinguishes the three
pads from each other.

| SDV | meaning |
|---|---|
| **13** | **branch: 1 ductile, 2 brittle** |
| 14 | the chip thickness the point was given |
| 15 | $d_c$ actually used |
| 19 | strain-gradient amplification |
| 12 | deletion flag |
| 21, 22 | the energy criterion's own accumulators |

### Before quoting a force

The Johnson-Cook constants for both WC-Co and SiC are **placeholders** except
$A$, which is derived from the JH-2 card's own quasi-static compressive
strength so the two branches meet at the transition. $B, n, C, m$ and
$D_1..D_5$ are defensible orders of magnitude and nothing more.

**The branch map is the result; the force magnitudes are not**, until those are
calibrated against nanoindentation or scratch data on your own material.
"""))

CELLS.append(code(_extra.DOWNLOAD))

# ---------------------------------------------------------------------------
nb = {"cells": CELLS,
      "metadata": {"colab": {"provenance": [], "toc_visible": True},
                   "kernelspec": {"display_name": "Python 3",
                                  "name": "python3"},
                   "language_info": {"name": "python"}},
      "nbformat": 4, "nbformat_minor": 0}

blob = payload()
CHUNK = 3000
chunks = ['    "%s"\n' % blob[i:i + CHUNK] for i in range(0, len(blob), CHUNK)]
spliced = 0
for cell in nb["cells"]:
    if cell["cell_type"] != "code" or "__PAYLOAD__\n" not in cell["source"]:
        continue
    i = cell["source"].index("__PAYLOAD__\n")
    cell["source"][i:i + 1] = chunks
    spliced += 1
assert spliced == 1, 'payload marker found %d times' % spliced

OUT = 'SEM_TO_ABAQUS_SAG.ipynb'
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
print('%d files embedded, payload %.0f KB' % (len(FILES), len(blob) / 1024))
for other in ('SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb',
              'SEM_TO_ABAQUS_MULTI_ABRASIVE.ipynb',
              'SEM_TO_ABAQUS_PRESENTATION.ipynb'):
    print('  %-42s %s' % (other,
                          'untouched' if os.path.exists(other) else 'MISSING'))
