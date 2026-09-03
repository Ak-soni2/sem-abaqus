"""The extra cells that bring SEM_TO_ABAQUS_SAG.ipynb to feature parity.

Kept in a separate module because _make_notebook_sag.py was already long, and
because these cells are the ones a reader of the paper will actually drive:
the CAD viewer, the mesh viewer, the abrasive-height and grinding-theory
readouts, the Plotly and glTF views, and the downloads.

The CAD viewer needs a PLACED plan -- bond, grits, workpiece, boundary
conditions -- which ``sagdeck.plan`` does not produce: it plans the compliant
two-scale model, not a rigid wheel. So cell A1 below builds a rigid-wheel plan
of the SAME tool geometry (125 mm diameter, the pad's own measured areal
density, the SAG depth of cut) purely so the viewer has something to show. That
is stated in the cell rather than hidden, because the viewed body is a
*visualisation* of the pad, not the deck that gets solved.
"""

CAD = '''
#@title A1 - A viewable model of the pad { display-mode: "form" }
#@markdown The CAD viewer draws a **placed** model: bond, grains, workpiece and
#@markdown every boundary condition the deck writes. `sagdeck.plan` does not
#@markdown produce one -- it plans the compliant two-scale model, where the
#@markdown "bond" is a hyperelastic ring and the grain count is in the hundreds
#@markdown of thousands.
#@markdown
#@markdown So this cell builds a rigid-wheel plan of the **same tool geometry**
#@markdown -- your wheel diameter, the pad's measured areal density, the SAG
#@markdown depth of cut -- so the viewer has real placed grains to show. It is a
#@markdown **visualisation of the pad**, not the deck that gets solved. The
#@markdown decks come from cell 6.
CAD_ARC_MM = 1.0        #@param {type:"number"}
CAD_WIDTH_MM = 0.30     #@param {type:"number"}
CAD_RIM_DEPTH_MM = 0.05 #@param {type:"number"}
need("PLAN SOLIDS", "cells 3 and 4")
from semgrit import materials as _materials
from semgrit.analysis import AnalysisParams
from semgrit.build_deck import DeckParams, plan_deck

_c = PLAN["contact"]
_dens = _c.active_grains / max(_c.spot_area_mm2, 1e-12)
CAD_PARAMS = DeckParams(
    name="sag_pad_view", diameter_mm=P.diameter_mm,
    include_bond=True, include_workpiece=True,
    sector_mode="arc", arc_length_mm=CAD_ARC_MM,
    rim_depth_mm=CAD_RIM_DEPTH_MM, width_mm=CAD_WIDTH_MM,
    grit_mode="areal_density", areal_density_per_mm2=_dens,
    wp_length_mm=CAD_ARC_MM * 0.2, wp_width_mm=CAD_WIDTH_MM * 0.7,
    wp_depth_mm=max(20.0 * PLAN["material"]["dc_nm"] * 1e-6, 0.005),
    wp_element_size_length_mm=CAD_ARC_MM / 100.0,
    wp_element_size_width_mm=CAD_WIDTH_MM / 100.0,
    wp_element_size_depth_mm=PLAN["micro"]["element_mm"],
    clearance_um=0.0, wp_position="centred",
    surface_speed_mm_s=_c.surface_speed_mm_s, cores=P.cores,
    analysis=AnalysisParams(
        enabled=True, depth_of_cut_um=_c.indentation_nm * 1e-3,
        material_model="hybrid",
        hybrid=_materials.hybrid_params(P.material, h_source=0, dc_form=2)))
_materials.apply(CAD_PARAMS, P.material)
CAD_PLAN = plan_deck(CAD_PARAMS, SOLIDS)
print("a viewable pad: %s grains placed on a %.0f mm tool"
      % (format(CAD_PLAN["n_grits"], ","), P.diameter_mm))
print("pad density   %.0f grains/mm2 (from the contact solution)" % _dens)
print("depth of cut  %.4f um (the per-grain indentation)"
      % (_c.indentation_nm * 1e-3))
print("")
print("This is for VIEWING. The solved decks come from cell 6.")
'''

CADVIEW = '''
#@title A2 - CAD viewer: the state-of-the-art one { display-mode: "form" }
#@markdown The same three.js viewer the main notebook uses, on the SAG pad.
#@markdown
#@markdown | | |
#@markdown |---|---|
#@markdown | **Shaded with edges** | feature edges over a lit surface |
#@markdown | **Wheel / Contact** | the whole 125 mm tool, or the grains on the work |
#@markdown | **Face / Axial** | straight at the pad, or down the tool axis |
#@markdown | **Section plane** | cut on any axis and drag through the model |
#@markdown | **Click a grain** | id, protrusion, height, width, volume, position |
#@markdown | **Shift-click twice** | distance and X Y Z, plus radial / along-arc / across-face |
#@markdown | **Parts tree** | show or hide the pad, the grains, the workpiece |
#@markdown | **Boundary conditions** | every symbol stands for a keyword the deck really writes |
#@markdown | **Drag block** (`G`) | drag the workpiece along the arc, shift-drag for standoff |
#@markdown | **Depth-of-cut band** | the valid window, shaded green |
#@markdown | **Colour the grains by** | protrusion, height, width, volume, or engages-the-block |
#@markdown | **Explode** | pull pad, grains and work apart along the radius |
#@markdown | **Cap the cut face** | a solid face instead of a hollow shell |
#@markdown | **Fullscreen**, **Save PNG**, **Keyboard** (`?`) | 12 shortcuts |
#@markdown
#@markdown No account, no upload. three.js loads from a CDN; the model is
#@markdown embedded in the page.
SHOW_CAD = True          #@param {type:"boolean"}
CAD_MODE = "whole wheel" #@param ["whole wheel", "wheel", "contact"]
CAD_HEIGHT = 720         #@param {type:"integer"}
CAD_MAX_INLINE_MB = 24.0 #@param {type:"number"}
need("CAD_PLAN", "cell A1")
from IPython.display import HTML, display
from semgrit.cadviewer import build as build_cad_view

if SHOW_CAD:
    _html, _meta, _info = build_cad_view(
        CAD_PLAN, os.path.join(WORK, "sag_pad.glb"), mode=CAD_MODE,
        max_grits=0, height=CAD_HEIGHT, max_inline_mb=CAD_MAX_INLINE_MB)
    print("%s: %s triangles, %d of %d grains drawn (%d in full detail)"
          % (CAD_MODE, format(_info["triangles"], ","), _meta["grits_drawn"],
             _meta["grits_total"], _meta["grits_full_detail"]))
    for _n in _meta.get("notes", []):
        print("note:", _n)
    display(HTML(_html))
else:
    print("set SHOW_CAD to draw the pad.")
'''

MESHVIEW = '''
#@title A3 - Mesh viewer: see what will actually be solved { display-mode: "form" }
#@markdown The CAD view above is the *geometry*. This is the **mesh** -- and the
#@markdown mesh is where the arguments are.
#@markdown
#@markdown | question | how you answer it here |
#@markdown |---|---|
#@markdown | Is $d_c$ actually resolved? | the element edges are drawn; count them through the surface band |
#@markdown | Can the compliant layer **bend**? | a layer with too few elements through its thickness only shears |
#@markdown | Is anything inverted? | inverted elements are **refused**, not drawn -- Abaqus reports this as a cryptic preprocessing failure with no element numbers |
#@markdown | Is the grading where it should be? | section the block and look at the depth transition |
#@markdown
#@markdown It is the *same viewer*, fed element geometry instead of solids, so
#@markdown it keeps section capping, explode, the measuring tool and every
#@markdown shortcut. The panel is retitled for a mesh -- "click an element face"
#@markdown rather than "click a grain".
SHOW_MESH = True       #@param {type:"boolean"}
MESH_PART = "all"      #@param ["all", "tool only", "workpiece only"]
MESH_EDGES = True      #@param {type:"boolean"}
MESH_HEIGHT = 700      #@param {type:"integer"}
need("PLAN", "cell 4")
from IPython.display import HTML, display
from semgrit import meshview as _mv
from semgrit.sagwrite import build_block, build_compliant_ring

if SHOW_MESH:
    _r_out = 0.5 * P.diameter_mm
    _r_in = _r_out - P.polyurethane.thickness_mm
    _sect = min(PLAN["macro"]["sector_deg"], 30.0)
    _mic = PLAN["micro"]
    _meshes = []
    if MESH_PART in ("all", "tool only"):
        _hub = build_compliant_ring(
            inner_r_mm=max(_r_in - 2.5, 1.0), outer_r_mm=_r_in,
            width_mm=P.width_mm, sector_deg=_sect,
            n_circ=28, n_rad=2, n_axial=6)
        _pu = build_compliant_ring(
            inner_r_mm=_r_in, outer_r_mm=_r_out, width_mm=P.width_mm,
            sector_deg=_sect, n_circ=28, n_rad=6, n_axial=6)
        _meshes += [
            dict(name="hub (rigid)", nodes=_hub[0], conn=_hub[1],
                 color=_mv.C_HUB),
            dict(name="polyurethane %0.1f mm" % P.polyurethane.thickness_mm,
                 nodes=_pu[0], conn=_pu[1], color=_mv.C_COMPLIANT)]
    if MESH_PART in ("all", "workpiece only"):
        _wp = build_block(
            length_mm=_mic["side_mm"], width_mm=_mic["side_mm"],
            depth_mm=_mic["depth_mm"],
            el_length_mm=_mic["element_inplane_mm"],
            el_width_mm=_mic["element_inplane_mm"],
            fine_depth_mm=_mic["element_mm"],
            band_mm=_mic["depth_mm"] * 0.5, growth=1.3,
            x0_mm=-0.5 * _mic["side_mm"], y0_mm=-0.5 * _mic["side_mm"])
        _meshes.append(dict(name="workpiece (MICRO, dc/%g)"
                            % P.elements_per_dc,
                            nodes=_wp[0], conn=_wp[1], color=_mv.C_WORK))

    _h, _m, _i = _mv.build(_meshes, os.path.join(WORK, "sag_mesh.glb"),
                           height=MESH_HEIGHT, edges=MESH_EDGES)
    print("%-34s %10s %10s %9s %s"
          % ("part", "elements", "min edge", "aspect", "inverted"))
    for _k, _v in _m["stats"].items():
        print("%-34s %10s %9.4f nm %8.1f:1 %8d"
              % (_k[:34], format(_v["elements"], ","),
                 _v["min_edge"] * 1e6, _v["aspect_max"], _v["inverted"]))
    print("")
    print("dc = %.1f nm, surface element %.2f nm -> %.1f elements across dc"
          % (PLAN["material"]["dc_nm"], _mic["element_mm"] * 1e6,
             PLAN["material"]["dc_nm"] / (_mic["element_mm"] * 1e6)))
    for _n in _m["notes"]:
        print("note:", _n)
    display(HTML(_h))
else:
    print("set SHOW_MESH to draw the mesh.")
'''

HEIGHTS = '''
#@title A4 - Abrasive heights, and what the pad can reach { display-mode: "form" }
#@markdown A grit cuts only as deep as it stands proud of its backing. On a
#@markdown rigid wheel that sets a hard ceiling on the depth of cut. On a SAG
#@markdown pad it matters for a different reason: the indentation is *tiny*
#@markdown against the grain, so the pad is nowhere near its geometric limit --
#@markdown and this cell shows by how much.
need("SOLIDS PLAN", "cells 3 and 4")
import numpy as _np

_h = _np.array([s.height_um for s in SOLIDS])
_c = PLAN["contact"]
_dc = PLAN["material"]["dc_nm"]
print("measured grain heights, %d solids" % len(_h))
for _q in (0, 5, 25, 50, 75, 95, 100):
    print("   %3d%%  %8.3f um" % (_q, _np.percentile(_h, _q)))
print("")
print("the pad's nominal grain size   %8.3f um" % P.grain_um)
print("mean measured height           %8.3f um" % _h.mean())
print("")
print("indentation this process makes %8.5f um  (%.3f nm)"
      % (_c.indentation_nm * 1e-3, _c.indentation_nm))
print("as a fraction of a mean grain  %8.2e" % (_c.indentation_nm * 1e-3
                                                / _h.mean()))
print("as a multiple of dc            %8.5f  (dc = %.1f nm)"
      % (_c.indentation_nm / _dc, _dc))
print("")
if _c.indentation_nm * 1e-3 < 0.01 * _h.mean():
    print("The grain is >100x deeper than the cut, so protrusion is NOT the")
    print("limit here -- which is exactly what makes SAG a finishing process")
    print("rather than a stock-removal one.")
else:
    print("The cut is a significant fraction of the grain height: check that")
    print("the pad is not being asked to cut deeper than it protrudes.")
'''

THEORY = '''
#@title A5 - Is this a real finishing regime? { display-mode: "form" }
#@markdown The deck can be geometrically perfect and still describe a process
#@markdown nobody would call grinding. These are the first questions a reviewer
#@markdown asks, and verifying the `.inp` answers none of them.
#@markdown
#@markdown **measured** rows are counted off the contact solution. **theory**
#@markdown rows are the textbook expressions for an equivalent traverse grind,
#@markdown so they need a work speed; with `WORK_SPEED_MM_MIN = 0` they are
#@markdown reported as not applicable rather than quietly computed from zero.
WORK_SPEED_MM_MIN = 15.0   #@param {type:"number"}
need("PLAN", "cell 4")
import math as _math

_c = PLAN["contact"]
_dc = PLAN["material"]["dc_nm"]
_R = 0.5 * P.diameter_mm
print("MEASURED, off the contact solution")
print("  normal load FN            %10.4f N" % _c.normal_load_n)
print("  tangential FT             %10.4f N" % (P.friction
                                                * _c.normal_load_n))
print("  spot area As              %10.2f mm2" % _c.spot_area_mm2)
print("  spot length Ls            %10.3f mm" % (2 * _c.semi_axis_a_mm))
print("  mean pressure             %10.5f MPa" % _c.mean_pressure_mpa)
print("  active grains             %10s" % format(int(_c.active_grains), ","))
print("  load per grain Fn         %10.4e N" % _c.load_per_grain_n)
print("  indentation d             %10.4f nm" % _c.indentation_nm)
print("  groove width              %10.1f nm" % _c.groove_width_nm)
print("  surface speed vs          %10.1f mm/s" % _c.surface_speed_mm_s)
print("  grain crossings / rev     %10s" % format(int(_c.grains_per_rev), ","))
print("  MRR                       %10.4f mm3/min" % _c.mrr_mm3_min)
print("")
_vw = float(WORK_SPEED_MM_MIN) / 60.0
if _vw > 0:
    print("THEORY, for an equivalent traverse grind at %.1f mm/min"
          % WORK_SPEED_MM_MIN)
    _ae = _c.indentation_nm * 1e-6
    print("  contact length sqrt(ae*de)%10.4f mm"
          % _math.sqrt(max(_ae, 0) * P.diameter_mm))
    print("  equivalent chip h_eq      %10.4e mm"
          % (_ae * _vw / max(_c.surface_speed_mm_s, 1e-9)))
    print("  speed ratio vs/vw         %10.0f"
          % (_c.surface_speed_mm_s / _vw))
    print("  removal rate Q'w          %10.4e mm3/s per mm" % (_ae * _vw))
else:
    print("THEORY: not applicable -- set WORK_SPEED_MM_MIN to compare with a")
    print("traverse grind. This is a plunge/spot configuration, and the")
    print("chip-thickness formulas need a work speed to mean anything.")
print("")
print("FINDINGS")
_bad = []
if _c.indentation_nm >= _dc:
    _bad.append("the indentation already exceeds dc, so removal is brittle "
                "from the first pass")
if _c.face_overrun > 1.0:
    _bad.append("the elliptical patch is %.1f%% wider than the %.0f mm face, "
                "so it is clipped by the wheel edges (%.1f%% of the nominal "
                "area is off the wheel)"
                % (100.0 * (_c.face_overrun - 1.0), P.width_mm,
                   100.0 * _c.area_clipped_fraction))
if not _c.density_measured:
    _bad.append("the pad density is interpolated, not measured for this "
                "grain size")
if PLAN["infeasible"]:
    _bad += list(PLAN["infeasible"])
if _bad:
    for _b in _bad:
        print("  - %s" % _b)
else:
    print("  nothing to flag: the regime is self-consistent.")
'''

PLOTLY = '''
#@title A6 - Quick 3-D scatter of the pad (Plotly) { display-mode: "form" }
#@markdown Every placed grain as a point, sized by protrusion. Cheaper than the
#@markdown CAD viewer and useful for seeing the *distribution* rather than the
#@markdown geometry -- whether the pad is uniform, whether the seeding clumped.
SHOW_SCATTER = True   #@param {type:"boolean"}
need("CAD_PLAN", "cell A1")
if SHOW_SCATTER:
    try:
        import plotly.graph_objects as _go
    except ImportError:
        import subprocess as _sp
        _sp.run([sys.executable, "-m", "pip", "-q", "install", "plotly"],
                check=True)
        import plotly.graph_objects as _go
    # The placement objects are on the model, not under plan["_place"] --
    # that key is a dict of per-plan arrays (baked vertices, frames, the
    # engaged set), which is a different thing entirely.
    _pl = CAD_PLAN["_model"].placements
    _x = [q.translation_mm[0] for q in _pl]
    _y = [q.translation_mm[1] for q in _pl]
    _z = [q.translation_mm[2] for q in _pl]
    _pr = [q.protrusion_mm * 1000.0 for q in _pl]
    _fig = _go.Figure(_go.Scatter3d(
        x=_x, y=_y, z=_z, mode="markers",
        marker=dict(size=3, color=_pr, colorscale="Viridis",
                    colorbar=dict(title="protrusion (um)"), opacity=0.85),
        text=["grain %d: %.2f um proud" % (i, p)
              for i, p in enumerate(_pr)]))
    _fig.update_layout(height=620, margin=dict(l=0, r=0, t=28, b=0),
                       title="%s grains on the pad, coloured by protrusion"
                             % format(len(_pl), ","),
                       scene=dict(aspectmode="data"))
    _fig.show()
else:
    print("set SHOW_SCATTER to draw it.")
'''

DOWNLOAD = '''
#@title A7 - Download everything { display-mode: "form" }
#@markdown Bundles the decks, the reports, the run scripts and the figures into
#@markdown one archive. On Colab it downloads; elsewhere it just says where the
#@markdown file is.
WHAT = "decks and reports"  #@param ["decks and reports", "everything in the output folder"]
need("MICRO", "cell 6")
import glob as _glob
import shutil as _shutil
import tarfile as _tf

_out = os.path.dirname(MICRO["path"]) or "."
_arc = os.path.join(WORK, "sag_bundle.tar.gz")
_pats = ["*.inp", "*.json", "*.csv", "*.for", "*.bat", "*.sh", "*.md",
         "*.png"] if WHAT == "decks and reports" else ["*"]
with _tf.open(_arc, "w:gz") as _t:
    _n = 0
    for _p in _pats:
        for _f in sorted(_glob.glob(os.path.join(_out, "**", _p),
                                    recursive=True)):
            if os.path.isfile(_f):
                _t.add(_f, arcname=os.path.relpath(_f, _out))
                _n += 1
print("%d file(s), %.1f MB -> %s" % (_n, os.path.getsize(_arc) / 1e6, _arc))
try:
    from google.colab import files as _files
    _files.download(_arc)
except Exception:
    print("(not Colab: copy the file from the path above)")
'''
