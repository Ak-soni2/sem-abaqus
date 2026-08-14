"""One-shot patch: add the transition visuals and the custom-trajectory cell.

Both notebook generators get the same two things:

* the trajectory figure -- tip depth against station with dc drawn across it,
  so the depth of cut the transition happens at can be read straight off the
  axis, plus the same paths from above;
* for the multi-abrasive notebook, a depth-of-cut sweep curve and a cell that
  loads a MEASURED trajectory from a CSV or an image and replays it.

Kept as a file rather than a shell heredoc because the inserted text contains
triple quotes of both kinds. Run once; it refuses to run twice.
"""
import io
import os
import sys

Q = "'''"
CO = "CELLS.append(code(" + Q
MARK1 = "# transition visuals (single abrasive)"
MARK2 = "# transition visuals + measured trajectory (multi abrasive)"


def edit(text, anchor, replacement, path):
    n = text.count(anchor)
    if n != 1:
        raise SystemExit("%s: anchor found %d times:\n%r"
                         % (path, n, anchor[:80]))
    return text.replace(anchor, replacement)


# --------------------------------------------------------------------------
# 1. the single-abrasive notebook: payload + a trajectory figure in B2
# --------------------------------------------------------------------------

NB1_FILES_OLD = """FILES = sorted(glob.glob('semgrit/*.py')) + [
    'verify_rigid_deck.py',"""
NB1_FILES_NEW = """FILES = (sorted(glob.glob('semgrit/*.py'))
         + sorted(glob.glob('semgrit_multi/*.py'))) + [
    'verify_rigid_deck.py',"""

NB1_PLOT_OLD = """    print(hybrid_summary(SA_FIELD, SA_DC, SA_PLAN["_wp"], SA_HP))
    print()
    print("-" * 78)"""
NB1_PLOT_NEW = """    print(hybrid_summary(SA_FIELD, SA_DC, SA_PLAN["_wp"], SA_HP))

    # @@MARK1@@
    # The picture the numbers above describe: where the grit goes, how deep it
    # is at each station, and where that crosses dc. The depth of cut the
    # transition happens at is read straight off the vertical axis.
    import dataclasses
    import matplotlib.pyplot as plt
    from semgrit.analysis import wheel_motion
    from semgrit_multi.plot import trajectory_figure

    _st = float((SA_PLAN.get("cost") or {}).get("step_time_s") or 0.0)
    _an = dataclasses.replace(SA_PARAMS.analysis,
                              depth_of_cut_um=float(SA_PLAN["depth_of_cut_um"]))
    _mot = wheel_motion(_an, SA_PLAN["_place"]["theta_c"],
                        SA_PARAMS.surface_speed_mm_s,
                        SA_PARAMS.outer_radius_mm, _st)
    try:
        _fig = trajectory_figure(SA_PLAN["_place"], _mot, SA_PLAN["_wp"],
                                 SA_DC, step_time_s=_st,
                                 rotation_reversed=bool(
                                     SA_PARAMS.analysis.rotation_reversed))
        plt.show()
    except ValueError as _exc:
        print("trajectory not drawn: %s" % _exc)

    print()
    print("-" * 78)"""


# --------------------------------------------------------------------------
# 2. the multi-abrasive notebook: trajectory + sweep figures, and a
#    measured-trajectory cell
# --------------------------------------------------------------------------

NB2_PLOT_OLD = """fig = preview_figure(MA_ENV, MA_DC, MA_PLAN["plan"]["_wp"], title=MA_NAME)
plt.show()
fig2 = field_slice_figure(MA_ENV, MA_DC, MA_PLAN["plan"]["_wp"])
plt.show()
print("Nothing written. When the split is what you want, run the next cell.")"""

NB2_PLOT_NEW = """# @@MARK2@@
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

print("Nothing written. When the split is what you want, run the next cell.")"""

NB2_SWEEPPARAM_OLD = """MA_FACET_SUBDIVISION = 2         #@param {type:"integer"}
print("settings captured")"""
NB2_SWEEPPARAM_NEW = """MA_FACET_SUBDIVISION = 2         #@param {type:"integer"}

#@markdown ### The depth-of-cut sweep curve
MA_SWEEP_DEPTHS = "0.10,0.20,0.30,0.40,0.60,0.80"  #@param {type:"string"}
#@markdown &nbsp;&nbsp;Depths of cut, in microns, to plot the ductile share
#@markdown against. Each one is a real sweep and costs a few seconds. Leave it
#@markdown empty to skip the curve.
print("settings captured")"""

# the measured-trajectory cell, inserted before the build cell
NB2_TRAJ_CELL = CO + '''
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
''' + Q + "))\n"

NB2_BUILD_OLD = """    MA_INFO = build_multi(MA_PARAMS, SOLIDS, MA_OUT)"""
NB2_BUILD_NEW = """    MA_INFO = build_multi(MA_PARAMS, SOLIDS, MA_OUT,
                          paths=(MA_PATHS if "MA_PATHS" in dir() else None))"""

NB2_ANCHOR_BUILD = CO + """
#@title 6 - Build it, inject the field, verify it { display-mode: "form" }"""


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))

    p1 = os.path.join(here, "_make_notebook.py")
    s1 = io.open(p1, encoding="utf-8").read()
    if MARK1 in s1:
        print("_make_notebook.py: already applied")
    else:
        s1 = edit(s1, NB1_FILES_OLD, NB1_FILES_NEW, p1)
        s1 = edit(s1, NB1_PLOT_OLD,
                  NB1_PLOT_NEW.replace("@@MARK1@@", MARK1), p1)
        io.open(p1, "w", encoding="utf-8", newline="\n").write(s1)
        print("_make_notebook.py: trajectory figure added to B2")

    p2 = os.path.join(here, "_make_notebook2.py")
    s2 = io.open(p2, encoding="utf-8").read()
    if MARK2 in s2:
        print("_make_notebook2.py: already applied")
        return 0
    s2 = edit(s2, NB2_SWEEPPARAM_OLD, NB2_SWEEPPARAM_NEW, p2)
    s2 = edit(s2, NB2_PLOT_OLD, NB2_PLOT_NEW.replace("@@MARK2@@", MARK2), p2)
    s2 = edit(s2, NB2_BUILD_OLD, NB2_BUILD_NEW, p2)
    s2 = edit(s2, NB2_ANCHOR_BUILD, NB2_TRAJ_CELL + "\n" + NB2_ANCHOR_BUILD, p2)
    # cell 5 now needs dataclasses and plan_multi at module scope
    s2 = edit(s2, "from semgrit_multi.build import MultiParams, plan_multi, summary_text",
              "import dataclasses\nfrom semgrit_multi.build import "
              "MultiParams, plan_multi, summary_text", p2)
    io.open(p2, "w", encoding="utf-8", newline="\n").write(s2)
    print("_make_notebook2.py: 4 figures + the measured-trajectory cell added")
    return 0


if __name__ == "__main__":
    sys.exit(main())
