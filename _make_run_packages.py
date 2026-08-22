"""Build the three ready-to-submit Abaqus packages into RUN_ME/.

    python _make_run_packages.py

Each package is a folder holding one ``.inp``, the subroutine it needs, the
report JSON, the post-processing script and a one-line command file. All three
share the same wheel, the same grain, the same seating and the same mesh, so any
difference between their results is the constitutive law and nothing else.

    1_single_abrasive     one grit, chip thickness from the four wedge
                          constants in the card       -> vumat_grind.for
    2_multi_abrasive      several grits, chip thickness swept per element and
                          carried as field variable 1 -> vumat_grind.for
    3_energy_criterion    the same deck as 2, switched to the local
                          W_p L_c >= PSI Kc^2/E rule  -> vumat_grind2.for

Every package is verified before it is written, by the gates that apply to it,
and the script refuses to leave a package behind that did not pass.
"""

from __future__ import annotations

import hashlib
import json
import os
import pickle
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LIBRARY = os.path.join(HERE, "WHEEL_FIXED", "1_measurements",
                       "grain_library.pkl")
OUT = os.path.join(HERE, "RUN_ME")
BUILD = os.path.join(HERE, "_runbuild")

# ---------------------------------------------------------------------------
# The configuration that makes the transition VISIBLE, and why each number.
#
# dc form 2 (Bifano, the calibrated one) gives 87.75 nm on this rock, where
# form 1 gives 5.3 nm. At 5.3 nm nothing is ductile at any depth of cut a mesh
# can resolve -- which is honest physics but shows nothing, so the packages use
# form 2 and say so.
#
# The depth of cut then has to straddle dc across the pass -- h ramps from zero
# to roughly ae along the scratch. It is per material, so it lives in
# MATERIAL_SETUP below along with the sweep it came from.
#
# Resolving a sub-micron cut needs a surface layer far finer than the 0.3 um the
# other decks use, so the depth is GRADED: fine elements at the face, then growing
# at 1.45 into the body. dt follows the smallest element and would have paid for
# that layer anyway; the grading is what stops the coarse body below costing 150
# layers of elements.
#
# THE DEPTH ELEMENT IS DERIVED FROM dc, PER MATERIAL, and that is the whole point
# of ELEMENTS_PER_DC. It used to be a flat 0.03 um for both materials, which is
# 2.9 elements through sandstone's dc and 1.76 through SiC's -- the SiC mesh had
# been sized for sandstone's dc and never revisited, so the ductile/brittle
# transition the model exists to resolve was spanned by less than two elements.
# The transition cannot be located to better than the element that straddles it.
#
# Axially 0.15 um in the groove lane, not the old flat 1.5 um. The old value came
# from "the stable increment already follows the depth element, so a fine axial
# mesh buys nothing" -- true about dt, and wrong about the mesh: 1.5 um against a
# 0.03 um depth element is an aspect ratio of 51:1, and a C3D8R that far from
# cubic integrates a bending-dominated chip badly no matter what dt is. What
# makes it affordable is grading the WIDTH too (width_band_mm): the groove is a
# few microns wide in a 20 um block, so only the lane needs the fine columns and
# the edges coarsen away at 1.35. Aspect ratio comes to ~9:1 sandstone, ~14:1 SiC.
#
# protrusion_std 0.015 is a well-dressed wheel -- grits at nearly one height.
# Without it only the tallest grit reaches the work at 0.4 um of infeed and the
# multi-abrasive case degenerates to the single one.
#
# 12 grits confined to a 0.1 mm arc window sit 8 um apart, so ten of them cross
# a 48 um block. Spread over the full 2 mm they would be 167 um apart and only
# one would ever touch it.
# ---------------------------------------------------------------------------
ARC_WINDOW_MM = 0.1
GRIT_COUNT = 12
DC_FORM = 2
ELEMENT_UM = 0.15
ELEMENT_AXIAL_UM = 0.15
# How many elements must span the critical depth of cut. Five is the smallest
# number that puts an element wholly inside the ductile zone, one wholly outside
# and one straddling; below about four the transition is an artefact of where the
# element boundary happens to fall. The depth element follows from it and from the
# material's own dc, so a new material is sized correctly without editing this.
ELEMENTS_PER_DC = 5.0
SURFACE_LAYER_UM = 0.35
DEPTH_GROWTH = 1.45
# The fine axial lane, centred on the groove. 6 um holds the widest measured grain
# (13.2 um across) at the depth it actually cuts, plus the plastic zone beside it.
WIDTH_BAND_UM = 6.0
WIDTH_GROWTH = 1.35


def element_depth_um(mat) -> float:
    """The depth element that puts ELEMENTS_PER_DC elements across dc."""
    return mat.dc_nm(DC_FORM) / 1000.0 / ELEMENTS_PER_DC
PROTRUSION_STD = 0.015
# 20 um, not 9. A 9 um block is 0.68 of ONE grain diameter (the largest measured
# grain is 13.2 um across), so grains overhang it and the groove is an edge
# chamfer rather than a confined cut. It also forced the placement sampler into
# its degenerate branch, which put all twelve grits at exactly two axial
# positions 374 nm from the free face -- two of them cutting, both in the same
# lane. 20 um holds ~1.5 grain diameters, engages many grits at a shallow
# infeed, and costs 54,080 elements against 24,960.
WP_WIDTH_MM = 0.020

# Per material: the depth of cut that straddles that material's dc, and where
# the packages go. Everything else above is shared, deliberately -- same wheel,
# same grains, same mesh -- so a difference between two materials' results is
# the material card and the depth of cut, and nothing structural.
#
# ae is NOT a fixed multiple of dc, and assuming it was is a trap worth naming.
# Sandstone wants 4.6x dc, SiC 6.8x, and the reason is that ae also has to be
# deep enough for enough ELEMENTS to be cut at all: the mesh is absolute and dc
# is not. SiC's dc is 53 nm against sandstone's 88, so scaling ae down by the
# same ratio leaves only 24 cut elements and 100% of them ductile -- a deck that
# runs, shows one regime, and looks like a modelling result. Both numbers below
# come from _sweep_depth_of_cut(), which is kept runnable for the next material:
#
#   sandstone  0.20 um -> 454 cut, 246 ductile, 208 brittle   (54% ductile)
#   SiC        0.20 um -> 454 cut, 146 ductile, 308 brittle   (32% ductile)
#
# Both materials now want the SAME ae, which is worth noticing: once the block is
# wide enough for the placement to spread properly, what sets the engagement is
# the grain-height distribution -- shared between the two materials -- and dc
# only decides where the split falls within the cut. On the old 9 um block the
# two wanted 0.40 and 0.36 because engagement was dominated by which of two
# stacked grits happened to reach.
#
# SiC is the slower of the two to run by 7x: 1.54 h on 8 cores against 0.22 h,
# 233,000 increments against 33,000. Its dilatational wave speed is 1.23e7 mm/s
# against sandstone's 1.76e6, so the stable increment is 7x smaller on the same
# mesh. Note which branch sets it -- 1.23e7 is the DUCTILE branch's
# sqrt(E(1-nu)/((1+nu)(1-2nu)rho)) with E = 450 GPa, which is stiffer than the
# JH-2 K1/G pair, and cost_model takes the faster of the two. The mesh is not
# the reason and refining it is not the fix.
MATERIAL_SETUP = {
    "sandstone":       dict(out="RUN_ME", depth_of_cut_um=0.20),
    "silicon_carbide": dict(out="RUN_ME_SIC", depth_of_cut_um=0.20),
}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(build_utc: str) -> dict:
    """What produced these decks, so it can be proved a year from now.

    Without this, RUN_ME is not reproducible from the deliverable: the report
    JSONs point into ``_runbuild_*`` trees that every rebuild deletes, and
    nothing records WHICH vumat_grind.for generated the constants that got
    published. Six lines, and it is the only thing that answers "is this the
    subroutine the paper used".
    """
    out = {"built_utc": build_utc, "python": sys.version.split()[0]}
    for f in ("vumat_grind.for", "vumat_grind2.for", "vumat_jh2.for"):
        p = os.path.join(HERE, f)
        if os.path.exists(p):
            out[f] = {"sha256": sha256(p), "bytes": os.path.getsize(p)}
    if os.path.exists(LIBRARY):
        out["grain_library"] = {"name": os.path.basename(LIBRARY),
                                "sha256": sha256(LIBRARY),
                                "bytes": os.path.getsize(LIBRARY)}
    for f in ("semgrit/hybrid.py", "semgrit/materials.py",
              "semgrit_multi/envelope.py", "semgrit_multi/swmode.py"):
        p = os.path.join(HERE, f)
        if os.path.exists(p):
            out[f] = {"sha256": sha256(p)[:16]}
    return out


def run(cmd, tag):
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=HERE)
    ok = r.returncode == 0
    last = [ln for ln in r.stdout.splitlines() if ln.strip()]
    print("   %-46s %s  %s" % (tag, "PASS" if ok else "FAIL",
                               last[-1][:60] if last else ""))
    if not ok:
        print(r.stdout[-3000:])
        print(r.stderr[-1500:])
    return ok


def command_file(folder, job, deck, subroutine, cores=8):
    """A .bat and a .sh so the command cannot be mistyped."""
    base = ("abaqus job=%s input=%s user=%s double=both"
            % (job, deck, subroutine))
    line = "%s cpus=%d interactive" % (base, cores)
    # Step 1 is a DATACHECK, and it is not politeness. It costs seconds, reads
    # the whole material card and every keyword, and it is exactly the stage the
    # only real submission this project has ever made died at. Finding a card
    # problem after an hour of solving is a wasted hour.
    rem = [
        "Submit this deck. THREE steps, in order.",
        "",
        "  0. abaqus verify -user_explicit",
        "     Does the Fortran toolchain work at all? Every gate in this",
        "     project compiled with gfortran, so this is the ONLY check that",
        "     Abaqus can build a user subroutine on this machine.",
        "  1. datacheck -- seconds. Reads the card and every keyword. The only",
        "     real submission this project has made died here, on *User",
        "     Material written four values to a line instead of eight.",
        "  2. the solve.",
        "",
        "double=both is REQUIRED, not a preference: the chip thickness is",
        "compared against a few nanometres on a 25 mm radius, a ratio of 1e-7,",
        "and single precision does not have the digits.",
        "",
        "LICENCE: cpus=%d needs int(5*%d^0.422) = %d Abaqus tokens, against 5"
        % (cores, cores, int(5 * cores ** 0.422)),
        "at cpus=1. Every wall clock in the README is the %d-core figure."
        % cores,
    ]
    # cd to the script's own folder: Abaqus resolves input= and user= against the
    # working directory, so RUN_ME\1_single_abrasive\run.bat invoked from the repo
    # root would otherwise look for the deck and the .for in the root and abort.
    # And stop on a failed datacheck -- without the guard a deck that died at
    # preprocessing still fires the solve line and burns the licence tokens.
    with open(os.path.join(folder, "run.bat"), "w", newline="") as fh:
        fh.write("@echo off\r\n")
        fh.write("cd /d \"%~dp0\"\r\n")
        for ln in rem:
            fh.write(("rem  " + ln).rstrip() + "\r\n")
        fh.write("abaqus verify -user_explicit\r\n")
        fh.write(base + " cpus=1 datacheck\r\n")
        fh.write("if errorlevel 1 exit /b 1\r\n")
        fh.write(line + "\r\n")
    with open(os.path.join(folder, "run.sh"), "w", newline="\n") as fh:
        fh.write("#!/bin/sh\nset -e\ncd \"$(dirname \"$0\")\"\n")
        for ln in rem:
            fh.write(("#  " + ln).rstrip() + "\n")
        fh.write("abaqus verify -user_explicit\n")
        fh.write(base + " cpus=1 datacheck\n")
        fh.write(line + "\n")
    return line


def main(material: str = "sandstone") -> int:
    if not os.path.exists(LIBRARY):
        print("missing grain library: %s" % LIBRARY)
        return 2
    from semgrit import materials
    from semgrit.build_deck import build_deck, hybrid_single_grit
    from semgrit_multi.build import MultiParams, build_multi
    from semgrit_multi.swmode import set_energy_mode

    mat = materials.get(material)
    setup = MATERIAL_SETUP[material]
    out = os.path.join(HERE, setup["out"])
    build = BUILD + "_" + material
    depth_of_cut_um = setup["depth_of_cut_um"]

    with open(LIBRARY, "rb") as fh:
        solids = pickle.load(fh)["solids"]
    shutil.rmtree(out, ignore_errors=True)
    shutil.rmtree(build, ignore_errors=True)
    os.makedirs(out)
    print("=" * 78)
    print("building three submittable packages from %d measured grains"
          % len(solids))
    print("material: %s" % mat.label)
    print("  dc = %.4f nm (form %d), depth of cut %.3f um = %.1f x dc"
          % (mat.dc_nm(DC_FORM), DC_FORM, depth_of_cut_um,
             depth_of_cut_um * 1000.0 / mat.dc_nm(DC_FORM)))
    print("=" * 78)
    import datetime
    build_utc = datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds")
    manifest = {"_provenance": provenance(build_utc), "_material": dict(
        key=mat.key, label=mat.label, inp_material=mat.inp_material,
        jh2=list(mat.jh2), density_kg_m3=mat.density_kg_m3,
        dc_nm=mat.dc_nm(DC_FORM), depth_of_cut_um=depth_of_cut_um,
        notes=list(mat.notes))}
    # Copy the grain library in beside the decks. Without it the packages are
    # not regenerable from the deliverable at all.
    if os.path.exists(LIBRARY):
        shutil.copy(LIBRARY, os.path.join(out, os.path.basename(LIBRARY)))

    # ------------------------------------------------------------------
    # 1. single abrasive, chip thickness from the card's four constants
    # ------------------------------------------------------------------
    print("\n1_single_abrasive")
    hp1 = mat.hybrid_params(h_source=0, dc_form=DC_FORM)
    p1 = hybrid_single_grit(
        hybrid=hp1, name="single_abrasive",
        wp_width_mm=WP_WIDTH_MM,
        wp_element_size_mm=ELEMENT_UM / 1000.0,
        wp_element_size_width_mm=ELEMENT_AXIAL_UM / 1000.0,
        wp_element_size_depth_mm=element_depth_um(mat) / 1000.0,
        wp_surface_layer_mm=SURFACE_LAYER_UM / 1000.0,
        wp_depth_growth=DEPTH_GROWTH,
        wp_width_band_mm=WIDTH_BAND_UM / 1000.0,
        wp_width_growth=WIDTH_GROWTH, protrusion_std=PROTRUSION_STD)
    p1.analysis.depth_of_cut_um = depth_of_cut_um
    # "Single abrasive" now means what it says: one grit and the workpiece, with no
    # bond rim. The rim was 2,812 of 2,928 rigid facets, all of them in general
    # contact, and it never reaches the work -- the grit stands 3.6 um proud of it.
    p1.include_bond = False
    # Moves the JH-2 card, its density and the *Material name together. Setting
    # only the hybrid params above would leave the brittle branch on sandstone.
    materials.apply(p1, material)
    d1 = os.path.join(build, "one")
    i1 = build_deck(p1, solids, d1)
    print("   %s  %.1f MB, %s C3D8R, dc = %.4f nm"
          % (os.path.basename(i1["path"]), i1["size_bytes"] / 1e6,
             format(i1["n_workpiece_elements"], ","), i1["hybrid"]["dc_nm"]))
    cf = i1["hybrid"]["chip_field"]
    print("   transition at u = %s mm of a %.3f mm scratch"
          % (cf["transition_u_mm"], p1.wp_length_mm))
    ok = (run([sys.executable, "verify_rigid_deck.py", i1["path"]],
              "verify_rigid_deck") and
          run([sys.executable, "verify_rigid_deck2.py", i1["path"]],
              "verify_rigid_deck2") and
          run([sys.executable, "verify_hybrid_deck.py", i1["path"]],
              "verify_hybrid_deck"))
    if not ok:
        print("   package 1 did not verify; not shipping it")
        return 1
    f1 = os.path.join(out, "1_single_abrasive")
    os.makedirs(f1)
    for src in (i1["path"],
                os.path.splitext(i1["path"])[0] + "_report.json",
                os.path.join(d1, "single_abrasive_placements.csv"),
                os.path.join(d1, "single_abrasive_postprocess_odb.py")):
        if os.path.exists(src):
            shutil.copy(src, f1)
    shutil.copy(os.path.join(HERE, "vumat_grind.for"), f1)
    cmd1 = command_file(f1, "single_abrasive", "single_abrasive.inp",
                        "vumat_grind.for", p1.cores)
    manifest["1_single_abrasive"] = dict(
        command=cmd1, subroutine="vumat_grind.for", n_props=56, n_depvar=20,
        dc_nm=i1["hybrid"]["dc_nm"], n_grits=i1["n_grits"],
        n_elements=i1["n_workpiece_elements"],
        chip_field=cf, h_source=0,
        what="chip thickness from the four wedge constants in the card")

    # ------------------------------------------------------------------
    # 2. several abrasives, chip thickness swept per element
    # ------------------------------------------------------------------
    print("\n2_multi_abrasive")
    hp2 = mat.hybrid_params(h_source=1, dc_form=DC_FORM)
    mp = MultiParams(name="multi_abrasive", grit_mode="count",
                     material=material,
                     grit_count=GRIT_COUNT, grit_arc_window_mm=ARC_WINDOW_MM,
                     hybrid=hp2, depth_of_cut_um=depth_of_cut_um,
                     wp_width_mm=WP_WIDTH_MM, element_um=ELEMENT_UM,
                     element_axial_um=ELEMENT_AXIAL_UM,
                     element_depth_um=element_depth_um(mat),
                     surface_layer_um=SURFACE_LAYER_UM,
                     depth_growth=DEPTH_GROWTH,
                     width_band_um=WIDTH_BAND_UM,
                     width_growth=WIDTH_GROWTH,
                     protrusion_std=PROTRUSION_STD)
    d2 = os.path.join(build, "many")
    i2 = build_multi(mp, solids, d2, log=lambda *a: None)
    sp = i2["split"]
    print("   %s  %.1f MB, %s C3D8R"
          % (os.path.basename(i2["path"]), i2["size_bytes"] / 1e6,
             format(i2["n_workpiece_elements"], ",")))
    print("   %d of %d grits cross the block; %s elements cut, %s ductile, "
          "%s brittle" % (i2["envelope"]["n_grits_engaged"], i2["n_grits"],
                          format(sp["n_cut"], ","),
                          format(sp["n_ductile_of_cut"], ","),
                          format(sp["n_brittle_of_cut"], ",")))
    ok = (run([sys.executable, "verify_rigid_deck.py", i2["plain_path"]],
              "verify_rigid_deck (un-injected geometry)") and
          run([sys.executable, "verify_rigid_deck2.py", i2["plain_path"]],
              "verify_rigid_deck2 (un-injected geometry)") and
          run([sys.executable, "verify_envelope.py", "--library", LIBRARY],
              "verify_envelope") and
          # On the INJECTED deck, which is the one that gets submitted. Until
          # this was added the two decks carrying the novel physics were read by
          # no deck-level gate at all -- and the gate, when finally pointed at
          # them, failed 3 of 33 because it assumed the single-grit h source.
          run([sys.executable, "verify_hybrid_deck.py", i2["path"]],
              "verify_hybrid_deck (injected field deck)"))
    if not ok:
        print("   package 2 did not verify; not shipping it")
        return 1
    f2 = os.path.join(out, "2_multi_abrasive")
    os.makedirs(f2)
    for src in (i2["path"],
                os.path.join(d2, "multi_abrasive_field_report.json"),
                os.path.join(d2, "multi_abrasive_placements.csv"),
                os.path.join(d2, "multi_abrasive_postprocess_odb.py"),
                os.path.join(d2, "multi_abrasive_h_elem.npy"),
                os.path.join(d2, "multi_abrasive_depth_removed.npy")):
        if os.path.exists(src):
            shutil.copy(src, f2)
    shutil.copy(os.path.join(HERE, "vumat_grind.for"), f2)
    cmd2 = command_file(f2, "multi_abrasive", "multi_abrasive_field.inp",
                        "vumat_grind.for", mp.cores)
    manifest["2_multi_abrasive"] = dict(
        command=cmd2, subroutine="vumat_grind.for", n_props=56, n_depvar=20,
        dc_nm=i2["dc_nm"], n_grits=i2["n_grits"],
        n_grits_engaged=i2["envelope"]["n_grits_engaged"],
        n_elements=i2["n_workpiece_elements"], split=sp, h_source=1,
        envelope=i2["envelope"],
        what="chip thickness swept per element, carried as field variable 1")

    # ------------------------------------------------------------------
    # 3. the same deck, switched to the local energy criterion
    # ------------------------------------------------------------------
    print("\n3_energy_criterion")
    f3 = os.path.join(out, "3_energy_criterion")
    os.makedirs(f3)
    e_deck = os.path.join(f3, "energy_criterion.inp")
    # SWMODE 2, not 1. With PSI derived from dc the local criterion needs a
    # plastic displacement of dc within one element, which at a 0.03 um element
    # and dc = 88 nm is a plastic strain of about 3 -- reachable in a chip but
    # not in the ploughing zone. SWMODE 1 alone would therefore start every
    # point ductile and might flip none of them, which shows nothing. SWMODE 2
    # keeps the geometric split visible AND lets the work criterion move it,
    # which is the comparison worth running.
    inf = set_energy_mode(
        i2["path"], e_deck, 2, 0.0,
        comment=["Same geometry, same grains, same mesh and the same swept",
                 "field as 2_multi_abrasive. Only the criterion differs, so a",
                 "difference in the result is the criterion.",
                 "SWMODE 2: the geometric split at t=0, plus points that flip",
                 "once their plastic work reaches the threshold. Watch SDV22:",
                 "it reaching 1 is what flips a point. If it never approaches",
                 "1, lower PSI on PROPS(58) and say what you lowered it to."])
    print("   %s  %.1f MB, constants 56 -> %d, *Depvar %d -> %d"
          % (os.path.basename(e_deck), inf["size_bytes"] / 1e6,
             inf["n_props"], inf["depvar_was"], inf["n_depvar"]))
    # Copy the report BEFORE gating: verify_hybrid_deck reads the report that
    # sits beside the deck, so gating first fails on a file the build was about
    # to write.
    for src in (os.path.join(d2, "multi_abrasive_field_report.json"),
                os.path.join(d2, "multi_abrasive_postprocess_odb.py")):
        if os.path.exists(src):
            dst_ = os.path.join(
                f3, os.path.basename(src).replace("multi_abrasive_field",
                                                  "energy_criterion")
                .replace("multi_abrasive", "energy_criterion"))
            shutil.copy(src, dst_)
            if dst_.endswith("_report.json"):
                # Same re-stamp as the ablation arms: the energy deck is the
                # multi deck with PROPS(56)/(58) rewritten, so it is a different
                # size and the copied report described the wrong file.
                with open(dst_) as fh:
                    _rep = json.load(fh)
                _rep["size_bytes"] = inf["size_bytes"]
                _rep["path"] = e_deck
                with open(dst_, "w") as fh:
                    json.dump(_rep, fh, indent=2, default=str)
    ok = (run([sys.executable, "verify_vumat_grind2.py"],
              "verify_vumat_grind2") and
          run([sys.executable, "verify_hybrid_deck.py", e_deck],
              "verify_hybrid_deck (energy deck, %d props / %d SDVs)"
              % (inf["n_props"], inf["n_depvar"])))
    if not ok:
        print("   package 3 did not verify; not shipping it")
        return 1
    shutil.copy(os.path.join(HERE, "vumat_grind2.for"), f3)
    cmd3 = command_file(f3, "energy_criterion", "energy_criterion.inp",
                        "vumat_grind2.for", mp.cores)
    manifest["3_energy_criterion"] = dict(
        command=cmd3, subroutine="vumat_grind2.for",
        n_props=inf["n_props"], n_depvar=inf["n_depvar"],
        swmode=2, psi="0 = derived from dc, threshold H*dc",
        dc_nm=i2["dc_nm"], h_source=1,
        what="geometric split plus the local W_p L_c >= PSI Kc^2/E rule")

    # ------------------------------------------------------------------
    # 4. the ablation arms: the same deck under each pure law
    # ------------------------------------------------------------------
    print(chr(10) + "4_ablation")
    from semgrit_multi.ablate import ARMS, write_arms
    f4 = os.path.join(out, "4_ablation")
    arms = write_arms(i2["path"], f4, cores=mp.cores)
    for name, a in sorted(arms.items()):
        shutil.copy(os.path.join(HERE, a["subroutine"]),
                    os.path.dirname(a["path"]))
        rep_src = os.path.join(d2, "multi_abrasive_field_report.json")
        if os.path.exists(rep_src):
            # Copy the multi report, but re-stamp the fields that are about THIS
            # file. An arm differs from the multi deck by the PROPS(56) token, so
            # a verbatim copy carries the wrong size and path, and
            # verify_rigid_deck2's "report matches the file on disk" check --
            # correctly -- fails on it.
            with open(rep_src) as fh:
                _rep = json.load(fh)
            _rep["size_bytes"] = a["size_bytes"]
            _rep["path"] = a["path"]
            with open(os.path.join(
                    os.path.dirname(a["path"]),
                    os.path.basename(a["path"]).replace(
                        ".inp", "_report.json")), "w") as fh:
                json.dump(_rep, fh, indent=2, default=str)
        print("   %-18s h_source %d -> %d, %.1f MB"
              % (name, a["h_source_was"], a["h_source"],
                 a["size_bytes"] / 1e6))
    ok = all(run([sys.executable, "verify_hybrid_deck.py", a["path"]],
                 "verify_hybrid_deck (%s)" % name)
             for name, a in sorted(arms.items()))
    if not ok:
        print("   the ablation arms did not verify; not shipping them")
        return 1
    manifest["4_ablation"] = {
        name: dict(command=a["command"], h_source=a["h_source"],
                   subroutine=a["subroutine"], why=a["why"],
                   folder=a.get("folder", name))
        for name, a in sorted(arms.items())}
    manifest["4_ablation"]["_what"] = (
        "The same deck as 2_multi_abrasive with PROPS(56) changed and nothing "
        "else. forced_brittle is bit-identical to vumat_jh2.for, so the hybrid "
        "arm sits bracketed between two known references instead of standing "
        "alone.")

    with open(os.path.join(out, "MANIFEST.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    _readme(out, manifest, mat, depth_of_cut_um,
            float((i1.get("cost") or {}).get("est_hours_at_requested_cores")
                  or 0.0))

    print("\n" + "=" * 78)
    for k in sorted(k for k in manifest if not k.startswith("_")):
        cmd = manifest[k].get("command")
        print("  %-20s %s" % (k, cmd if cmd else "(3 ablation arms)"))
    print("=" * 78)
    print("%s/ is ready. Every package verified before it was written."
          % os.path.basename(out))
    return 0


def _readme(out, man, mat, ae_um, est_hours=0.0):
    L = []
    a = L.append
    a("# %s - three submittable Abaqus packages" % os.path.basename(out))
    a("")
    a("**Workpiece material: %s**" % mat.label)
    a("")
    a("| | |")
    a("|---|---|")
    a("| critical depth of cut `dc` | %.4f nm (form %d) |"
      % (mat.dc_nm(DC_FORM), DC_FORM))
    a("| depth of cut `ae` | %.3f um = %.1f x dc |"
      % (ae_um, ae_um * 1000.0 / mat.dc_nm(DC_FORM)))
    a("| JH-2 `K1, G, HEL, PHEL, T` | %.6g, %.6g, %.6g, %.6g, %.6g MPa |"
      % (mat.jh2[0], mat.jh2[1], mat.jh2[2], mat.jh2[3], mat.jh2[4]))
    a("| density | %.6g kg/m3 |" % mat.density_kg_m3)
    a("| hardness `H`, toughness `Kc` | %.6g MPa, %.6g MPa*sqrt(m) |"
      % (mat.dc["hardness_mpa"], mat.dc["kic_mpa_sqrt_m"]))
    a("| estimated wall clock, 8 cores | %.2f h per package |" % est_hours)
    a("")
    for n in mat.notes:
        a("> %s" % n)
    a("")
    a("Same wheel, same grain library and same mesh throughout. Packages 2, 3")
    a("and the three arms in 4 also share the same **seating** and the same")
    a("swept field, so a difference between those is the criterion or the law")
    a("and nothing else. That is the comparison to quote.")
    a("")
    a("**Package 1 is not directly comparable to them.** It is one grit, not")
    a("twelve, so its seating, its ground radius and its peak chip thickness")
    a("all differ -- `MANIFEST.json` in this folder carries the numbers. Treat")
    a("it as a smoke test, and as the closed-form check the swept field was")
    a("validated against, rather than as an arm of the comparison.")
    a("")
    a("| folder | what | subroutine |")
    a("|---|---|---|")
    keys = sorted(k for k in man if not k.startswith("_"))
    for k in keys:
        if "what" not in man[k]:
            a("| `%s` | %s | (three arms) |" % (k, man[k].get("_what", "")))
            continue
        a("| `%s` | %s | `%s` |" % (k, man[k]["what"], man[k]["subroutine"]))
    a("")
    a("## Running them")
    a("")
    a("Each folder has `run.bat` (Windows) and `run.sh`. Or by hand:")
    a("")
    a("```")
    for k in keys:
        a("cd %s" % k)
        if "command" in man[k]:
            a(man[k]["command"])
        else:
            for sub in sorted(x for x in man[k] if not x.startswith("_")):
                a("cd %s && %s && cd .."
                  % (man[k][sub].get("folder", sub), man[k][sub]["command"]))
        a("")
    a("```")
    a("")
    a("`double=both` is **required**, not a preference: the chip thickness is")
    a("compared against a threshold of a few nanometres on a 25 mm radius, a")
    a("ratio of 1e-7, and single precision does not have the digits.")
    a("")
    a("The `.for` file must keep its name. A filename with a space or a")
    a("bracket makes Abaqus read part of it as a separate argument and abort.")
    a("")
    a("## What to plot afterwards")
    a("")
    a("| SDV | |")
    a("|---|---|")
    a("| **13** | the branch: 1 ductile, 2 brittle. This is the picture. |")
    a("| **14** | the chip thickness that point was given |")
    a("| **15** | dc |")
    a("| **19** | the strain-gradient amplification, 1 = no size effect |")
    a("| **21, 22** | plastic work and the energy ratio (package 3 only); "
      "SDV22 reaching 1 is what flips a point |")
    a("| 1, 2, 12 | damage, equivalent plastic strain, STATUS |")
    a("")
    a("`RF` and `RM` at `A_WHEEL_REF` are the grinding force. The")
    a("`*_postprocess_odb.py` in each folder reads exactly those:")
    a("")
    a("```")
    a("abaqus python <name>_postprocess_odb.py <job>.odb")
    a("```")
    a("")
    a("That writes the CSVs and the summary JSON. It also tries to draw the")
    a("figures, but **Abaqus' bundled Python usually has no matplotlib**, so in")
    a("practice it writes data and no pictures -- which is why the first six runs")
    a("of this project ended up documented by photographs of a screen. Draw them")
    a("with the host Python instead, which needs no Abaqus and no re-run:")
    a("")
    a("```")
    a("python REPOST/plots.py <folder holding the CSVs>")
    a("```")
    a("")
    a("It also writes `compare_all.png` across every job it finds, and marks on")
    a("the figures anything that should stop a number being quoted -- artificial")
    a("energy over 5% of internal, and byte-identical duplicate datasets.")
    a("")
    a("## What the three should show")
    a("")
    a("Package 1 puts the transition at one station along one scratch, from a")
    a("wedge written into the card. Package 2 computes the chip thickness for")
    a("every element from the real grit trajectories, so several grits and")
    a("their overlapping grooves are handled properly. Package 3 reaches the")
    a("transition from the material point's own plastic work instead of from")
    a("geometry at all. **If 2 and 3 disagree, that disagreement is a result**")
    a("-- it is the difference between a geometric and an energetic reading of")
    a("the same transition.")
    a("")
    a("## Read this before quoting a number")
    a("")
    a("* The **Johnson-Cook constants are placeholders**. `A` is tied to the")
    a("  JH-2 card's own quasi-static compressive strength so the two branches")
    a("  meet at the transition, but `B, n, C, m` and `D1..D5` are")
    a("  order-of-magnitude values. The decks say so in their own headers.")
    a("* `lambda_c` belongs to whichever `dc` form was used. The two published")
    a("  forms differ by `(E/H)^1.5`, about 17x on this rock, and the energy")
    a("  criterion is a third member of the family again.")
    a("* `PSI` in package 3 is **mesh-dependent** by construction. Quote the")
    a("  element size with it.")
    a("* Check `SDV13` against the split the build reported. If they differ,")
    a("  the field did not reach the material points.")
    with open(os.path.join(out, "README.md"), "w", encoding="utf-8",
              newline="\n") as fh:
        fh.write("\n".join(L) + "\n")


def _sweep_depth_of_cut(material="sandstone", values=None):
    """How the ductile/brittle split moves with the depth of cut.

        python _make_run_packages.py --sweep silicon_carbide

    This is how the ``depth_of_cut_um`` in MATERIAL_SETUP was chosen, and it is
    kept runnable so a new material can be tuned the same way instead of by
    guessing a multiple of dc. Prints a table; writes nothing.
    """
    from semgrit import materials
    from semgrit_multi.build import MultiParams, plan_multi

    mat = materials.get(material)
    dc_nm = mat.dc_nm(DC_FORM)
    if values is None:
        values = [round(f * dc_nm / 1000.0, 4) for f in (2, 3, 4.5, 6, 9)]
    with open(LIBRARY, "rb") as fh:
        solids = pickle.load(fh)["solids"]
    print("%s: dc = %.4f nm (form %d)" % (mat.label, dc_nm, DC_FORM))
    print("%9s %7s %9s %9s %9s %7s" % ("ae [um]", "ae/dc", "cut", "ductile",
                                       "brittle", "%duct"))
    for ae in values:
        mp = MultiParams(name="sweep", grit_mode="count", material=material,
                         grit_count=GRIT_COUNT,
                         grit_arc_window_mm=ARC_WINDOW_MM,
                         hybrid=mat.hybrid_params(h_source=1, dc_form=DC_FORM),
                         depth_of_cut_um=ae, wp_width_mm=WP_WIDTH_MM,
                         element_um=ELEMENT_UM,
                         element_axial_um=ELEMENT_AXIAL_UM,
                         element_depth_um=element_depth_um(mat),
                         surface_layer_um=SURFACE_LAYER_UM,
                         depth_growth=DEPTH_GROWTH,
                         width_band_um=WIDTH_BAND_UM,
                         width_growth=WIDTH_GROWTH,
                         protrusion_std=PROTRUSION_STD)
        sp = plan_multi(mp, solids, log=lambda *a: None)["split"]
        n = max(sp["n_cut"], 1)
        print("%9.4f %7.1f %9s %9s %9s %6.1f%%"
              % (ae, ae * 1000.0 / dc_nm, format(sp["n_cut"], ","),
                 format(sp["n_ductile_of_cut"], ","),
                 format(sp["n_brittle_of_cut"], ","),
                 100.0 * sp["n_ductile_of_cut"] / n))
    return 0


if __name__ == "__main__":
    argv = sys.argv[1:]
    if argv and argv[0] == "--sweep":
        sys.exit(_sweep_depth_of_cut(argv[1] if len(argv) > 1
                                     else "sandstone"))
    if argv and argv[0] == "--all":
        sys.exit(max(main(m) for m in MATERIAL_SETUP))
    sys.exit(main(argv[0] if argv else "sandstone"))
