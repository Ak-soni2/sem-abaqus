"""Post-process a semgrit grinding .odb: forces, energies, material removed.

    abaqus python THIS_FILE.py <job>.odb [<name>_report.json]

Writes <job>_forces.csv, <job>_energy.csv, <job>_summary.json next to the .odb, and
PNGs if matplotlib is importable from Abaqus' Python.

Sign convention, stated because it is easy to get backwards: the reported normal force
is the component pushing the wheel *away* from the workpiece, so it is positive while
the wheel is cutting. The tangential force is taken along the wheel's surface-velocity
direction, and is also computed independently as -RM3 / R.

Both come from rigid-body equilibrium at the reference node. A contact force F acting on
the wheel at radius R gives a moment about the axis of M_z = R*F_t; the velocity BC
reports the reaction, so RF = -F and RM3 = -R*F_t. Hence F_t = -(RF . e_t) and, equally,
F_t = -RM3/R. The two should therefore agree in sign as well as magnitude, and a
disagreement means the resolution angle is wrong. (Sign derived, not measured -- there is
no Abaqus on the machine this was written on. Confirm on the first real run by checking
that the tangential force opposes the surface velocity.)

Caveat kept in view: the reaction at the reference node also carries the rigid body's
own inertia. For a sector rim this is of order 1e-7 N against grinding forces in the
milli-newton range, so it is negligible here, but it is not exactly zero.
"""
from __future__ import print_function

import json
import math
import os
import sys

try:
    from odbAccess import openOdb
except ImportError:
    print("this script must be run with Abaqus' Python:")
    print("    abaqus python %s <job>.odb" % os.path.basename(__file__))
    raise


def find_report(odb_path, given):
    if given and os.path.exists(given):
        return given
    here = os.path.dirname(os.path.abspath(odb_path))
    base = os.path.splitext(os.path.basename(odb_path))[0]
    for cand in (base + "_report.json", base + "_cae_report.json"):
        p = os.path.join(here, cand)
        if os.path.exists(p):
            return p
    for f in sorted(os.listdir(here or ".")):
        if f.endswith("_report.json"):
            return os.path.join(here, f)
    return None


def series(hr, name):
    """One history output as two lists, or None if the odb does not have it."""
    if name not in hr.historyOutputs:
        return None
    d = hr.historyOutputs[name].data
    return [float(p[0]) for p in d], [float(p[1]) for p in d]


def pick_regions(step):
    """(reference-node region, whole-model energy region)."""
    node_r, en_r = None, None
    for key in step.historyRegions.keys():
        hr = step.historyRegions[key]
        names = set(hr.historyOutputs.keys())
        if node_r is None and "RF1" in names:
            node_r = hr
        if en_r is None and "ALLIE" in names:
            en_r = hr
    return node_r, en_r


def write_csv(path, header, columns):
    n = max([len(c) for c in columns]) if columns else 0
    fh = open(path, "w")
    try:
        fh.write(",".join(header) + "\n")
        for i in range(n):
            row = []
            for c in columns:
                row.append("%.9g" % c[i] if i < len(c) else "")
            fh.write(",".join(row) + "\n")
    finally:
        fh.close()


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    odb_path = argv[1]
    rep_path = find_report(odb_path, argv[2] if len(argv) > 2 else None)
    rep = {}
    if rep_path:
        fh = open(rep_path)
        try:
            rep = json.load(fh)
        finally:
            fh.close()
        print("report     : %s" % os.path.basename(rep_path))
    else:
        print("report     : none found - forces will not be resolved into "
              "normal/tangential")

    theta = math.radians(float(rep.get("theta_workpiece_deg", 0.0)))
    R = float(rep.get("outer_radius_mm", 0.0)) or None
    params = rep.get("params", {}) or {}
    analysis = params.get("analysis") or {}

    odb = openOdb(odb_path, readOnly=True)
    try:
        step_names = list(odb.steps.keys())
        if not step_names:
            print("the odb has no steps - the job did not get as far as solving")
            return 1
        step = odb.steps[step_names[-1]]
        print("step       : %s, %d frames" % (step_names[-1], len(step.frames)))

        node_r, en_r = pick_regions(step)
        out = {"odb": os.path.basename(odb_path), "step": step_names[-1],
               "frames": len(step.frames)}

        # ---------------- forces ----------------
        forces = None
        if node_r is None:
            print("forces     : no RF history at the wheel reference node. The deck "
                  "must request it:")
            print("             *Output, history / *Node Output, nset=A_WHEEL_REF / "
                  "RF1, RF2, RF3, RM3")
        else:
            t, rf1 = series(node_r, "RF1")
            _, rf2 = series(node_r, "RF2")
            rf3 = (series(node_r, "RF3") or (None, [0.0] * len(t)))[1]
            rm3 = (series(node_r, "RM3") or (None, [0.0] * len(t)))[1]
            ct, st = math.cos(theta), math.sin(theta)
            fn, ft, fmag, ftq = [], [], [], []
            for i in range(len(t)):
                # e_r points from the axis out through the contact; e_t is the
                # direction the surface travels away from.
                r_comp = rf1[i] * ct + rf2[i] * st
                t_comp = -rf1[i] * st + rf2[i] * ct
                fn.append(-r_comp)
                ft.append(-t_comp)
                fmag.append(math.sqrt(rf1[i] ** 2 + rf2[i] ** 2 + rf3[i] ** 2))
                # -RM3/R, not +RM3/R. For a force F on the wheel at radius R the
                # moment about the axis is M_z = R*F_t, and the velocity BC reports the
                # *reaction*, so RM3 = -R*F_t and RF = -F. That makes both estimates
                # F_t = -(RF . e_t) and F_t = -RM3/R. The script used +RM3/R, so its
                # two "independent" estimates were guaranteed to come out equal and
                # opposite, and the docstring told you that meant your angle was wrong.
                ftq.append((-rm3[i] / R) if R else 0.0)
            base = os.path.splitext(odb_path)[0]
            write_csv(base + "_forces.csv",
                      ["time_s", "RF1_N", "RF2_N", "RF3_N", "RM3_Nmm",
                       "F_normal_N", "F_tangential_N", "F_tangential_from_torque_N",
                       "F_magnitude_N"],
                      [t, rf1, rf2, rf3, rm3, fn, ft, ftq, fmag])
            peak = max(fmag) if fmag else 0.0
            # First contact: the first time the force leaves the numerical floor. A
            # job that ran to completion with this never triggering is the failure
            # mode that looks like success.
            floor = max(peak * 0.02, 1e-9)
            first = None
            for i in range(len(t)):
                if fmag[i] > floor:
                    first = t[i]
                    break
            i_pk = fmag.index(peak) if fmag else 0
            # Mean over the ENGAGED window, not over the whole step. The infeed
            # ramps linearly from zero, so the early samples are a near-zero load
            # the grit was not carrying; averaging them in halves the reported
            # force for a reason that has nothing to do with the material.
            # first_contact_s was already computed here and then never used.
            i0 = 0
            if first is not None:
                for i in range(len(t)):
                    if t[i] >= first:
                        i0 = i
                        break
            fn_e, ft_e = fn[i0:], ft[i0:]
            forces = {"peak_magnitude_N": peak,
                      "peak_time_s": t[i_pk] if t else 0.0,
                      "peak_normal_N": max(fn) if fn else 0.0,
                      "peak_tangential_N": max(ft) if ft else 0.0,
                      "mean_normal_N_engaged":
                          (sum(fn_e) / len(fn_e)) if fn_e else 0.0,
                      "mean_tangential_N_engaged":
                          (sum(ft_e) / len(ft_e)) if ft_e else 0.0,
                      "mean_normal_N_whole_step":
                          (sum(fn) / len(fn)) if fn else 0.0,
                      "mean_tangential_N_whole_step":
                          (sum(ft) / len(ft)) if ft else 0.0,
                      "engaged_samples": len(fn_e),
                      "total_samples": len(fn),
                      "first_contact_s": first,
                      "csv": os.path.basename(base + "_forces.csv")}
            # Backwards-compatible aliases: the engaged window is the number to
            # quote, so these now point at it rather than at the whole step.
            forces["mean_normal_N"] = forces["mean_normal_N_engaged"]
            forces["mean_tangential_N"] = forces["mean_tangential_N_engaged"]
            if forces["mean_normal_N"]:
                forces["force_ratio_Ft_over_Fn"] = (forces["mean_tangential_N"]
                                                    / forces["mean_normal_N"])
            out["forces"] = forces
            print("forces     : peak |F| %.6g N at t = %.6g s" % (peak, t[i_pk]))
            if first is None:
                print("             *** THE WHEEL NEVER TOUCHED THE WORK ***")
                print("             The job ran to completion with no contact force. "
                      "Increase the depth")
                print("             of cut, or reduce the standoff, and resubmit.")
            else:
                print("             first contact at t = %.6g s (%.1f%% into the step)"
                      % (first, 100.0 * first / t[-1] if t and t[-1] else 0.0))
                print("             mean over the ENGAGED window: Fn %.6g N, Ft %.6g N"
                      % (forces["mean_normal_N"], forces["mean_tangential_N"]))

            # ---------------- force against instantaneous chip thickness ----
            # The linear infeed ramp is usually called a loading artefact, and it
            # is -- but it also means one pass sweeps the chip thickness from zero
            # to its final value, so a single run already contains the size-effect
            # curve u(h). h at time tau is the wedge the card itself carries:
            #     u   = -v_s * (t_end - tau)      station of the grit, mm
            #     h(u) = H0 + HG*u - u^2/(2*RTIP)
            # with H0, HG, RTIP the four constants written into *User Material.
            # Without this the curve is discarded by two sum()/len() calls.
            hyb = rep.get("hybrid") or {}
            cf = hyb.get("chip_field") or {}
            if cf.get("rtip_mm") and t:
                h0 = float(cf["h0_mm"])
                hg = float(cf["hg"])
                rtip = float(cf["rtip_mm"])
                dc_mm = float(hyb.get("dc_mm") or 0.0)
                vs_mm = float(analysis.get("surface_speed_mm_s")
                              or params.get("surface_speed_mm_s") or 0.0)
                t_end = t[-1]
                us, hs = [], []
                for tau in t:
                    u = -vs_mm * (t_end - tau)
                    h = h0 + hg * u - u * u / (2.0 * rtip)
                    us.append(u)
                    hs.append(h if h > 0.0 else 0.0)
                write_csv(base + "_force_vs_h.csv",
                          ["time_s", "station_u_mm", "h_mm", "h_over_dc",
                           "F_normal_N", "F_tangential_N"],
                          [t, us, hs,
                           [(h / dc_mm if dc_mm else 0.0) for h in hs], fn, ft])
                out["force_vs_h"] = {
                    "csv": os.path.basename(base + "_force_vs_h.csv"),
                    "dc_mm": dc_mm,
                    "h_range_mm": [min(hs), max(hs)],
                    "h_over_dc_range": ([min(hs) / dc_mm, max(hs) / dc_mm]
                                        if dc_mm else None),
                    "note": ("h is the PRESCRIBED wedge from the card, not a "
                             "measurement; it is exact for the rigid wheel."),
                }
                print("size effect: h swept %.4f -> %.4f nm (%.2f -> %.2f dc) in "
                      "one pass" % (min(hs) * 1e6, max(hs) * 1e6,
                                    (min(hs) / dc_mm) if dc_mm else 0.0,
                                    (max(hs) / dc_mm) if dc_mm else 0.0))
                print("             wrote %s"
                      % os.path.basename(base + "_force_vs_h.csv"))

        # ---------------- energies ----------------
        energy = None
        if en_r is None:
            print("energy     : no ALLIE history; the deck needs "
                  "*Output, history, variable=PRESELECT")
        else:
            names = ["ALLIE", "ALLKE", "ALLAE", "ALLSE", "ALLPD", "ALLDMD",
                     "ALLWK", "ALLVD", "ETOTAL"]
            cols, hdr, got = [], ["time_s"], {}
            tt = None
            for n in names:
                s = series(en_r, n)
                if s is None:
                    continue
                if tt is None:
                    tt = s[0]
                    cols.append(tt)
                hdr.append(n + "_mJ")
                cols.append(s[1])
                got[n] = s[1]
            if tt is not None:
                base = os.path.splitext(odb_path)[0]
                write_csv(base + "_energy.csv", hdr, cols)
                ie = got.get("ALLIE", [0.0])[-1]
                ae = got.get("ALLAE", [0.0])[-1]
                ke = got.get("ALLKE", [0.0])[-1]
                energy = {"ALLIE_final": ie, "ALLAE_final": ae, "ALLKE_final": ke,
                          "artificial_fraction": (ae / ie) if ie else None,
                          "kinetic_fraction": (ke / ie) if ie else None,
                          "csv": os.path.basename(base + "_energy.csv")}
                out["energy"] = energy
                print("energy     : ALLIE %.6g, ALLAE %.6g, ALLKE %.6g" % (ie, ae, ke))
                if ie:
                    af, kf = ae / ie, ke / ie
                    print("             artificial/internal %.1f%%  %s"
                          % (100 * af,
                             "ok" if af < 0.05 else
                             "HIGH - hourglassing is carrying the load, refine the "
                             "mesh or change the hourglass control"))
                    print("             kinetic/internal    %.1f%%  %s"
                          % (100 * kf,
                             "ok" if kf < 0.10 else
                             "HIGH - mass scaling is adding inertia, lower "
                             "MASS_SCALING"))

        # ---------------- the branch map, SDV13 ----------------
        # This is the figure the whole model exists to produce, and until now
        # nothing read it: the only "SDV" in this script was an error message.
        # SDV13 = 1 ductile, 2 brittle; 14 = the h that point was given;
        # 15 = dc; 19 = the SGE amplification; 21/22 = plastic work and the
        # energy ratio (package 3 only, which carries 22 SDVs not 20).
        sdv = None
        try:
            last = step.frames[-1]
            names = [k for k in last.fieldOutputs.keys()]
            want = [("SDV13", "branch"), ("SDV14", "h_mm"), ("SDV15", "dc_mm"),
                    ("SDV19", "sge_factor"), ("SDV21", "plastic_work"),
                    ("SDV22", "energy_ratio")]
            have = [(k, lab) for k, lab in want if k in names]
            if not have:
                print("branch map : no SDV fields in the .odb. The deck needs")
                print("             *Element Output ... SDV  (it requests it, so a")
                print("             missing SDV here means the VUMAT never ran)")
            else:
                alive = {}
                if "STATUS" in names:
                    for v in last.fieldOutputs["STATUS"].values:
                        alive[v.elementLabel] = float(v.data) >= 0.5
                cols, labels = [], []
                labs_by_el = {}
                for k, lab in have:
                    d = {}
                    for v in last.fieldOutputs[k].values:
                        d[v.elementLabel] = float(v.data)
                    labs_by_el[lab] = d
                    labels.append(lab)
                els = sorted(labs_by_el[labels[0]].keys())
                cols = [els] + [[labs_by_el[l].get(e, 0.0) for e in els]
                                for l in labels]
                cols.append([(1 if alive.get(e, True) else 0) for e in els])
                write_csv(base + "_sdv.csv",
                          ["element"] + labels + ["alive"], cols)
                br = labs_by_el.get("branch") or {}
                n_d = n_b = n_other = 0
                for e in els:
                    if not alive.get(e, True):
                        continue
                    b = int(round(br.get(e, 0.0)))
                    if b == 1:
                        n_d += 1
                    elif b == 2:
                        n_b += 1
                    else:
                        n_other += 1
                sdv = {"csv": os.path.basename(base + "_sdv.csv"),
                       "fields": labels,
                       "n_ductile_alive": n_d, "n_brittle_alive": n_b,
                       "n_unset_alive": n_other,
                       "ductile_fraction_alive":
                           (float(n_d) / (n_d + n_b)) if (n_d + n_b) else None}
                # The claim the README tells you to check by hand: SDV13 against
                # the split the build predicted. Do it here instead.
                pred = rep.get("split") or {}
                if pred.get("n_ductile_law") is not None:
                    sdv["predicted_n_ductile_law"] = pred["n_ductile_law"]
                    sdv["predicted_n_brittle_law"] = pred["n_brittle_law"]
                    got_b = n_b + sum(1 for e in els if not alive.get(e, True)
                                      and int(round(br.get(e, 0.0))) == 2)
                    sdv["brittle_including_deleted"] = got_b
                out["sdv"] = sdv
                print("branch map : %s ductile, %s brittle, %s unset (of the "
                      "elements still alive)" % (n_d, n_b, n_other))
                if pred.get("n_ductile_law") is not None:
                    print("             build predicted %s ductile / %s brittle "
                          "before any deletion"
                          % (pred["n_ductile_law"], pred["n_brittle_law"]))
                if n_d and not n_b:
                    print("             *** EVERY LIVE POINT IS DUCTILE ***  If this")
                    print("             is a field-carrying deck (PROPS(56)=1), the")
                    print("             most likely cause is that field variable 1")
                    print("             never reached the VUMAT, which makes hloc 0")
                    print("             and 0 < dc everywhere. Check SDV14: it should")
                    print("             be the injected h, not zero.")
                if n_other:
                    print("             %s live elements have SDV13 unset -- the "
                          "VUMAT never ran on them" % n_other)
                print("             wrote %s" % os.path.basename(base + "_sdv.csv"))
        except Exception as exc:                      # noqa: BLE001
            print("branch map : could not be computed: %s" % exc)

        # ---------------- material removed ----------------
        removal = None
        try:
            last = step.frames[-1]
            if "STATUS" in last.fieldOutputs:
                st_f = last.fieldOutputs["STATUS"]
                total = len(st_f.values)
                dead = 0
                for v in st_f.values:
                    if float(v.data) < 0.5:
                        dead += 1
                p = params
                wpv = None
                if p.get("wp_length_mm"):
                    wpv = (float(p["wp_length_mm"]) * float(p["wp_width_mm"])
                           * float(p["wp_depth_mm"]))
                n_el = int(rep.get("n_workpiece_elements") or total or 0)
                # The SURFACE-LAYER element volume, not the mean element volume.
                #
                # The depth mesh is GRADED -- 0.03 um layers near the surface
                # growing into the body -- so the mean element is 7.7x larger
                # than any element that can be deleted. Using the mean made the
                # removed volume 7.7x too big and therefore the specific energy
                # 7.7x too SMALL, which turns a size-effect number (hundreds of
                # J/mm3) into an ordinary macro-grinding one (tens) and erases
                # the effect the model exists to show.
                #
                # Valid while the cut stays inside the uniform surface layer,
                # which is asserted below rather than assumed.
                cost = rep.get("cost") or {}
                surf_v = None
                try:
                    lx = float(cost["element_size_cutting_mm"])
                    lz = float(cost["element_size_axial_mm"])
                    ly = float(cost["depth_layer_min_mm"])
                    surf_v = lx * ly * lz
                except (KeyError, TypeError, ValueError):
                    surf_v = None
                mean_v = (wpv / n_el) if (wpv and n_el) else None
                elem_v, basis = (surf_v, "surface-layer element")                      if surf_v else (mean_v, "mean element (GRADED MESH: this "
                                            "over-estimates the volume)")
                # Does the cut actually stay in the uniform layer?
                env = rep.get("envelope") or {}
                deep_um = float(env.get("max_depth_removed_um") or 0.0)
                layer_um = float(env.get("surface_layer_um")
                                 or (1000.0 * float(cost.get(
                                     "element_size_depth_mm") or 0.0)) or 0.0)
                layer_ok = (deep_um <= layer_um) if (deep_um and layer_um) else None
                if layer_ok is False:
                    print("removal    : WARNING the deepest swept cut (%.4f um) "
                          "leaves the uniform" % deep_um)
                    print("             surface layer (%.4f um), so the element "
                          "volume is a lower" % layer_um)
                    print("             bound and the specific energy an upper one.")
                removal = {"elements_total": total, "elements_deleted": dead,
                           "deleted_fraction": (float(dead) / total) if total else 0.0,
                           "element_volume_mm3": elem_v,
                           "element_volume_basis": basis,
                           "mean_element_volume_mm3": mean_v,
                           "removed_volume_mm3_approx":
                               (dead * elem_v) if elem_v else None,
                           "cut_inside_uniform_layer": layer_ok,
                           "volume_is_approximate": True}
                mean_v = elem_v
                out["removal"] = removal
                print("removal    : %d of %d elements deleted (%.3f%%)"
                      % (dead, total, 100.0 * dead / total if total else 0.0))
                if mean_v:
                    print("             ~%.6g mm3 removed (%s)"
                          % (dead * mean_v, removal["element_volume_basis"]))
                if dead == 0:
                    print("             nothing was deleted: either nothing cut, or "
                          "the VUMAT never set")
                    print("             the deletion flag SDV%s to 0."
                          % analysis.get("n_depvar", 12))
            else:
                print("removal    : no STATUS field, so element deletion cannot be "
                      "counted")
        except Exception as exc:                      # noqa: BLE001 - report, do not die
            print("removal    : could not be computed: %s" % exc)

        # ---------------- subsurface damage, from fields already present ----
        # Nothing in this project measured subsurface damage, crack depth or
        # surface integrity -- which is the entire reason dc matters
        # industrially. Both proxies are READ from fields the deck already
        # requests; neither changes the model.
        #
        #   depth of the deepest DELETED element in each (i,j) column  -> the
        #     groove that was actually cut
        #   depth of the deepest element with damage above a threshold -> how
        #     far the damage field reached below that groove
        #
        # CAVEAT that must travel with the number: under the geometric switch h
        # depends only on the tangential station, so the whole column below a
        # ductile station runs Johnson-Cook forever. "No subsurface damage under
        # the ductile zone" is therefore imposed by the switch, not predicted by
        # it. Say so wherever this is plotted.
        try:
            divs = (rep.get("cost") or {}).get("element_divisions")
            last = step.frames[-1]
            if divs and len(divs) == 3 and "STATUS" in last.fieldOutputs:
                nl, nw, nd = int(divs[0]), int(divs[1]), int(divs[2])
                st = {}
                for v in last.fieldOutputs["STATUS"].values:
                    st[v.elementLabel] = float(v.data)
                dmg = {}
                if "SDV1" in last.fieldOutputs:
                    for v in last.fieldOutputs["SDV1"].values:
                        dmg[v.elementLabel] = float(v.data)
                # Depth of each layer, from the graded mesh if the report has it.
                cost = rep.get("cost") or {}
                dz = float(cost.get("depth_layer_min_mm") or 0.0)
                d_thresh = 0.2
                cut_k, dam_k = [], []
                for i in range(nl):
                    for j in range(nw):
                        deep_cut = -1
                        deep_dam = -1
                        for k in range(nd):
                            e = ((i * nw) + j) * nd + k + 1
                            if st.get(e, 1.0) < 0.5:
                                deep_cut = k
                            if dmg.get(e, 0.0) > d_thresh:
                                deep_dam = k
                        if deep_cut >= 0:
                            cut_k.append(deep_cut + 1)
                        if deep_dam >= 0:
                            dam_k.append(deep_dam + 1)
                sub = {"damage_threshold": d_thresh,
                       "columns_with_material_removed": len(cut_k),
                       "columns_with_damage": len(dam_k),
                       "max_removed_layers": max(cut_k) if cut_k else 0,
                       "max_damaged_layers": max(dam_k) if dam_k else 0,
                       "surface_layer_mm": dz,
                       "caveat": ("under the geometric switch the whole column "
                                  "below a ductile station is Johnson-Cook by "
                                  "construction, so an absence of subsurface "
                                  "damage there is imposed, not predicted")}
                if dz > 0:
                    sub["max_removed_depth_um_approx"] = (
                        (max(cut_k) if cut_k else 0) * dz * 1000.0)
                    sub["max_damaged_depth_um_approx"] = (
                        (max(dam_k) if dam_k else 0) * dz * 1000.0)
                out["subsurface"] = sub
                print("subsurface : %d columns cut, deepest %d layers; %d "
                      "columns damaged (D > %.2f), deepest %d layers"
                      % (len(cut_k), sub["max_removed_layers"], len(dam_k),
                         d_thresh, sub["max_damaged_layers"]))
                if dz > 0:
                    print("             ~%.4f um removed, damage reaching "
                          "~%.4f um (uniform-layer estimate)"
                          % (sub.get("max_removed_depth_um_approx", 0.0),
                             sub.get("max_damaged_depth_um_approx", 0.0)))
                if not dmg:
                    print("             no SDV1 field, so only the removed "
                          "depth could be measured")
        except Exception as exc:                      # noqa: BLE001
            print("subsurface : could not be computed: %s" % exc)

        # specific energy, if both halves are available
        if forces and removal and removal.get("removed_volume_mm3_approx"):
            vs = float(analysis.get("surface_speed_mm_s")
                       or params.get("surface_speed_mm_s") or 0.0)
            tt_end = step.frames[-1].frameValue
            vol = removal["removed_volume_mm3_approx"]
            if vs and tt_end and vol > 0:
                # u = Ft * vs / Q, with Q the volumetric removal rate
                u = forces["mean_tangential_N"] * vs / (vol / tt_end)
                out["specific_energy_J_mm3"] = u / 1000.0
                out["specific_energy_basis"] = removal["element_volume_basis"]
                print("specific   : %.4g J/mm3  (engaged-window Ft, %s)"
                      % (u / 1000.0, removal["element_volume_basis"]))
                mv = removal.get("mean_element_volume_mm3")
                ev = removal.get("element_volume_mm3")
                if mv and ev and abs(mv - ev) > 1e-30:
                    print("             the old mean-element denominator would have "
                          "reported %.4g J/mm3" % (u / 1000.0 * ev / mv))

        base = os.path.splitext(odb_path)[0]
        fh = open(base + "_summary.json", "w")
        try:
            json.dump(out, fh, indent=2)
        finally:
            fh.close()
        print("wrote      : %s" % os.path.basename(base + "_summary.json"))

        # ---------------- plots, if we can ----------------
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            # Units are scaled at the axis (us, mN, nm). Leaving SI here puts a
            # "1e-6" offset in the corner and nanometre data in the fifth decimal.
            if forces:
                fig, ax = plt.subplots(figsize=(9, 4.4))
                ax.plot([x * 1e6 for x in t], [y * 1e3 for y in fn],
                        color="#0072B2", lw=1.2, label="normal")
                ax.plot([x * 1e6 for x in t], [y * 1e3 for y in ft],
                        color="#D55E00", lw=1.2, label="tangential")
                _fc = forces.get("first_contact_s")
                if _fc:
                    ax.axvline(_fc * 1e6, color="#666666", ls=":",
                               lw=1.2, label="first contact")
                ax.set_xlabel("time (us)")
                ax.set_ylabel("force on the wheel (mN)")
                ax.set_title("%s -- grinding force" % os.path.basename(base))
                ax.legend()
                ax.grid(alpha=0.3)
                fig.tight_layout()
                fig.savefig(base + "_forces.png", dpi=200)
                print("wrote      : %s" % os.path.basename(base + "_forces.png"))
            if out.get("force_vs_h") and t:
                fig, ax = plt.subplots(figsize=(9, 4.4))
                # h in nm: dc is tens of nanometres, so millimetres put every
                # interesting point in the fifth decimal -- and the dc label was
                # already quoting nm, so the axis and the legend disagreed.
                hn = [x * 1e6 for x in hs]
                ax.plot(hn, [y * 1e3 for y in fn], ".", ms=2, color="#0072B2",
                        label="normal")
                ax.plot(hn, [y * 1e3 for y in ft], ".", ms=2, color="#D55E00",
                        label="tangential")
                if dc_mm:
                    d_nm = dc_mm * 1e6
                    ax.axvline(d_nm, color="k", ls="--", lw=1.2,
                               label="dc = %.1f nm" % d_nm)
                    ax.axvspan(0, d_nm, color="#0072B2", alpha=0.07)
                    ax.text(0.02, 0.95, "ductile  h < dc", transform=ax.transAxes,
                            fontsize=9, color="#0072B2", va="top")
                    ax.text(0.98, 0.95, "brittle  h > dc", transform=ax.transAxes,
                            fontsize=9, color="#D55E00", va="top", ha="right")
                ax.set_xlabel("undeformed chip thickness h (nm)")
                ax.set_ylabel("force on the wheel (mN)")
                ax.set_title("%s -- force against chip thickness (the size "
                             "effect)" % os.path.basename(base))
                ax.legend()
                ax.grid(alpha=0.3)
                fig.tight_layout()
                fig.savefig(base + "_force_vs_h.png", dpi=200)
                print("wrote      : %s"
                      % os.path.basename(base + "_force_vs_h.png"))
            if sdv and sdv.get("n_ductile_alive") is not None:
                divs = (rep.get("cost") or {}).get("element_divisions")
                if divs and len(divs) == 3:
                    nl, nw, nd = int(divs[0]), int(divs[1]), int(divs[2])
                    grid = [[0.0] * nl for _ in range(nd)]
                    br = labs_by_el.get("branch") or {}
                    # Elements were written i (along) outer, j (axial), k (depth)
                    # inner, so element (i,j,k) is ((i*nw)+j)*nd + k + 1. Collapse
                    # the axial index by taking the mid lane.
                    jmid = nw // 2
                    for i in range(nl):
                        for k in range(nd):
                            e = ((i * nw) + jmid) * nd + k + 1
                            v = br.get(e, 0.0)
                            if not alive.get(e, True):
                                v = -1.0
                            grid[k][i] = v
                    # FOUR DISCRETE CATEGORIES, so a discrete colormap. coolwarm
                    # is continuous and diverging: it rendered deleted and unset
                    # as near-identical pale blues and ductile and brittle as
                    # near-identical pale reds -- the one distinction the figure
                    # exists to show. Blue/vermillion also survives colour
                    # blindness, which the green/red pair used elsewhere does not.
                    from matplotlib.colors import BoundaryNorm, ListedColormap
                    cmap = ListedColormap(["#BBBBBB", "#FFFFFF",
                                           "#0072B2", "#D55E00"])
                    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
                    # Axes in microns, not element indices: the reader cannot
                    # judge "how deep did the brittle zone reach" from an index.
                    cost = rep.get("cost") or {}
                    ex = (cost.get("element_size_cutting_mm") or 0) * 1000.0 * nl
                    ey = (cost.get("depth_layer_min_mm") or 0) * 1000.0 * nd
                    extent = ([0, ex, ey, 0] if ex and ey else None)
                    fig, ax = plt.subplots(figsize=(11, 3.6))
                    im = ax.imshow(grid, aspect="auto", origin="upper",
                                   cmap=cmap, norm=norm, extent=extent,
                                   interpolation="nearest")
                    ax.set_xlabel("station along the scratch (%s)"
                                  % ("um" if extent else "element"))
                    ax.set_ylabel("depth (%s)" % ("um" if extent else "element"))
                    ax.set_title("%s -- SDV13 ductile/brittle branch, mid axial "
                                 "lane" % os.path.basename(base))
                    cb = fig.colorbar(im, ax=ax, ticks=[-1, 0, 1, 2], shrink=0.9)
                    cb.ax.set_yticklabels(["deleted", "unset", "ductile",
                                           "brittle"])
                    # State the sampling bias ON the figure. This is the last
                    # frame, by which time the elements that went brittle and
                    # failed have been deleted -- so the surviving brittle count
                    # is a lower bound, and a reader cannot know that otherwise.
                    ax.text(0.005, -0.30, "last frame only: elements that went "
                            "brittle and failed now read 'deleted', so the "
                            "brittle fraction here is a LOWER bound "
                            "(REPOST/hotspot.py walks every frame)",
                            transform=ax.transAxes, fontsize=8, color="#B00020")
                    fig.tight_layout()
                    fig.savefig(base + "_branch_map.png", dpi=200)
                    print("wrote      : %s"
                          % os.path.basename(base + "_branch_map.png"))
            if energy and tt is not None:
                fig, ax = plt.subplots(figsize=(9, 4.6))
                _c = {"ALLIE": "#0072B2", "ALLKE": "#009E73",
                      "ALLAE": "#D55E00", "ALLPD": "#CC79A7",
                      "ALLDMD": "#56B4E9"}
                for n in ("ALLIE", "ALLKE", "ALLAE", "ALLPD", "ALLDMD"):
                    if n in got and any(abs(v) > 0 for v in got[n]):
                        ax.plot([x * 1e6 for x in tt], [abs(v) for v in got[n]],
                                color=_c.get(n), lw=1.2, label=n)
                # LOG, because these span five decades: ALLKE/ALLIE has run as
                # high as 56,000 on this project's own results, and on a linear
                # axis that is one flat line with everything else at zero -- a
                # figure that conveys nothing at all.
                ax.set_yscale("log")
                ax.set_xlabel("time (us)")
                ax.set_ylabel("energy (mJ), log scale")
                ax.set_title("%s -- energy balance" % os.path.basename(base))
                # The two quality bars, drawn. Both were computed, printed to the
                # console and then left off the figure, on runs that fail them.
                _e = out.get("energy") or {}
                _msg = []
                if _e.get("artificial_fraction") is not None:
                    _msg.append("ALLAE/ALLIE = %.0f%% (bar 5%%)"
                                % (100 * _e["artificial_fraction"]))
                if _e.get("kinetic_fraction") is not None:
                    _msg.append("ALLKE/ALLIE = %.3gx (bar 0.1)"
                                % _e["kinetic_fraction"])
                if _msg:
                    _bad = (_e.get("artificial_fraction") or 0) > 0.05
                    ax.text(0.01, 0.97, "\n".join(_msg), transform=ax.transAxes,
                            va="top", fontsize=9,
                            color="#B00020" if _bad else "#666666")
                ax.legend(ncol=2, fontsize=9)
                ax.grid(alpha=0.3)
                fig.tight_layout()
                fig.savefig(base + "_energy.png", dpi=200)
                print("wrote      : %s" % os.path.basename(base + "_energy.png"))
        except ImportError:
            # This is the NORMAL path, not an edge case: Abaqus' bundled Python
            # generally has no matplotlib, which is why no run of this project has
            # ever produced a PNG from here. Say what to run instead, or the CSVs
            # sit unplotted and the results get documented by phone photographs.
            print("plots      : matplotlib not available in this Abaqus Python.")
            print("             CSVs are written; draw the figures with the host")
            print("             Python, which does not need Abaqus at all:")
            print("                 python REPOST/plots.py <dir with the CSVs>")
    finally:
        odb.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
