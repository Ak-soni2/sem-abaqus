"""Where this model sits against classical grinding theory.

The deck can be geometrically perfect and still describe a process nobody would call
grinding -- a depth of cut no grain reaches, a mesh too coarse to carry a chip, a
contact zone shorter than the block. Those are the first questions a reviewer asks, and
none of them is answered by verifying the ``.inp``.

Two columns, deliberately:

**measured** rows are counted off the geometry the deck actually contains -- grain
density from the grains that were placed, active grains from the ones that reach the
work at this infeed, mesh resolution from the elements that were written. They carry no
assumptions.

**theory** rows are the textbook expressions. They assume a *traverse* grind at some
work speed ``v_w``. This model is a rotating wheel with radial infeed against a fixed
block, which is a plunge configuration, so ``v_w`` has to be supplied to say what
traverse case the numbers correspond to. With ``work_speed_mm_s = 0`` those rows are
reported as not applicable rather than quietly computed from a speed of zero.

Where a quantity appears in both columns, agreement is the claim worth making; a
disagreement is the finding.

Symbols follow Malkin & Guo, *Grinding Technology*: ``a_e`` depth of cut, ``d_e``
equivalent diameter, ``l_c`` geometric contact length, ``C`` active grain density,
``r`` chip width-to-thickness ratio, ``h_max`` maximum undeformed chip thickness,
``h_eq`` equivalent chip thickness.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

# Ratio of chip width to chip thickness. Not measurable here -- it depends on the
# grain's cutting-edge geometry -- and the literature spans about 5 to 20, so it is an
# input with a stated default rather than a constant.
DEFAULT_SHAPE_FACTOR = 10.0

# What counts as an ordinary surface-grinding regime. Only used to decide whether to
# raise a note; nothing is clamped.
SANE_WHEEL_SPEED_M_S = (10.0, 80.0)


def report(plan: dict, *, work_speed_mm_s: float = 0.0,
           shape_factor: float = DEFAULT_SHAPE_FACTOR) -> dict:
    """Everything in one dict: ``measured``, ``theory``, ``mesh`` and ``notes``."""
    wp = plan.get("workpiece")
    cost = plan.get("cost") or {}
    R = plan["outer_radius_mm"]
    d_e = 2.0 * R                                   # surface grinding, flat work
    v_s = float(cost.get("surface_speed_mm_s") or 0.0)
    ae_um = float(plan.get("depth_of_cut_um") or 0.0)
    ae = ae_um / 1000.0
    notes: list[str] = []

    # ---- measured off the deck ------------------------------------------
    swept = np.asarray(plan.get("swept_clearances_um") or [], dtype=float)
    active = int((swept <= ae_um).sum()) if (swept.size and ae_um) else 0
    band_a = float(plan.get("grit_band_arc_mm") or 0.0)
    band_w = float(plan.get("grit_band_width_mm") or 0.0)
    band_area = band_a * band_w
    C_static = (plan["n_grits"] / band_area) if band_area > 0 else 0.0
    prot = plan.get("protrusion_um") or {}
    gh = plan.get("grain_height_um") or {}
    b_mm = wp["width_mm"] if wp else 0.0

    # The arc over which a grain is within the depth of cut of the flat ground face:
    # a chord of the rim, so 2*sqrt(2 R ae) to first order. This is the model's own
    # contact length, computed from its radius and its infeed and nothing else.
    l_c_geom = 2.0 * math.sqrt(max(2.0 * R * ae - ae * ae, 0.0)) if ae > 0 else 0.0
    # Active grain density: grains that actually reach the work, per unit of the area
    # they sweep. The swept band is the block plus its travel, times the block width.
    swept_area = ((wp["length_mm"] + (cost.get("travel_mm") or 0.0)) * b_mm
                  if wp else 0.0)
    C_active = (active / swept_area) if swept_area > 0 else 0.0

    measured = {
        "wheel_diameter_mm": d_e,
        "wheel_speed_m_s": v_s / 1000.0,
        "wheel_rpm": (v_s / R) * 30.0 / math.pi if R else 0.0,
        "depth_of_cut_um": ae_um,
        "contact_length_mm": l_c_geom,
        "grains_placed": plan["n_grits"],
        "static_grain_density_per_mm2": C_static,
        "mean_grain_spacing_um": (1000.0 / math.sqrt(C_static)) if C_static > 0 else 0.0,
        "grains_under_the_block": int(plan.get("n_grits_under_block") or 0),
        "grains_in_the_swept_band": int(swept.size),
        "active_grains": active,
        "active_grain_density_per_mm2": C_active,
        "mean_protrusion_um": float(prot.get("mean") or 0.0),
        "max_protrusion_um": float(prot.get("max") or 0.0),
        "mean_grain_height_um": float(gh.get("mean") or 0.0),
        "grain_contact_time_us": (l_c_geom / v_s * 1e6) if v_s > 0 else 0.0,
        "step_time_us": float(cost.get("step_time_s") or 0.0) * 1e6,
        "wheel_travel_mm": float(cost.get("travel_mm") or 0.0),
    }

    # ---- classical, for an equivalent traverse grind ---------------------
    v_w = float(work_speed_mm_s or 0.0)
    theory: dict = {"work_speed_mm_s": v_w, "shape_factor_r": shape_factor,
                    "applicable": v_w > 0}
    if ae > 0:
        theory["contact_length_mm"] = math.sqrt(ae * d_e)
    if v_w > 0 and v_s > 0 and ae > 0:
        C_use = C_active if C_active > 0 else C_static
        theory["speed_ratio_vs_over_vw"] = v_s / v_w
        theory["equivalent_chip_thickness_um"] = ae * v_w / v_s * 1000.0
        if C_use > 0:
            # Malkin's h_max for a plunge-equivalent traverse grind.
            h_max = math.sqrt((4.0 * v_w) / (v_s * C_use * shape_factor)
                              * math.sqrt(ae / d_e))
            theory["max_chip_thickness_um"] = h_max * 1000.0
        theory["removal_rate_mm3_s_per_mm"] = ae * v_w
        if b_mm:
            theory["removal_rate_mm3_s"] = ae * v_w * b_mm
    else:
        notes.append(
            "work speed not given, so the theory column is geometry only: this deck is "
            "a plunge (fixed block, radial infeed), and chip-thickness formulas need a "
            "traverse speed. Set WORK_SPEED_MM_S to the case you want to compare with.")

    # ---- is the mesh able to carry the answer ----------------------------
    el = plan.get("element_um") or (0.0, 0.0, 0.0, 0.0)
    h_cut, h_ax, h_dep_lo = el[0], el[1], el[2]
    ref = theory.get("max_chip_thickness_um") or ae_um
    mesh = {
        "element_cutting_um": h_cut, "element_axial_um": h_ax,
        "element_depth_um": h_dep_lo,
        "elements_across_the_cut": (ae_um / h_dep_lo) if h_dep_lo > 0 else 0.0,
        "elements_across_a_chip": (ref / h_dep_lo) if h_dep_lo > 0 else 0.0,
        "reference_chip_um": ref,
        "n_elements": int(plan.get("n_workpiece_elements") or 0),
        "stable_dt_s": float(cost.get("stable_dt_s") or 0.0),
    }
    gw = plan.get("grain_width_um") or {}
    mesh["elements_across_a_grain"] = ((gw.get("mean") or 0.0) / h_cut) if h_cut else 0.0

    # ---- the findings ----------------------------------------------------
    if ae_um <= 0:
        notes.append("no depth of cut set, so nothing here describes a cut: enable the "
                     "run-ready analysis or set DEPTH_OF_CUT_UM.")
    elif active == 0:
        notes.append("NO grain reaches the work at %.3f um of infeed. The wheel would "
                     "turn for the whole step and cut nothing." % ae_um)
    elif active < 3:
        notes.append("only %d grain(s) cut at this infeed, so the result is one "
                     "grain's scratch rather than a grinding average." % active)
    if wp and l_c_geom > wp["length_mm"]:
        notes.append(
            "the geometric contact arc is %.3f mm but the block is only %.3f mm long, "
            "so this models a slice through the contact zone, not the whole of it. "
            "That is a normal micro-scale choice -- state it, do not scale forces to "
            "the full contact without accounting for it."
            % (l_c_geom, wp["length_mm"]))
    if mesh["elements_across_a_chip"] and mesh["elements_across_a_chip"] < 3:
        notes.append(
            "only %.1f elements through the reference chip thickness (%.3f um over a "
            "%.3f um element). Below about 3 the chip is not resolved and the force "
            "will read low; refine WP_ELEM_DEPTH_MM."
            % (mesh["elements_across_a_chip"], ref, h_dep_lo))
    if mesh["elements_across_a_grain"] and mesh["elements_across_a_grain"] < 4:
        notes.append(
            "only %.1f elements across a mean-width grain, so the groove shape is "
            "carried by very few elements." % mesh["elements_across_a_grain"])
    vs_ms = v_s / 1000.0
    if vs_ms and not (SANE_WHEEL_SPEED_M_S[0] <= vs_ms <= SANE_WHEEL_SPEED_M_S[1]):
        notes.append("wheel speed %.1f m/s is outside the usual %.0f-%.0f m/s for "
                     "surface grinding." % (vs_ms, *SANE_WHEEL_SPEED_M_S))
    if ae_um and prot.get("max") and ae_um > prot["max"]:
        notes.append("the depth of cut exceeds even the tallest grain's protrusion, so "
                     "the bond would reach the work.")

    # ---- the OTHER threshold, which nothing here used to compute ----------
    #
    # A blunt edge stops cutting and starts ploughing below a minimum chip
    # thickness set by the edge radius and the friction, not by the material's
    # fracture behaviour (Son/Lim's form, from the stagnation point on a
    # round edge):
    #
    #     h_min / r_e = 1 - cos(pi/4 - beta/2),     beta = arctan(mu)
    #
    # This matters because it is a SECOND transition at a similar depth. On the
    # shipped decks r_e = 0.35 um and mu = 0.2 give h_min = 0.0795 um against
    # dc = 0.08775 um -- 9.4% apart. Two mechanisms with thresholds that close
    # cannot be told apart in a force trace, and only one of them is this
    # project's subject. Nothing in the repo computed h_min before, so the
    # confound was invisible.
    edge = dict(edge_radius_um=0.0, friction=0.0, h_min_um=0.0,
                h_min_over_dc=None, ratio_note="")
    r_e = float(plan.get("edge_radius_um") or 0.0)
    if r_e <= 0:
        pr = plan.get("_params")
        r_e = float(getattr(pr, "edge_radius_um", 0.0) or 0.0)
    mu = 0.0
    an = getattr(plan.get("_params"), "analysis", None)
    if an is not None:
        mu = float(getattr(an, "friction", 0.0) or 0.0)
    if r_e > 0:
        beta = math.atan(mu) if mu > 0 else 0.0
        h_min = r_e * (1.0 - math.cos(math.pi / 4.0 - beta / 2.0))
        edge.update(edge_radius_um=r_e, friction=mu, h_min_um=h_min)
        dc_um = 0.0
        hy = plan.get("hybrid") or {}
        if hy.get("dc_nm"):
            dc_um = float(hy["dc_nm"]) / 1000.0
        elif an is not None and getattr(an, "hybrid", None) is not None:
            dc_um = an.hybrid.critical_depth_mm() * 1000.0
        if dc_um > 0:
            edge["h_min_over_dc"] = h_min / dc_um
            edge["dc_um"] = dc_um
            r = h_min / dc_um
            if 0.5 <= r <= 2.0:
                notes.append(
                    "h_min from the edge radius is %.4f um and dc is %.4f um -- "
                    "within a factor of %.2f. Below h ~ h_min a blunt edge "
                    "ploughs instead of cutting, so ANY force feature near this "
                    "depth has two possible causes and the constitutive switch "
                    "is only one of them. Separate them: rebuild once at a much "
                    "smaller edge radius (0.05 um gives h_min = %.4f um) leaving "
                    "dc where it is."
                    % (h_min, dc_um, max(r, 1.0 / r) if r else 0.0,
                       0.05 * (1.0 - math.cos(math.pi / 4.0
                                              - math.atan(mu) / 2.0))))
        if ae_um and h_min > 0 and ae_um / h_min < 3.0:
            notes.append(
                "the whole cut sits at h/r_e = %.3f or below, i.e. at or under "
                "the ploughing-to-cutting threshold. Expect a high specific "
                "energy for edge-geometry reasons alone."
                % (ae_um / r_e))

    return {"measured": measured, "theory": theory, "mesh": mesh,
            "edge": edge, "notes": notes}


def format_report(rep: dict) -> str:
    """The report as a table, for a notebook cell or a log."""
    m, t, h = rep["measured"], rep["theory"], rep["mesh"]
    L: list[str] = []
    a = L.append

    def row(label, val, unit="", fmt="%.4g"):
        a("  %-34s %14s %s" % (label, (fmt % val) if val is not None else "n/a", unit))

    a("PROCESS, as the model has it")
    row("wheel diameter", m["wheel_diameter_mm"], "mm")
    row("wheel speed", m["wheel_speed_m_s"], "m/s")
    row("wheel speed", m["wheel_rpm"], "rpm", "%.0f")
    row("depth of cut  a_e", m["depth_of_cut_um"], "um")
    row("contact arc within a_e", m["contact_length_mm"], "mm")
    row("grain contact time", m["grain_contact_time_us"], "us")
    row("step time", m["step_time_us"], "us")
    row("wheel travel in the step", m["wheel_travel_mm"], "mm")

    a("")
    a("GRAINS, counted off the deck")
    row("placed on the wheel", m["grains_placed"], "", "%d")
    row("static density  C", m["static_grain_density_per_mm2"], "/mm2", "%.0f")
    row("mean spacing", m["mean_grain_spacing_um"], "um")
    row("under the block at t=0", m["grains_under_the_block"], "", "%d")
    row("in the swept band", m["grains_in_the_swept_band"], "", "%d")
    row("ACTIVE at this infeed", m["active_grains"], "", "%d")
    row("active density", m["active_grain_density_per_mm2"], "/mm2", "%.0f")
    row("mean protrusion", m["mean_protrusion_um"], "um")
    row("tallest protrusion", m["max_protrusion_um"], "um")

    a("")
    if t.get("applicable"):
        a("CLASSICAL, for an equivalent traverse grind at v_w = %.1f mm/s"
          % t["work_speed_mm_s"])
        row("speed ratio  v_s / v_w", t.get("speed_ratio_vs_over_vw"))
        row("contact length  l_c = sqrt(a_e d_e)", t.get("contact_length_mm"), "mm")
        row("equivalent chip thickness h_eq", t.get("equivalent_chip_thickness_um"),
            "um")
        row("max chip thickness  h_max", t.get("max_chip_thickness_um"), "um")
        row("removal rate", t.get("removal_rate_mm3_s"), "mm3/s")
        a("  (h_max uses the ACTIVE grain density above and r = %.0f)"
          % t["shape_factor_r"])
    else:
        a("CLASSICAL  -- needs a work speed; this deck is a plunge, so nothing here")
        a("            is computed from a traverse that was never specified.")
        row("contact length  l_c = sqrt(a_e d_e)", t.get("contact_length_mm"), "mm")

    a("")
    a("MESH, is it able to carry the answer")
    row("element, cutting direction", h["element_cutting_um"], "um")
    row("element, axial", h["element_axial_um"], "um")
    row("element, into the depth", h["element_depth_um"], "um")
    row("elements across the cut", h["elements_across_the_cut"], "", "%.1f")
    row("elements across a chip", h["elements_across_a_chip"], "", "%.1f")
    row("elements across a grain", h["elements_across_a_grain"], "", "%.1f")
    row("stable increment", h["stable_dt_s"], "s", "%.3e")

    if rep["notes"]:
        a("")
        a("WHAT TO LOOK AT")
        for n in rep["notes"]:
            a("  * " + n)
    return "\n".join(L)
