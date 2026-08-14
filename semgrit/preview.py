"""See the model before you build it.

A 200 MB deck takes minutes to write and longer to open, so getting a dimension wrong
is expensive. This draws the whole assembly from the *same* placement code the writer
uses -- so what you see is what the deck will contain -- without writing anything.

Five panels, because the model spans four orders of magnitude and no single view can
show a 50 mm wheel and a 3 um grit at once:

  1  the wheel sector and the workpiece, to scale, in the XY plane
  2  the dressed band unwrapped, every grit footprint and the block outline on it
  3  the contact, in microns: the rim line, grit profiles, the block, the depth of cut
  4  how many grits engage as the depth of cut increases
  5  the numbers -- dimensions, counts, dt, run time, and whatever the build warned about
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np


def _wedge(ax, r0, r1, t0, t1, **kw):
    th = np.linspace(t0, t1, 200)
    x = np.r_[r1 * np.cos(th), r0 * np.cos(th[::-1])]
    y = np.r_[r1 * np.sin(th), r0 * np.sin(th[::-1])]
    ax.fill(x, y, **kw)


def preview(plan: dict, figsize=(15, 9)):
    """Draw the plan produced by ``build_deck(..., dry_run=True)``."""
    import matplotlib
    import matplotlib.pyplot as plt

    R = plan["outer_radius_mm"]
    r0 = R - plan["rim_depth_mm"]
    sec = math.radians(plan["sector_deg"])
    full = plan["full_wheel"]
    thc = math.radians(plan["theta_workpiece_deg"])
    frames = plan["frames"]          # grit vertices in the (radial, tangential, axial) frame
    wp = plan["workpiece"]
    r_ground = plan["ground_radius_mm"]

    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(2, 3)

    # ---------------- 1. the wheel, to scale ----------------
    ax = fig.add_subplot(gs[0, 0])
    # Always draw the complete wheel as a dashed outline, then highlight the sector on
    # it. A 4.6 deg sector with a 12 um rim is a hairline on its own and tells you
    # nothing; against the whole wheel it reads immediately as a small piece of one.
    _c = np.linspace(0, 2 * math.pi, 400)
    ax.plot(R * np.cos(_c), R * np.sin(_c), ls="--", lw=0.8, color="#9aa0a6")
    ax.plot(r0 * np.cos(_c), r0 * np.sin(_c), ls="--", lw=0.6, color="#c2c6cb")
    t0, t1 = (0.0, 2 * math.pi) if full else (0.0, sec)
    _wedge(ax, r0, R, t0, t1, facecolor="#c9ccd1", edgecolor="#4a4e54", lw=1.2)
    if not full:
        for t in (t0, t1):
            ax.plot([0, R * math.cos(t)], [0, R * math.sin(t)], lw=0.6, ls=":",
                    color="#9aa0a6")
    if wp is not None:
        # the block, drawn in its own tangent frame then rotated into place
        hb, d = wp["length_mm"] / 2.0, wp["depth_mm"]
        c, s = math.cos(thc), math.sin(thc)
        pts = np.array([[r_ground, -hb], [r_ground + d, -hb],
                        [r_ground + d, hb], [r_ground, hb]])
        ax.fill(pts[:, 0] * c - pts[:, 1] * s, pts[:, 0] * s + pts[:, 1] * c,
                facecolor="#2f6fb5", edgecolor="#12365e", lw=1.0, zorder=3)
    ax.plot([0], [0], "+", color="#4a4e54", ms=9)
    if wp is not None:
        # The block is microns across on a 50 mm wheel, so it is a single pixel here.
        # Point at it rather than pretending it is visible.
        ax.annotate("workpiece\n(%.0f x %.0f um)"
                    % (wp["length_mm"] * 1000, wp["width_mm"] * 1000),
                    xy=(r_ground * math.cos(thc), r_ground * math.sin(thc)),
                    xytext=(0.55, 0.10), textcoords="axes fraction", fontsize=7,
                    color="#12365e",
                    arrowprops=dict(arrowstyle="->", color="#2f6fb5", lw=1.2))
    ax.set_aspect("equal")
    ax.set_title("wheel %.0f mm dia, %s\n(dashed = the whole wheel it is cut from)"
                 % (2 * R, "full" if full else "%.1f deg sector" % plan["sector_deg"]),
                 fontsize=9)
    ax.set_xlabel("x (mm)", fontsize=8)
    ax.set_ylabel("y (mm)", fontsize=8)
    ax.tick_params(labelsize=7)

    # ---------------- 2. the dressed band, unwrapped ----------------
    ax = fig.add_subplot(gs[0, 1:])
    if frames:
        b = np.array([f[:, 1].mean() for f in frames])
        z = np.array([f[:, 2].mean() for f in frames])
        rad = np.array([0.5 * (f[:, 1].max() - f[:, 1].min()) for f in frames])
        # True protrusion is the distance from the axis, not the projection onto the
        # tangent frame: at 750 um along the arc the curvature alone drops that
        # projection 11 um, which would read as a deeply buried grain.
        prot = np.array([float(np.hypot(f[:, 0], f[:, 1]).max()) for f in frames]) - R
        sc = ax.scatter(b * 1000, z * 1000, s=(rad * 2000) ** 2 * 0.6,
                        c=prot * 1000, cmap="viridis", alpha=0.75,
                        edgecolors="none")
        cb = fig.colorbar(sc, ax=ax, pad=0.01)
        cb.set_label("protrusion above the bond (um)", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    if wp is not None:
        hb, hz = wp["length_mm"] / 2.0, wp["width_mm"] / 2.0
        ax.add_patch(matplotlib.patches.Rectangle(
            (-hb * 1000, -hz * 1000), wp["length_mm"] * 1000, wp["width_mm"] * 1000,
            fill=False, ec="#2f6fb5", lw=2.0, zorder=5, label="workpiece footprint"))
        sw = plan.get("sweep_mm") or 0.0
        if sw:
            ax.add_patch(matplotlib.patches.Rectangle(
                ((-hb - sw) * 1000, -hz * 1000), (wp["length_mm"] + sw) * 1000,
                wp["width_mm"] * 1000, fill=False, ec="#b5432f", lw=1.2, ls="--",
                zorder=4, label="swept during the step"))
        ax.legend(fontsize=7, loc="upper right")
    _ar = (plan["grit_band_arc_mm"] / max(plan["grit_band_width_mm"], 1e-9))
    _eq = _ar <= 12.0
    ax.set_aspect("equal" if _eq else "auto")
    _note = ("" if _eq else
             "\n(band is %.0fx longer than wide - axes NOT to the same scale)" % _ar)
    ax.set_title("the dressed band, unwrapped: %d grits, each dot one grain%s"
                 % (len(frames), _note), fontsize=9)
    ax.set_xlabel("along the arc, from the block centre (um)", fontsize=8)
    ax.set_ylabel("across the face (um)", fontsize=8)
    ax.tick_params(labelsize=7)

    # ---------------- 3. the contact, in microns ----------------
    ax = fig.add_subplot(gs[1, 0:2])
    if wp is not None and frames:
        hb, hz = wp["length_mm"] / 2.0, wp["width_mm"] / 2.0
        near = [f for f in frames if abs(f[:, 1].mean()) < hb
                and abs(f[:, 2].mean()) < hz]
        bb = np.linspace(-hb, hb, 400)
        ax.plot(bb * 1000, (np.sqrt(np.maximum(R ** 2 - bb ** 2, 0.0)) - r_ground) * 1000,
                color="#4a4e54", lw=1.4, label="bond rim")
        for f in near[:400]:
            ax.plot(f[:, 1] * 1000, (f[:, 0] - r_ground) * 1000, ".",
                    ms=1.2, color="#2a9d5c")
        ax.axhline(0.0, color="#2f6fb5", lw=2.0, label="workpiece ground face")
        ae = plan.get("depth_of_cut_um") or 0.0
        if ae:
            ax.axhline(-ae, color="#b5432f", lw=1.4, ls="--",
                       label="after %.2f um infeed" % ae)
        ax.set_ylim(-max(plan["bond_clearance_um"] * 1.4, 1.0),
                    max(plan["bond_clearance_um"] * 0.4, 0.5))
        ax.legend(fontsize=7, loc="lower right")
        _n_cut = (sum(1 for f in near
                      if (f[:, 0].max() - r_ground) * 1000 >= -ae) if ae else 0)
        ax.set_title("AT t=0, under the block: %d grits, %d of them cut at this infeed"
                     % (len(near), _n_cut), fontsize=9)
        ax.set_xlabel("along the arc (um)", fontsize=8)
        ax.set_ylabel("height relative to the ground face (um)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.25)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "no workpiece in this deck", ha="center", fontsize=10)

    # ---------------- 4. engagement vs depth of cut ----------------
    ax = fig.add_subplot(gs[1, 2])
    clr = plan.get("swept_clearances_um")
    if clr is not None and len(clr):
        cs = np.sort(np.asarray(clr))
        ax.plot(cs, np.arange(1, len(cs) + 1), color="#2a9d5c", lw=1.6)
        ax.axvline(plan["bond_clearance_um"], color="#b5432f", ls=":", lw=1.4)
        ax.text(plan["bond_clearance_um"], len(cs) * 0.5, " bond hits", rotation=90,
                fontsize=7, color="#b5432f", va="center")
        ae = plan.get("depth_of_cut_um") or 0.0
        if ae:
            n = int((cs <= ae).sum())
            ax.axvline(ae, color="#2f6fb5", lw=1.6)
            ax.text(ae, len(cs) * 0.85, " %d grits cut" % n, fontsize=7,
                    color="#2f6fb5")
        ax.set_title("OVER THE WHOLE PASS: grits engaged vs depth of cut", fontsize=9)
        ax.set_xlabel("depth of cut (um)", fontsize=8)
        ax.set_ylabel("grits cutting", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.25)
    else:
        ax.axis("off")

    fig.suptitle(plan.get("title", "grinding model preview"), fontsize=11)
    return fig


def summary_text(plan: dict) -> str:
    """The numbers a user needs to sanity-check before committing to a build."""
    L = []
    a = L.append
    R = plan["outer_radius_mm"]
    a("WHEEL      %.1f mm diameter, %s, rim %.3f mm, %.3f mm thick"
      % (2 * R, "full wheel" if plan["full_wheel"]
         else "%.3f deg sector" % plan["sector_deg"],
         plan["rim_depth_mm"], plan["width_mm"]))
    if not plan["full_wheel"]:
        arc = plan["arc_length_mm"]
        sag = arc * arc / (8 * R)
        a("           arc %.3f mm, sagitta %.1f um = %.1f%% of the chord -> %s"
          % (arc, sag * 1000, 100 * sag / arc,
             "reads as an arc" if (plan["sector_deg"] >= 10 or sag > plan["rim_depth_mm"])
             else "WILL LOOK FLAT"))
    a("GRITS      %s placed, %.0f/mm2, dressed band %.3f x %.3f mm"
      % (format(plan["n_grits"], ","), plan["areal_density"],
         plan["grit_band_arc_mm"], plan["grit_band_width_mm"]))
    if plan.get("protrusion_um"):
        pa, gh = plan["protrusion_um"], plan.get("grain_height_um") or {}
        a("           stand proud of the bond by %.3f to %.3f um (mean %.3f, "
          "median %.3f)" % (pa["min"], pa["max"], pa["mean"], pa["median"]))
        if gh.get("n"):
            a("           grain height as measured %.3f to %.3f um (mean %.3f)"
              % (gh["min"], gh["max"], gh["mean"]))
    if plan["workpiece"] is not None:
        wp = plan["workpiece"]
        a("WORKPIECE  %g x %g x %g mm at theta = %.3f deg, %s C3D8R"
          % (wp["length_mm"], wp["width_mm"], wp["depth_mm"],
             plan["theta_workpiece_deg"], format(plan["n_workpiece_elements"], ",")))
        pos = plan.get("wp_position", "centred")
        gr = plan.get("grit_theta_range_deg") or (0.0, 0.0)
        a("           placed '%s': spans theta %.3f (entry) to %.3f deg, over %d "
          "grain%s" % (pos, plan.get("wp_entry_theta_deg", 0.0),
                       plan.get("wp_exit_theta_deg", 0.0),
                       plan.get("n_grits_under_block", 0),
                       "" if plan.get("n_grits_under_block") == 1 else "s"))
        grr = plan.get("grit_theta_reachable_deg") or gr
        a("           grit occupies theta %.3f to %.3f deg (%.3f to %.3f within the "
          "block's width)" % (gr[0], gr[1], grr[0], grr[1]))
        a("           the surface travels toward decreasing theta, so grains arrive "
          "from the high-theta end")
        if plan.get("wp_relocated"):
            a("           NOTE the footprint you asked for held no grit, so the block "
              "moved to the tallest grain it can reach")
        pu = plan.get("protrusion_under_block_um") or {}
        if pu.get("n"):
            a("           grains under it stand %.3f to %.3f um proud; standoff set "
              "to %.3f um" % (pu["min"], pu["max"], plan.get("standoff_um", 0.0)))
        a("           mesh %.2f cutting x %.2f axial x %.2f-%.2f depth um"
          % (plan["element_um"][0], plan["element_um"][1],
             plan["element_um"][2], plan["element_um"][3]))
        a("           ground face r = %.6f mm, tangent, bond clearance %.3f um"
          % (plan["ground_radius_mm"], plan["bond_clearance_um"]))
        clr = plan.get("swept_clearances_um")
        ceil = plan.get("depth_ceiling_um") or plan["bond_clearance_um"]
        floor = plan.get("first_contact_um")
        if clr is not None and len(clr):
            cs = np.sort(np.asarray(clr))
            q = lambda k: cs[min(k, len(cs) - 1)]
            a("ENGAGEMENT ae for 1 / 10 / 50 grits: %.2f / %.2f / %.2f um"
              % (q(0), q(9), q(49)))
            a("           so a usable depth of cut is between %.3f and %.3f um"
              % (max(floor or 0.0, 0.0), ceil))
        ae = plan.get("depth_of_cut_um") or 0.0
        if ae and ceil:
            if ae >= ceil:
                a("")
                a("*** DEPTH OF CUT TOO DEEP ***  ae = %.3f um but the face-to-bond gap"
                  % ae)
                a("    is only %.3f um, so the bond would hit the workpiece and the "
                  "build" % ceil)
                a("    will refuse. Set DEPTH_OF_CUT_UM to %.3f um or less."
                  % (0.85 * ceil))
            elif floor is not None and ae <= floor:
                a("")
                a("*** DEPTH OF CUT NEVER REACHES THE WORK ***  ae = %.3f um, but with"
                  % ae)
                a("    a %.3f um standoff the nearest grain in the swept band is %.3f "
                  "um" % (plan.get("standoff_um", 0.0), floor))
                a("    clear. The wheel would turn for the whole step and touch "
                  "nothing.")
                a("    Cut deeper than %.3f um, or lower CLEARANCE_UM." % floor)
            elif ae < 0.2 * ceil:
                a("NOTE       ae = %.3f um is only %.0f%% of the %.3f um gap; few "
                  "grits will reach the work" % (ae, 100 * ae / ceil, ceil))
    if plan.get("cost"):
        c = plan["cost"]
        a("RUN        dt %.3e s, %s increments, %.1f h on 8 cores / %.1f h on 4"
          % (c["stable_dt_s"], format(int(c["increments"]), ","),
             c["est_hours"]["8"], c["est_hours"]["4"]))
    a("SIZE       roughly %.0f MB of .inp" % plan["estimated_mb"])
    for w in plan.get("warnings", []) + plan.get("notes", []):
        a("  note: %s" % w)
    return "\n".join(L)
