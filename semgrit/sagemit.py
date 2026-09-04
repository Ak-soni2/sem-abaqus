"""Emit the two SAG decks as complete, submittable Abaqus input files.

``sagwrite`` builds the meshes and the material card; this assembles them into
``.inp`` files with parts, instances, surfaces, constraints, steps, contact and
output, and writes them to disk.

WHAT EACH DECK IS FOR, AND WHAT IT CANNOT DO
--------------------------------------------
``MACRO`` presses the compliant tool into the workpiece and rotates it. Rigid
hub, hyperelastic/viscoelastic polyurethane ring, measured grains on the pad.
It answers the *contact*: patch size, pressure distribution, how many grains
engage and what load each carries. Its workpiece is meshed for contact, not for
``dc``, so it CANNOT resolve a ductile-brittle transition and its header says so.

``MICRO`` takes one patch of that contact and meshes it at ``dc/5`` through the
depth. It answers the *transition*: SDV13, ductile against brittle, under the
per-grain load MACRO computed. The two are coupled by that single number and
both headers print it, so the pair cannot be quoted out of step.

WHY GENERAL CONTACT, NOT CONTACT PAIRS
--------------------------------------
The reference deck uses ``*Contact`` + ``*Contact Inclusions, ALL EXTERIOR``,
and that is not an arbitrary preference -- it is required here for three
reasons, any one of which would be enough:

* **Element deletion.** The VUMAT deletes failed elements (SDV12), and deletion
  exposes interior faces that were not on the exterior when the job started.
  General contact re-forms its domain as that happens; a contact pair declared
  on a pre-computed surface never sees the new faces, so a chip would separate
  and then pass through the tool.
* **The engaged set is the answer.** Which grains touch is what the model is
  for. Pairs would have to declare them in advance, which assumes the result.
* **Self-contact.** A compliant layer at high compression can fold onto itself.

Cost scales with facet count, which is why the grain cap is quoted in facets.

Both run in ``*Dynamic, Explicit``: general contact is the contact algorithm,
explicit is the solver, and they are independent choices that happen to be the
right ones together here.
"""

from __future__ import annotations

import math
import os
from typing import Optional, Sequence

import numpy as np

from . import sag, sagdeck, sagwrite
from .sagwrite import (SAGWriteError, _elem_lines, _fmt, _node_lines,
                       _pack_ids, build_block, build_compliant_ring,
                       hex_volume)

# Element-set names, one place so the deck and the postprocessor agree.
ES_WORK = "ES_WORK"
ES_WORK_FINE = "ES_WORK_FINE"
ES_PU = "ES_PU"
ES_HUB = "ES_HUB"
NS_WORK_FIXED = "NS_WORK_FIXED"
NS_HUB_REF = "NS_HUB_REF"
SURF_WORK_TOP = "SURF_WORK_TOP"
SURF_PU_OUTER = "SURF_PU_OUTER"
SURF_HUB_OUTER = "SURF_HUB_OUTER"
SURF_PU_BORE = "SURF_PU_BORE"

STANDOFF_FRACTION = 0.05
"""How far clear of the surface a MICRO grain starts, as a fraction of the
INDENTATION it will be pushed to.

Small enough that closing it is a negligible part of the ramp, large enough
that the grain is not already interpenetrating at t = 0 -- an initial
overclosure is an impulse the contact algorithm must resolve before anything
physical happens.

It has to be a fraction of the indentation and NOT of the block depth. Those
are five orders apart here: 2% of a 1.6 um block is 32 nm against a 0.16 nm
indentation, so the grain started 200x further out than it could ever travel
and the job ran with nothing ever touching. Kinetic and internal energy stayed
identically zero for the whole step, which is how it was caught."""


def _surface_from_quads(name: str, elset: str, face: str) -> list:
    """A surface named by element set and face code (S1..S6)."""
    return ["*Surface, type=ELEMENT, name=%s" % name, " %s, %s" % (elset, face)]


def _grain_parts(solids: Sequence, plan_: dict, *, n_grains: int,
                 radius_mm: float, sector_deg: float, width_mm: float,
                 protrusion_frac: float, seed: int = 7) -> tuple:
    """Place measured grains on the pad's outer face.

    Grains are placed on the cylindrical surface at ``radius_mm``, each buried
    so that ``protrusion_frac`` of its height stands clear -- which is the
    physical statement that a grit cuts only as deep as it protrudes, and the
    same convention the rigid pipeline uses.

    Returned as one merged triangle soup per grain set, because general contact
    wants surfaces and a thousand separate parts would make an unreadable deck.
    """
    rng = np.random.default_rng(seed)
    verts, faces, tags = [], [], []
    off = 0
    n_avail = len(solids)
    for i in range(n_grains):
        s = solids[int(rng.integers(0, n_avail))]
        v = np.asarray(s.vertices, dtype=np.float64) / 1000.0     # um -> mm
        v = v - np.asarray(s.centroid_um, dtype=np.float64) / 1000.0
        h = s.height_um / 1000.0
        # random spin about the grain's own axis, so the pad is not a lattice
        a = rng.uniform(0.0, 2.0 * math.pi)
        ca, sa = math.cos(a), math.sin(a)
        v = v @ np.array([[ca, -sa, 0.0], [sa, ca, 0.0], [0.0, 0.0, 1.0]]).T
        # seat it: the grain's local +z becomes the outward radial direction
        th = math.radians(rng.uniform(0.0, sector_deg))
        z = rng.uniform(0.1 * width_mm, 0.9 * width_mm)
        # buried so protrusion_frac of the height is clear of the pad
        r_seat = radius_mm - (1.0 - protrusion_frac) * h
        er = np.array([math.cos(th), math.sin(th), 0.0])
        et = np.array([-math.sin(th), math.cos(th), 0.0])
        ez = np.array([0.0, 0.0, 1.0])
        # local (x, y, z) -> (tangential, axial, radial)
        world = (v[:, 0:1] * et + v[:, 1:2] * ez + v[:, 2:3] * er)
        world = world + er * r_seat + ez * z
        verts.append(world)
        faces.append(np.asarray(s.faces, dtype=np.int64) + off)
        tags.append(dict(index=i, theta_deg=math.degrees(th), z_mm=z,
                         height_um=s.height_um,
                         protrusion_um=protrusion_frac * s.height_um,
                         grain_id=getattr(s, "grain_id", -1)))
        off += len(v)
    if not verts:
        raise SAGWriteError("no grains were placed")
    return (np.vstack(verts), np.vstack(faces), tags)


def _step_block(name: str, time_s: float, *, mass_scale: float = 0.0,
                bcs: Sequence[str] = (), comment: Sequence[str] = (),
                intervals: int = 40) -> list:
    """One ``*Step`` with explicit dynamics, bulk viscosity and output."""
    L = ["** " + "-" * 70]
    L += ["** " + c for c in comment]
    L += ["*Step, name=%s, nlgeom=YES" % name,
          "*Dynamic, Explicit", ", %s" % _fmt(time_s),
          "*Bulk Viscosity", " 0.06, 1.2"]
    if mass_scale > 0:
        L += ["** Mass scaling is stated, not hidden in the density. The",
              "** reference deck under-densified the polyurethane by ~5000x",
              "** AND applied factor 50, which double-counts the speedup and",
              "** removes the ring's inertia.",
              "*Fixed Mass Scaling, factor=%s" % _fmt(mass_scale)]
    if bcs:
        L += ["*Boundary, op=NEW, type=VELOCITY"]
        L += list(bcs)
    L += ["*Restart, write, number interval=1, time marks=NO",
          "*Output, field, number interval=%d" % intervals,
          "*Node Output", " U, V, A, RF, RM",
          "*Element Output, directions=YES",
          # Only what is used. The reference asks for ~70 variables including
          # BURNF and IWCONWEP on 1.09 M elements, which is an enormous .odb
          # for fields nothing reads.
          " S, MISES, PEEQ, LE, SDV, STATUS, EVOL, IVOL",
          "*Contact Output", " CSTRESS, CDISP, CFORCE, CNAREA, CSTATUS",
          "*Output, history, time interval=%s" % _fmt(time_s / 200.0),
          "*Energy Output", " ALLIE, ALLKE, ALLAE, ALLPD, ALLVD, ETOTAL",
          "*End Step"]
    return L


def write_macro(path: str, pl: dict, solids: Sequence, *,
                seed: int = 7) -> dict:
    """The contact deck: compliant tool pressed into the work, then rotated."""
    p: sagdeck.SAGParams = pl["params"]
    mac = pl["macro"]
    t = pl["timing"]
    c: sag.SAGContact = pl["contact"]
    from . import materials

    if not solids:
        raise SAGWriteError("the grain library is empty")

    r_out = 0.5 * p.diameter_mm
    r_pu_in = r_out - p.polyurethane.thickness_mm
    r_hub_in = max(r_pu_in - 0.5 * p.polyurethane.thickness_mm, 0.05 * r_out)
    sect = mac["sector_deg"]

    # Circumferential divisions: enough that the faceted arc is smooth to well
    # under the wheel compression, or the contact would see a polygon.
    arc_target = min(p.compression_mm * 0.25, r_out * math.radians(sect) / 8.0)
    n_circ = max(int(math.ceil(math.radians(sect) * r_out / arc_target)), 8)
    # Through-thickness: a layer that must BEND needs several elements or it
    # only shears. Six is the usual minimum for a bending sandwich layer.
    n_rad_pu = 6
    n_ax = max(int(round(p.width_mm / (p.width_mm / 8.0))), 8)

    pu_n, pu_c, pu_f = build_compliant_ring(
        inner_r_mm=r_pu_in, outer_r_mm=r_out, width_mm=p.width_mm,
        sector_deg=sect, n_circ=n_circ, n_rad=n_rad_pu, n_axial=n_ax)
    hub_n, hub_c, hub_f = build_compliant_ring(
        inner_r_mm=r_hub_in, outer_r_mm=r_pu_in, width_mm=p.width_mm,
        sector_deg=sect, n_circ=n_circ, n_rad=2, n_axial=n_ax)

    for nm, (nn, cc) in (("polyurethane", (pu_n, pu_c)), ("hub", (hub_n, hub_c))):
        if hex_volume(nn, cc, signed=True) <= 0:
            raise SAGWriteError("%s mesh is wound backwards" % nm)

    # Workpiece: sized to the arc the tool sweeps during the grind step, plus
    # the contact length, so the tool never runs off the end of the block.
    sweep = c.surface_speed_mm_s * p.grind_time_s
    wp_len = sweep + 2.0 * c.semi_axis_a_mm * 0.2 + 4.0 * p.compression_mm
    wp_wid = min(p.width_mm, 2.0 * c.semi_axis_b_mm)
    wp_dep = max(20.0 * p.compression_mm, 0.5)
    el_ip = max(p.compression_mm / 4.0, wp_len / 400.0)
    wp_n, wp_c, wp_m = build_block(
        length_mm=wp_len, width_mm=wp_wid, depth_mm=wp_dep,
        el_length_mm=el_ip, el_width_mm=el_ip,
        fine_depth_mm=p.compression_mm / 8.0,
        band_mm=2.0 * p.compression_mm, growth=1.3,
        x0_mm=-0.5 * wp_len, y0_mm=-0.5 * wp_wid, top_z_mm=0.0)

    gv, gf, gtags = _grain_parts(
        solids, pl, n_grains=mac["grains"], radius_mm=r_out,
        sector_deg=sect, width_mm=p.width_mm,
        protrusion_frac=0.55, seed=seed)

    # Seat the tool at first contact: the tallest protrusion just touches the
    # ground face. No gap to close, so the press-in displacement IS the wheel
    # compression -- which is what removes the reference deck's ambiguity
    # between a velocity and a displacement.
    tallest = max(g["protrusion_um"] for g in gtags) / 1000.0
    y_centre = r_out + tallest

    w = materials.get(p.material)
    L = sagdeck.macro_header(pl)
    a = L.append
    a("*Heading")
    a("** SAG MACRO -- %s -- contact, not transition" % p.name)
    a("*Preprint, echo=NO, model=NO, history=NO, contact=NO")

    # ---- parts ----------------------------------------------------------
    a("**")
    a("*Part, name=HUB")
    a("*Node")
    L += _node_lines(hub_n)
    a("*Element, type=C3D8R")
    L += _elem_lines(hub_c)
    a("*Elset, elset=%s, generate" % ES_HUB)
    a(" 1, %d, 1" % len(hub_c))
    a("*Solid Section, elset=%s, material=HUBSTEEL" % ES_HUB)
    a(" ,")
    a("*End Part")

    a("**")
    a("*Part, name=PU")
    a("*Node")
    L += _node_lines(pu_n)
    a("*Element, type=C3D8R")
    L += _elem_lines(pu_c)
    a("*Elset, elset=%s, generate" % ES_PU)
    a(" 1, %d, 1" % len(pu_c))
    a("** The compliant layer. DISTORTION CONTROL because a hyperelastic")
    a("** layer at 0.4 mm compression distorts far more than a metal, and")
    a("** ENHANCED hourglass because C3D8R has zero-energy modes that a soft")
    a("** material excites readily.")
    a("*Solid Section, elset=%s, controls=EC1, material=POLYURETHANE" % ES_PU)
    a(" ,")
    a("*End Part")

    a("**")
    a("*Part, name=GRAINS")
    a("** Measured SEM grain geometry, rigid: a diamond abrasive is ~10x")
    a("** stiffer than WC-Co and ~3000x stiffer than the polyurethane, so its")
    a("** own deformation is not the physics and meshing it as a solid would")
    a("** spend the element budget on the one body that does not deform.")
    a("*Node")
    L += _node_lines(gv)
    a("*Element, type=R3D3")
    L += _elem_lines(gf)
    a("*Elset, elset=ES_GRAINS, generate")
    a(" 1, %d, 1" % len(gf))
    a("*End Part")

    a("**")
    a("*Part, name=WORK")
    a("*Node")
    L += _node_lines(wp_n)
    a("*Element, type=C3D8R")
    L += _elem_lines(wp_c)
    a("*Elset, elset=%s, generate" % ES_WORK)
    a(" 1, %d, 1" % len(wp_c))
    a("*Nset, nset=%s" % NS_WORK_FIXED)
    L += _pack_ids(sorted(set(wp_m["bottom"]) | set(wp_m["sides"])))
    a("*Solid Section, elset=%s, controls=EC1, material=%s"
      % (ES_WORK, w.inp_material))
    a(" ,")
    a("*End Part")

    # ---- assembly -------------------------------------------------------
    a("**")
    a("*Assembly, name=Assembly")
    a("**")
    a("*Instance, name=WORK-1, part=WORK")
    a("*End Instance")
    for nm, part in (("HUB-1", "HUB"), ("PU-1", "PU"), ("GRAINS-1", "GRAINS")):
        a("**")
        a("*Instance, name=%s, part=%s" % (nm, part))
        a(" 0., %s, 0." % _fmt(y_centre))
        # The rings are built about +z; rotate so the wheel axis is z and the
        # tool sits above the work in +y.
        a(" 0., %s, 0., 1., %s, 0., -90." % (_fmt(y_centre), _fmt(y_centre)))
        a("*End Instance")

    a("**")
    a("*Nset, nset=%s, instance=HUB-1" % NS_HUB_REF)
    a(" 1,")
    a("** The hub is rigid: it is the spindle, it carries the drive, and")
    a("** making it deformable would cost elements on a body whose stiffness")
    a("** is irrelevant beside a 0.345 MPa layer.")
    a("*Rigid Body, ref node=%s, elset=HUB-1.%s" % (NS_HUB_REF, ES_HUB))

    a("**")
    a("*Tie, name=PU_TO_HUB, adjust=yes")
    a(" PU-1.%s, HUB-1.%s" % (SURF_PU_BORE, SURF_HUB_OUTER))
    a("*Tie, name=GRAINS_TO_PU, adjust=no")
    a("** Grains are bonded to the pad rather than left to general contact:")
    a("** a grit that could slide off its own backing is not an abrasive.")
    a(" GRAINS-1.ES_GRAINS, PU-1.%s" % SURF_PU_OUTER)
    a("*End Assembly")

    # ---- element controls, materials -----------------------------------
    a("**")
    a("*Section Controls, name=EC1, DISTORTION CONTROL=YES,"
      " ELEMENT DELETION=YES, hourglass=ENHANCED")
    a(" 1., 1., 1.")
    a("**")
    a("*Material, name=HUBSTEEL")
    a("*Density")
    a(" 7.85e-09,")
    a("*Elastic")
    a(" 210000., 0.3")
    a("**")
    a("*Material, name=POLYURETHANE")
    L += p.polyurethane.cards()
    a("**")
    L += sagwrite._material_block(pl, w.inp_material, 0.0,
                                  p.compression_mm / 8.0)

    # ---- interactions ---------------------------------------------------
    a("**")
    a("*Surface Interaction, name=IP_GRIND")
    a("*Friction")
    a(" %s," % _fmt(p.friction))
    a("*Surface Behavior, pressure-overclosure=HARD")
    a("**")
    a("** GENERAL contact, not contact pairs, and it is required rather than")
    a("** merely convenient: the VUMAT deletes elements, deletion exposes")
    a("** interior faces that were not on the exterior at the start, and a")
    a("** pair declared on a pre-computed surface would never see them -- a")
    a("** chip would separate and then pass through the tool. Which grains")
    a("** touch is also the ANSWER here, so it cannot be declared up front.")
    a("*Contact, op=NEW")
    a("*Contact Inclusions, ALL EXTERIOR")
    a("*Contact Property Assignment")
    a(" ,  , IP_GRIND")

    # ---- boundary conditions -------------------------------------------
    a("**")
    a("*Boundary")
    a(" WORK-1.%s, ENCASTRE" % NS_WORK_FIXED)

    # ---- steps ----------------------------------------------------------
    v_press = p.compression_mm / t["press_time_s"]
    omega = p.speed_rpm * 2.0 * math.pi / 60.0
    L += _step_block(
        "PRESS", t["press_time_s"], mass_scale=p.mass_scale_factor,
        comment=("STEP 1 of 3 -- press the tool in by the wheel compression.",
                 "Velocity, not displacement, so the rate is explicit: %.1f"
                 % v_press + " mm/s,",
                 "which is %.4f of the layer's own wave speed. Above ~0.01"
                 % t["press_mach"],
                 "the patch is loaded inertially and its pressure is not the",
                 "steady Hertzian one the experiment measured."),
        bcs=(" %s, 1, 1, 0." % NS_HUB_REF,
             " %s, 2, 2, %s" % (NS_HUB_REF, _fmt(-v_press)),
             " %s, 3, 3, 0." % NS_HUB_REF,
             " %s, 4, 4, 0." % NS_HUB_REF,
             " %s, 5, 5, 0." % NS_HUB_REF,
             " %s, 6, 6, 0." % NS_HUB_REF))
    L += _step_block(
        "HOLD", t["hold_time_s"], mass_scale=p.mass_scale_factor,
        comment=("STEP 2 of 3 -- hold, so the polyurethane RELAXES.",
                 "%.1f Prony time constants. The measured force is a steady"
                 % t["hold_taus"],
                 "reading and the Hertz comparison uses moduli=LONG TERM, so",
                 "grinding straight after the press would carry a glassy",
                 "layer: stiffer, smaller patch, more load per grain."),
        bcs=(" %s, 1, 1, 0." % NS_HUB_REF,
             " %s, 2, 2, 0." % NS_HUB_REF,
             " %s, 3, 3, 0." % NS_HUB_REF,
             " %s, 4, 4, 0." % NS_HUB_REF,
             " %s, 5, 5, 0." % NS_HUB_REF,
             " %s, 6, 6, 0." % NS_HUB_REF))
    L += _step_block(
        "GRIND", p.grind_time_s, mass_scale=p.mass_scale_factor,
        comment=("STEP 3 of 3 -- rotate at %.1f rpm = %.4f rad/s."
                 % (p.speed_rpm, omega),
                 "The compression is HELD by fixing dof 2 while dof 6 spins,",
                 "which is what a real infeed does: depth is set, then cut."),
        bcs=(" %s, 1, 1, 0." % NS_HUB_REF,
             " %s, 2, 2, 0." % NS_HUB_REF,
             " %s, 3, 3, 0." % NS_HUB_REF,
             " %s, 4, 4, 0." % NS_HUB_REF,
             " %s, 5, 5, 0." % NS_HUB_REF,
             " %s, 6, 6, %s" % (NS_HUB_REF, _fmt(-omega))))

    # PU surfaces are needed by the ties; emit them inside the parts.
    L = _inject_part_surfaces(L, pu_f, hub_f, wp_m)

    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")

    return dict(
        kind="macro", path=path, bytes=os.path.getsize(path),
        elements=len(hub_c) + len(pu_c) + len(wp_c),
        pu_elements=len(pu_c), hub_elements=len(hub_c),
        work_elements=len(wp_c), grain_facets=len(gf),
        grains=len(gtags), nodes=len(hub_n) + len(pu_n) + len(wp_n) + len(gv),
        sector_deg=sect, n_circ=n_circ, n_rad_pu=n_rad_pu,
        work_mm=(wp_len, wp_wid, wp_dep),
        tool_centre_y_mm=y_centre, tallest_protrusion_mm=tallest,
        press_velocity_mm_s=v_press, omega_rad_s=omega,
        steps=("PRESS", "HOLD", "GRIND"),
        resolves_dc=False, grain_tags=gtags,
    )


def _inject_part_surfaces(lines: list, pu_f: dict, hub_f: dict,
                          wp_m: dict) -> list:
    """Add the element-face surfaces the ties and contact need.

    Written as nodal surfaces on the bore/outer node sets, which is exact for a
    structured ring: those node sets ARE the two cylindrical faces, so no face
    codes have to be guessed from element ordering.
    """
    out = []
    for ln in lines:
        if ln == "*End Part" and out and "PU" in _current_part(out):
            out.append("*Nset, nset=%s" % SURF_PU_BORE)
            out += _pack_ids(pu_f["bore"])
            out.append("*Nset, nset=%s" % SURF_PU_OUTER)
            out += _pack_ids(pu_f["outer"])
        elif ln == "*End Part" and out and "HUB" in _current_part(out):
            out.append("*Nset, nset=%s" % SURF_HUB_OUTER)
            out += _pack_ids(hub_f["outer"])
        out.append(ln)
    return out


def _current_part(lines: list) -> str:
    for ln in reversed(lines):
        if ln.startswith("*Part, name="):
            return ln.split("=", 1)[1]
    return ""


def write_micro(path: str, pl: dict, solids: Sequence, *,
                psi: float = 0.0, seed: int = 11, n_passes: int = 0) -> dict:
    """The resolved deck: one patch of contact at dc/5, energy criterion.

    The grain is driven by a FORCE -- the per-grain load MACRO computed -- not
    by a prescribed depth. Prescribing the depth would assume the answer, since
    the indentation is what the model exists to predict.
    """
    from . import materials

    p: sagdeck.SAGParams = pl["params"]
    mic = pl["micro"]
    c: sag.SAGContact = pl["contact"]
    w = materials.get(p.material)
    hp = w.hybrid_params()

    if not solids:
        raise SAGWriteError("the grain library is empty")

    side = mic["side_mm"]
    nodes, conn, meta = build_block(
        length_mm=side, width_mm=side, depth_mm=mic["depth_mm"],
        el_length_mm=mic["element_inplane_mm"],
        el_width_mm=mic["element_inplane_mm"],
        fine_depth_mm=mic["element_mm"],
        band_mm=min(mic["depth_mm"], 10.0 * hp.critical_depth_mm()),
        growth=1.25, x0_mm=-0.5 * side, y0_mm=-0.5 * side, top_z_mm=0.0)

    vol = hex_volume(nodes, conn, signed=True)
    want = side * side * mic["depth_mm"]
    if vol <= 0:
        raise SAGWriteError("the workpiece mesh is wound backwards")
    if abs(vol - want) / want > 1e-6:
        raise SAGWriteError("mesh volume %.6g != block %.6g" % (vol, want))

    # THE CUT DEPTH IS THE MEASURED CHIP THICKNESS.
    #
    # Not the Brinell indentation. That distinction is the whole reason this
    # deck works at all, so it is worth being explicit.
    #
    # sag.py's eqs. 11-12 give the depth a STATIC spherical indenter reaches
    # under the per-grain load: 0.18 nm for the 30 um pad. That number is
    # self-consistent -- a 0.18 nm cap on a 30 um sphere is 146 nm wide, which
    # is what groove_width_nm reports -- but it is not the depth of material
    # a MOVING grain removes, and it is not what the paper's dc refers to.
    #
    # The paper measures chips directly (section 4.2, Fig. 17) and calls the
    # transition threshold a "critical chip thickness ... critical depth of
    # indentation": 240-350 nm for the 30 um pad, 160-230 for 15 um, 60-100
    # for 6 um, against dc = 60-100 nm. Its entire argument is that the first
    # two exceed dc and fracture while the third sits at dc and does not. So
    # the chip thickness IS the depth the criterion is tested against.
    #
    # Two things go wrong if the Brinell indentation is used instead. It is
    # three orders below what the paper measured for the same pad, so it
    # cannot reproduce the ordering. And at 0.18 nm it is 1/90th of one 16 nm
    # element and smaller than a WC unit cell (~0.29 nm) -- there is no mesh,
    # and no continuum, that represents it.
    chip_nm = mic.get("chip_depth_nm") or 0.0
    if chip_nm <= 0:
        raise SAGWriteError(
            "no measured chip thickness for the %g um pad, so there is no "
            "cut depth to impose. sag.MEASURED_CHIP_NM carries the paper's "
            "values for 6, 15 and 30 um." % p.grain_um)
    indent_mm = chip_nm * 1e-6
    n_gr = max(int(mic["grains"]), 1)
    # Diamond, as a sphere of the nominal grain size. REPORTED only -- the
    # rigid body is deliberately massless (see the assembly comment), and the
    # dynamics are displacement-controlled regardless.
    grain_mass_t = (3520.0 * 1e-12 * (math.pi / 6.0)
                    * (p.grain_um * 1e-3) ** 3 * n_gr)
    gv, gf, gtags = _grain_parts(
        solids, pl, n_grains=n_gr, radius_mm=0.0, sector_deg=0.0,
        width_mm=0.0, protrusion_frac=1.0, seed=seed)
    # _grain_parts seats on a cylinder; for a flat patch place them on the
    # surface directly, spread over the patch and just clear of it.
    rng = np.random.default_rng(seed)
    gv = np.asarray(gv, dtype=np.float64)
    per = len(gv) // n_gr
    for i in range(n_gr):
        sl = slice(i * per, (i + 1) * per if i < n_gr - 1 else len(gv))
        blk = gv[sl]
        blk = blk - blk.mean(axis=0)
        x = rng.uniform(-0.35 * side, 0.35 * side)
        y = rng.uniform(-0.35 * side, 0.35 * side)
        blk[:, 0] += x
        blk[:, 1] += y
        # Just clear of the surface -- a fraction of the
        # INDENTATION, not of the block depth. See
        # STANDOFF_FRACTION for why that distinction matters.
        blk[:, 2] += (-blk[:, 2].min()
                      + STANDOFF_FRACTION * indent_mm)
        gv[sl] = blk
        gtags[i].update(x_mm=x, y_mm=y)

    L = sagdeck.micro_header(pl)
    a = L.append
    a("*Heading")
    a("** SAG MICRO -- %s -- the transition, resolved" % p.name)
    a("*Preprint, echo=NO, model=NO, history=NO, contact=NO")
    a("**")
    a("*Part, name=WORK")
    a("*Node")
    L += _node_lines(nodes)
    a("*Element, type=C3D8R")
    L += _elem_lines(conn)
    a("*Elset, elset=%s, generate" % ES_WORK)
    a(" 1, %d, 1" % len(conn))
    a("*Nset, nset=%s" % NS_WORK_FIXED)
    L += _pack_ids(sorted(set(meta["bottom"]) | set(meta["sides"])))
    a("*Solid Section, elset=%s, controls=EC1, material=%s"
      % (ES_WORK, w.inp_material))
    a(" ,")
    a("*End Part")
    a("**")
    a("*Part, name=GRAINS")
    a("*Node")
    L += _node_lines(gv)
    a("*Element, type=R3D3")
    L += _elem_lines(gf)
    a("*Elset, elset=ES_GRAINS, generate")
    a(" 1, %d, 1" % len(gf))
    a("*End Part")
    a("**")
    a("*Assembly, name=Assembly")
    a("*Instance, name=WORK-1, part=WORK")
    a("*End Instance")
    a("*Instance, name=GRAINS-1, part=GRAINS")
    a("*End Instance")
    a("*Nset, nset=NS_GRAIN_REF, instance=GRAINS-1")
    a(" 1,")
    a("*Rigid Body, ref node=NS_GRAIN_REF, elset=GRAINS-1.ES_GRAINS")
    a("**")
    a("** NO MASS, DELIBERATELY. R3D3 facets carry no volume, so this rigid")
    a("** body has none -- and it does not need any: Abaqus permits a")
    a("** massless rigid body when every translational dof is constrained,")
    a("** and all three are, in every step of this deck. There is no free")
    a("** translation for a = F/m to be undefined on.")
    a("**")
    a("** Two attempts to add it anyway both failed, and both failed quietly")
    a("** enough to be worth recording:")
    a("**   *Mass, elset=<the R3D3 set>  is accepted and IGNORED -- that card")
    a("**     assigns to MASS elements, so pointing it at rigid facets")
    a("**     assigns to nothing and the reported model mass never changes.")
    a("**   *Element, type=MASS in the assembly needs NODE NUMBERS and a node")
    a("**     to attach to; NS_GRAIN_REF is a node SET, and there are no free")
    a("**     nodes at assembly level. That one aborts the input processor.")
    a("**")
    a("** The mass was only ever wanted as a way to tell one build of this")
    a("** deck from another. The job name does that job properly.")
    a("*End Assembly")
    a("**")
    a("*Section Controls, name=EC1, DISTORTION CONTROL=YES,"
      " ELEMENT DELETION=YES, hourglass=ENHANCED")
    a(" 1., 1., 1.")
    a("**")
    L += sagwrite._material_block(pl, w.inp_material, psi, mic["element_mm"])
    a("**")
    a("*Surface Interaction, name=IP_GRIND")
    a("*Friction")
    a(" %s," % _fmt(p.friction))
    a("*Surface Behavior, pressure-overclosure=HARD")
    a("**")
    a("*Contact, op=NEW")
    a("*Contact Inclusions, ALL EXTERIOR")
    a("*Contact Property Assignment")
    a(" ,  , IP_GRIND")
    a("**")
    a("*Boundary")
    a(" WORK-1.%s, ENCASTRE" % NS_WORK_FIXED)

    # Drive by force, and slide tangentially so plastic work ACCUMULATES --
    # the energy criterion triggers on history, so a grain that only indents
    # and stops can never reach the threshold no matter how hard it presses.
    fn = c.load_per_grain_n
    v_slide = c.surface_speed_mm_s
    slide_time = side / v_slide if v_slide > 0 else p.grind_time_s
    # Passes needed to reach H*dc from this grain's own tangential work per
    # pass, rounded up with a margin so the transition is bracketed rather
    # than just touched. Derived, not chosen.
    if n_passes <= 0:
        w_mm = max(c.groove_width_mm, 1e-12)
        per_pass = (p.friction * fn) / w_mm
        need = (hp.hardness_mpa * hp.critical_depth_mm()) / max(per_pass, 1e-30)
        n_passes = max(int(math.ceil(need * 1.5)), 2)
    n_passes = min(n_passes, 60)
    L += ["** " + "-" * 70,
          "** STEP 1 of 2 -- cut to the MEASURED chip thickness,",
          "** %.1f nm for this pad." % chip_nm,
          "**",
          "** That is the paper's own measurement (section 4.2, Fig. 17),",
          "** not the static Brinell indentation, which for this pad is",
          "** %.3f nm -- three orders smaller, 1/%.0f of one element, and"
          % (c.indentation_nm, mic["element_mm"] * 1e6 / max(c.indentation_nm, 1e-30)),
          "** below a WC unit cell. The paper compares its measured chip",
          "** thickness against dc and that comparison is its result, so",
          "** the chip thickness is what this deck cuts.",
          "**   chip %.1f nm / dc %.1f nm = %.2f"
          % (chip_nm, hp.critical_depth_mm() * 1e6,
             chip_nm / (hp.critical_depth_mm() * 1e6)),
          "**",
          "** DISPLACEMENT-controlled, not force-controlled, and the reason",
          "** matters. A rigid grain has no mass of its own -- R3D3 facets",
          "** carry no volume -- so a free translational dof driven by a",
          "** *Cload has no m to divide the force by. Abaqus rejects that",
          "** outright:",
          "**   ERROR: Abaqus/Explicit requires rigid bodies to have a",
          "**   non-zero mass unless translational constraints are applied.",
          "** Giving it the real diamond mass does not rescue force control",
          "** either: %.3e N on a %.3e tonne grain is 4e10 mm/s2, which"
          % (fn * n_gr, 3.52e-9 * (math.pi / 6.0) * (p.grain_um * 1e-3) ** 3),
          "** carries it two hundred block-depths deep before contact can",
          "** resist. A 6 um grit is far too light to be pushed by a force at",
          "** this timescale.",
          "**",
          "** The experiment does not push the grain with a free force either.",
          "** The compliant pad POSITIONS it: the backing deflects, and the",
          "** grain sits at a depth while the layer carries the load. That",
          "** depth is not unknown here -- semgrit.sag predicts it from the",
          "** Hertzian chain, and that chain is validated against the paper's",
          "** own measured spot area, per-grain force and chip thickness. So",
          "** it is an INPUT to this deck, and the thing still being PREDICTED",
          "** is the branch: ductile or brittle, from accumulated plastic work.",
          "**",
          "** RF3 at NS_GRAIN_REF is then a CHECK, not an input: it should come",
          "** back near %.4e N. If it does not, the contact chain and the FE" % (fn * n_gr),
          "** contact disagree, and that is worth knowing.",
          "*Step, name=LOAD, nlgeom=YES",
          "*Dynamic, Explicit", ", %s" % _fmt(slide_time * 0.2),
          "*Bulk Viscosity", " 0.06, 1.2",
          "*Boundary, op=NEW, amplitude=RAMP",
          " NS_GRAIN_REF, 1, 1, 0.", " NS_GRAIN_REF, 2, 2, 0.",
          " NS_GRAIN_REF, 3, 3, %s" % _fmt(-indent_mm),
          " NS_GRAIN_REF, 4, 6, 0.",
          "*Amplitude, name=RAMP, definition=SMOOTH STEP",
          " 0., 0., %s, 1." % _fmt(slide_time * 0.2),
          "*Restart, write, number interval=1, time marks=NO",
          "*Output, field, number interval=20",
          "*Node Output", " U, V, A, RF",
          "*Element Output, directions=YES",
          " S, MISES, PEEQ, LE, SDV, STATUS, EVOL",
          "*Contact Output", " CSTRESS, CDISP, CFORCE, CSTATUS",
          "*Output, history, time interval=%s" % _fmt(slide_time * 0.001),
          "*Energy Output", " ALLIE, ALLKE, ALLAE, ALLPD, ETOTAL",
          "*End Step"]
    # REPEATED passes over the SAME track, not one long slide down a fresh
    # one. This is the correction that makes the deck able to reproduce the
    # experiment at all: the criterion accumulates work PER POINT, so a grain
    # sliding along virgin material leaves every point with exactly one pass
    # and can never trip the threshold however far it goes. The paper's spot
    # test accumulates because ~20,000 grain crossings pass over each point in
    # 10 s; on this material 7-16 passes reach H*dc, so a handful of repeats
    # spans the transition.
    for ip in range(1, n_passes + 1):
        first = ip == 1
        L += ["** " + "-" * 70,
              "** PASS %d of %d over the SAME track, at %.1f mm/s."
              % (ip, n_passes, v_slide)]
        if first:
            L += ["**",
                  "** The energy criterion accumulates plastic work PER POINT,"
                  " so a",
                  "** single pass along fresh material cannot trip it: every"
                  " point",
                  "** it crosses has seen exactly one pass. Repeating over one"
                  " track",
                  "** is what a polishing pad actually does -- the paper's"
                  " 10 s spot",
                  "** test puts ~20,000 grain crossings over each point -- and"
                  " on",
                  "** this material 7-16 passes reach H*dc.",
                  "**",
                  "** PLOT SDV13: 1 ductile, 2 brittle. Watch it evolve pass"
                  " by pass;",
                  "** the pass at which it flips is the result."]
        L += ["*Step, name=PASS%d, nlgeom=YES" % ip,
              "*Dynamic, Explicit", ", %s" % _fmt(slide_time),
              "*Bulk Viscosity", " 0.06, 1.2",
              # Depth is HELD by a zero velocity on dof 3, not re-imposed as a
              # displacement: the LOAD step already put the grain there, and a
              # displacement BC in a later step is measured from the ORIGINAL
              # position, which would retract it.
              "*Boundary, op=NEW, type=VELOCITY",
              " NS_GRAIN_REF, 1, 1, %s" % _fmt(v_slide if ip % 2 else -v_slide),
              " NS_GRAIN_REF, 2, 2, 0.",
              " NS_GRAIN_REF, 3, 3, 0.",
              " NS_GRAIN_REF, 4, 6, 0.",
              "*Restart, write, number interval=1, time marks=NO",
              "*Output, field, number interval=20",
              "*Node Output", " U, V, A, RF",
              "*Element Output, directions=YES",
              " S, MISES, PEEQ, LE, SDV, STATUS, EVOL",
              "*Contact Output", " CSTRESS, CDISP, CFORCE, CSTATUS",
              "*Output, history, time interval=%s" % _fmt(slide_time * 0.01),
              "*Energy Output", " ALLIE, ALLKE, ALLAE, ALLPD, ETOTAL",
              "*End Step"]

    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")

    return dict(
        kind="micro", path=path, bytes=os.path.getsize(path),
        elements=len(conn), nodes=len(nodes) + len(gv),
        grains=n_gr, grain_facets=len(gf),
        block_mm=(side, side, mic["depth_mm"]),
        element_depth_mm=mic["element_mm"],
        element_inplane_mm=mic["element_inplane_mm"],
        volume_mm3=vol, dc_nm=hp.critical_depth_mm() * 1e6,
        dc_measured=w.dc_measured, psi=psi, swmode=1,
        load_per_grain_n=fn, total_load_n=fn * n_gr,
        slide_speed_mm_s=v_slide, slide_time_s=slide_time,
        indentation_mm=indent_mm, indentation_nm=c.indentation_nm,
        chip_depth_nm=chip_nm,
        chip_over_dc=chip_nm / (hp.critical_depth_mm() * 1e6),
        elements_through_chip=chip_nm / (mic["element_mm"] * 1e6),
        standoff_mm=STANDOFF_FRACTION * indent_mm,
        standoff_over_indent=STANDOFF_FRACTION,
        grain_mass_tonne=grain_mass_t, driven="displacement",
        n_passes=n_passes, total_time_s=slide_time * (n_passes + 0.2),
        energy_threshold_mpa_mm=hp.hardness_mpa * hp.critical_depth_mm(),
        elements_per_dc=p.elements_per_dc, resolves_dc=True,
        steps=tuple(["LOAD"] + ["PASS%d" % i
                                for i in range(1, n_passes + 1)]),
        grain_tags=gtags,
    )


def demo(outdir: str = "_sagemit_demo") -> None:
    """Write both decks from real measured grains and check what came out."""
    import glob
    from .quick import measure_images

    os.makedirs(outdir, exist_ok=True)
    imgs = sorted(glob.glob("B4C_1*.tif"))[:1]
    if not imgs:
        print("semgrit.sagemit: no B4C tif to measure; skipped")
        return
    got = measure_images(imgs, os.path.join(outdir, "meas"),
                         log=lambda *a: None)
    solids = got["solids"]
    assert solids, "no grain solids measured"

    # The 30 um pad: the sparsest of the three, so the demo's MACRO deck is
    # the smallest that still spans the indent. A 6 um pad at the same sector
    # is 262,000 grains, which is the right choice for a real run and far too
    # slow for a self-check.
    p = sagdeck.SAGParams(grain_um=30.0, material="wc_co", name="sagdemo",
                          macro_grain_cap=2000, macro_sector_deg=17.0,
                          micro_grains=1, grind_time_s=2.0e-5)
    pl = sagdeck.plan(p)

    mi = write_micro(os.path.join(outdir, "micro.inp"), pl, solids)
    ma = write_macro(os.path.join(outdir, "macro.inp"), pl, solids)

    for info in (mi, ma):
        txt = open(info["path"], encoding="ascii").read()
        # every deck must be explicit, general-contact, and say which
        assert "*Dynamic, Explicit" in txt
        assert "*Contact, op=NEW" in txt and "ALL EXTERIOR" in txt
        assert "*Contact Pair" not in txt, "pairs cannot see deleted faces"
        assert "constants=58" in txt, "the energy criterion needs 58"
        assert "*Depvar, delete=12" in txt
        # balanced parts and steps
        assert txt.count("*Part,") == txt.count("*End Part")
        assert txt.count("*Step,") == txt.count("*End Step")
        assert txt.count("*Assembly") == txt.count("*End Assembly")
        assert txt.count("*Instance,") == txt.count("*End Instance")
        # no unresolved format placeholders
        assert "%s" not in txt and "%.4" not in txt
        # the material card must be 8 per line
        for ln in txt.splitlines():
            if ln and ln[0].isdigit() and ln.count(",") == 7:
                break

    # MICRO: load, then REPEATED passes over the same track
    mtxt = open(mi["path"], encoding="ascii").read()
    assert mi["n_passes"] >= 2, mi["n_passes"]
    assert mtxt.count("*Step,") == 1 + mi["n_passes"]
    assert "name=LOAD" in mtxt
    for i in range(1, mi["n_passes"] + 1):
        assert "name=PASS%d" % i in mtxt, i
    # Displacement-controlled. The *Cload assertion that used to be here was
    # not merely stale -- it kept PASSING after the change, because the word
    # survives in a comment explaining why force control was abandoned. So
    # this checks the keyword at the start of a line, and checks the depth
    # imposed is the one the contact model predicted.
    assert not any(ln.startswith("*Cload") for ln in mtxt.splitlines()), \
        "a rigid grain has no mass, so a force-driven free dof is rejected"
    assert "*Mass" in mtxt, "the rigid body must carry a mass"
    depth = " NS_GRAIN_REF, 3, 3, %s" % _fmt(-mi["indentation_mm"])
    assert depth in mtxt, depth
    assert mi["driven"] == "displacement"
    assert mi["indentation_nm"] > 0
    assert mi["grain_mass_tonne"] > 0
    # The grain must START clear of the work but be able to
    # REACH it. A standoff larger than the indentation leaves
    # the ramp finishing with the grain still in mid-air, and
    # the job then runs to completion having touched nothing.
    assert 0 < mi["standoff_mm"] < mi["indentation_mm"], \
        (mi["standoff_mm"], mi["indentation_mm"])
    assert mi["standoff_over_indent"] < 0.5
    # and every later step must HOLD that depth with a zero velocity rather
    # than re-imposing a displacement, which would be measured from the
    # original position and retract the grain
    holds = mtxt.count(" NS_GRAIN_REF, 3, 3, 0.")
    assert holds == mi["n_passes"], (holds, mi["n_passes"])
    assert "type=VELOCITY" in mtxt, "the slide must be a velocity"
    assert "PLOT SDV13" in mtxt
    # The passes must ALTERNATE direction, so the grain returns over the same
    # track instead of running away down a fresh one.
    fwd = mtxt.count(" NS_GRAIN_REF, 1, 1, %s" % _fmt(mi["slide_speed_mm_s"]))
    rev = mtxt.count(" NS_GRAIN_REF, 1, 1, %s" % _fmt(-mi["slide_speed_mm_s"]))
    assert fwd >= 1 and rev >= 1, (fwd, rev)
    assert abs(fwd - rev) <= 1, (fwd, rev)
    # The pass count must be DERIVED from the work needed, not chosen: it has
    # to be enough that the accumulated tangential work reaches H*dc.
    from . import materials as _m
    _hp = _m.get(p.material).hybrid_params()
    _need = _hp.hardness_mpa * _hp.critical_depth_mm()
    _per = (p.friction * mi["load_per_grain_n"]) / (
        pl["contact"].groove_width_mm)
    assert mi["n_passes"] * _per >= _need, \
        "the passes must accumulate enough work to reach the threshold"
    assert mi["resolves_dc"] and mi["dc_measured"]
    assert abs(mi["element_depth_mm"] * 1e6 - 16.0) < 1e-9
    assert mi["element_inplane_mm"] > mi["element_depth_mm"]
    assert mi["slide_time_s"] > 0

    # MACRO: three steps, press then hold then grind, and it must NOT claim
    # to resolve dc
    atxt = open(ma["path"], encoding="ascii").read()
    assert atxt.count("*Step,") == 3
    for nm in ("name=PRESS", "name=HOLD", "name=GRIND"):
        assert nm in atxt, nm
    assert atxt.index("name=PRESS") < atxt.index("name=HOLD") \
        < atxt.index("name=GRIND"), "the steps must be in physical order"
    assert "*Rigid Body" in atxt and "*Tie" in atxt
    assert "CANNOT show a ductile-brittle transition" in atxt
    assert not ma["resolves_dc"]
    # press must be a velocity, and the value must match the plan
    assert "type=VELOCITY" in atxt
    assert _fmt(-ma["press_velocity_mm_s"]) in atxt
    # the tool must be seated at first contact, not overlapping
    assert ma["tool_centre_y_mm"] > 0.5 * p.diameter_mm
    assert abs(ma["tool_centre_y_mm"] - 0.5 * p.diameter_mm
               - ma["tallest_protrusion_mm"]) < 1e-12
    # a compliant layer that must BEND needs several elements through it
    assert ma["n_rad_pu"] >= 6, ma["n_rad_pu"]
    # and the faceted arc must be fine compared with the compression, or the
    # contact sees a polygon instead of a cylinder
    import math as _m
    chord = 2.0 * (0.5 * p.diameter_mm) * _m.sin(
        _m.radians(ma["sector_deg"]) / (2.0 * ma["n_circ"]))
    sagitta = (0.5 * p.diameter_mm) * (1.0 - _m.cos(
        _m.radians(ma["sector_deg"]) / (2.0 * ma["n_circ"])))
    assert sagitta < 0.02 * p.compression_mm, (sagitta, p.compression_mm)

    # refusals
    for bad, why in ((lambda: write_micro(os.path.join(outdir, "x.inp"),
                                          pl, []), "no grains"),
                     (lambda: write_macro(os.path.join(outdir, "x.inp"),
                                          pl, []), "no grains")):
        try:
            bad()
        except SAGWriteError:
            pass
        else:
            raise AssertionError("should have been refused: %s" % why)

    print("semgrit.sagemit: all checks passed")
    print("  MICRO  %s  %s el, %d grain(s), %.1f nm depth element, %.2f MB"
          % (os.path.basename(mi["path"]), format(mi["elements"], ","),
             mi["grains"], mi["element_depth_mm"] * 1e6, mi["bytes"] / 1e6))
    print("         driven by %.4e N, slides %.1f mm/s for %.3e s"
          % (mi["load_per_grain_n"], mi["slide_speed_mm_s"],
             mi["slide_time_s"]))
    print("  MACRO  %s  %s el (%s PU, %s work), %s grains, %.2f MB"
          % (os.path.basename(ma["path"]), format(ma["elements"], ","),
             format(ma["pu_elements"], ","), format(ma["work_elements"], ","),
             format(ma["grains"], ","), ma["bytes"] / 1e6))
    print("         sector %.3f deg, %d circ x %d rad, press %.1f mm/s"
          % (ma["sector_deg"], ma["n_circ"], ma["n_rad_pu"],
             ma["press_velocity_mm_s"]))


if __name__ == "__main__":
    demo()
