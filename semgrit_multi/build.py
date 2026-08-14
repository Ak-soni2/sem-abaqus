"""One call from a grain library to a multi-abrasive ductile/brittle deck.

    build_multi(MultiParams(...), solids, outdir)

Three steps, each of which already exists and is already gated:

1. ``semgrit.build_deck`` writes the deck exactly as it does today, with
   ``HybridParams(h_source=1)`` so the material card asks for the chip
   thickness from field variable 1 instead of from four hard-coded constants.
2. ``semgrit_multi.envelope`` sweeps the grit trajectories and returns the
   undeformed chip thickness for every element.
3. ``semgrit_multi.fieldinject`` adds that field to the deck.

Nothing in ``semgrit`` or in ``vumat_grind.for`` is modified by any of it.

The single-grit deck stays available and stays the reference: run this with one
grit and the swept field must reproduce the closed-form wedge that
``verify_hybrid_deck.py`` already validates. ``verify_envelope.py`` checks
exactly that, which is what lets the general engine be trusted at 700 grits.
"""

from __future__ import annotations

import dataclasses
import json
import math
import os
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from semgrit.analysis import AnalysisParams, wheel_motion
from semgrit.build_deck import DeckParams, PRESETS, build_deck, plan_deck
from semgrit.hybrid import HybridParams, kic_from_mpa_sqrt_m
from semgrit.rigid_wheel import place_workpiece

from .envelope import EnvelopeParams, nodal_field, sweep_envelope
from .fieldinject import inject_field


class MultiError(RuntimeError):
    pass


@dataclass
class MultiParams:
    """A multi-abrasive hybrid deck. Everything else is inherited."""

    name: str = "multi_abrasive_hybrid"

    # ---- the wheel and how many abrasives -------------------------------
    diameter_mm: float = 50.0
    arc_length_mm: float = 2.0
    rim_depth_mm: float = 0.012
    width_mm: float = 0.030
    grit_mode: str = "count"
    """'single', 'count', 'areal_density' or 'concentration'. The whole point of
    this package is that anything other than 'single' now works."""
    grit_count: int = 12
    areal_density_per_mm2: float = 5000.0
    concentration: float = 100.0
    grit_arc_window_mm: float = 0.0
    grit_width_window_mm: float = 0.0
    seed: int = 20260731

    # ---- the workpiece ---------------------------------------------------
    wp_length_mm: float = 0.048
    wp_width_mm: float = 0.015
    wp_depth_mm: float = 0.006
    element_um: float = 0.30
    element_axial_um: float = 0.0
    element_depth_um: float = 0.0
    # Graded depth mesh, straight through to DeckParams. Resolving a
    # nanometre-scale cut needs a very fine surface layer, and grading is what
    # makes that affordable: dt follows the smallest element, so the fine layer
    # costs increments it was going to cost anyway, while the coarse body below
    # costs almost no elements.
    surface_layer_um: float = 0.0
    depth_growth: float = 1.3
    max_depth_element_um: float = 0.0
    protrusion_std: float = 0.12
    """Spread of grit protrusion. Small values mean a well-dressed wheel where
    many grits stand at nearly the same height -- which is what lets SEVERAL of
    them cut at a shallow depth of cut instead of only the tallest."""
    standoff_um: float = 0.0
    depth_of_cut_um: float = 0.0
    wp_position: str = "centred"
    """'centred' matches the validated single-grit deck, which is what makes the
    single-grit regression against the closed form apples-to-apples. 'under the
    tallest grit' centres the block on one grain, which halves the arc that
    grain sweeps across it."""
    surface_speed_m_s: float = 30.0
    cores: int = 8

    # ---- the material ----------------------------------------------------
    material: str = "sandstone"
    """A key from :data:`semgrit.materials.MATERIALS`. Sets the JH-2 card, its
    density and the ``*Material`` name together. If ``hybrid`` is left None the
    ductile constants come from the same entry, so the two branches cannot end
    up describing two different materials."""
    hybrid: Optional[HybridParams] = None
    envelope: Optional[EnvelopeParams] = None

    # ---- output ----------------------------------------------------------
    keep_geometry_deck: bool = True
    """Keep the un-injected deck as well. It is the same model without the
    field, so a diff of the two is exactly the field and nothing else."""

    def deck_params(self) -> DeckParams:
        from semgrit import materials
        mat = materials.get(self.material)
        hp = self.hybrid or mat.hybrid_params()
        if hp.h_source != 1:
            hp = dataclasses.replace(hp, h_source=1)
        from semgrit.hybrid import HYBRID_DEPVAR
        an = AnalysisParams(
            enabled=True, material_model="hybrid", hybrid=hp,
            jh2_constants=mat.jh2, jh2_density_kg_m3=mat.density_kg_m3,
            n_depvar=HYBRID_DEPVAR, element_deletion=True,
            depth_of_cut_um=self.depth_of_cut_um)
        return DeckParams(
            wp_material=mat.inp_material,
            name=self.name, sector_mode="arc",
            arc_length_mm=self.arc_length_mm, diameter_mm=self.diameter_mm,
            rim_depth_mm=self.rim_depth_mm, width_mm=self.width_mm,
            grit_mode=self.grit_mode, grit_count=self.grit_count,
            areal_density_per_mm2=self.areal_density_per_mm2,
            concentration=self.concentration,
            grit_arc_window_mm=self.grit_arc_window_mm,
            grit_width_window_mm=self.grit_width_window_mm,
            seed=self.seed,
            wp_length_mm=self.wp_length_mm, wp_width_mm=self.wp_width_mm,
            wp_depth_mm=self.wp_depth_mm,
            wp_element_size_mm=self.element_um / 1000.0,
            wp_element_size_width_mm=self.element_axial_um / 1000.0,
            wp_element_size_depth_mm=self.element_depth_um / 1000.0,
            wp_surface_layer_mm=self.surface_layer_um / 1000.0,
            wp_depth_growth=self.depth_growth,
            wp_max_depth_element_mm=self.max_depth_element_um / 1000.0,
            protrusion_std=self.protrusion_std,
            clearance_um=self.standoff_um, wp_position=self.wp_position,
            surface_speed_mm_s=self.surface_speed_m_s * 1000.0,
            cores=self.cores, analysis=an, also_write_cae_deck=False)


def plan_multi(p: MultiParams, solids: Sequence, *, paths=None,
               log=print) -> dict:
    """What the deck would contain, and how the switch would split it.

    Writes nothing. Runs the same sweep the build will, so the split reported
    here is the split the deck gets.
    """
    dp = p.deck_params()
    plan = plan_deck(dp, solids)
    hp = dp.analysis.hybrid
    dc = hp.critical_depth_mm()
    step_time = float((plan.get("cost") or {}).get("step_time_s") or 0.0)
    an = dataclasses.replace(dp.analysis,
                             depth_of_cut_um=float(plan["depth_of_cut_um"]))
    motion = wheel_motion(an, plan["_place"]["theta_c"], dp.surface_speed_mm_s,
                          dp.outer_radius_mm, step_time)
    env = sweep_envelope(plan["_place"], motion, plan["_wp"],
                         step_time_s=step_time,
                         rotation_reversed=bool(an.rotation_reversed),
                         params=p.envelope, paths=paths, log=log)
    return {"plan": plan, "envelope": env, "motion": motion, "dc_mm": dc,
            "split": env.split(dc), "deck_params": dp, "step_time_s": step_time}


def build_multi(p: MultiParams, solids: Sequence, outdir: str,
                *, paths=None, log=print) -> dict:
    """Build, sweep, inject. Returns a summary including the field statistics.

    ``paths`` is forwarded to :func:`semgrit_multi.envelope.sweep_envelope`:
    ``{grit_index: (n, 4) array}`` replaces that grit's ideal arc with a
    measured one. See :mod:`semgrit_multi.trajectory`.
    """
    os.makedirs(outdir, exist_ok=True)
    dp = p.deck_params()
    hp = dp.analysis.hybrid
    dc = hp.critical_depth_mm()

    log("1/3  writing the deck (%s grits) ..." % (
        p.grit_count if p.grit_mode == "count" else p.grit_mode))
    info = build_deck(dp, solids, outdir)
    log("     %s: %.1f MB, %s grits, %s C3D8R"
        % (os.path.basename(info["path"]), info["size_bytes"] / 1e6,
           format(info["n_grits"], ","),
           format(info["n_workpiece_elements"], ",")))

    log("2/3  sweeping the grit envelope ...")
    # Rebuild the placement the writer used. place_workpiece is deterministic
    # and takes only the model and the block, so this is the same seating, not
    # a re-derivation of it.
    wp = _workpiece_of(dp)
    model = _model_of(dp, solids)
    place = place_workpiece(model, wp, dp.clearance_um, dp.wp_position,
                            dp.wp_position_deg,
                            not bool(dp.analysis.rotation_reversed))
    step_time = float(info["cost"]["step_time_s"])
    motion = info["motion"]
    if motion is None:
        raise MultiError("the deck carries no motion, so nothing sweeps")
    env = sweep_envelope(place, motion, wp, step_time_s=step_time,
                         rotation_reversed=bool(dp.analysis.rotation_reversed),
                         params=p.envelope, paths=paths, log=log)
    split = env.split(dc)
    log("     dc = %.4f nm -> of %s elements CUT, %s ductile (%.1f%%), %s brittle"
        % (dc * 1e6, format(split["n_cut"], ","),
           format(split["n_ductile_of_cut"], ","),
           100.0 * split["ductile_fraction_of_cut"],
           format(split["n_brittle_of_cut"], ",")))

    log("3/3  injecting the field ...")
    nodal = nodal_field(env, wp)
    plain = info["path"]
    final = os.path.join(outdir, p.name + "_field.inp")
    inj = inject_field(
        plain, final, nodal,
        comment=["swept from %d of %d grits, %s surface points"
                 % (env.n_grits_engaged, len(place["frames"]),
                    format(env.stats["n_surface_points_swept"], ",")),
                 "dc = %.6e mm; of %d elements cut, %d ductile and %d brittle"
                 % (dc, split["n_cut"], split["n_ductile_of_cut"],
                    split["n_brittle_of_cut"]),
                 "deepest groove %.4f um, mean %.4f um"
                 % (env.stats["max_depth_removed_um"],
                    env.stats["mean_depth_removed_um"])])
    log("     %s: %.1f MB (+%.1f MB of field), %s values"
        % (os.path.basename(final), inj["size_bytes"] / 1e6,
           (inj["size_bytes"] - info["size_bytes"]) / 1e6,
           format(inj["n_values"], ",")))
    if not p.keep_geometry_deck:
        os.remove(plain)

    out = dict(info)
    out.update(path=final, plain_path=(plain if p.keep_geometry_deck else None),
               size_bytes=inj["size_bytes"],
               envelope=dict(env.stats), split=split, injected=inj,
               dc_mm=dc, dc_nm=dc * 1e6,
               grit_order=env.grit_order,
               per_grit_h_um={str(k): [a * 1000.0, b * 1000.0]
                              for k, (a, b) in env.per_grit_h.items()})
    with open(os.path.join(outdir, p.name + "_multi_report.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    # verify_hybrid_deck.py looks for "<deck>_report.json" beside the deck it is
    # given, and the injected deck is "<name>_field.inp". Writing the report
    # under that name too lets the existing gate check the injected deck without
    # the gate having to learn anything about this package.
    with open(os.path.join(outdir, p.name + "_field_report.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    # The element field itself, for plotting and for the gate.
    np.save(os.path.join(outdir, p.name + "_h_elem.npy"), env.h_elem)
    np.save(os.path.join(outdir, p.name + "_depth_removed.npy"),
            env.depth_removed)
    return out


def _workpiece_of(dp: DeckParams):
    from semgrit.build_deck import workpiece_of
    return workpiece_of(dp)


def _model_of(dp: DeckParams, solids: Sequence):
    """The placed wheel model, by the writer's own code path.

    ``build_deck(..., dry_run=True)`` stops after deciding what the model is and
    before writing anything, and returns that model. Using it here means the
    sweep cannot be looking at a different grit layout from the deck.
    """
    _info, model = build_deck(dp, solids, None, return_model=True,
                              dry_run=True)
    return model


def summary_text(res: dict, wp_length_mm: float) -> str:
    """The split, as a block of text for a notebook cell."""
    s = res["split"]
    e = res["envelope"]
    L: list[str] = []
    a = L.append
    a("SWEPT CHIP-THICKNESS FIELD")
    a("  grits crossing the block   : %d of %d"
      % (e["n_grits_engaged"], e["n_grits_total"]))
    a("  surface points swept       : %s per grit crossing"
      % format(e["n_surface_points_swept"], ","))
    a("  station resolution         : %.4f um" % e["station_step_um"])
    a("  time samples per grit      : %s to %s, depth blur %.2f nm"
      % (format(e["time_samples_min"], ","), format(e["time_samples_max"], ","),
         e["depth_resolution_nm"]))
    a("  chip thickness where cut   : %.4f to %.4f nm, mean %.4f nm"
      % (e["h_cut_min_um"] * 1000.0, e["h_cut_max_um"] * 1000.0,
         e["h_cut_mean_um"] * 1000.0))
    a("  groove depth               : max %.4f um, mean %.4f um"
      % (e["max_depth_removed_um"], e["mean_depth_removed_um"]))
    a("  material removed           : %.6e mm3" % e["removed_volume_mm3"])
    a("  elements through the cut   : %.2f  (surface layer %.4f um)"
      % (e["elements_through_deepest_cut"], e["surface_layer_um"]))
    if e["elements_through_deepest_cut"] < 1.0:
        a("    the deepest cut is thinner than one element: the mesh cannot")
        a("    hold this chip. Refine the depth direction or cut deeper.")
    elif e["elements_through_deepest_cut"] < 3.0:
        a("    under 3 elements through the cut: the chip is barely resolved")
        a("    and the force will read low.")
    a("")
    a("THE SWITCH")
    a("  dc                         : %.4f nm" % res["dc_nm"])
    a("  elements the grits cut     : %s of %s"
      % (format(s["n_cut"], ","), format(s["n_elements"], ",")))
    a("    of those, ductile        : %s  (%.1f%%)"
      % (format(s["n_ductile_of_cut"], ","),
         100.0 * s["ductile_fraction_of_cut"]))
    a("    of those, brittle        : %s" % format(s["n_brittle_of_cut"], ","))
    a("  never cut (subsurface etc) : %s, given the nearest cut station's h"
      % format(s["n_never_cut"], ","))
    a("  law each element will run  : %s ductile, %s brittle"
      % (format(s["n_ductile_law"], ","), format(s["n_brittle_law"], ",")))
    if s["n_cut"] == 0:
        a("  NOTHING IS CUT. No grit removes material from any element, so the")
        a("  switch has nothing to act on. Either the infeed never reaches the")
        a("  work, or the cut is far thinner than one element.")
    elif s["n_ductile_of_cut"] == 0:
        a("  every cut element is BRITTLE: no cut is thinner than dc.")
        a("  Cut shallower, raise dc, or use the Bifano form of dc.")
    elif s["n_brittle_of_cut"] == 0:
        a("  every cut element is DUCTILE: nothing will fracture.")
        a("  Cut deeper or lower dc.")
    else:
        a("  Both regimes present. Plot SDV13 to see where they divide.")
    return "\n".join(L)
