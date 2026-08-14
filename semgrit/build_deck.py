"""One entry point from a measured grain library to a verified Abaqus deck.

``build_deck(DeckParams(...), solids, outdir)`` writes the ``.inp``, the CAE loader
script, the placement table, the report JSON and optionally STEP/STL CAD, and returns
a summary dict. Everything a user should be able to choose lives in
:class:`DeckParams`; nothing else in the pipeline needs editing.

This exists because the knowledge of *how* to assemble a correct deck -- inset the
grit band so grains do not hang off the bond, seat the workpiece tangent to the
tallest grit that can actually reach it, keep the whole wheel one rigid body -- was
spread across two hand-written build scripts. A notebook cannot ask a user to edit
those, so the logic is collected here and the scripts become presets.

Wheel extent can be given three ways, whichever the user thinks in:
``sector_mode='angle'`` (degrees), ``'arc'`` (arc length in mm) or ``'full'`` (360).

Grit population can be given four ways: ``grit_mode='concentration'`` (C-number),
``'areal_density'`` (grains/mm^2), ``'count'`` (exactly N grains) or ``'single'``
(one chosen grain, for a single-grit scratch test).
"""

from __future__ import annotations

import csv
import dataclasses
import json
import math
import os
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from .abaqus import write_cae_import_script
from .analysis import AnalysisParams
from .rigid_wheel import write_rigid_wheel_inp
from .wheel import (UM_PER_MM, GrainPopulationSpec, WheelModel, WheelSpec,
                    build_rim_mesh, build_wheel, check_grain_overlaps)
from .wheel_workpiece import WorkpieceBlock, rotate_placements_about_z

# Throughput used for the wall-clock estimate: element-increments per second per
# core for C3D8R under a scalar-loop VUMAT. Deliberately pessimistic -- being wrong
# in this direction only means the job finishes early.
RATE_PER_CORE = 3.0e5
PARALLEL_EFFICIENCY = 0.70
CONTACT_OVERHEAD = 1.6


class DeckError(RuntimeError):
    pass


@dataclass
class DeckParams:
    """Every user-facing knob. Millimetres, tonnes, seconds, MPa, N throughout."""

    # ---- wheel ----------------------------------------------------------
    diameter_mm: float = 50.0
    sector_mode: str = "arc"              # 'angle' | 'arc' | 'full'
    sector_deg: float = 30.0              # used when sector_mode='angle'
    arc_length_mm: float = 2.0            # used when sector_mode='arc'
    rim_depth_mm: float = 0.012
    width_mm: float = 0.030
    shell_circumferential_divisions: int = 200
    shell_axial_divisions: int = 6
    shell_radial_divisions: int = 1
    bond_density_kg_m3: float = 2700.0

    # ---- grits ----------------------------------------------------------
    grit_mode: str = "concentration"      # 'concentration'|'areal_density'|'count'|'single'
    concentration: float = 100.0
    areal_density_per_mm2: float = 5000.0
    grit_count: int = 500
    single_grain_index: int = -1           # -1 = largest grain in the library
    single_grit_offset_mm: float = 0.015   # tangential offset from the block centre
    grit_arc_window_mm: float = 0.0
    """Confine grits to this arc length, centred mid-arc. 0 = the whole sector.

    A big sector at a realistic density needs an unusable number of grains -- 13 mm of
    arc at 3500/mm^2 is 22,000 of them and 300 MB. Only the arc the workpiece sweeps
    can ever touch it, so dressing a window and leaving the rest of the rim bare costs
    nothing physically and keeps the deck openable."""
    grit_width_window_mm: float = 0.0
    """Dress only this much of the wheel's face, centred. 0 = the full width.

    Lets the slice be as thick as it needs to look like a real chunk of wheel while the
    grains stay in the band the workpiece actually runs in. A 3 mm face dressed at
    5000/mm^2 would be 30,000 grains for no benefit -- the workpiece is 0.25 mm wide."""
    inset_grit_band: bool = True
    protrusion_mean: float = 0.55
    protrusion_std: float = 0.12
    protrusion_min: float = 0.25
    protrusion_max: float = 0.85
    max_tilt_deg: float = 35.0
    spacing_factor: float = 1.05
    seed: int = 20260731

    # ---- workpiece ------------------------------------------------------
    include_workpiece: bool = True
    wp_length_mm: float = 0.048
    wp_width_mm: float = 0.015
    wp_depth_mm: float = 0.006
    wp_element_size_mm: float = 0.0003
    # Per-direction mesh size. 0 means "use wp_element_size_mm for this direction".
    # The element type is fixed at C3D8R; only the size is a user choice.
    wp_element_size_length_mm: float = 0.0
    wp_element_size_width_mm: float = 0.0
    wp_element_size_depth_mm: float = 0.0
    # Graded depth mesh: fine at the ground face, coarsening into the body. With
    # wp_surface_layer_mm = 0 the mesh is uniform, as before.
    wp_surface_layer_mm: float = 0.0
    wp_depth_growth: float = 1.3
    wp_max_depth_element_mm: float = 0.0
    wp_material: str = "STONE"
    wp_density_kg_m3: float = 2650.0
    wp_youngs_modulus_mpa: float = 50_000.0
    wp_poisson_ratio: float = 0.25
    clearance_um: float = 0.0
    """Radial standoff, in microns, between the block's ground face and the tallest
    grain underneath it. 0 = touching, with zero initial overclosure. Choose it against
    the grain heights the plan reports -- a standoff larger than the depth of cut means
    the infeed never closes the gap and nothing is ground."""
    wp_position: str = "centred"
    """Where along the arc the block sits: 'centred', 'first grit at entry',
    'under the tallest grit' or 'custom angle'. See rigid_wheel.WP_POSITIONS."""
    wp_position_deg: float = 0.0

    # ---- kinematics used for the run estimate ---------------------------
    surface_speed_mm_s: float = 30_000.0
    travel_mm: float = 0.0                # 0 = workpiece length + travel_margin_mm
    travel_margin_mm: float = 0.006
    cores: int = 8

    # ---- run-ready analysis (step, BCs, contact, material, output) -------
    analysis: Optional[AnalysisParams] = None
    """None or disabled = geometry only, to be finished in CAE. Enabled = a deck you
    can submit from the terminal with no GUI at all."""

    # ---- output ---------------------------------------------------------
    name: str = "wheel"
    also_write_cae_deck: bool = False
    """With a run-ready analysis, also emit the geometry-only twin plus its CAE loader,
    so the same grit layout can be opened and assembled by hand as well as submitted
    from the terminal. Written from the same placed model, so the two decks describe
    exactly the same wheel."""
    write_step: bool = False              # assembled wheel, for SOLIDWORKS
    write_stl: bool = False
    step_max_grains: int = 0              # 0 = every grain
    stl_max_grains: int = 0
    write_grain_stls: bool = False        # one STL per measured grain
    write_grains_step: bool = False       # all grains laid out on a grid, one STEP
    grains_step_max: int = 200

    # ---------------------------------------------------------------------
    @property
    def outer_radius_mm(self) -> float:
        return self.diameter_mm / 2.0

    def resolved_sector_deg(self) -> float:
        if self.sector_mode == "full":
            return 360.0
        if self.sector_mode == "angle":
            return float(self.sector_deg)
        if self.sector_mode == "arc":
            deg = math.degrees(self.arc_length_mm / self.outer_radius_mm)
            if deg > 360.0:
                raise DeckError(
                    "arc_length_mm=%g exceeds the circumference of a %g mm wheel "
                    "(%.3f mm); use sector_mode='full'"
                    % (self.arc_length_mm, self.diameter_mm,
                       math.pi * self.diameter_mm))
            return deg
        raise DeckError("sector_mode must be 'angle', 'arc' or 'full', not %r"
                        % self.sector_mode)

    def validate(self) -> None:
        if self.diameter_mm <= 0:
            raise DeckError("diameter_mm must be positive")
        if not 0 < self.rim_depth_mm < self.outer_radius_mm:
            raise DeckError("rim_depth_mm must be in (0, radius)")
        if self.width_mm <= 0:
            raise DeckError("width_mm must be positive")
        if self.grit_mode not in ("concentration", "areal_density", "count", "single"):
            raise DeckError("grit_mode must be concentration|areal_density|count|single")
        if self.grit_mode == "count" and self.grit_count < 1:
            raise DeckError("grit_count must be at least 1")
        if self.include_workpiece:
            for nm in ("wp_length_mm", "wp_width_mm", "wp_depth_mm",
                       "wp_element_size_mm"):
                if getattr(self, nm) <= 0:
                    raise DeckError("%s must be positive" % nm)
            h = self.wp_element_size_mm
            for nm in ("wp_length_mm", "wp_width_mm", "wp_depth_mm"):
                if getattr(self, nm) < h:
                    raise DeckError("%s is smaller than one element (%g mm)" % (nm, h))
            # Element ASPECT RATIO. C3D8R with enhanced hourglass control is
            # tolerant but not indifferent, and the shipped surface brick is
            # 0.3 x 1.5 x 0.03 um -- 50:1 axially. That is the reason the groove
            # comes out 1.7 elements wide, which is what makes lateral pile-up
            # and lateral cracking, the two signatures that separate the
            # regimes, unrepresentable. Warn rather than raise: a coarse axial
            # mesh is a legitimate choice when nothing axial is being measured.
            sizes = [s for s in (self.wp_element_size_mm,
                                 self.wp_element_size_width_mm or
                                 self.wp_element_size_mm,
                                 self.wp_element_size_depth_mm or
                                 self.wp_element_size_mm) if s > 0]
            if sizes:
                ar = max(sizes) / min(sizes)
                if ar > 10.0:
                    object.__setattr__(self, "_aspect_warning", (
                        "workpiece element aspect ratio is %.0f:1 (%g / %g mm). "
                        "Above about 10:1 a C3D8R groove is only a couple of "
                        "elements wide and lateral flow cannot be resolved."
                        % (ar, max(sizes), min(sizes))))
        self.resolved_sector_deg()


def workpiece_of(p: DeckParams) -> Optional[WorkpieceBlock]:
    if not p.include_workpiece:
        return None
    return WorkpieceBlock(
        length_mm=p.wp_length_mm, width_mm=p.wp_width_mm, depth_mm=p.wp_depth_mm,
        element_size_mm=p.wp_element_size_mm, material=p.wp_material,
        density_kg_m3=p.wp_density_kg_m3,
        youngs_modulus_mpa=p.wp_youngs_modulus_mpa, poisson_ratio=p.wp_poisson_ratio,
        element_size_length_mm=p.wp_element_size_length_mm or None,
        element_size_width_mm=p.wp_element_size_width_mm or None,
        element_size_depth_mm=p.wp_element_size_depth_mm or None,
        surface_layer_mm=p.wp_surface_layer_mm,
        depth_growth=p.wp_depth_growth,
        max_depth_element_mm=p.wp_max_depth_element_mm,
    )


def _bond_spec(p: DeckParams, sector_deg: float) -> WheelSpec:
    return WheelSpec(
        diameter_mm=p.diameter_mm, width_mm=p.width_mm, sector_deg=sector_deg,
        rim_depth_mm=p.rim_depth_mm,
        radial_divisions=max(int(p.shell_radial_divisions), 1),
        axial_divisions=max(int(p.shell_axial_divisions), 1),
        circumferential_divisions_per_deg=(
            max(int(p.shell_circumferential_divisions), 3) / sector_deg),
    )


def _empty_model(spec: WheelSpec, pop: GrainPopulationSpec) -> WheelModel:
    """Bond mesh only, no grits yet.

    The rim hex mesh is built even though the Abaqus deck emits a rigid *surface*:
    ``write_wheel_step`` needs the solid to hand a closed bond body to CAD.
    """
    nodes, hexes, sets = build_rim_mesh(spec, with_element_sets=True)
    return WheelModel(spec=spec, population=pop, body_nodes=nodes, body_hexes=hexes,
                      placements=[], shapes=[], requested_grains=0,
                      achieved_grains=0, node_sets=sets, stats={}, warnings=[])


def _population(p: DeckParams, solids: Sequence, sector_deg: float,
                inset_mm: float, band_area_mm2: float) -> GrainPopulationSpec:
    common = dict(
        protrusion_mean=p.protrusion_mean, protrusion_std=p.protrusion_std,
        protrusion_min=p.protrusion_min, protrusion_max=p.protrusion_max,
        max_tilt_deg=p.max_tilt_deg, spacing_factor=p.spacing_factor, seed=p.seed,
    )
    if p.grit_mode == "concentration":
        return GrainPopulationSpec(concentration=p.concentration, **common)
    if p.grit_mode == "areal_density":
        return GrainPopulationSpec(areal_density_per_mm2=p.areal_density_per_mm2,
                                   **common)
    # 'count': turn the requested number into the density that produces it.
    return GrainPopulationSpec(
        areal_density_per_mm2=p.grit_count / max(band_area_mm2, 1e-12), **common)


def cost_model(p: DeckParams, wp, travel: float, step_time: float,
               an, n_wp_elements: int) -> dict:
    """Stable increment, increment count and wall clock, from the material and
    mesh the deck will actually run with. Shared by the writer and the planner so
    a preview cannot quote a different number from the build."""
    if wp is None:
        return {}
    R = p.outer_radius_mm
    if an is not None and an.enabled and an.material_model in ("jh2", "hybrid"):
        # The stable increment has to come from the material the deck will really
        # run with. JH-2 props 1 and 2 are the bulk and shear moduli, so the
        # dilatational speed is sqrt((K + 4G/3)/rho). Using the placeholder
        # elasticity instead overstates the wave speed by 2.7x, and so overstates
        # the increment count and the run time by the same factor.
        K, G = an.jh2_constants[0], an.jh2_constants[1]
        rho = an.jh2_density_kg_m3 * 1e-12
        c = math.sqrt((K + 4.0 * G / 3.0) / rho)
        if an.material_model == "hybrid":
            # The hybrid card carries TWO elasticities, one per branch, and any
            # element may be running either. dt has to satisfy the stiffer of
            # them, so take the faster wave. With the placeholder constants the
            # two coincide by construction -- E = 6500, nu = 0.21 is exactly
            # K1 = 3735.6, G = 2686 -- but a calibrated JC set need not.
            rho = an.hybrid.density_kg_m3 * 1e-12
            nu_j, e_j = an.hybrid.poisson, an.hybrid.youngs_mpa
            c_j = math.sqrt(e_j * (1 - nu_j)
                            / ((1 + nu_j) * (1 - 2 * nu_j) * rho))
            c = max(math.sqrt((K + 4.0 * G / 3.0) / rho), c_j)
    else:
        nu, E = p.wp_poisson_ratio, p.wp_youngs_modulus_mpa
        rho = p.wp_density_kg_m3 * 1e-12
        c = math.sqrt(E * (1 - nu) / ((1 + nu) * (1 - 2 * nu) * rho))
    # The smallest *achieved* dimension, not the requested one: rounding the block
    # to a whole number of elements can shrink it, and with per-direction sizes the
    # governing direction need not be the one the user thought of.
    hl, hw, hd = wp.element_sizes()
    h_min = wp.min_element_size()
    scale = (an.mass_scaling_factor if an is not None and an.enabled else 1.0)
    dt = h_min / c * math.sqrt(scale)
    t_step = step_time or (travel / p.surface_speed_mm_s)
    inc = t_step / dt
    el_inc = n_wp_elements * inc
    per_s = RATE_PER_CORE * max(p.cores, 1) * PARALLEL_EFFICIENCY / CONTACT_OVERHEAD
    cost = dict(dilatational_speed_mm_s=c, stable_dt_s=dt, step_time_s=t_step,
                mass_scaling_factor=scale,
                element_size_cutting_mm=hl, element_size_axial_mm=hw,
                element_size_depth_mm=hd, governing_element_size_mm=h_min,
                depth_layer_min_mm=wp.depth_layer_range()[0],
                depth_layer_max_mm=wp.depth_layer_range()[1],
                element_divisions=list(wp.divisions()),
                increments=inc, element_increments=el_inc, travel_mm=travel,
                surface_speed_mm_s=p.surface_speed_mm_s,
                omega_rad_s=p.surface_speed_mm_s / R,
                rpm=p.surface_speed_mm_s / R * 30.0 / math.pi,
                est_hours={str(n): el_inc / (RATE_PER_CORE * n
                                             * PARALLEL_EFFICIENCY
                                             / CONTACT_OVERHEAD) / 3600.0
                           for n in (1, 2, 4, 8, 16, 32)},
                est_hours_at_requested_cores=el_inc / per_s / 3600.0)
    return cost


def _auto_depth(p, clearance_um) -> float:
    """The depth of cut the build will use: as asked, or automatic when 0.

    Automatic means close the standoff first, then cut 85% of the way through the
    grain protrusion. Ignoring the standoff would pick a depth that never reaches the
    work as soon as the block is parked clear of the grits.
    """
    ae = p.analysis.depth_of_cut_um
    if ae > 0 or not clearance_um:
        return ae
    return round(p.clearance_um + 0.85 * clearance_um, 4)


def plan_deck(p: DeckParams, solids: Sequence) -> dict:
    """Everything the deck *would* contain, without writing a byte.

    Uses the same placement code as the writer -- ``bake_grit`` for the grit vertices
    and ``ground_radius`` for where the block sits -- so a preview drawn from this
    cannot disagree with the deck that follows.
    """
    import math as _m

    from .rigid_wheel import place_workpiece

    info, model = build_deck(p, solids, None, return_model=True, dry_run=True)
    wp = workpiece_of(p)
    R = p.outer_radius_mm
    sector = p.resolved_sector_deg()
    thc = 0.0 if abs(sector - 360.0) < 1e-9 else _m.radians(sector) / 2.0

    place = place_workpiece(
        model, wp, p.clearance_um, p.wp_position, p.wp_position_deg,
        not bool(getattr(p.analysis, "rotation_reversed", False)))
    baked, frames = place["baked"], place["frames"]
    thc = place["theta_c"]
    r_ground = place["r_ground"]
    clr_um = place["clearance_um"] or 0.0

    swept = []
    if wp is not None:
        travel = info.get("cost", {}).get("travel_mm") or 0.0
        # The same band the writer calls "able to engage": half the block plus the
        # distance it travels. Using a different one here would let the preview quote
        # an engagement the deck does not have.
        hb = wp.length_mm / 2.0 + travel
        hz = wp.width_mm / 2.0
        for f in frames:
            sel = (np.abs(f[:, 1]) <= hb) & (np.abs(f[:, 2]) <= hz)
            if sel.any():
                swept.append((r_ground - f[sel, 0].max()) * 1000.0)

    el = (0.0, 0.0, 0.0, 0.0)
    if wp is not None:
        hl, hw, _ = wp.element_sizes()
        lo, hi = wp.depth_layer_range()
        el = (hl * 1000, hw * 1000, lo * 1000, hi * 1000)

    # How tall the abrasive actually stands. You cannot choose a sensible standoff
    # without it: a gap wider than the tallest grain guarantees nothing ever touches,
    # and a depth of cut deeper than it drives the bond into the work.
    def _stats(a):
        a = np.asarray(list(a), dtype=float)
        if not a.size:
            return dict(n=0, min=0.0, max=0.0, mean=0.0, median=0.0)
        return dict(n=int(a.size), min=float(a.min()), max=float(a.max()),
                    mean=float(a.mean()), median=float(np.median(a)))

    prot = np.asarray(place["protrusion_um"], dtype=float)
    under = []
    if wp is not None:
        hb, hz = wp.length_mm / 2.0, wp.width_mm / 2.0
        under = [i for i, f in enumerate(frames)
                 if bool(((np.abs(f[:, 1]) <= hb) & (np.abs(f[:, 2]) <= hz)).any())]
    # Angular span of the block and of the dressed grit, so it is clear which grains
    # the pass will actually meet.
    # The deck's own ES_GRITS_ENGAGE, from the shared helper the writer uses.
    from .rigid_wheel import engaging_grits
    _engage_idx, _engage_win = engaging_grits(
        frames, wp, (info.get("cost", {}) or {}).get("travel_mm") or 0.0)
    th = [float(np.arctan2(v[:, 1], v[:, 0]).mean()) for v in baked]
    # Placement keys off the grains the block can reach across the face, so report
    # that span too -- otherwise the printed grit range can sit outside the block's
    # entry edge and read like a contradiction.
    reach_i = place.get("reachable") or list(range(len(baked)))
    th_reach = [float(np.arctan2(baked[i][:, 1], baked[i][:, 0]).max())
                for i in reach_i]
    half_ang = _m.atan2(wp.length_mm / 2.0, r_ground) if wp is not None else 0.0
    _entry_way = -1.0 if getattr(p.analysis, "rotation_reversed", False) else 1.0

    return dict(
        outer_radius_mm=R, sector_deg=sector, full_wheel=abs(sector - 360.0) < 1e-9,
        rim_depth_mm=p.rim_depth_mm, width_mm=p.width_mm,
        arc_length_mm=R * _m.radians(sector),
        n_grits=len(model.placements), frames=frames,
        areal_density=info.get("achieved_areal_density_per_mm2") or 0.0,
        grit_band_arc_mm=info.get("grit_band_arc_mm") or 0.0,
        grit_band_width_mm=info.get("grit_band_width_mm") or 0.0,
        theta_workpiece_deg=_m.degrees(thc),
        ground_radius_mm=r_ground, bond_clearance_um=clr_um or 0.0,
        swept_clearances_um=swept,
        # The window a depth of cut has to land in: deeper than the nearest grain in
        # the swept band or nothing is touched, shallower than the face-to-bond gap or
        # the rim itself hits the work.
        first_contact_um=(min(swept) if swept else None),
        depth_ceiling_um=(clr_um or 0.0) + p.clearance_um,
        # where the block sits, and what it will meet there
        wp_position=p.wp_position, wp_position_deg=p.wp_position_deg,
        wp_relocated=place["relocated"],
        wp_requested_theta_deg=place["requested_theta_deg"],
        # The entry edge is the one the grains reach first, so it follows the rotation
        # sense. Reporting it as the high-theta edge regardless would contradict a
        # reversed wheel -- and the whole point of naming an "entry" is that it tells
        # you where the pass begins.
        wp_entry_theta_deg=_m.degrees(thc + _entry_way * half_ang),
        wp_exit_theta_deg=_m.degrees(thc - _entry_way * half_ang),
        rotation_reversed=bool(getattr(p.analysis, "rotation_reversed", False)),
        grit_theta_range_deg=((_m.degrees(min(th)), _m.degrees(max(th)))
                              if th else (0.0, 0.0)),
        grit_theta_reachable_deg=((_m.degrees(min(th_reach)),
                                   _m.degrees(max(th_reach)))
                                  if th_reach else (0.0, 0.0)),
        n_grits_reachable=len(reach_i),
        n_grits_under_block=len(under),
        # abrasive heights above the bond, the numbers a standoff is chosen against
        protrusion_um=_stats(prot),
        protrusion_under_block_um=_stats(prot[under] if len(under) else []),
        grain_height_um=_stats(model.shapes[q.shape_index].height_um
                               for q in model.placements),
        grain_width_um=_stats(max(model.shapes[q.shape_index].extent_um()[:2])
                              for q in model.placements),
        standoff_um=p.clearance_um,
        sweep_mm=info.get("cost", {}).get("travel_mm") or 0.0,
        depth_of_cut_um=(_auto_depth(p, clr_um)
                         if p.analysis is not None and p.analysis.enabled else 0.0),
        workpiece=(dict(length_mm=wp.length_mm, width_mm=wp.width_mm,
                        depth_mm=wp.depth_mm) if wp is not None else None),
        n_workpiece_elements=(wp.n_elements() if wp is not None else 0),
        element_um=el, cost=info.get("cost"),
        estimated_mb=info["estimated_mb"],
        warnings=info.get("warnings", []), notes=info.get("notes", []),
        # handed to the 3-D viewer so it draws the very same geometry, the very same
        # contact surface, and the very same boundary conditions
        _model=model, _wp=wp, _place=place, _params=p,
        _engage=_engage_idx, engage_window_mm=_engage_win,
        title="%s  -  %s" % (p.name, "run-ready" if (p.analysis is not None
                                                     and p.analysis.enabled)
                             else "geometry only"),
    )


def build_deck(p: DeckParams, solids: Sequence, outdir: str,
               return_model: bool = False, dry_run: bool = False):
    """Build and write one deck. Returns the report dict (also saved as JSON).

    With ``return_model=True`` returns ``(info, model)`` instead, so a caller can do
    geometry on the placed grits -- e.g. work out which of them a given depth of cut
    actually engages, which needs the real baked vertices and cannot be had from the
    placement centres alone.
    """
    p.validate()
    if not solids:
        raise DeckError("the grain library is empty")
    if not dry_run:
        os.makedirs(outdir, exist_ok=True)

    sector_deg = p.resolved_sector_deg()
    R = p.outer_radius_mm
    bond_spec = _bond_spec(p, sector_deg)
    full = bond_spec.is_full_circle
    theta_c_deg = 0.0 if full else sector_deg / 2.0

    warnings: list[str] = []
    notes: list[str] = []

    # Grit centres are sampled over a band inset from the bond by one max grain
    # radius. The sampler places centres, not bodies, so sampling the full band
    # leaves grains at the edges hanging past the sector cut faces and the wheel's
    # side faces. A full wheel needs no arc inset -- it has no cut faces -- but
    # still needs the axial one.
    max_r_mm = max(s.bounding_radius_um for s in solids) / UM_PER_MM
    inset = 1.02 * max_r_mm if p.inset_grit_band else 0.0
    face = (min(p.grit_width_window_mm, p.width_mm) if p.grit_width_window_mm > 0
            else p.width_mm)
    grit_width = face - 2 * inset
    if grit_width <= 0:
        raise DeckError(
            "dressed face %g mm cannot hold the largest grain (%.4f mm across) with "
            "an inset band; widen it or set inset_grit_band=False"
            % (face, 2 * max_r_mm))
    full_arc = R * math.radians(sector_deg)
    window = min(p.grit_arc_window_mm, full_arc) if p.grit_arc_window_mm > 0 else full_arc
    if full and window >= full_arc:
        grit_sector_deg, arc_offset_deg = 360.0, 0.0
    else:
        grit_arc = window - 2 * inset
        if grit_arc <= 0:
            raise DeckError(
                "grit band %.4f mm cannot hold the largest grain (%.4f mm across) with "
                "an inset band; widen it or set inset_grit_band=False"
                % (window, 2 * max_r_mm))
        grit_sector_deg = math.degrees(grit_arc / R)
        # Centre the window in the sector, then step past the inset. With no window
        # this reduces to the inset alone, so existing decks are unaffected.
        arc_offset_deg = math.degrees((0.5 * (full_arc - window) + inset) / R)


    if p.grit_mode == "single":
        model, grain_index = _single_grit_model(p, solids, bond_spec, theta_c_deg)
        notes.append(
            "single-grit deck: rotation direction decides whether the grit traverses "
            "the block. It starts at b = %+.4f mm, so a wheel turning toward "
            "%s theta drags it across." % (p.single_grit_offset_mm,
                                           "decreasing" if p.single_grit_offset_mm > 0
                                           else "increasing"))
    else:
        grain_index = None
        grit_spec = WheelSpec(
            diameter_mm=p.diameter_mm, width_mm=grit_width,
            sector_deg=grit_sector_deg, rim_depth_mm=p.rim_depth_mm,
            radial_divisions=bond_spec.radial_divisions,
            axial_divisions=bond_spec.axial_divisions,
            circumferential_divisions_per_deg=(
                bond_spec.circumferential_divisions_per_deg),
        )
        pop = _population(p, solids, sector_deg, inset, grit_spec.surface_area_mm2)
        grits = build_wheel(grit_spec, solids, pop)
        warnings.extend(grits.warnings)
        placed = (grits.placements if full else
                  rotate_placements_about_z(grits.placements, arc_offset_deg))
        for i, pl in enumerate(placed, start=1):
            pl.placement_id = i
        model = _empty_model(bond_spec, pop)
        model = dataclasses.replace(
            model, placements=placed, shapes=grits.shapes,
            requested_grains=grits.requested_grains, achieved_grains=len(placed),
            stats=grits.stats)

    if not model.placements:
        raise DeckError("no grits were placed; lower the density or widen the band")

    wp = workpiece_of(p)
    travel = p.travel_mm if p.travel_mm > 0 else (
        (wp.length_mm + p.travel_margin_mm) if wp else 0.0)

    if dry_run:
        # Everything above decides what the model is; everything below writes it out.
        _wp = workpiece_of(p)
        _tr = p.travel_mm if p.travel_mm > 0 else (
            (_wp.length_mm + p.travel_margin_mm) if _wp else 0.0)
        _an = p.analysis
        _st = ((_an.step_time_s or (_tr / p.surface_speed_mm_s))
               if _an is not None and _an.enabled else 0.0)
        _nf = sum(len(model.shapes[q.shape_index].faces) for q in model.placements)
        _nv = sum(len(model.shapes[q.shape_index].vertices) for q in model.placements)
        _nel = _wp.n_elements() if _wp else 0
        _nwn = 0
        if _wp is not None:
            _a, _b, _c = _wp.divisions()
            _nwn = (_a + 1) * (_b + 1) * (_c + 1)
        # Line lengths the writer actually emits: coordinates at %.12e are ~65 chars a
        # node, an R3D3 facet ~28, a C3D8R hex ~62. Node sets add roughly a fifth.
        _mb = (_nf * 28 + _nv * 65 + _nel * 62 + _nwn * 65) * 1.2 / 1e6
        return (dict(cost=cost_model(p, _wp, _tr, _st, _an, _nel),
                     achieved_areal_density_per_mm2=model.stats.get(
                         "achieved_areal_density_per_mm2"),
                     grit_band_arc_mm=R * math.radians(grit_sector_deg),
                     grit_band_width_mm=grit_width,
                     estimated_mb=_mb,
                     warnings=warnings, notes=notes), model)

    an = p.analysis
    if an is not None and an.enabled:
        an.validate()
        if wp is None:
            raise DeckError("a run-ready deck needs a workpiece to grind")
    step_time = 0.0
    if an is not None and an.enabled:
        step_time = an.step_time_s or (travel / p.surface_speed_mm_s)

    if an is not None and an.enabled and an.depth_of_cut_um <= 0:
        # A safe depth is not knowable until the grits are placed and the block is
        # seated on the tallest of them, so 0 means "choose one that works". Resolved
        # before writing, from the same placement the writer will use. Without this a
        # deck with no infeed runs to completion and grinds nothing -- the most
        # expensive way there is to discover a missing number.
        from .rigid_wheel import place_workpiece
        _clr = place_workpiece(
            model, wp, p.clearance_um, p.wp_position, p.wp_position_deg,
            not bool(getattr(an, "rotation_reversed", False)))["clearance_um"]
        an.depth_of_cut_um = round(p.clearance_um + 0.85 * _clr, 4)
        notes.append("depth of cut auto-set to %.4f um: %.4f um to close the standoff "
                     "plus 85%% of this wheel's %.4f um grain protrusion"
                     % (an.depth_of_cut_um, p.clearance_um, _clr))

    inp = os.path.join(outdir, p.name + ".inp")
    info = write_rigid_wheel_inp(
        inp, model, wp, clearance_um=p.clearance_um,
        wp_position=p.wp_position, wp_position_deg=p.wp_position_deg,
        bond_density_kg_m3=p.bond_density_kg_m3,
        engage_window_mm=travel or None, model_name=p.name,
        analysis=an, step_time_s=step_time,
        surface_speed_mm_s=p.surface_speed_mm_s)

    # The depth of cut can only be judged once the grits are placed and the ground
    # face has been seated on the tallest of them. Feeding deeper than the bond
    # clearance drives the rim itself into the workpiece, which is not grinding and
    # usually just explodes the contact.
    if an is not None and an.enabled:
        # The bond arrives at the work after the standoff *and* the protrusion have
        # been fed through, so the ceiling is the whole gap from face to bond.
        clr = (info["workpiece_ground_radius_mm"] - p.outer_radius_mm) * 1000.0
        # The other end of the same question: with a standoff, the infeed has to close
        # the gap before anything can cut. A job that runs to completion having touched
        # nothing looks like a success until you plot the contact force.
        floor = info.get("min_engaging_infeed_um")
        if floor is not None and an.depth_of_cut_um <= floor:
            raise DeckError(
                "depth of cut %.3f um never reaches the work: with a %.3f um standoff "
                "the nearest grain in the swept band is %.3f um clear, so the wheel "
                "would turn for the whole step without touching. Cut deeper than "
                "%.3f um, or reduce the standoff (clearance_um)."
                % (an.depth_of_cut_um, p.clearance_um, floor, floor))
        if an.depth_of_cut_um >= clr:
            raise DeckError(
                "depth of cut %.3f um exceeds the bond-rim clearance %.3f um for this "
                "wheel, so the bond would hit the workpiece. Use at most %.3f um, or "
                "raise the grit protrusion."
                % (an.depth_of_cut_um, clr, 0.85 * clr))
        if an.depth_of_cut_um < 0.2 * clr:
            notes.append(
                "depth of cut %.3f um is only %.0f%% of the %.3f um clearance; few "
                "grits will reach the work"
                % (an.depth_of_cut_um, 100 * an.depth_of_cut_um / clr, clr))

    # The geometry-only twin: same wheel, no history, for assembling in CAE by hand.
    cae_inp = inp
    if p.also_write_cae_deck and an is not None and an.enabled:
        cae_inp = os.path.join(outdir, p.name + "_cae.inp")
        cae_info = write_rigid_wheel_inp(
            cae_inp, model, wp, clearance_um=p.clearance_um,
            wp_position=p.wp_position, wp_position_deg=p.wp_position_deg,
            bond_density_kg_m3=p.bond_density_kg_m3,
            engage_window_mm=travel or None, model_name=p.name + "_cae",
            analysis=None)
        info["cae_deck"] = cae_inp
        # Each deck needs its own report: the verifiers cross-check the file against
        # it, and size_bytes alone differs between the twin and the run-ready deck.
        cae_info.update(resolved_sector_deg=sector_deg,
                        grit_band_arc_mm=R * math.radians(grit_sector_deg),
                        params=dataclasses.asdict(p))
        with open(os.path.join(outdir, p.name + "_cae_report.json"), "w") as fh:
            json.dump(cae_info, fh, indent=2, default=str)
    write_cae_import_script(os.path.join(outdir, p.name + "_import_into_cae.py"),
                           cae_inp, model_name=p.name)
    _write_placements(os.path.join(outdir, p.name + "_placements.csv"), model)
    # Ships with every run-ready deck, because the job is only half the work: without
    # it there is nothing to turn the .odb into a force, an energy balance or a
    # removed volume.
    if an is not None and an.enabled:
        from .odbpost import write_odb_postprocess_script
        info["postprocess_script"] = write_odb_postprocess_script(
            os.path.join(outdir, p.name + "_postprocess_odb.py"))

    cost = cost_model(p, wp, travel, step_time, an,
                      info["n_workpiece_elements"])

    # ---- CAD --------------------------------------------------------------
    cad: dict = {}
    if p.write_step:
        from .step import StepExportOptions, check_step_solids, write_wheel_step
        sp = os.path.join(outdir, p.name + ".step")
        cad["step"] = write_wheel_step(
            sp, model, StepExportOptions(max_grains=p.step_max_grains, name=p.name))
        cad["step_audit"] = check_step_solids(sp)
    if p.write_stl:
        from .step import wheel_triangles, write_binary_stl
        tri = wheel_triangles(model, max_grains=(p.stl_max_grains or None))
        cad["stl"] = write_binary_stl(os.path.join(outdir, p.name + ".stl"), tri)
    if p.write_grains_step:
        # The measured grains themselves, laid out on a grid rather than at their
        # wheel positions: this is the file for inspecting or measuring one grit.
        from .step import check_step_solids, write_grains_step
        gp = os.path.join(outdir, p.name + "_grains.step")
        cad["grains_step"] = write_grains_step(gp, solids,
                                               max_grains=p.grains_step_max,
                                               laid_out=True)
        cad["grains_step_audit"] = check_step_solids(gp)
    if p.write_grain_stls:
        from .abaqus import write_grain_stl
        d = os.path.join(outdir, "grits_stl")
        os.makedirs(d, exist_ok=True)
        for sol in solids:
            write_grain_stl(os.path.join(d, "grain_%d.stl" % sol.grain_id), sol)
        cad["grain_stls"] = {"dir": d, "count": len(solids)}

    overlaps = check_grain_overlaps(model)
    if overlaps["n_overlapping"]:
        notes.append(
            "%d grit bounding-sphere pairs intersect. The spheres are conservative and "
            "all grits are one rigid body, so this affects appearance only, not contact."
            % overlaps["n_overlapping"])

    info.update(
        params=dataclasses.asdict(p),
        resolved_sector_deg=sector_deg, theta_workpiece_deg_check=theta_c_deg,
        grit_band_inset_mm=inset, grit_band_width_mm=grit_width,
        grit_band_sector_deg=grit_sector_deg,
        grit_band_arc_mm=R * math.radians(grit_sector_deg),
        grit_arc_window_mm=window, grit_face_window_mm=face,
        single_grain_index=grain_index,
        requested_grains=model.requested_grains,
        achieved_areal_density_per_mm2=model.stats.get(
            "achieved_areal_density_per_mm2"),
        grit_overlaps=overlaps, cost=cost, cad=cad,
        warnings=warnings, notes=notes,
    )
    with open(os.path.join(outdir, p.name + "_report.json"), "w") as fh:
        json.dump(info, fh, indent=2, default=str)
    return (info, model) if return_model else info


def _single_grit_model(p: DeckParams, solids: Sequence, bond_spec: WheelSpec,
                       theta_c_deg: float) -> tuple[WheelModel, int]:
    """One chosen grain, seated by the same code path as a full population.

    A patch is sized so the density asks for exactly one grain and the library is cut
    down to the chosen grain, so the sampler cannot pick anything else. Going through
    ``build_wheel`` rather than hand-building the placement keeps the seating
    identical: tip outward, random spin, bounded tilt, protrusion from the same
    truncated normal, then the radius solved so the furthest vertex sits exactly that
    far above the *curved* rim.
    """
    idx = p.single_grain_index
    if idx < 0:
        idx = int(np.argmax([s.mesh_volume_um3 for s in solids]))
    if idx >= len(solids):
        raise DeckError("single_grain_index %d is outside the library of %d grains"
                        % (idx, len(solids)))
    grain = solids[idx]
    R = p.outer_radius_mm

    patch = 0.02
    seed_spec = WheelSpec(diameter_mm=p.diameter_mm, width_mm=patch,
                          sector_deg=math.degrees(patch / R),
                          rim_depth_mm=p.rim_depth_mm, radial_divisions=1,
                          axial_divisions=3,
                          circumferential_divisions_per_deg=100.0)
    pop = GrainPopulationSpec(
        areal_density_per_mm2=1.0 / (patch * patch),
        protrusion_mean=p.protrusion_mean, protrusion_std=p.protrusion_std,
        protrusion_min=p.protrusion_min, protrusion_max=p.protrusion_max,
        max_tilt_deg=p.max_tilt_deg, spacing_factor=p.spacing_factor, seed=p.seed)
    one = build_wheel(seed_spec, [grain], pop)
    if len(one.placements) != 1:
        raise DeckError("expected exactly 1 grit on the seeding patch, got %d"
                        % len(one.placements))

    target = theta_c_deg + math.degrees(p.single_grit_offset_mm / R)
    placed = rotate_placements_about_z(one.placements,
                                       target - one.placements[0].theta_deg)
    pl = placed[0]
    pl.placement_id = 1
    pl.translation_mm[2] = 0.0     # a pure Z shift cannot disturb the radial seating
    pl.axial_mm = 0.0

    model = _empty_model(bond_spec, pop)
    return dataclasses.replace(model, placements=[pl], shapes=[grain],
                               requested_grains=1, achieved_grains=1,
                               stats=one.stats), idx


def _write_placements(path: str, model: WheelModel) -> None:
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["placement_id", "shape_index", "x_mm", "y_mm", "z_mm",
                     "theta_deg", "radius_mm", "protrusion_um", "rot_axis_x",
                     "rot_axis_y", "rot_axis_z", "rot_angle_deg"])
        for q in model.placements:
            wr.writerow([q.placement_id, q.shape_index,
                         "%.9f" % q.translation_mm[0], "%.9f" % q.translation_mm[1],
                         "%.9f" % q.translation_mm[2], "%.6f" % q.theta_deg,
                         "%.6f" % q.radius_mm, "%.4f" % (q.protrusion_mm * 1000),
                         "%.6f" % q.rotation_axis[0], "%.6f" % q.rotation_axis[1],
                         "%.6f" % q.rotation_axis[2],
                         "%.4f" % q.rotation_angle_deg])


# --------------------------------------------------------------------------
# Presets: the two decks that were validated against Abaqus by hand.
# --------------------------------------------------------------------------

PRESETS = {
    "final_712_grit": DeckParams(
        name="wheel_rigid_2mm", sector_mode="arc", arc_length_mm=2.0,
        diameter_mm=50.0, rim_depth_mm=0.012, width_mm=0.030,
        grit_mode="concentration", concentration=100.0, seed=20260731,
        wp_length_mm=0.048, wp_width_mm=0.015, wp_depth_mm=0.006,
        wp_element_size_mm=0.0003),
    "single_grit": DeckParams(
        name="wheel_single_grit", sector_mode="arc", arc_length_mm=2.0,
        diameter_mm=50.0, rim_depth_mm=0.012, width_mm=0.030,
        grit_mode="single", single_grain_index=-1, single_grit_offset_mm=0.015,
        seed=20260731, wp_length_mm=0.048, wp_width_mm=0.015, wp_depth_mm=0.006,
        wp_element_size_mm=0.0003),
}


def hybrid_single_grit(hybrid=None, **over) -> DeckParams:
    """The single-grit deck, run with the hybrid ductile/brittle law.

    Geometrically identical to ``PRESETS['single_grit']`` -- same wheel, same
    grain, same seating -- so a hybrid run and a pure-JH-2 run differ in the
    constitutive law and in nothing else. That is what makes the comparison
    worth doing: any difference in the chip is the switch, not the mesh.

    Not in ``PRESETS`` because it needs a :class:`semgrit.hybrid.HybridParams`
    and because ``_check_presets.py`` compares that dict against two decks
    Abaqus has already accepted.
    """
    from .hybrid import HYBRID_DEPVAR, HybridParams

    an = AnalysisParams(
        enabled=True, material_model="hybrid",
        hybrid=hybrid if hybrid is not None else HybridParams(enabled=True),
        n_depvar=HYBRID_DEPVAR, element_deletion=True,
        depth_of_cut_um=0.0)
    base = dataclasses.replace(PRESETS["single_grit"],
                               name="hybrid_single_grit", analysis=an,
                               also_write_cae_deck=False)
    return dataclasses.replace(base, **over) if over else base
