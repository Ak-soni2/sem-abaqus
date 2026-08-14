"""Abaqus ``.inp`` export of the grinding-wheel model.

Scalability
-----------
Each distinct grain *shape* is written once as a ``*Part``; every grain on the
wheel is a ``*Instance`` of one of those parts carrying its own translation and
rotation. A 30 deg sector with 20,000 grains built from 96 measured shapes
therefore stores 97 meshes, not 20,000 -- the deck stays small enough to open.

Instance transform convention
-----------------------------
``*Instance`` accepts one translation line and one rotation line::

    *Instance, name=G-1, part=GRAIN-3
      cx, cy, cz
      cx, cy, cz, cx+ax, cy+ay, cz+az, angle
    *End Instance

Abaqus applies the **translation first, then the rotation** about the axis given
in assembly coordinates. Writing the axis so that it passes through the grain's
own final centre therefore leaves the centre fixed and simply spins the grain in
place. This is the same pattern Abaqus/CAE emits when you translate and then
rotate an instance.

Because a wrong assumption here would silently misplace every grain,
:func:`semgrit.verify.verify_inp_roundtrip` re-reads the written file and
reconstructs the positions under this convention. ``grain_parts="baked"`` is
available as an escape hatch: it writes pre-rotated coordinates and no instance
transform at all, so placement cannot depend on the convention.

Boolean unions are deliberately not performed. Grains must stay separate bodies
for contact or tie constraints; fusing them into the bond would erase the
interfaces the analysis is about.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Optional, Sequence, TextIO

import numpy as np

from .grain3d import GrainSolid
from .wheel import UM_PER_MM, GrainPlacement, WheelModel, _rotation_matrix

# mm / N / tonne / MPa / s. Density in tonne/mm^3 = (kg/m^3) * 1e-12.
KGM3_TO_TONNE_MM3 = 1e-12


@dataclass
class MaterialSpec:
    name: str
    youngs_modulus_mpa: float
    poisson_ratio: float
    density_kg_m3: float
    note: str = ""


# Representative room-temperature properties. Stated so the deck runs; confirm
# against your own material data before drawing conclusions from results.
MATERIALS = {
    "diamond": MaterialSpec("DIAMOND", 1_050_000.0, 0.07, 3520.0, "typical synthetic diamond"),
    "b4c": MaterialSpec("B4C", 450_000.0, 0.17, 2520.0, "typical hot-pressed boron carbide"),
    "cbn": MaterialSpec("CBN", 680_000.0, 0.12, 3480.0, "typical cubic boron nitride"),
    "vitrified_bond": MaterialSpec(
        "VITRIFIED_BOND", 60_000.0, 0.22, 2500.0, "typical vitrified bond"
    ),
    "resin_bond": MaterialSpec("RESIN_BOND", 8_000.0, 0.30, 1800.0, "typical resin bond"),
    "metal_bond": MaterialSpec("METAL_BOND", 110_000.0, 0.30, 8400.0, "typical bronze bond"),
}


@dataclass
class AbaqusExportOptions:
    """Controls how the model is written."""

    grain_element: str = "C3D4"
    """'C3D4' deformable tetrahedra, or 'R3D3' discrete-rigid surface triangles.
    R3D3 is much cheaper and is usually the right choice: the abrasive is far
    harder than the workpiece, so grain deformation rarely matters."""

    grain_parts: str = "shared"
    """'shared' = one part per measured shape plus instance transforms (compact).
    'baked' = one part per placement with pre-rotated coordinates (larger, but
    independent of the instance-transform convention)."""

    grain_material: str = "diamond"
    bond_material: str = "vitrified_bond"

    include_body: bool = True
    include_step_template: bool = True
    include_contact: bool = True

    wheel_part_name: str = "WHEEL_BOND"
    model_name: str = "GRINDING_WHEEL"

    float_fmt: str = "{: .8e}"

    def __post_init__(self) -> None:
        if self.grain_element not in ("C3D4", "R3D3"):
            raise ValueError("grain_element must be 'C3D4' or 'R3D3'")
        if self.grain_parts not in ("shared", "baked"):
            raise ValueError("grain_parts must be 'shared' or 'baked'")
        if self.grain_material not in MATERIALS:
            raise ValueError(f"unknown grain material {self.grain_material!r}")
        if self.bond_material not in MATERIALS:
            raise ValueError(f"unknown bond material {self.bond_material!r}")


def _fmt_row(fh: TextIO, first: int, values: Sequence[float], fmt: str) -> None:
    fh.write(f"{first:d}, " + ", ".join(fmt.format(v) for v in values) + "\n")


def _write_int_set(fh: TextIO, ids: Sequence[int], per_line: int = 8) -> None:
    """Abaqus tolerates at most 16 entries per line; 8 keeps it readable."""
    ids = list(ids)
    for i in range(0, len(ids), per_line):
        fh.write(", ".join(str(v) for v in ids[i : i + per_line]) + "\n")


def _grain_part_name(index: int) -> str:
    return f"GRAIN-{index + 1}"


def _write_grain_part(
    fh: TextIO,
    name: str,
    vertices_mm: np.ndarray,
    solid: GrainSolid,
    opts: AbaqusExportOptions,
) -> None:
    """One grain shape as an Abaqus part, authored centred on its own centroid."""
    fh.write(f"*Part, name={name}\n")
    fh.write("*Node\n")
    for i, v in enumerate(vertices_mm, start=1):
        _fmt_row(fh, i, v, opts.float_fmt)

    if opts.grain_element == "C3D4":
        fh.write("*Element, type=C3D4\n")
        for e, tet in enumerate(solid.tets, start=1):
            fh.write(f"{e}, {tet[0] + 1}, {tet[1] + 1}, {tet[2] + 1}, {tet[3] + 1}\n")
        n_el = len(solid.tets)
        fh.write(f"*Elset, elset=GRAIN_ALL, generate\n1, {n_el}, 1\n")
        fh.write("*Solid Section, elset=GRAIN_ALL, material="
                 f"{MATERIALS[opts.grain_material].name}\n,\n")
    else:
        fh.write("*Element, type=R3D3\n")
        for e, tri in enumerate(solid.faces, start=1):
            fh.write(f"{e}, {tri[0] + 1}, {tri[1] + 1}, {tri[2] + 1}\n")
        n_el = len(solid.faces)
        fh.write(f"*Elset, elset=GRAIN_ALL, generate\n1, {n_el}, 1\n")
        # A discrete rigid body needs a reference node; put it at the centroid,
        # which is the part origin by construction.
        ref = len(vertices_mm) + 1
        fh.write("*Node\n")
        _fmt_row(fh, ref, np.zeros(3), opts.float_fmt)
        fh.write(f"*Nset, nset=GRAIN_REF\n{ref},\n")
        fh.write(f"*Rigid Body, ref node=GRAIN_REF, elset=GRAIN_ALL\n")

    fh.write(f"*Nset, nset=GRAIN_NODES, generate\n1, {len(vertices_mm)}, 1\n")
    if opts.grain_element == "R3D3":
        fh.write("*Surface, type=ELEMENT, name=GRAIN_SURF\nGRAIN_ALL, SPOS\n")
    fh.write("*End Part\n")


def _write_body_part(fh: TextIO, model: WheelModel, opts: AbaqusExportOptions) -> None:
    fh.write(f"*Part, name={opts.wheel_part_name}\n")
    fh.write("*Node\n")
    for i, v in enumerate(model.body_nodes, start=1):
        _fmt_row(fh, i, v, opts.float_fmt)

    fh.write("*Element, type=C3D8\n")
    for e, h in enumerate(model.body_hexes, start=1):
        fh.write(f"{e}, " + ", ".join(str(int(n) + 1) for n in h) + "\n")

    n_el = len(model.body_hexes)
    fh.write(f"*Elset, elset=BOND_ALL, generate\n1, {n_el}, 1\n")
    fh.write(
        f"*Solid Section, elset=BOND_ALL, material={MATERIALS[opts.bond_material].name}\n,\n"
    )

    for set_name, ids in model.node_sets.items():
        if set_name.startswith("_"):
            continue
        fh.write(f"*Nset, nset={set_name}\n")
        _write_int_set(fh, [int(i) + 1 for i in ids])

    # Element-based surfaces. Face numbering follows the C3D8 node ordering used
    # in build_rim_mesh: S4 outer, S6 bore, S3/S5 sector cuts, S1/S2 axial.
    face_map = {
        "_EL_OUTER": ("BOND_OUTER_SURF", "S4"),
        "_EL_BORE": ("BOND_BORE_SURF", "S6"),
        "_EL_ZMIN": ("BOND_ZMIN_SURF", "S1"),
        "_EL_ZMAX": ("BOND_ZMAX_SURF", "S2"),
        "_EL_SECTOR_START": ("BOND_SECTOR_START_SURF", "S3"),
        "_EL_SECTOR_END": ("BOND_SECTOR_END_SURF", "S5"),
    }
    for key, (surf_name, face) in face_map.items():
        ids = model.node_sets.get(key)
        if ids is None or len(ids) == 0:
            continue
        elset = f"_{surf_name}_EL"
        fh.write(f"*Elset, elset={elset}\n")
        _write_int_set(fh, [int(i) + 1 for i in ids])
        fh.write(f"*Surface, type=ELEMENT, name={surf_name}\n{elset}, {face}\n")

    fh.write("*End Part\n")


def _write_instance(
    fh: TextIO,
    inst_name: str,
    part_name: str,
    placement: Optional[GrainPlacement],
    opts: AbaqusExportOptions,
) -> None:
    fh.write(f"*Instance, name={inst_name}, part={part_name}\n")
    if placement is not None:
        c = placement.translation_mm
        f = opts.float_fmt
        fh.write(", ".join(f.format(v) for v in c) + "\n")
        if abs(placement.rotation_angle_deg) > 1e-12:
            a = placement.rotation_axis
            p2 = c + a
            fh.write(
                ", ".join(f.format(v) for v in c)
                + ", "
                + ", ".join(f.format(v) for v in p2)
                + f", {placement.rotation_angle_deg: .8e}\n"
            )
    fh.write("*End Instance\n")


def write_inp(
    path: str,
    model: WheelModel,
    opts: Optional[AbaqusExportOptions] = None,
    provenance: Optional[dict] = None,
) -> dict:
    """Write the model as an Abaqus input deck. Returns a summary dict."""
    opts = opts or AbaqusExportOptions()
    spec = model.spec

    # Grain vertices are authored in mm, centred on each grain's centroid, so an
    # instance translation places the centroid directly.
    shape_vertices = [
        (s.vertices - s.centroid_um) / UM_PER_MM for s in model.shapes
    ]

    used_shapes = sorted({p.shape_index for p in model.placements})

    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("*Heading\n")
        fh.write(f"** {opts.model_name}: SEM-measured abrasive grains on a "
                 f"{spec.sector_deg:g} deg wheel sector\n")
        fh.write("** Units: mm, N, tonne, MPa, s. Wheel axis = Z. "
                 "Sector spans theta=0..%g deg from +X.\n" % spec.sector_deg)
        fh.write("** Generated by semgrit. Grain geometry is derived from SEM "
                 "measurements; see the provenance block below.\n")
        if provenance:
            for k, v in provenance.items():
                fh.write(f"** {k}: {v}\n")
        fh.write(f"** wheel diameter (mm): {spec.diameter_mm:g}\n")
        fh.write(f"** wheel width (mm): {spec.width_mm:g}\n")
        fh.write(f"** sector (deg): {spec.sector_deg:g}\n")
        fh.write(f"** rim depth (mm): {spec.outer_radius_mm - spec.inner_radius_mm:g}\n")
        fh.write(f"** grains placed: {len(model.placements)}\n")
        fh.write(f"** distinct grain shapes: {len(used_shapes)}\n")
        fh.write(f"** grain element type: {opts.grain_element}\n")
        fh.write("** Instance convention: translation applied first, then "
                 "rotation about an axis through the translated centre.\n")
        for w in model.warnings:
            fh.write(f"** WARNING: {w}\n")

        # ---------------- parts ----------------
        if opts.grain_parts == "shared":
            for idx in used_shapes:
                _write_grain_part(
                    fh, _grain_part_name(idx), shape_vertices[idx],
                    model.shapes[idx], opts,
                )
        else:
            for p in model.placements:
                v = shape_vertices[p.shape_index]
                r = _rotation_matrix(p.rotation_axis, math.radians(p.rotation_angle_deg))
                _write_grain_part(
                    fh, f"GRAIN-P{p.placement_id}", v @ r.T,
                    model.shapes[p.shape_index], opts,
                )

        if opts.include_body:
            _write_body_part(fh, model, opts)

        # ---------------- assembly ----------------
        fh.write("*Assembly, name=ASSEMBLY\n")
        if opts.include_body:
            _write_instance(fh, "BOND-1", opts.wheel_part_name, None, opts)

        for p in model.placements:
            if opts.grain_parts == "shared":
                _write_instance(
                    fh, f"G-{p.placement_id}", _grain_part_name(p.shape_index), p, opts
                )
            else:
                # Coordinates already rotated; translate only.
                fh.write(f"*Instance, name=G-{p.placement_id}, "
                         f"part=GRAIN-P{p.placement_id}\n")
                fh.write(", ".join(opts.float_fmt.format(v) for v in p.translation_mm) + "\n")
                fh.write("*End Instance\n")

        # Assembly-level grouping so one statement can address every grain.
        # Instance-qualified member names are used because a set declared with
        # instance=... would only ever cover that single instance.
        if model.placements and opts.grain_element == "R3D3":
            names = [f"G-{p.placement_id}.GRAIN_REF" for p in model.placements]
            fh.write("** Reference node of every rigid grain. A discrete rigid body\n")
            fh.write("** has no stiffness, so each one must be constrained, tied or\n")
            fh.write("** driven -- otherwise the step is singular.\n")
            fh.write("*Nset, nset=ALL_GRAIN_REF\n")
            for i in range(0, len(names), 8):
                fh.write(", ".join(names[i : i + 8]) + "\n")

        if opts.include_body:
            for set_name in model.node_sets:
                if set_name.startswith("_"):
                    continue
                fh.write(f"*Nset, nset=A_{set_name}, instance=BOND-1\n{set_name},\n")

            if not spec.is_full_circle:
                fh.write("** ------------------------------------------------------\n")
                fh.write("** Sector model: the two cut faces are meshed identically,\n")
                fh.write("** so entry k of A_SECTOR_FACE_START pairs with entry k of\n")
                fh.write("** A_SECTOR_FACE_END. Apply cyclic symmetry with *Equation\n")
                fh.write("** on those node pairs, or clamp both faces for a local\n")
                fh.write("** grinding study.\n")
                fh.write("** ------------------------------------------------------\n")

        fh.write("*End Assembly\n")

        # ---------------- materials ----------------
        for key in {opts.grain_material, opts.bond_material}:
            m = MATERIALS[key]
            fh.write(f"** {m.note}\n")
            fh.write(f"*Material, name={m.name}\n")
            fh.write(f"*Density\n{m.density_kg_m3 * KGM3_TO_TONNE_MM3: .8e},\n")
            fh.write(f"*Elastic\n{m.youngs_modulus_mpa: .8e}, {m.poisson_ratio: .4f}\n")

        if opts.include_contact:
            fh.write("*Surface Interaction, name=GRAIN_BOND_INT\n")
            fh.write("*Friction\n0.2,\n")
            # General contact belongs to the *initial* step, i.e. the model data
            # section before the first *Step. Abaqus/Standard rejects it inside a
            # step with "General Contact (Std) can only be defined in the initial
            # step", and the contact definition is then silently dropped on import.
            fh.write("** General contact, defined in the initial step as Abaqus/Standard\n")
            fh.write("** requires. Applies to all exterior faces, which avoids having to\n")
            fh.write("** enumerate contact pairs across tens of thousands of grains.\n")
            fh.write("*Contact\n")
            fh.write("*Contact Inclusions, ALL EXTERIOR\n")
            fh.write("*Contact Property Assignment\n ,  , GRAIN_BOND_INT\n")

        if opts.include_step_template:
            _write_step_template(fh, model, opts)

    return {
        "path": path,
        "n_grain_parts": (
            len(used_shapes) if opts.grain_parts == "shared" else len(model.placements)
        ),
        "n_instances": len(model.placements) + (1 if opts.include_body else 0),
        "n_grain_elements": (
            sum(model.shapes[p.shape_index].n_tets for p in model.placements)
            if opts.grain_element == "C3D4"
            else sum(len(model.shapes[p.shape_index].faces) for p in model.placements)
        ),
        "n_body_elements": len(model.body_hexes) if opts.include_body else 0,
        "grain_element": opts.grain_element,
        "grain_parts": opts.grain_parts,
    }


def _write_step_template(fh: TextIO, model: WheelModel, opts: AbaqusExportOptions) -> None:
    """A minimal runnable step. Loads and BCs are the analyst's to set."""
    spec = model.spec
    fh.write("** ==========================================================\n")
    fh.write("** TEMPLATE STEP -- edit before using for real analysis.\n")
    fh.write("** Boundary conditions, loads and step type depend entirely on\n")
    fh.write("** the grinding process being modelled; what follows only makes\n")
    fh.write("** the deck syntactically complete and importable.\n")
    fh.write("** ==========================================================\n")
    if opts.include_body:
        fh.write("*Boundary\n")
        if spec.inner_radius_mm > 0:
            fh.write("A_WHEEL_BORE, 1, 3\n")
        else:
            fh.write("A_WHEEL_ZMIN, 1, 3\n")
    if opts.grain_element == "R3D3" and model.placements:
        fh.write("** Rigid grains pinned to the bond for import validation. Replace\n")
        fh.write("** with a *Tie, *Coupling or prescribed motion for real analysis.\n")
        fh.write("*Boundary\nALL_GRAIN_REF, 1, 6\n")
    fh.write("*Step, name=Step-1, nlgeom=YES\n")
    fh.write("*Static\n0.1, 1.0, 1e-08, 0.1\n")
    fh.write("*Output, field\n*Node Output\nU, RF\n")
    fh.write("*Element Output, directions=YES\nS, E\n")
    fh.write("*Output, history, variable=PRESELECT\n")
    fh.write("*End Step\n")


# --------------------------------------------------------------------------
# Auxiliary exports
# --------------------------------------------------------------------------

def write_cae_import_script(path: str, inp_path: str, model_name: str = "WheelModel") -> str:
    """Write a script that loads an ``.inp`` into Abaqus/CAE as a *complete* model.

    Needed because CAE's GUI has no reliable menu route for this:

    * ``File -> Import -> Part`` reads the ``*Part`` definitions but ignores the
      ``*Assembly`` block, so every grain arrives as a part and none are placed --
      the wheel looks bare even though the deck is correct.
    * ``File -> Import -> Model`` only offers ``*.cae`` in its file filter.

    ``mdb.ModelFromInputFile`` is the API those menus wrap, and it does read the
    assembly. Run the generated file from CAE via **File -> Run Script**.
    """
    inp_name = os.path.basename(inp_path)
    body = f'''"""Load a semgrit wheel deck into Abaqus/CAE as a complete model.

Run from Abaqus/CAE:   File -> Run Script...   and pick this file.

Why a script: CAE's File -> Import -> Part reads only the *Part blocks and skips
the *Assembly, so the grains arrive unplaced; and File -> Import -> Model does not
accept .inp. ModelFromInputFile reads both.

The .inp is located relative to this script, so the pair can be copied to any
folder or machine without editing paths.
"""
import os

from abaqus import mdb, session
from abaqusConstants import *   # noqa: F401,F403  (CAE scripts rely on these names)

INP_NAME = "{inp_name}"
MODEL = "{model_name}"

# Extra places to look, in case the script and the deck get separated.
SEARCH_ROOTS = [r"D:\\temp", r"C:\\temp"]


def find_inp():
    cands = []
    try:
        cands.append(os.path.dirname(os.path.abspath(__file__)))
    except NameError:          # some CAE versions do not set __file__
        pass
    cands.append(os.getcwd())
    for root in SEARCH_ROOTS:
        cands.append(root)
        if os.path.isdir(root):
            try:
                for entry in os.listdir(root):
                    sub = os.path.join(root, entry)
                    if os.path.isdir(sub):
                        cands.append(sub)
                        for entry2 in os.listdir(sub):
                            sub2 = os.path.join(sub, entry2)
                            if os.path.isdir(sub2):
                                cands.append(sub2)
            except OSError:
                pass
    seen = set()
    for c in cands:
        if not c or c in seen:
            continue
        seen.add(c)
        p = os.path.join(c, INP_NAME)
        if os.path.isfile(p):
            return p
    raise IOError(
        "could not find %s. Looked in:\\n  %s\\nPut this script next to the .inp, "
        "or edit SEARCH_ROOTS at the top." % (INP_NAME, "\\n  ".join(sorted(seen)))
    )


INP = find_inp()

if MODEL in mdb.models:
    del mdb.models[MODEL]

print("reading " + INP)
mdb.ModelFromInputFile(name=MODEL, inputFileName=INP)

m = mdb.models[MODEL]
a = m.rootAssembly
print("parts      : %d" % len(m.parts))
print("instances  : %d" % len(a.instances))
grits = len([k for k in a.instances.keys() if k.startswith("G-")])
print("grit instances placed on the wheel : %d" % grits)

# Show the assembled model rather than a single part.
try:
    vp = session.viewports[session.currentViewportName]
    vp.setValues(displayedObject=a)
    vp.view.fitView()
except Exception as exc:              # viewport naming differs between versions
    print("could not set the viewport automatically: %s" % exc)
    print("switch Module to Assembly by hand instead")

print("")
print("If 'grit instances placed' matches the grain count in the report JSON,")
print("every grit is in the model. Zoom onto the outer diameter to see them:")
print("they are only a few microns across.")
'''
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write(body)
    return path


def write_placement_csv(path: str, model: WheelModel) -> None:
    """Per-grain placement table.

    Everything needed to rebuild the assembly outside Abaqus, or to drive a DEM
    or analytical kinematic model of the same wheel.
    """
    import csv

    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "placement_id", "shape_index", "source_image", "grain_id",
                "x_mm", "y_mm", "z_mm",
                "rot_axis_x", "rot_axis_y", "rot_axis_z", "rot_angle_deg",
                "theta_deg", "axial_mm", "radius_mm", "protrusion_mm",
                "bounding_radius_mm", "grain_height_um", "grain_volume_um3",
                "equivalent_diameter_um", "feret_max_um", "feret_min_um",
            ]
        )
        for p in model.placements:
            s = model.shapes[p.shape_index]
            m = s.measurement
            w.writerow(
                [
                    p.placement_id, p.shape_index, s.source_image, s.grain_id,
                    f"{p.translation_mm[0]:.6f}", f"{p.translation_mm[1]:.6f}",
                    f"{p.translation_mm[2]:.6f}",
                    f"{p.rotation_axis[0]:.6f}", f"{p.rotation_axis[1]:.6f}",
                    f"{p.rotation_axis[2]:.6f}", f"{p.rotation_angle_deg:.6f}",
                    f"{p.theta_deg:.4f}", f"{p.axial_mm:.6f}", f"{p.radius_mm:.6f}",
                    f"{p.protrusion_mm:.6f}", f"{p.bounding_radius_mm:.6f}",
                    f"{s.height_um:.4f}", f"{s.mesh_volume_um3:.4f}",
                    f"{m.equivalent_diameter_um:.4f}" if m else "",
                    f"{m.feret_max_um:.4f}" if m else "",
                    f"{m.feret_min_um:.4f}" if m else "",
                ]
            )


def write_grain_stl(path: str, solid: GrainSolid, scale: float = 1.0) -> None:
    """ASCII STL of one grain, for quick visual inspection."""
    v = solid.vertices * scale
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write(f"solid grain_{solid.grain_id}\n")
        for tri in solid.faces:
            p = v[tri]
            n = np.cross(p[1] - p[0], p[2] - p[0])
            ln = np.linalg.norm(n)
            n = n / ln if ln > 0 else np.zeros(3)
            fh.write(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n")
            fh.write("    outer loop\n")
            for q in p:
                fh.write(f"      vertex {q[0]:.6e} {q[1]:.6e} {q[2]:.6e}\n")
            fh.write("    endloop\n  endfacet\n")
        fh.write(f"endsolid grain_{solid.grain_id}\n")
