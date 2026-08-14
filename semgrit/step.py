"""STEP (ISO 10303-21) export, for opening the model in SOLIDWORKS / Fusion / Inventor.

Writes a **faceted B-rep** (`FACETED_BREP` with planar `FACE_SURFACE`s bounded by
`POLY_LOOP`s), which is the canonical STEP representation for polyhedral solids and
what SOLIDWORKS imports as genuine solid bodies rather than as a graphics mesh.

Why STEP rather than STL: SOLIDWORKS treats an STL as a mesh body, so you cannot
readily measure, section or apply features to it. A STEP `FACETED_BREP` arrives as a
real solid body per grain in a multibody part.

Why faceted rather than analytic surfaces: the grains genuinely *are* polyhedra --
they come from a measured outline lofted between planar rings, so every face is
exactly planar and nothing is lost. The wheel rim is faceted too, at the same
angular resolution as its hex mesh, so the CAD and the FE model describe the same
body.

Scale limits
------------
A STEP solid costs roughly 6 entities per triangle plus one point per vertex, and
CAD kernels are far slower per body than an FE solver. Ten thousand grains is fine
in Abaqus and hopeless in SOLIDWORKS. Exports are therefore capped
(``max_grains``) and the caller is told what was dropped -- a small sector is the
intended use for CAD.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

import numpy as np

# Written into the STEP header; the length unit the file declares.
STEP_LENGTH_UNIT = "MILLI"


class StepWriter:
    """Streams faceted solids into a STEP part file.

    Entity ids are assigned sequentially. Vertices are deduplicated per solid, so
    a shared corner is one CARTESIAN_POINT referenced by every face using it.

    Geometry is written **straight to disk** as each solid is added, rather than
    accumulated in memory. Exporting every grain of a full wheel means tens of
    millions of faces and a multi-gigabyte file; buffering that as Python strings
    would need several GB of RAM and fail. STEP does not care about entity order,
    so the product/context boilerplate is appended at the end once the solid ids
    are known.

    Use as a context manager, or call :meth:`finalize` when done.
    """

    def __init__(
        self,
        path: str,
        name: str = "semgrit_model",
        timestamp: str = "2026-01-01T00:00:00",
    ) -> None:
        self.name = name
        self.path = path
        self._next = 1
        self._solids: list[int] = []
        self._solid_names: list[str] = []
        self._closed = False
        self._fh = open(path, "w", encoding="ascii", newline="\n", buffering=1 << 20)
        self._fh.write(
            "ISO-10303-21;\n"
            "HEADER;\n"
            "FILE_DESCRIPTION(('faceted brep of SEM-measured abrasive grains'),'2;1');\n"
            f"FILE_NAME('{_safe(os.path.basename(path))}','{timestamp}',('semgrit'),"
            f"(''),'semgrit','semgrit','');\n"
            "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n"
            "ENDSEC;\n"
            "DATA;\n"
        )

    def __enter__(self) -> "StepWriter":
        return self

    def __exit__(self, *exc) -> None:
        if not self._closed:
            self.finalize()

    # ISO 10303-21 limits a physical line to 32768 characters. A single entity may
    # be continued across lines because the terminator is the semicolon, so long
    # parameter lists are wrapped. Without this, the top-level
    # FACETED_BREP_SHAPE_REPRESENTATION -- which lists every solid in the file --
    # reached 199510 characters on a 20944-grain export and would be rejected.
    _MAX_LINE = 200

    # -- low level ---------------------------------------------------------
    def _add(self, body: str) -> int:
        i = self._next
        self._next += 1
        text = f"#{i}={body};"
        if len(text) <= self._MAX_LINE:
            self._fh.write(text + "\n")
        else:
            self._fh.write(self._wrap(text))
        return i

    @classmethod
    def _wrap(cls, text: str) -> str:
        """Break an over-long entity onto continuation lines at commas."""
        out: list[str] = []
        line = ""
        for piece in text.split(","):
            candidate = piece if not line else line + "," + piece
            if len(candidate) > cls._MAX_LINE and line:
                out.append(line + ",")
                line = piece
            else:
                line = candidate
        if line:
            out.append(line)
        return "\n".join(out) + "\n"

    @staticmethod
    def _f(v: float) -> str:
        """STEP reals must contain a '.' -- ``1`` is invalid, ``1.`` is not.

        12 significant digits, not the more common 9. Grains sit ~50 mm from the
        axis but have ~5 um features, so coordinates are large numbers with small
        differences; at 9 digits the quantisation showed up as a 5e-5 relative
        error in the reconstructed grain volumes.
        """
        if v == 0:
            return "0."
        s = f"{v:.12G}"
        if "E" in s:
            mant, exp = s.split("E")
            if "." not in mant:
                mant += "."
            return f"{mant}E{int(exp)}"
        if "." not in s:
            s += "."
        return s

    def _point(self, p: Sequence[float]) -> int:
        return self._add(
            f"CARTESIAN_POINT('',({self._f(p[0])},{self._f(p[1])},{self._f(p[2])}))"
        )

    def _direction(self, d: Sequence[float]) -> int:
        return self._add(
            f"DIRECTION('',({self._f(d[0])},{self._f(d[1])},{self._f(d[2])}))"
        )

    # -- solids ------------------------------------------------------------
    def add_faceted_solid(
        self, vertices: np.ndarray, faces: Iterable[Sequence[int]], name: str
    ) -> Optional[int]:
        """Add one closed polyhedron. ``faces`` are outward-oriented index loops."""
        verts = np.asarray(vertices, dtype=np.float64)
        face_list = [list(map(int, f)) for f in faces]
        if len(verts) < 4 or not face_list:
            return None

        # One CARTESIAN_POINT per vertex, reused by every face that touches it.
        pids = [self._point(v) for v in verts]

        face_ids: list[int] = []
        for loop in face_list:
            if len(loop) < 3:
                continue
            pts = verts[loop]
            # Plane through the face, oriented along the face normal so the STEP
            # surface normal agrees with the loop winding.
            n = np.cross(pts[1] - pts[0], pts[2] - pts[0])
            ln = float(np.linalg.norm(n))
            if ln <= 0.0:
                continue  # degenerate; skip rather than emit an invalid PLANE
            n = n / ln
            ref = pts[1] - pts[0]
            rl = float(np.linalg.norm(ref))
            if rl <= 0.0:
                continue
            ref = ref / rl

            # Reuse the loop's own first vertex as the plane origin instead of
            # emitting a duplicate CARTESIAN_POINT: one fewer entity per face,
            # which is ~12% off the file size at tens of millions of faces.
            origin = pids[loop[0]]
            axis = self._direction(n)
            refd = self._direction(ref)
            placement = self._add(f"AXIS2_PLACEMENT_3D('',#{origin},#{axis},#{refd})")
            plane = self._add(f"PLANE('',#{placement})")
            poly = self._add(
                "POLY_LOOP('',(" + ",".join(f"#{pids[i]}" for i in loop) + "))"
            )
            bound = self._add(f"FACE_OUTER_BOUND('',#{poly},.T.)")
            face_ids.append(self._add(f"FACE_SURFACE('',(#{bound}),#{plane},.T.)"))

        if not face_ids:
            return None
        shell = self._add(
            "CLOSED_SHELL('',(" + ",".join(f"#{i}" for i in face_ids) + "))"
        )
        solid = self._add(f"FACETED_BREP('{_safe(name)}',#{shell})")
        self._solids.append(solid)
        self._solid_names.append(name)
        return solid

    # -- output ------------------------------------------------------------
    def finalize(self) -> dict:
        """Append the product/context boilerplate and close the file."""
        if self._closed:
            raise ValueError("writer already finalized")
        if not self._solids:
            self._fh.close()
            self._closed = True
            raise ValueError("no solids to write")

        # The product/context boilerplate must come after the geometry because ids
        # are sequential and it has to reference the solids.
        unit_len = self._add(
            f"(NAMED_UNIT(*)SI_UNIT(.{STEP_LENGTH_UNIT}.,.METRE.)LENGTH_UNIT())"
        )
        unit_ang = self._add("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
        unit_sol = self._add("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())")
        tol = self._add(
            f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#{unit_len},"
            f"'distance_accuracy_value','confusion accuracy')"
        )
        geo_ctx = self._add(
            f"(GEOMETRIC_REPRESENTATION_CONTEXT(3)"
            f"GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{tol}))"
            f"GLOBAL_UNIT_ASSIGNED_CONTEXT((#{unit_len},#{unit_ang},#{unit_sol}))"
            f"REPRESENTATION_CONTEXT('',''))"
        )
        app_ctx = self._add("APPLICATION_CONTEXT('automotive design')")
        self._add(
            f"APPLICATION_PROTOCOL_DEFINITION('international standard',"
            f"'automotive_design',2000,#{app_ctx})"
        )
        prod_ctx = self._add(f"PRODUCT_CONTEXT('',#{app_ctx},'mechanical')")
        pdef_ctx = self._add(f"PRODUCT_DEFINITION_CONTEXT('part definition',#{app_ctx},'design')")
        product = self._add(
            f"PRODUCT('{_safe(self.name)}','{_safe(self.name)}','',(#{prod_ctx}))"
        )
        formation = self._add(f"PRODUCT_DEFINITION_FORMATION('','',#{product})")
        pdef = self._add(f"PRODUCT_DEFINITION('design','',#{formation},#{pdef_ctx})")
        pshape = self._add(f"PRODUCT_DEFINITION_SHAPE('','',#{pdef})")
        rep = self._add(
            f"FACETED_BREP_SHAPE_REPRESENTATION('{_safe(self.name)}',("
            + ",".join(f"#{i}" for i in self._solids)
            + f"),#{geo_ctx})"
        )
        self._add(f"SHAPE_DEFINITION_REPRESENTATION(#{pshape},#{rep})")

        self._fh.write("ENDSEC;\nEND-ISO-10303-21;\n")
        self._fh.close()
        self._closed = True
        return {
            "path": self.path,
            "n_solids": len(self._solids),
            "n_entities": self._next - 1,
            "size_bytes": os.path.getsize(self.path),
        }

    # Kept so existing call sites read naturally; the path is fixed at construction.
    def write(self, path: Optional[str] = None) -> dict:
        if path is not None and os.path.abspath(path) != os.path.abspath(self.path):
            raise ValueError(
                f"StepWriter streams to {self.path!r}; pass the path to the "
                f"constructor instead"
            )
        return self.finalize()


def estimate_step_size_bytes(n_solids: int, faces_per_solid: float) -> int:
    """Rough STEP size for a faceted export. Measured at ~430 bytes per face."""
    return int(n_solids * faces_per_solid * 430)


def _safe(s: str) -> str:
    """Escape a STEP string literal."""
    return str(s).replace("\\", "").replace("'", "''")


# --------------------------------------------------------------------------
# Hex-mesh boundary extraction
# --------------------------------------------------------------------------

# Outward-oriented faces of a positively-oriented C3D8 hex with nodes 0-3 on the
# -z face (wound so their normal is +z) and 4-7 above them. Verified against the
# rim mesh, whose hexes all have positive volume.
_HEX_FACES = (
    (0, 3, 2, 1),   # -z
    (4, 5, 6, 7),   # +z
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)


def surface_from_hexes(hexes: np.ndarray) -> np.ndarray:
    """Outward-oriented boundary quads of a hex mesh.

    A quad on the boundary belongs to exactly one element; interior quads appear
    twice and cancel.
    """
    quads = np.vstack([hexes[:, list(f)] for f in _HEX_FACES])
    keys = np.sort(quads, axis=1)
    _, inv, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    return quads[counts[inv] == 1]


def signed_volume(vertices: np.ndarray, faces: Sequence[Sequence[int]]) -> float:
    """Volume enclosed by a closed polygonal surface, via the divergence theorem.

    Negative means the faces are wound inward.
    """
    total = 0.0
    for loop in faces:
        p = vertices[list(loop)]
        for k in range(1, len(loop) - 1):
            a, b, c = p[0], p[k], p[k + 1]
            total += float(np.dot(np.cross(b - a, c - a), a)) / 6.0
    return total


# --------------------------------------------------------------------------
# High-level exports
# --------------------------------------------------------------------------

@dataclass
class StepExportOptions:
    max_grains: int = 200
    """Cap on grain bodies. **0 or negative means export every grain.**

    A cap exists because CAD kernels are far slower per solid body than an FE
    solver is per element: a STEP B-rep costs roughly 430 bytes per planar face,
    so a full wheel runs to gigabytes and will not open. Set it to 0 when you want
    the complete model regardless."""
    include_body: bool = True
    include_grains: bool = True
    name: str = "grinding_wheel"

    @property
    def unlimited(self) -> bool:
        return self.max_grains <= 0


def write_wheel_step(
    path: str, model, options: Optional[StepExportOptions] = None
) -> dict:
    """Export a wheel model (bond rim + grains) as a STEP part.

    Grains are placed at their assembly positions, so the STEP describes the same
    geometry as the Abaqus deck rather than a pile of loose parts.
    """
    from .wheel import UM_PER_MM, _rotation_matrix

    options = options or StepExportOptions()
    w = StepWriter(path, options.name)
    warnings: list[str] = []

    if options.include_body:
        quads = surface_from_hexes(model.body_hexes)
        vol = signed_volume(model.body_nodes, quads)
        loops = [list(q) for q in quads]
        if vol < 0:
            loops = [l[::-1] for l in loops]
            vol = -vol
        w.add_faceted_solid(model.body_nodes, loops, "BOND_RIM")

    n_grains = 0
    if options.include_grains and model.placements:
        placements = model.placements
        if not options.unlimited and len(placements) > options.max_grains:
            # Take an angular spread rather than the first N, so the exported
            # patch is representative instead of a clump at theta=0.
            idx = np.linspace(0, len(placements) - 1, options.max_grains).astype(int)
            idx = np.unique(idx)
            warnings.append(
                f"exported {len(idx)} of {len(placements)} grains to STEP "
                f"(max_grains={options.max_grains}); grains sampled evenly across "
                f"the sector. The Abaqus .inp still contains all "
                f"{len(placements)}."
            )
            placements = [placements[i] for i in idx]

        for p in placements:
            s = model.shapes[p.shape_index]
            rot = _rotation_matrix(p.rotation_axis, np.radians(p.rotation_angle_deg))
            v = (s.vertices - s.centroid_um) / UM_PER_MM @ rot.T + p.translation_mm
            if w.add_faceted_solid(v, s.faces, f"GRAIN_{p.placement_id}") is not None:
                n_grains += 1

    info = w.finalize()
    info["n_grain_bodies"] = n_grains
    info["warnings"] = warnings
    return info


def write_grains_step(
    path: str, solids: Sequence, max_grains: int = 200, laid_out: bool = True
) -> dict:
    """Export individual grain solids as a multibody STEP part.

    ``laid_out`` spreads the grains on a grid instead of leaving them at their
    measured positions in the micrograph, which is what you want for inspecting or
    measuring single grains in CAD.
    """
    w = StepWriter(path, "sem_grains")
    warnings: list[str] = []
    use = list(solids)
    if 0 < max_grains < len(use):
        warnings.append(f"exported {max_grains} of {len(use)} grain solids")
        use = use[:max_grains]

    if laid_out:
        pitch = 1.25 * max(
            (float(max(s.extent_um()[0], s.extent_um()[1])) for s in use), default=1.0
        )
        cols = max(int(np.ceil(np.sqrt(len(use)))), 1)

    n = 0
    for k, s in enumerate(use):
        v = s.vertices.copy()
        if laid_out:
            v = v - s.centroid_um
            v[:, 0] += (k % cols) * pitch
            v[:, 1] += (k // cols) * pitch
            v[:, 2] -= v[:, 2].min()
        base = os.path.splitext(os.path.basename(s.source_image))[0]
        if w.add_faceted_solid(v, s.faces, f"{base}_g{s.grain_id}") is not None:
            n += 1
    info = w.finalize()
    info["n_grain_bodies"] = n
    info["warnings"] = warnings
    info["units"] = "micrometres" if not laid_out else "micrometres"
    return info


# --------------------------------------------------------------------------
# Reading back (verification)
# --------------------------------------------------------------------------

def read_step_faceted_solids(path: str) -> list[dict]:
    """Parse the faceted solids out of a STEP file written by :class:`StepWriter`.

    Deliberately independent of the writer: it resolves the entity graph from the
    text, so it can confirm that what landed on disk really describes the intended
    closed solids. Returns one dict per solid with ``name``, ``vertices``,
    ``faces`` and ``volume``.

    Not a general STEP parser -- it understands the faceted subset this module
    emits (FACETED_BREP / CLOSED_SHELL / FACE_SURFACE / FACE_OUTER_BOUND /
    POLY_LOOP / CARTESIAN_POINT).
    """
    import re

    with open(path, "r", encoding="ascii", errors="replace") as fh:
        text = fh.read()

    entities: dict[str, str] = {}
    for m in re.finditer(r"#(\d+)\s*=\s*(.*?);", text, re.S):
        entities[m.group(1)] = m.group(2).strip()

    def refs(body: str) -> list[str]:
        return re.findall(r"#(\d+)", body)

    points: dict[str, tuple[float, float, float]] = {}
    poly_loops: dict[str, list[str]] = {}
    outer_bounds: dict[str, str] = {}
    face_surfaces: dict[str, str] = {}
    closed_shells: dict[str, list[str]] = {}
    breps: dict[str, tuple[str, str]] = {}

    for eid, body in entities.items():
        if body.startswith("CARTESIAN_POINT"):
            nums = re.findall(r"[-+]?\d*\.\d*(?:E[-+]?\d+)?|[-+]?\d+\.", body)
            vals = [float(n) for n in nums[:3]]
            if len(vals) == 3:
                points[eid] = (vals[0], vals[1], vals[2])
        elif body.startswith("POLY_LOOP"):
            poly_loops[eid] = refs(body)
        elif body.startswith("FACE_OUTER_BOUND"):
            r = refs(body)
            if r:
                outer_bounds[eid] = r[0]
        elif body.startswith("FACE_SURFACE"):
            r = refs(body)
            if r:
                face_surfaces[eid] = r[0]   # first ref is the bound
        elif body.startswith("CLOSED_SHELL"):
            closed_shells[eid] = refs(body)
        elif body.startswith("FACETED_BREP"):
            name = re.match(r"FACETED_BREP\s*\(\s*'([^']*)'", body)
            r = refs(body)
            if r:
                breps[eid] = (name.group(1) if name else "", r[0])

    out: list[dict] = []
    for _, (name, shell_id) in sorted(breps.items(), key=lambda kv: int(kv[0])):
        loops: list[list[str]] = []
        for face_id in closed_shells.get(shell_id, []):
            bound_id = face_surfaces.get(face_id)
            if bound_id is None:
                continue
            loop_id = outer_bounds.get(bound_id)
            if loop_id is None:
                continue
            pts = [p for p in poly_loops.get(loop_id, []) if p in points]
            if len(pts) >= 3:
                loops.append(pts)
        if not loops:
            continue
        unique = sorted({p for lp in loops for p in lp}, key=int)
        index = {p: k for k, p in enumerate(unique)}
        verts = np.array([points[p] for p in unique], dtype=np.float64)
        faces = [[index[p] for p in lp] for lp in loops]
        out.append(
            {
                "name": name,
                "vertices": verts,
                "faces": faces,
                "volume": signed_volume(verts, faces),
            }
        )
    return out


ISO_10303_MAX_LINE = 32768


def check_step_solids(path: str) -> dict:
    """Structural audit of a STEP file: closure, orientation, finite coordinates."""
    solids = read_step_faceted_solids(path)
    issues: list[str] = []

    # ISO 10303-21 caps a physical line at 32768 characters. Exceeding it can make
    # a CAD translator reject or truncate the file, and it is easy to do by
    # accident: the representation entity lists every solid in the model.
    with open(path, "r", encoding="ascii", errors="replace") as fh:
        max_line = max((len(l.rstrip("\n")) for l in fh), default=0)
    if max_line > ISO_10303_MAX_LINE:
        issues.append(
            f"longest line is {max_line} chars, over the ISO 10303-21 limit of "
            f"{ISO_10303_MAX_LINE}"
        )
    n_open = 0
    n_negative = 0
    for s in solids:
        # Every directed edge must appear exactly once, and its reverse exactly
        # once, for a closed consistently-oriented surface.
        seen: dict[tuple[int, int], int] = {}
        for loop in s["faces"]:
            for k in range(len(loop)):
                e = (loop[k], loop[(k + 1) % len(loop)])
                seen[e] = seen.get(e, 0) + 1
        bad = sum(1 for e, c in seen.items() if c != 1 or seen.get((e[1], e[0]), 0) != 1)
        if bad:
            n_open += 1
        if s["volume"] <= 0:
            n_negative += 1
        if not np.isfinite(s["vertices"]).all():
            issues.append(f"{s['name']}: non-finite coordinates")

    if n_open:
        issues.append(f"{n_open} solids are not closed / consistently oriented")
    if n_negative:
        issues.append(f"{n_negative} solids have non-positive volume (inward normals)")

    return {
        "path": path,
        "ok": not issues,
        "issues": issues,
        "n_solids": len(solids),
        "n_faces": int(sum(len(s["faces"]) for s in solids)),
        "total_volume": float(sum(s["volume"] for s in solids)),
        "min_volume": float(min((s["volume"] for s in solids), default=0.0)),
        "max_line_chars": int(max_line),
    }


# --------------------------------------------------------------------------
# STL (fallback / visualisation)
# --------------------------------------------------------------------------

def write_binary_stl(path: str, triangles: np.ndarray) -> dict:
    """Binary STL from an (N, 3, 3) array of triangle vertices.

    Kept as a fallback: an STL opens anywhere and stays manageable at grain counts
    where a STEP B-rep would not, at the cost of arriving as a mesh body rather
    than a solid.
    """
    import struct

    tris = np.asarray(triangles, dtype=np.float32).reshape(-1, 3, 3)
    with open(path, "wb") as fh:
        fh.write(b"semgrit binary STL".ljust(80, b"\0"))
        fh.write(struct.pack("<I", len(tris)))
        for t in tris:
            n = np.cross(t[1] - t[0], t[2] - t[0])
            ln = np.linalg.norm(n)
            n = n / ln if ln > 0 else np.zeros(3, np.float32)
            fh.write(struct.pack("<3f", *n.astype(np.float32)))
            for v in t:
                fh.write(struct.pack("<3f", *v))
            fh.write(struct.pack("<H", 0))
    return {"path": path, "n_triangles": int(len(tris)), "size_bytes": os.path.getsize(path)}


def wheel_triangles(model, max_grains: Optional[int] = None) -> np.ndarray:
    """All triangles of a wheel model in assembly coordinates, for STL export."""
    from .wheel import UM_PER_MM, _rotation_matrix

    out: list[np.ndarray] = []
    quads = surface_from_hexes(model.body_hexes)
    for q in quads:
        p = model.body_nodes[list(q)]
        out.append(p[[0, 1, 2]])
        out.append(p[[0, 2, 3]])

    placements = model.placements
    if max_grains is not None and len(placements) > max_grains:
        idx = np.unique(np.linspace(0, len(placements) - 1, max_grains).astype(int))
        placements = [placements[i] for i in idx]
    for p in placements:
        s = model.shapes[p.shape_index]
        rot = _rotation_matrix(p.rotation_axis, np.radians(p.rotation_angle_deg))
        v = (s.vertices - s.centroid_um) / UM_PER_MM @ rot.T + p.translation_mm
        out.extend(v[f] for f in s.faces)
    return np.asarray(out, dtype=np.float64)
