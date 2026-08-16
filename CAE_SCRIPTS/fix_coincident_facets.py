"""Repair the "facet has coincident nodes" error in an already-set-up CAE model.

Run from Abaqus/CAE:   File -> Run Script...

Fixes this without rebuilding, so materials, sections, the step, contact, BCs and the
VUMAT assignment all survive:

  ***ERROR: A facet that is part of a general contact surface has coincident
            nodes. ... See nodes 64 and 65 on element 125 of instance G-1992.

Cause: a few rigid-grain facets were generated with two nodes ~1e-7 mm apart. Their
area is ~5e-12 mm^2 -- zero to any physical purpose -- but Abaqus will not put such a
facet in a general contact surface.

Method: delete those facets. Their area is zero, so no coordinate changes and no
geometry is lost. Abaqus was already going to discard them ("These collapsed faces
will be ignored, creating a seam"); it just refused to start first. Deletion is
preferred over merging nodes, because merging *moves* nodes and perturbs the
neighbouring facets.

Editing a *part* repairs every instance of it, so a handful of edits fixes thousands
of grits.

SAFETY MEASURES, each guarding a specific way this could go wrong:

  1. deleteUnreferencedNodes=OFF. Every rigid grain part carries a reference node
     that belongs to no element (*Rigid Body, ref node=GRAIN_REF). Deleting
     unreferenced nodes would remove it and break the rigid body definition.
  2. Only surface facet elements (R3D3/R3D4/S3/S3R/S4/S4R) are considered. Solids are
     skipped entirely: for a C3D8 hex, walking the nodes cyclically crosses body
     diagonals rather than edges, so the distances would be meaningless.
  3. Node coordinates come from element.getNodes(), not from connectivity indices.
     Abaqus returns internal indices for native meshes but node *labels* for orphan
     meshes; indexing part.nodes[label] would silently return the wrong node.
  4. A part is skipped if the fix would remove more than MAX_FRACTION of its facets.
     That is the backstop: if anything about the geometry read is misinterpreted, the
     script refuses to mass-delete and says so instead.
  5. DRY_RUN lets you see the counts before anything is modified.
"""

import abaqusConstants as AC
from abaqus import mdb

# Facets with any edge shorter than this (model length units; mm here) are removed.
# Observed bad edges were <= 5.4e-7 mm and the next genuine edge was 1.4e-6 mm, so
# anything in that gap works. 1e-5 mm matches the tolerance the corrected generator
# enforces and also clears a few marginal facets.
TOLERANCE = 1.0e-5

# Refuse to delete more than this share of a part's facets.
MAX_FRACTION = 0.25

# Set True to report only, changing nothing.
DRY_RUN = False

MODELS_TO_DO = []          # empty = every model

FACET_TYPES = ("R3D3", "R3D4", "S3", "S3R", "S4", "S4R", "SFM3D3", "SFM3D4")


def shortest_edge(coords):
    """Shortest edge of a triangle or quad, walking the perimeter."""
    n = len(coords)
    best = None
    for i in range(n):
        a = coords[i]
        b = coords[(i + 1) % n]
        d = ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
        if best is None or d < best:
            best = d
    return best


def facet_coords(elem):
    """Node coordinates of a facet element, or None if they cannot be read.

    Uses getNodes() so the label-versus-index ambiguity in .connectivity cannot
    silently pick the wrong nodes.
    """
    try:
        return [n.coordinates for n in elem.getNodes()]
    except Exception:
        return None


def scan_part(part):
    """Returns (bad_labels, n_facets, worst_edge, skipped_types, read_failures)."""
    try:
        elems = part.elements
    except Exception:
        return [], 0, None, set(), 0
    if not len(elems):
        return [], 0, None, set(), 0

    bad = []
    n_facets = 0
    worst = None
    skipped = set()
    failures = 0

    for e in elems:
        etype = str(e.type)
        if etype not in FACET_TYPES:
            skipped.add(etype)
            continue
        n_facets += 1
        coords = facet_coords(e)
        if coords is None or len(coords) < 3:
            failures += 1
            continue
        repeated = len(set(tuple(c) for c in coords)) < len(coords)
        s = shortest_edge(coords)
        if worst is None or s < worst:
            worst = s
        if repeated or s < TOLERANCE:
            bad.append(e.label)
    return bad, n_facets, worst, skipped, failures


def main():
    names = MODELS_TO_DO if MODELS_TO_DO else sorted(mdb.models.keys())
    print("=" * 76)
    print("Removing zero-area contact facets   tolerance=%.1e   %s" % (
        TOLERANCE, "DRY RUN - nothing will change" if DRY_RUN else "LIVE"))
    print("=" * 76)

    for mname in names:
        if mname not in mdb.models:
            continue
        m = mdb.models[mname]
        if not m.parts:
            continue
        print("")
        print("MODEL %r  (%d parts)" % (mname, len(m.parts)))

        total_bad = total_facets = n_deleted = 0
        edited = []
        refused = []
        all_skipped = set()
        total_failures = 0
        worst_before = None

        for pname in sorted(m.parts.keys()):
            part = m.parts[pname]
            bad, n_facets, worst, skipped, failures = scan_part(part)
            all_skipped |= skipped
            total_failures += failures
            total_facets += n_facets
            if worst is not None and (worst_before is None or worst < worst_before):
                worst_before = worst
            if not bad:
                continue
            total_bad += len(bad)

            # Backstop: never mass-delete.
            if n_facets and float(len(bad)) / n_facets > MAX_FRACTION:
                refused.append((pname, len(bad), n_facets))
                continue

            if DRY_RUN:
                edited.append((pname, len(bad)))
                continue
            try:
                arr = part.elements.sequenceFromLabels(bad)
                # OFF, not ON: the rigid-body reference node is referenced by no
                # element and must survive.
                part.deleteElement(elements=arr, deleteUnreferencedNodes=AC.OFF)
                n_deleted += len(bad)
                edited.append((pname, len(bad)))
            except Exception as exc:
                refused.append((pname, len(bad), n_facets))
                print("    could not edit %s: %s" % (pname, exc))

        print("  facet elements examined     : %d" % total_facets)
        if all_skipped:
            print("  non-facet types skipped     : %s" % ", ".join(sorted(all_skipped)))
        if total_failures:
            print("  facets whose nodes could not be read : %d  <-- investigate before"
                  " trusting the result" % total_failures)
        print("  degenerate facets found     : %d" % total_bad)
        print("  facets deleted              : %d" % n_deleted)
        if edited:
            print("  parts edited                : %d" % len(edited))
            for pname, n in edited[:15]:
                print("      %-16s %d" % (pname, n))
            if len(edited) > 15:
                print("      ... and %d more" % (len(edited) - 15))
        if refused:
            print("  parts REFUSED (over %.0f%% of facets, or an error) : %d" % (
                100 * MAX_FRACTION, len(refused)))
            for pname, n, tot in refused:
                print("      %-16s %d of %d" % (pname, n, tot))
            print("      nothing was deleted from those. Do not raise MAX_FRACTION")
            print("      without checking why the count is so high.")
        if worst_before is not None:
            print("  shortest facet edge before  : %.4e" % worst_before)

        if n_deleted and not DRY_RUN:
            worst_after = None
            for pname in sorted(m.parts.keys()):
                _, _, w, _, _ = scan_part(m.parts[pname])
                if w is not None and (worst_after is None or w < worst_after):
                    worst_after = w
            if worst_after is not None:
                print("  shortest facet edge after   : %.4e   %s" % (
                    worst_after,
                    "OK" if worst_after >= TOLERANCE else "<-- STILL SHORT"))
            try:
                m.rootAssembly.regenerate()
                print("  assembly regenerated")
            except Exception as exc:
                print("  assembly regenerate failed: %s" % exc)
                print("  do Assembly -> Regenerate by hand before submitting")

    print("")
    print("=" * 76)
    if DRY_RUN:
        print("DRY RUN complete - nothing was modified. Set DRY_RUN = False to apply.")
    else:
        print("Done. Save the model, then resubmit the job.")
    print("=" * 76)


main()
