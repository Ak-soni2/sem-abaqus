---
name: wheel-all-discrete-rigid
description: In the Abaqus grinding model everything must be discrete rigid except the workpiece — the bond too, not just the grits
metadata:
  type: feedback
---

In the grinding-wheel decks, make **every part discrete rigid except the workpiece**. That includes the **bond body**, which I had been writing as deformable C3D8R while only the grits were R3D3. Stated by the user on 2026-07-30, for all future decks.

**Why:** the workpiece is the only thing whose stress and deformation is of interest. A rigid wheel is cheaper and simpler in three ways that matter:
- The bond stops contributing to the Explicit stable time increment. In the D50 deck the bond's 4.4 µm elements gave dt = 9.0e-10 s; once rigid, only the workpiece governs dt.
- Bond elements leave the deformable element count entirely.
- Most importantly, bond + grits can be tied into **one rigid body with a single reference node**, so wheel rotation is one velocity BC on one node instead of constraining 75+ separate rigid grit bodies.

**How to apply:** the bond must become a **surface** mesh (R3D4 quads on the rim faces, or R3D3 triangles), not a solid — `*Rigid Body` needs a surface element set plus a reference node, and `*Solid Section` disappears. Put the reference node on the wheel axis at the origin so a rotational velocity about Z drives the whole wheel directly. Keep the workpiece as C3D8R with its real material. See [[semgrit-package-layout]] and [[sem-grinding-wheel-abaqus-goal]].

**Implemented in `semgrit/rigid_wheel.py`** on 2026-07-31. Two design points there were not obvious in advance and would otherwise be re-derived from scratch:

1. The grits must be **merged into the same part** as the bond, with the instance transform baked in as `v_local @ R.T + t`. A `*Rigid Body` takes one element set, and an element set cannot span parts — so keeping the grits as instances (which is what keeps the file small) makes a single rigid body impossible.
2. Give the body inertia with `*Element, type=MASS` / `type=ROTARYI` on the reference node, **not** `density=`/`thickness=` on `*Rigid Body`. A rigid-element thickness offsets the contact surface by t/2, which destroys exact grit-to-workpiece tangency.
