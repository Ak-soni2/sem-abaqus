---
name: sem-grinding-wheel-abaqus-goal
description: Ultimate aim of the SEM grinding-wheel project is an Abaqus-importable model; includes the angular-sector requirement
metadata:
  type: project
---

The SEM → grit-detection → CAD pipeline exists to produce a **grinding wheel model that is imported into Abaqus and simulated** (as of 2026-07-28). The STL/mesh stage is a means to that end, not the deliverable.

Requirements stated by the user:
- Must emit file(s) directly importable into Abaqus, carrying a **large number of micro-scale grits**.
- Must support an **angular sector** parameter: 360° (full wheel), 180° (half), 30°, 25°, etc. — generate only that wedge.
- Grit detection and the overall pipeline are expected to be improved substantially, not just patched.
- User asked that changes be **verified multiple times and checked completely** before reporting done.

**Why:** it dictates the output format and the geometry representation. Boolean-unioning grains into the wheel would be actively wrong — Abaqus needs grains as separate bodies to define contact/tie constraints.

**How to apply:** target native Abaqus `.inp` with `*Part` defined once per unique grain shape plus many `*Instance` blocks (translate + axis-angle rotate), which keeps the deck small at high grit counts. Prefer R3D3 discrete-rigid surface meshes for grains when they are treated as rigid, C3D4/C3D8 when deformable. For sector models, emit node sets on the two cut faces so cyclic-symmetry BCs can be applied. See [[sem-wheel-scale-bug]].
