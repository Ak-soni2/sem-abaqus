"""Presentation deck: a 30 deg slice of a D50 wheel that reads as a wheel on sight.

Design choices, and why
-----------------------
**A 30 deg sector, not a few degrees.** What makes a shape read as "cut from a wheel"
is the two radial faces visibly converging. At 30 deg they are unmistakable, and the
outer face bows 0.86 mm across a 13.09 mm chord (6.5%), which the eye sees as a curve.
Every earlier deck used a fraction of a degree and had to rely on a thin rim to look
curved, which is why it kept reading as a rectangle.

**A deep rim: bore at r = 10 mm.** The bond is a rigid *surface*, so radial depth costs
a few hundred quads and nothing in the solver. A 30 deg slice running from r = 10 to
r = 25 mm is a pie slice of a wheel rather than a sliver of its rim -- the single
biggest visual win available, and it is free.

**A 3 mm thick slice, dressed over a 0.6 mm band.** Thickness is what stops the slice
looking like a thin washer segment, and for a rigid surface it costs a few hundred
quads. Dressing the whole 3 mm face at 5000/mm^2 would be 30,000 grains for no benefit,
because the workpiece is 0.25 mm wide -- so the band is 0.6 mm, centred.

**Grits over a 1.5 mm arc window, not the whole 13 mm arc.** At true C100 the whole arc would
need ~65,000 grains and hundreds of MB. The window is 5x wider than anything that
touches the work, so the dressed band is plainly visible at intermediate zoom, and at
the zoom where individual grits resolve you cannot see the rest of the arc anyway.

**A 0.25 mm thick block with a graded depth mesh.** The chip is removed *into the depth*,
so the depth direction is what resolves chip thickness -- and at a uniform 5 um a 2.5 um
cut sat entirely inside one element. Grading fixes both complaints at once: 1 um layers
over the top 5 um where the chip forms, coarsening to 20 um in the body. The block is
0.25 mm thick instead of 0.1 mm, so it is a solid rectangular specimen rather than a
plate that would bend under the grinding load, and it is plainly visible when zoomed.
Grading costs nothing in time: dt follows the *smallest* element, and 1 um in the depth
just matches the 1 um cutting size that already set it.

**Only the workpiece is deformable**, as always: one rigid body, one reference node.

Budget: ~12 h on 4 cores, ~6 h on 8, which is the stated allowance.
"""
import json
import math
import os
import pickle
import sys

sys.path.insert(0, '.')

import numpy as np

from semgrit.build_deck import DeckParams, build_deck

OUT = 'PRESENT'
R = 25.0                 # D50 wheel
SECTOR_DEG = 30.0        # the two cut faces converge at 30 deg -> obviously a sector
BORE_R = 10.0            # rim depth 15 mm: a pie slice, not a sliver
WIDTH = 3.0              # axial thickness: a solid chunk, not a thin plate
FACE = 0.6               # dressed band across the face (workpiece is 0.25 mm wide)
GRIT_ARC = 1.5           # dressed band along the arc, mm
GRIT_DENSITY = 5000.0    # grains/mm2; the bounding-sphere spacing caps this near 5200

lib = pickle.load(open('WHEEL_FIXED/1_measurements/grain_library.pkl', 'rb'))
solids = lib['solids']

params = DeckParams(
    name='wheel_present_30deg',
    diameter_mm=2 * R, sector_mode='angle', sector_deg=SECTOR_DEG,
    rim_depth_mm=R - BORE_R, width_mm=WIDTH,
    shell_circumferential_divisions=300, shell_axial_divisions=30,
    shell_radial_divisions=12,
    grit_mode='areal_density', areal_density_per_mm2=GRIT_DENSITY,
    grit_arc_window_mm=GRIT_ARC, grit_width_window_mm=FACE, seed=20260801,
    # workpiece: a real specimen, and the only deformable body
    wp_length_mm=0.40, wp_width_mm=0.25, wp_depth_mm=0.25,
    wp_element_size_mm=0.001,          # cutting direction
    wp_element_size_width_mm=0.0025,   # axial
    wp_element_size_depth_mm=0.001,    # the fine surface layer, where the chip forms
    wp_surface_layer_mm=0.005,         # 5 um of 1 um layers at the ground face
    wp_depth_growth=1.3,               # then coarsen into the body
    wp_max_depth_element_mm=0.020,     # capped so deep elements do not become slivers
    travel_mm=0.18, surface_speed_mm_s=30_000.0, cores=8,
)

info, model = build_deck(params, solids, OUT, return_model=True)
c = info['cost']

# --- what depth of cut actually engages, inside the swept band -----------------
# From the real baked grit vertices, not from centre+protrusion: the furthest point of
# a tilted grain is not at its centre angle, and that approximation gave clearances
# wrong by more than a micron - enough to recommend a useless depth of cut.
from semgrit.rigid_wheel import bake_grit, quantise

thc = math.radians(info['theta_workpiece_deg'])
e_r = np.array([math.cos(thc), math.sin(thc), 0.0])
e_t = np.array([-math.sin(thc), math.cos(thc), 0.0])
basis = np.column_stack([e_r, e_t, np.array([0.0, 0.0, 1.0])])
half_b = params.travel_mm / 2.0 + 0.01
half_z = params.wp_width_mm / 2.0
clr = []
for pl in model.placements:
    v = quantise(bake_grit(model, pl)) @ basis
    sel = (np.abs(v[:, 1]) <= half_b) & (np.abs(v[:, 2]) <= half_z)
    if sel.any():
        clr.append((info['workpiece_ground_radius_mm'] - v[sel, 0].max()) * 1000.0)
swept = np.sort(np.array(clr)) if clr else np.array([np.nan])
band = np.array(clr) if clr else np.array([])

print('WHEEL   D%g, %.1f deg sector, r = %.1f .. %.1f mm (rim %.1f mm), %g mm thick'
      % (2 * R, info['resolved_sector_deg'], BORE_R, R, info['rim_depth_mm'], WIDTH))
print('        dressed band %.2f mm arc x %.2f mm of the %g mm face'
      % (info['grit_band_arc_mm'], info['grit_band_width_mm'], WIDTH))
print('        arc %.3f mm, sagitta %.0f um = %.1f%% of the chord -> visibly curved'
      % (info['arc_length_mm'], info['sagitta_um'],
         100 * info['sagitta_um'] / 1000 / info['arc_length_mm']))
print('        ONE rigid body: %s shell quads + %s grit facets, ref node %d'
      % (format(info['n_bond_shell_quads'], ','),
         format(info['n_grit_facets'], ','), info['wheel_ref_node']))
print('GRITS   %s placed of %s requested (%.0f/mm2 achieved, %.0f%% of true C100)'
      % (format(info['n_grits'], ','), format(info['requested_grains'], ','),
         info['achieved_areal_density_per_mm2'],
         100 * info['achieved_areal_density_per_mm2'] / 28254))
print('        %s lie under the block' % format(info['n_grits_engaging'], ','))
print('WP      %g x %g x %g mm -> %s C3D8R, %d x %d x %d  (only deformable part)'
      % (params.wp_length_mm, params.wp_width_mm, params.wp_depth_mm,
         format(info['n_workpiece_elements'], ','), *c['element_divisions']))
print('        element %.2f cutting x %.2f axial x depth %.2f (surface) -> %.2f (back) um'
      % (c['element_size_cutting_mm'] * 1000, c['element_size_axial_mm'] * 1000,
         c['depth_layer_min_mm'] * 1000, c['depth_layer_max_mm'] * 1000))
print('        %.2f um sets dt; %.1f elements through a %.1f um cut'
      % (c['governing_element_size_mm'] * 1000,
         2.5 / (c['depth_layer_min_mm'] * 1000), 2.5))
print('        ground face r = %.6f mm, tangent to placement %s, penetration 0'
      % (info['workpiece_ground_radius_mm'], info['governing_grit_placement_id']))
print('        bond rim clearance %.3f um  (hard ceiling on the depth of cut)'
      % info['max_engaging_protrusion_um'])
if not np.isnan(swept[0]):
    q = lambda k: swept[min(k, len(swept) - 1)]
    print('DEPTH   in the %.0f um swept band, %d grits; the tallest sits %.3f um down'
          % (params.travel_mm * 1000, len(swept), swept[0]))
    print('        ae to engage 1 / 10 / 50 grits: %.2f / %.2f / %.2f um'
          % (q(0), q(9), q(49)))
print('RUN     dt = %.3e s, omega = %.1f rad/s (%.0f rpm), travel %.3f mm'
      % (c['stable_dt_s'], c['omega_rad_s'], c['rpm'], c['travel_mm']))
print('        step %.4e s = %s increments, %.2e element-increments'
      % (c['step_time_s'], format(int(c['increments']), ','),
         c['element_increments']))
print('        estimate ' + ', '.join('%s core %.1f h' % (k, v) for k, v in
      sorted(c['est_hours'].items(), key=lambda kv: int(kv[0]))))
for m in info['warnings'] + info['notes']:
    print('  note: %s' % m)
print()
print('FILE    %s  (%.1f MB)' % (info['path'], info['size_bytes'] / 1e6))
