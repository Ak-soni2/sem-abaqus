"""wheel_rigid_2mm with a full-size 0.2 x 0.2 x 0.2 mm workpiece.

The wheel is byte-for-byte the one that was accepted: D50, 2.000 mm arc, 12 um rim,
30 um width, C100 grits, seed 20260731, one discrete rigid body. Only the workpiece
changes.

Why the mesh is anisotropic
---------------------------
0.2 mm cubed at the original 0.3 um would be 297 million elements -- roughly 93 days on
8 cores. The three directions do not need the same resolution:

* **cutting 1.0 um** - the chip forms here, so this is the direction to spend on. It
  also sets the stable increment, so it is the expensive one.
* **axial 2.5 um** - 12 elements across the 30 um wheel width is enough to shape the
  groove, and this direction costs elements without touching the increment.
* **depth 5.0 um** - the cut is ~2 um deep into a 200 um block, so below the first few
  microns the mesh is carrying nothing.

640,000 elements, ~1.5 h on 8 cores. A *uniform* mesh of the same element count would
be 2.5 um, which puts less than one element through the chip and cannot fracture at
all; the anisotropy buys 3x the cutting-direction resolution for the same money.

The pass length stays 54 um, as in the accepted deck. A bigger specimen is for removing
edge effects, not for grinding all of it -- letting travel scale with the block would
quadruple the run for no new physics.
"""
import json
import os
import pickle
import sys

sys.path.insert(0, '.')

from semgrit.build_deck import DeckParams, build_deck

OUT = 'WP200'
lib = pickle.load(open('WHEEL_FIXED/1_measurements/grain_library.pkl', 'rb'))
solids = lib['solids']

params = DeckParams(
    name='wheel_rigid_2mm_wp200',
    # ---- wheel: identical to the accepted wheel_rigid_2mm ----
    diameter_mm=50.0, sector_mode='arc', arc_length_mm=2.0,
    rim_depth_mm=0.012, width_mm=0.030,
    shell_circumferential_divisions=200, shell_axial_divisions=6,
    shell_radial_divisions=1,
    grit_mode='concentration', concentration=100.0, seed=20260731,
    # ---- workpiece: the requested full size ----
    wp_length_mm=0.2, wp_width_mm=0.2, wp_depth_mm=0.2,
    wp_element_size_mm=0.001,          # cutting direction
    wp_element_size_width_mm=0.0025,   # axial
    wp_element_size_depth_mm=0.005,    # into the depth
    # ---- same grinding pass as before ----
    travel_mm=0.054, surface_speed_mm_s=30_000.0, cores=8,
)

info = build_deck(params, solids, OUT)
c = info['cost']
R = info['outer_radius_mm']

print('WHEEL   D%g, %.4f deg = %.4f mm arc, rim %.4f mm, width %g mm  (unchanged)'
      % (2 * R, info['resolved_sector_deg'], info['arc_length_mm'],
         info['rim_depth_mm'], params.width_mm))
print('        sagitta %.2f um = %.0f%% of the rim depth'
      % (info['sagitta_um'], 100 * info['sagitta_um'] / 1000 / info['rim_depth_mm']))
print('        ONE rigid body: %s shell quads + %s grit facets, ref node %d'
      % (format(info['n_bond_shell_quads'], ','),
         format(info['n_grit_facets'], ','), info['wheel_ref_node']))
print('GRITS   %s placed of %s requested at C100 (%.0f/mm2)'
      % (format(info['n_grits'], ','), format(info['requested_grains'], ','),
         info['achieved_areal_density_per_mm2']))
print('        %d lie under the block; tallest reaching protrusion %.4f um'
      % (info['n_grits_engaging'], info['max_engaging_protrusion_um']))
print('WP      %g x %g x %g mm -> %s C3D8R, %d x %d x %d'
      % (params.wp_length_mm, params.wp_width_mm, params.wp_depth_mm,
         format(info['n_workpiece_elements'], ','), *c['element_divisions']))
print('        element %.4f cutting x %.4f axial x %.4f depth um; %.4f um sets dt'
      % (c['element_size_cutting_mm'] * 1000, c['element_size_axial_mm'] * 1000,
         c['element_size_depth_mm'] * 1000, c['governing_element_size_mm'] * 1000))
print('        ground face r = %.6f mm, tangent to placement %s, penetration 0'
      % (info['workpiece_ground_radius_mm'], info['governing_grit_placement_id']))
print('RUN     dt = %.3e s, omega = %.1f rad/s (%.0f rpm), travel %.4f mm'
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
