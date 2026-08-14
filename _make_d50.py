"""Build the D50 / 10 deg deck: wide visibly-curved bond, grits only where they bite."""
import dataclasses
import json
import math
import os
import pickle
import sys

sys.path.insert(0, '.')

from semgrit.abaqus import write_cae_import_script
from semgrit.wheel import GrainPopulationSpec, WheelSpec, build_wheel
from semgrit.wheel_workpiece import (
    WorkpieceBlock,
    rotate_placements_about_z,
    write_wheel_workpiece_inp,
)

R = 25.0                 # outer radius (mm) -> D50 wheel
SECTOR_DEG = 10.0        # bond arc
GRIT_ARC_MM = 0.12       # grits confined to this arc length
WIDTH = 0.03             # axial width, matched to the workpiece

lib = pickle.load(open('WHEEL_FIXED/1_measurements/grain_library.pkl', 'rb'))
solids = lib['solids']

wp = WorkpieceBlock(length_mm=0.05, width_mm=WIDTH, depth_mm=0.015,
                    element_size_mm=0.0015)

# --- bond: the full 10 deg arc, coarse circumferentially (it is only there to hold
#     the grits and to show the curvature, so it does not need fine elements) ---
bond_spec = WheelSpec(diameter_mm=2 * R, width_mm=WIDTH, sector_deg=SECTOR_DEG,
                      rim_depth_mm=0.02, radial_divisions=2, axial_divisions=3,
                      circumferential_divisions_per_deg=6.0)
bond = build_wheel(bond_spec, solids,
                   GrainPopulationSpec(areal_density_per_mm2=1.0, seed=1))

# --- grits: a small window at true C100, then rotated to the middle of the arc ---
grit_spec = WheelSpec(diameter_mm=2 * R, width_mm=WIDTH,
                      sector_deg=math.degrees(GRIT_ARC_MM / R),
                      rim_depth_mm=0.02, radial_divisions=2, axial_divisions=3,
                      circumferential_divisions_per_deg=400.0)
grits = build_wheel(grit_spec, solids,
                    GrainPopulationSpec(concentration=100.0, seed=20260730))

offset = SECTOR_DEG / 2.0 - grit_spec.sector_deg / 2.0
placed = rotate_placements_about_z(grits.placements, offset)
for i, p in enumerate(placed, start=1):
    p.placement_id = i

# bond geometry + grit placements
model = dataclasses.replace(bond, placements=placed, shapes=grits.shapes)

arc = R * math.radians(SECTOR_DEG)
sag = arc * arc / (8 * R)
print('bond   : R=%.1f mm (D%.0f), %.1f deg = %.3f mm arc x %g mm wide'
      % (R, 2 * R, SECTOR_DEG, arc, WIDTH))
print('         sagitta %.1f um vs 20 um rim depth -> %.0f%% bow, clearly curved'
      % (sag * 1000, 100 * sag * 1000 / 20))
print('         %d C3D8R hexes' % len(bond.body_hexes))
print('grits  : %d over a %.3f mm window at %.0f/mm2 (true C100), centred at %.2f deg'
      % (len(placed), GRIT_ARC_MM, grits.stats['achieved_areal_density_per_mm2'],
         SECTOR_DEG / 2))

os.makedirs('D50_MODEL', exist_ok=True)
path = 'D50_MODEL/wheel_D50_10deg.inp'
info = write_wheel_workpiece_inp(path, model, wp, clearance_um=0.0)
write_cae_import_script('D50_MODEL/wheel_D50_10deg_import_into_cae.py', path,
                        model_name='wheel_D50_10deg')

E, rho = 50000.0, 2.65e-9
dt = wp.element_size_mm / math.sqrt(E / rho)
T = (GRIT_ARC_MM - wp.length_mm) / 30000.0
info.update(stable_dt_s=dt, suggested_step_time_s=T, increments=T / dt,
            total_elements=info['n_grit_facets'] + info['n_bond_elements']
            + info['n_workpiece_elements'],
            sector_deg=SECTOR_DEG, outer_radius_mm=R,
            arc_length_mm=arc, sagitta_um=sag * 1000,
            omega_rad_s_at_30ms=30000.0 / R)
json.dump(info, open('D50_MODEL/wheel_D50_10deg_report.json', 'w'), indent=2, default=str)

print()
for k in ('size_bytes', 'n_grits', 'n_grit_parts', 'n_grit_facets', 'n_bond_elements',
          'n_workpiece_elements', 'total_elements', 'tallest_grit_tip_radius_mm',
          'workpiece_ground_radius_mm', 'bond_to_workpiece_gap_um',
          'stable_dt_s', 'suggested_step_time_s', 'increments', 'omega_rad_s_at_30ms'):
    v = info[k]
    if isinstance(v, int):
        print('  %-30s %s' % (k, format(v, ',')))
    elif isinstance(v, float):
        print('  %-30s %.6g' % (k, v))
    else:
        print('  %-30s %s' % (k, v))
