"""Build the final deck: 2 mm visibly-curved arc, ALL rigid except the workpiece.

Sizing rationale
----------------
* Arc 2.0 mm on a D50 wheel gives a sagitta of 20 um. Against a 12 um rim depth that
  is a 167% bow, so the sector reads unambiguously as an arc rather than a rectangle.
* Only the workpiece is deformable, so it alone sets both the element count and the
  stable time increment. 0.3 um elements put ~8 elements through the deepest cut a
  grit can take, which is what a brittle (JH-2) chip needs in order to form at all --
  the earlier 1.5 um mesh was coarser than the chip itself.
* Grits cover the whole arc at the measured concentration, so fresh grits keep
  arriving as the wheel turns, and the deck still looks like a wheel.
"""
import csv
import dataclasses
import json
import math
import os
import pickle
import sys

sys.path.insert(0, '.')

from semgrit.abaqus import write_cae_import_script
from semgrit.rigid_wheel import write_rigid_wheel_inp
from semgrit.wheel import (UM_PER_MM, GrainPopulationSpec, WheelSpec, build_wheel,
                           check_grain_overlaps)
from semgrit.wheel_workpiece import WorkpieceBlock, rotate_placements_about_z

OUT = 'FINAL_RIGID'
R = 25.0                 # outer radius, mm  (D50 wheel, unchanged)
ARC_MM = 2.0             # <-- the requested 1-2 mm arc
RIM_MM = 0.012           # rim depth: sagitta/rim = 167%, clearly curved
WHEEL_W = 0.030          # axial width

SECTOR_DEG = math.degrees(ARC_MM / R)

# Workpiece: the only deformable body in the model.
WP = WorkpieceBlock(length_mm=0.048, width_mm=0.015, depth_mm=0.006,
                    element_size_mm=0.0003)

V_SURF = 30_000.0        # wheel surface speed, mm/s (30 m/s)
TRAVEL = WP.length_mm + 0.006    # a full pass over the block, plus run-in

lib = pickle.load(open('WHEEL_FIXED/1_measurements/grain_library.pkl', 'rb'))
solids = lib['solids']

bond_spec = WheelSpec(diameter_mm=2 * R, width_mm=WHEEL_W, sector_deg=SECTOR_DEG,
                      rim_depth_mm=RIM_MM, radial_divisions=1, axial_divisions=6,
                      circumferential_divisions_per_deg=200.0 / SECTOR_DEG)
bond = build_wheel(bond_spec, solids,
                   GrainPopulationSpec(areal_density_per_mm2=1.0, seed=1))

# Grit *centres* are sampled over a band inset from the bond by one max grain radius,
# because the sampler places centres, not bodies: sampling the full band leaves the
# grains at the edges hanging past the sector cut faces and the wheel's side faces.
INSET = 1.02 * max(s.bounding_radius_um for s in solids) / UM_PER_MM
grit_spec = WheelSpec(diameter_mm=2 * R, width_mm=WHEEL_W - 2 * INSET,
                      sector_deg=math.degrees((ARC_MM - 2 * INSET) / R),
                      rim_depth_mm=RIM_MM, radial_divisions=1, axial_divisions=6,
                      circumferential_divisions_per_deg=200.0 / SECTOR_DEG)
grits = build_wheel(grit_spec, solids,
                    GrainPopulationSpec(concentration=100.0, seed=20260731))
placed = rotate_placements_about_z(grits.placements, math.degrees(INSET / R))
for i, p in enumerate(placed, start=1):
    p.placement_id = i

# Bond geometry from bond_spec, grits from the inset band.
model = dataclasses.replace(bond, placements=placed, shapes=grits.shapes,
                            requested_grains=grits.requested_grains,
                            stats=grits.stats)

for msg in grits.warnings:
    print('  warn: %s' % msg)
print('  grit band inset by %.6f mm on every side of the bond' % INSET)

os.makedirs(OUT, exist_ok=True)
path = os.path.join(OUT, 'wheel_rigid_2mm.inp')
info = write_rigid_wheel_inp(path, model, WP, clearance_um=0.0,
                             engage_window_mm=TRAVEL,
                             model_name='wheel_rigid_2mm')
write_cae_import_script(os.path.join(OUT, 'wheel_rigid_2mm_import_into_cae.py'),
                        path, model_name='wheel_rigid_2mm')

with open(os.path.join(OUT, 'wheel_rigid_2mm_placements.csv'), 'w',
          newline='') as fh:
    wr = csv.writer(fh)
    wr.writerow(['placement_id', 'shape_index', 'x_mm', 'y_mm', 'z_mm',
                 'theta_deg', 'radius_mm', 'protrusion_um'])
    for p in model.placements:
        wr.writerow([p.placement_id, p.shape_index, '%.9f' % p.translation_mm[0],
                     '%.9f' % p.translation_mm[1], '%.9f' % p.translation_mm[2],
                     '%.6f' % p.theta_deg, '%.6f' % p.radius_mm,
                     '%.4f' % (p.protrusion_mm * 1000)])

# ---- cost model -----------------------------------------------------------
# Dilatational wave speed, which is what Abaqus/Explicit uses for the stable
# increment of a solid element -- not sqrt(E/rho).
E, nu, rho = WP.youngs_modulus_mpa, WP.poisson_ratio, WP.density_kg_m3 * 1e-12
c = math.sqrt(E * (1 - nu) / ((1 + nu) * (1 - 2 * nu) * rho))
dt = WP.element_size_mm / c
T = TRAVEL / V_SURF
inc = T / dt
el_inc = info['n_workpiece_elements'] * inc

# Conservative throughput: 3e5 element-increments/s/core for C3D8R under a
# scalar-loop VUMAT, 70% parallel efficiency, and general contact taken to cost
# a further 60%. Anything faster than this only shortens the run.
RATE_PER_CORE = 3.0e5
overlap = load = 0.0
est = {n: el_inc / (RATE_PER_CORE * n * 0.70 / 1.6) / 3600.0 for n in (4, 8, 16)}

overlaps = check_grain_overlaps(model)

info.update(
    workpiece_element_um=WP.element_size_mm * 1000,
    dilatational_speed_mm_s=c,
    stable_dt_s=dt,
    surface_speed_mm_s=V_SURF,
    omega_rad_s=V_SURF / R,
    travel_mm=TRAVEL,
    step_time_s=T,
    increments=inc,
    element_increments=el_inc,
    est_hours_4core=est[4], est_hours_8core=est[8], est_hours_16core=est[16],
    grit_overlaps=overlaps,
    achieved_areal_density_per_mm2=model.stats['achieved_areal_density_per_mm2'],
    requested_grains=model.requested_grains,
    warnings=model.warnings,
)
json.dump(info, open(os.path.join(OUT, 'wheel_rigid_2mm_report.json'), 'w'),
          indent=2, default=str)

print()
print('WHEEL  R=%.1f mm (D%.0f), %.4f deg = %.3f mm arc, rim %.4f mm, width %g mm'
      % (R, 2 * R, SECTOR_DEG, ARC_MM, RIM_MM, WHEEL_W))
print('       sagitta %.2f um = %.0f%% of the rim depth  -> reads as an arc'
      % (info['sagitta_um'], 100 * info['sagitta_um'] / 1000 / RIM_MM))
print('       ONE rigid body: %d R3D4 shell quads + %d R3D3 grit facets, ref node %d'
      % (info['n_bond_shell_quads'], info['n_grit_facets'], info['wheel_ref_node']))
print('GRITS  %d placed of %d requested at C100, %.0f/mm2 achieved; %d can engage'
      % (info['n_grits'], info['requested_grains'],
         info['achieved_areal_density_per_mm2'], info['n_grits_engaging']))
print('       overlaps: %d of %d bounding-sphere pairs'
      % (overlaps['n_overlapping'], overlaps['n_pairs_checked']))
print('       max protrusion able to engage: %.4f um'
      % info['max_engaging_protrusion_um'])
print('WP     %g x %g x %g mm, %g um elements -> %s C3D8R  (only deformable part)'
      % (WP.length_mm, WP.width_mm, WP.depth_mm, WP.element_size_mm * 1000,
         format(info['n_workpiece_elements'], ',')))
print('       ground face r = %.6f mm, tangent to placement %d, penetration 0'
      % (info['workpiece_ground_radius_mm'], info['governing_grit_placement_id']))
print('RUN    dt = %.3e s, omega = %.1f rad/s, travel %.3f mm -> step %.4e s'
      % (dt, V_SURF / R, TRAVEL, T))
print('       %s increments, %.2e element-increments' % (format(int(inc), ','), el_inc))
print('       estimate: %.1f h on 4 cores, %.1f h on 8, %.1f h on 16'
      % (est[4], est[8], est[16]))
print('FILE   %s  (%.2f MB)' % (path, info['size_bytes'] / 1e6))
