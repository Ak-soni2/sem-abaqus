"""Single-grit version of FINAL_RIGID: identical settings, one big grit instead of 712.

Everything is copied from _make_final_rigid.py -- same D50 wheel, same 2.0 mm arc, same
12 um rim, same 30 um width, same 48 x 15 x 6 um workpiece at 0.3 um, same 30 m/s, same
54 um pass, still one discrete rigid body with the workpiece as the only deformable
part. The only change is the grit population.

Two decisions the 712-grit deck did not have to make:

* **Which grain.** Index 0, the largest in the library: 193.9 um^3, 6.54 um tall,
  10.35 um across, aspect ratio 1.05, 116 facets. Chunky and equiaxed rather than a
  sliver, and its 10.35 um footprint still sits inside the 15 um wide block.
* **Where.** At the trailing (high-theta) end of the block footprint, not at its
  centre. With 712 grits, rotation direction did not matter because grits covered the
  whole arc; with one it decides whether the grit traverses the workpiece or leaves it
  immediately. Placed at b = +15 um, a wheel turning toward decreasing theta
  (VR3 < 0) drags the grit 39 um across the 48 um block and out the far edge.
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
from semgrit.wheel import GrainPopulationSpec, WheelSpec, build_wheel
from semgrit.wheel_workpiece import WorkpieceBlock, rotate_placements_about_z

OUT = 'SINGLE_GRIT'
R = 25.0                 # identical to FINAL_RIGID
ARC_MM = 2.0
RIM_MM = 0.012
WHEEL_W = 0.030
SECTOR_DEG = math.degrees(ARC_MM / R)

GRAIN = 0                # largest grain in the library
START_B_MM = 0.015       # grit centre, tangential offset from the block centre

WP = WorkpieceBlock(length_mm=0.048, width_mm=0.015, depth_mm=0.006,
                    element_size_mm=0.0003)

V_SURF = 30_000.0
TRAVEL = WP.length_mm + 0.006

lib = pickle.load(open('WHEEL_FIXED/1_measurements/grain_library.pkl', 'rb'))
solids = lib['solids']
grain = solids[GRAIN]

# --- bond: unchanged --------------------------------------------------------
bond_spec = WheelSpec(diameter_mm=2 * R, width_mm=WHEEL_W, sector_deg=SECTOR_DEG,
                      rim_depth_mm=RIM_MM, radial_divisions=1, axial_divisions=6,
                      circumferential_divisions_per_deg=200.0 / SECTOR_DEG)
bond = build_wheel(bond_spec, solids,
                   GrainPopulationSpec(areal_density_per_mm2=1.0, seed=1))

# --- one grit, seated by the same code path as the 712-grit deck ------------
# A patch sized so the density asks for exactly one grain, with the library cut down
# to the chosen grain so the sampler cannot pick anything else. Going through
# build_wheel rather than hand-building the placement keeps the seating identical:
# tip outward, random spin, tilt up to 35 deg, protrusion drawn from the same
# truncated normal, then the radius solved so the furthest vertex sits exactly that
# far above the bond surface on the *curved* rim.
PATCH = 0.02
seed_spec = WheelSpec(diameter_mm=2 * R, width_mm=PATCH,
                      sector_deg=math.degrees(PATCH / R), rim_depth_mm=RIM_MM,
                      radial_divisions=1, axial_divisions=3,
                      circumferential_divisions_per_deg=100.0)
one = build_wheel(seed_spec, [grain],
                  GrainPopulationSpec(areal_density_per_mm2=1.0 / (PATCH * PATCH),
                                      seed=20260731))
if len(one.placements) != 1:
    raise SystemExit('expected exactly 1 grit, got %d' % len(one.placements))

theta_target = SECTOR_DEG / 2.0 + math.degrees(START_B_MM / R)
placed = rotate_placements_about_z(one.placements,
                                   theta_target - one.placements[0].theta_deg)
p = placed[0]
p.placement_id = 1
p.translation_mm[2] = 0.0        # mid-width; a pure Z shift cannot disturb the seating
p.axial_mm = 0.0

model = dataclasses.replace(bond, placements=[p], shapes=[grain],
                            requested_grains=1, achieved_grains=1)

os.makedirs(OUT, exist_ok=True)
path = os.path.join(OUT, 'wheel_single_grit.inp')
info = write_rigid_wheel_inp(path, model, WP, clearance_um=0.0,
                             engage_window_mm=TRAVEL,
                             model_name='wheel_single_grit')
write_cae_import_script(os.path.join(OUT, 'wheel_single_grit_import_into_cae.py'),
                        path, model_name='wheel_single_grit')

with open(os.path.join(OUT, 'wheel_single_grit_placements.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh)
    wr.writerow(['placement_id', 'shape_index', 'x_mm', 'y_mm', 'z_mm', 'theta_deg',
                 'radius_mm', 'protrusion_um', 'rot_axis_x', 'rot_axis_y',
                 'rot_axis_z', 'rot_angle_deg'])
    wr.writerow([p.placement_id, GRAIN, '%.9f' % p.translation_mm[0],
                 '%.9f' % p.translation_mm[1], '%.9f' % p.translation_mm[2],
                 '%.6f' % p.theta_deg, '%.6f' % p.radius_mm,
                 '%.4f' % (p.protrusion_mm * 1000),
                 '%.6f' % p.rotation_axis[0], '%.6f' % p.rotation_axis[1],
                 '%.6f' % p.rotation_axis[2], '%.4f' % p.rotation_angle_deg])

# --- cost model: identical formula to FINAL_RIGID, for comparability --------
E, nu, rho = WP.youngs_modulus_mpa, WP.poisson_ratio, WP.density_kg_m3 * 1e-12
c = math.sqrt(E * (1 - nu) / ((1 + nu) * (1 - 2 * nu) * rho))
dt = WP.element_size_mm / c
T = TRAVEL / V_SURF
inc = T / dt
el_inc = info['n_workpiece_elements'] * inc
est = {n: el_inc / (3.0e5 * n * 0.70 / 1.6) / 3600.0 for n in (4, 8, 16)}

# Where the grit is, in the frame the workpiece is built in.
theta_c = math.radians(SECTOR_DEG / 2.0)
e_t = (-math.sin(theta_c), math.cos(theta_c), 0.0)
b0 = sum(a * b for a, b in zip(p.translation_mm, e_t))

info.update(
    grain_index=GRAIN, grain_volume_um3=grain.mesh_volume_um3,
    grain_height_um=grain.height_um, grain_footprint_um=max(grain.extent_um()[:2]),
    grain_facets=len(grain.faces), grit_protrusion_um=p.protrusion_mm * 1000,
    grit_tilt_deg=p.rotation_angle_deg, grit_theta_deg=p.theta_deg,
    grit_b_offset_mm=b0, cut_length_mm=b0 + WP.length_mm / 2.0,
    workpiece_element_um=WP.element_size_mm * 1000,
    dilatational_speed_mm_s=c, stable_dt_s=dt, surface_speed_mm_s=V_SURF,
    omega_rad_s=V_SURF / R, travel_mm=TRAVEL, step_time_s=T, increments=inc,
    element_increments=el_inc, est_hours_4core=est[4], est_hours_8core=est[8],
    est_hours_16core=est[16],
)
json.dump(info, open(os.path.join(OUT, 'wheel_single_grit_report.json'), 'w'),
          indent=2, default=str)

print('WHEEL  identical to FINAL_RIGID: R=%.1f mm, %.4f deg = %.3f mm arc, rim %.4f, '
      'width %g' % (R, SECTOR_DEG, ARC_MM, RIM_MM, WHEEL_W))
print('       sagitta %.2f um = %.0f%% of the rim depth' % (
    info['sagitta_um'], 100 * info['sagitta_um'] / 1000 / RIM_MM))
print('       ONE rigid body: %d R3D4 shell quads + %d R3D3 grit facets, ref node %d'
      % (info['n_bond_shell_quads'], info['n_grit_facets'], info['wheel_ref_node']))
print('GRIT   library grain %d: %.1f um3, %.2f um tall, %.2f um across, %d facets'
      % (GRAIN, grain.mesh_volume_um3, grain.height_um,
         max(grain.extent_um()[:2]), len(grain.faces)))
print('       protrusion %.4f um, tilt %.2f deg, at theta %.6f deg (b = %+.4f mm)'
      % (info['grit_protrusion_um'], p.rotation_angle_deg, p.theta_deg, b0))
print('       cut length before it leaves the block: %.4f mm of the %.3f mm block'
      % (info['cut_length_mm'], WP.length_mm))
print('WP     %g x %g x %g mm, %g um elements -> %s C3D8R  (only deformable part)'
      % (WP.length_mm, WP.width_mm, WP.depth_mm, WP.element_size_mm * 1000,
         format(info['n_workpiece_elements'], ',')))
print('       ground face r = %.6f mm, tangent, penetration 0; bond clearance %.4f um'
      % (info['workpiece_ground_radius_mm'], info['max_engaging_protrusion_um']))
print('RUN    dt = %.3e s, omega = %.1f rad/s, travel %.3f mm -> step %.4e s'
      % (dt, V_SURF / R, TRAVEL, T))
print('       %s increments, %.2e element-increments' % (format(int(inc), ','), el_inc))
print('       estimate: %.1f h on 4 cores, %.1f h on 8, %.1f h on 16  (same formula '
      'as the 712-grit deck;' % (est[4], est[8], est[16]))
print('       with 116 contact facets instead of 82,656 it should come in well under)')
print('FILE   %s  (%.2f MB)' % (path, info['size_bytes'] / 1e6))
