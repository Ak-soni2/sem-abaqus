"""Pass A - is the front half of the pipeline right? Image -> grits -> measurements.

Independent of semgrit wherever it matters: the TIFF tags are parsed from the raw bytes
here, the scale bar is re-measured with plain numpy, and every descriptor is recomputed
from the label mask. Where this file and semgrit agree, two implementations agree.

  A1  calibration    pixel size from the raw TIFF bytes vs what the pipeline used
  A2  scale bar      bar re-measured in its own bbox; bar x pixel size must equal its
                     printed label. This is the check that would have caught the
                     original 15-30x calibration bug on sight.
  A3  grits found    labels partition the foreground, none in the databar, all disjoint
  A4  boundaries     detected edges sit on real image gradient, not on nothing
  A5  measurements   every descriptor recomputed from the mask, in pixels x pixel size
  A6  scale once     with the segmentation held fixed, doubling the pixel size must
                     double every length, quadruple every area, and leave ratios alone
  A7  size chain     measured grain um -> lofted solid um -> grit in the .inp, in mm

Deliberately not verified: the *search* heuristics that locate the databar and the bar.
A2 re-measures the bar inside the bbox the pipeline reports, so it validates the
measurement and its arithmetic, not the hunt for it. The independent guarantee that the
calibration is right comes from A1, which never touches semgrit.

Exits non-zero on any failure.
"""
import dataclasses
import math
import os
import re
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FAIL = []
NOTE = []


def chk(name, ok, detail=''):
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', name,
                           (': ' + detail) if detail else ''))
    if not ok:
        FAIL.append(name)


# --------------------------------------------------------------------------
# A1  Read AP_IMAGE_PIXEL_SIZE straight out of the file, with our own TIFF walk.
# --------------------------------------------------------------------------
def tiff_zeiss_text(path):
    with open(path, 'rb') as fh:
        data = fh.read()
    if data[:2] == b'II':
        end = '<'
    elif data[:2] == b'MM':
        end = '>'
    else:
        return None
    magic, first = struct.unpack(end + 'HI', data[2:8])
    if magic != 42:
        return None
    off, out, seen = first, [], set()
    while off and off not in seen and off + 2 <= len(data):
        seen.add(off)
        n, = struct.unpack(end + 'H', data[off:off + 2])
        for i in range(n):
            p = off + 2 + i * 12
            tag, typ, cnt = struct.unpack(end + 'HHI', data[p:p + 8])
            if tag not in (34118, 34119):
                continue
            size = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 7: 1}.get(typ, 1) * cnt
            if size <= 4:
                blob = data[p + 8:p + 8 + size]
            else:
                vo, = struct.unpack(end + 'I', data[p + 8:p + 12])
                blob = data[vo:vo + size]
            for enc in ('utf-16-le', 'latin-1'):
                out.append(blob.decode(enc, errors='ignore'))
        off, = struct.unpack(end + 'I', data[off + 2 + n * 12:off + 6 + n * 12])
    return '\n'.join(out) if out else None


_UNIT_UM = {'pm': 1e-6, 'nm': 1e-3, 'um': 1.0, '\xb5m': 1.0, 'mm': 1e3, 'm': 1e6}


def independent_pixel_size_um(path):
    txt = tiff_zeiss_text(path)
    if not txt:
        return None
    m = re.search(r'Image Pixel Size\s*=?\s*([0-9.]+)\s*([pnu\xb5m]{1,2})', txt)
    if not m:
        m = re.search(r'AP_IMAGE_PIXEL_SIZE[^0-9-]{0,40}?([0-9.]+)\s*([pnu\xb5m]{1,2})',
                      txt, re.S)
    return float(m.group(1)) * _UNIT_UM.get(m.group(2), 1.0) if m else None


# --------------------------------------------------------------------------
# A2
# --------------------------------------------------------------------------
def bar_length_in_bbox(full, bbox):
    """Longest bright horizontal run inside the reported bar bbox."""
    # Both axes must be cropped. Taking the full row instead lets the databar's own
    # wide bright band satisfy an "overlaps the bar's x-range" test and the bar comes
    # back 20x too long -- which is precisely the failure mode this file exists to
    # catch, so it must not be reproduced by the check itself.
    y0, x0, y1, x1 = bbox
    sub = full[y0:y1 + 1, x0:x1 + 1]
    thr = 0.5 * (float(sub.max()) + float(sub.min()))
    best = 0
    for r in range(sub.shape[0]):
        idx = np.flatnonzero(sub[r] >= thr)
        if idx.size < 5:
            continue
        brk = np.flatnonzero(np.diff(idx) > 1)
        for a, b in zip(idx[np.r_[0, brk + 1]], idx[np.r_[brk, idx.size - 1]]):
            best = max(best, int(b - a + 1))
    return best


def snap_125(v):
    """Nearest 1-2-5 value, the only lengths a microscope prints on a scale bar."""
    if v <= 0:
        return v
    e = math.floor(math.log10(v))
    m = v / 10.0 ** e
    return min((1.0, 2.0, 5.0, 10.0), key=lambda c: abs(c - m)) * 10.0 ** e


def max_caliper(pts):
    """Largest distance between any two points -- rotation invariant."""
    p = np.asarray(pts, dtype=np.float64)
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=-1)
    return float(d.max())


def max_caliper_2d(pts, n=180):
    p = np.asarray(pts, dtype=np.float64)[:, :2]
    best = 0.0
    for a in np.linspace(0.0, math.pi, n, endpoint=False):
        proj = p @ np.array([math.cos(a), math.sin(a)])
        best = max(best, float(proj.max() - proj.min()))
    return best


# --------------------------------------------------------------------------
def main():
    from semgrit.build_deck import DeckParams, build_deck
    from semgrit.grain3d import HeightModel, LoftProfile, build_grain_library
    from semgrit.measure import measure_all
    from semgrit.metrology import load_sem_image
    from semgrit.segment import SegmentationParams, segment_grains

    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    quick = '--quick' in sys.argv[1:]
    images = args or sorted(f for f in os.listdir('.')
                            if f.lower().endswith(('.tif', '.tiff')))
    images = [f for f in images if os.path.exists(f)]
    if not images:
        print('no images given and none found in this directory')
        return 2

    print('=' * 78)
    print('PASS A   image -> grits -> measurements    (%d images)' % len(images))
    print('=' * 78)
    print('A1/A2  CALIBRATION  (raw-TIFF parse, and bar x pixel size vs its label)')
    print('  %-15s %10s %10s %8s %7s %8s %8s %7s'
          % ('image', 'raw TIFF', 'pipeline', 'match', 'bbox px', 'bar->um', 'label',
             'err'))
    seg_params = SegmentationParams()
    sems, bad_cal, bad_len, bad_label = {}, 0, 0, 0
    for f in images:
        sem = load_sem_image(f)
        sems[f] = sem
        mine = independent_pixel_size_um(f)
        ok_cal = mine is not None and abs(mine - sem.pixel_size_um) <= 1e-12
        bad_cal += (not ok_cal)
        b = sem.scale_bar
        if b is None:
            print('  %-15s %10.6f %10.6f %8s %7s %8s %8s %7s'
                  % (f[:15], mine or float('nan'), sem.pixel_size_um,
                     'yes' if ok_cal else 'NO', '-', '-', '-', '-'))
            continue
        # The bbox is the bar's *outer* extent; the pipeline measures tick centre to
        # tick centre, which is what the printed label refers to, so it is legitimately
        # a couple of pixels shorter. Require it to be essentially the bar and no wider
        # than the bar's own bounding box -- the label arithmetic below is what decides
        # whether the convention is the right one.
        mine_px = bar_length_in_bbox(sem.full_intensity, b.bbox)
        if not (0.80 * mine_px <= b.length_px <= mine_px + 1e-9):
            bad_len += 1
        implied = b.length_px * sem.pixel_size_um     # bar length in metadata microns
        snapped = snap_125(implied)
        err = abs(implied - snapped) / snapped
        if err > 0.02 or abs(snapped - b.snapped_label_um) > 1e-9:
            bad_label += 1
        print('  %-15s %10.6f %10.6f %8s %7d %8.4f %8.3f %6.2f%%'
              % (f[:15], mine or float('nan'), sem.pixel_size_um,
                 'yes' if ok_cal else 'NO', mine_px, implied, snapped, 100 * err))

    chk('pixel size matches an independent raw-TIFF parse on every image',
        bad_cal == 0, '%d of %d disagree' % (bad_cal, len(images)))
    chk('bar length is consistent with its own bbox (tick centres vs outer edge)',
        bad_len == 0, '%d of %d outside 80-100%% of the bbox width' % (bad_len, len(images)))
    chk('bar length x metadata pixel size equals the printed 1-2-5 label',
        bad_label == 0, '%d of %d inconsistent' % (bad_label, len(images)))
    ag = [s.scalebar_agreement for s in sems.values()
          if s.scalebar_agreement is not None]
    # scalebar_agreement is a *relative difference*, so 0 is perfect, not 1.
    chk('metadata and drawn bar agree within 5% on every image',
        bool(ag) and max(abs(a) for a in ag) <= 0.05,
        'worst %.2f%% over %d images' % (100 * max(abs(a) for a in ag), len(ag)))
    chk('no image silently fell back to the bar for calibration',
        all(s.pixel_size_source == 'metadata' for s in sems.values()),
        ' '.join(sorted({s.pixel_size_source for s in sems.values()})))
    labels = sorted({s.scale_bar.snapped_label_um for s in sems.values()
                     if s.scale_bar})
    if len(sems) > 1:
        # Only meaningful over a set: it demonstrates the label is read per image
        # rather than hardcoded. A single image has exactly one label by definition,
        # so asserting variety there would fail on a perfectly good picture.
        chk('bar labels are read per image, not assumed constant', len(labels) > 1,
            'labels found: %s um  (a hardcoded 2.0 would misread the 1 um ones)'
            % labels)
    else:
        print('  [info] bar label read from this image: %s um (variety is only '
              'checkable over a set)' % labels)

    # ------------------------------------------------------------------
    subset = [s for s in ('DIAMOND_11.tif', 'B4C_15.tif', 'DIAMOND_14.tif')
              if s in sems] or images[:2]
    print()
    print('A3/A4/A5  GRITS AND MEASUREMENTS')
    segs = {}
    for f in subset:
        sem = sems[f]
        seg = segment_grains(sem, seg_params)
        segs[f] = seg
        grains = measure_all(seg, sem)
        lab = seg.labels
        px = sem.pixel_size_um
        print('  %s: %d grains, %d seeds, %.6f um/px' % (f, len(grains), seg.n_seeds, px))

        ids = np.unique(lab)
        ids = ids[ids > 0]
        chk('%s: one region per label, all disjoint' % f,
            len(ids) == len(grains) and int((lab > 0).sum()) ==
            sum(int((lab == i).sum()) for i in ids),
            '%d labels, %d grains' % (len(ids), len(grains)))
        chk('%s: segmentation covers the micrograph only, not the databar' % f,
            lab.shape == sem.intensity.shape,
            'labels %s, micrograph %s, databar cut at row %d'
            % (lab.shape, sem.intensity.shape, sem.databar_top))

        gy, gx = np.gradient(sem.intensity.astype(np.float64))
        gmag = np.hypot(gx, gy)
        inner = np.zeros(lab.shape, bool)
        edge = np.zeros(lab.shape, bool)
        for i in ids:
            m = lab == i
            er = m.copy()
            er[1:, :] &= m[:-1, :]
            er[:-1, :] &= m[1:, :]
            er[:, 1:] &= m[:, :-1]
            er[:, :-1] &= m[:, 1:]
            edge |= m & ~er
            inner |= er
        ratio = (gmag[edge].mean() / gmag[inner].mean()) if inner.any() else 0.0
        chk('%s: grain boundaries sit on real image gradient' % f, ratio > 2.0,
            'boundary/interior gradient = %.2f' % ratio)

        worst = {}
        for g in grains[:60]:
            m = lab == g.label
            n = int(m.sum())
            ys, xs = np.nonzero(m)
            rec = {'area_um2': n * px * px,
                   'equivalent_diameter_um': 2.0 * math.sqrt(n * px * px / math.pi),
                   'bbox_width_um': (xs.max() - xs.min() + 1) * px,
                   'bbox_height_um': (ys.max() - ys.min() + 1) * px,
                   'centroid_x_um': xs.mean() * px,
                   'centroid_y_um': ys.mean() * px}
            for k, v in rec.items():
                got = getattr(g, k, None)
                if got is not None and v != 0:
                    worst[k] = max(worst.get(k, 0.0), abs(got - v) / abs(v))
        for k, v in sorted(worst.items()):
            chk('%s: %s recomputed from the mask' % (f, k), v < 1e-9,
                'worst relative difference %.2e' % v)

        wf = 0.0
        for g in grains[:40]:
            ys, xs = np.nonzero(lab == g.label)
            wf = max(wf, abs(g.feret_max_um - max_caliper_2d(
                np.column_stack([xs, ys]).astype(float) * px)) / g.feret_max_um)
        chk('%s: max Feret recomputed by brute-force projection' % f, wf < 0.02,
            'worst relative difference %.2e over 40 grains' % wf)

    # ------------------------------------------------------------------
    print()
    print('A6  SCALE APPLIED EXACTLY ONCE  (same segmentation, pixel size x2)')
    f = subset[0]
    sem = sems[f]
    seg = segs[f]
    K = 2.0
    ga = measure_all(seg, sem)
    gb = measure_all(seg, dataclasses.replace(sem, pixel_size_um=sem.pixel_size_um * K))
    chk('the same grains are measured', len(ga) == len(gb),
        '%d vs %d' % (len(ga), len(gb)))
    groups = {1: ['equivalent_diameter_um', 'feret_max_um', 'feret_min_um',
                  'perimeter_um', 'bbox_width_um', 'bbox_height_um',
                  'centroid_x_um', 'centroid_y_um'],
              2: ['area_um2', 'convex_area_um2'],
              0: ['solidity', 'circularity', 'aspect_ratio', 'convexity']}
    for power, names in groups.items():
        worst, which = 0.0, ''
        for x, y in zip(ga, gb):
            for k in names:
                u, v = getattr(x, k, None), getattr(y, k, None)
                if not u or v is None:
                    continue
                e = abs(v - u * K ** power) / abs(u * K ** power)
                if e > worst:
                    worst, which = e, k
        chk('%s scale as pixel_size**%d' % (
            {1: 'lengths', 2: 'areas', 0: 'dimensionless ratios'}[power], power),
            worst < 1e-12, 'worst %.2e (%s)' % (worst, which))

    # ------------------------------------------------------------------
    print()
    print('A7  SIZE CHAIN  measured grain um -> lofted solid um -> grit in the .inp mm')
    solids, _ = build_grain_library(grains_for := measure_all(segs[f], sems[f]),
                                    segs[f], sems[f],
                                    height_model=HeightModel(seed=1),
                                    profile=LoftProfile(), simplify_um=0.10,
                                    max_vertices=64, interior_only=True)
    by_id = {g.grain_id: g for g in grains_for}
    worst_w, ratios = 0.0, []
    for s in solids:
        g = by_id.get(s.grain_id)
        if g is None:
            continue
        # The full measured outline sits at 42% height (scale 1.0), so the widest
        # cross-section of the solid is the measured outline itself.
        worst_w = max(worst_w, abs(max_caliper_2d(s.vertices - s.centroid_um)
                                   - g.feret_max_um) / g.feret_max_um)
        ratios.append(s.height_um / g.feret_min_um)
    chk('the solid keeps the measured grain width', worst_w < 0.08,
        'worst %.2f%% over %d solids (outline simplified at 0.10 um)'
        % (100 * worst_w, len(solids)))
    r = np.array(ratios)
    chk('modelled height stays inside the declared 0.45-0.95 x minFeret band',
        bool((r >= 0.45 - 1e-9).all() and (r <= 0.95 + 1e-9).all()),
        'range %.3f..%.3f' % (r.min(), r.max()))
    chk('modelled height averages the declared 0.70 x minFeret',
        abs(r.mean() - 0.70) < 0.03, 'mean %.4f over %d solids' % (r.mean(), len(r)))
    NOTE.append('height is the one modelled quantity: an SEM gives no depth, so it is '
                '%.2f x min Feret with scatter, not a measurement' % 0.70)

    if quick:
        print('  (--quick: skipping the deck build; the deck verifiers cover that)')
        print()
        for n in NOTE:
            print('  note: %s' % n)
        print()
        print('=' * 78)
        print('PASS A TOTAL: %d failure(s)%s'
              % (len(FAIL), '' if not FAIL else '  -> ' + str(FAIL)))
        print('=' * 78)
        return 1 if FAIL else 0

    idx = int(np.argmax([s.mesh_volume_um3 for s in solids]))
    info = build_deck(DeckParams(
        name='chain', sector_mode='arc', arc_length_mm=2.0, grit_mode='single',
        single_grain_index=idx, wp_element_size_mm=0.002, seed=3), solids, '_pipeA')
    txt = open(info['path'], encoding='ascii').read()
    part = txt.split('*Part, name=WHEEL')[1].split('*End Part')[0]
    nodes = {}
    for ln in part.split('*Element')[0].strip().split('\n'):
        v = ln.split(',')
        if len(v) >= 4:
            nodes[int(v[0])] = [float(x) for x in v[1:4]]
    tri_block = part.split('*Element, type=R3D3')[1].split('*Elset')[0]
    gids = sorted({int(x) for ln in tri_block.strip().split('\n')
                   for x in ln.split(',')[1:4]})
    deck_um = max_caliper([nodes[i] for i in gids]) * 1000.0
    lib_um = max_caliper(solids[idx].vertices)
    chk('the grit in the .inp is the measured grain, size preserved into mm',
        abs(deck_um - lib_um) / lib_um < 1e-6,
        'deck %.6f um vs library %.6f um (max caliper, rotation invariant)'
        % (deck_um, lib_um))
    chk('microns became millimetres exactly once (no 1000x slip)',
        0.001 < lib_um / 1000.0 < 0.1,
        'largest grain %.3f um = %.6f mm' % (lib_um, lib_um / 1000.0))

    print()
    for n in NOTE:
        print('  note: %s' % n)
    print()
    print('=' * 78)
    print('PASS A TOTAL: %d failure(s)%s'
          % (len(FAIL), '' if not FAIL else '  -> ' + str(FAIL)))
    print('=' * 78)
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
