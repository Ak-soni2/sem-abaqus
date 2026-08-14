"""The settings a viewer may edit, and the one function that applies them.

The browser is allowed to collect numbers. It is not allowed to decide anything. Every
edit made in the CAD viewer comes back here as a flat dict of named scalars, and this
module turns it into a :class:`DeckParams` -- validated, and by exactly the same path
whether it arrived through a Colab kernel callback, a pasted JSON string or a downloaded
``viewer_settings.json``.

That matters more than it sounds. The seating of the block is decided by
``rigid_wheel.ground_radius``, which clips every grit facet to the block footprint; the
engaging set is decided by ``rigid_wheel.engaging_grits``. Neither can be reproduced in
JavaScript, and a project that has twice shipped two implementations agreeing on the same
wrong answer should not try. So the viewer previews an edit by moving a box, and says so;
the numbers that depend on seating are recomputed here and only here.

Two tiers, because they cost differently:

``live``
    Moves or resizes geometry that already exists. The viewer can show it immediately by
    transforming the drawn box, and nothing has to be re-placed.
``rebuild``
    Changes which grains exist or where they sit, so the grain packing has to run again.
    ``plan_deck`` takes 0.13 s at 60 grains and 6.6 s at 4,000, which is why these are
    committed deliberately rather than on every keystroke.
"""

from __future__ import annotations

import dataclasses
import json
from typing import Optional

from .analysis import AnalysisParams
from .build_deck import DeckError, DeckParams
from .rigid_wheel import WP_POSITIONS


@dataclasses.dataclass(frozen=True)
class Field:
    """One editable quantity: where it lives, what it means, how far it may go.

    ``widget`` is the notebook form variable this corresponds to, and ``widget_scale``
    the factor from this field's unit to the widget's. Both live here, next to the field,
    because the mapping has two traps that a separate lookup table would get wrong:
    ``width_mm`` is called ``WHEEL_WIDTH_MM``, and ``surface_speed_mm_s`` is presented as
    ``SURFACE_SPEED_M_S`` in metres per second.
    """

    key: str
    label: str
    unit: str
    tier: str                     # 'live' | 'rebuild'
    lo: Optional[float] = None
    hi: Optional[float] = None
    on_analysis: bool = False     # lives on AnalysisParams, not DeckParams
    choices: tuple = ()
    note: str = ""
    widget: str = ""
    widget_scale: float = 1.0
    seating: bool = False
    """True if changing this re-seats the block, so anything derived from the seating
    must be shown as stale until Python has recomputed it."""


# Ordered as a user reads them: the block first, because that is what they came to move.
FIELDS = (
    Field("wp_length_mm", "block length (along the cut)", "mm", "live", 1e-4, 50.0,
          widget="WP_LENGTH_MM", seating=True),
    Field("wp_width_mm", "block width (across the face)", "mm", "live", 1e-4, 50.0,
          widget="WP_WIDTH_MM", seating=True),
    Field("wp_depth_mm", "block depth", "mm", "live", 1e-4, 50.0,
          widget="WP_DEPTH_MM"),
    Field("wp_position", "where on the arc", "", "rebuild", choices=WP_POSITIONS,
          note="dragging the block sets this to 'custom angle'",
          widget="WP_POSITION", seating=True),
    Field("wp_position_deg", "custom angle", "deg", "live", -360.0, 360.0,
          note="only used by 'custom angle'; dragging writes this",
          widget="WP_POSITION_DEG", seating=True),
    Field("clearance_um", "standoff", "um", "live", 0.0, 1000.0,
          note="0 = the tallest grain under the block is tangent to it. Lifting it "
               "shifts the whole depth-of-cut window by the same amount",
          widget="CLEARANCE_UM"),
    Field("depth_of_cut_um", "depth of cut", "um", "live", 0.0, 1000.0,
          on_analysis=True,
          note="must exceed the first-contact infeed and stay under the "
               "face-to-bond gap",
          widget="DEPTH_OF_CUT_UM"),
    Field("surface_speed_mm_s", "wheel surface speed", "mm/s", "live", 1.0, 200_000.0,
          widget="SURFACE_SPEED_M_S", widget_scale=1e-3),
    Field("rotation_reversed", "reverse the rotation", "", "live",
          note="flips which end of the arc the grains arrive from",
          widget="ROTATION_REVERSED", on_analysis=True, seating=True),
    Field("diameter_mm", "wheel diameter", "mm", "rebuild", 0.1, 2000.0,
          widget="DIAMETER_MM", seating=True),
    Field("arc_length_mm", "slice arc length", "mm", "rebuild", 1e-4, 5000.0,
          note="in force when the wheel extent is given as an arc",
          widget="ARC_LENGTH_MM", seating=True),
    Field("sector_deg", "sector angle", "deg", "rebuild", 0.01, 360.0,
          note="in force when the wheel extent is given as an angle; editing it "
               "switches the extent to 'angle'",
          widget="SECTOR_DEG", seating=True),
    Field("rim_depth_mm", "rim depth", "mm", "rebuild", 1e-4, 500.0,
          widget="RIM_DEPTH_MM"),
    Field("width_mm", "wheel face width", "mm", "rebuild", 1e-4, 500.0,
          widget="WHEEL_WIDTH_MM", seating=True),
    Field("grit_mode", "how the grits are counted", "", "rebuild",
          choices=("concentration", "areal_density", "count", "single"),
          widget="GRIT_MODE", seating=True),
    Field("grit_count", "number of grits", "", "rebuild", 1, 500_000,
          widget="GRIT_COUNT", seating=True),
    Field("concentration", "concentration (C-number)", "", "rebuild", 1.0, 400.0,
          widget="CONCENTRATION", seating=True),
    Field("areal_density_per_mm2", "grits per mm2", "1/mm2", "rebuild", 1.0, 200_000.0,
          widget="AREAL_DENSITY_PER_MM2", seating=True),
    Field("friction", "friction coefficient", "", "live", 0.0, 2.0, on_analysis=True,
          widget="FRICTION"),
)
BY_KEY = {f.key: f for f in FIELDS}
LIVE = tuple(f.key for f in FIELDS if f.tier == "live")
REBUILD = tuple(f.key for f in FIELDS if f.tier == "rebuild")

# Not a DeckParams field: the deck expresses the sense through the sign of VR3, which
# analysis.wheel_motion derives. Reversing is offered as a flag and applied there.
SYNTHETIC = ("rotation_reversed",)


def settings_from_params(p: DeckParams) -> dict:
    """The editable subset of a DeckParams, as a flat JSON-safe dict."""
    an = p.analysis
    out = {}
    for f in FIELDS:
        if f.key in SYNTHETIC:
            out[f.key] = bool(getattr(an, "rotation_reversed", False)) if an else False
        elif f.on_analysis:
            out[f.key] = getattr(an, f.key) if an is not None else None
        else:
            out[f.key] = getattr(p, f.key)
    # sector_mode is not editable directly -- editing the angle implies 'angle' -- but
    # the viewer needs to know which one is in force to label the field honestly.
    #
    # This is why arc_length_mm had to become editable: on an arc-mode wheel the raw
    # sector_deg is whatever the dataclass default happens to be, so a panel showing it
    # reported 30 degrees for a 2.29 degree arc. The resolved angle is reported beside
    # it so the two can never be confused again.
    out["_sector_mode"] = p.sector_mode
    out["_sector_resolved_deg"] = p.resolved_sector_deg()
    out["_arc_resolved_mm"] = p.outer_radius_mm * (
        p.resolved_sector_deg() * 3.141592653589793 / 180.0)
    return out


def param_block(settings: dict, base: DeckParams, only_changed: bool = True) -> str:
    """The edited settings as notebook form assignments, ready to paste.

    The form widgets cannot be written to from Python, so after an edit is applied the
    widgets still show their old values. This closes that gap by hand: it prints exactly
    the lines to paste back, in the widgets' own names and units -- including the two
    that differ (``WHEEL_WIDTH_MM``, and ``SURFACE_SPEED_M_S`` in m/s).
    """
    before = settings_from_params(base)
    lines = []
    for f in FIELDS:
        if not f.widget or f.key not in settings:
            continue
        v = settings[f.key]
        if only_changed and f.key in before and before[f.key] == v:
            continue
        if isinstance(v, bool):
            txt = "True" if v else "False"
        elif f.choices:
            txt = '"%s"' % v
        else:
            scaled = float(v) * f.widget_scale
            txt = ("%d" % round(scaled)) if f.key == "grit_count" else "%g" % scaled
        lines.append("%s = %s" % (f.widget, txt))
    # Editing the angle or the arc also decides which one is in force.
    if any(l.startswith("SECTOR_DEG ") for l in lines):
        lines.append('SECTOR_MODE = "angle"')
    elif any(l.startswith("ARC_LENGTH_MM ") for l in lines):
        lines.append('SECTOR_MODE = "arc"')
    return "\n".join(lines)


def params_from_settings(settings: dict, base: DeckParams) -> DeckParams:
    """Apply an edited settings dict to ``base``. The only way edits become a deck.

    Unknown keys are refused rather than ignored: a typo that silently did nothing would
    be worse than an error, because the viewer would show one model and the deck would
    contain another.
    """
    if not isinstance(settings, dict):
        raise DeckError("settings must be a JSON object of named values")
    junk = [k for k in settings
            if not k.startswith("_") and k not in BY_KEY]
    if junk:
        raise DeckError("settings contains fields that are not editable: %s"
                        % ", ".join(sorted(junk)))

    # Only *changes* are applied. Re-submitting the settings unchanged has to be a
    # no-op, and it was not: every settings dict carries sector_deg, and applying it
    # switched sector_mode from 'arc' to 'angle' -- quietly redefining the wheel on an
    # edit the user never made.
    before = settings_from_params(base)
    deck_kw, an_kw = {}, {}
    for k, v in settings.items():
        if k.startswith("_"):
            continue
        f = BY_KEY[k]
        if k in before and before[k] == v:
            continue
        if f.choices:
            if v not in f.choices:
                raise DeckError("%s must be one of %s, not %r"
                                % (k, ", ".join(map(str, f.choices)), v))
        elif isinstance(v, bool):
            pass
        elif v is not None:
            try:
                v = float(v)
            except (TypeError, ValueError):
                raise DeckError("%s must be a number, not %r" % (k, v))
            if f.lo is not None and v < f.lo:
                raise DeckError("%s = %g is below the %g %s minimum"
                                % (k, v, f.lo, f.unit))
            if f.hi is not None and v > f.hi:
                raise DeckError("%s = %g is above the %g %s maximum"
                                % (k, v, f.hi, f.unit))
            if k in ("grit_count",):
                v = int(round(v))
        (an_kw if (f.on_analysis or k in SYNTHETIC) else deck_kw)[k] = v

    # The wheel's extent is given one way or the other, never both. Editing the angle
    # puts it in 'angle' mode and editing the arc puts it in 'arc' mode; asking for both
    # at once has no single answer, so it is refused rather than resolved by whichever
    # branch happens to run last.
    if "sector_deg" in deck_kw and "arc_length_mm" in deck_kw:
        raise DeckError(
            "the wheel extent can be set as an angle or as an arc length, not both: "
            "sector_deg = %g and arc_length_mm = %g were both changed. Edit one."
            % (deck_kw["sector_deg"], deck_kw["arc_length_mm"]))
    if "sector_deg" in deck_kw:
        deck_kw["sector_mode"] = "angle"
    elif "arc_length_mm" in deck_kw:
        deck_kw["sector_mode"] = "arc"

    an = base.analysis
    if an_kw:
        if an is None:
            raise DeckError("this deck has no analysis, so %s cannot be set. Enable the "
                            "run-ready analysis first." % ", ".join(sorted(an_kw)))
        an = dataclasses.replace(an, **an_kw)
    out = dataclasses.replace(base, analysis=an, **deck_kw)
    validate_edit(out)
    return out


def validate_edit(p: DeckParams) -> None:
    """Geometric sanity an edit can break, checked before anything is written.

    These are the states a viewer makes easy to reach by dragging: a block longer than
    the slice it sits on, or wider than the face it runs across.
    """
    p.validate()
    if not p.include_workpiece:
        return
    arc = p.outer_radius_mm * (p.resolved_sector_deg() * 3.141592653589793 / 180.0)
    if p.wp_length_mm > arc:
        raise DeckError(
            "the block is %.4f mm long but the wheel slice is only %.4f mm of arc, so "
            "it would hang off both ends. Lengthen the slice or shorten the block."
            % (p.wp_length_mm, arc))
    if p.wp_width_mm > p.width_mm:
        raise DeckError(
            "the block is %.4f mm wide but the wheel face is only %.4f mm, so it would "
            "overhang the sides." % (p.wp_width_mm, p.width_mm))


def apply(settings: dict, base: DeckParams, solids) -> dict:
    """Edited settings -> a re-planned model. What both commit paths call.

    Returns ``{"params", "plan", "settings", "changed", "tier"}``. ``tier`` is 'live' if
    nothing that re-places grains was touched, which is the caller's cue that it may
    reuse a cached measurement rather than re-measuring.
    """
    from .build_deck import plan_deck

    before = settings_from_params(base)
    p = params_from_settings(settings, base)
    after = settings_from_params(p)
    changed = [k for k in after if not k.startswith("_") and after[k] != before[k]]
    tier = "rebuild" if any(k in REBUILD for k in changed) else "live"
    return {"params": p, "plan": plan_deck(p, solids),
            "settings": after, "changed": changed, "tier": tier}


class CommitReply:
    """What a viewer Apply callback must return, in a form the browser can read.

    Colab does not hand a callback's return value to the page as-is. It runs it through
    IPython's display formatter and the page receives the resulting mimetype bundle. A
    plain ``dict`` formats to ``text/plain`` **only** -- a Python repr, not JSON -- so a
    viewer reading ``application/json`` saw no reply at all. That is not a cosmetic bug:
    the viewer then had nothing to distinguish success from failure, and reported every
    Apply as a refusal with no reason, including the ones that had already been applied.

    So this carries the payload under *both* mimetypes, exactly:

    ``application/json``
        via ``_repr_json_``, which is what the formatter asks for.
    ``text/plain``
        via ``__repr__``, as :data:`SENTINEL` followed by the same payload as JSON text,
        so a runtime that only forwards plain text is still read exactly rather than
        guessed at.

    Both are the same ``json.dumps``-able dict, so the two can never disagree.
    """

    SENTINEL = "CADREPLY "

    def __init__(self, **payload):
        self.payload = dict(payload)

    def _repr_json_(self):
        return dict(self.payload)

    def __repr__(self):
        return self.SENTINEL + json.dumps(self.payload, sort_keys=True, default=str)

    def __getitem__(self, key):
        return self.payload[key]

    def get(self, key, default=None):
        return self.payload.get(key, default)


def commit_reply(ok: bool, message: str = "", error: str = "") -> CommitReply:
    """Build the Apply reply. ``error`` is never allowed to be empty on a refusal.

    A refusal the user cannot act on is barely better than a silent one, so if something
    raised without a message this says that, rather than leaving the viewer to print
    'no reason given'.
    """
    if not ok and not error:
        error = ("Python refused the edit but raised no message. Re-run the rebuild cell "
                 "with these settings to see the full traceback.")
    return CommitReply(ok=bool(ok), message=message, error=error)


def load(path_or_text: str) -> dict:
    """Read a settings dict from a file path or from pasted JSON text."""
    txt = path_or_text
    try:
        import os
        if os.path.exists(path_or_text):
            with open(path_or_text, encoding="utf-8") as fh:
                txt = fh.read()
    except (OSError, ValueError):
        pass
    try:
        got = json.loads(txt)
    except ValueError as exc:
        raise DeckError("that is not valid JSON: %s" % exc)
    # The viewer exports {"settings": {...}, "meta": {...}}; a bare dict is fine too.
    return got.get("settings", got) if isinstance(got, dict) else got
