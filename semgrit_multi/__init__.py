"""Multi-abrasive extension: the chip thickness computed, not prescribed.

``semgrit`` and ``vumat_grind.for`` are untouched by everything in here. This
package adds two things on top of them:

``envelope``
    A swept-envelope engine. The wheel is a discrete rigid body driven by a
    prescribed velocity boundary condition, so every grit's path over the step
    is known in closed form *before* the run. Sweeping that envelope in time
    order gives the undeformed chip thickness for every element of the
    workpiece, for any number of grits, including the fact that a later grit
    cuts into the groove an earlier one left.

``fieldinject``
    Writes that field into an already-built deck as
    ``*Initial Conditions, type=FIELD, variable=1``. ``vumat_grind.for``
    already reads field variable 1 when PROPS(56) = 1, so nothing in the
    verified subroutine has to change.

``build``
    The two of them plus ``semgrit.build_deck`` in one call.

Why this is better than the four constants the single-grit deck passes: those
constants describe one wedge. Two grits need two, seven hundred need seven
hundred, and none of them can express one grit cutting into another's groove.
The envelope engine has no such limit and needs no new Fortran.
"""

from .envelope import (ChipEnvelope, EnvelopeError, EnvelopeParams,
                       nodal_field, sweep_envelope)
from .fieldinject import InjectError, inject_field

__all__ = ["ChipEnvelope", "EnvelopeError", "EnvelopeParams", "InjectError",
           "inject_field", "nodal_field", "sweep_envelope"]
