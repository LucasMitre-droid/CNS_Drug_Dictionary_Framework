"""
The exposure window, under a single-dose or a compounding half-life.

THE SINGLE-DOSE MODEL

The published framework ends a drug's effective exposure window k half-lives
after the recorded stop date:

    e_EF = e + k * h

with k = 5, at which point 1/32 of the peak remains. This treats the course as
one dose. It is wrong in a specific, one-directional way: it ignores that
repeated dosing accumulates, so a drug given every three weeks for a year does
not start decaying from a single-dose peak on its last day.

THE COMPOUNDING MODEL

Dosing every tau days drives concentration to a plateau

    R = 1 / (1 - (1/2) ** (tau / h))

in units of one single-dose peak. Decay after the last dose therefore starts
from R * D rather than D, and the window becomes

    e_EF = e + (k + log2 R) * h

log2(R) does not depend on k. Compounding is an axis orthogonal to the choice
of k, not a rescaling of it -- which is why the two must be swept separately
rather than folded into a single tuned parameter.

R is a property of the REGIMEN, the drug's half-life measured against its
dosing interval, not of the patient. A per-agent value therefore generalises to
patients whose administration dates were never abstracted. Whether tau exceeds
the half-life decides everything: an antibody with h = 27 d given q20d reaches
R = 2.48 and gains five weeks of window, while an oral agent with h = 0.03 d
given daily reaches R = 1.000 and gains nothing at all. A meaningful fraction of
any real dictionary cannot accumulate under its own schedule.

THE COMBINATION GUARD

A reconstructed cadence belongs to the therapy RECORD, not to each component of
a combination. Handing one row's cadence to every component is safe for
"nivolumab + ipilimumab", where both are three-weekly antibodies, and produces
nonsense for "tucatinib + trastuzumab + capecitabine", where the oral daily
cadence of two components would be inherited by an antibody actually dosed
q21d. Left unguarded that single substitution inflates R by an order of
magnitude and adds months to one agent's window.

So a cadence is attributed to components only when they sit on the same dosing
scale, max(h) / min(h) <= TAU_SHARE_RATIO. Rejected records are returned, never
silently dropped.

THE COVERAGE MODEL

Where administration dates are known, presence can be evaluated pointwise
instead of as a window. Concentration in single-dose-peak units is

    C(t) = SUM_i (1/2) ** ((t - a_i) / h)      over all doses a_i <= t

and the agent is present while C(t) >= (1/2) ** k. For a single dose this is
exactly the single-dose rule, since

    (1/2) ** (delta / h) >= (1/2) ** k   <=>   delta <= k * h

which is what makes a k grid comparable across the two models: at any k they
mean the same thing for a one-shot drug and diverge only through accumulation.
C is a sum of non-negative terms, so superposition presence always contains box
presence -- the model can lengthen exposure and can never shorten it.
"""

from __future__ import annotations

import math
from collections import Counter
from datetime import date, timedelta

# Components of a combination may share a reconstructed cadence only if their
# half-lives sit within this ratio of each other.
TAU_SHARE_RATIO = 3.0

# Relative slack so that an exact tie, delta == k * h, counts as present. The
# two accumulation models must agree exactly at that boundary; without the
# tolerance, floating-point representation decides which side it lands on.
TIE_TOL = 1e-12

BOX = "box"
COMPOUNDING = "compounding"
MODELS = (BOX, COMPOUNDING)


def accumulation_factor(half_life_days: float,
                        tau_days: float | None,
                        n_doses: int | None = None) -> float:
    """
    Steady-state accumulation factor R, in single-dose-peak units.

    R = 1 / (1 - (1/2) ** (tau / h)), capped at what a finite course of
    n_doses can actually reach. A drug given twice has not arrived at the
    plateau of the infinite series and should not be credited with it.

    Returns 1.0 -- no accumulation, identical to the single-dose model --
    whenever the inputs cannot support the calculation.
    """
    if not tau_days or float(tau_days) <= 0:
        return 1.0
    h = float(half_life_days)
    if h <= 0:
        return 1.0
    if n_doses is not None and int(n_doses) < 2:
        return 1.0

    tau = float(tau_days)
    plateau = 1.0 / (1.0 - 0.5 ** (tau / h))
    if n_doses is None:
        return float(plateau)
    reached = sum(0.5 ** ((i * tau) / h) for i in range(int(n_doses)))
    return float(min(plateau, reached))


def extra_half_lives(accumulation: float) -> float:
    """log2 R -- the half-lives compounding adds to the window, independent of k."""
    r = float(accumulation)
    return math.log2(r) if r > 0 else 0.0


def window_half_lives(k: float, accumulation: float = 1.0,
                      model: str = BOX) -> float:
    """Total half-lives from the recorded stop date to the end of the window."""
    _check_model(model)
    if model == BOX:
        return float(k)
    return float(k) + extra_half_lives(accumulation)


def effective_window_end(end_date, half_life_days: float | None, k: float = 5.0,
                         accumulation: float = 1.0, model: str = BOX):
    """
    e_EF, the biologically meaningful end of exposure.

    Returns the recorded end date unchanged when no half-life is available.
    That is the conservative fallback: an agent with an unknown half-life is
    credited with no persistence rather than a guessed one, so it can only ever
    be under-counted as exposure, never over-counted.
    """
    if end_date is None or half_life_days is None:
        return end_date
    h = float(half_life_days)
    if h <= 0:
        return end_date
    return end_date + timedelta(
        days=window_half_lives(k, accumulation, model) * h)


def presence_threshold(k: float) -> float:
    """(1/2) ** k -- the fraction of one single-dose peak that counts as present."""
    return 0.5 ** float(k)


def concentration(admins, t: date, half_life_days: float) -> float:
    """
    C(t), summed over every dose given on or before t, in single-dose peaks.

    Evaluated by direct summation rather than by a day-to-day recurrence. A
    recurrence multiplies by (1/2) ** (1/h) once per day and its rounding
    drifts downward, so at an exact tie -- where the two accumulation models
    must agree by construction -- it lands a hair below threshold and silently
    drops a day of presence. Direct summation is exact there. The cost is
    O(doses) per evaluation, which is negligible at any realistic course length.
    """
    h = float(half_life_days)
    if h <= 0:
        return 0.0
    return sum(0.5 ** ((t - a).days / h) for a in admins if a <= t)


def is_present(admins, t: date, half_life_days: float, k: float,
               model: str = BOX) -> bool:
    """Whether the agent counts as on board on day t."""
    _check_model(model)
    h = float(half_life_days)
    if h <= 0 or not admins:
        return False

    if model == BOX:
        persist = float(k) * h
        return any(a <= t and (t - a).days <= persist for a in admins)

    thr = presence_threshold(k)
    return concentration(admins, t, h) >= thr * (1.0 - TIE_TOL)


def reconstruct_tau(admins) -> float | None:
    """
    Dosing interval from a course's administration dates: the modal gap between
    consecutive doses.

    The mode, not the mean, because a real course contains held doses and
    rescheduled visits. One six-week gap for toxicity should not turn a
    three-weekly antibody into a six-weekly one.
    """
    dates = sorted(set(admins or []))
    if len(dates) < 2:
        return None
    gaps = [(b - a).days for a, b in zip(dates, dates[1:]) if (b - a).days > 0]
    if not gaps:
        return None
    return float(Counter(gaps).most_common(1)[0][0])


def tau_attributable(component_half_lives,
                     ratio: float = TAU_SHARE_RATIO) -> tuple[bool, str]:
    """
    Whether one record's cadence may be attributed to each of its components.

    Returns the verdict and the reason, so a rejection is always explainable in
    the audit output rather than appearing as an unexplained missing value.
    """
    usable = [float(h) for h in (component_half_lives or [])
              if h is not None and float(h) > 0]
    if len(usable) <= 1:
        return True, ("single agent" if len(component_half_lives or []) <= 1
                      else "one component with a known half-life")
    spread = max(usable) / min(usable)
    ok = spread <= float(ratio)
    return ok, (f"combination, half-life spread {spread:.1f}x "
                f"{'<=' if ok else '>'} {float(ratio):g}x")


def _check_model(model: str) -> None:
    if model not in MODELS:
        raise ValueError(
            f"unknown accumulation model {model!r}; expected one of {MODELS}")
