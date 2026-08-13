"""
The two concurrency rules, and the split of adjuvant exposure.

Two defensible definitions of "concurrent" applied to the same cohort will
disagree, and the disagreement is not small. They answer different questions:

  WINDOW    Did the drug's effective exposure window overlap the radiotherapy
            interval at all? Needs only a therapy start and stop date, so it
            covers every patient in a dataset. Resolves time to the course,
            not the dose.

  COVERAGE  Was the drug on board for at least half the days of radiotherapy?
            Needs reconstructed administration dates, so it covers only the
            subset of patients whose records were abstracted to that depth.
            Resolves time to the dose.

WINDOW is the published rule. COVERAGE is stricter by construction: a drug
given once in the first week of a six-week course satisfies WINDOW and fails
COVERAGE. Neither is a better-powered version of the other and their outputs
are not interchangeable. Report which one produced a result, and prefer to
report both.

Both rules take k, the number of half-lives an agent keeps counting as present.
Both accept either accumulation model. That gives a grid, and the grid is the
point: an exposure effect that appears at one cell and vanishes at the
neighbouring one is a property of the definition, not of the biology. Read the
distribution across the grid, never the cell with the smallest p-value.

ADJUVANT EXPOSURE HAS TWO SOURCES, AND THEY ARE NOT THE SAME CLAIM

  persistence  the drug was stopped on or before the last fraction, but its
               effective window still extends past it. The patient took nothing
               new; pharmacokinetics carried the exposure forward.

  initiation   the drug was started after the last fraction. This is a
               clinical decision to begin adjuvant therapy.

Collapsing them loses the distinction between continued exposure and true
adjuvant treatment, which is usually the distinction an adjuvant-therapy
analysis is trying to make. They are returned separately and their union is
returned as well.

A NOTE ON tau

The published framework writes tau for the post-radiotherapy attribution
window. The compounding model writes tau for a drug's dosing interval. They are
unrelated quantities. This module calls the first `attribution_window_days` and
the second never appears here; it lives in `pharmacokinetics`.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, timedelta

from . import pharmacokinetics as pk
from . import scoring

WINDOW = "window"
COVERAGE = "coverage"
RULES = (WINDOW, COVERAGE)

DEFAULT_COVERAGE_THRESHOLD_PCT = 50.0


@dataclass
class Timing:
    """
    One agent's classification against one radiotherapy interval.

    `unresolved` marks a record the rules could not decide -- a missing stop
    date, an unreadable radiotherapy interval. Such a record must not be
    reported as unexposed. Exposure variables are aggregated by disjunction, so
    a False that really means "unknown" is absorbed silently by any other
    record and biases exposure downward with no trace in the output. Unresolved
    records surface as missing values and are counted separately.
    """
    cns_active: bool
    concurrent: bool
    adjuvant_persistence: bool
    adjuvant_initiation: bool
    category: str
    effective_window_end: date | None = None
    covered_days: int | None = None
    window_days: int | None = None
    coverage_pct: float | None = None
    flag: str = ""
    unresolved: bool = False

    @property
    def adjuvant(self) -> bool:
        return self.adjuvant_persistence or self.adjuvant_initiation


def classify_window(cns_class, start, end, rt_start: date, rt_end: date,
                    half_life_days: float | None, k: float = 5.0,
                    accumulation: float = 1.0, model: str = pk.BOX,
                    activity_threshold: int = scoring.DEFAULT_ACTIVITY_THRESHOLD,
                    attribution_window_days: float | None = None) -> Timing:
    """
    The overlap rule.

        unrelated    e_EF < t0
        concurrent   active and s <= t1 and e_EF >= t0
        persistence  active and s <= t1 and e_EF >  t1
        initiation   active and s >  t1 and s <= t1 + attribution window

    `attribution_window_days = None` leaves initiation unbounded: any CNS-active
    therapy started after the last fraction counts, however late. Bound it when
    the analysis needs adjuvant exposure attributable to the radiotherapy
    course rather than to the patient's subsequent disease history.
    """
    active = scoring.is_cns_active(cns_class, activity_threshold)

    if end is None:
        return Timing(active, False, False, False, "manual_review",
                      None, flag="missing_end_date", unresolved=True)

    e_ef = pk.effective_window_end(end, half_life_days, k, accumulation, model)
    flag = "" if half_life_days is not None else "no_half_life_window_not_extended"

    if e_ef is not None and e_ef < rt_start:
        return Timing(active, False, False, False, "unrelated", e_ef, flag=flag)

    concurrent = bool(active and start is not None
                      and start <= rt_end and e_ef >= rt_start)
    persistence = bool(active and start is not None
                       and start <= rt_end and e_ef > rt_end)
    initiation = bool(active and start is not None and start > rt_end)
    if initiation and attribution_window_days is not None:
        initiation = start <= rt_end + timedelta(days=float(attribution_window_days))

    return Timing(active, concurrent, persistence, initiation,
                  _category(active, concurrent, persistence or initiation),
                  e_ef, flag=flag)


def classify_coverage(cns_class, admins, rt_start: date, rt_end: date,
                      half_life_days: float | None, k: float = 5.0,
                      model: str = pk.BOX,
                      threshold_pct: float = DEFAULT_COVERAGE_THRESHOLD_PCT,
                      activity_threshold: int = scoring.DEFAULT_ACTIVITY_THRESHOLD
                      ) -> Timing:
    """
    The coverage rule: concurrent when the agent is on board on at least
    `threshold_pct` of radiotherapy days.

    Presence is evaluated once per calendar day, at the day's maximum. That is
    a deliberate coarsening. At true clock resolution an oral agent with a
    45-minute half-life given daily is present about 9% of the time, which
    would make every short-half-life oral drug non-concurrent by construction
    and turn the rule into a proxy for route of administration.

    Adjuvant exposure is not derivable from this rule -- it reads only the
    radiotherapy interval -- so use `classify_window` for the adjuvant terms.
    Those two questions are independent, and gating adjuvant persistence on a
    concurrency test would discard exactly the courses that persistence exists
    to capture.
    """
    active = scoring.is_cns_active(cns_class, activity_threshold)
    window_days = (rt_end - rt_start).days + 1
    if window_days <= 0:
        return Timing(active, False, False, False, "manual_review",
                      flag="invalid_rt_interval", unresolved=True)

    if half_life_days is None or not admins:
        return Timing(active, False, False, False,
                      "non_cns" if not active else "not_concurrent",
                      covered_days=0, window_days=window_days, coverage_pct=0.0,
                      flag="no_administration_dates" if not admins
                           else "no_half_life")

    dates = sorted(set(admins))
    covered = sum(
        1 for i in range(window_days)
        if pk.is_present(dates, rt_start + timedelta(days=i),
                         half_life_days, k, model))
    pct = round(100.0 * covered / window_days, 1)
    concurrent = bool(active and pct >= float(threshold_pct))

    flag = ("borderline_within_5pct_of_threshold"
            if abs(pct - float(threshold_pct)) < 5.0 else "")

    return Timing(active, concurrent, False, False,
                  _category(active, concurrent, False),
                  covered_days=covered, window_days=window_days,
                  coverage_pct=pct, flag=flag)


def _category(active: bool, concurrent: bool, adjuvant: bool) -> str:
    if not active:
        return "non_cns"
    if concurrent and adjuvant:
        return "concurrent_and_adjuvant"
    if concurrent:
        return "concurrent"
    if adjuvant:
        return "adjuvant"
    return "not_concurrent"
