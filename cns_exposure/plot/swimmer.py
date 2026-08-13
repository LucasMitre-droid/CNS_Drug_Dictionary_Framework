"""
The swimmer plot, with the concentration curve drawn on it.

A presence bar states that an agent was on board. It is the OUTPUT of the
pharmacokinetic model, with the model itself invisible, so a reader cannot tell
a drug that barely cleared the threshold from one sitting an order of magnitude
above it, and cannot see accumulation happening at all. This figure draws the
model underneath the answer:

    C(t) = Σᵢ (1/2)^((t − aᵢ)/h)      over every dose aᵢ ≤ t

in units of one single-dose peak — the same direct summation the coverage rule
uses to decide presence. Each administration steps the curve up, it decays
between doses, under repeated dosing it climbs to the plateau
R = 1/(1 − (1/2)^(τ/h)), and after the last dose it decays from R·D rather than
from D. That last clause is the entire content of the compounding model, and on
this figure it is visible rather than asserted.

THE LANE BASELINE IS THE THRESHOLD

C is drawn on a log₂ axis mapped so the lane's dashed baseline sits exactly at

    C = (1/2)^k          the presence threshold at this k

and the top of the lane at C = 2^LOG2_TOP. "Curve above the baseline" and
"presence bar drawn" are therefore the same statement, and the point where the
tail crosses the baseline is the moment the drug stops counting as present.

A linear axis cannot show this. At k = 5 the threshold is 0.031 while an agent
dosed daily against a 4.5-day half-life plateaus near 7 — a span of more than
two orders of magnitude, which pins the entire decision boundary to the axis
floor. The floor is set FLOOR_PAD log₂ units *below* the threshold so a
decaying tail is seen crossing it rather than being clipped flat against the
bottom of the lane.

Days before the first dose carry C = 0, which has no logarithm. They are masked
out of the curve rather than pinned to the floor, so the line begins at the
first administration instead of running in along the bottom of the lane and
implying a measurement that was never made.

Lanes are about three times the height of a plain swimmer so the sawtooth is
legible; the figure is correspondingly tall.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from .. import pharmacokinetics as pk
from . import palette

# -- lane geometry ----------------------------------------------------------
LOG2_TOP = 3.0        # top of a lane = 8x a single-dose peak
FLOOR_PAD = 2.5       # log2 units drawn BELOW threshold, so the crossing shows
CURVE_LO = -0.30      # lane-relative y of the axis floor
CURVE_HI = 0.34       # lane-relative y of C = 2^LOG2_TOP

AX_LEFT, AX_RIGHT = 0.215, 0.965
GROUP_X = -0.225
HEADING_X = AX_LEFT + GROUP_X * (AX_RIGHT - AX_LEFT)

LEAD_IN_DAYS = 90.0
GRID_STEP_DAYS = 0.25
MIN_GRID_STEP_DAYS = 0.002   # ~3 minutes; bounds the array for very short t½
MAX_MARKS_DRAWN = 60
GROUP_LABEL_CHARS = 26


@dataclass
class SwimmerRecord:
    """One patient-drug record: a single lane."""
    subject: str
    agent: str
    rt_start: date
    rt_end: date
    half_life_days: float
    administrations: list[date] = field(default_factory=list)
    called: bool = False
    coverage_pct: float = 0.0
    covered_days: int = 0

    @property
    def window_days(self) -> int:
        return (self.rt_end - self.rt_start).days + 1

    @property
    def offsets(self) -> np.ndarray:
        """Administration days relative to the first fraction."""
        return np.array([(a - self.rt_start).days for a in self.administrations],
                        dtype=float)


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------
def concentration_curve(offsets: np.ndarray, t: np.ndarray,
                        half_life: float) -> np.ndarray:
    """
    C(t) over a grid, in units of one single-dose peak.

    The vectorised twin of `pharmacokinetics.concentration`, evaluated by direct
    summation over an outer difference rather than by a running recurrence. A
    recurrence multiplies by a decay factor once per step and drifts downward,
    which matters precisely at the threshold crossing this figure is drawn to
    show. The two implementations agree to floating-point tolerance; the scalar
    one remains authoritative for classification.
    """
    if offsets.size == 0 or half_life <= 0:
        return np.zeros_like(t)
    dt = t[:, None] - offsets[None, :]
    return np.where(dt >= 0, np.power(0.5, np.abs(dt) / half_life), 0.0).sum(axis=1)


def _span(k: float) -> float:
    return LOG2_TOP + float(k) + FLOOR_PAD


def _frac_to_y(frac, y: float):
    return y + CURVE_LO + np.clip(frac, 0.0, 1.0) * (CURVE_HI - CURVE_LO)


def to_lane_y(c: np.ndarray, y: float, k: float) -> np.ndarray:
    """
    Map C onto a lane with log2 spacing.

    Anything at or below the floor is NaN and simply is not drawn. Clipping it
    to the floor instead would run a flat line the full width of the lane, which
    reads as a drug sitting at a low but constant concentration -- a claim about
    the data rather than about the plotted range. The floor is FLOOR_PAD log2
    units below the threshold, so the crossing this figure exists to show is
    still drawn in full before the line ends.

    C = 0 before the first dose has no logarithm and is masked by the same rule.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        l2 = np.log2(np.where(c > 0, c, np.nan))
    frac = (l2 + float(k) + FLOOR_PAD) / _span(k)
    out = _frac_to_y(frac, y)
    return np.where(np.isnan(frac) | (frac <= 0.0), np.nan, out)


def threshold_y(y: float, k: float) -> float:
    """Lane y of the presence threshold C = (1/2)^k."""
    return float(_frac_to_y(FLOOR_PAD / _span(k), y))


def presence_spans(record: SwimmerRecord, k: float,
                   model: str = pk.COMPOUNDING) -> list[tuple[float, float]]:
    """
    Merged intervals during which the agent is above threshold.

    Because every dose decays with the same half-life, the summed concentration
    after dose i is a single exponential reseeded there, so presence following
    it runs until `aᵢ + h·(k + log₂Sᵢ)` with `Sᵢ` the concentration just after
    that dose. The single-dose model is the special case `Sᵢ = 1`, a flat k
    half-lives — which is exactly why the two models coincide for a drug given
    once and diverge only through accumulation.
    """
    offsets = sorted(record.offsets.tolist())
    h = float(record.half_life_days)
    if not offsets or h <= 0:
        return []

    out: list[list[float]] = []
    for i, start in enumerate(offsets):
        if model == pk.BOX:
            extra = 0.0
        else:
            s_i = sum(0.5 ** ((start - earlier) / h) for earlier in offsets[:i + 1])
            extra = float(np.log2(s_i)) if s_i > 0 else 0.0
        end = start + h * (float(k) + extra)
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return [(a, b) for a, b in out]


# ---------------------------------------------------------------------------
# building lanes from the framework's own inputs
# ---------------------------------------------------------------------------
def build_records(therapy, rt, administrations, drug_dictionary,
                  k: float = 5.0, model: str = pk.COMPOUNDING,
                  coverage_threshold_pct: float = 50.0,
                  activity_threshold: int = 1) -> list[SwimmerRecord]:
    """
    Turn the framework's input tables into lanes.

    A lane needs administration dates, so this covers only the records the
    coverage rule can see. Patients without them are absent from the figure
    rather than drawn empty: a blank lane reads as "no drug", which is a
    different claim from "not abstracted to dose level". The count of records
    dropped for that reason is worth reporting in any caption.

    The concurrency call is taken from `rules.classify_coverage`, not recomputed
    here, so a green lane is by construction a record the analysis counts.
    """
    from ..classify import _admin_index, _rt_index, _to_date
    from ..rules import classify_coverage

    rt_by_patient = _rt_index(rt)
    admins_by_key = _admin_index(administrations, drug_dictionary)

    out: list[SwimmerRecord] = []
    for _, row in therapy.iterrows():
        pid = str(row["patient_id"]).strip()
        interval = rt_by_patient.get(pid)
        if interval is None:
            continue
        rt_start, rt_end = interval

        for entry in drug_dictionary.resolve_regimen(row["drug_or_regimen"]):
            # A combination scored as a unit has no administration record of its
            # own -- dose dates are recorded per component -- so its lane is the
            # union of its components' dates. The regimen is present whenever
            # any component is, which is the same rule the combination's
            # half-life follows in taking max() over components.
            keys = (entry.components if entry.is_combination
                    else (drug_dictionary.resolve_name(entry.name),))
            admins = sorted({d for key in keys
                             for d in admins_by_key.get((pid, key), [])})
            if not admins or entry.half_life_days is None:
                continue
            timing = classify_coverage(
                entry.cns_class, admins, rt_start, rt_end, entry.half_life_days,
                k=k, model=model, threshold_pct=coverage_threshold_pct,
                activity_threshold=activity_threshold)
            out.append(SwimmerRecord(
                subject=pid, agent=entry.name,
                rt_start=rt_start, rt_end=rt_end,
                half_life_days=float(entry.half_life_days),
                administrations=admins,
                called=bool(timing.concurrent),
                coverage_pct=float(timing.coverage_pct or 0.0),
                covered_days=int(timing.covered_days or 0)))
    return out


def to_frame(records: list[SwimmerRecord]):
    """Every number on the figure, as a table, so a lane can be traced."""
    import pandas as pd

    return pd.DataFrame([{
        "subject": r.subject,
        "agent": r.agent,
        "half_life_days": r.half_life_days,
        "rt_start": r.rt_start,
        "rt_end": r.rt_end,
        "rt_window_days": r.window_days,
        "n_administrations": len(r.administrations),
        "first_admin_rel_days": (min(r.offsets) if len(r.offsets) else None),
        "last_admin_rel_days": (max(r.offsets) if len(r.offsets) else None),
        "covered_days": r.covered_days,
        "coverage_pct": r.coverage_pct,
        "called_concurrent": int(r.called),
    } for r in records])


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------
def _lane_grid(xmin: float, xmax: float, half_life: float) -> np.ndarray:
    """
    Sampling grid for one lane, refined to the agent's own half-life.

    A fixed grid under-samples short-half-life agents badly. Temozolomide at
    h = 0.075 d rises and falls inside two hours, so on a quarter-day grid most
    of its doses fall between samples and the lane shows two spikes where the
    record holds five. Sampling at h/6 keeps every dose visible.

    The step is floored so that a very short half-life over a long axis cannot
    generate an unbounded array; at the floor the spikes are drawn narrow but
    are all present.
    """
    step = max(min(GRID_STEP_DAYS, float(half_life) / 6.0), MIN_GRID_STEP_DAYS)
    return np.arange(xmin, xmax + step, step)


def _agent_heading(name: str) -> tuple[str, float]:
    """Wrap a regimen one component per line and pick a font size.

    Components are wrapped rather than truncated: "tucatinib + trastuzumab +
    capecitabine" must not read as a two-drug regimen.
    """
    parts = [p.strip() for p in str(name).strip().split("+")]
    wrapped = "\n+ ".join(p[:GROUP_LABEL_CHARS] for p in parts)
    return wrapped, 7.4 if len(parts) == 1 else 6.6


def _blocks(records: list[SwimmerRecord]) -> list[tuple[str, int, int, int]]:
    """Contiguous same-agent runs as (agent, y_low, y_high, rank), top-down."""
    n = len(records)
    out, rank, start = [], 0, 0
    for i in range(1, n + 1):
        if i == n or records[i].agent != records[start].agent:
            out.append((records[start].agent, n - 1 - (i - 1), n - 1 - start, rank))
            rank, start = rank + 1, i
    return out


def _draw_lane(ax, axb, rec: SwimmerRecord, y: int, band: bool,
               xmin: float, xmax: float, k: float, model: str) -> None:
    window = float(rec.window_days)

    if band:
        ax.add_patch(Rectangle((xmin, y - 0.5), xmax - xmin, 1.0,
                               facecolor=palette.ROW_BAND, edgecolor="none",
                               zorder=0))
    ax.add_patch(Rectangle((0, y - 0.46), min(window, xmax), 0.92,
                           facecolor=palette.RT_BAND, edgecolor="none", zorder=1))

    # presence bars, thinned so the curve reads on top of them
    for s, e in presence_spans(rec, k, model):
        s, e = max(s, xmin), min(e, xmax)
        if e <= s:
            continue
        for a, b, inside in ((s, min(e, 0.0), False),
                             (max(s, 0.0), min(e, window), True),
                             (max(s, window), e, False)):
            if b <= a:
                continue
            ax.add_patch(Rectangle((a, y - 0.40), b - a, 0.13, zorder=3,
                                   facecolor=palette.fill_for(rec.called, inside),
                                   edgecolor="none"))

    # the concentration curve
    offsets = rec.offsets
    h = float(rec.half_life_days)
    if offsets.size and h > 0:
        grid = _lane_grid(xmin, xmax, h)
        c = concentration_curve(offsets, grid, h)
        yy = to_lane_y(c, y, k)
        thr_y = threshold_y(y, k)
        ax.plot([xmin, xmax], [thr_y, thr_y], color=palette.THRESHOLD, lw=0.55,
                ls=(0, (3, 2.5)), zorder=4, alpha=0.8)
        above = np.nan_to_num(c) >= pk.presence_threshold(k)
        ax.fill_between(grid, thr_y, yy, where=above, color=palette.CURVE,
                        alpha=0.13, zorder=4, linewidth=0)
        ax.plot(grid, yy, color=palette.CURVE, lw=0.85, zorder=5,
                solid_joinstyle="round")

    marks = [t for t in offsets.tolist() if xmin <= t <= xmax]
    if 0 < len(marks) <= MAX_MARKS_DRAWN:
        ax.vlines(marks, y - 0.46, y - 0.40, color=palette.INK, lw=0.5,
                  alpha=0.7, zorder=6)

    axb.barh(y, rec.coverage_pct, height=0.42,
             color=palette.CALLED if rec.called else palette.NOT_CALLED, zorder=3)
    axb.text(min(rec.coverage_pct, 100) + 2.5, y,
             f"{palette.STATUS_GLYPH[rec.called]} {rec.coverage_pct:.0f}%  "
             f"({rec.covered_days}/{rec.window_days}d)",
             va="center", ha="left", fontsize=6.4,
             color=palette.INK if rec.called else palette.INK_2, zorder=6)


def render(records: list[SwimmerRecord], path: str | Path, k: float = 5.0,
           model: str = pk.COMPOUNDING, title: str = "",
           coverage_threshold_pct: float = 50.0,
           deidentify: bool = False) -> Path:
    """Draw one group's swimmer plot with concentration curves and write it."""
    if not records:
        raise ValueError("no records to draw")
    pk._check_model(model)
    palette.apply_rc()

    records = sorted(records, key=lambda r: (r.agent.strip().lower(),
                                             -r.coverage_pct, r.subject))
    if deidentify:
        labels_by_subject = {
            s: f"P{i + 1:02d}"
            for i, s in enumerate(sorted({r.subject for r in records}))
        }
    else:
        labels_by_subject = {r.subject: r.subject for r in records}

    n = len(records)
    xmin = -LEAD_IN_DAYS
    xmax = max(float(np.percentile([r.window_days for r in records], 95)) + 30,
               60.0)

    height_in = 0.55 * n + 4.2
    fig = plt.figure(figsize=(16.0, height_in))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.45, 1.0], wspace=0.055,
                          left=AX_LEFT, right=AX_RIGHT,
                          top=1 - 2.9 / height_in, bottom=1.35 / height_in)
    ax = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1], sharey=ax)

    blocks = _blocks(records)
    rank_of = {y: rank for _a, lo, hi, rank in blocks for y in range(lo, hi + 1)}

    for i, rec in enumerate(records):
        y = n - 1 - i
        _draw_lane(ax, axb, rec, y, rank_of[y] % 2 == 0, xmin, xmax, k, model)

    ax.set_yticks(range(n))
    ax.set_yticklabels([labels_by_subject[r.subject] for r in records][::-1],
                       fontsize=6.2, color=palette.INK_2)

    for agent, y_lo, y_hi, _rank in blocks:
        text, size = _agent_heading(agent)
        ax.text(GROUP_X, 0.5 * (y_lo + y_hi), text,
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=size, weight="bold", color=palette.INK,
                clip_on=False, linespacing=0.98)
        if y_lo > 0:
            for axis in (ax, axb):
                axis.axhline(y_lo - 0.5, color=palette.GRID, lw=0.8, zorder=2)

    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(xmin, xmax)
    ax.axvline(0, color=palette.INK_2, lw=0.9, zorder=7)
    ax.set_xlabel("Days relative to RT start", fontsize=9)
    ax.grid(axis="x", color=palette.GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=8)

    axb.axvline(coverage_threshold_pct, color=palette.INK_2, lw=1.0,
                ls=(0, (4, 3)), zorder=2)
    axb.set_xlim(0, 185)
    axb.set_xticks([0, 50, 100])
    axb.set_xlabel("% of RT window with drug present", fontsize=9)
    axb.grid(axis="x", color=palette.GRID, lw=0.6, zorder=0)
    axb.set_axisbelow(True)
    for s in ("top", "right", "left"):
        axb.spines[s].set_visible(False)
    axb.tick_params(axis="y", length=0, labelleft=False)
    axb.tick_params(axis="x", labelsize=8)

    n_called = sum(1 for r in records if r.called)
    label = "compounding" if model == pk.COMPOUNDING else "single-dose"
    fig.suptitle(title or f"Drug concentration through the RT window "
                          f"({label} PK, k = {k:g})",
                 x=HEADING_X, ha="left", fontsize=15, weight="bold",
                 y=1 - 0.6 / height_in)
    fig.text(HEADING_X, 1 - 1.25 / height_in,
             "Blue curve: C(t) = Σ over doses of (1/2)^((t − dose)/t½), in units "
             "of one single-dose peak — the same summation the coverage rule "
             "uses. Each dose steps it up; it accumulates toward the plateau R "
             "and decays from R×D after the last dose.",
             ha="left", fontsize=9, color=palette.INK_2)
    fig.text(HEADING_X, 1 - 1.72 / height_in,
             f"Dashed line is the presence threshold C = (1/2)^{k:g} = "
             f"{pk.presence_threshold(k):.4g}; lane top is 8× a single-dose peak "
             f"and the floor is {FLOOR_PAD:g} log₂ units below threshold, so the "
             f"tail is seen crossing it. Curve above the dashed line = presence "
             f"bar drawn.",
             ha="left", fontsize=9, color=palette.INK_2)
    fig.text(HEADING_X, 1 - 2.19 / height_in,
             f"{n_called}/{n} drug records counted concurrent "
             f"(≥ {coverage_threshold_pct:.0f}% of RT-window days covered). "
             f"One lane per patient-drug record, grouped by agent.",
             ha="left", fontsize=9, color=palette.MUTED)

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=palette.RT_BAND, label="RT window"),
        Line2D([0], [0], color=palette.CURVE, lw=1.4, label=f"C(t), {label}"),
        Line2D([0], [0], color=palette.THRESHOLD, lw=1.0, ls=(0, (3, 2.5)),
               label=f"threshold (1/2)^{k:g}"),
        Rectangle((0, 0), 1, 1, facecolor=palette.CALLED,
                  label="present — counted concurrent"),
        Rectangle((0, 0), 1, 1, facecolor=palette.NOT_CALLED,
                  label="present — not counted"),
        Line2D([0], [0], color=palette.INK, lw=0.8, label="administration"),
    ]
    fig.legend(handles=handles, loc="upper right",
               bbox_to_anchor=(AX_RIGHT, 1 - 0.5 / height_in), ncol=3,
               frameon=False, fontsize=9, handlelength=1.6, columnspacing=1.4,
               labelcolor=palette.INK_2)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path
