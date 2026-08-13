"""
CNS activity scoring.

A drug's CNS activity class is computed from three declared evidence
components rather than read off a curated list:

    CNS_SCORE = P_iORR + P_consensus + P_IT

P_iORR carries the most weight because published intracranial objective
response rate is the most direct clinical evidence that an agent works inside
the brain. P_consensus admits agents the field treats as CNS-relevant when the
response literature is sparse or histology-specific. P_IT exists because
intrathecal delivery bypasses the blood-brain barrier by design, which is a
qualitatively different basis for CNS activity than systemic penetration, and
without it agents such as methotrexate are under-classified.

The score maps to an ordinal class 0-3. Note that the map is not the identity:
scores 2 and 3 both give class 2, so an agent reaching class 3 needs either a
>= 60% iORR plus consensus endorsement, or intrathecal use plus supporting
evidence. No single component can produce class 3 alone.

The class is separate from the evidence grade G, which records how strong the
underlying literature is. A drug may hold class 2 on a grade C base, or class 1
on a grade D base. The grade never modifies the class threshold or the
effective-window equation; it is carried so that a reader can see which
classifications rest on prospective CNS trials and which rest on extrapolation.
"""

from __future__ import annotations

# iORR tier lower bounds, in percent, highest first. The strongest available
# peer-reviewed intracranial response estimate selected under the dictionary's
# prespecified evidence hierarchy is the value scored.
IORR_TIERS: tuple[tuple[float, int, str], ...] = (
    (60.0, 3, ">=60%"),
    (40.0, 2, "40-59%"),
    (20.0, 1, "20-39%"),
    (0.0, 0, "<20%"),
)

TIER_POINTS: dict[str, int] = {label: pts for _, pts, label in IORR_TIERS}

CONSENSUS_POINTS = 1
INTRATHECAL_POINTS = 2

CLASS_LABELS: dict[int, str] = {
    0: "No CNS activity",
    1: "Minimal CNS activity",
    2: "Moderate CNS activity",
    3: "High CNS activity",
}

EVIDENCE_GRADES: dict[str, str] = {
    "A": "Prospective CNS-specific trial or prespecified CNS endpoint",
    "B": "Prospective subgroup with standardized intracranial response assessment",
    "C": "Retrospective cohort with defined CNS response criteria and adequate "
         "sample size",
    "D": "Case series, post hoc subgroup, permissive trial eligibility, or "
         "consensus extrapolation",
    "E": "Mechanistic or pharmacokinetic plausibility without clinical "
         "intracranial response evidence",
}

# The threshold at which a class counts as CNS-active for exposure purposes.
# Declare it explicitly per study and keep it fixed for the primary analysis:
# it is the single most consequential free parameter in the framework, and
# moving it between >= 1 and >= 2 changes which drugs enter the exposure
# variable at all, before any timing rule is applied.
DEFAULT_ACTIVITY_THRESHOLD = 1


def iorr_tier(iorr_percent: float | None) -> str | None:
    """Tier label for a raw intracranial objective response rate in percent."""
    if iorr_percent is None:
        return None
    for lower, _, label in IORR_TIERS:
        if float(iorr_percent) >= lower:
            return label
    return None


def iorr_points(iorr_percent: float | None = None,
                tier: str | None = None) -> int:
    """
    Points from intracranial response, given either a raw percentage or a
    pre-assigned tier label. A raw percentage takes precedence.

    An agent with no published iORR scores 0 here; it can still reach class 1
    through consensus endorsement alone.
    """
    if iorr_percent is not None:
        for lower, pts, _ in IORR_TIERS:
            if float(iorr_percent) >= lower:
                return pts
        return 0
    if tier is not None:
        return TIER_POINTS.get(str(tier).strip(), 0)
    return 0


def consensus_points(endorsed: bool | str | None) -> int:
    """
    Points from consensus endorsement by NCCN, RANO or ESMO for CNS use, or
    from routine permission in dedicated CNS trials.

    Mere eligibility of patients with stable brain metastases on a systemic
    trial does not qualify. The guideline or trial must identify the agent as
    CNS-active, CNS-directed, or clinically relevant for intracranial disease.
    """
    return CONSENSUS_POINTS if _truthy(endorsed) else 0


def intrathecal_points(intrathecal: bool | str | None) -> int:
    """Points from a clinically relevant CNS use that is intrathecal."""
    return INTRATHECAL_POINTS if _truthy(intrathecal) else 0


def cns_score(iorr_percent: float | None = None,
              tier: str | None = None,
              consensus: bool | str | None = False,
              intrathecal: bool | str | None = False) -> int:
    """Total CNS evidence score."""
    return (iorr_points(iorr_percent, tier)
            + consensus_points(consensus)
            + intrathecal_points(intrathecal))


def cns_class(score: int) -> int:
    """
    Ordinal CNS activity class from the total evidence score.

        >= 4  -> 3    2 or 3 -> 2    1 -> 1    0 -> 0
    """
    s = int(score)
    if s >= 4:
        return 3
    if s >= 2:
        return 2
    if s >= 1:
        return 1
    return 0


def class_label(cns_class_value: int) -> str:
    return CLASS_LABELS.get(int(cns_class_value), "Unknown")


def is_cns_active(cns_class_value,
                  threshold: int = DEFAULT_ACTIVITY_THRESHOLD) -> bool:
    """
    Whether a class counts as CNS-active at the declared threshold.

    An unparseable class is treated as not CNS-active, which is the
    conservative direction: an agent is never counted as exposure on the
    strength of a value nobody could read.
    """
    try:
        return int(float(str(cns_class_value).strip())) >= int(threshold)
    except (TypeError, ValueError):
        return False


def _truthy(v) -> bool:
    if v is None:
        return False
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    return s in {"1", "y", "yes", "true", "t"}
