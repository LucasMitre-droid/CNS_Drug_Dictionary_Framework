"""
End-to-end classification of therapy records into exposure variables.

Input is three tables, two required:

  therapy_records   patient_id, drug_or_regimen, start_date, end_date
  rt_records        patient_id, rt_start_date, rt_end_date
  administrations   patient_id, drug_or_regimen, admin_date       (optional)

Administration dates are needed only by the coverage rule. Their absence is not
an error: patients without them keep their window-rule classification and are
marked so, because dropping them would silently restrict the cohort to the
patients whose charts happened to be abstracted in more depth, which is not a
clinical criterion and is unlikely to be independent of anything being measured.

Output is one row per resolved agent per patient. Aggregation to the patient,
lesion or treatment-course level is left to the caller, since the correct unit
is a property of the study design and not of the exposure definition;
`aggregate` handles the common patient-level case.

Dates are parsed leniently and anything unparseable becomes missing and is
flagged. A record that cannot be read is reported, never guessed at.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import pandas as pd

from . import pharmacokinetics as pk
from . import rules, scoring
from .dictionary import DrugDictionary, normalize

THERAPY_COLUMNS = ("patient_id", "drug_or_regimen", "start_date", "end_date")
RT_COLUMNS = ("patient_id", "rt_start_date", "rt_end_date")
ADMIN_COLUMNS = ("patient_id", "drug_or_regimen", "admin_date")


def classify(therapy: pd.DataFrame,
             rt: pd.DataFrame,
             drug_dictionary: DrugDictionary,
             administrations: pd.DataFrame | None = None,
             k: float = 5.0,
             rule: str = rules.WINDOW,
             model: str = pk.BOX,
             activity_threshold: int = scoring.DEFAULT_ACTIVITY_THRESHOLD,
             attribution_window_days: float | None = None,
             coverage_threshold_pct: float = rules.DEFAULT_COVERAGE_THRESHOLD_PCT,
             ) -> pd.DataFrame:
    """
    Classify every therapy record under one cell of the definition grid.

    A cell is (rule, k, accumulation model, activity threshold). Changing any
    of the four changes the exposure variable, so all four are recorded on
    every output row: a classified table that does not carry its own
    definition cannot be compared to another one.
    """
    _require(therapy, THERAPY_COLUMNS, "therapy records")
    _require(rt, RT_COLUMNS, "RT records")
    if rule not in rules.RULES:
        raise ValueError(f"unknown rule {rule!r}; expected one of {rules.RULES}")
    pk._check_model(model)

    rt_by_patient = _rt_index(rt)
    admins_by_key = _admin_index(administrations, drug_dictionary)

    out: list[dict] = []
    for _, row in therapy.iterrows():
        pid = str(row["patient_id"]).strip()
        regimen = row["drug_or_regimen"]
        interval = rt_by_patient.get(pid)

        if interval is None:
            out.append(_row(pid, str(regimen), None, None, None,
                            flag="no_rt_record", rule=rule, k=k, model=model,
                            activity_threshold=activity_threshold))
            continue

        rt_start, rt_end = interval
        s, e = _to_date(row["start_date"]), _to_date(row["end_date"])

        for entry in drug_dictionary.resolve_regimen(regimen):
            accumulation = entry.accumulation if model == pk.COMPOUNDING else 1.0

            timing = rules.classify_window(
                entry.cns_class, s, e, rt_start, rt_end, entry.half_life_days,
                k=k, accumulation=accumulation, model=model,
                activity_threshold=activity_threshold,
                attribution_window_days=attribution_window_days)

            coverage_timing = None
            if rule == rules.COVERAGE:
                admins = admins_by_key.get((pid, normalize(entry.name)), [])
                coverage_timing = rules.classify_coverage(
                    entry.cns_class, admins, rt_start, rt_end,
                    entry.half_life_days, k=k, model=model,
                    threshold_pct=coverage_threshold_pct,
                    activity_threshold=activity_threshold)

            out.append(_row(
                pid, str(regimen), entry, timing, coverage_timing,
                rt_start=rt_start, rt_end=rt_end, start=s, end=e,
                rule=rule, k=k, model=model,
                activity_threshold=activity_threshold))

    return _typed(pd.DataFrame(out, columns=_OUTPUT_COLUMNS))


def aggregate(classified: pd.DataFrame, by: str = "patient_id") -> pd.DataFrame:
    """
    Collapse agent-level rows to one row per analysis unit.

    Exposure is the disjunction: a unit is exposed if any of its agents is. The
    contributing drug names are carried through, because an exposure flag whose
    provenance cannot be inspected is not auditable, and reviewers ask.

    A unit with no exposed agent but at least one unresolved record is reported
    as 0 with `n_unresolved > 0`, not as a clean 0. The disjunction cannot
    distinguish the two on its own -- an unknown absorbed into a False looks
    exactly like a negative -- so the count travels alongside. Treat a unit
    carrying unresolved records as missing, not unexposed, unless the records
    have been reviewed.
    """
    columns = [by, "cCNS_aD", "aCNS_aD", "aCNS_aD_persistence",
               "aCNS_aD_initiation", "concurrent_drugs", "adjuvant_drugs",
               "n_records", "n_unresolved", "n_flagged"]
    if classified.empty:
        return pd.DataFrame(columns=columns)

    def names(g, col):
        picked = g.loc[g[col] == 1, "drug"].dropna().astype(str)
        return "; ".join(sorted(set(picked)))

    rows = []
    for key, g in classified.groupby(by, sort=True):
        unresolved = (int(g.unresolved.fillna(0).sum())
                      if "unresolved" in g.columns else 0)
        rows.append({
            by: key,
            "cCNS_aD": int((g.cCNS_aD == 1).any()),
            "aCNS_aD": int((g.aCNS_aD == 1).any()),
            "aCNS_aD_persistence": int((g.aCNS_aD_persistence == 1).any()),
            "aCNS_aD_initiation": int((g.aCNS_aD_initiation == 1).any()),
            "concurrent_drugs": names(g, "cCNS_aD"),
            "adjuvant_drugs": names(g, "aCNS_aD"),
            "n_records": len(g),
            "n_unresolved": unresolved,
            "n_flagged": int((g.flag.fillna("") != "").sum()),
        })
    return pd.DataFrame(rows, columns=columns)


# ---------------------------------------------------------------------------
# input handling
# ---------------------------------------------------------------------------
def read_csv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str)


def _require(df: pd.DataFrame, columns, label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} missing required column(s): {', '.join(missing)}. "
            f"Expected {', '.join(columns)}.")


def _rt_index(rt: pd.DataFrame) -> dict[str, tuple[date, date]]:
    out: dict[str, tuple[date, date]] = {}
    for _, row in rt.iterrows():
        pid = str(row["patient_id"]).strip()
        s, e = _to_date(row["rt_start_date"]), _to_date(row["rt_end_date"])
        if s is not None and e is not None and e >= s:
            out[pid] = (s, e)
    return out


def _admin_index(administrations: pd.DataFrame | None,
                 drug_dictionary: DrugDictionary
                 ) -> dict[tuple[str, str], list[date]]:
    """
    (patient, canonical drug) -> administration dates.

    Keyed on the canonical name rather than the recorded string so that a
    course written as "KEYTRUDA" in one table and "pembrolizumab" in the other
    still joins.
    """
    if administrations is None or administrations.empty:
        return {}
    _require(administrations, ADMIN_COLUMNS, "administration records")

    out: dict[tuple[str, str], list[date]] = {}
    for _, row in administrations.iterrows():
        d = _to_date(row["admin_date"])
        if d is None:
            continue
        pid = str(row["patient_id"]).strip()
        for token in drug_dictionary.split_regimen(row["drug_or_regimen"]):
            out.setdefault((pid, token), []).append(d)
    return {k: sorted(set(v)) for k, v in out.items()}


def _to_date(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v).strip()
    if not s or s.lower() in {"nan", "nat", "none", "n/a", "na", "-", "?"}:
        return None
    parsed = pd.to_datetime(s, errors="coerce")
    return None if pd.isna(parsed) else parsed.date()


_NUMERIC_COLUMNS = (
    "matched", "is_combination", "cns_class", "half_life_days", "tau_days",
    "R", "extra_half_lives", "covered_days", "rt_window_days", "coverage_pct",
    "cCNS_aD", "aCNS_aD", "aCNS_aD_persistence", "aCNS_aD_initiation",
    "unresolved", "k", "activity_threshold",
)


def _typed(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pin numeric columns to float even when a cell produced none of them.

    A cell run under the window rule never fills `coverage_pct`, so that column
    arrives all-missing and would otherwise be inferred as object. Frames from
    different cells then disagree on dtype and cannot be concatenated for the
    sweep without pandas guessing. Pinning here keeps every cell's output
    shaped identically, which is what makes cells comparable at all.
    """
    for col in _NUMERIC_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


_OUTPUT_COLUMNS = [
    "patient_id", "regimen", "drug", "matched", "is_combination",
    "cns_class", "class_label", "evidence_grade", "half_life_days",
    "tau_days", "R", "extra_half_lives",
    "start_date", "end_date", "rt_start_date", "rt_end_date",
    "effective_window_end", "covered_days", "rt_window_days", "coverage_pct",
    "cCNS_aD", "aCNS_aD", "aCNS_aD_persistence", "aCNS_aD_initiation",
    "unresolved", "category", "flag",
    "rule", "k", "accumulation_model", "activity_threshold",
]


def _row(pid, regimen, entry, timing, coverage_timing, rt_start=None,
         rt_end=None, start=None, end=None, flag="", rule="", k=None,
         model="", activity_threshold=None) -> dict:
    concurrent = coverage_timing if coverage_timing is not None else timing

    flags = [f for f in (flag,
                         timing.flag if timing else "",
                         coverage_timing.flag if coverage_timing else "") if f]

    return {
        "patient_id": pid,
        "regimen": regimen,
        "drug": entry.name if entry else None,
        "matched": None if entry is None else int(entry.matched),
        "is_combination": None if entry is None else int(entry.is_combination),
        "cns_class": entry.cns_class if entry else None,
        "class_label": scoring.class_label(entry.cns_class) if entry else None,
        "evidence_grade": entry.evidence_grade if entry else None,
        "half_life_days": entry.half_life_days if entry else None,
        "tau_days": entry.tau_days if entry else None,
        "R": round(entry.accumulation, 4) if entry else None,
        "extra_half_lives": round(entry.extra_half_lives, 4) if entry else None,
        "start_date": start,
        "end_date": end,
        "rt_start_date": rt_start,
        "rt_end_date": rt_end,
        "effective_window_end": timing.effective_window_end if timing else None,
        "covered_days": coverage_timing.covered_days if coverage_timing else None,
        "rt_window_days": coverage_timing.window_days if coverage_timing else None,
        "coverage_pct": coverage_timing.coverage_pct if coverage_timing else None,
        # Unresolved records emit missing values, not zeros. See rules.Timing.
        "cCNS_aD": (None if concurrent is None or concurrent.unresolved
                    else int(concurrent.concurrent)),
        "aCNS_aD": (None if timing is None or timing.unresolved
                    else int(timing.adjuvant)),
        "aCNS_aD_persistence": (None if timing is None or timing.unresolved
                                else int(timing.adjuvant_persistence)),
        "aCNS_aD_initiation": (None if timing is None or timing.unresolved
                               else int(timing.adjuvant_initiation)),
        "unresolved": int(bool(
            (timing is not None and timing.unresolved)
            or (coverage_timing is not None and coverage_timing.unresolved)
            or timing is None)),
        "category": concurrent.category if concurrent else None,
        "flag": "; ".join(flags),
        "rule": rule,
        "k": k,
        "accumulation_model": model,
        "activity_threshold": activity_threshold,
    }
