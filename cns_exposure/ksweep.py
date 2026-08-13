"""
Sweeping the definition grid.

There is no single correct k, and there is no single correct concurrency rule.
Choosing one and reporting only that result hides the fact that a different
defensible choice would have produced a different exposure variable, and
possibly a different conclusion. The sweep runs the whole grid and reports the
distribution.

The grid has three axes:

  rule    window | coverage          two questions, not two precisions
  k       1..5                       how many half-lives count as present
  model   box | compounding          single-dose or steady-state tail

`k` and `model` are orthogonal: compounding adds log2(R) half-lives to the
window, and log2(R) does not depend on k. Sweeping k does not explore
compounding and cannot substitute for it, which is why the model is a separate
axis rather than a k offset.

WHAT TO READ

`summary` gives exposure counts per cell. `stability` gives, for each analysis
unit, the fraction of cells calling it exposed. A unit exposed in 10 of 10
cells is exposed under any defensible definition; one exposed in 3 of 10 is an
artefact of the definition and should be visible as such.

The honest headline from a sweep is the range across cells, not the cell that
happens to give the cleanest number. A sweep used to select a cell is a
multiple-comparisons machine and its output is not confirmatory.
"""

from __future__ import annotations

import pandas as pd

from . import pharmacokinetics as pk
from . import rules as rule_module
from . import scoring
from .classify import aggregate, classify
from .dictionary import DrugDictionary

DEFAULT_KS = (1, 2, 3, 4, 5)
DEFAULT_RULES = (rule_module.WINDOW, rule_module.COVERAGE)
DEFAULT_MODELS = (pk.COMPOUNDING,)


def sweep(therapy: pd.DataFrame,
          rt: pd.DataFrame,
          drug_dictionary: DrugDictionary,
          administrations: pd.DataFrame | None = None,
          ks=DEFAULT_KS,
          rules=DEFAULT_RULES,
          models=DEFAULT_MODELS,
          activity_threshold: int = scoring.DEFAULT_ACTIVITY_THRESHOLD,
          unit: str = "patient_id",
          **kwargs) -> dict[str, pd.DataFrame]:
    """
    Run every cell of the grid.

    Returns three frames:
      long       every classified row from every cell, tagged with `cell`
      summary    one row per cell: exposure counts and rates
      stability  one row per analysis unit: how often it is called exposed
    """
    frames = []
    for rule in rules:
        for model in models:
            for k in ks:
                out = classify(therapy, rt, drug_dictionary,
                               administrations=administrations,
                               k=k, rule=rule, model=model,
                               activity_threshold=activity_threshold, **kwargs)
                out.insert(0, "cell", cell_name(rule, k, model))
                frames.append(out)

    populated = [f for f in frames if not f.empty]
    long = (pd.concat(populated, ignore_index=True) if populated
            else pd.DataFrame(columns=["cell"]))
    return {"long": long,
            "summary": summarize(long, unit=unit),
            "stability": stability(long, unit=unit)}


def cell_name(rule: str, k: float, model: str) -> str:
    suffix = "" if model == pk.COMPOUNDING else f"_{model}"
    return f"{rule}_k{int(k) if float(k).is_integer() else k}{suffix}"


def summarize(long: pd.DataFrame, unit: str = "patient_id") -> pd.DataFrame:
    """Exposure counts per cell, at the analysis unit."""
    if long.empty:
        return pd.DataFrame()

    rows = []
    for cell, g in long.groupby("cell", sort=True):
        agg = aggregate(g, by=unit)
        n = len(agg)
        first = g.iloc[0]
        rows.append({
            "cell": cell,
            "rule": first["rule"],
            "k": first["k"],
            "accumulation_model": first["accumulation_model"],
            "activity_threshold": first["activity_threshold"],
            "n_units": n,
            "n_records": len(g),
            "cCNS_aD": int(agg.cCNS_aD.sum()) if n else 0,
            "aCNS_aD": int(agg.aCNS_aD.sum()) if n else 0,
            "aCNS_aD_persistence": int(agg.aCNS_aD_persistence.sum()) if n else 0,
            "aCNS_aD_initiation": int(agg.aCNS_aD_initiation.sum()) if n else 0,
            "cCNS_aD_pct": round(100.0 * agg.cCNS_aD.mean(), 1) if n else None,
            "aCNS_aD_pct": round(100.0 * agg.aCNS_aD.mean(), 1) if n else None,
            "n_units_unresolved": int((agg.n_unresolved > 0).sum()) if n else 0,
            "n_unmatched_drugs": int((g.matched == 0).sum()),
            "n_flagged": int((g.flag.fillna("") != "").sum()),
        })

    out = pd.DataFrame(rows)
    return out.sort_values(["rule", "accumulation_model", "k"]).reset_index(drop=True)


def stability(long: pd.DataFrame, unit: str = "patient_id") -> pd.DataFrame:
    """
    How often each analysis unit is called exposed across the grid.

    `verdict` is the reading a reviewer needs:
      always      exposed in every cell; the call does not depend on the rule
      never       exposed in no cell
      unstable    exposed in some cells only; the call IS the definition
    """
    if long.empty:
        return pd.DataFrame()

    cells = sorted(long["cell"].unique())
    per_cell = []
    for cell in cells:
        agg = aggregate(long[long.cell == cell], by=unit)
        agg = agg[[unit, "cCNS_aD", "aCNS_aD"]].copy()
        agg["cell"] = cell
        per_cell.append(agg)

    stacked = pd.concat(per_cell, ignore_index=True)
    grouped = stacked.groupby(unit, sort=True)

    out = pd.DataFrame({
        unit: [k for k, _ in grouped],
        "n_cells": [len(g) for _, g in grouped],
        "cCNS_aD_cells": [int(g.cCNS_aD.sum()) for _, g in grouped],
        "aCNS_aD_cells": [int(g.aCNS_aD.sum()) for _, g in grouped],
    })
    out["cCNS_aD_fraction"] = (out.cCNS_aD_cells / out.n_cells).round(3)
    out["aCNS_aD_fraction"] = (out.aCNS_aD_cells / out.n_cells).round(3)
    out["verdict"] = [
        "always" if c == n else "never" if c == 0 else "unstable"
        for c, n in zip(out.cCNS_aD_cells, out.n_cells)
    ]
    return out


def compare_models(therapy: pd.DataFrame,
                   rt: pd.DataFrame,
                   drug_dictionary: DrugDictionary,
                   administrations: pd.DataFrame | None = None,
                   ks=DEFAULT_KS,
                   rule: str = rule_module.WINDOW,
                   **kwargs) -> pd.DataFrame:
    """
    Single-dose against compounding at matched k and rule.

    Two different quantities are reported and they behave differently.

    `days_added` is how much longer the effective window runs under
    compounding. It equals log2(R) * h and is INVARIANT IN k: the same column
    repeats down the table. That invariance is the claim -- compounding is a
    fixed extension of the window, so sweeping k does not explore it and cannot
    substitute for it.

    `cCNS_aD_delta` is how many analysis units that extension actually
    reclassifies. It is NOT invariant in k and is not monotone in it either.
    Whether a longer window changes a classification depends on where its end
    lands relative to the radiotherapy interval, so the extension can move
    several units at one k, none at the next, and nothing at all once k alone
    already carries the window past the interval. A zero in this column means
    the reclassification saturated at that k, not that compounding did nothing.
    """
    res = sweep(therapy, rt, drug_dictionary, administrations=administrations,
                ks=ks, rules=(rule,), models=(pk.BOX, pk.COMPOUNDING), **kwargs)
    s = res["summary"]
    long = res["long"]

    box = s[s.accumulation_model == pk.BOX].set_index("k")
    comp = s[s.accumulation_model == pk.COMPOUNDING].set_index("k")
    shared = sorted(set(box.index) & set(comp.index))

    added = _days_added(long, rule, shared)

    return pd.DataFrame({
        "k": shared,
        "rule": rule,
        "cCNS_aD_box": [int(box.at[k, "cCNS_aD"]) for k in shared],
        "cCNS_aD_compounding": [int(comp.at[k, "cCNS_aD"]) for k in shared],
        "cCNS_aD_delta": [int(comp.at[k, "cCNS_aD"] - box.at[k, "cCNS_aD"])
                          for k in shared],
        "aCNS_aD_box": [int(box.at[k, "aCNS_aD"]) for k in shared],
        "aCNS_aD_compounding": [int(comp.at[k, "aCNS_aD"]) for k in shared],
        "aCNS_aD_delta": [int(comp.at[k, "aCNS_aD"] - box.at[k, "aCNS_aD"])
                          for k in shared],
        "mean_days_added": [added[k][0] for k in shared],
        "max_days_added": [added[k][1] for k in shared],
    })


def _days_added(long: pd.DataFrame, rule: str, ks) -> dict:
    """Mean and max window extension per k, over records the model can move."""
    out = {}
    for k in ks:
        sel = long[(long.rule == rule) & (long.k == k)]
        b = sel[sel.accumulation_model == pk.BOX]
        c = sel[sel.accumulation_model == pk.COMPOUNDING]
        merged = b[["patient_id", "drug", "effective_window_end"]].merge(
            c[["patient_id", "drug", "effective_window_end"]],
            on=["patient_id", "drug"], suffixes=("_b", "_c"))
        deltas = [(y - x).days
                  for x, y in zip(merged.effective_window_end_b,
                                  merged.effective_window_end_c)
                  if x is not None and y is not None]
        moved = [d for d in deltas if d > 0]
        out[k] = (round(sum(moved) / len(moved), 1) if moved else 0.0,
                  max(deltas) if deltas else 0)
    return out
