"""Command line interface. See `python -m cns_exposure --help`."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from . import pharmacokinetics as pk
from . import resources
from . import rules as rule_module
from . import ksweep, scoring
from .classify import aggregate, classify, read_csv
from .dictionary import DrugDictionary

EXAMPLES = resources.prefer_repo_path("examples")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="cns_exposure",
        description="Classify CNS-active systemic therapy exposure around "
                    "radiotherapy.")
    sub = parser.add_subparsers(dest="command", required=True)

    # -- classify ----------------------------------------------------------
    c = sub.add_parser("classify", help="classify therapy records in one cell")
    _io_args(c)
    _grid_args(c)
    c.add_argument("--out", help="output CSV (default: stdout)")
    c.add_argument("--aggregate", action="store_true",
                   help="collapse to one row per patient")

    # -- sweep -------------------------------------------------------------
    s = sub.add_parser("sweep", help="run the whole definition grid")
    _io_args(s)
    s.add_argument("--k", type=float, nargs="+", default=list(ksweep.DEFAULT_KS))
    s.add_argument("--rule", nargs="+", choices=rule_module.RULES,
                   default=list(ksweep.DEFAULT_RULES))
    s.add_argument("--model", nargs="+", choices=pk.MODELS,
                   default=list(ksweep.DEFAULT_MODELS))
    s.add_argument("--activity-threshold", type=int,
                   default=scoring.DEFAULT_ACTIVITY_THRESHOLD)
    s.add_argument("--out-dir", help="write long/summary/stability CSVs here")

    # -- compare -----------------------------------------------------------
    m = sub.add_parser("compare-models",
                       help="single-dose against compounding at matched k")
    _io_args(m)
    m.add_argument("--k", type=float, nargs="+", default=list(ksweep.DEFAULT_KS))
    m.add_argument("--rule", choices=rule_module.RULES,
                   default=rule_module.WINDOW)

    # -- swimmer -----------------------------------------------------------
    w = sub.add_parser("swimmer",
                       help="swimmer plot with the concentration curve drawn on it")
    _io_args(w)
    w.add_argument("--k", type=float, default=5.0)
    w.add_argument("--model", choices=pk.MODELS, default=pk.COMPOUNDING,
                   help="accumulation model (default: compounding)")
    w.add_argument("--activity-threshold", type=int,
                   default=scoring.DEFAULT_ACTIVITY_THRESHOLD)
    w.add_argument("--out", default="swimmer.png", help="output PNG")
    w.add_argument("--title", default="")
    w.add_argument("--deidentify", action="store_true",
                   help="replace patient identifiers with stable P01-style labels")
    w.add_argument("--lanes-csv", help="write every number on the figure here")

    # -- validate ----------------------------------------------------------
    v = sub.add_parser("validate", help="check the drug dictionary")
    v.add_argument("--dictionary")
    v.add_argument("--out", help="output CSV (default: stdout)")

    # -- dictionary --------------------------------------------------------
    d = sub.add_parser("dictionary", help="print the resolved dictionary")
    d.add_argument("--dictionary")
    d.add_argument("--out")

    # -- demo --------------------------------------------------------------
    sub.add_parser("demo", help="run the worked example shipped in examples/")

    args = parser.parse_args(argv)
    return _dispatch(args)


def _io_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--therapy", required=True, help="therapy records CSV")
    p.add_argument("--rt", required=True, help="RT records CSV")
    p.add_argument("--administrations",
                   help="administration dates CSV; required by the coverage rule")
    p.add_argument("--dictionary", help="dictionary directory")


def _grid_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--k", type=float, default=5.0,
                   help="half-lives counted as present (default: 5)")
    p.add_argument("--rule", choices=rule_module.RULES,
                   default=rule_module.WINDOW)
    p.add_argument("--model", choices=pk.MODELS, default=pk.BOX,
                   help="accumulation model (default: box, the published rule)")
    p.add_argument("--activity-threshold", type=int,
                   default=scoring.DEFAULT_ACTIVITY_THRESHOLD,
                   help="lowest CNS class counting as active (default: 1)")
    p.add_argument("--attribution-window", type=float,
                   help="days after the last fraction in which a newly started "
                        "drug counts as adjuvant initiation (default: unbounded)")


def _dispatch(args) -> int:
    if args.command == "demo":
        return _demo()

    dictionary = DrugDictionary.load(getattr(args, "dictionary", None))

    if args.command == "validate":
        issues = dictionary.validate()
        _emit(issues, getattr(args, "out", None))
        if issues.empty:
            print("No issues.", file=sys.stderr)
        else:
            counts = issues.severity.value_counts().to_dict()
            print(f"{len(issues)} issue(s): "
                  + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())),
                  file=sys.stderr)
        return 1 if (not issues.empty and (issues.severity == "error").any()) else 0

    if args.command == "dictionary":
        _emit(dictionary.summary(), getattr(args, "out", None))
        return 0

    therapy = read_csv(args.therapy)
    rt = read_csv(args.rt)
    admins = read_csv(args.administrations) if args.administrations else None

    if args.command == "classify":
        if args.rule == rule_module.COVERAGE and admins is None:
            print("The coverage rule needs --administrations.", file=sys.stderr)
            return 2
        out = classify(therapy, rt, dictionary, administrations=admins,
                       k=args.k, rule=args.rule, model=args.model,
                       activity_threshold=args.activity_threshold,
                       attribution_window_days=args.attribution_window)
        if args.aggregate:
            out = aggregate(out)
        _emit(out, args.out)
        return 0

    if args.command == "sweep":
        if rule_module.COVERAGE in args.rule and admins is None:
            print("The coverage rule needs --administrations.", file=sys.stderr)
            return 2
        res = ksweep.sweep(therapy, rt, dictionary, administrations=admins,
                           ks=args.k, rules=args.rule, models=args.model,
                           activity_threshold=args.activity_threshold)
        if args.out_dir:
            out_dir = Path(args.out_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            for name, frame in res.items():
                frame.to_csv(out_dir / f"{name}.csv", index=False)
                print(f"wrote {out_dir / f'{name}.csv'}", file=sys.stderr)
        else:
            _print(res["summary"])
        return 0

    if args.command == "swimmer":
        if admins is None:
            print("The swimmer plot needs --administrations: it draws one lane "
                  "per dose-level record.", file=sys.stderr)
            return 2
        from .plot import build_records, render, to_frame

        records = build_records(therapy, rt, admins, dictionary, k=args.k,
                                model=args.model,
                                activity_threshold=args.activity_threshold)
        if not records:
            print("No record has both administration dates and a half-life; "
                  "nothing to draw.", file=sys.stderr)
            return 2
        out = render(records, args.out, k=args.k, model=args.model,
                     title=args.title, deidentify=args.deidentify)
        print(f"wrote {out} ({len(records)} lanes)", file=sys.stderr)
        if args.lanes_csv:
            to_frame(records).to_csv(args.lanes_csv, index=False)
            print(f"wrote {args.lanes_csv}", file=sys.stderr)
        return 0

    if args.command == "compare-models":
        out = ksweep.compare_models(therapy, rt, dictionary,
                                    administrations=admins,
                                    ks=args.k, rule=args.rule)
        _print(out)
        return 0

    return 2


def _demo() -> int:
    dictionary = DrugDictionary.load()
    therapy = read_csv(EXAMPLES / "therapy_records.csv")
    rt = read_csv(EXAMPLES / "rt_records.csv")
    admins = read_csv(EXAMPLES / "administrations.csv")

    _header("DICTIONARY VALIDATION")
    issues = dictionary.validate()
    if issues.empty:
        print("  no issues")
    else:
        _print(issues)

    _header("PUBLISHED RULE  --  window, k = 5, single-dose half-life")
    published = classify(therapy, rt, dictionary, k=5.0,
                         rule=rule_module.WINDOW, model=pk.BOX)
    _print(published[["patient_id", "drug", "cns_class", "half_life_days",
                      "effective_window_end", "cCNS_aD", "aCNS_aD",
                      "category", "flag"]])

    _header("SAME CELL, COMPOUNDING HALF-LIFE")
    compounded = classify(therapy, rt, dictionary, k=5.0,
                          rule=rule_module.WINDOW, model=pk.COMPOUNDING)
    merged = published[["patient_id", "drug", "effective_window_end"]].merge(
        compounded[["patient_id", "drug", "R", "extra_half_lives",
                    "effective_window_end"]],
        on=["patient_id", "drug"], suffixes=("_single_dose", "_compounding"))
    merged["days_added"] = [
        (b - a).days if a is not None and b is not None else None
        for a, b in zip(merged.effective_window_end_single_dose,
                        merged.effective_window_end_compounding)]
    _print(merged)

    _header("SINGLE-DOSE AGAINST COMPOUNDING, ACROSS k")
    _print(ksweep.compare_models(therapy, rt, dictionary, rule=rule_module.WINDOW))
    print("\n  mean_days_added is constant down the table: the extension is")
    print("  log2(R) * half-life, which does not depend on k. Sweeping k does")
    print("  not explore compounding, which is why it is a separate axis.")
    print("  cCNS_aD_delta is not constant, and is not monotone either -- it")
    print("  counts how many units that fixed extension happens to reclassify,")
    print("  which depends on where each window end falls relative to RT.")

    _header("THE FULL GRID")
    res = ksweep.sweep(therapy, rt, dictionary, administrations=admins)
    _print(res["summary"][["cell", "rule", "k", "n_units", "cCNS_aD",
                           "cCNS_aD_pct", "aCNS_aD", "aCNS_aD_pct"]])

    _header("STABILITY ACROSS THE GRID")
    _print(res["stability"])
    unstable = (res["stability"].verdict == "unstable").sum()
    print(f"\n  {unstable} of {len(res['stability'])} patients change "
          f"concurrency call across the grid.")
    print("  For those patients the exposure variable is a property of the")
    print("  definition, not of the treatment they received.")
    return 0


def _header(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def _print(df: pd.DataFrame) -> None:
    with pd.option_context("display.width", 200, "display.max_columns", 40,
                           "display.max_colwidth", 30):
        print(df.to_string(index=False))


def _emit(df: pd.DataFrame, out: str | None) -> None:
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
        print(f"wrote {out} ({len(df)} rows)", file=sys.stderr)
    else:
        df.to_csv(sys.stdout, index=False)
