# Worked example

Ten fictional cases. **No patient data is in this directory.** The records are
the synthetic dataset published with the preprint, chosen so that each case
exercises one behaviour of the framework.

```bash
python -m cns_exposure demo
```

| File | Contents |
|---|---|
| `therapy_records.csv` | `patient_id`, `drug_or_regimen`, `start_date`, `end_date` |
| `rt_records.csv` | `patient_id`, `rt_start_date`, `rt_end_date` |
| `administrations.csv` | dose dates, expanded from each course at its agent's `τ`; needed by the `coverage` rule |
| `expected_classification.csv` | published output of the original cell: `window`, `k = 5`, single-dose half-life |
| `case_notes.csv` | what each case is designed to test |

## The cases

| Case | Tests |
|---|---|
| PT001 | long half-life, stopped 61 days before RT; the effective window bridges the gap |
| PT002 | short half-life, cleared before RT starts → unrelated |
| PT003 | drug spans the RT interval and extends past it |
| PT004 | started after the last fraction → adjuvant initiation, not persistence |
| PT005 | combination override: `nivolumab + ipilimumab` → class 3, `h = max(25, 15.4)` |
| PT006 | class 0 agent present throughout RT → not exposure |
| PT007 | missing end date → flagged for review, never guessed |
| PT008 | drug absent from the dictionary → unmatched, class 0, visible in output |
| PT009 | single-fraction RT, `t₀ = t₁` |
| PT010 | started during RT and continued past it |

`expected_classification.csv` reproduces the published cell only. It is not a
fixture for the compounding model or the coverage rule, both of which are
extensions and deliberately give different answers — that difference is what
the demo shows.

## The swimmer plot

```bash
pip install -r ../requirements-plot.txt
python -m cns_exposure swimmer \
    --therapy examples/therapy_records.csv \
    --rt examples/rt_records.csv \
    --administrations examples/administrations.csv \
    --k 2 --model compounding --deidentify --out swimmer.png
```

Eight of the ten cases draw a lane. PT007 has no end date and PT008 is not in
the dictionary, so neither has a half-life to propagate and both are absent —
a blank lane would read as "no drug", which is a different claim from "not
classifiable".

PT010 (osimertinib, dosed daily against a 2-day half-life) shows the sawtooth
clearly: each dose steps the curve up and it accumulates toward `R = 3.41`.
PT003 and PT005 are three-weekly antibodies whose half-life is close to their
dosing interval, so their curves climb in shallow steps instead. PT005 is the
combination case — its lane unions the administration dates of both components,
since a combination scored as a unit has no dose record of its own.

## What the demo demonstrates

Under the published rule, 5 of 10 cases are concurrent. Switching to the
compounding half-life at the same `k = 5` adds a mean of 21.8 days to the
windows it can move, and up to 29. Across the ten-cell grid, three of the ten
cases change their concurrency call — for those, the exposure variable is a
property of the definition rather than of the treatment given.
