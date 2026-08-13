# CNS Drug Dictionary Framework

Reproducible concurrent and adjuvant CNS-active drug exposure variables from
ordinary clinical therapy records.

Clinical oncology datasets store systemic therapy as a regimen label, a start
date and a stop date. Studies that then adjust for "concurrent CNS-active
systemic therapy" rarely state what made a drug CNS-active, or what made it
concurrent. Both choices are usually a curated list and a calendar cutoff, and
neither travels between studies: one paper's concurrent is same-day, another's
is within 90 days, and the drug lists behind them are not the same list.

This repository implements a framework that makes both choices explicit and
auditable. CNS activity comes from declared evidence components rather than
from a list. Concurrency comes from a pharmacokinetic exposure window rather
than from a calendar cutoff. And because the window still has free parameters,
the tooling sweeps them and reports the distribution instead of a single number.

It is the reference implementation for the preprint below, extended with the
compounding half-life model, the coverage rule, and the sweep.

## What the framework does

```
regimen string ─┬─ normalize, split on delimiters, resolve aliases
                ├─ dictionary lookup      → CNS class c, half-life h, grade G
                ├─ combination override    → joint class, h = max(components)
                ├─ effective window end    → e_EF = e + (k + log₂R)·h
                └─ timing rule vs RT       → cCNS-aD, aCNS-aD
```

A drug's CNS activity class is the sum of three declared evidence components:

| Component | Points | Basis |
|---|---|---|
| Intracranial ORR | 0–3 | ≥60% → 3, 40–59% → 2, 20–39% → 1, else 0 |
| Consensus endorsement | 0–1 | NCCN, RANO or ESMO for CNS use, or routine use in dedicated CNS trials |
| Intrathecal route | 0–2 | clinically relevant CNS use is intrathecal |

The total maps to a class: `≥4 → 3`, `2–3 → 2`, `1 → 1`, `0 → 0`. The map is
not the identity, so no single component reaches class 3 alone. A separate
evidence grade `A–E` records how strong the underlying literature is; it never
modifies the class or the window.

## What is new here

The preprint defines the exposure window as `e_EF = e + 5h` — the recorded stop
date plus five half-lives. That treats a course of therapy as a single dose. It
is wrong in one direction: a drug given q21d for a year does not begin decaying
from a single-dose peak on its last day.

**Compounding half-life.** Dosing every `τ` days drives concentration to a
plateau `R = 1 / (1 − (1/2)^(τ/h))` in single-dose-peak units, so decay starts
from `R·D` and the window becomes

```
e_EF = e + (k + log₂R)·h
```

`log₂R` does not depend on `k`. Compounding is an axis **orthogonal** to the
choice of `k`, not a rescaling of it — which is why sweeping `k` does not
explore it and the two must be swept separately. Whether `τ` exceeds `h`
decides everything: an antibody at `h = 27 d` given q20d reaches `R = 2.49` and
gains five weeks of window, while an oral agent at `h = 0.03 d` given daily
reaches `R = 1.000` and gains nothing. In the shipped dosing table, 9 of 37
agents cannot accumulate under their own schedule at all.

**A second concurrency rule.** Where administration dates exist, presence can
be evaluated per dose rather than per course, from the superposed concentration

```
C(t) = Σᵢ (1/2)^((t − aᵢ)/h)     present while C(t) ≥ (1/2)^k
```

For a single dose this reduces **exactly** to the window rule, since
`(1/2)^(δ/h) ≥ (1/2)^k ⟺ δ ≤ k·h`. That degeneracy is what makes a `k` grid
comparable across the two models: at any `k` they mean the same thing for a
one-shot drug and diverge only through accumulation.

**A sweep instead of a chosen cell.** `k`, the rule and the accumulation model
form a grid. Running one cell and reporting it hides that a different
defensible cell would have produced a different exposure variable. The sweep
runs all of them and reports, per analysis unit, the fraction of cells calling
it exposed.

**A dictionary validator.** The dictionary is a declared scientific input, so
the validator reports and never repairs. It checks that evidence components sum
to the recorded total, that the class follows from the score, that iORR tiers
score consistently, that CNS-active drugs carry a grade, that aliases resolve,
that combination components exist, and that each `R` is reproducible from its
own `τ` and `h`.

## The two rules

| | `window` | `coverage` |
|---|---|---|
| Concurrent when | effective window intersects the RT interval | agent on board on ≥50% of RT days |
| Needs | therapy start and stop dates | reconstructed administration dates |
| Covers | every patient in the dataset | only patients abstracted to dose level |
| Resolves | the course | the dose |
| Answers | was drug plausibly on board | how much drug was actually on board |

`window` is the published rule. `coverage` is stricter by construction: a drug
given once in week one of a six-week course satisfies `window` and fails
`coverage`. **They are different questions, not two precisions of one question,
and their outputs are not interchangeable.** Report which one produced a
result, and prefer to report both.

Adjuvant exposure is derived only from `window`, and is returned split:

- **persistence** — therapy stopped on or before the last fraction, but its
  effective window extends past it. Nothing new was taken.
- **initiation** — therapy started after the last fraction. A clinical decision.

Collapsing them loses the distinction between continued exposure and true
adjuvant treatment, which is usually the distinction the analysis is after.

## Install

```bash
pip install -r requirements.txt          # core: pandas only
pip install -r requirements-plot.txt     # adds matplotlib + numpy for figures
```

Python 3.10 or later.

## Use

Run the worked example, which classifies the published synthetic cases under
the original rule, the compounding rule, and the full grid:

```bash
python -m cns_exposure demo
```

Classify your own records under one cell:

```bash
python -m cns_exposure classify \
    --therapy therapy_records.csv \
    --rt rt_records.csv \
    --k 5 --rule window --model compounding \
    --out classified.csv
```

Sweep the grid:

```bash
python -m cns_exposure sweep \
    --therapy therapy_records.csv \
    --rt rt_records.csv \
    --administrations administrations.csv \
    --out-dir sweep/
```

writing `summary.csv` (exposure counts per cell), `stability.csv` (per unit,
how often it is called exposed) and `long.csv` (every classified row).

Isolate the effect of compounding at matched `k`:

```bash
python -m cns_exposure compare-models --therapy ... --rt ...
```

Check the dictionary:

```bash
python -m cns_exposure validate
```

## The swimmer plot

A presence bar states that an agent was on board. It is the *output* of the
pharmacokinetic model, with the model itself invisible — you cannot tell a drug
that barely cleared the threshold from one sitting an order of magnitude above
it, and you cannot see accumulation happening at all. The swimmer plot draws the
model underneath the answer, one lane per patient-drug record:

```bash
pip install -r requirements-plot.txt

python -m cns_exposure swimmer \
    --therapy therapy_records.csv \
    --rt rt_records.csv \
    --administrations administrations.csv \
    --k 2 --model compounding \
    --deidentify \
    --out swimmer.png --lanes-csv lanes.csv
```

The blue sawtooth is `C(t)`, in units of one single-dose peak. Each
administration steps it up, it decays between doses, it climbs toward the
plateau `R`, and after the last dose it decays from `R·D`. **The lane's dashed
baseline is the presence threshold `(1/2)^k`**, so "curve above the baseline"
and "presence bar drawn" are the same statement, and the point where the tail
crosses the baseline is the moment the drug stops counting as present.

The vertical axis is log₂ within each lane. A linear axis cannot show this: at
`k = 5` the threshold is 0.031 while an agent dosed daily against a 4.5-day
half-life plateaus near 7, which pins the whole decision boundary to the axis
floor. The floor sits 2.5 log₂ units *below* the threshold so a decaying tail is
seen crossing it; below that the line simply ends rather than running flat, which
would read as a constant low concentration that was never measured.

`--lanes-csv` writes every number on the figure so a lane can be traced back.
`--deidentify` replaces patient identifiers with stable `P01`-style labels
assigned in sorted order — **figures generated without it print whatever
`patient_id` you supplied and are not publishable.**

Expect the two accumulation models to look similar *on this figure* and differ
sharply under `compare-models`. Compounding lengthens the tail after the last
dose; inside the radiotherapy window the drug is usually still being given, so
the coverage rule already counts those days and accumulation has nothing to add.
That is a documented negative result, not a defect — and it is the reason
compounding and the coverage rule are reported as separate axes.

Plotting needs matplotlib and numpy. The core classifier depends only on pandas,
so a pipeline producing exposure variables and no figures does not carry a
plotting stack; importing `cns_exposure.plot` without them raises with the
install command rather than a bare `ModuleNotFoundError`.

From Python:

```python
from cns_exposure import DrugDictionary, classify, sweep

d = DrugDictionary.load()
out = classify(therapy, rt, d, k=5, rule="window", model="compounding")
grid = sweep(therapy, rt, d, administrations=admins)
grid["stability"].query("verdict == 'unstable'")
```

## Inputs

| File | Columns |
|---|---|
| `therapy_records.csv` | `patient_id`, `drug_or_regimen`, `start_date`, `end_date` |
| `rt_records.csv` | `patient_id`, `rt_start_date`, `rt_end_date` |
| `administrations.csv` | `patient_id`, `drug_or_regimen`, `admin_date` — `coverage` rule only |

Dates in ISO 8601. A record the rules cannot decide — a blank `end_date`, an
unreadable RT interval — returns **missing** exposure values, not zeros, and is
counted in `n_unresolved`. Exposure is aggregated by disjunction, so a `False`
that really means "unknown" is absorbed by any other record and biases exposure
downward leaving no trace; treat a unit carrying unresolved records as missing
rather than unexposed until they are reviewed. A drug absent from the dictionary
is returned as an unmatched class-0 row rather than dropped, so an unrecognised
name is visible in the output instead of quietly reducing the exposure count.
Patients without administration dates keep
their `window` classification and are marked, because dropping them would
restrict the cohort to patients whose charts happened to be abstracted more
deeply — not a clinical criterion, and unlikely to be independent of anything
being measured.

Output is one row per resolved agent per patient, carrying its full definition
cell. Aggregate to patient, lesion or treatment-course level with `aggregate()`;
the correct unit is a property of the study design, not of the exposure
definition.

## The dictionary

Five tables in [`dictionary/`](dictionary/), documented in
[`dictionary/README.md`](dictionary/README.md): 114 scored agents with
half-lives and evidence grades, 114 aliases, half-life provenance, 6
combination overrides, and dosing intervals with accumulation factors for 37
agents.

Point `--dictionary` at your own directory to substitute one. The framework
assumes the dictionary is **declared before analysis and not tuned afterwards
against outcome results**; nothing in the code can enforce that.

## Limitations

- **The dictionary reflects the literature at its review date.** CNS activity
  evidence moves. Re-review before reuse and record the version.
- **`τ` in the shipped dosing table came from one derivation cohort**, and one
  entry (trastuzumab) rests on a single non-standard course. Substituting label
  schedules is usually better. See [`dictionary/README.md`](dictionary/README.md).
- **The sweep is not a selection procedure.** Using it to find the cell with
  the cleanest result is a multiple-comparisons machine, and nothing chosen
  that way is confirmatory. Report the distribution.
- **Day-level presence in the `coverage` rule is a deliberate coarsening.** At
  clock resolution a 45-minute-half-life oral agent given daily is present
  about 9% of the time, which would make every short-half-life oral drug
  non-concurrent by construction and turn the rule into a proxy for route.
- **This is an exposure classifier, not a causal design.** Implementation
  audits and endpoint-specific sensitivity analyses remain necessary before use.

## Data

No patient data is in this repository and none should be added to it. The
worked example in [`examples/`](examples/) is the fictional dataset published
with the preprint; the dictionary tables are literature-derived reference
values. The package reads clinical inputs from paths you supply.

## Citing

Mitre LP, Drapkin B, Dohopolski M. A Drug-Specific, Half-Life-Adjusted
Framework for Classifying CNS-Active Systemic Therapy Exposure During and After
Radiotherapy. medRxiv. 2026. doi:10.64898/2026.06.11.26354463

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff). The compounding
model, the coverage rule and the sweep are extensions to the preprint; see
[`CHANGELOG.md`](CHANGELOG.md).

## Licence

CC BY 4.0 — see [LICENSE](LICENSE). Attribution required; commercial and
derivative use permitted.
