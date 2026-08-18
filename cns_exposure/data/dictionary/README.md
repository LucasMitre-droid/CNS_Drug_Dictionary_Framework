# Reference tables

Five CSVs. Together they are the framework's declared input: fix them before
analysis, record the version, and do not tune them afterwards against outcome
results. Nothing in the code can enforce that discipline.

None of these files contains patient data. `S1`–`S4` are literature-derived and
published as the preprint supplement. `S5` is a per-agent regimen property; see
the provenance note below.

## S1_drug_dictionary.csv

112 agents, one row each. The scored dictionary.

| Column | Meaning |
|---|---|
| `Canonical_Name` | the name every lookup resolves to |
| `Route` | O oral, I intravenous, SC subcutaneous, M intramuscular |
| `iORR_Tier_Derived` | `<20%`, `20-39%`, `40-59%`, `>=60%` |
| `Intrathecal_Flag` | clinically relevant CNS use is intrathecal |
| `pts_iORR`, `pts_Consensus`, `pts_IT` | evidence component points |
| `CNS_Score_Total` | their sum |
| `CNS_Class`, `CNS_Class_Label` | 0–3, from the score |
| `Evidence_Grade_G` | A–E, strength of the supporting literature |
| `Half_Life_Value`, `Half_Life_Unit` | terminal elimination half-life, days |
| `Dictionary_Version`, `Date_Last_Reviewed` | provenance |
| `Manual_Notes` | curation rationale, with PMIDs |

Read `Manual_Notes` before trusting a borderline agent. Several rows record a
correction against an earlier version, an iORR that turned out to come from a
combination rather than monotherapy, or a class retained on a consensus point
alone.

Class distribution: 56 class 0, 19 class 1, 26 class 2, 11 class 3. Nine
agents carry no half-life; for those the effective window falls back to the
recorded stop date, so they can only ever be under-counted as exposure. All
nine are class 0 — see the known issues below.

## S2_alias_normalization.csv

Brand names, abbreviations and spelling variants → canonical name. Lookup also
strips hyphens and spaces, so `trastuzumab-deruxtecan`, `trastuzumab
deruxtecan` and `trastuzumabderuxtecan` resolve without an alias row.

## S3_halflife_sources.csv

Half-life provenance, one row per agent: value, source type, citation, and the
DailyMed URL for FDA-label values. Where a range is reported, the selected
value is the one in `S1`.

## S4_combination_overrides.csv

Six regimens whose combination intracranial response evidence exceeds any of
their components alone. A combination matched here is scored as a unit and
takes `h = max(components)`, the component governing how long the regimen
persists. Combinations are matched longest-first, so a three-agent override is
tested before a two-agent override contained in it.

## S5_dosing_intervals.csv

37 agents: dosing interval `τ`, half-life `h`, accumulation factor `R`, the
half-lives it adds (`log₂R`), and the resulting days added to the window.

```
R = 1 / (1 − (1/2)^(τ/h))
```

`R` here is exactly that formula applied to the `Tau_Days` and `Half_Life_Days`
in the same row, so the table is reproducible from its own columns. The
finite-course cap — a drug given twice has not reached the plateau of an
infinite series — depends on how many doses a particular course contained and
therefore cannot live in a per-agent table. It is available in the API as
`accumulation_factor(h, tau, n_doses)`.

`Tau_Confidence` is `high` when three or more courses contributed to the
agent's interval and `low` otherwise.

### Provenance, and why you may want to replace this table

`τ` was reconstructed as the modal interval between consecutive
administrations in the framework's derivation cohort, then taken per agent as a
median. `R` is a property of the **regimen** — a drug's half-life measured
against its dosing schedule — not of a patient, which is what allows a
per-agent value to generalise to patients whose administration dates were never
abstracted.

Two cautions:

- **trastuzumab** rests on a single course at `τ = 10 d`, which is not a
  standard schedule (q21d IV, or q7d weekly). It carries the largest tail
  extension in the table at 61 days. The label-schedule alternative `τ = 21`
  gives `R = 2.34` and 38 days; `R = 1.0` excludes it. Override before use.
- A reconstructed cadence belongs to a therapy **record**, not to each
  component of a combination. Attributing one row's cadence to every component
  is safe for `nivolumab + ipilimumab`, both three-weekly antibodies, and
  produces nonsense for `tucatinib + trastuzumab + capecitabine`, where an oral
  daily cadence would be inherited by an antibody actually dosed q21d — one
  substitution that inflates `R` by an order of magnitude. The API guards this
  with `tau_attributable(half_lives)`, which attributes a cadence only when
  `max(h)/min(h) ≤ 3`.

Substituting label schedules is usually the better input, and rebuilding `τ`
from your own administration records is better still. `reconstruct_tau(dates)`
does the reconstruction; `tau_attributable(...)` applies the guard.

## Known issues in the shipped dictionary

`python -m cns_exposure validate` reports these. They are left in place rather
than silently corrected, because the dictionary is a record of what was
declared and a validator that rewrites its input destroys that record.

`v3` reports **no errors**. Two warnings stand, and neither changes the class of
any CNS-active agent.

| Severity | Item | Issue |
|---|---|---|
| warning | 9 agents | no usable half-life. Every one is class 0, so none contributes to CNS-active exposure and the gap cannot move an exposure call. Two have no value to record: `domvanalimab` is investigational with no FDA label, and `quavonlimab` is unresolved manual review. The other seven — `cabazitaxel`, `fludarabine`, `hyaluronidase`, `leuprorelin`, `mitomycin`, `relatlimab`, `tocilizumab` — name a specific FDA label in `S3` whose value was never transcribed into `S1`. That is a transcription gap, not an evidence gap, and it is the one item here worth closing. |
| warning | `tucatinib + trastuzumab + capecitabine` | override class 2 does not exceed its highest component class, so the override changes nothing. Deliberate: `tucatinib` already reaches class 2 alone, and the row exists so that the HER2CLIMB regimen is matched and labelled as a unit rather than silently decomposing into components. The warning is the validator correctly reporting a redundancy that is there on purpose. |

Two encodings of a missing half-life appear in `S1` and they are not
interchangeable: `-` means no such value exists — the agent is investigational,
or its half-life is not applicable — while an empty cell means the value was not
recorded. Both coerce to missing, so the classifier treats them alike; the
distinction survives only for a reader deciding whether a gap is worth chasing.

### Resolved in v3

Earlier versions of this table carried three errors. All three are corrected in
the shipped `v3` files, and the corrections are recorded here rather than
dropped, because a known-issues list that quietly loses entries is no more of a
record than a validator that rewrites its input.

| Item | Was | Now |
|---|---|---|
| `durvalumab` | `CNS_Score_Total` 1 with all three component columns 0 | `pts_Consensus = 1`, which is the consensus point `Manual_Notes` always described. Class 1, unchanged. |
| `irinotecan` | same pattern, with `Manual_Notes` recording an open question over `pts_iORR` | question closed against the source: the 60%/28%/37% figures were bevacizumab + irinotecan combination data, and solo iORR is 6.3% (PMID 19066728). All components 0, total 0, **class 0** — irinotecan is out of the CNS-active set. |
| `krazati` | alias resolved to `adagrasib`, which had no `S1` row | alias removed. Neither name appears in any table, so `krazati` now resolves as an unmatched class-0 row carrying its own name, which is visible in the output. Adding a scored `adagrasib` row remains open work. |

## Substituting your own

Point `--dictionary` at a directory holding files with these names. Only
`S1_drug_dictionary.csv` is required; the rest degrade gracefully — without
`S2` no aliases resolve, without `S4` combinations decompose to components,
without `S5` every `R` is 1.0 and the compounding model reduces to the
published single-dose model.

Run `python -m cns_exposure validate --dictionary <dir>` before using one.
