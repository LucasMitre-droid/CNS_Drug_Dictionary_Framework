# Changelog

## 4.0.0

Extends the framework as published in the preprint
(doi:10.64898/2026.06.11.26354463). The published behaviour is preserved as the
default: `--rule window --model box --k 5` reproduces the original classifier.

### Added

- **Compounding half-life.** `e_EF = e + (k + log₂R)·h`, with
  `R = 1 / (1 − (1/2)^(τ/h))`. `log₂R` is independent of `k`, so compounding is
  an axis orthogonal to the `k` sweep rather than a rescaling of it. Select
  with `--model compounding`.
- **Coverage rule.** Dose-level concurrency from superposed concentration
  `C(t) = Σᵢ (1/2)^((t − aᵢ)/h)`, concurrent when the agent is on board on at
  least 50% of radiotherapy days. Requires administration dates. Reduces
  exactly to the window rule for a single dose, which is what makes the two
  comparable at a shared `k`.
- **Definition grid sweep.** `sweep()` over rule × `k` × accumulation model,
  with per-unit `stability` reporting how often each analysis unit is called
  exposed across cells.
- **`compare-models`.** Single-dose against compounding at matched `k`,
  reporting both the window extension (invariant in `k`) and the count of
  reclassified units (not invariant, and not monotone).
- **Swimmer plot with the concentration curve drawn on it**
  (`cns_exposure.plot`, `python -m cns_exposure swimmer`). One lane per
  patient-drug record, with `C(t)` drawn on a per-lane log₂ axis mapped so the
  lane baseline is exactly the presence threshold `(1/2)^k` — making "curve
  above baseline" and "presence bar drawn" the same statement. Sampling is
  refined per lane to the agent's own half-life, so short-half-life agents are
  not aliased away by a fixed grid. Values below the plotted floor are masked
  rather than clipped flat. Carries `--deidentify` and a `--lanes-csv` audit
  export. matplotlib and numpy are optional extras; the core stays pandas-only.
- **Dictionary validator.** Checks score-component sums, class-from-score,
  iORR tier consistency, evidence grades, alias targets, combination
  components, and that each `R` is reproducible from its own `τ` and `h`.
  Reports; never repairs.
- **`S5_dosing_intervals.csv`.** Dosing interval and accumulation factor for 37
  agents.
- **Combination cadence guard.** `tau_attributable()` attributes a record's
  reconstructed cadence to its components only when their half-lives sit within
  a factor of 3, preventing an oral daily cadence from being inherited by a
  three-weekly antibody.
- **Configurable activity threshold.** `--activity-threshold` selects the
  lowest CNS class counting as active. Default 1, matching the publication.
- **Bounded adjuvant attribution.** `--attribution-window` limits how long
  after the last fraction a newly started drug counts as adjuvant initiation.
  Unbounded by default.

### Changed

- Adjuvant exposure is returned split into `aCNS_aD_persistence` (the window
  carried past the last fraction) and `aCNS_aD_initiation` (therapy started
  after it), alongside their union. The published implementation returned only
  the union.
- Records the rules cannot decide return missing exposure values and are
  counted in a new `unresolved` column, propagated to `n_unresolved` at the
  analysis unit. Exposure aggregates by disjunction, so a zero standing in for
  an unknown is absorbed by any other record and biases exposure downward
  without trace. The published implementation already emitted a missing value
  for the manual-review case; this generalises the behaviour and makes it
  countable.
- Unmatched drugs are returned as class-0 rows carrying the unrecognised name
  rather than being dropped, so a name the dictionary does not know is visible
  instead of quietly reducing the exposure count.
- Every output row carries its definition cell — rule, `k`, accumulation model,
  activity threshold — since a classified table without its own definition
  cannot be compared to another one.
- Name normalization strips hyphens and non-breaking spaces, so published
  spelling variants of the same agent resolve to one entry.
- `DrugEntry` carries `components` for a combination scored as a unit, so
  callers joining against per-drug tables — the swimmer plot's administration
  dates, for one — can still reach the agents the combination stands for.
- Superposed concentration is evaluated by direct summation rather than by a
  day-to-day recurrence, whose downward rounding drift drops a day of presence
  at exact ties — the boundary at which the two accumulation models must agree
  by construction.

### Fixed in the dictionary

Three errors the validator reported against the published dictionary are
corrected in the shipped `v3` tables. Each was corrected in the data, against
its source, rather than by relaxing the check that caught it.

- **`durvalumab`** recorded a total score of 1 with all three component point
  columns 0. `pts_Consensus` now reads 1, which is the consensus point
  `Manual_Notes` described all along. Class 1, unchanged.
- **`irinotecan`** recorded the same pattern, with an open question over
  `pts_iORR`. Closed against the source: the iORR figures in the earlier table
  were bevacizumab + irinotecan combination data, and solo iORR is 6.3%
  (PMID 19066728). Scored 0 on every component, **class 0** — irinotecan leaves
  the CNS-active set, which is the one correction here that changes an exposure
  variable.
- **`krazati`** aliased to `adagrasib`, which had no row in S1. The alias is
  removed, so the name now surfaces as an unmatched class-0 row instead of
  resolving to a canonical name the dictionary cannot score. A scored
  `adagrasib` row is still to be added.

`validate` reports no errors against `v3`. The two remaining warnings — nine
class-0 agents with no half-life, and one deliberately redundant combination
override — are documented in `dictionary/README.md`.

### Removed

- **`NVL-655` and `WTXX-124`** are dropped from `S1` and `S3`, taking the
  dictionary from 114 scored agents to 112 and class 0 from 58 to 56. Both were
  investigational agents carried under a trial code rather than an INN, with no
  FDA label and no half-life. Neither had an alias in `S2` or a dosing row in
  `S5`, so nothing else referenced them. Exposure to either name is unaffected
  in class terms — both scored 0 — but it now resolves as an **unmatched**
  class-0 row carrying the unrecognised name rather than as a scored one. That
  is a real change in what the dictionary claims: a scored 0 asserts the agent
  was reviewed and found non-CNS-active, while an unmatched row asserts only
  that the dictionary does not know the name.

## 3.0.0

Dictionary version as published: 114 scored agents, single-dose exposure window
`e_EF = e + 5h`, overlap-based timing rules.
