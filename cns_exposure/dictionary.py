"""
The drug dictionary: loading, normalization, regimen decomposition, validation.

The dictionary is the framework's declared input. It must be fixed before
analysis and not tuned afterwards against outcome results; every rule here
assumes that discipline and none of them can enforce it.

Four tables, and one derived from the study's own records:

  S1  canonical name -> evidence points, CNS class, evidence grade, half-life
  S2  alias, brand name or abbreviation -> canonical name
  S3  half-life provenance, one row per agent
  S4  combination regimens whose joint CNS class exceeds their components'
  S5  dosing interval and accumulation factor per agent

S5 is the only table a new study is likely to need to rebuild, because tau
depends on how the agent was actually given. It ships with values from the
derivation cohort so the compounding model runs out of the box, but a
label-schedule table or a table reconstructed from your own administration
records is the better input. See `dictionary/README.md`.

REGIMEN DECOMPOSITION

Combination regimens are split into components, because CNS activity and
half-life are properties of a drug and not of a regimen label. A regimen
pairing one CNS-active targeted agent with a non-CNS-active cytotoxic backbone
must not collapse into a single unanalysed string, or the distinction the
framework exists to draw is lost at the first step.

Decomposition is checked against S4 first. A handful of combinations have
intracranial response evidence exceeding any of their components alone, and for
those the combination is scored as a unit and takes the longest component
half-life, which is the component that governs how long the regimen persists.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import pharmacokinetics as pk
from . import resources
from . import scoring

DICTIONARY_DIR = resources.prefer_repo_path("dictionary")

FILES = {
    "drugs": "S1_drug_dictionary.csv",
    "aliases": "S2_alias_normalization.csv",
    "halflife_sources": "S3_halflife_sources.csv",
    "combinations": "S4_combination_overrides.csv",
    "dosing": "S5_dosing_intervals.csv",
}

# Tokens that join components of a regimen. Matched case-insensitively.
DELIMITERS = re.compile(r"[/+,;&]|\bwith\b|\bplus\b|\band\b", re.IGNORECASE)


@dataclass(frozen=True)
class DrugEntry:
    """One resolved agent, or one combination scored as a unit."""
    name: str
    cns_class: int
    half_life_days: float | None
    evidence_grade: str | None = None
    tau_days: float | None = None
    accumulation: float = 1.0
    is_combination: bool = False
    matched: bool = True
    components: tuple[str, ...] = ()
    """Canonical component names. Empty for a single agent; populated for a
    combination scored as a unit, so callers joining on per-drug tables can
    still reach the components the combination stands for."""

    @property
    def extra_half_lives(self) -> float:
        return pk.extra_half_lives(self.accumulation)

    def is_cns_active(self, threshold: int = scoring.DEFAULT_ACTIVITY_THRESHOLD
                      ) -> bool:
        return scoring.is_cns_active(self.cns_class, threshold)


@dataclass
class DrugDictionary:
    entries: dict[str, DrugEntry] = field(default_factory=dict)
    aliases: dict[str, str] = field(default_factory=dict)
    combinations: list[dict] = field(default_factory=list)
    source_dir: Path | None = None

    # -- loading -----------------------------------------------------------
    @classmethod
    def load(cls, directory: str | Path | None = None) -> "DrugDictionary":
        d = Path(directory) if directory else DICTIONARY_DIR
        missing = [f for f in (FILES["drugs"],) if not (d / f).exists()]
        if missing:
            raise FileNotFoundError(
                f"{d / FILES['drugs']} not found. Point --dictionary at a "
                f"directory holding {', '.join(FILES.values())}.")

        drugs = pd.read_csv(d / FILES["drugs"])
        dosing = _read_optional(d / FILES["dosing"])
        aliases = _read_optional(d / FILES["aliases"])
        combos = _read_optional(d / FILES["combinations"])

        tau_by_name, r_by_name = {}, {}
        if dosing is not None:
            for _, row in dosing.iterrows():
                key = normalize(row.get("Canonical_Name"))
                if not key:
                    continue
                tau_by_name[key] = _num(row.get("Tau_Days"))
                r_by_name[key] = _num(row.get("R")) or 1.0

        entries: dict[str, DrugEntry] = {}
        for _, row in drugs.iterrows():
            key = normalize(row.get("Canonical_Name"))
            if not key or key in entries:
                continue
            entries[key] = DrugEntry(
                name=str(row.get("Canonical_Name")).strip(),
                cns_class=int(_num(row.get("CNS_Class")) or 0),
                half_life_days=_num(row.get("Half_Life_Value")),
                evidence_grade=_text(row.get("Evidence_Grade_G")),
                tau_days=tau_by_name.get(key),
                accumulation=r_by_name.get(key, 1.0),
            )

        alias_map: dict[str, str] = {}
        if aliases is not None:
            for _, row in aliases.iterrows():
                a, c = normalize(row.get("Alias")), normalize(row.get("Canonical_Name"))
                if a and c and a != c:
                    alias_map[a] = c

        combo_list: list[dict] = []
        if combos is not None:
            for _, row in combos.iterrows():
                comps = [normalize(p) for p in
                         str(row.get("Components", "")).split(";") if normalize(p)]
                if not comps:
                    continue
                combo_list.append({
                    "label": str(row.get("Combo_Label", "")).strip(),
                    "components": comps,
                    "cns_class": int(_num(row.get("Override_Class")) or 0),
                })
        # Longest combinations first, so a three-agent override is tested
        # before a two-agent override that is a subset of it.
        combo_list.sort(key=lambda c: len(c["components"]), reverse=True)

        return cls(entries=entries, aliases=alias_map,
                   combinations=combo_list, source_dir=d)

    # -- lookup ------------------------------------------------------------
    def resolve_name(self, token: str) -> str:
        key = normalize(token)
        return self.aliases.get(key, key)

    def get(self, token: str) -> DrugEntry | None:
        return self.entries.get(self.resolve_name(token))

    def split_regimen(self, regimen: str) -> list[str]:
        """Component tokens of a regimen string, canonicalised, order preserved."""
        parts = DELIMITERS.split(str(regimen or ""))
        out, seen = [], set()
        for p in parts:
            name = self.resolve_name(p)
            if name and name not in seen:
                seen.add(name)
                out.append(name)
        return out

    def resolve_regimen(self, regimen: str) -> list[DrugEntry]:
        """
        A regimen string as a list of scored entries.

        A combination listed in S4 resolves to one entry carrying the override
        class and the longest component half-life. Anything else resolves to
        one entry per component. Components not in the dictionary are returned
        as unmatched class-0 entries rather than dropped, so that a name the
        dictionary does not know is visible in the output instead of quietly
        reducing the exposure count.
        """
        tokens = self.split_regimen(regimen)
        if not tokens:
            return []

        present = set(tokens)
        for combo in self.combinations:
            if set(combo["components"]).issubset(present):
                members = [self.entries.get(c) for c in combo["components"]]
                halves = [m.half_life_days for m in members
                          if m is not None and m.half_life_days is not None]
                accs = [m.accumulation for m in members if m is not None]
                taus = [m.tau_days for m in members
                        if m is not None and m.tau_days is not None]
                combo_entry = DrugEntry(
                    name=combo["label"],
                    cns_class=combo["cns_class"],
                    half_life_days=max(halves) if halves else None,
                    evidence_grade=None,
                    tau_days=max(taus) if taus else None,
                    accumulation=max(accs) if accs else 1.0,
                    is_combination=True,
                    components=tuple(combo["components"]),
                )
                rest = [t for t in tokens if t not in set(combo["components"])]
                return [combo_entry] + [self._entry_for(t) for t in rest]

        return [self._entry_for(t) for t in tokens]

    def _entry_for(self, token: str) -> DrugEntry:
        e = self.entries.get(token)
        if e is not None:
            return e
        return DrugEntry(name=token, cns_class=0, half_life_days=None,
                         matched=False)

    # -- validation --------------------------------------------------------
    def validate(self) -> pd.DataFrame:
        """
        Check the dictionary against the framework's own rules.

        Reports rather than repairs. A dictionary is a declared scientific
        input; silently rewriting one to satisfy a consistency check would
        destroy the record of what was actually declared.
        """
        issues: list[dict] = []

        def add(severity, name, check, detail):
            issues.append({"severity": severity, "drug": name,
                           "check": check, "detail": detail})

        drugs = _read_optional(
            (self.source_dir or DICTIONARY_DIR) / FILES["drugs"])
        if drugs is not None:
            for _, row in drugs.iterrows():
                name = str(row.get("Canonical_Name", "")).strip()
                if not name:
                    continue
                pts_iorr = _num(row.get("pts_iORR")) or 0
                pts_cons = _num(row.get("pts_Consensus")) or 0
                pts_it = _num(row.get("pts_IT")) or 0
                total = _num(row.get("CNS_Score_Total"))
                klass = _num(row.get("CNS_Class"))
                grade = _text(row.get("Evidence_Grade_G"))

                component_sum = pts_iorr + pts_cons + pts_it
                if total is not None and component_sum != total:
                    add("error", name, "score_components_sum",
                        f"pts_iORR {pts_iorr:g} + pts_Consensus {pts_cons:g} + "
                        f"pts_IT {pts_it:g} = {component_sum:g}, but "
                        f"CNS_Score_Total is {total:g}")

                if total is not None and klass is not None:
                    expected = scoring.cns_class(total)
                    if int(klass) != expected:
                        add("error", name, "class_from_score",
                            f"score {total:g} maps to class {expected}, "
                            f"recorded class is {int(klass)}")

                tier = _text(row.get("iORR_Tier_Derived"))
                if tier is not None:
                    expected_pts = scoring.iorr_points(tier=tier)
                    if expected_pts != pts_iorr:
                        add("error", name, "iorr_tier_points",
                            f"tier {tier} scores {expected_pts}, "
                            f"pts_iORR is {pts_iorr:g}")

                if klass is not None and int(klass) >= 1:
                    if grade is None:
                        add("warning", name, "missing_evidence_grade",
                            f"class {int(klass)} carries no evidence grade")
                    elif grade not in scoring.EVIDENCE_GRADES:
                        add("error", name, "unknown_evidence_grade",
                            f"grade {grade!r} is not one of "
                            f"{'/'.join(scoring.EVIDENCE_GRADES)}")

                if _num(row.get("Half_Life_Value")) is None:
                    add("warning", name, "missing_half_life",
                        "no half-life: the effective window falls back to the "
                        "recorded end date, so this agent can only be "
                        "under-counted as exposure")

        for combo in self.combinations:
            unknown = [c for c in combo["components"] if c not in self.entries]
            if unknown:
                add("error", combo["label"], "combination_component_unknown",
                    f"component(s) not in the dictionary: {', '.join(unknown)}")
            members = [self.entries[c] for c in combo["components"]
                       if c in self.entries]
            if members and combo["cns_class"] <= max(m.cns_class for m in members):
                add("warning", combo["label"], "combination_override_redundant",
                    f"override class {combo['cns_class']} does not exceed the "
                    f"highest component class "
                    f"{max(m.cns_class for m in members)}")
            if not any(m.half_life_days is not None for m in members):
                add("warning", combo["label"], "combination_no_half_life",
                    "no component carries a half-life")

        for alias, canonical in sorted(self.aliases.items()):
            if canonical not in self.entries:
                add("error", alias, "alias_target_unknown",
                    f"alias resolves to {canonical!r}, which is not in S1")

        for key, entry in sorted(self.entries.items()):
            if entry.tau_days and entry.half_life_days:
                expected = pk.accumulation_factor(
                    entry.half_life_days, entry.tau_days)
                if abs(expected - entry.accumulation) > 1e-3:
                    add("error", entry.name, "accumulation_factor",
                        f"tau {entry.tau_days:g} and half-life "
                        f"{entry.half_life_days:g} give R = {expected:.4f}, "
                        f"recorded R is {entry.accumulation:.4f}")

        return pd.DataFrame(issues, columns=["severity", "drug", "check", "detail"])

    def summary(self) -> pd.DataFrame:
        rows = []
        for _, e in sorted(self.entries.items()):
            rows.append({
                "drug": e.name,
                "cns_class": e.cns_class,
                "class_label": scoring.class_label(e.cns_class),
                "evidence_grade": e.evidence_grade,
                "half_life_days": e.half_life_days,
                "tau_days": e.tau_days,
                "R": round(e.accumulation, 4),
                "extra_half_lives": round(e.extra_half_lives, 4),
            })
        return pd.DataFrame(rows)


def normalize(value) -> str:
    """
    Canonical form of a drug name: lowercase, whitespace collapsed, hyphens and
    non-breaking spaces removed.

    Hyphens go because published names disagree with each other about them:
    trastuzumab-deruxtecan, trastuzumab deruxtecan and trastuzumabderuxtecan
    are one agent, and a lookup that treats them as three will report a real
    exposure as an unmatched drug.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    s = str(value).replace(" ", " ").strip().lower()
    s = re.sub(r"\s+", " ", s)
    return s.replace("-", "").replace(" ", "")


def _read_optional(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path) if Path(path).exists() else None


def _num(v) -> float | None:
    n = pd.to_numeric(v, errors="coerce")
    return None if pd.isna(n) else float(n)


def _text(v) -> str | None:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v).strip()
    return s or None
