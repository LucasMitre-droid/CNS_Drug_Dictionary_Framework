"""
Reproducible CNS-active systemic therapy exposure variables around radiotherapy.

    from cns_exposure import DrugDictionary, classify, sweep

    d = DrugDictionary.load()
    out = classify(therapy, rt, d, k=5, rule="window", model="compounding")
"""

from .classify import aggregate, classify, read_csv
from .dictionary import DrugDictionary, DrugEntry, normalize
from .ksweep import compare_models, stability, summarize, sweep
from .pharmacokinetics import (
    BOX,
    COMPOUNDING,
    accumulation_factor,
    concentration,
    effective_window_end,
    extra_half_lives,
    is_present,
    reconstruct_tau,
    tau_attributable,
)
from .rules import COVERAGE, WINDOW, classify_coverage, classify_window
from .scoring import cns_class, cns_score, is_cns_active

__version__ = "4.0.0"

__all__ = [
    "DrugDictionary", "DrugEntry", "normalize",
    "classify", "aggregate", "read_csv",
    "sweep", "summarize", "stability", "compare_models",
    "accumulation_factor", "extra_half_lives", "effective_window_end",
    "concentration", "is_present", "reconstruct_tau", "tau_attributable",
    "BOX", "COMPOUNDING",
    "classify_window", "classify_coverage", "WINDOW", "COVERAGE",
    "cns_score", "cns_class", "is_cns_active",
    "__version__",
]
