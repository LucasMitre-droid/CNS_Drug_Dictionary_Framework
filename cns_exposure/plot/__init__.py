"""
Figures. Optional: needs matplotlib and numpy, which the core does not.

    pip install -r requirements-plot.txt
    # or
    pip install "cns-exposure[plot]"

The core classifier depends only on pandas, so a pipeline that produces
exposure variables and no figures does not carry a plotting stack. Importing
this subpackage without matplotlib raises with that instruction rather than a
bare ModuleNotFoundError.
"""

from __future__ import annotations

_MISSING = (
    "Figures need matplotlib and numpy, which are optional extras.\n"
    "  pip install -r requirements-plot.txt\n"
    '  or: pip install "cns-exposure[plot]"'
)

try:
    from .swimmer import (
        SwimmerRecord,
        build_records,
        concentration_curve,
        presence_spans,
        render,
        to_frame,
    )
except ModuleNotFoundError as exc:  # pragma: no cover - depends on install
    raise ModuleNotFoundError(_MISSING) from exc

__all__ = [
    "SwimmerRecord", "build_records", "render", "to_frame",
    "concentration_curve", "presence_spans",
]
