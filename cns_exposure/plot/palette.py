"""
Colour and typography.

Status is carried three ways -- by colour, by a filled or hollow glyph, and by a
printed percentage on every lane -- so no reading of the figure depends on hue
alone.

The concurrent / not-concurrent pair measures ΔE 20.2 under normal vision and
ΔE 10.7 under simulated deuteranopia, both at contrast ≥ 3:1 against the
surface. The neutral is deliberately low-chroma: it encodes "no signal" rather
than a second category. Green against red would measure ΔE 4.1 under
deuteranopia, so reference lines are drawn in neutral ink instead.

The concentration curve gets its own hue, distinct from both status colours,
because it is the model rather than the verdict.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

CALLED = "#0ca30c"
CALLED_FAINT = "#b6e3b6"
NOT_CALLED = "#898781"
NOT_CALLED_FAINT = "#dedcd6"

RT_BAND = "#cde2fb"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
ROW_BAND = "#f4f3ef"

CURVE = "#1f4e79"
THRESHOLD = "#8a6d3b"

STATUS_GLYPH = {True: "●", False: "○"}


def fill_for(called: bool, inside_window: bool) -> str:
    """Bar colour for a span, given its status and whether it lies inside RT."""
    if inside_window:
        return CALLED if called else NOT_CALLED
    return CALLED_FAINT if called else NOT_CALLED_FAINT


def apply_rc() -> None:
    """Install the figure-wide style. Idempotent."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Segoe UI", "Arial"],
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "text.color": INK,
        "axes.labelcolor": INK_2,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.titlecolor": INK,
        "savefig.facecolor": SURFACE,
    })
