from __future__ import annotations

from contextlib import ExitStack
from importlib import resources
from pathlib import Path


_PACKAGE = "cns_exposure"
_ROOT = Path(__file__).resolve().parent.parent


def repo_root() -> Path:
    return _ROOT


def bundled_path(*parts: str) -> Path:
    """Materialize a bundled resource path for the duration of the process."""
    traversable = resources.files(_PACKAGE).joinpath(*parts)
    if not traversable.exists():
        raise FileNotFoundError(f"Bundled resource not found: {'/'.join(parts)}")
    return _RESOURCE_STACK.enter_context(resources.as_file(traversable))


def prefer_repo_path(*parts: str) -> Path:
    """Use the repo copy when present, otherwise fall back to installed data."""
    candidate = _ROOT.joinpath(*parts)
    return candidate if candidate.exists() else bundled_path(*parts)


_RESOURCE_STACK = ExitStack()
