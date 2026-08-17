"""Shared validation and provenance helpers for the command-line workflow."""

from __future__ import annotations

import json
import platform
from collections.abc import Iterable
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse


def get_counts(adata, layer: str = "counts"):
    """Return the requested count matrix and reject invalid count-like data."""
    matrix = adata.layers.get(layer, adata.X)
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if values.size and (not np.isfinite(values).all() or values.min() < 0):
        raise ValueError("Counts must be finite and non-negative.")
    if values.size and not np.allclose(values, np.rint(values), atol=1e-6):
        raise ValueError(
            f"The {layer!r} matrix is not integer-like. Supply raw UMI counts, not normalized data."
        )
    return matrix


def require_columns(adata, columns: Iterable[str]) -> None:
    missing = [column for column in columns if column not in adata.obs]
    if missing:
        raise KeyError(f"Missing required adata.obs columns: {', '.join(missing)}")


def package_versions(packages: Iterable[str]) -> dict[str, str]:
    found: dict[str, str] = {"python": platform.python_version()}
    for package in packages:
        try:
            found[package] = version(package)
        except PackageNotFoundError:
            found[package] = "not-installed"
    return found


def record_provenance(adata, stage: str, parameters: dict[str, Any], packages: Iterable[str]) -> None:
    pipeline = adata.uns.setdefault("scrna_workflow", {})
    existing = pipeline.get("runs", [])
    if isinstance(existing, np.ndarray):
        runs = existing.astype(str).tolist()
    elif isinstance(existing, str):
        runs = [existing]
    else:
        runs = list(existing)
    runs.append(
        json.dumps(
            {
                "stage": stage,
                "completed_at_utc": datetime.now(timezone.utc).isoformat(),
                "parameters": _json_safe(parameters),
                "versions": package_versions(packages),
            },
            sort_keys=True,
        )
    )
    pipeline["runs"] = runs


def write_json(path: str | Path, payload: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(_json_safe(payload), indent=2, sort_keys=True) + "\n")


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, Path):
        return str(value)
    return value
