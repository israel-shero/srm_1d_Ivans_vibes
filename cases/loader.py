"""Load registered validation cases without guessing file schemas."""

from __future__ import annotations

import copy
import json
from pathlib import Path

from srm_1d import plotting


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = Path(__file__).with_name("registry.json")


class CaseNotFoundError(KeyError):
    """Raised when a requested case ID is not registered."""


class MeasurementUnavailableError(RuntimeError):
    """Raised when registered measurement data are unsafe or unavailable."""


def load_registry(path: str | Path = DEFAULT_REGISTRY) -> dict:
    """Return the parsed case registry."""
    registry_path = Path(path)
    with registry_path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def get_case(case_id: str, registry: dict | None = None) -> dict:
    """Return one case by stable ID."""
    source = registry if registry is not None else load_registry()
    for case in source["cases"]:
        if case["id"] == case_id:
            return copy.deepcopy(case)
    raise CaseNotFoundError(case_id)


def _file_for_role(case: dict, role: str) -> Path:
    for entry in case["files"]:
        if entry["role"] == role and entry.get("available", True):
            return PROJECT_ROOT / entry["path"]
    raise MeasurementUnavailableError(
        f"{case['id']}: no available file registered for role {role!r}"
    )


def load_measurement(case_id: str, registry: dict | None = None) -> dict:
    """Load a registry-approved processed trace in seconds and MPa.

    Raw, missing, and absent measurements deliberately raise instead of being
    calibrated, aligned, or scaled implicitly.
    """
    case = get_case(case_id, registry)
    measurement = case["measurement"]
    kind = measurement["kind"]

    if kind == "csv":
        path = _file_for_role(case, "measurement_csv")
        return plotting.load_experimental_csv(
            str(path),
            time_col=measurement["time_column"],
            pressure_col=measurement["pressure_column"],
            delimiter=measurement["delimiter"],
            skip_header=measurement["skip_header"],
            pressure_unit=measurement["pressure_unit"],
            label=case["label"],
        )

    if kind == "python_symbol":
        symbol = measurement["symbol"]
        try:
            trace = getattr(plotting, symbol)
        except AttributeError as exc:
            raise MeasurementUnavailableError(
                f"{case_id}: plotting symbol {symbol!r} is missing"
            ) from exc
        return {
            "time": trace["time"].copy(),
            "pressure": trace["pressure"].copy(),
            "label": trace.get("label", case["label"]),
        }

    reason = {
        "raw_csv": "raw trace is not approved for direct loading",
        "missing": "measurement file is missing",
        "none": "no measured trace is registered",
    }.get(kind, f"unsupported measurement kind {kind!r}")
    raise MeasurementUnavailableError(f"{case_id}: {reason}")
