#!/usr/bin/env python3
"""Validate case-registry structure, paths, checksums, and CSV row counts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = PROJECT_ROOT / "cases" / "registry.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def csv_data_rows(path: Path, skip_header: int) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = sum(1 for _ in csv.reader(stream))
    return max(0, rows - skip_header)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    args = parser.parse_args()

    registry_path = args.registry.resolve()
    try:
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: cannot read registry {registry_path}: {exc}")
        return 2

    errors: list[str] = []
    warnings: list[str] = []
    known_readiness = set(registry.get("readiness_values", []))
    known_flags = set(registry.get("quality_flags", {}))
    seen_ids: set[str] = set()

    for case in registry.get("cases", []):
        case_id = case.get("id")
        if not isinstance(case_id, str) or not case_id:
            errors.append("case has missing or invalid id")
            continue
        if case_id in seen_ids:
            errors.append(f"{case_id}: duplicate case id")
        seen_ids.add(case_id)

        if case.get("readiness") not in known_readiness:
            errors.append(f"{case_id}: unknown readiness {case.get('readiness')!r}")

        for flag in case.get("quality_flags", []):
            if flag not in known_flags:
                errors.append(f"{case_id}: unknown quality flag {flag!r}")

        roles: set[str] = set()
        for entry in case.get("files", []):
            role = entry.get("role")
            if not isinstance(role, str) or not role:
                errors.append(f"{case_id}: file entry has no role")
                continue
            if role in roles:
                errors.append(f"{case_id}: duplicate file role {role!r}")
            roles.add(role)

            relative = entry.get("path")
            if not isinstance(relative, str) or not relative:
                errors.append(f"{case_id}/{role}: missing path")
                continue
            path = PROJECT_ROOT / relative
            available = entry.get("available", True)

            if not available:
                if path.exists():
                    errors.append(
                        f"{case_id}/{role}: registry says unavailable but file exists: "
                        f"{relative}"
                    )
                continue

            if not path.is_file():
                errors.append(f"{case_id}/{role}: missing file {relative}")
                continue

            expected_hash = entry.get("sha256")
            if not isinstance(expected_hash, str) or len(expected_hash) != 64:
                errors.append(f"{case_id}/{role}: missing or invalid SHA-256")
            else:
                actual_hash = sha256(path)
                if actual_hash != expected_hash:
                    errors.append(
                        f"{case_id}/{role}: checksum changed for {relative}: "
                        f"{actual_hash}"
                    )

            if "data_rows" in entry:
                skip_header = int(entry.get("skip_header", 0))
                actual_rows = csv_data_rows(path, skip_header)
                if actual_rows != entry["data_rows"]:
                    errors.append(
                        f"{case_id}/{role}: expected {entry['data_rows']} data rows, "
                        f"found {actual_rows}"
                    )

        if "motor" not in roles:
            errors.append(f"{case_id}: no motor file registered")
        if "runner" not in roles:
            warnings.append(f"{case_id}: no runner file registered")

    print(f"registry: {registry_path}")
    print(f"cases: {len(seen_ids)}")
    for warning in warnings:
        print(f"WARNING: {warning}")
    for error in errors:
        print(f"ERROR: {error}")

    if errors:
        print(f"FAILED: {len(errors)} error(s), {len(warnings)} warning(s)")
        return 1
    print(f"PASS: case registry verified ({len(warnings)} warning(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
