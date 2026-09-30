#!/usr/bin/env python3
"""Verify the pinned srm_1d/openMotor development workspace."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LOCK = PROJECT_ROOT / "requirements" / "workspace-lock.json"


def run(command: list[str], cwd: Path) -> tuple[int, str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return completed.returncode, completed.stdout.strip()


def normalized_remote(value: str) -> str:
    """Normalize HTTPS and SSH GitHub remotes for stable comparison."""
    value = value.strip().removesuffix(".git")
    if value.startswith("git@github.com:"):
        return value.removeprefix("git@github.com:").lower()
    parsed = urlparse(value)
    if parsed.hostname == "github.com":
        return parsed.path.strip("/").lower()
    return value.lower()


def check_repository(
    name: str,
    config: dict[str, object],
    strict_clean: bool,
) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    repo = (PROJECT_ROOT / str(config["path"])).resolve()

    if not repo.is_dir():
        return [f"{name}: missing checkout at {repo}"], warnings

    code, actual_commit = run(["git", "rev-parse", "HEAD"], repo)
    if code != 0:
        errors.append(f"{name}: not a Git checkout at {repo}: {actual_commit}")
        return errors, warnings

    expected_commit = str(config["commit"])
    pin_mode = str(config.get("pin_mode", "exact"))
    if pin_mode == "exact":
        if actual_commit != expected_commit:
            errors.append(
                f"{name}: commit {actual_commit} does not match {expected_commit}"
            )
    elif pin_mode == "contains":
        code, _ = run(
            ["git", "merge-base", "--is-ancestor", expected_commit, actual_commit],
            repo,
        )
        if code != 0:
            errors.append(
                f"{name}: HEAD {actual_commit} does not contain baseline "
                f"{expected_commit}"
            )
    else:
        errors.append(f"{name}: unsupported pin_mode {pin_mode!r}")

    code, actual_origin = run(["git", "remote", "get-url", "origin"], repo)
    if code != 0:
        errors.append(f"{name}: cannot read origin remote: {actual_origin}")
    elif normalized_remote(actual_origin) != normalized_remote(str(config["origin"])):
        errors.append(
            f"{name}: origin {actual_origin!r} does not match {config['origin']!r}"
        )

    for required_path in config.get("required_paths", []):
        if not (repo / str(required_path)).is_file():
            errors.append(f"{name}: missing required file {required_path}")

    code, status = run(["git", "status", "--short"], repo)
    if code != 0:
        errors.append(f"{name}: cannot read worktree status: {status}")
    elif status:
        message = f"{name}: worktree has local changes:\n{status}"
        (errors if strict_clean else warnings).append(message)

    return errors, warnings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lock",
        type=Path,
        default=DEFAULT_LOCK,
        help="workspace lock JSON (default: requirements/workspace-lock.json)",
    )
    parser.add_argument(
        "--strict-clean",
        action="store_true",
        help="treat local repository changes as verification errors",
    )
    args = parser.parse_args()

    lock_path = args.lock.resolve()
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: cannot read workspace lock {lock_path}: {exc}")
        return 2

    errors: list[str] = []
    warnings: list[str] = []

    expected_python = str(lock["python"]["major_minor"])
    actual_python = f"{sys.version_info.major}.{sys.version_info.minor}"
    if actual_python != expected_python:
        errors.append(
            f"python: running {actual_python}, expected {expected_python} "
            f"({sys.executable})"
        )

    for name, config in lock["repositories"].items():
        repo_errors, repo_warnings = check_repository(
            name, config, args.strict_clean
        )
        errors.extend(repo_errors)
        warnings.extend(repo_warnings)

    code, pip_output = run([sys.executable, "-m", "pip", "check"], PROJECT_ROOT)
    if code != 0:
        errors.append(f"pip check failed:\n{pip_output}")

    print(f"workspace lock: {lock_path}")
    print(f"python: {sys.version.split()[0]} ({sys.executable})")
    for warning in warnings:
        print(f"WARNING: {warning}")
    for error in errors:
        print(f"ERROR: {error}")

    if errors:
        print(f"FAILED: {len(errors)} error(s), {len(warnings)} warning(s)")
        return 1

    print(f"PASS: pinned workspace verified ({len(warnings)} warning(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
