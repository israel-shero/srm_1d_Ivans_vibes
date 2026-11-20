"""Reproducibility identity helpers shared by verification runners."""

from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import platform
import sys
from pathlib import Path

SCRIPT_ROOT = Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _package_files(root: Path) -> list[Path]:
    """Return source and runtime-data files, excluding generated/cache files."""
    return sorted(
        path for path in root.rglob("*")
        if path.is_file()
        and ".git" not in path.parts
        and "__pycache__" not in path.parts
        and path.suffix not in {".pyc", ".pyo", ".nbc", ".nbi"}
    )


def _motorlib_root() -> Path:
    spec = importlib.util.find_spec("motorlib")
    if spec is None or spec.origin is None:
        raise RuntimeError("cannot locate imported motorlib package")
    package_dir = Path(spec.origin).resolve().parent
    # Fingerprint the imported checkout rather than assuming a sibling path.
    return package_dir


def source_fingerprints() -> dict[str, str]:
    """Hash relevant runner, simulator, case, motorlib, and runtime-data files."""
    candidates: dict[str, Path] = {}
    motorlib_root = _motorlib_root()
    roots = (("srm_1d", SCRIPT_ROOT / "srm_1d"), ("motorlib", motorlib_root))
    mathlib_root = motorlib_root.parent / "mathlib"
    if mathlib_root.is_dir():
        roots += (("mathlib", mathlib_root),)
    for label, root in roots:
        for path in _package_files(root):
            candidates[f"{label}/{path.relative_to(root).as_posix()}"] = path

    for directory, suffixes in (
        (SCRIPT_ROOT / "scripts", {".py"}),
        (SCRIPT_ROOT / "cases", {".py", ".json"}),
    ):
        for path in sorted(directory.iterdir()):
            if path.is_file() and path.suffix in suffixes:
                candidates[f"{directory.name}/{path.name}"] = path

    manifests = (
        ("srm_1d/pyproject.toml", SCRIPT_ROOT / "pyproject.toml"),
        ("srm_1d/requirements", SCRIPT_ROOT / "requirements"),
        ("openmotor/requirements.txt", motorlib_root.parent / "requirements.txt"),
    )
    for prefix, item in manifests:
        paths = _package_files(item) if item.is_dir() else ([item] if item.is_file() else [])
        for path in paths:
            relative = path.relative_to(item).as_posix() if item.is_dir() else path.name
            candidates[f"manifest/{prefix}/{relative}"] = path

    return {name: _sha256(path) for name, path in sorted(candidates.items())}


def runtime_versions() -> dict[str, str]:
    """Return Python and direct numerical/runtime dependency versions."""
    versions = {"python": sys.version}
    versions["platform"] = sys.platform
    versions["platform_detail"] = platform.platform()
    versions["machine"] = platform.machine()
    for distribution in (
        "numpy", "scipy", "numba", "PyYAML", "matplotlib", "scikit-fmm",
        "scikit-image", "llvmlite",
    ):
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = "not-installed"
    return versions


def effective_options(
    base_options: dict, variants: tuple[tuple[str, dict], ...],
) -> dict[str, dict]:
    """Resolve final keyword options for each named run variant."""
    return {
        name: {**base_options, **overrides}
        for name, overrides in variants
    }


# Private aliases retain the naming convention used by existing runners.
_source_fingerprints = source_fingerprints
_runtime_versions = runtime_versions
