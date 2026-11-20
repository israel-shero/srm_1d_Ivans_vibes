import json

import numpy as np
import pytest

import scripts.run_limiting_case_verification as limiting
from scripts.run_limiting_case_verification import (
    _closure_limiting_cases,
    _relative_change,
    _steady_interval_metrics,
    run_limiting_cases,
)
from scripts.verification_identity import (
    effective_options,
    runtime_versions,
    source_fingerprints,
)
from srm_1d.openmotor_adapter import convert_propellant, load_ric


def test_steady_interval_metrics_identifies_stationary_pressure():
    result = {
        "time": np.linspace(0.0, 2.0, 21),
        "P_head": np.full(21, 8.0e6),
    }

    metrics = _steady_interval_metrics(result, 0.6, 1.8)

    assert metrics["status"] == "evaluated"
    assert metrics["mean_pressure_mpa"] == pytest.approx(8.0)
    assert metrics["pressure_cv"] == 0.0
    assert metrics["pressure_slope_percent_per_s"] == pytest.approx(0.0)
    assert metrics["screen"]["passed"]


def test_steady_interval_metrics_rejects_large_pressure_slope():
    time = np.linspace(0.0, 2.0, 21)
    result = {
        "time": time,
        "P_head": 8.0e6 + 1.0e6 * time,
    }

    metrics = _steady_interval_metrics(result, 0.6, 1.8)

    assert metrics["pressure_slope_percent_per_s"] > 5.0
    assert not metrics["screen"]["passed"]


def test_steady_interval_requires_two_samples():
    result = {
        "time": np.array([0.0, 1.0]),
        "P_head": np.array([1.0e5, 1.0e5]),
    }

    metrics = _steady_interval_metrics(result, 0.2, 0.8)

    assert metrics["status"] == "not_evaluated"


def test_relative_change_uses_reference_magnitude():
    assert _relative_change(10.0, 8.0) == pytest.approx(-20.0)
    assert _relative_change(0.0, 1.0) is None


def test_closure_limiting_cases_use_converted_propellant():
    motor = load_ric("motors/machbusterNew.ric")
    propellant = convert_propellant(motor["propellant"])

    checks = _closure_limiting_cases(propellant)

    assert checks["zero_crossflow"]["normal_rate_preserved"]
    assert checks["finite_crossflow"]["normal_rate_preserved"]
    assert checks["finite_crossflow"]["erosive_rate_m_s"] > 0.0
    mukunda_paul = checks["finite_crossflow"]["mukunda_paul"]
    assert mukunda_paul["dimensionless_g"] > 35.0
    assert mukunda_paul["eta"] > 1.0
    assert mukunda_paul["erosive_rate_m_s"] > 0.0


def test_resume_rejects_changed_motor_checksum(tmp_path):
    manifest = tmp_path / "limiting_cases.json"
    manifest.write_text(json.dumps({
        "provenance": {
            "motor_sha256": "wrong",
            "baseline_config_sha256": "wrong",
        },
        "runs": {},
    }))

    with pytest.raises(ValueError, match="motor checksum"):
        run_limiting_cases(resume_path=manifest)


def _complete_resume_manifest(path):
    baseline = json.loads(limiting.CONFIG_PATH.read_text())["cases"]["chunc"]
    options = dict(baseline["run_options"])
    options.update({
        "target_propellant_cells": 100,
        "cfl_target": 0.3,
        "print_interval": 0.0,
        "snapshot_interval": 0.1,
    })
    variants = (
        ("baseline", {}),
        ("erosive_disabled", {"diagnostic_disable_erosive": True}),
        ("mukunda_paul", {"erosive_closure": "mukunda_paul"}),
    )
    effective = effective_options(options, variants)
    source = source_fingerprints()
    source_digest = limiting.hashlib.sha256(
        json.dumps(source, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    run_metrics = {
        "peak_pressure_mpa": 1.0,
        "peak_time_s": 1.0,
        "total_impulse_ns": 1.0,
        "performance_burn_time_s": 1.0,
    }
    identity = {
        "source_fingerprint_sha256": source_digest,
        "runtime_versions": runtime_versions(),
    }
    path.write_text(json.dumps({
        "schema_version": 2,
        "configuration": options,
        "provenance": {
            "motor_sha256": limiting._sha256(limiting.PROJECT_ROOT / baseline["motor_path"]),
            "baseline_config_sha256": limiting._sha256(limiting.CONFIG_PATH),
            "source_sha256": source,
            "runtime_versions": runtime_versions(),
        },
        "effective_options": effective,
        "runs": {
            name: {
                "overrides": overrides,
                "identity": {**identity, "effective_options": effective[name]},
                "metrics": dict(run_metrics),
            }
            for name, overrides in variants
        },
    }))


def test_resume_rejects_changed_source_fingerprint(tmp_path, monkeypatch):
    manifest = tmp_path / "limiting_cases.json"
    _complete_resume_manifest(manifest)
    changed = dict(source_fingerprints())
    changed["srm_1d/solver.py"] = "changed"
    monkeypatch.setattr(limiting, "_source_fingerprints", lambda: changed)

    original = manifest.read_bytes()
    with pytest.raises(ValueError, match="source checksum identity"):
        run_limiting_cases(resume_path=manifest)
    assert manifest.read_bytes() == original


def test_resume_rejects_changed_baseline_configuration(tmp_path, monkeypatch):
    manifest = tmp_path / "limiting_cases.json"
    _complete_resume_manifest(manifest)
    changed_config = tmp_path / "baseline_configs.json"
    changed_config.write_text(limiting.CONFIG_PATH.read_text() + "\n")
    monkeypatch.setattr(limiting, "CONFIG_PATH", changed_config)

    original = manifest.read_bytes()
    with pytest.raises(ValueError, match="baseline-config checksum"):
        run_limiting_cases(resume_path=manifest)
    assert manifest.read_bytes() == original


def test_resume_rejects_changed_effective_options_and_run_overrides(tmp_path):
    manifest = tmp_path / "limiting_cases.json"
    _complete_resume_manifest(manifest)
    payload = json.loads(manifest.read_text())
    payload["runs"]["mukunda_paul"]["overrides"] = {"erosive_closure": "ma"}
    manifest.write_text(json.dumps(payload))

    original = manifest.read_bytes()
    with pytest.raises(ValueError, match="overrides do not match"):
        run_limiting_cases(resume_path=manifest)
    assert manifest.read_bytes() == original


@pytest.mark.parametrize(
    ("section", "match"),
    (("effective_options", "effective options"), ("configuration", "configuration"),
     ("runtime_versions", "runtime versions")),
)
def test_resume_rejects_changed_saved_settings_without_manifest_mutation(
    tmp_path, section, match,
):
    manifest = tmp_path / "limiting_cases.json"
    _complete_resume_manifest(manifest)
    payload = json.loads(manifest.read_text())
    if section == "effective_options":
        payload[section]["baseline"]["cfl_target"] = 0.01
    elif section == "configuration":
        payload[section]["cfl_target"] = 0.01
    else:
        payload["provenance"][section]["numpy"] = "0.0"
    manifest.write_text(json.dumps(payload))
    original = manifest.read_bytes()

    with pytest.raises(ValueError, match=match):
        run_limiting_cases(resume_path=manifest)
    assert manifest.read_bytes() == original


def test_unchanged_resume_skips_completed_simulations(tmp_path, monkeypatch):
    manifest = tmp_path / "limiting_cases.json"
    _complete_resume_manifest(manifest)
    monkeypatch.setattr(
        limiting, "run_from_ric",
        lambda *args, **kwargs: pytest.fail("completed runs should be skipped"),
    )

    output = run_limiting_cases(resume_path=manifest)

    assert output == tmp_path
    assert json.loads(manifest.read_text())["status"] == "complete"
