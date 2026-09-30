import importlib
import json

import pytest

from cases.baseline_runner import DEFAULT_CONFIGS, run_baseline


@pytest.fixture(scope="module")
def baseline_configs():
    return json.loads(DEFAULT_CONFIGS.read_text(encoding="utf-8"))["cases"]


def test_hasegawa_config_matches_documented_example(baseline_configs):
    example = importlib.import_module("examples.hasegawa_motor_a")
    config = baseline_configs["hasegawa_a"]

    assert config["motor_path"] == "motors/hasegawa_a.ric"
    assert config["run_options"] == {
        "roughness": 32e-6,
        "kappa": 0.44,
        "pyrogen": "bpnv",
        "pyrogen_mass": None,
        "T_ignition": 756.0,
        "P_cutoff": 0.05e6,
        "snapshot_interval": 0.2,
        "print_interval": 0.2,
    }
    assert config["comparison"]["time_offset_s"] == pytest.approx(
        example.EXPERIMENTAL_TIME_OFFSET
    )


def test_chunc_config_matches_documented_example(baseline_configs):
    example = importlib.import_module("examples.run_chunc")
    config = baseline_configs["chunc"]

    assert config["motor_path"] == "motors/machbusterNew.ric"
    assert config["run_options"] == {
        "pyrogen": example.PYROGEN,
        "pyrogen_mass": example.PYROGEN_MASS,
        "pyrogen_throat_area": example.PYROGEN_THROAT_AREA,
        "pyrogen_volume": example.PYROGEN_VOLUME,
        "roughness": example.ROUGHNESS,
        "kappa": example.KAPPA,
        "T_ignition": example.T_IGNITION,
        "k_solid": example.K_SOLID,
        "P_cutoff": example.P_CUTOFF,
        "t_max": example.TIME_MAX,
        "cfl_target": example.CFL_TARGET,
        "snapshot_interval": example.SNAPSHOT_INTERVAL,
        "print_interval": example.PRINT_INTERVAL,
    }
    assert config["plots"]["flow_snapshot_time_s"] == pytest.approx(
        example.FLOW_SNAPSHOT_TIME
    )


@pytest.mark.parametrize("case_id", ["hasegawa_a", "chunc"])
def test_dry_run_writes_provenance_without_reproduction_claim(
    case_id, tmp_path, monkeypatch
):
    monkeypatch.setenv(
        "SRM1D_OPENMOTOR_PATH",
        str(tmp_path / "intentionally-absent-openmotor"),
    )
    output = run_baseline(case_id, output_root=tmp_path, dry_run=True)
    manifest = json.loads((output / "run_manifest.json").read_text())

    assert manifest["status"] == "dry_run_complete"
    assert manifest["case_id"] == case_id
    assert manifest["claims"]["historical_research_results_reproduced"] is False
    assert manifest["repositories"]["srm_1d"]["commit"]
    assert manifest["repositories"]["openMotor"]["available"] is False
    assert manifest["active_inputs"]["transport"]["mode"] == "embedded_ric"
    assert manifest["configuration"]["resolved_simulation_options"][
        "P_ambient"
    ] > 0
