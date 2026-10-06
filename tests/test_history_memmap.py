import numpy as np
import pytest

from srm_1d import run_simulation
from srm_1d.igniter_plenum import PyrogenChamber
from srm_1d.nozzle import Nozzle
from srm_1d.propellant import Pyrogen
from tests._motor_fixtures import hasegawa_propellant_1, single_cylinder_geo


def _inputs():
    geo = single_cylinder_geo(
        D_bore=0.030, D_outer=0.060, length=0.120,
        target_propellant_cells=12,
    )
    propellant = hasegawa_propellant_1()
    nozzle = Nozzle(D_throat=0.010, D_exit=0.020, efficiency=0.95)
    pyrogen = Pyrogen(
        name="history-memmap-test", a=3.0e-5, n=0.5, rho=1700.0,
        T_flame=2800.0, M=0.030, gamma=1.25, impetus_W=5000.0,
        heat_flux_cal_cm2_s=69.4,
    )
    chamber = PyrogenChamber(
        pyrogen=pyrogen, m_pyrogen_initial=0.003, A_burn_initial=5.0e-4,
        A_throat=2.0e-5, V_plenum=3.0e-6, burn_law="end_burning",
    )
    return geo, propellant, nozzle, chamber


def _run(path=None):
    return run_simulation(
        *_inputs(), T_ignition=294.0, t_max=0.001, P_cutoff=1.0,
        history_capacity=64, snapshot_interval=0.001, dt_max=1.0e-5,
        burn_update_interval=1, verbose=False, history_memmap_path=path,
    )


def _memmap_base(array):
    base = array
    while base is not None:
        if isinstance(base, np.memmap):
            return base
        base = base.base
    return None


def test_disk_backed_scalar_histories_match_memory_results(tmp_path):
    expected = _run()
    path = tmp_path / "scalar-history.bin"
    actual = _run(path)

    scalar_history_keys = [
        key for key, value in expected.items()
        if isinstance(value, np.ndarray)
        and value.ndim == 1
        and value.size == expected["summary"]["steps"]
    ]
    assert len(scalar_history_keys) == 44
    for key in scalar_history_keys:
        array = actual[key]
        assert isinstance(array, np.ndarray)
        assert _memmap_base(array) is not None, key
        np.testing.assert_array_equal(array, expected[key])
    assert isinstance(actual["snapshots"][0]["P"], np.ndarray)
    assert _memmap_base(actual["snapshots"][0]["P"]) is None

    memory_summary = expected["summary"].copy()
    memmap_summary = actual["summary"].copy()
    for field in (
        "wall_time", "history_storage", "history_memmap_path",
    ):
        memory_summary.pop(field)
        memmap_summary.pop(field)
    np.testing.assert_equal(memmap_summary, memory_summary)
    assert expected["summary"]["history_storage"] == "memory"
    assert expected["summary"]["history_memmap_path"] is None
    summary = actual["summary"]
    expected_bytes = 42 * 64 * np.dtype(np.float64).itemsize
    assert summary["history_storage"] == "memmap"
    assert summary["history_memmap_path"] == str(path.resolve())
    assert summary["history_storage_bytes"] == expected_bytes
    assert path.is_file()
    assert path.stat().st_size == expected_bytes


def test_memmap_refuses_existing_file(tmp_path):
    path = tmp_path / "already-there.bin"
    path.write_bytes(b"owned by caller")

    with pytest.raises(FileExistsError):
        _run(path)

    assert path.read_bytes() == b"owned by caller"


def test_memmap_requires_existing_parent(tmp_path):
    path = tmp_path / "missing" / "history.bin"

    with pytest.raises(FileNotFoundError, match="parent directory does not exist"):
        _run(path)
