# srm_1d

A **1D transient finite-volume solid rocket motor internal-ballistics
simulator** with the Ma et al. (2020) erosive-burning model. A Numba-JIT
compiled PISO time loop (≈45–90k steps/s) resolves the axial pressure,
velocity, temperature, and per-cell grain regression through the full
firing — ignition transient, plateau, erosive lift, and tail-off.

## Features

- **PISO transient core** — pressure–velocity coupling, TDMA, adaptive CFL.
- **Ma 2020 erosive burning** — Haaland → Gnielinski → bisection closure.
- **Pyrogen ignition** — Goodman integral solid-heating subsolver, choked
  igniter plenum, and uncontained (`head_basket`/`aft_basket`) topologies.
- **N-species bore gas** with per-cell transport (frozen / effective).
- **openMotor integration** — loads `.ric` motors and registers as an
  openMotor solver plugin (`srm_1d.srm1d_plugin`).
- **FMM grain regression** via a bridge to a local openMotor checkout.

## Install

For the tested Python 3.10 Intel-macOS environment, use the checked-in
constraints so pip does not select incompatible current releases:

```bash
python3.10 -m venv ../.venv-srm
../.venv-srm/bin/python -m pip install \
  pip==26.0.1 setuptools==70.3.0 wheel==0.46.3
../.venv-srm/bin/python -m pip install \
  -c requirements/py310-macos-x86_64.lock.txt -e ".[fmm,dev]"
../.venv-srm/bin/python -m pip check
```

The generic editable install remains `pip install -e ".[fmm,dev]"`, but it is
not a reproducible environment. See `requirements/README.md` for the validated
versions, compatible openMotor setup, and repository-pin verifier.

Requires Python >= 3.10 (developed and currently locked on 3.10). Core deps:
numpy, scipy, numba, pyyaml, matplotlib.

## Quick start

Run from the repo root:

```bash
python -m examples.hasegawa_motor_a   # runs motors/hasegawa_a.ric
python -m pytest tests/               # test suite
python -m cases.baseline_runner --case hasegawa_a --dry-run
python scripts/run_numerical_verification.py  # Chunc startup grid/CFL study
python scripts/run_numerical_verification.py --profile full
python scripts/run_numerical_verification.py --profile full \
  --single-cells 100 --single-cfl 0.075 --history-capacity 6500000
```

```python
from srm_1d import run_from_ric

result, perf, nozzle, geo, prop = run_from_ric(
    "motors/hasegawa_a.ric", pyrogen="bpnv"
)
print(result["summary"]["P_peak"], perf["total_impulse"])
```

## Layout

| Path | Role |
|------|------|
| `srm_1d/` | the importable package (ships): solver core, `tools/`, `pyrogens/` |
| `motors/` | motor data — `<motor>.ric` (transport embedded per-tab) |
| `examples/` | runnable studies (`python -m examples.<name>`) |
| `cases/` | checksummed validation-case registry and safe measurement loaders |
| `tests/` | pytest suite |
| `docs/` | design packages and development narrative |

See [`srm_1d/README.md`](srm_1d/README.md) for the full public API and
validated parameters, and [`CLAUDE.md`](CLAUDE.md) for an orientation map.

## License

MIT.
