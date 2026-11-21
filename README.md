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
python scripts/run_limiting_case_verification.py  # Full-duration limiting cases
```

```python
from srm_1d import run_from_ric

result, perf, nozzle, geo, prop = run_from_ric(
    "motors/hasegawa_a.ric", pyrogen="bpnv"
)
print(result["summary"]["P_peak"], perf["total_impulse"])
```

## Propellant catalog

The package includes the supplied openMotor `propellants.yaml` catalog
(77 propellant names). The file is kept in its native YAML format and is
installed with the package. The catalog does **not** modify solver defaults
or existing motor definitions.

List entries or inspect one entry:

```bash
python -m srm_1d.propellant_catalog
python -m srm_1d.propellant_catalog --show "RCS - Super Thunder (H13ST Version)"
```

Python usage:

```python
from srm_1d.openmotor_adapter import load_ric
from srm_1d.propellant_catalog import (
    list_propellants, get_propellant,
    missing_transport_fields, with_catalog_propellant,
)

names = list_propellants()
selected = get_propellant(names[0])
print(missing_transport_fields(selected))
motor = with_catalog_propellant(load_ric("motors/example_bates.ric"), names[0])
```

The returned motor dictionary can be passed to the existing
`ric_to_sim_args` adapter **only when all required transport
properties are supplied**. Many library entries have zero sentinel
values for viscosity, gas conductivity, or heat capacity. They are
**not simulation-ready** until those physical properties are provided
from appropriate sources; the importer deliberately does not estimate
or backfill them. This is a reference-data import, not a validation
of the catalog's values or a GUI library synchronization.

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
