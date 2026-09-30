# Validated Python Environment

`py310-macos-x86_64.lock.txt` is the tested dependency lock for the current
Intel-macOS development machine. It resolves the two installation failures seen
with unconstrained dependencies:

- current Numba/llvmlite releases attempted unsupported local builds;
- `scikit-image 0.25.2` required SciPy >= 1.11.4 while the openMotor numerical
  stack uses SciPy 1.10.1.

The resolution keeps the established numerical core and selects
`scikit-image 0.21.0`, which provides a CPython 3.10 Intel-macOS wheel and is
compatible with SciPy 1.10.1.

## Create the environment

From the `srm_1d` repository root:

```bash
/usr/local/opt/python@3.10/bin/python3.10 -m venv ../.venv-srm
../.venv-srm/bin/python -m pip install \
  pip==26.0.1 setuptools==70.3.0 wheel==0.46.3
../.venv-srm/bin/python -m pip install \
  -c requirements/py310-macos-x86_64.lock.txt -e ".[fmm,dev]"
../.venv-srm/bin/python -m pip check
```

The compatible GUI/plugin checkout is:

```text
https://github.com/eJar0k/openMotor.git
branch: staging
validated commit: 92e425f7832ed3f3e926ff985117649e3c77881a
```

Install that checkout into the same environment with its requirements
constrained by this lock, then build its extension and generated UI modules.
The openMotor setup currently writes its normal log under the user Library, so
the commands may need normal host permissions rather than a restricted sandbox.

```bash
../.venv-srm/bin/python -m pip install \
  -c requirements/py310-macos-x86_64.lock.txt \
  Cython pyqt-distutils PyQt6 ezdxf platformdirs decorator docopt Sphinx
../.venv-srm/bin/python -m pip install \
  --no-build-isolation --no-deps -e ../openMotor-compatible

cd ../openMotor-compatible
PATH="../.venv-srm/bin:/usr/bin:/bin" MPLCONFIGDIR=/private/tmp \
  ../.venv-srm/bin/python setup.py build_ui
cd ../srm_1d
```

Set `SRM1D_OPENMOTOR_PATH` for tests and examples:

```bash
export SRM1D_OPENMOTOR_PATH="$(cd ../openMotor-compatible && pwd)"
MPLCONFIGDIR=/private/tmp ../.venv-srm/bin/python -m pytest tests -q
MPLCONFIGDIR=/private/tmp ../.venv-srm/bin/python \
  -m examples.hasegawa_motor_a
```

## Validation record

Validated on 2026-09-29 with CPython 3.10.10:

- `pip check`: no broken requirements;
- imports: all numerical, FMM, Qt, `srm_1d`, `motorlib.solvers`, and
  `motorlib.taper` imports passed;
- test suite: 408 passed, 95 warnings in 755.99 seconds;
- Hasegawa A: all reported physics/performance metrics matched the established
  baseline, including peak pressure 6.14 MPa, first burnout 3.678 s, and total
  impulse 22,657.8 N s.

The generated PNG hashes differ from the earlier unlocked environment because
the rendering dependencies are now pinned. This is not recorded as numerical
identity of plot files; the solver metrics matched at their reported precision.

This lock is platform-specific. Other Python versions, Apple Silicon, Windows,
and Linux require separately tested locks rather than silently reusing this one.

## Verify repository pins

`workspace-lock.json` pins the compatible source pair independently of the
Python package lock. From the `srm_1d` repository root, run:

```bash
../.venv-srm/bin/python scripts/verify_workspace.py
```

The verifier checks Python 3.10, both origin URLs, the exact openMotor commit,
that the current `srm_1d` branch contains the validated baseline commit, the
fork-only solver/taper APIs, and `pip check`. The ancestry rule avoids a
self-referential lock that would fail as soon as this repository advances.
Local modifications are warnings so normal development remains usable. Release
or archival checks can require pristine trees:

```bash
../.venv-srm/bin/python scripts/verify_workspace.py --strict-clean
```

The script never fetches, checks out, resets, or modifies either repository.
