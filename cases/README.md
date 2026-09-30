# Validation Case Registry

`registry.json` is the machine-readable source/data inventory for baseline and
validation work. It records case readiness, exact input checksums, measurement
schemas, and limitations. A `candidate` label means the inputs can be prepared
for comparison; it does not mean a repository research result was reproduced.

Validate metadata and source files from the repository root:

```bash
../.venv-srm/bin/python scripts/validate_case_registry.py
```

Load a processed trace without repeating per-file column assumptions:

```python
from cases import load_measurement

trace = load_measurement("chunc")
print(trace["time"], trace["pressure"])  # seconds, MPa
```

The loader refuses raw, missing, and absent measurements. In particular, it
will not silently tare BALLSStick, invent the missing ST2.0 trace, or treat a
motor with no paired data as validated.

## Recorded baseline runner

`baseline_configs.json` freezes the arguments from the documented Hasegawa A
and Chunc examples without changing solver defaults. Prepare a provenance
manifest without running the solver:

```bash
MPLCONFIGDIR=/private/tmp \
SRM1D_OPENMOTOR_PATH=../openMotor-compatible \
../.venv-srm/bin/python -m cases.baseline_runner \
  --case hasegawa_a --dry-run
```

Run either registered baseline and save metrics, a sampled CSV, and plots:

```bash
MPLCONFIGDIR=/private/tmp \
SRM1D_OPENMOTOR_PATH=../openMotor-compatible \
../.venv-srm/bin/python -m cases.baseline_runner \
  --case hasegawa_a --plots
```

Each run writes `run_manifest.json` before solver execution and records a
failure there if execution aborts. A completed manifest captures repository
state, dependency versions, input hashes, requested and resolved options,
runtime motor properties, and output hashes. Its reproduction claim remains
false; only metrics created in that run directory are established by the run.
