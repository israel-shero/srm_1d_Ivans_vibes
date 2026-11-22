#!/usr/bin/env python3
"""Inspect SimulationChannels attributes."""
from srm_1d.openmotor_adapter import run_from_ric
from pathlib import Path
import sys

motor_path = Path('motors/zerox.ric')
result, perf, nozzle, geo, prop = run_from_ric(
    str(motor_path), t_max=0.1, roughness=32e-6, kappa=0.44, pyrogen='bpnv',
    pyrogen_mass=None, T_ignition=756.0, P_cutoff=0.01e6,
    snapshot_interval=0.2, print_interval=0.0,
)
print('result type:', type(result).__name__)
attrs = [a for a in dir(result) if not a.startswith('_')]
print('all attrs:', attrs)
# Check commonly expected ones
for name in ['t', 'time', 'times', 'P', 'pressure', 'P_head', 'head_pressure', 'channels']:
    if hasattr(result, name):
        val = getattr(result, name)
        print(f'  {name}: type={type(val).__name__}', end='')
        if hasattr(val, '__len__'):
            print(f' len={len(val)}', end='')
        print()
