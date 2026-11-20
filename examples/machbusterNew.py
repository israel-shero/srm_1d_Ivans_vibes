"""
machbusterNew.py -- run the machbusterNew .ric through the v0.7.0 workflow.

This mirrors `BALLSstick.py` and is intentionally simple: load the .ric,
auto-discover the sibling transport YAML, use the built-in BPNV pyrogen,
and write plots under `artifacts/machbusterNew/`.

Usage:
    python -m examples.machbusterNew_test
"""

from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.plotting import plot_pressure, plot_flow_snapshot, plot_summary
from srm_1d.run_artifacts import artifact_dir

# Legacy plotting overlay only: its (0, 0) point is not in the primary
# Chunc CSV/workbook and has no independently documented event-time meaning.
CHUNC_EXPERIMENTAL = {
    'label': 'Experimental (Chunc / machbusterNew)',
    'time': np.array([
        0.0, 0.01, 0.044, 0.085, 0.126, 0.167, 0.207, 0.248, 0.289, 0.33,
        0.37, 0.411, 0.452, 0.493, 0.533, 0.574, 0.615, 0.655, 0.696,
        0.737, 0.778, 0.818, 0.859, 0.9, 0.941, 0.981, 1.022, 1.063,
        1.104, 1.144, 1.185, 1.226, 1.267, 1.307, 1.348, 1.389, 1.429,
        1.47, 1.511, 1.552, 1.592, 1.633, 1.674, 1.715, 1.755, 1.796,
        1.837, 1.878, 1.918, 1.959, 2, 2.041, 2.081, 2.122, 2.163, 2.203, 2.244, 2.285, 2.326

    ]),
    'pressure': np.array([
        0.0, 8.466, 8.735, 8.807, 8.83, 8.863, 8.815, 8.808, 8.801, 8.834,
        8.834, 8.883, 8.841, 8.808, 8.789, 8.76, 8.709, 8.67, 8.617, 8.57,
        8.61, 8.602, 8.592, 8.52, 8.412, 8.106, 7.254, 6.251, 5.644, 5.241,
        4.852, 4.273, 3.859, 3.498, 3.112, 2.699, 2.44, 2.199, 1.955, 1.69,
        1.492, 1.308, 1.149, 0.979, 0.802, 0.669, 0.553, 0.45, 0.355, 0.288,
        0.237, 0.193, 0.15, 0.125, 0.103, 0.081, 0.062, 0.049, 0.041


    ]),
}


CASE_NAME = 'machbusterNew'
MOTOR_PATH = Path(__file__).resolve().parents[1] / 'motors' / 'machbusterNew.ric'


def main():
    # Per-run artifact dir matches the run_template.py convention.
    OUTPUT_DIR = artifact_dir(CASE_NAME)

    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        # v0.7.5 cross-motor re-LHS optimum (shared fired-motor fit;
        # docs/v0_7_5/RESULT.md). Chunc stays the worst fit (high-L/D
        # QS-erosive limitation, IGNITION_SPIKE_CLOSEOUT.md).
        roughness=32e-6,
        kappa=0.44,
        pyrogen='bpnv',
        pyrogen_mass=None,
        T_ignition=756.0,
        P_cutoff=0.01e6,
        snapshot_interval=0.2,
        print_interval=0.2,
        verbose=True,
    )

    summary = result['summary']
    print(
        f"machbusterNew: P_peak={summary['P_peak'] / 1e6:.2f} MPa, "
        f"t_burn={summary['t_burn']:.3f} s, "
        f"impulse={perf['total_impulse']:.1f} N*s, "
        f"designation={perf['motor_designation']}"
    )

    plot_pressure(
        result,
        title="machbusterNew Pressure Trace",
        experimental=CHUNC_EXPERIMENTAL,
        save_path=str(OUTPUT_DIR / "pressure.png"),
    )

    plot_flow_snapshot(
        result,
        t_target=0.01,
        save_path=str(OUTPUT_DIR / "flow.png"),
    )

    plot_summary(
        result,
        performance=perf,
        title="machbusterNew Simulation Summary",
        experimental=CHUNC_EXPERIMENTAL,
        save_path=str(OUTPUT_DIR / "summary.png"),
    )

    plt.close('all')
    print(f"Plots saved under {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
