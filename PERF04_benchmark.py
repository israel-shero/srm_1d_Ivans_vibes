"""PERF-04: Batch throughput, parallelism and resource budgets benchmark."""
from __future__ import annotations

import time
from pathlib import Path
from typing import Dict

import numpy as np

# Absolute paths for this worktree
_REPO_ROOT = Path("/Users/israelshero/Development/MotorDevelopment/srm_1d")
_MOTOR_PATH = _REPO_ROOT / "motors" / "hasegawa_a.ric"

sys.path.insert(0, str(_REPO_ROOT))

from srm_1d.openmotor_adapter import run_from_ric


def run_single() -> dict:
    result, perf = run_from_ric(
        str(_MOTOR_PATH),
        gas_props={"mu": 8.842e-5, "k": 0.3685, "Cp": 2060.0},
        kappa=0.45, t_max=10.0, P_cutoff=0.05e6, cfl_target=0.3,
    )
    return {"result": result, "perf": perf}


def bench_serial(n: int = 5) -> Dict[str, float]:
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        run_single()
        times.append(time.perf_counter() - t0)
    return {
        "mean_ms": float(np.mean(times) * 1000),
        "median_ms": float(np.median(times) * 1000),
        "std_ms": float(np.std(times) * 1000),
        "min_ms": float(np.min(times) * 1000),
        "max_ms": float(np.max(times) * 1000),
    }


def bench_pool(executor, max_workers: int, n: int = 3) -> Dict[str, float]:
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        futures = [executor.submit(run_single) for _ in range(max_workers)]
        for f in futures:
            f.result()
        times.append(time.perf_counter() - t0)
    return {
        "mean_ms": float(np.mean(times) * 1000),
        "median_ms": float(np.median(times) * 1000),
        "std_ms": float(np.std(times) * 1000),
    }


def main():
    print("PERF-04 Benchmark Suite")
    print("=" * 50)
    s = bench_serial(5)
    print(f"Serial: mean={s['mean_ms']:.1f}ms median={s['median_ms']:.1f}ms")
    from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
    for w in [1, 2, 4]:
        p = bench_pool(ProcessPoolExecutor(max_workers=w), w, 3)
        print(f"Process w={w}: mean={p['mean_ms']:.1f}ms median={p['median_ms']:.1f}ms")
    for w in [1, 2, 4, 8]:
        t = bench_pool(ThreadPoolExecutor(max_workers=w), w, 3)
        print(f"Thread w={w}: mean={t['mean_ms']:.1f}ms median={s['median_ms']:.1f}ms")
    print("\nDeterminism check:")
    outs = [run_single() for _ in range(3)]
    peaks = [r["result"]["summary"]["P_peak"] for r in outs]
    burns = [float(r["result"]["summary"]["t_burn"]) for r in outs]
    print(f"  P_peak: {peaks} consistent={len(set(round(p, 6) for p in peaks)) == 1}")
    print(f"  t_burn: {burns} consistent={len(set(round(t, 6) for t in burns)) == 1}")


if __name__ == "__main__":
    main()