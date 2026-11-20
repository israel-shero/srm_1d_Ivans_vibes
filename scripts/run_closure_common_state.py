#!/usr/bin/env python3
"""Common-state closure comparison: Ma vs Mukunda-Paul on identical local states.

Scope (plan T06, Section 8):
  * Run a Ma baseline simulation on Chunc / machbusterNew and save flow snapshots.
  * Replay BOTH closures on the identical saved local states (P, Re, D_hyd, x,
    is_burning) so only the constitutive burn-rate response differs — the flow
    that drives each closure is held fixed. This isolates driver effects from
    closure response effects.
  * Report normal rate, total rate, relative augmentation (eta), and the
    Ma/MP divergence across the reached regime (startup spike through plateau).

This is a *common-state* (constitutive) comparison only. It does NOT run fully
coupled simulations with each closure; that is a separate step and would let
each closure alter the flow that drives it.

Outputs a JSON manifest under artifacts/verification_chunc_common_state/.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from scripts.verification_identity import (
    effective_options as _effective_options,
    runtime_versions as _runtime_versions,
    source_fingerprints as _source_fingerprints,
)
from srm_1d.burn_rate import (
    burn_rate_cell,
    mukunda_paul_enhancement,
)
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.run_artifacts import artifact_dir, verify_run_health

MOTOR_PATH = SCRIPT_ROOT / "motors" / "machbusterNew.ric"
ROUGHNESS = 32e-6
KAPPA = 0.44
T_MAX_S = 1.3
SNAPSHOT_INTERVAL_S = 0.02


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_value(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args], cwd=SCRIPT_ROOT, check=False,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else ""


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
def main() -> int:
    source_sha256 = _sha256(SCRIPT_ROOT / "srm_1d" / "burn_rate.py")
    source_identity_sha256 = _sha256(MOTOR_PATH)
    runtime_versions = _runtime_versions()

    # --- 1. Ma baseline simulation, snapshots saved -----------------------
    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        roughness=ROUGHNESS,
        kappa=KAPPA,
        pyrogen="bpnv",
        pyrogen_mass=None,
        T_ignition=756.0,
        P_cutoff=0.01e6,
        t_max=T_MAX_S,
        snapshot_interval=SNAPSHOT_INTERVAL_S,
        print_interval=0.0,
        verbose=False,
    )
    verify_run_health(result, motor_name="chunc-common-state", min_t_burn_s=0.4,
                      raise_on_fail=True)
    snaps = result["snapshots"]

    # --- 2. Material / correlation properties (shared by both closures) ---
    tab = prop.representative_tab()
    tab_min_p, tab_max_p, tab_a, tab_n = prop.tab_arrays()
    n_tabs = len(prop.tabs)
    mu_gas = float(prop.mu_gas)
    k_thermal = float(prop.k_gas)
    cp_gas = float(prop.Cp_gas)
    pr = cp_gas * mu_gas / k_thermal
    t_flame = float(tab.T_flame)
    t_surface = float(prop.T_surface)
    rho_p = float(prop.rho_propellant)
    cps = float(prop.Cps)
    t_initial = float(prop.T_initial)

    # --- 3. Replay both closures on identical saved local states ----------
    runs = []
    for snap in snaps:
        t = float(snap["t"])
        g = snap["is_grain"] & snap["is_burning"]
        if not np.any(g):
            continue
        P = np.asarray(snap["P"][g], dtype=float)
        rho = np.asarray(snap["rho"][g], dtype=float)
        u = np.asarray(snap["u"][g], dtype=float)
        d_hyd = np.asarray(snap["D_port"][g], dtype=float)
        x = np.asarray(snap["x"][g], dtype=float)

        # Local Reynolds number from the saved flow state (cell-centered |u|).
        re = rho * np.abs(u) * d_hyd / mu_gas

        N = len(P)
        r0 = np.zeros(N)
        r_ma_total = np.zeros(N)
        r_ma_erosive = np.zeros(N)
        r_mp_total = np.zeros(N)
        r_mp_erosive = np.zeros(N)
        eta_mp = np.zeros(N)
        g_val = np.zeros(N)
        for i in range(N):
            rt_m, re_m = burn_rate_cell(
                P[i], re[i], d_hyd[i], x[i], ROUGHNESS,
                pr, k_thermal, cp_gas, t_flame, t_surface,
                rho_p, cps, t_initial,
                tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
                KAPPA, mu_gas, 0,  # erosive_closure = Ma
            )
            r_ma_total[i], r_ma_erosive[i] = rt_m, re_m
            # Saint-Robert normal rate is total minus Ma erosive increment.
            r0[i] = r_ma_total[i] - r_ma_erosive[i]

            # Mukunda-Paul constitutive response on the identical state.
            eta, gv = mukunda_paul_enhancement(
                re[i] * mu_gas / d_hyd[i],  # mass flux G = Re*mu/D
                r0[i], d_hyd[i], rho_p, mu_gas,
            )
            r_mp_total[i] = r0[i] * eta
            r_mp_erosive[i] = max(r_mp_total[i] - r0[i], 0.0)
            eta_mp[i] = eta
            g_val[i] = gv

        eta_ma = np.where(r0 > 0, r_ma_total / np.maximum(r0, 1e-12), 1.0)
        eta_mp_safe = np.where(r0 > 0, r_mp_total / np.maximum(r0, 1e-12), 1.0)

        ipk = int(np.argmax(r_ma_erosive))
        runs.append({
            "t_s": t,
            "peak_pressure_pa": float(np.max(P)),
            "n_burning_cells": int(np.count_nonzero(g)),
            "re_range": [float(re.min()), float(re.max())],
            "pressure_range_mpa": [float(P.min() / 1e6), float(P.max() / 1e6)],
            "peak_erosive_cell": {
                "re": float(re[ipk]),
                "g": float(g_val[ipk]),
                "r0_mm_s": float(r0[ipk] * 1e3),
                "eta_ma": float(eta_ma[ipk]),
                "eta_mp": float(eta_mp_safe[ipk]),
                "ma_over_mp": float(eta_ma[ipk] / max(eta_mp_safe[ipk], 1e-12)),
                "r_erosive_ma_mm_s": float(r_ma_erosive[ipk] * 1e3),
                "r_erosive_mp_mm_s": float(r_mp_erosive[ipk] * 1e3),
            },
            "domain_max": {
                "eta_ma": float(eta_ma.max()),
                "eta_mp": float(eta_mp_safe.max()),
                "ma_over_mp": float(
                    eta_ma.max() / max(eta_mp_safe.max(), 1e-12)),
            },
            "cells_g_above_threshold": int(np.count_nonzero(g_val > 35.0)),
            "frac_cells_ma_over_mp_1pct": float(
                np.count_nonzero(eta_ma > 1.01 * eta_mp_safe) / N),
        })

    # --- 4. Regime aggregation -------------------------------------------
    times = np.array([r["t_s"] for r in runs])
    # Startup/spike: before plateau (~<0.5 s); plateau: >=0.5 s
    startup = [r for r in runs if r["t_s"] < 0.5]
    plateau = [r for r in runs if r["t_s"] >= 0.5]

    def agg(group):
        if not group:
            return {}
        return {
            "n_snapshots": len(group),
            "t_range_s": [group[0]["t_s"], group[-1]["t_s"]],
            "peak_eta_ma": max(r["domain_max"]["eta_ma"] for r in group),
            "peak_eta_mp": max(r["domain_max"]["eta_mp"] for r in group),
            "peak_ma_over_mp": max(r["domain_max"]["ma_over_mp"] for r in group),
            "peak_re": max(r["re_range"][1] for r in group),
            "mean_frac_cells_ma_over_mp_1pct": float(
                np.mean([r["frac_cells_ma_over_mp_1pct"] for r in group])),
        }

    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "study": "chunc_common_state_closure_comparison",
        "scope": {
            "note": "Common-state (constitutive) comparison only; no coupled runs.",
            "motor": str(MOTOR_PATH.relative_to(SCRIPT_ROOT)),
            "motor_sha256": source_identity_sha256,
            "burn_rate_source_sha256": source_sha256,
            "closure_baseline": "ma",
            "closures_compared": ["ma", "mukunda_paul"],
            "t_max_s": T_MAX_S,
            "snapshot_interval_s": SNAPSHOT_INTERVAL_S,
            "n_snapshots_replayed": len(runs),
        },
        "identity": {
            "source_fingerprint_sha256": _source_fingerprints(),
            "runtime_versions": runtime_versions,
            "git_commit": _git_value("rev-parse", "HEAD"),
            "git_status": _git_value("status", "--short"),
        },
        "regime": {
            "startup_spike_t_lt_0p5_s": agg(startup),
            "plateau_t_ge_0p5_s": agg(plateau),
        },
        "snapshots": runs,
    }
    out = artifact_dir("verification_chunc_common_state")
    _write_json(out / "comparison.json", manifest)
    print(f"Common-state closure comparison: {out / 'comparison.json'}")
    print(f"  snapshots replayed: {len(runs)}")
    print(f"  startup spike:   peak eta_Ma={manifest['regime']['startup_spike_t_lt_0p5_s']['peak_eta_ma']:.2f} "
          f"peak eta_MP={manifest['regime']['startup_spike_t_lt_0p5_s']['peak_eta_mp']:.2f} "
          f"max Ma/MP={manifest['regime']['startup_spike_t_lt_0p5_s']['peak_ma_over_mp']:.2f}")
    print(f"  plateau:         peak eta_Ma={manifest['regime']['plateau_t_ge_0p5_s']['peak_eta_ma']:.2f} "
          f"peak eta_MP={manifest['regime']['plateau_t_ge_0p5_s']['peak_eta_mp']:.2f} "
          f"max Ma/MP={manifest['regime']['plateau_t_ge_0p5_s']['peak_ma_over_mp']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

