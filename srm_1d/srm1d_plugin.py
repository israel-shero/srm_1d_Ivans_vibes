"""
srm1d_plugin.py — v0.8.0 Phase 5: srm_1d as an openMotor solver plugin (D1).

Registers srm_1d's 1-D PISO transient solver against openMotor's solver
registry (``motorlib.solvers``) so the GUI / simulation manager can select
it alongside the built-in quasi-steady solver. Importing this module
performs the registration.

The solver consumes an openMotor ``Motor`` (via ``Motor.getDict()`` — the
same ``{nozzle, propellant, grains, config}`` shape the ``.ric`` adapter
already converts, now carrying per-tab transport from Phase 3), runs
``run_simulation``, and maps the resulting channels back into an openMotor
``SimulationResult`` so the GUI renders it natively.

Scope: the scalar per-step channels (time, pressure, force, exit pressure,
throat change, Kn) and the per-grain multi-value channels (mass, massFlow,
massFlux, regression, web, machNumber) are populated — the latter aggregated
from srm_1d's per-cell axial snapshots via the cell→segment map and
interpolated onto the per-step time grid (Phase 6 task 3). The igniter is
read from the motor's own ``data.igniter`` block (Phase 6 task 4), falling
back to the built-in BPNV pyrogen for motors that predate the block. The GUI
solver picker (D6) is Phase 6 task 1, verified live in the openMotor app.
"""

import math as _math
import threading
import time as _time

import numpy as np

from .fmm_grain import _setup_openmotor_path

SOLVER_NAME = 'srm_1d-transient'

# The transient solver emits one sample per solver step (10^5+ for a typical
# burn). openMotor's graph widget re-converts (per-point, in Python) and
# re-plots the whole series on every channel switch, so handing it the raw
# resolution causes multi-second hangs. We decimate the per-step channels to
# this many points for the GUI SimulationResult (a display artifact), while
# always preserving the peak-pressure / peak-thrust / endpoint samples so
# peak P, burn time, and impulse stay faithful. srm_1d's own full-resolution
# channels are untouched (analysis / CSV use those).
GUI_MAX_POINTS = 5000


def _om():
    """Lazily import the openMotor modules this plugin builds on."""
    _setup_openmotor_path()
    from motorlib.solvers import SolverPlugin, register_solver  # type: ignore
    from motorlib.simResult import SimulationResult  # type: ignore
    return SolverPlugin, register_solver, SimulationResult


def _transient_config_schema():
    """Build the srm_1d transient solver's run-parameter schema as an
    openMotor ``PropertyCollection``. Property keys match ``run_simulation``
    keyword arguments so the GUI-collected dict can be passed straight through
    as the solver ``config`` (→ ``simulate_motor`` overrides). Defaults mirror
    ``run_simulation``'s. ``ambPressure`` is intentionally absent — it stays in
    the shared global config (the adapter reads it from the motor)."""
    _setup_openmotor_path()
    from motorlib.properties import (  # type: ignore
        PropertyCollection, FloatProperty,
    )
    schema = PropertyCollection()
    schema.props['t_max'] = FloatProperty('Max Simulation Time', 's', 0.01, 600.0)
    schema.props['P_cutoff'] = FloatProperty('Pressure Cutoff', 'Pa', 1.0e3, 5.0e6)
    schema.props['cfl_target'] = FloatProperty('CFL Target', '', 0.01, 1.0)
    # Roughness is edited in micrometers (its own unit category, default um),
    # converted to metres at the run boundary (see ROUGHNESS_KEY handling).
    schema.props['roughness'] = FloatProperty('Surface Roughness', 'um', 0.0, 1000.0)
    # kappa is the Gnielinski (Ma Eq. 9) temperature-ratio exponent in the
    # turbulent Nusselt correction (T_gas/T_surface)^kappa; recommended ~0.45
    # for heated gas. It is NOT an erosive-burning coefficient.
    schema.props['kappa'] = FloatProperty('Gnielinski Temp-Ratio Exponent', '', 0.0, 2.0)
    schema.props['T_ignition'] = FloatProperty('Ignition Temperature', 'K', 300.0, 2000.0)
    # Seed defaults (FloatProperty clamps to range; values mirror run_simulation,
    # roughness expressed in um). NOTE: transport variant is a per-propellant
    # property (propellant editor / .ric), not a run-control param.
    schema.setProperties({
        't_max': 10.0, 'P_cutoff': 0.5e6, 'cfl_target': 0.3,
        'roughness': 50.0, 'kappa': 0.45, 'T_ignition': 850.0,
    })
    return schema


def simulate_motor(motor, igniter_pyrogen=None, callback=None, **sim_overrides):
    """Run srm_1d on an openMotor ``Motor`` and return an openMotor
    ``SimulationResult``. Reusable headlessly without the registry.

    The igniter is read from the motor's own ``data.igniter`` block (Phase 4
    self-describing motors, now carried by openMotor's ``Motor`` — Phase 6
    task 4). Passing ``igniter_pyrogen='<name>'`` forces a named library
    pyrogen with auto chamber sizing instead, overriding the motor's block.
    Motors that predate the block fall back to the BPNV default.

    ``callback`` (the openMotor progress callback) is driven live: the @njit
    time loop publishes a 0..1 progress metric into a shared array while it
    runs (GIL released via ``nogil``), a poller thread here forwards it to
    ``callback``, and a truthy callback return requests a cooperative cancel.
    """
    from .openmotor_adapter import (
        ric_to_sim_args, build_pyrogen_chamber, load_pyrogen, load_igniter,
    )
    from .simulation import run_simulation
    from .nozzle import compute_motor_performance

    motor_dict = motor.getDict()
    # Transport is read from the propellant's per-tab schema (gas_props=None);
    # a sentinel/missing transport hard-faults (D7), consistent with the .ric path.
    args = ric_to_sim_args(motor_dict, gas_props=None, **sim_overrides)
    geo = args.pop('geo')
    prop = args.pop('propellant')
    nozzle = args['nozzle']

    # Igniter resolution: explicit override > motor's own block > BPNV default.
    if igniter_pyrogen is not None:
        pyro, sizing = load_pyrogen(igniter_pyrogen), {}
    else:
        loaded = load_igniter(motor_dict)
        if loaded is not None:
            pyro, sizing = loaded
        else:
            pyro, sizing = load_pyrogen('bpnv'), {}
    args['pyrogen_chamber'] = build_pyrogen_chamber(pyro, geo, nozzle, **sizing)
    args.setdefault('verbose', False)

    # v0.8.x station-viz: snapshots are the source of BOTH the per-cell axial
    # payload and the per-grain channels, so the run_simulation default 0.2 s
    # interval is far too coarse for the GUI — it heavily interpolates the
    # per-grain channels and under-samples the axial fields, making the station
    # traces jagged and eating the ignition transient. Target ~2000 frames over
    # the run (>= 5 ms) so the axial pressure tracks the full-resolution head
    # pressure. Short runs naturally get fewer frames.
    t_max_eff = float(args.get('t_max', 10.0))
    args.setdefault('snapshot_interval', max(0.005, t_max_eff / 2000.0))

    if callback is None:
        result = run_simulation(geo, prop, **args)
    else:
        result = _run_with_progress(run_simulation, geo, prop, args, callback)

    perf = compute_motor_performance(result, nozzle, prop)

    # openMotor's post-run consumers — motor-stats ``getPortRatio`` and the
    # grain burnback cross-section in ``resultsWidget`` — call grain-geometry
    # methods (``getPortArea`` → ``getFaceArea``) that need each grain's FMM
    # regression map / ``faceArea``, which are populated ONLY by
    # ``simulationSetup``. srm_1d runs its own FMM bridge and never calls it,
    # so the motor's own grain objects keep ``faceArea = None`` and crash
    # ``getFaceArea`` for FMM-class grains (Finocyl/Star/Custom/...). Mirror
    # openMotor's solver setup loop. (Base ``Grain.simulationSetup`` is a
    # no-op, so this is safe for every grain type.)
    for grain in motor.grains:
        grain.simulationSetup(motor.config)

    return _result_to_om_simresult(motor, result, perf, geo, prop)


def _run_with_progress(run_simulation, geo, prop, args, callback,
                       poll_interval=0.05):
    """Run ``run_simulation`` in a worker thread while polling its shared
    ``progress_state`` array and forwarding progress to ``callback``. A truthy
    callback return sets the cancel flag the @njit loop reads each step.

    The loop releases the GIL (``nogil=True``), so this thread runs
    concurrently and sees progress update live. Returns the worker's result;
    re-raises any exception the worker hit."""
    progress_state = np.zeros(2, dtype=np.float64)
    args['progress_state'] = progress_state
    holder = {}

    def _worker():
        try:
            holder['result'] = run_simulation(geo, prop, **args)
        except BaseException as exc:  # propagate to the caller's thread
            holder['error'] = exc

    # Progress display has three regimes, because no single signal predicts the
    # wall-clock cost of the whole run:
    #   • WARMUP — before the @njit loop publishes anything (JIT compile + the
    #     first geometry/burn setup), creep slowly to WARMUP_CAP so the bar
    #     shows life instead of sitting frozen at 0.
    #   • BURN (phys < TAIL_GATE) — follow the physics fraction. The in-loop
    #     metric is now web-depth-weighted (total web burned / total web), so it
    #     tracks the BULK burn and no longer snaps to the gate the instant a
    #     thin-web OD-taper tip cell burns through.
    #   • TAIL (phys >= TAIL_GATE) — the end-of-burn / depressurization. Keep
    #     following the (still-moving) physics tail, but add an ASYMPTOTIC
    #     wall-clock floor toward TAIL_CAP so the bar keeps creeping even if the
    #     physics value flattens — always moving, never a hard stall at a fixed
    #     %, decelerating over a horizon scaled by the burn-to-gate time.
    #     Everything is capped at TAIL_CAP while the worker runs; completion
    #     snaps to 1.0.
    WARMUP_CAP = 0.05
    WARMUP_DUR = 8.0          # ~JIT-compile time; creep 0 -> WARMUP_CAP over it
    TAIL_GATE = 0.9
    TAIL_CAP = 0.995          # asymptote; only true completion fills the bar
    displayed = 0.0
    tail_start = None
    tail_tau = None
    start = _time.monotonic()

    worker = threading.Thread(target=_worker, daemon=True)
    worker.start()
    while worker.is_alive():
        phys = float(progress_state[0])
        now = _time.monotonic()
        if phys <= 0.0:
            # Setup / JIT-compile: indeterminate slow creep (bounded).
            target = min(WARMUP_CAP, (now - start) / WARMUP_DUR * WARMUP_CAP)
        else:
            target = phys                  # follow the (representative) burn
            if phys >= TAIL_GATE:
                if tail_start is None:     # entering the tail — calibrate
                    tail_start = now
                    tail_tau = max(1.5, 0.5 * (now - start))
                # Asymptotic anti-stall floor — decelerates but never stops.
                floor = TAIL_GATE + (TAIL_CAP - TAIL_GATE) * (
                    1.0 - _math.exp(-(now - tail_start) / tail_tau))
                if floor > target:
                    target = floor
            if target > TAIL_CAP:
                target = TAIL_CAP          # only true completion fills the bar
        if target > displayed:
            displayed = target             # monotonic
        if callback(displayed):
            progress_state[1] = 1.0         # request cooperative cancel
        _time.sleep(poll_interval)
    worker.join()
    if progress_state[1] <= 0.5:
        callback(1.0)                       # completed → fill the bar

    if 'error' in holder:
        raise holder['error']
    return holder['result']


def _per_grain_series(result, geo, prop):
    """Aggregate srm_1d's per-cell axial snapshot fields to per-grain series,
    sampled on the snapshot time base (Phase 6 task 3).

    Returns ``(snap_times, series, vol_loading)`` where ``series`` maps each
    openMotor multi-value channel name to an ``(n_frames, n_grains)`` array
    and ``vol_loading`` is the ``(n_frames,)`` scalar volume-loading series:

    - ``mass``       — solid propellant mass per grain, integrated from the
      per-cell port diameter: ρ·dx·Σ (A_outer − A_port).
    - ``regression`` / ``web`` — per-grain radial regression / remaining web
      (already computed per-grain in ``result['grains']``).
    - ``machNumber`` — peak core Mach over each grain's cells.
    - ``massFlow``   — cumulative downstream mass-loss rate (forward→aft sum
      of each grain's −d(mass)/dt), matching openMotor's convention.
    - ``massFlux``   — ``massFlow`` divided by the grain's aft port area.
    - ``vol_loading`` — 100 · (solid propellant volume / grain bounding
      volume), the per-step scalar openMotor reports as Volume Loading.

    Returns ``None`` if the result carries no snapshots (nothing to map).
    """
    snapshots = result['snapshots'] if 'snapshots' in result else None
    grains = result['grains'] if 'grains' in result else None
    if not snapshots or not grains:
        return None

    ga = geo.compile_geometry_arrays()
    cell_seg = np.asarray(ga['cell_segment_id'])
    n_grains = len(grains)
    n_frames = len(snapshots)
    dx = geo.dx
    A_outer = np.pi / 4.0 * geo.D_outer ** 2
    rho = float(prop.rho_propellant)

    # Per-grain cell masks + aft-most cell index (for the aft port area).
    seg_cells = [np.flatnonzero(cell_seg == k) for k in range(n_grains)]
    seg_aft = [int(cells[-1]) if cells.size else -1 for cells in seg_cells]

    snap_times = np.array([snap['t'] for snap in snapshots])
    mass = np.zeros((n_frames, n_grains))
    mach = np.zeros((n_frames, n_grains))
    aft_port_area = np.zeros((n_frames, n_grains))
    for s, snap in enumerate(snapshots):
        d_port = np.asarray(snap['D_port'])
        mach_field = np.asarray(snap['Mach'])
        a_port = np.pi / 4.0 * d_port ** 2
        solid = np.maximum(A_outer - a_port, 0.0)
        for k in range(n_grains):
            cells = seg_cells[k]
            if cells.size:
                mass[s, k] = rho * dx * float(np.sum(solid[cells]))
                mach[s, k] = float(np.max(mach_field[cells]))
            if seg_aft[k] >= 0:
                aft_port_area[s, k] = float(a_port[seg_aft[k]])

    regression = np.stack(
        [np.asarray(grains[k]['regression']) for k in range(n_grains)], axis=1)
    web = np.stack(
        [np.asarray(grains[k]['web']) for k in range(n_grains)], axis=1)

    # massFlow: each grain's mass-loss rate, accumulated forward→aft so a
    # grain reports the total flow passing through it (openMotor convention).
    dmass = np.zeros_like(mass)
    if n_frames > 1:
        dt = np.diff(snap_times)
        dt[dt <= 0.0] = np.nan  # guard against duplicate snapshot times
        rate = np.maximum(-(np.diff(mass, axis=0)) / dt[:, None], 0.0)
        dmass[1:] = np.nan_to_num(rate)
    mass_flow = np.cumsum(dmass, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        mass_flux = np.where(aft_port_area > 0.0,
                             mass_flow / aft_port_area, 0.0)

    # Volume loading: solid propellant volume / grain bounding volume. The
    # bounding volume is the grain cells' cylinder (A_outer · dx · n_cells),
    # matching openMotor's "fraction of the chamber occupied by propellant".
    n_grain_cells = int(np.count_nonzero(cell_seg >= 0))
    bounding_volume = A_outer * dx * max(n_grain_cells, 1)
    solid_volume = mass.sum(axis=1) / rho
    vol_loading = 100.0 * solid_volume / bounding_volume

    return snap_times, {
        'mass': mass, 'regression': regression, 'web': web,
        'machNumber': mach, 'massFlow': mass_flow, 'massFlux': mass_flux,
    }, vol_loading


def _decimate_indices(n, p_head, thrust, max_points):
    """Return a sorted, unique index array selecting at most ~``max_points``
    samples from a length-``n`` series, always including the first, last, and
    peak-pressure / peak-thrust samples so the decimated trace preserves peak
    P, burn time, and (closely) impulse. Returns all indices when ``n`` is
    already under the cap."""
    if n <= max_points:
        return np.arange(n)
    stride = int(np.ceil(n / max_points))
    keep = set(range(0, n, stride))
    keep.update((0, n - 1))
    if p_head is not None and p_head.size:
        keep.add(int(np.argmax(p_head)))
    if thrust is not None and thrust.size:
        keep.add(int(np.argmax(thrust)))
    return np.array(sorted(keep))


def _axial_payload_for_gui(result):
    """Build the per-station axial payload the GUI station panel consumes,
    as **plain, GUI-friendly structures** (numpy arrays + dicts; no srm_1d
    types) attached to the SimulationResult as ``sr.srm1d_axial``.

    This is the v0.8.x station-viz data contract (design phase 3): the
    capability-gated station panel slices these per-cell field matrices on
    demand instead of routing through openMotor's fixed per-grain channels.
    The default fore/mid/aft station model is embedded so the GUI needs no
    srm_1d import to populate its selector.

    Returns ``None`` when the result lacks the axial contract (e.g. results
    produced before v0.8.x), so the GUI cleanly falls back to per-grain mode.
    """
    if 'cell_segment_id' not in result or not result.get('snapshots'):
        return None
    from .station_viz import build_axial_payload, default_stations
    # Carry the snapshots at (nearly) full resolution — the snapshot interval
    # (set in simulate_motor) already bounds the frame count (~2000), so a high
    # cap avoids a second decimation stage that would re-introduce jaggedness.
    payload = build_axial_payload(result, max_frames=4000)
    if payload is None:
        return None
    stations = default_stations(payload.cell_segment_id, payload.x_cell)
    return {
        'snap_times': payload.snap_times,
        'x_cell': payload.x_cell,
        'cell_segment_id': payload.cell_segment_id,
        'fields': payload.fields,
        # Roadmap #2 longitudinal motor-slice geometry (constant per run).
        'dx': payload.dx,
        'D_outer': payload.D_outer,
        'cell_D_outer': payload.cell_D_outer,   # per-cell OD-tapered casing
        'cell_wall_web': payload.cell_wall_web,
        'seg_geom': payload.seg_geom,   # per-segment end-face burnback
        'stations': [
            {
                'grain': s.grain, 'cell_index': s.cell_index,
                'position_m': s.position_m, 'active': s.active,
                'role': s.role, 'label': s.label,
            }
            for s in stations
        ],
    }


def _result_to_om_simresult(motor, result, perf, geo=None, prop=None):
    """Map srm_1d results (SimulationChannels) + performance into an
    openMotor SimulationResult. Populates the scalar per-step channels and —
    when ``geo``/``prop`` are supplied — the per-grain multi-value channels
    (interpolated from the snapshot time base onto the per-step time grid).

    Also attaches ``sr.srm1d_axial`` (the per-station axial payload) when the
    result carries the v0.8.x axial contract, so the GUI's station panel can
    slice per-cell fields without re-deriving them."""
    _, _, SimulationResult = _om()
    sr = SimulationResult(motor)

    time = np.asarray(result['time'])
    p_head = np.asarray(result['P_head'])
    p_exit = np.asarray(result['P_exit'])
    d_throat = np.asarray(result['D_throat'])
    thrust = np.asarray(perf['thrust'])
    kn = np.asarray(result['Kn']) if 'Kn' in result else None
    d0 = float(d_throat[0]) if d_throat.size else 0.0

    # Per-grain channels live on the sparser snapshot time base; interpolate
    # each grain's series onto the per-step time grid so every channel shares
    # the 'time' channel's length (an openMotor SimulationResult invariant).
    grain_step = None
    vol_loading_step = None
    if geo is not None and prop is not None:
        pg = _per_grain_series(result, geo, prop)
        if pg is not None:
            snap_times, series, vol_loading = pg
            grain_step = {
                name: np.column_stack([
                    np.interp(time, snap_times, arr[:, k])
                    for k in range(arr.shape[1])
                ]) for name, arr in series.items()
            }
            vol_loading_step = np.interp(time, snap_times, vol_loading)

    # Decimate to a GUI-friendly sample count, keeping ballistics-critical
    # samples (peak pressure, peak thrust, first/last).
    idx = _decimate_indices(time.size, p_head, thrust, GUI_MAX_POINTS)

    for i in idx:
        sr.channels['time'].addData(float(time[i]))
        sr.channels['pressure'].addData(float(p_head[i]))
        sr.channels['exitPressure'].addData(float(p_exit[i]))
        sr.channels['force'].addData(float(thrust[i]))
        sr.channels['dThroat'].addData(float(d_throat[i] - d0))
        if kn is not None:
            sr.channels['kn'].addData(float(kn[i]))
        if vol_loading_step is not None:
            sr.channels['volumeLoading'].addData(float(vol_loading_step[i]))
        if grain_step is not None:
            for name, arr in grain_step.items():
                sr.channels[name].addData(tuple(float(v) for v in arr[i]))

    summary = result['summary'] if 'summary' in result else {}
    # Run-health watchdog: surface failed runs instead of silently plotting
    # nothing. A run is UNHEALTHY when it (a) numerically collapsed (code 4),
    # (b) was canceled (code 5), or (c) reached t_max (code 0) having never
    # pressurized — i.e. it never ignited (peak pressure near ambient). Note
    # ``t_burn`` is just the run duration, so it can't detect non-ignition;
    # peak pressure is the right signal. The ~0.3 MPa floor (≈3× ambient) is
    # far below any real solid-motor chamber pressure.
    term_code = summary.get('termination_code', 0)
    p_peak_mpa = float(summary.get('P_peak', 0.0) or 0.0) / 1.0e6
    canceled = term_code == 5
    collapsed = term_code == 4
    no_ignition = (term_code == 0) and (p_peak_mpa < 0.3)
    sr.success = not (canceled or collapsed or no_ignition)
    if not sr.success and not canceled:
        from motorlib.simResult import SimAlert, SimAlertLevel, SimAlertType
        why = ("the motor never ignited" if no_ignition
               else "the solution collapsed")
        sr.addAlert(SimAlert(
            SimAlertLevel.ERROR, SimAlertType.VALUE,
            "srm_1d: no usable trace -- {} (P_peak {:.2f} MPa). Check the "
            "igniter / propellant / geometry.".format(why, p_peak_mpa),
            'srm_1d'))

    # v0.8.x station-viz: attach the per-station axial payload (capability-
    # gated side attribute; the GUI station panel consumes it, QS ignores it).
    axial = _axial_payload_for_gui(result)
    if axial is not None:
        sr.srm1d_axial = axial
    return sr


def register():
    """Register the srm_1d transient solver with openMotor's registry."""
    SolverPlugin, register_solver, _ = _om()

    class Srm1dTransientSolver(SolverPlugin):
        name = SOLVER_NAME
        capabilities = {
            'transient': True,
            'axial_fields': True,
            'needs_transport': True,
            'igniter': True,  # uses the pyrogen igniter (drives the GUI Igniter row)
        }

        def simulate(self, motor, config=None, callback=None):
            overrides = dict(config) if isinstance(config, dict) else {}
            # The config screen edits roughness in micrometers (its own unit
            # category); run_simulation wants metres. Convert at this boundary
            # so simulate_motor / run_simulation stay SI.
            if 'roughness' in overrides:
                overrides['roughness'] = overrides['roughness'] * 1.0e-6
            return simulate_motor(motor, callback=callback, **overrides)

        def get_config_schema(self):
            return _transient_config_schema()

        def upgrade_solver_config(self, values):
            if not isinstance(values, dict):
                return values, []
            new = dict(values)
            upgraded = []
            if 'T_ignition' in new and abs(new['T_ignition'] - 850.0) < 1e-4:
                new['T_ignition'] = 756.0
                upgraded.append('T_ignition')
            if 'kappa' in new and abs(new['kappa'] - 0.45) < 1e-4:
                new['kappa'] = 0.44
                upgraded.append('kappa')
            if 'roughness' in new and abs(new['roughness'] - 50.0) < 1e-4:
                new['roughness'] = 32.0
                upgraded.append('roughness')
            return new, upgraded

        def upgrade_igniter_block(self, block):
            if not isinstance(block, dict) or 'pyrogen' not in block:
                return block, False
            pyro = block['pyrogen']
            if not isinstance(pyro, dict) or pyro.get('name') != 'BPNV':
                return block, False
            if 'gas_mass_fraction' in pyro:
                return block, False
            import copy
            import math
            legacy_numeric = {
                'a': 2.0e-5, 'n': 0.50, 'rho': 1700.0,
                'T_flame': 2800.0, 'M': 0.030, 'gamma': 1.25,
                'impetus_W': 5000.0, 'heat_flux_cal_cm2_s': 69.4,
                'kappa_jet': 8.0,
                'particle_diameter_m': 0.0032, 'particle_LD_ratio': 1.0,
                'pellet_emissivity': 0.7, 'radiation_absorption_length_m': 1.0,
            }
            legacy_str = {
                'form': 'pellets', 'heat_delivery_mode': 'demar',
            }
            for k, v in legacy_str.items():
                if pyro.get(k) != v:
                    return block, False
            for k, v in legacy_numeric.items():
                val = pyro.get(k)
                if val is None or not math.isclose(float(val), v, rel_tol=1e-3, abs_tol=1e-5):
                    return block, False
            from .openmotor_adapter import _pyrogen_to_block, load_pyrogen
            current = _pyrogen_to_block(load_pyrogen('bpnv'))
            new_block = copy.deepcopy(block)
            new_block['pyrogen'] = current
            return new_block, True

        def get_default_igniter(self):
            from .openmotor_adapter import default_igniter_block
            return default_igniter_block('bpnv')

        def get_igniter_materials(self):
            from .openmotor_adapter import _pyrogen_to_block, load_pyrogen
            return [
                _pyrogen_to_block(load_pyrogen(name))
                for name in ('bpnv', 'mtv', 'mtv_fast', 'thermite_cual')
            ]

        def get_igniter_material_schema(self):
            _setup_openmotor_path()
            from motorlib.properties import (
                PropertyCollection, FloatProperty, StringProperty,
                EnumProperty, BooleanProperty,
            )
            from .openmotor_adapter import load_pyrogen

            schema = PropertyCollection()
            schema.props['name'] = StringProperty('Pyrogen Name')
            schema.props['a'] = FloatProperty('Burn Rate Coefficient', '', 0.0, 1.0)
            schema.props['n'] = FloatProperty('Burn Rate Exponent', '', 0.0, 2.0)
            schema.props['rho'] = FloatProperty('Density', 'kg/m^3', 1.0, 10000.0)
            schema.props['T_flame'] = FloatProperty('Flame Temperature', 'K', 300.0, 6000.0)
            schema.props['M'] = FloatProperty('Molar Mass', 'kg/mol', 0.001, 0.5)
            schema.props['gamma'] = FloatProperty('Specific Heat Ratio', '', 1.0, 3.0)
            schema.props['impetus_W'] = FloatProperty('Impetus', 'W', 0.0, 50000.0)
            schema.props['heat_flux_cal_cm2_s'] = FloatProperty(
                'Heat Flux', 'cal/(cm^2*s)', 0.0, 1000.0)
            schema.props['particle_diameter_m'] = FloatProperty(
                'Particle Diameter', 'm', 1e-6, 0.1)
            schema.props['particle_LD_ratio'] = FloatProperty(
                'Particle L/D', '', 0.1, 100.0)
            schema.props['form'] = EnumProperty(
                'Physical Form', ['powder', 'pellets', 'chunks', 'solid'])
            schema.props['heat_delivery_mode'] = EnumProperty(
                'Heat Delivery', ['demar', 'radiation', 'none'])
            schema.props['pellet_emissivity'] = FloatProperty(
                'Pellet Emissivity', '', 0.0, 1.0)
            schema.props['radiation_absorption_length_m'] = FloatProperty(
                'Radiation Absorption Length', 'm', 1e-3, 100.0)
            schema.props['kappa_jet'] = FloatProperty(
                'Jet Decay Factor', '', 0.0, 100.0)
            schema.props['gas_mass_fraction'] = FloatProperty(
                'Gas Mass Fraction', '', 0.0, 1.0)
            schema.setProperties({
                'name': 'BPNV',
                'a': 2e-5,
                'n': 0.5,
                'rho': 1700.0,
                'T_flame': 2800.0,
                'M': 0.03,
                'gamma': 1.25,
                'impetus_W': 5000.0,
                'heat_flux_cal_cm2_s': 69.4,
                'particle_diameter_m': 0.0032,
                'particle_LD_ratio': 1.0,
                'form': 'pellets',
                'heat_delivery_mode': 'demar',
                'pellet_emissivity': 0.7,
                'radiation_absorption_length_m': 1.0,
                'kappa_jet': 8.0,
                'gas_mass_fraction': 1.0,
            })
            return schema

        def get_igniter_chamber_schema(self):
            _setup_openmotor_path()
            from motorlib.properties import (
                PropertyCollection, FloatProperty, EnumProperty,
            )
            schema = PropertyCollection()
            schema.props['mass'] = FloatProperty('Pyrogen Mass', 'kg', -1.0, 100.0)
            schema.props['throat_area'] = FloatProperty(
                'Vent Throat Area', 'm^2', -1.0, 1.0)
            schema.props['volume'] = FloatProperty(
                'Chamber Volume', 'm^3', -1.0, 1.0)
            schema.props['burn_area'] = FloatProperty(
                'Pyrogen Burn Area', 'm^2', -1.0, 100.0)
            schema.props['cartridge_length_m'] = FloatProperty(
                'Cartridge Length', 'm', -1.0, 100.0)
            schema.props['basket_fill_fraction'] = FloatProperty(
                'Basket Fill Fraction', '', 0.0, 1.0)
            schema.props['pellet_packing_fraction'] = FloatProperty(
                'Pellet Packing Fraction', '', 0.0, 1.0)
            schema.props['injection_topology'] = EnumProperty(
                'Injection Topology',
                ['forward_plenum', 'head_basket', 'aft_basket'])
            schema.props['burn_law'] = EnumProperty(
                'Burn Law', ['0d'])
            schema.setProperties({
                'mass': -1.0,
                'throat_area': -1.0,
                'volume': -1.0,
                'burn_area': -1.0,
                'cartridge_length_m': -1.0,
                'basket_fill_fraction': 0.5,
                'pellet_packing_fraction': 0.6,
                'injection_topology': 'forward_plenum',
                'burn_law': '0d',
            })
            return schema

    return register_solver(Srm1dTransientSolver())


# Register on import so `import srm_1d.srm1d_plugin` makes the solver
# selectable via motorlib.solvers.get_solver('srm_1d-transient').
try:
    register()
except Exception:  # pragma: no cover - openMotor checkout may be absent
    pass
