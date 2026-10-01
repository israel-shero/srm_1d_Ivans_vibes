"""
channels.py — v0.8.0 Phase 1 result channel model
==================================================

Aligns srm_1d's simulation results to openMotor's ``LogChannel`` shape so
the frontend-integration work (the ``motorlib`` plugin in Phase 5, the
unit-aware ``plot_channels`` in Phase 2) can consume a single channel
object instead of the historical ``run_simulation()`` results dict.

Two channel types (DESIGN.md D5):

- :class:`Channel` — a scalar-per-step or per-grain-per-step series. Mirrors
  ``motorlib.simResult.LogChannel`` (``name``, ``unit``, ``getData(unit)``,
  ``getMax``/``getMin``/``getLast``/``getAverage``). openMotor's list-type
  channels are per-GRAIN; this type covers both the scalar and per-grain
  cases.
- :class:`AxialChannel` — the srm_1d-specific extension: a per-CELL axial
  field over time (``time × N_cells``) plus the axial coordinate
  ``x_cells``. openMotor has no per-cell concept; the GUI ignores these
  until the visualization work (v0.8.x) consumes them.

:func:`build_channels` maps the ``run_simulation()`` dict into a
:class:`SimulationChannels` container **without touching the numerics** —
the data is re-shaped, never recomputed, so results are byte-for-byte
identical to the dict. The dict remains the source of truth for this phase;
consumers migrate onto channels incrementally (Phase 2 plotting first).

Unit conversion defers to openMotor's ``motorlib.units.convert`` (a
GUI-free, dependency-light table), imported lazily so the core stays
importable without the openMotor checkout for raw (un-converted) access.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np


# ----------------------------------------------------------------------
# Lazy unit conversion — defers to openMotor's motorlib.units table.
# ----------------------------------------------------------------------
def _convert(value, origin_unit: str, dest_unit: str):
    """Convert ``value`` from ``origin_unit`` to ``dest_unit``.

    Identity when the units match (no openMotor needed). Otherwise defers
    to ``motorlib.units.convert`` via the same checkout bootstrap that
    ``fmm_grain`` uses. ``value`` may be a scalar or a numpy array; the
    conversion is a pure ratio so it broadcasts.
    """
    if origin_unit == dest_unit:
        return value
    from .fmm_grain import _setup_openmotor_path
    _setup_openmotor_path()
    from motorlib import units  # type: ignore
    ratio = units.getConversion(origin_unit, dest_unit)
    return value * ratio


class Channel:
    """A scalar or per-grain time-series, aligned to openMotor ``LogChannel``.

    Parameters
    ----------
    name : str
        Human-readable channel name (e.g. ``"Head Pressure"``).
    unit : str
        openMotor unit string the data is stored in (e.g. ``"Pa"``). May be
        a unit outside openMotor's conversion table (``"K"``, ``"W"``,
        ``""``); those support only identity ``getData``.
    data : np.ndarray
        1-D ``(n_steps,)`` for a scalar channel, or 2-D ``(n_steps,
        n_grains)`` for a per-grain channel.
    per_grain : bool
        Whether ``data`` is per-grain (2-D). Mirrors ``LogChannel``'s
        ``valueType in (list, tuple)`` distinction.
    """

    def __init__(self, name: str, unit: str, data: np.ndarray,
                 per_grain: bool = False):
        self.name = name
        self.unit = unit
        self.per_grain = per_grain
        self.data = np.asarray(data)

    def getData(self, unit: Optional[str] = None) -> np.ndarray:
        """Return the channel data, converted to ``unit`` if given."""
        if unit is None or unit == self.unit:
            return self.data
        return _convert(self.data, self.unit, unit)

    def getPoint(self, i: int):
        """Return the datapoint (scalar or per-grain row) at index ``i``."""
        return self.data[i]

    def getLast(self):
        return self.data[-1]

    def getMax(self) -> float:
        """Largest single value in the channel (matches LogChannel)."""
        return float(np.max(self.data))

    def getMin(self) -> float:
        return float(np.min(self.data))

    def getAverage(self) -> float:
        if self.per_grain:
            raise NotImplementedError("Average not supported for per-grain channels")
        return float(np.mean(self.data))

    def __len__(self) -> int:
        return self.data.shape[0]

    def __repr__(self) -> str:
        shape = "x".join(str(d) for d in self.data.shape)
        kind = "per-grain" if self.per_grain else "scalar"
        return f"Channel({self.name!r}, unit={self.unit!r}, {kind}, [{shape}])"


class AxialChannel:
    """A per-cell axial field over time — srm_1d's extension to the channel
    model (DESIGN.md D5).

    Parameters
    ----------
    name : str
        Human-readable field name (e.g. ``"Bore Pressure"``).
    unit : str
        openMotor unit string the data is stored in.
    times : np.ndarray
        ``(n_frames,)`` sample times [s] (srm_1d's snapshot times — sparser
        than the solver step series).
    data : np.ndarray
        ``(n_frames, n_cells)`` field values.
    x_cells : np.ndarray
        ``(n_cells,)`` axial cell-center coordinate [m].
    """

    def __init__(self, name: str, unit: str, times: np.ndarray,
                 data: np.ndarray, x_cells: np.ndarray):
        self.name = name
        self.unit = unit
        self.times = np.asarray(times)
        self.data = np.asarray(data)
        self.x_cells = np.asarray(x_cells)
        if self.data.ndim != 2:
            raise ValueError(
                f"AxialChannel {name!r} data must be 2-D (frames x cells); "
                f"got shape {self.data.shape}")
        if self.data.shape[0] != self.times.shape[0]:
            raise ValueError(
                f"AxialChannel {name!r}: frames ({self.data.shape[0]}) != "
                f"times ({self.times.shape[0]})")
        if self.data.shape[1] != self.x_cells.shape[0]:
            raise ValueError(
                f"AxialChannel {name!r}: cells ({self.data.shape[1]}) != "
                f"x_cells ({self.x_cells.shape[0]})")

    @property
    def n_frames(self) -> int:
        return self.data.shape[0]

    @property
    def n_cells(self) -> int:
        return self.data.shape[1]

    def getData(self, unit: Optional[str] = None) -> np.ndarray:
        """Return the ``(n_frames, n_cells)`` field, converted if ``unit`` given."""
        if unit is None or unit == self.unit:
            return self.data
        return _convert(self.data, self.unit, unit)

    def getFrame(self, i: int, unit: Optional[str] = None) -> np.ndarray:
        """Return the axial profile at frame index ``i`` (one value per cell)."""
        return self.getData(unit)[i]

    def getCell(self, j: int, unit: Optional[str] = None) -> np.ndarray:
        """Return the time history at cell index ``j`` (one value per frame)."""
        return self.getData(unit)[:, j]

    def getMax(self) -> float:
        return float(np.max(self.data))

    def getMin(self) -> float:
        return float(np.min(self.data))

    def __repr__(self) -> str:
        return (f"AxialChannel({self.name!r}, unit={self.unit!r}, "
                f"{self.n_frames} frames x {self.n_cells} cells)")


class SimulationChannels:
    """Container for a simulation's channels, the channel-model analogue of
    ``motorlib.simResult.SimulationResult``.

    Attributes
    ----------
    channels : dict[str, Channel]
        Scalar + per-grain time series keyed by dict name (``"P_head"``, …).
    axial : dict[str, AxialChannel]
        Per-cell axial fields keyed by field name (``"P"``, ``"T"``, …).
    summary : dict
        The scalar summary block, passed through unchanged.
    extras : dict
        Non-channel result entries retained for back-compat / diagnostics
        (final-state arrays, species registry, per-cell ignition time).
    """

    def __init__(self):
        self.channels: Dict[str, Channel] = {}
        self.axial: Dict[str, AxialChannel] = {}
        self.summary: dict = {}
        self.extras: dict = {}
        # The original run_simulation results dict. SimulationChannels is a
        # drop-in for it: item access / iteration proxy to ``raw`` so legacy
        # ``result['P_head']`` / ``result['summary']`` code keeps working,
        # while ``.channels`` / ``.axial`` expose the unit-aware channel API.
        self.raw: dict = {}

    # --- legacy mapping interface (drop-in for the results dict) ---
    def __getitem__(self, key):
        return self.raw[key]

    def __setitem__(self, key, value):
        self.raw[key] = value

    def __contains__(self, key) -> bool:
        return key in self.raw

    def __iter__(self):
        return iter(self.raw)

    def __len__(self) -> int:
        return len(self.raw)

    def get(self, key, default=None):
        return self.raw.get(key, default)

    def keys(self):
        return self.raw.keys()

    def items(self):
        return self.raw.items()

    def values(self):
        return self.raw.values()

    # --- channel API ---
    def channel(self, name: str) -> Channel:
        """Return the unit-aware :class:`Channel` for ``name``."""
        return self.channels[name]

    def add(self, name: str, channel: Channel) -> Channel:
        """Insert a scalar/per-grain channel (e.g. derived thrust/Isp from
        ``compute_motor_performance``) so it can flow through the generic
        plotting path. Returns the channel for chaining."""
        self.channels[name] = channel
        return channel

    def names(self) -> List[str]:
        return list(self.channels.keys())

    def axialNames(self) -> List[str]:
        return list(self.axial.keys())

    def __repr__(self) -> str:
        return (f"SimulationChannels({len(self.channels)} channels, "
                f"{len(self.axial)} axial fields)")


# ----------------------------------------------------------------------
# Dict -> channel mapping. Units are openMotor unit strings where a
# conversion exists; srm_1d-only quantities (power, enthalpy, per-length
# source rates) carry a descriptive unit string with identity conversion
# only — they exist for labelling and future viz, not GUI unit switching.
# ----------------------------------------------------------------------

# Per-step scalar series: dict key -> (display name, unit).
_SCALAR_UNITS = {
    'time': ('Time', 's'),
    'P_head': ('Head Pressure', 'Pa'),
    'P_exit': ('Nozzle Exit Pressure', 'Pa'),
    'D_throat': ('Throat Diameter', 'm'),
    'Kn': ('Kn', ''),
    'massflow': ('Nozzle Mass Flow', 'kg/s'),
    'P_ig': ('Igniter Pressure', 'Pa'),
    'T_ig': ('Igniter Temperature', 'K'),
    'mdot_ig': ('Igniter Mass Flow', 'kg/s'),
    'm_pyrogen': ('Pyrogen Mass Remaining', 'kg'),
    'gas_sensible_energy_before': ('Gas Sensible Energy (pre-step)', 'J'),
    'gas_sensible_energy': ('Gas Sensible Energy', 'J'),
    'gas_sensible_dE_dt': ('Gas Sensible dE/dt', 'W'),
    'normal_sidewall_thermal_power': ('Sidewall Thermal Power (normal)', 'W'),
    'erosive_sidewall_thermal_power': ('Sidewall Thermal Power (erosive)', 'W'),
    'endface_thermal_power': ('Endface Thermal Power', 'W'),
    'pyrogen_gas_thermal_power': ('Pyrogen Gas Thermal Power', 'W'),
    'convective_scalar_flux_power': ('Convective Scalar Flux Power', 'W'),
    'nozzle_scalar_flux_power': ('Nozzle Scalar Flux Power', 'W'),
    'clipping_correction_power': ('Clipping Correction Power', 'W'),
    'pyrogen_enthalpy_power': ('Pyrogen Enthalpy Power', 'W'),
    'pyrogen_surface_heat_power': ('Pyrogen Surface Heat Power', 'W'),
    'gas_surface_heat_sink_power': ('Pyrogen-to-Wall Gas Debit', 'W'),
    'radiation_heat_power': ('Radiation Heat Power', 'W'),
    'radiation_sink_power': ('Radiation Sink Power', 'W'),
    'convective_wall_heat_sink_power': ('Convective Wall Heat Sink Power', 'W'),
    'nozzle_enthalpy_power': ('Nozzle Enthalpy Power', 'W'),
    'thermal_source_power': ('Thermal Source Power', 'W'),
    'energy_residual': ('Energy Residual', 'W'),
    'pyrogen_momentum_expected': ('Pyrogen Momentum (expected)', 'N'),
    'pyrogen_momentum_deposited': ('Pyrogen Momentum (deposited)', 'N'),
    'pyrogen_momentum_residual': ('Pyrogen Momentum Residual', 'N'),
    'dt': ('Timestep', 's'),
    'n_burning': ('Burning Cell Count', ''),
    'n_ignited': ('Ignited Cell Count', ''),
    'radiation_emitter_count': ('Radiation Emitter Count', ''),
    'radiation_receiver_count': ('Radiation Receiver Count', ''),
    'min_gas_temperature': ('Min Gas Temperature', 'K'),
    'max_gas_temperature': ('Max Gas Temperature', 'K'),
    'min_surface_temperature': ('Min Surface Temperature', 'K'),
    'max_surface_temperature': ('Max Surface Temperature', 'K'),
    'min_pressure': ('Min Pressure', 'Pa'),
    'max_pressure': ('Max Pressure', 'Pa'),
    'max_mach': ('Max Mach Number', ''),
}

# Axial snapshot fields: snapshot-dict key -> (display name, unit).
_AXIAL_UNITS = {
    'P': ('Bore Pressure', 'Pa'),
    'u': ('Face Velocity', 'm/s'),
    'Mach': ('Mach Number', ''),
    'T': ('Gas Temperature', 'K'),
    'r_total': ('Total Burn Rate', 'm/s'),
    'r_erosive': ('Erosive Burn Rate', 'm/s'),
    'D_port': ('Port Diameter', 'm'),
    'C_burn': ('Burn Perimeter', 'm'),
    'endface_msource': ('Endface Mass Source', 'kg/s'),
    'T_surf': ('Surface Temperature', 'K'),
    'mass_source': ('Mass Source', 'kg/(m*s)'),
    'thermal_source': ('Thermal Source', 'W/m'),
    'momentum_source': ('Momentum Source', 'N/m'),
    'pyrogen_surface_heat_flux': ('Pyrogen Surface Heat Flux', 'W/m^2'),
    'radiation_heat_flux': ('Radiation Heat Flux', 'W/m^2'),
    'is_burning': ('Is Burning', ''),
    'is_grain': ('Is Grain', ''),
}

# Per-grain history fields (over snapshot times): key -> (display name, unit).
_GRAIN_UNITS = {
    'regression': ('Regression Depth', 'm'),
    'web': ('Web', 'm'),
}

# Result-dict entries deliberately retained as ``extras`` rather than
# channelised (final-state arrays, registries, per-cell scalars).
_EXTRA_KEYS = (
    'ignition_time_by_cell', 'P_ambient',
    'Y_species_final', 'species_params', 'species_names',
    'rho_final', 'A_port_final',
    'gamma_mix_final', 'Cp_mix_final', 'R_mix_final', 'M_mix_final',
)


def as_channels(result) -> SimulationChannels:
    """Normalize a result to :class:`SimulationChannels`.

    Accepts either a ``run_simulation()`` results dict (re-shaped via
    :func:`build_channels` — pure, no recompute) or an existing
    :class:`SimulationChannels` (returned as-is). Lets consumer code accept
    both during the v0.8.0 migration: dict callers stay byte-for-byte
    identical, channel callers pass through.
    """
    if isinstance(result, SimulationChannels):
        return result
    return build_channels(result)


def build_channels(results: dict) -> SimulationChannels:
    """Map a ``run_simulation()`` results dict into a SimulationChannels.

    Pure re-shaping: every value comes straight from ``results`` (no
    recomputation), so the channels carry byte-for-byte identical data.
    Unknown keys are ignored (forward-compatible with new diagnostics);
    missing known keys are skipped (back-compatible with older dicts).
    """
    sc = SimulationChannels()
    sc.raw = results  # legacy drop-in access (item/iteration proxy)

    # --- Scalar per-step channels ---
    for key, (disp, unit) in _SCALAR_UNITS.items():
        if key in results and results[key] is not None:
            sc.channels[key] = Channel(disp, unit, np.asarray(results[key]))

    # --- Axial fields, pivoted from the snapshot list-of-dicts ---
    snapshots = results.get('snapshots') or []
    if snapshots:
        times = np.array([snap['t'] for snap in snapshots])
        x_cells = np.asarray(snapshots[0]['x'])
        for key, (disp, unit) in _AXIAL_UNITS.items():
            if key not in snapshots[0]:
                continue
            field = np.array([np.asarray(snap[key]) for snap in snapshots])
            sc.axial[key] = AxialChannel(disp, unit, times, field, x_cells)

    # --- Per-grain history channels, over the snapshot time base ---
    grains = results.get('grains') or []
    if grains:
        for key, (disp, unit) in _GRAIN_UNITS.items():
            if key not in grains[0]:
                continue
            # (n_frames, n_grains): transpose the per-grain (n_frames,) arrays.
            stacked = np.array([np.asarray(g[key]) for g in grains]).T
            sc.channels[key] = Channel(disp, unit, stacked, per_grain=True)

    # --- Summary + extras passthrough ---
    sc.summary = results.get('summary', {})
    for key in _EXTRA_KEYS:
        if key in results:
            sc.extras[key] = results[key]

    return sc
