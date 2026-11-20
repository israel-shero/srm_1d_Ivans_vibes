"""
plotting.py — Visualization Utilities
=======================================

Quick plotting for simulation results. This is a convenience module,
not part of the solver — it has no effect on numerics.

Usage:

    from srm_1d.openmotor_adapter import run_from_ric
    from srm_1d.nozzle import Nozzle, compute_motor_performance
    from srm_1d.plotting import plot_pressure, plot_thrust, plot_flow_snapshot

    result, perf, nozzle, geo, prop = run_from_ric(
        "srm_1d/motors/hasegawa_a.ric", pyrogen="bpnv"
    )

    # Pressure trace with experimental overlay
    plot_pressure(result, title="Motor A",
                  experimental=HASEGAWA_MOTOR_A_EXPERIMENTAL)

    # Thrust (requires nozzle post-processing)
    nozzle = Nozzle(D_throat=0.034, D_exit=0.050)
    perf = compute_motor_performance(result, nozzle, prop)
    plot_thrust(result, perf, title="Motor A")

    # Flow field snapshot at a specific time
    plot_flow_snapshot(result, t_target=2.0)
"""

import numpy as np

from srm_1d.run_artifacts import save_figure
from srm_1d.channels import (
    Channel, SimulationChannels, as_channels,
)

try:
    import matplotlib.pyplot as plt
    import matplotlib
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ================================================================
# Generic, unit-aware channel plotting (v0.8.0 Phase 2)
# ================================================================
# plot_channels is the generic engine the bespoke figures below are
# expressed through. It reads the channel model (any dict result is
# normalized via as_channels) and converts each channel to its display
# unit at draw time via Channel.getData(unit) — units are carried on the
# channels, conversion happens only at the boundary.

def plot_channels(result, names, *, display_units=None, time_name='time',
                  time_unit='s', overlay=False, labels=None, colors=None,
                  axes=None, title=None, save_path=None,
                  ylim_bottom_zero=True):
    """Plot one or more scalar channels vs a time channel, unit-aware.

    Parameters
    ----------
    result : SimulationChannels | dict | mapping of name->Channel
        Channel source. A run_simulation dict is re-shaped via
        ``as_channels``; a plain ``{name: Channel}`` mapping (e.g. wrapped
        performance) is accepted directly.
    names : str | sequence of str
        Channel name(s) to plot.
    display_units : dict, optional
        Per-channel display unit override (``{name: unit}``); defaults to
        each channel's stored unit. Conversion via ``Channel.getData``.
    time_name : str
        Name of the time channel (default ``'time'``).
    time_unit : str
        Display unit for the time axis (default ``'s'``).
    overlay : bool
        True → all channels on one axis; False (default) → one stacked
        subplot per channel sharing the x-axis.
    labels, colors : dict, optional
        Per-channel legend label / line color overrides.
    axes : matplotlib Axes or array, optional
        Draw onto these instead of creating a figure.
    title, save_path : str, optional
    ylim_bottom_zero : bool
        Clamp each y-axis bottom to 0 (default True).

    Returns
    -------
    fig, axes
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    sc = _normalize_channel_source(result)
    if isinstance(names, str):
        names = [names]
    names = list(names)
    display_units = dict(display_units or {})
    labels = dict(labels or {})
    colors = dict(colors or {})

    t = sc.channels[time_name].getData(time_unit)

    own_fig = axes is None
    if overlay:
        if own_fig:
            fig, ax = plt.subplots(figsize=(10, 6))
            axes_list = [ax]
        else:
            ax = axes if not isinstance(axes, (list, np.ndarray)) else axes[0]
            axes_list = [ax]
            fig = ax.get_figure()
        for name in names:
            ch = sc.channels[name]
            unit = display_units.get(name, ch.unit)
            y = ch.getData(unit)
            ax.plot(t, y, label=labels.get(name, ch.name),
                    color=colors.get(name))
        ax.set_xlabel(f'Time [{time_unit}]', fontsize=12)
        ax.legend(fontsize=10, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0)
        if ylim_bottom_zero:
            ax.set_ylim(bottom=0)
        if title:
            ax.set_title(title, fontsize=14)
        result_axes = axes_list
    else:
        n = len(names)
        if own_fig:
            fig, axarr = plt.subplots(n, 1, figsize=(10, 3.2 * n + 0.6),
                                      sharex=True, squeeze=False)
            axarr = axarr[:, 0]
        else:
            axarr = np.atleast_1d(axes)
            fig = axarr[0].get_figure()
        for i, name in enumerate(names):
            ch = sc.channels[name]
            unit = display_units.get(name, ch.unit)
            y = ch.getData(unit)
            ax = axarr[i]
            ax.plot(t, y, color=colors.get(name, 'b'), linewidth=2,
                    label=labels.get(name, ch.name))
            unit_sfx = f' [{unit}]' if unit else ''
            ax.set_ylabel(f'{ch.name}{unit_sfx}', fontsize=11)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(left=0)
            if ylim_bottom_zero:
                ax.set_ylim(bottom=0)
        axarr[-1].set_xlabel(f'Time [{time_unit}]', fontsize=12)
        if title:
            axarr[0].set_title(title, fontsize=14)
        result_axes = axarr

    if own_fig:
        plt.tight_layout()
    if save_path:
        save_figure(fig, save_path)
    return fig, result_axes


def _normalize_channel_source(result):
    """Accept a SimulationChannels, a results dict, or a {name: Channel}
    mapping; return a SimulationChannels."""
    if isinstance(result, SimulationChannels):
        return result
    if isinstance(result, dict) and result and all(
            isinstance(v, Channel) for v in result.values()):
        sc = SimulationChannels()
        sc.channels = dict(result)
        return sc
    return as_channels(result)


def _frames_from_channels(sc):
    """Reconstruct a list of per-snapshot field dicts from the axial
    channels, so the axial/flow plotters read the channel model rather
    than the raw results dict. Mirrors the historical snapshot dict
    structure (``'t'``, ``'x'``, and each field) for unchanged rendering.
    Returns ``[]`` when there are no axial channels.
    """
    if not sc.axial:
        return []
    any_ax = next(iter(sc.axial.values()))
    times = any_ax.times
    x_cells = any_ax.x_cells
    frames = []
    for i in range(any_ax.n_frames):
        frame = {'t': float(times[i]), 'x': x_cells}
        for key, ch in sc.axial.items():
            frame[key] = ch.getFrame(i)
        frames.append(frame)
    return frames


# ================================================================
# Experimental Data
# ================================================================
# Paste experimental data here as dict with 'time' and 'pressure'
# arrays. Pressure in MPa, time in seconds.
#
# To add your own data: either define a new dict following the same
# format, or use load_experimental_csv() to read from a file.

HASEGAWA_MOTOR_A_EXPERIMENTAL = {
    'label': 'Experimental (Hasegawa)',
    'time': np.array([
        0.000, 0.045, 0.053, 0.064, 0.068, 0.092, 0.116, 0.147, 0.206,
        0.358, 0.541, 0.953, 1.074, 1.214, 1.545, 1.822, 2.133,
        2.258, 2.348, 2.422, 2.555, 2.688, 2.754, 2.891, 3.036,
        3.168, 3.328, 3.481, 3.633, 3.852, 4.113, 4.297, 4.430,
        4.609, 4.675, 4.738,
    ]),
    'pressure': np.array([
        0.143, 2.299, 5.634, 6.219, 6.436, 5.922, 5.489, 5.297, 5.088,
        4.983, 5.014, 5.284, 5.388, 5.411, 5.578, 5.705, 5.887,
        5.951, 5.990, 5.854, 5.348, 4.762, 4.441, 3.935, 3.550,
        3.253, 2.923, 2.602, 2.256, 1.894, 1.219, 0.778, 0.464,
        0.191, 0.110, 0.110,

    ]),
}


# Zerox static-fire (Risky Batman V3, forward-Finocyl + aft-BATES).
# 58-point digitized pressure trace (MPa). time_offset shifts the
# experimental ignition event to align with sim t=0.
ZEROX_EXPERIMENTAL = {
    'label': 'Experimental (Zerox)',
    'time': np.array([
        0.000, 0.155, 0.292, 0.430, 0.574, 0.719, 0.874, 1.025, 1.179,
        1.316, 1.478, 1.610, 1.754, 1.898, 2.032, 2.170, 2.310, 2.446,
        2.588, 2.730, 2.867, 3.002, 3.147, 3.298, 3.438, 3.572, 3.711,
        3.852, 4.002, 4.143, 4.278, 4.417, 4.558, 4.700, 4.842, 4.974,
        5.115, 5.266, 5.410, 5.559, 5.694, 5.832, 5.985, 6.135, 6.276,
        6.426, 6.560, 6.717, 6.849, 6.991, 7.132, 7.284, 7.425, 7.559,
        7.701, 7.834, 7.973, 8.113,
    ]),
    'pressure': np.array([
        0.0000, 0.0239, 0.4591, 3.9938, 3.8630, 3.6684, 3.5838, 3.4890,
        3.4492, 3.3936, 3.3551, 3.2957, 3.2805, 3.2319, 3.2073, 3.1725,
        3.0986, 3.0645, 3.0449, 2.9944, 2.9640, 2.9394, 2.8680, 2.8017,
        2.6892, 2.6406, 2.5509, 2.4397, 2.3771, 2.2855, 2.2476, 2.2047,
        2.1472, 2.1042, 2.0562, 1.9899, 1.9779, 1.9507, 1.8686, 1.8029,
        1.7043, 1.5546, 1.3101, 1.1705, 0.9961, 0.8060, 0.6777, 0.4794,
        0.4017, 0.3650, 0.2905, 0.2197, 0.1679, 0.1313, 0.1073, 0.0510,
        0.0599, 0.0302,
    ]),
    'time_offset': -0.3,
}


# Legacy Chunc / machbusterNew plotting overlay (59-point digitized). Its
# (time=0, pressure=0) point is absent from the 346-point source, and its
# time-zero/baseline are undocumented; do not use it as protocol v0.1 data.
# Originally mis-labeled as "Hasegawa"
# in srm_1d/examples/ISP_Super_Loki.py via a copy-paste, then
# mis-promoted to ISP_SUPER_LOKI_EXPERIMENTAL in v0.7.3 Phase A.3
# (2026-05-24). User flagged the provenance error 2026-05-25: this is
# actually Chunc data (matches srm_1d/examples/machbusterNew.py's
# CHUNC_EXPERIMENTAL). We do NOT have a verified Super Loki
# experimental dataset in the repo.
CHUNC_EXPERIMENTAL = {
    'label': 'Experimental (Chunc / machbusterNew)',
    'time': np.array([
        0.0, 0.01, 0.044, 0.085, 0.126, 0.167, 0.207, 0.248, 0.289, 0.33,
        0.37, 0.411, 0.452, 0.493, 0.533, 0.574, 0.615, 0.655, 0.696,
        0.737, 0.778, 0.818, 0.859, 0.9, 0.941, 0.981, 1.022, 1.063,
        1.104, 1.144, 1.185, 1.226, 1.267, 1.307, 1.348, 1.389, 1.429,
        1.47, 1.511, 1.552, 1.592, 1.633, 1.674, 1.715, 1.755, 1.796,
        1.837, 1.878, 1.918, 1.959, 2.0, 2.041, 2.081, 2.122, 2.163,
        2.203, 2.244, 2.285, 2.326,
    ]),
    'pressure': np.array([
        0.0, 8.466, 8.735, 8.807, 8.83, 8.863, 8.815, 8.808, 8.801, 8.834,
        8.834, 8.883, 8.841, 8.808, 8.789, 8.76, 8.709, 8.67, 8.617, 8.57,
        8.61, 8.602, 8.592, 8.52, 8.412, 8.106, 7.254, 6.251, 5.644, 5.241,
        4.852, 4.273, 3.859, 3.498, 3.112, 2.699, 2.44, 2.199, 1.955, 1.69,
        1.492, 1.308, 1.149, 0.979, 0.802, 0.669, 0.553, 0.45, 0.355, 0.288,
        0.237, 0.193, 0.15, 0.125, 0.103, 0.081, 0.062, 0.049, 0.041,
    ]),
}


def load_experimental_csv(filepath, time_col=0, pressure_col=1,
                          delimiter=',', skip_header=1,
                          pressure_unit='MPa', label=None):
    """
    Load experimental data from a CSV file.

    Parameters
    ----------
    filepath : str
        Path to CSV file.
    time_col : int
        Column index for time data (0-based). Default: 0.
    pressure_col : int
        Column index for pressure data (0-based). Default: 1.
    delimiter : str
        Column delimiter. Default: ','.
    skip_header : int
        Number of header rows to skip. Default: 1.
    pressure_unit : str
        Unit of pressure in the file. One of 'Pa', 'kPa', 'MPa', 'psi'.
        Data will be converted to MPa for plotting.
    label : str or None
        Legend label. If None, uses the filename.

    Returns
    -------
    dict with 'time', 'pressure' (in MPa), and 'label'.
    """
    data = np.loadtxt(filepath, delimiter=delimiter, skiprows=skip_header)
    time = data[:, time_col]
    pressure = data[:, pressure_col]

    # Convert to MPa
    conversion = {
        'Pa': 1e-6,
        'kPa': 1e-3,
        'MPa': 1.0,
        'psi': 0.00689476,
    }
    if pressure_unit not in conversion:
        raise ValueError(
            f"Unknown pressure unit '{pressure_unit}'. "
            f"Use one of: {list(conversion.keys())}"
        )
    pressure = pressure * conversion[pressure_unit]

    if label is None:
        label = filepath.split('/')[-1].split('\\')[-1]

    return {
        'time': time,
        'pressure': pressure,
        'label': label,
    }


# ================================================================
# Pressure Plot
# ================================================================

def plot_pressure(result, title="Head-End Pressure",
                  experimental=None, time_offset=0.0,
                  save_path=None, ax=None, n_head_cells=1):
    """
    Plot head-end pressure vs time.

    Parameters
    ----------
    result : dict
        Output from run_simulation. Must have 'time' and 'P_head'.
    title : str
        Plot title.
    experimental : dict or list of dicts, optional
        Experimental data to overlay. Each dict must have 'time'
        (seconds), 'pressure' (MPa), and optionally 'label'.
    time_offset : float
        Seconds added to all experimental time axes. Use to align
        the experimental ignition event with the simulation t=0.
        For per-dataset offsets, pre-shift the 'time' arrays before
        passing.
    save_path : str or None
        If provided, save figure to this path.
    ax : matplotlib Axes or None
        If provided, plot on this axes. Otherwise create a new figure.
    n_head_cells : int
        Number of head-end cells to plot. 1 = just P[0] (default).
        Values > 1 show pressure at cells 0..n-1 to reveal axial
        variation during the ignition transient.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes (None, None if no matplotlib).
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    sc = as_channels(result)

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.get_figure()

    if n_head_cells <= 1:
        # Single head-end trace, drawn via the generic unit-aware path
        # (Pa -> MPa carried on the channel).
        plot_channels(sc, 'P_head', display_units={'P_head': 'MPa'},
                      overlay=True, labels={'P_head': '1D PISO (this work)'},
                      colors={'P_head': 'b'}, axes=ax)
    else:
        # Multi-cell traces from the bore-pressure axial channel
        P_axial = sc.axial.get('P')
        if P_axial is not None:
            snap_t = P_axial.times
            cell_colors = plt.cm.viridis(np.linspace(0, 0.7, n_head_cells))
            for ci in range(n_head_cells):
                P_cell = P_axial.getCell(ci, 'MPa')
                x_mm = P_axial.x_cells[ci] * 1000
                ax.plot(snap_t, P_cell, '-', color=cell_colors[ci],
                        linewidth=1.5, label=f'Cell {ci} (x={x_mm:.0f}mm)')
        # Always show P_head as the primary trace
        t = sc['time']
        P = sc.channels['P_head'].getData('MPa')
        ax.plot(t, P, 'b-', linewidth=2, alpha=0.4, label='P[0] (full res)')

    # Experimental overlay(s)
    if experimental is not None:
        if isinstance(experimental, dict):
            experimental = [experimental]
        colors = ['k', 'r', 'g', 'orange', 'purple']
        markers = ['o', 's', '^', 'D', 'v']
        for i, expt in enumerate(experimental):
            color = colors[i % len(colors)]
            marker = markers[i % len(markers)]
            label = expt.get('label', f'Experimental {i+1}')
            ax.plot(expt['time'] + time_offset, expt['pressure'],
                    color=color, linewidth=2, marker=marker,
                    markersize=3, markevery=max(1, len(expt['time'])//20),
                    label=label)

    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Head-End Pressure [MPa]', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=11, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)

    if own_fig:
        plt.tight_layout()
    if save_path:
        save_figure(fig, save_path)

    return fig, ax


# ================================================================
# Thrust Plot
# ================================================================

def plot_thrust(result, performance, title="Thrust",
                save_path=None, ax=None):
    """
    Plot thrust and Isp vs time.

    Parameters
    ----------
    result : dict
        Output from run_simulation.
    performance : dict
        Output from compute_motor_performance. Must have 'thrust' and 'Isp'.
    title : str
        Plot title.
    save_path : str or None
        If provided, save figure to this path.
    ax : matplotlib Axes or None
        If provided, plot on this axes (thrust only). Otherwise create
        a two-panel figure with thrust and Isp.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes.
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    # Wrap the performance series as channels so thrust/Isp flow through
    # the generic plot_channels path (unit-aware, same as every other
    # scalar trace).
    sc = as_channels(result)
    pc = SimulationChannels()
    pc.add('time', sc.channels['time'])
    pc.add('thrust', Channel('Thrust', 'N', np.asarray(performance['thrust'])))
    pc.add('Isp', Channel('Isp', 's', np.asarray(performance['Isp'])))

    own_fig = ax is None
    if own_fig:
        fig, axes = plot_channels(pc, ['thrust', 'Isp'],
                                  colors={'thrust': 'b', 'Isp': 'r'},
                                  title=title)
        ax_thrust = axes[0]
    else:
        fig, _ = plot_channels(pc, 'thrust', overlay=True, axes=ax,
                               colors={'thrust': 'b'}, title=title)
        ax_thrust = ax
        ax_thrust.set_ylabel('Thrust [N]', fontsize=12)
        ax_thrust.get_legend().remove()
        axes = [ax]

    # Add designation annotation
    desig = performance.get('motor_designation', '')
    impulse = performance.get('total_impulse', 0)
    if desig:
        ax_thrust.annotate(
            f"{desig}\n{impulse:.0f} N·s",
            xy=(0.98, 0.95), xycoords='axes fraction',
            ha='right', va='top', fontsize=12,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8),
        )

    if save_path:
        save_figure(fig, save_path)

    return fig, axes


# ================================================================
# Flow Snapshot Plot
# ================================================================

def plot_flow_snapshot(result, t_target=None, snap_index=None,
                      title=None, save_path=None):
    """
    Plot a flow field snapshot (P, Mach, burn rate, port diameter).

    Parameters
    ----------
    result : dict
        Output from run_simulation. Must have 'snapshots'.
    t_target : float or None
        Target simulation time [s]. Uses the nearest snapshot.
    snap_index : int or None
        Snapshot index directly. Overrides t_target.
    title : str or None
        Plot title. Auto-generated if None.
    save_path : str or None
        If provided, save figure to this path.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes (2x2 grid).
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    snapshots = _frames_from_channels(as_channels(result))
    if not snapshots:
        print("No snapshots in result dict.")
        return None, None

    # Select snapshot
    if snap_index is not None:
        snap_index = max(0, min(snap_index, len(snapshots) - 1))
    elif t_target is not None:
        snap_times = [s['t'] for s in snapshots]
        snap_index = min(range(len(snap_times)),
                         key=lambda i: abs(snap_times[i] - t_target))
    else:
        # Default: ~30% through the burn
        snap_index = len(snapshots) // 3

    snap = snapshots[snap_index]
    x_mm = snap['x'] * 1000

    if title is None:
        title = f"Flow Snapshot at t = {snap['t']:.3f} s"

    # v0.7.3 Phase A: 3x2 grid adds u_cell (velocity, with sign-band) and
    # T panels for aft_basket / head_basket diagnostics. Velocity sign
    # reversal under aft_basket topology shows up cleanly in the new
    # u_cell panel.
    fig, axes = plt.subplots(3, 2, figsize=(13, 11))

    # Pressure
    axes[0, 0].plot(x_mm, snap['P'] / 1e6, 'b-', linewidth=1.5)
    axes[0, 0].set_ylabel('Pressure [MPa]')
    axes[0, 0].set_title('Pressure')
    axes[0, 0].grid(True, alpha=0.3)

    # Mach
    axes[0, 1].plot(x_mm, snap['Mach'], 'r-', linewidth=1.5)
    axes[0, 1].set_ylabel('Mach Number')
    axes[0, 1].set_title('Mach Number')
    axes[0, 1].grid(True, alpha=0.3)

    # Velocity (cell-centered, signed). Positive = downstream, negative =
    # upstream (back-firing aft_basket diagnostic). Color the line by
    # sign-band so reversed-flow regions are visually obvious.
    if 'u' in snap:
        u = snap['u']
        axes[1, 0].axhline(0.0, color='gray', linewidth=0.8, alpha=0.6)
        axes[1, 0].plot(x_mm, u, 'k-', linewidth=1.5)
        # Fill positive region green, negative region red
        axes[1, 0].fill_between(x_mm, u, 0.0, where=(u > 0), color='tab:green',
                                alpha=0.25, interpolate=True, label='downstream')
        axes[1, 0].fill_between(x_mm, u, 0.0, where=(u < 0), color='tab:red',
                                alpha=0.25, interpolate=True, label='upstream')
        axes[1, 0].set_ylabel('u_cell [m/s]')
        axes[1, 0].set_title('Cell Velocity (sign-banded)')
        axes[1, 0].legend(fontsize=8, loc='best')
        axes[1, 0].grid(True, alpha=0.3)
    else:
        axes[1, 0].text(0.5, 0.5, 'u_cell not in snapshot',
                        transform=axes[1, 0].transAxes, ha='center')

    # Temperature
    if 'T' in snap:
        axes[1, 1].plot(x_mm, snap['T'], 'tab:orange', linewidth=1.5)
        axes[1, 1].set_ylabel('Gas Temperature [K]')
        axes[1, 1].set_title('Gas Temperature')
        axes[1, 1].grid(True, alpha=0.3)
    else:
        axes[1, 1].text(0.5, 0.5, 'T not in snapshot',
                        transform=axes[1, 1].transAxes, ha='center')

    # Burn rate
    axes[2, 0].plot(x_mm, snap['r_total'] * 1000, 'b-', linewidth=1.5,
                    label='Total')
    axes[2, 0].plot(x_mm, snap['r_erosive'] * 1000, 'r--', linewidth=1.5,
                    label='Erosive')
    r_normal = (snap['r_total'] - snap['r_erosive']) * 1000
    axes[2, 0].plot(x_mm, r_normal, 'g:', linewidth=1.5, label='Normal')
    axes[2, 0].set_xlabel('Position [mm]')
    axes[2, 0].set_ylabel('Burn Rate [mm/s]')
    axes[2, 0].set_title('Burn Rate')
    axes[2, 0].legend(fontsize=9)
    axes[2, 0].grid(True, alpha=0.3)

    # Show end-face mass injection locations if available
    if 'endface_msource' in snap:
        ef = snap['endface_msource']
        ef_mask = ef > 0
        if np.any(ef_mask):
            ax2 = axes[2, 0].twinx()
            ax2.bar(x_mm[ef_mask], ef[ef_mask], width=x_mm[1]-x_mm[0],
                    alpha=0.25, color='orange', label='End-face source')
            ax2.set_ylabel('End-face [kg/(m·s)]', fontsize=8, color='orange')
            ax2.tick_params(axis='y', labelcolor='orange', labelsize=7)

    # Port diameter
    axes[2, 1].plot(x_mm, snap['D_port'] * 1000, 'k-', linewidth=1.5)
    axes[2, 1].set_xlabel('Position [mm]')
    axes[2, 1].set_ylabel('Port Diameter [mm]')
    axes[2, 1].set_title('Port Diameter')
    axes[2, 1].grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=14)
    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path)

    return fig, axes


# ================================================================
# Multi-Snapshot Subplot Grid (v0.7.3 Phase A)
# ================================================================

_FIELD_LABELS = {
    'P':         ('Pressure [MPa]',     lambda v: v / 1e6),
    'Mach':      ('Mach',               lambda v: v),
    'u':    ('u_cell [m/s]',       lambda v: v),
    'T':         ('Gas Temp [K]',       lambda v: v),
    'T_surf':    ('Surface Temp [K]',   lambda v: v),
    'r_total':   ('r_total [mm/s]',     lambda v: v * 1000),
    'r_erosive': ('r_erosive [mm/s]',   lambda v: v * 1000),
    'D_port':    ('Port D [mm]',        lambda v: v * 1000),
    'mass_source':    ('mass_source [kg/m/s]',    lambda v: v),
    'thermal_source': ('thermal_source [W/m]',    lambda v: v),
    'is_burning':     ('is_burning',              lambda v: v),
    'pyrogen_surface_heat_flux': ('pyro_heat_flux [W/m²]', lambda v: v),
    'radiation_heat_flux':       ('rad_heat_flux [W/m²]', lambda v: v),
}


def plot_flow_snapshots(result, t_targets, fields=('P', 'Mach', 'u', 'T'),
                        title=None, save_path=None):
    """v0.7.3 Phase A — multi-time snapshot subplot grid.

    Renders one row per ``t_target`` (using nearest captured snapshot)
    and one column per field. Each cell shows the field's spatial
    profile (x vs value) at that snapshot time. Useful for comparing
    flow evolution across ignition, plateau, and burnout phases in a
    single figure.

    Parameters
    ----------
    result : dict
        Output from ``run_simulation``. Must have ``'snapshots'``.
    t_targets : iterable of float
        Target times [s]; rendered in the supplied order.
    fields : iterable of str
        Snapshot keys to plot per row. Defaults to
        ``('P', 'Mach', 'u', 'T')``. Any of the keys in
        ``_FIELD_LABELS`` may be requested; unknown keys are skipped
        with a one-line note in the panel.
    title : str or None
        Overall figure title; auto-generated if None.
    save_path : str or None
        If provided, save figure via ``save_figure``; otherwise the
        caller is expected to ``plt.show()`` or further mutate the
        figure.

    Returns
    -------
    fig, axes : matplotlib Figure and 2D axes array
        ``axes`` is always 2D ``(len(t_targets), len(fields))`` for
        consistent indexing even when there's only one row/column.

    Velocity sign-banding (``u_cell``) is applied automatically per
    panel for at-a-glance reverse-flow diagnosis under aft_basket
    topology.
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    snapshots = _frames_from_channels(as_channels(result))
    if not snapshots:
        print("No snapshots in result dict.")
        return None, None

    t_targets = list(t_targets)
    fields = list(fields)
    n_rows = len(t_targets)
    n_cols = len(fields)
    if n_rows == 0 or n_cols == 0:
        return None, None

    snap_times = np.array([s['t'] for s in snapshots])

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.2 * n_cols + 0.6, 2.4 * n_rows + 0.8),
        squeeze=False,
    )

    for r, t_tgt in enumerate(t_targets):
        idx = int(np.argmin(np.abs(snap_times - t_tgt)))
        snap = snapshots[idx]
        x_mm = snap['x'] * 1000
        actual_t = snap['t']
        for c, field in enumerate(fields):
            ax = axes[r, c]
            if field not in _FIELD_LABELS:
                ax.text(0.5, 0.5, f"unknown field\n'{field}'",
                        ha='center', va='center', transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            if field not in snap:
                ax.text(0.5, 0.5, f"'{field}' not in snapshot",
                        ha='center', va='center', transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            ylabel, scale = _FIELD_LABELS[field]
            y = scale(np.asarray(snap[field]))
            if field == 'u':
                ax.axhline(0.0, color='gray', linewidth=0.7, alpha=0.6)
                ax.plot(x_mm, y, 'k-', linewidth=1.2)
                ax.fill_between(x_mm, y, 0.0, where=(y > 0),
                                color='tab:green', alpha=0.25,
                                interpolate=True)
                ax.fill_between(x_mm, y, 0.0, where=(y < 0),
                                color='tab:red', alpha=0.25,
                                interpolate=True)
            elif field == 'is_burning':
                ax.fill_between(x_mm, 0, y, step='mid',
                                color='tab:orange', alpha=0.5)
                ax.set_ylim(-0.05, 1.05)
            else:
                ax.plot(x_mm, y, linewidth=1.2)
            if c == 0:
                ax.set_ylabel(f"t={actual_t:.3f}s\n{ylabel}", fontsize=9)
            else:
                ax.set_ylabel(ylabel, fontsize=9)
            if r == 0:
                ax.set_title(field, fontsize=10)
            if r == n_rows - 1:
                ax.set_xlabel('x [mm]', fontsize=9)
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='both', labelsize=8)

    if title is None:
        title = (f"Flow Snapshots — "
                 f"{n_rows} time slice{'s' if n_rows > 1 else ''}")
    fig.suptitle(title, fontsize=13)
    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path)

    return fig, axes


# ================================================================
# x-t Field Heatmap (v0.7.3 Phase A)
# ================================================================

def plot_field_heatmap(result, fields=('P', 'u', 'T', 'is_burning'),
                       cmap=None, title=None, save_path=None,
                       t_max=None):
    """v0.7.3 Phase A — render snapshot fields as 2D x-t heatmaps.

    For each requested field, builds a pcolormesh with x on the
    horizontal axis and snapshot-time on the vertical axis (origin
    at lower-left, time increasing upward). Color encodes field
    value. Useful for visualizing back→front ignition cascades
    (``is_burning`` shows diagonal stripes), reverse-flow regions
    (``u_cell`` shows a sign-band), and pressure waves (``P``).

    Parameters
    ----------
    result : dict
        Output from ``run_simulation``. Must have ``'snapshots'``.
    fields : iterable of str
        Snapshot keys to render; one panel per field. Defaults to
        ``('P', 'u', 'T', 'is_burning')``.
    cmap : dict or None
        Optional per-field colormap override mapping field name to
        a matplotlib colormap name. Defaults:
        ``{'P': 'viridis', 'u': 'RdBu_r', 'T': 'inferno',
        'is_burning': 'Oranges'}``.
    title : str or None
        Figure title; auto-generated if None.
    save_path : str or None
        If provided, save figure via ``save_figure``.
    t_max : float or None
        If provided, clip the time axis to ``t <= t_max`` for a
        zoomed-in view of the ignition transient (set to ~0.5s for
        spike-shape diagnostics).

    Returns
    -------
    fig, axes : matplotlib Figure and axes (1D array len(fields))
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    snapshots = _frames_from_channels(as_channels(result))
    if not snapshots:
        print("No snapshots in result dict.")
        return None, None

    fields = list(fields)
    n_cols = len(fields)
    if n_cols == 0:
        return None, None

    snap_times = np.array([s['t'] for s in snapshots])
    if t_max is not None:
        keep_mask = snap_times <= float(t_max)
        if not np.any(keep_mask):
            print(f"No snapshots at t<={t_max}s.")
            return None, None
        snap_indices = np.where(keep_mask)[0]
    else:
        snap_indices = np.arange(len(snapshots))

    x = snapshots[0]['x']
    x_mm = x * 1000
    t_axis = snap_times[snap_indices]

    default_cmaps = {
        'P': 'viridis', 'u': 'RdBu_r', 'T': 'inferno',
        'is_burning': 'Oranges', 'Mach': 'plasma',
        'r_total': 'viridis', 'mass_source': 'cividis',
        'thermal_source': 'magma',
    }
    if cmap is None:
        cmap = {}
    cmap_for = {**default_cmaps, **dict(cmap)}

    fig, axes = plt.subplots(
        1, n_cols, figsize=(3.6 * n_cols + 0.8, 4.5), squeeze=False,
    )
    axes = axes[0]

    for c, field in enumerate(fields):
        ax = axes[c]
        if field not in _FIELD_LABELS or field not in snapshots[0]:
            ax.text(0.5, 0.5, f"'{field}'\nunavailable",
                    ha='center', va='center', transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        ylabel, scale = _FIELD_LABELS[field]
        # Build (n_t, n_x) matrix
        Z = np.array([scale(np.asarray(snapshots[i][field]))
                      for i in snap_indices])
        cmap_name = cmap_for.get(field, 'viridis')
        # Symmetric colorbar for signed fields
        if field == 'u':
            zmax = float(np.max(np.abs(Z))) if Z.size else 1.0
            zmax = zmax if zmax > 0 else 1.0
            im = ax.pcolormesh(x_mm, t_axis, Z, cmap=cmap_name,
                               vmin=-zmax, vmax=zmax, shading='auto')
        else:
            im = ax.pcolormesh(x_mm, t_axis, Z, cmap=cmap_name,
                               shading='auto')
        ax.set_title(f"{field} — {ylabel}", fontsize=10)
        ax.set_xlabel('x [mm]', fontsize=9)
        if c == 0:
            ax.set_ylabel('t [s]', fontsize=9)
        ax.tick_params(axis='both', labelsize=8)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=8)

    if title is None:
        title = "Field Heatmap (x vs t)"
    fig.suptitle(title, fontsize=13)
    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path)

    return fig, axes


# ================================================================
# Combined Summary Plot
# ================================================================

def plot_summary(result, performance=None, experimental=None,
                 time_offset=0.0,
                 title="Simulation Summary", save_path=None):
    """
    Combined plot: pressure, thrust (if available), and a flow snapshot.

    Parameters
    ----------
    result : dict
        Output from run_simulation.
    performance : dict or None
        Output from compute_motor_performance. If None, thrust panel
        is replaced with exit pressure.
    experimental : dict or list of dicts, optional
        Experimental pressure data for overlay.
    time_offset : float
        Seconds added to all experimental time axes (alignment knob).
    title : str
        Overall figure title.
    save_path : str or None
        If provided, save figure to this path.

    Returns
    -------
    fig, axes : matplotlib Figure and Axes.
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — cannot plot.")
        return None, None

    sc = as_channels(result)
    t = sc['time']

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (0,0) Pressure trace
    plot_pressure(sc, title="Head-End Pressure",
                  experimental=experimental, time_offset=time_offset,
                  ax=axes[0, 0])

    # (0,1) Thrust or exit pressure
    if performance is not None:
        axes[0, 1].plot(t, performance['thrust'], 'b-', linewidth=2)
        axes[0, 1].set_ylabel('Thrust [N]', fontsize=11)
        axes[0, 1].set_title('Thrust', fontsize=12)
        desig = performance.get('motor_designation', '')
        if desig:
            axes[0, 1].annotate(
                desig, xy=(0.98, 0.95), xycoords='axes fraction',
                ha='right', va='top', fontsize=12,
                bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='wheat', alpha=0.8),
            )
    else:
        axes[0, 1].plot(t, sc.channels['P_exit'].getData('MPa'), 'r-', linewidth=2)
        axes[0, 1].set_ylabel('Exit Pressure [MPa]', fontsize=11)
        axes[0, 1].set_title('Nozzle-End Pressure', fontsize=12)
    axes[0, 1].set_xlabel('Time [s]', fontsize=11)
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_ylim(bottom=0)
    axes[0, 1].set_xlim(left=0)

    # (1,0) and (1,1) Flow snapshot at ~30% burn
    snapshots = _frames_from_channels(sc)
    if snapshots:
        snap = snapshots[len(snapshots) // 3]
        x_mm = snap['x'] * 1000

        axes[1, 0].plot(x_mm, snap['r_total'] * 1000, 'b-', lw=1.5,
                        label='Total')
        axes[1, 0].plot(x_mm, snap['r_erosive'] * 1000, 'r--', lw=1.5,
                        label='Erosive')
        axes[1, 0].set_xlabel('Position [mm]', fontsize=11)
        axes[1, 0].set_ylabel('Burn Rate [mm/s]', fontsize=11)
        axes[1, 0].set_title(f'Burn Rate at t={snap["t"]:.2f}s', fontsize=12)
        axes[1, 0].legend(fontsize=9)
        axes[1, 0].grid(True, alpha=0.3)

        axes[1, 1].plot(x_mm, snap['Mach'], 'r-', lw=1.5)
        axes[1, 1].set_xlabel('Position [mm]', fontsize=11)
        axes[1, 1].set_ylabel('Mach Number', fontsize=11)
        axes[1, 1].set_title(f'Mach at t={snap["t"]:.2f}s', fontsize=12)
        axes[1, 1].grid(True, alpha=0.3)
    else:
        for a in [axes[1, 0], axes[1, 1]]:
            a.text(0.5, 0.5, 'No snapshots', ha='center', va='center',
                   transform=a.transAxes, fontsize=14, color='gray')

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path)

    return fig, axes


# ================================================================
# Comparison Plot (1D PISO vs openMotor or experimental)
# ================================================================

def plot_comparison(result, perf=None, reference=None,
                    ref_label="openMotor", title="1D PISO vs Reference",
                    save_path=None):
    """
    Plot pressure and thrust comparison between our sim and a reference.

    Parameters
    ----------
    result : dict
        Output from run_simulation.
    perf : dict or None
        Output from compute_motor_performance.
    reference : dict or None
        Reference data with at least 'time' and 'pressure' arrays.
        Can be from load_openmotor_csv() or load_experimental_csv().
        If it has 'force', thrust is also compared.
        Pressure should be in Pa (or if < 100, assumed MPa).
    ref_label : str
        Legend label for the reference data.
    title : str
        Plot title.
    save_path : str or None
        If provided, save figure.

    Returns
    -------
    fig, axes
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available.")
        return None, None

    has_thrust = perf is not None and reference is not None and 'force' in reference
    n_rows = 2 if has_thrust else 1
    fig, axes_raw = plt.subplots(n_rows, 1, figsize=(10, 5 * n_rows), sharex=True)
    if n_rows == 1:
        axes = [axes_raw]
    else:
        axes = list(axes_raw)

    sc = as_channels(result)
    t_sim = sc['time']
    P_sim = sc.channels['P_head'].getData('MPa')

    # Pressure panel
    axes[0].plot(t_sim, P_sim, 'b-', linewidth=2, label='1D PISO')

    if reference is not None:
        t_ref = reference['time']
        P_ref = reference['pressure']
        if np.max(P_ref) > 100:  # Assume Pa, convert to MPa
            P_ref = P_ref / 1e6
        axes[0].plot(t_ref, P_ref, 'r--', linewidth=2, label=ref_label)

    axes[0].set_ylabel('Pressure [MPa]', fontsize=12)
    axes[0].set_title(title, fontsize=14)
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(bottom=0)

    # Thrust panel
    if has_thrust:
        axes[1].plot(t_sim, perf['thrust'], 'b-', linewidth=2, label='1D PISO')
        t_ref = reference['time']
        F_ref = reference['force']
        axes[1].plot(t_ref, F_ref, 'r--', linewidth=2, label=ref_label)
        axes[1].set_ylabel('Thrust [N]', fontsize=12)
        axes[1].legend(fontsize=11)
        axes[1].grid(True, alpha=0.3)
        axes[1].set_ylim(bottom=0)

    axes[-1].set_xlabel('Time [s]', fontsize=12)
    axes[-1].set_xlim(left=0)
    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path)

    return fig, axes


# ================================================================
# Per-Grain Regression Plot
# ================================================================

def plot_grain_regression(grain_metrics, geo, title="Per-Grain Regression",
                          save_path=None):
    """
    Plot per-grain regression depth and web remaining vs time.

    Parameters
    ----------
    grain_metrics : dict
        Output from compute_grain_metrics().
    geo : MotorGeometry
        The geometry (for segment labels).
    title : str
        Plot title.
    save_path : str or None
        If provided, save figure.

    Returns
    -------
    fig, axes
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available.")
        return None, None

    t = grain_metrics['snap_times']
    reg = grain_metrics['regression'] * 1000  # to mm
    web = grain_metrics['web'] * 1000
    n_grains = reg.shape[1]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    colors = plt.cm.tab10(np.linspace(0, 1, max(n_grains, 1)))

    for k in range(n_grains):
        axes[0].plot(t, reg[:, k], color=colors[k], linewidth=1.5,
                     label=f'Grain {k}')
        axes[1].plot(t, web[:, k], color=colors[k], linewidth=1.5,
                     label=f'Grain {k}')

    axes[0].set_xlabel('Time [s]', fontsize=11)
    axes[0].set_ylabel('Regression [mm]', fontsize=11)
    axes[0].set_title('Radial Regression', fontsize=12)
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel('Time [s]', fontsize=11)
    axes[1].set_ylabel('Web Remaining [mm]', fontsize=11)
    axes[1].set_title('Web Remaining', fontsize=12)
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)
    axes[1].set_ylim(bottom=0)

    fig.suptitle(title, fontsize=14)
    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path)

    return fig, axes
