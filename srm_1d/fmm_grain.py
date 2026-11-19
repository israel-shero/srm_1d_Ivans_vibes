"""
fmm_grain.py — Fast Marching Method grain support (openMotor bridge)
=====================================================================

Bridges srm_1d to openMotor's FmmGrain machinery (skfmm-based regression
maps) for arbitrary cross-sections: Finocyl, Star, Moonburner, X, D, C,
RodTube, Custom.

WHAT THIS MODULE DOES
    1. Lazily adds the local openMotor checkout to sys.path so
       `import motorlib` and `import mathlib` work.
    2. Runs openMotor's FmmGrain pipeline (initGeometry → generateCoreMap →
       generateRegressionMap) on a populated grain instance.
    3. Samples the resulting regression map at fine resolution and packs
       it into an `FmmTable` of (regression_depth, perimeter, port_area)
       arrays suitable for srm_1d's @njit hot loop.

WHAT THIS MODULE DOES NOT
    - Reimplement skfmm or openMotor's FmmGrain — we call upstream.
    - Run openMotor's solver. We use the cross-section processor only.
    - Modify the time loop or geometry compilation directly. That's
      Phase 2, in grain_geometry.py.

DESIGN NOTE
    FMM tables are per-grain (one set per FmmTable). Per-cell regression
    in the hot loop is preserved — the Ma erosive model produces axially
    varying burn rates, and per-cell regression captures that variation
    while the heavy data (the regressionMap and sampled tables) stays
    shared at the grain level.
"""

import sys
import types
from dataclasses import dataclass
from pathlib import Path
import numpy as np

try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def wrapper(func):
            return func
        return wrapper


# Cached so we only do path manipulation once per process
_OPENMOTOR_PATH_CONFIGURED = False


# ================================================================
# Pure-Python (Numba-JIT) replacement for openMotor's Cython
# `mathlib._find_perimeter_cy._get_perimeter`. Avoids requiring the
# user to compile a Cython extension (MSVC build tools on Windows,
# etc.). The algorithm is a verbatim port of openMotor's marching-
# squares perimeter computation; only the language differs.
# ================================================================

@njit(cache=True)
def _marching_squares_perimeter(arr, level):
    """
    Marching-squares perimeter of the iso-contour at `level` in a 2D
    float64 array. Mirrors openMotor's
    `mathlib._find_perimeter_cy._get_perimeter`, including the
    "skip squares outside the motor radius" optimization.

    Returns the perimeter in array-index distance units (one unit =
    one pixel-width); the caller (openMotor's getCorePerimeter)
    multiplies by `mapToLength` to get meters.
    """
    perimeter = 0.0
    rows, cols = arr.shape[0], arr.shape[1]
    half_r = rows / 2.0
    half_c = cols / 2.0
    radius_thresh = half_r - 3.0

    for r0 in range(3, rows - 4):
        for c0 in range(3, cols - 4):
            r1 = r0 + 1
            c1 = c0 + 1

            ul = arr[r0, c0]
            ur = arr[r0, c1]
            ll = arr[r1, c0]
            lr = arr[r1, c1]

            sc = 0
            if ul > level:
                sc += 1
            if ur > level:
                sc += 2
            if ll > level:
                sc += 4
            if lr > level:
                sc += 8

            if sc == 0 or sc == 15:
                continue

            # Skip squares outside the motor radius (cylindrical mask).
            dist_r = r0 + 0.5 - half_r
            dist_c = c0 + 0.5 - half_c
            if radius_thresh < (dist_r * dist_r + dist_c * dist_c) ** 0.5:
                continue

            # _get_fraction (linear interp): zero-safe when both ends equal
            top = 0.0 if ul == ur else (level - ul) / (ur - ul)
            bottom = 0.0 if ll == lr else (level - ll) / (lr - ll)
            left = 0.0 if ll == ul else (level - ll) / (ul - ll)
            right = 0.0 if lr == ur else (level - lr) / (ur - lr)

            add = 0.0
            if sc == 1 or sc == 14:
                add = (top * top + (1.0 - left) ** 2) ** 0.5
            elif sc == 2 or sc == 13:
                add = ((1.0 - top) ** 2 + (1.0 - right) ** 2) ** 0.5
            elif sc == 3 or sc == 12:
                add = ((right - left) ** 2 + 1.0) ** 0.5
            elif sc == 4 or sc == 11:
                add = (left * left + bottom * bottom) ** 0.5
            elif sc == 5 or sc == 10:
                add = ((top - bottom) ** 2 + 1.0) ** 0.5
            elif sc == 6:
                add = (((1.0 - top) ** 2 + (1.0 - right) ** 2) ** 0.5
                       + (left * left + bottom * bottom) ** 0.5)
            elif sc == 7 or sc == 8:
                add = ((1.0 - bottom) ** 2 + right * right) ** 0.5
            elif sc == 9:
                add = ((top * top + (1.0 - left) ** 2) ** 0.5
                       + ((1.0 - bottom) ** 2 + right * right) ** 0.5)
            perimeter += add

    return perimeter


def _get_perimeter_replacement(image, level, vertex_connect_high, returning_contours):
    """
    Python-level wrapper matching the Cython
    `_get_perimeter(image, level, vertex_connect_high, returning_contours)`
    signature. Returns (perimeter, segments) where segments is always
    an empty list — openMotor's `getCorePerimeter` only uses [0]
    (the perimeter), and we don't need contour reassembly for the
    solver path. If `returning_contours` is requested (e.g. for
    plotting via openMotor's `getRegressionData`), this raises.
    """
    if returning_contours:
        raise NotImplementedError(
            "srm_1d's pure-Python find_perimeter shim doesn't support "
            "returning_contours=True. Build openMotor's Cython extension "
            "(`python setup.py build_ext --inplace` in the openMotor dir) "
            "if you need contour visualizations."
        )
    arr = np.ascontiguousarray(np.asarray(image, dtype=np.float64))
    return float(_marching_squares_perimeter(arr, float(level))), []


def _install_fake_find_perimeter_cy():
    """
    Inject a fake `mathlib._find_perimeter_cy` module exposing our
    Numba `_get_perimeter` so openMotor's mathlib package can import
    without requiring the Cython build. Idempotent.
    """
    mod_name = 'mathlib._find_perimeter_cy'
    if mod_name in sys.modules:
        return
    fake = types.ModuleType(mod_name)
    fake._get_perimeter = _get_perimeter_replacement
    sys.modules[mod_name] = fake


def _find_openmotor_root():
    """
    Walk upward from this file looking for an openMotor checkout. Returns the
    checkout path (which contains `motorlib/` and `mathlib/`), or None if not
    found.

    Resolution order:
        1. `SRM1D_OPENMOTOR_PATH` environment variable (if set, must point
           directly at the inner openMotor/openMotor/ directory)
        2. Walk upward from this file, checking sibling checkout names used by
           supported layouts, including `openMotor`, `openMotor-compatible`,
           and nested `openMotor/openMotor` layouts.

    Returning None lets the caller raise a clear error message including
    the directories that were checked.
    """
    import os

    env_path = os.environ.get('SRM1D_OPENMOTOR_PATH')
    if env_path:
        candidate = Path(env_path).resolve()
        if (candidate / 'motorlib').is_dir():
            return candidate
        # Fall through to upward search if env path is wrong, with a hint
        # in the eventual error message.

    here = Path(__file__).resolve().parent
    for ancestor in [here, *here.parents]:
        candidates = (
            ancestor / 'openMotor' / 'openMotor',
            ancestor / 'openMotor-compatible',
            ancestor / 'openMotor',
        )
        for candidate in candidates:
            if (candidate / 'motorlib').is_dir():
                return candidate
        candidate = ancestor
        if (candidate / 'motorlib').is_dir():
            return candidate
    return None


def _setup_openmotor_path():
    """
    Add the local openMotor checkout to sys.path so `motorlib` and
    `mathlib` are importable. Inject the find_perimeter shim before
    `mathlib` imports its missing Cython dependency. Idempotent.

    Locates openMotor by walking upward from this file (see
    `_find_openmotor_root`). Set `SRM1D_OPENMOTOR_PATH` env var to
    override.
    """
    global _OPENMOTOR_PATH_CONFIGURED
    if _OPENMOTOR_PATH_CONFIGURED:
        return

    om_root = _find_openmotor_root()
    if om_root is None:
        here = Path(__file__).resolve().parent
        searched = [str(here)] + [str(p) for p in here.parents]
        raise ImportError(
            "openMotor checkout not found. Searched for "
            "`openMotor/openMotor/motorlib/` under each of:\n  "
            + "\n  ".join(searched)
            + "\nClone it adjacent to this project:\n"
            "    git clone https://github.com/reilleya/openMotor\n"
            "Or set SRM1D_OPENMOTOR_PATH to the inner openMotor/openMotor/ dir."
        )

    om_str = str(om_root)
    if om_str not in sys.path:
        sys.path.insert(0, om_str)

    # Inject the Cython-replacement BEFORE mathlib tries to import it.
    _install_fake_find_perimeter_cy()

    # Verify both required packages import.
    try:
        import motorlib  # noqa: F401
        import mathlib   # noqa: F401
    except ImportError as e:
        raise ImportError(
            f"Found openMotor at {om_root}, but `import motorlib`/`import mathlib` "
            f"failed: {e}\n"
            "Make sure scikit-fmm is installed: pip install scikit-fmm"
        ) from e

    _OPENMOTOR_PATH_CONFIGURED = True


# ================================================================
# FmmTable: tabulated grain shape data ready for the hot loop
# ================================================================

@dataclass
class FmmTable:
    """
    Pre-tabulated FMM grain shape data — perimeter and port area as
    1D functions of regression depth. Built once at simulation start
    by sampling an openMotor FmmGrain's regression map.

    Attributes
    ----------
    reg_depth : ndarray (T,)
        Regression depths sampled, [m]. Uniform grid from 0 to wall_web.
    perimeter : ndarray (T,)
        Burning core perimeter at each regression depth, [m].
        Used as `C_burn` in the srm_1d solver.
    port_area : ndarray (T,)
        Open cross-sectional area at each regression depth, [m²].
        Used as `A_port`. Includes the casting-tube fully-open value
        once regression > wall_web (after burnout).
    wall_web : float
        Maximum regression depth for this grain, [m]. Beyond this
        the grain is fully consumed (port_area = casting-tube area).
    grain_outer_diameter : float
        Casting tube ID, [m]. (openMotor: `props['diameter']`.)
    grain_length : float
        Axial length, [m]. (openMotor: `props['length']`.)
    inhibited_fwd : bool
    inhibited_aft : bool
        End-face inhibition flags. (Translated from openMotor's
        `inhibitedEnds` enum: Top→fwd, Bottom→aft.)
    geom_name : str
        openMotor `geomName` of the source grain (e.g. 'Finocyl').
        Informational only.
    """
    reg_depth: np.ndarray
    perimeter: np.ndarray
    port_area: np.ndarray
    wall_web: float
    grain_outer_diameter: float
    grain_length: float
    inhibited_fwd: bool
    inhibited_aft: bool
    geom_name: str

    @property
    def n_samples(self) -> int:
        return len(self.reg_depth)

    @property
    def initial_perimeter(self) -> float:
        """C_burn at regress=0."""
        return float(self.perimeter[0])

    @property
    def initial_port_area(self) -> float:
        """A_port at regress=0."""
        return float(self.port_area[0])


# Maps openMotor's 'inhibitedEnds' enum to (inhibit_fwd, inhibit_aft)
# (Top = forward / head end, Bottom = aft / nozzle end.)
_INHIBIT_MAP = {
    'Neither': (False, False),
    'Top':     (True,  False),
    'Bottom':  (False, True),
    'Both':    (True,  True),
}


def from_openmotor(om_grain, map_dim: int = 1001) -> FmmTable:
    """
    Run openMotor's FmmGrain pipeline on the supplied grain instance and
    extract its perimeter / port-area tables as numpy arrays.

    Parameters
    ----------
    om_grain : openMotor FmmGrain (e.g. Finocyl, Star, Moonburner)
        Properties (diameter, length, inhibitedEnds, and grain-specific
        ones like coreDiameter/numFins/finWidth/...) must already be set.
    map_dim : int
        FMM regression-map resolution. openMotor default is 1001. Higher
        is more accurate but quadratically slower (skfmm.distance is
        O(mapDim²)). Sim-level config knob.

    Returns
    -------
    FmmTable

    Notes
    -----
    The regression map is sampled at one entry per FMM pixel-width in
    regression depth (matches openMotor's internal sampling density),
    plus a final sample slightly beyond burnout to make the table cover
    [0, wall_web] inclusively for safe interpolation at the boundary.
    """
    _setup_openmotor_path()

    # Late import — only available after _setup_openmotor_path()
    from motorlib.grain import FmmGrain
    if not isinstance(om_grain, FmmGrain):
        raise TypeError(
            f"Expected an openMotor FmmGrain subclass, got {type(om_grain).__name__}"
        )

    # The full simulationSetup wants a config object; do its three steps directly.
    om_grain.initGeometry(map_dim)
    om_grain.generateCoreMap()
    om_grain.generateRegressionMap()

    wall_web = float(om_grain.wallWeb)
    if wall_web <= 0.0:
        raise ValueError(
            f"FmmGrain {type(om_grain).__name__} returned wallWeb={wall_web}; "
            "check grain properties (likely zero core or zero web)."
        )

    # Sample at ~one point per FMM pixel-width of regression. openMotor's
    # internal sampling uses int(maxDist*mapDim)+2 entries; we mirror.
    n_samples = max(int(wall_web * map_dim / om_grain.props['diameter'].getValue()) + 2,
                    64)
    reg_depth = np.linspace(0.0, wall_web, n_samples)

    # Sample perimeter and port area directly from openMotor's API.
    perimeter = np.empty(n_samples)
    port_area = np.empty(n_samples)
    casting_area = np.pi / 4.0 * om_grain.props['diameter'].getValue() ** 2
    for i, r in enumerate(reg_depth):
        if r >= wall_web:
            # At/past burnout, perimeter and face_area both go to 0; the
            # cell is fully consumed and the port equals the casting tube.
            perimeter[i] = 0.0
            port_area[i] = casting_area
        else:
            perimeter[i] = float(om_grain.getCorePerimeter(r))
            port_area[i] = float(om_grain.getPortArea(r))

    inh_str = om_grain.props['inhibitedEnds'].getValue()
    if inh_str not in _INHIBIT_MAP:
        raise ValueError(
            f"Unknown inhibitedEnds value {inh_str!r}. "
            f"Expected one of {list(_INHIBIT_MAP)}."
        )
    inh_fwd, inh_aft = _INHIBIT_MAP[inh_str]

    return FmmTable(
        reg_depth=reg_depth,
        perimeter=perimeter,
        port_area=port_area,
        wall_web=wall_web,
        grain_outer_diameter=float(om_grain.props['diameter'].getValue()),
        grain_length=float(om_grain.props['length'].getValue()),
        inhibited_fwd=inh_fwd,
        inhibited_aft=inh_aft,
        geom_name=getattr(type(om_grain), 'geomName', None) or type(om_grain).__name__,
    )


def from_ric_grain(ric_grain: dict, map_dim: int = 1001) -> FmmTable:
    """
    Convenience wrapper: take a .ric grain dict (with 'type' and
    'properties' keys) and return an FmmTable. Dispatches on the grain
    type to instantiate the right openMotor class.
    """
    _setup_openmotor_path()

    gtype = ric_grain['type']
    om_grain = _instantiate_openmotor_grain(gtype)
    om_grain.setProperties(ric_grain['properties'])
    return from_openmotor(om_grain, map_dim=map_dim)


# Lazy registry of openMotor FmmGrain subclasses by geomName.
_OM_GRAIN_CLASSES = None


def _instantiate_openmotor_grain(geom_name: str):
    """Look up an openMotor grain class by its `geomName` and return a
    fresh instance. Raises ValueError for unknown / non-FMM types."""
    global _OM_GRAIN_CLASSES
    if _OM_GRAIN_CLASSES is None:
        _setup_openmotor_path()
        from motorlib.grain import FmmGrain
        from motorlib.grains.finocyl import Finocyl
        from motorlib.grains.star import StarGrain
        from motorlib.grains.moonBurner import MoonBurner
        from motorlib.grains.cGrain import CGrain
        from motorlib.grains.dGrain import DGrain
        from motorlib.grains.xCore import XCore
        from motorlib.grains.custom import CustomGrain

        # NOTE: RodTubeGrain is a PerforatedGrain (analytic), not FmmGrain;
        # it can be added later by extending from_openmotor to accept
        # PerforatedGrain and skipping the FMM setup steps for analytic types.
        registry = {}
        for cls in (Finocyl, StarGrain, MoonBurner, CGrain, DGrain,
                    XCore, CustomGrain):
            if issubclass(cls, FmmGrain):
                registry[cls.geomName] = cls
        _OM_GRAIN_CLASSES = registry

    if geom_name not in _OM_GRAIN_CLASSES:
        raise ValueError(
            f"Grain type {geom_name!r} is not a supported FMM grain. "
            f"Supported FMM types: {sorted(_OM_GRAIN_CLASSES)}"
        )
    return _OM_GRAIN_CLASSES[geom_name]()


# ================================================================
# Parametric axial tapers
# ================================================================
#
# A tapered grain is a single grain whose cross-section varies along its
# axis (e.g. a finocyl whose fins grow from 0.25" to 0.5" tip-to-root).
# Rather than hand-author many short stepped segments, we define the
# taper by a start and end cross-section (or an arbitrary list of control
# stations) and let srm_1d build a smooth axial stack of *real* FMM
# regression tables — one per axial station — interpolated between them.
#
# The taper is authored as an UNRESOLVED `TaperSpec` (no FMM solves yet).
# The actual number of stations is decided AFTER the geometry snapper
# fixes the segment's cell count (see grain_geometry.build_snapped_geometry),
# so station density tracks the mesh rather than a fixed integer. `map_dim`
# is the cross-sectional (radial) FMM resolution and is independent of the
# axial mesh.

# Properties that cannot be linearly interpolated (integer-valued count).
# Booleans / strings (invertedFins, inhibitedEnds) are caught generically
# below by the "non-numeric" branch.
_NON_INTERPOLABLE_KEYS = {'numFins'}


def _interpolate_props(props_a: dict, props_b: dict, t: float) -> dict:
    """
    Linearly interpolate between two openMotor grain property dicts.

    Float-valued dimensions (coreDiameter, finLength, finWidth, diameter,
    length, ...) blend as `a + (b - a)*t`. Integer counts (`numFins`),
    booleans (`invertedFins`) and strings (`inhibitedEnds`) cannot be
    interpolated — they MUST be equal in both endpoints, else a
    `ValueError` is raised. Pass dimensions as floats (e.g. 0.0, not 0).
    """
    if set(props_a) != set(props_b):
        raise ValueError(
            "taper endpoints must define the same property keys; got "
            f"{sorted(props_a)} vs {sorted(props_b)}"
        )
    out = {}
    for k, va in props_a.items():
        vb = props_b[k]
        numeric = (
            isinstance(va, (int, float)) and not isinstance(va, bool)
            and isinstance(vb, (int, float)) and not isinstance(vb, bool)
        )
        if numeric and k not in _NON_INTERPOLABLE_KEYS:
            out[k] = va + (vb - va) * t
        else:
            if va != vb:
                raise ValueError(
                    f"taper property {k!r} cannot be interpolated "
                    f"(integer count, boolean, or string): {va!r} != {vb!r}. "
                    "Hold it constant across the taper."
                )
            out[k] = va
    return out


@dataclass
class TaperSpec:
    """
    Unresolved definition of an axially-tapered FMM grain.

    Attributes
    ----------
    grain_type : str
        openMotor geomName ('Finocyl', 'Star', ...).
    control_stations : list of (float, dict)
        Sorted `[(frac, props), ...]` with `frac` in [0, 1] (0 = forward
        / head end, 1 = aft / nozzle end). Two points = a linear taper;
        more points define a piecewise-linear (curved / multi-stage)
        profile. Cross-sections at intermediate fracs are obtained by
        interpolating between bracketing control points.
    map_dim : int
        Cross-sectional FMM resolution per station (radial, NOT axial).
    max_stations : int
        Upper bound on the number of FMM solves along the axis. The
        resolved station count is `min(segment_cells, max_stations)`.
    od_ends : list of dict or None
        OD / end-taper entries (openMotor `taper['od']['ends']` schema). When
        set, each station's `props['diameter']` is reduced to the local
        casting diameter (`motorlib.taper.od_diameter_at`) BEFORE the FMM
        runs, so the regression map / port area / wall_web are clipped to the
        shrinking casing (no FMM-internal change — mirrors the QS expander).
        Forces the per-station path even with one control point (constant
        cross-section, varying casing diameter).
    grain_length : float
        Axial length [m] used for the OD-taper end-region mapping. Set to the
        SNAPPED segment length by `build_snapped_geometry` so the transient
        and analytic `cell_D_outer` use an identical end-region extent.
    """
    grain_type: str
    control_stations: list
    map_dim: int = 1001
    max_stations: int = 32
    od_ends: object = None        # list[dict] OD/end-taper entries (None = none)
    grain_length: float = 0.0     # snapped axial length for the OD mapping

    def __post_init__(self):
        if not self.control_stations:
            raise ValueError("TaperSpec needs at least one control station.")
        # Sort by fraction and validate the range.
        self.control_stations = sorted(
            ((float(f), p) for f, p in self.control_stations),
            key=lambda fp: fp[0],
        )
        for f, _ in self.control_stations:
            if not (0.0 <= f <= 1.0):
                raise ValueError(
                    f"taper control-station fraction {f} outside [0, 1]."
                )


def linear_taper(grain_type: str, props_fwd: dict, props_aft: dict,
                 map_dim: int = 1001, max_stations: int = 32,
                 od_ends=None, grain_length: float = 0.0) -> TaperSpec:
    """Convenience: a two-point (forward → aft) linear taper."""
    return TaperSpec(
        grain_type=grain_type,
        control_stations=[(0.0, dict(props_fwd)), (1.0, dict(props_aft))],
        map_dim=map_dim,
        max_stations=max_stations,
        od_ends=od_ends,
        grain_length=grain_length,
    )


def taper_profile(grain_type: str, control_stations, map_dim: int = 1001,
                  max_stations: int = 32, od_ends=None,
                  grain_length: float = 0.0) -> TaperSpec:
    """Convenience: an arbitrary piecewise-linear taper from control points."""
    return TaperSpec(
        grain_type=grain_type,
        control_stations=[(float(f), dict(p)) for f, p in control_stations],
        map_dim=map_dim,
        max_stations=max_stations,
        od_ends=od_ends,
        grain_length=grain_length,
    )


def od_ends_from_taper(taper_def):
    """The OD / end-taper entries of a `.ric` `taper` block ([] when absent or
    disabled). Thin re-export of `motorlib.taper.od_ends_from_taper` that sets
    up the openMotor path first, so the adapter has one import surface."""
    if not isinstance(taper_def, dict):
        return []
    _setup_openmotor_path()
    from motorlib.taper import od_ends_from_taper as _od
    return _od(taper_def)


def taper_spec_from_props(grain_type: str, base_props: dict, taper_def: dict,
                          map_dim: int = 1001, max_stations: int = 32) -> 'TaperSpec':
    """
    Build a `TaperSpec` from an openMotor `.ric` grain's properties + its
    `taper` definition block (the solver-agnostic schema written by
    openMotor's `TaperProperty`). The grain's normal properties are the
    forward (frac 0) cross-section; each `bore.controlStations` entry gives
    the overrides at its `frac`. The `od` sub-block (end taper) is carried as
    `od_ends` so each per-station FMM table is clipped to the local casting
    diameter. Mirrors `motorlib.taper`'s control-point construction so QS and
    transient read identical geometry.
    """
    base = {k: v for k, v in base_props.items() if k != 'taper'}
    bore = (taper_def.get('bore', {}) or {}) if isinstance(taper_def, dict) else {}
    stations = sorted(bore.get('controlStations', []),
                      key=lambda s: float(s['frac']))
    control = [(0.0, dict(base))]
    for station in stations:
        overrides = station.get('props', {}) or {}
        control.append((float(station['frac']), {**base, **overrides}))
    od_ends = od_ends_from_taper(taper_def)
    grain_length = float(base.get('length', 0.0))
    return taper_profile(grain_type, control, map_dim=map_dim,
                         max_stations=max_stations,
                         od_ends=(od_ends or None), grain_length=grain_length)


def _interp_control(control_stations: list, frac: float) -> dict:
    """Cross-section property dict at axial `frac` from sorted control points."""
    if len(control_stations) == 1:
        return dict(control_stations[0][1])
    if frac <= control_stations[0][0]:
        return dict(control_stations[0][1])
    if frac >= control_stations[-1][0]:
        return dict(control_stations[-1][1])
    for j in range(len(control_stations) - 1):
        f0, p0 = control_stations[j]
        f1, p1 = control_stations[j + 1]
        if f0 <= frac <= f1:
            t = 0.0 if f1 == f0 else (frac - f0) / (f1 - f0)
            return _interpolate_props(p0, p1, t)
    return dict(control_stations[-1][1])  # unreachable; defensive


def _closed_fmm_table(props: dict, grain_type: str) -> FmmTable:
    """A degenerate, fully-closed FMM table for a near-closed OD-tip station.

    When an OD / end taper shrinks the casing to (near) the core — e.g. a true
    hemispherical dome with `endDiameter = 0` — the cross-section is essentially
    consumed (web → 0) and openMotor's FMM pipeline can't characterize it
    (`savgol_filter` needs ≥ 31 face-area samples, which a sub-pixel web can't
    provide). Such a tip carries negligible propellant, so model it as CLOSED:
    zero burning perimeter, port == casting (no propellant), and a tiny web so
    any cell mapping here is already burnt out. This mirrors the QS expander,
    which dodges the tip entirely by sampling slice centers."""
    god = float(props.get('diameter', 0.0))
    casting = np.pi / 4.0 * god * god
    w = 1.0e-4  # tiny positive web (cells mapping here contribute no burn)
    inh_fwd, inh_aft = _INHIBIT_MAP.get(props.get('inhibitedEnds', 'Neither'),
                                        (False, False))
    return FmmTable(
        reg_depth=np.array([0.0, w]),
        perimeter=np.array([0.0, 0.0]),
        port_area=np.array([casting, casting]),
        wall_web=w,
        grain_outer_diameter=god,
        grain_length=float(props.get('length', 0.0)),
        inhibited_fwd=inh_fwd,
        inhibited_aft=inh_aft,
        geom_name=grain_type,
    )


def resolve_taper(taper: TaperSpec, n_stations: int):
    """
    Build the real per-station FMM tables for a taper.

    Runs openMotor's FMM pipeline once per station (deduplicated when
    consecutive interpolated cross-sections are identical, e.g. a
    degenerate fwd == aft taper). Stations are placed at
    `np.linspace(0, 1, n_stations)` so the endpoints reproduce the exact
    forward/aft cross-sections. For an OD / end taper that closes a tip to
    ~the core, the degenerate tip station(s) become a `_closed_fmm_table`
    (no burn, no propellant) instead of crashing openMotor's FMM pipeline.

    Parameters
    ----------
    taper : TaperSpec
    n_stations : int
        Number of axial stations (>= 1). Decided by the caller from the
        segment's snapped cell count.

    Returns
    -------
    (tables, station_frac) : (list[FmmTable], np.ndarray)
        `tables[m]` is the cross-section at `station_frac[m]`. The list may
        contain repeated table references where stations are identical.
    """
    _setup_openmotor_path()

    n = max(1, int(n_stations))
    if n == 1:
        fracs = np.array([0.5])
    else:
        fracs = np.linspace(0.0, 1.0, n)

    # OD / end taper: each station's casting diameter is reduced to the local
    # value BEFORE the FMM runs (the mask then clips the cross-section). Same
    # analytic profile the QS expander and the transient cell_D_outer use.
    od_ends = getattr(taper, 'od_ends', None)
    od_diameter_at = None
    if od_ends:
        from motorlib.taper import od_diameter_at as _od_at
        od_diameter_at = _od_at

    tables = []
    cache = {}  # interpolated-props signature -> FmmTable (avoid resolves)
    for f in fracs:
        props = _interp_control(taper.control_stations, float(f))
        if od_diameter_at is not None:
            full_d = props['diameter']
            d = od_diameter_at(float(f), taper.grain_length, full_d, od_ends)
            # Keep a tiny positive web at a near-closed tip so the FMM doesn't
            # choke (mirrors the QS expander's min-diameter clamp).
            min_d = props.get('coreDiameter', 0.0) + 2.0e-4
            props['diameter'] = max(d, min_d)
        sig = tuple(
            (k, round(v, 12) if isinstance(v, float) else v)
            for k, v in sorted(props.items())
        )
        tab = cache.get(sig)
        if tab is None:
            ric = {'type': taper.grain_type, 'properties': props}
            if od_ends:
                # A near-closed OD tip can give the FMM too thin a web; fall
                # back to a closed table rather than crash the solver.
                try:
                    tab = from_ric_grain(ric, map_dim=taper.map_dim)
                except (ValueError, ZeroDivisionError):
                    tab = _closed_fmm_table(props, taper.grain_type)
            else:
                tab = from_ric_grain(ric, map_dim=taper.map_dim)
            cache[sig] = tab
        tables.append(tab)
    return tables, fracs


# ================================================================
# Numba-compiled lookup helper
# ================================================================

@njit(cache=True)
def fmm_table_lookup(regress, reg_arr, val_arr, n_samples):
    """
    O(1) linear interpolation into an FMM table sampled on a UNIFORM
    regression-depth grid (which is what `from_openmotor` produces via
    np.linspace). Used inside the @njit hot loop for per-cell
    perimeter/port-area lookups.

    Parameters
    ----------
    regress : float
        Regression depth at this cell, [m].
    reg_arr : ndarray (n_samples,)
        Uniform regression-depth grid from 0 to wall_web.
    val_arr : ndarray (n_samples,)
        Tabulated values to interpolate (perimeter or port_area).
    n_samples : int
        Length of reg_arr / val_arr.

    Returns
    -------
    float
        Linearly interpolated value. Clamped at the table endpoints.
    """
    if regress <= reg_arr[0]:
        return val_arr[0]
    if regress >= reg_arr[n_samples - 1]:
        return val_arr[n_samples - 1]

    # Uniform grid → O(1) index calc
    dr = reg_arr[1] - reg_arr[0]
    idx_f = (regress - reg_arr[0]) / dr
    idx_lo = int(idx_f)
    if idx_lo >= n_samples - 1:
        return val_arr[n_samples - 1]
    f = idx_f - idx_lo
    v0 = val_arr[idx_lo]
    v1 = val_arr[idx_lo + 1]
    return v0 + f * (v1 - v0)
