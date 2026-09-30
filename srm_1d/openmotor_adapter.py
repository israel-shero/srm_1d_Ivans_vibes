"""
openmotor_adapter.py — openMotor .ric File Adapter
=====================================================

Reads openMotor .ric files (YAML format), converts motor definitions
into srm_1d data structures, runs simulations, and exports results
as openMotor-compatible CSV.

Standalone — does not require motorlib to be installed.

Usage:
    from srm_1d.openmotor_adapter import run_from_ric

    result, perf = run_from_ric(
        "my_motor.ric",
        gas_props={'mu': 8.842e-5, 'k': 0.3685, 'Cp': 2060.0},
    )

    # Or step-by-step:
    from srm_1d.openmotor_adapter import load_ric, ric_to_sim_args, result_to_csv

    motor = load_ric("my_motor.ric")
    args = ric_to_sim_args(motor, gas_props={...})
    nozzle = args['nozzle']
    geo = args.pop('geo')
    prop = args.pop('propellant')
    result = run_simulation(geo, prop, **args)
    csv_str = result_to_csv(result, perf)

Unit conversions from openMotor internal units:
    erosionCoeff: m/(s·Pa) → μm/(s·MPa)  (multiply by 1e12)
    slagCoeff:    (m·Pa)/s → (m·MPa)/s  (divide by 1e6)
    propellant.m: g/mol → kg/mol  (divide by 1000)
    propellant.a: m/s per Pa^n  (same as ours, no conversion)
"""

import os
import warnings
import numpy as np

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

from .propellant import Propellant, PropellantTab, Pyrogen
from .grain_geometry import build_snapped_geometry
from .igniter_plenum import PyrogenChamber, sutton_pyrogen_mass
from .nozzle import Nozzle, compute_motor_performance, print_performance_summary
from .simulation import run_simulation


# ================================================================
# .ric file reader
# ================================================================

def load_ric(filepath):
    """
    Load an openMotor .ric file.

    Parameters
    ----------
    filepath : str
        Path to the .ric file.

    Returns
    -------
    dict with keys 'grains', 'nozzle', 'propellant', 'config',
    'version'. Raw openMotor format, not yet converted.
    """
    if not HAS_YAML:
        raise ImportError(
            "PyYAML is required to read .ric files. "
            "Install with: pip install pyyaml"
        )

    with open(filepath, 'r', encoding='utf-8') as f:
        # openMotor .ric files contain Python-specific YAML tags
        # (!!python/object/apply, !!python/tuple) that safe_load
        # rejects. Use a custom loader that converts these to plain
        # Python objects — we only need the data dict.
        class _RicLoader(yaml.SafeLoader):
            pass
        # Handle any unknown tag by returning its constructed value
        _RicLoader.add_multi_constructor(
            'tag:yaml.org,2002:python/',
            lambda loader, suffix, node: loader.construct_mapping(node)
            if isinstance(node, yaml.MappingNode)
            else loader.construct_sequence(node)
            if isinstance(node, yaml.SequenceNode)
            else loader.construct_scalar(node),
        )
        raw = yaml.load(f, Loader=_RicLoader)

    if 'data' not in raw:
        raise ValueError(f"Invalid .ric file: missing 'data' key in {filepath}")

    data = raw['data']
    required = ['grains', 'nozzle', 'propellant', 'config']
    for key in required:
        if key not in data:
            raise ValueError(f"Invalid .ric file: missing '{key}' in {filepath}")

    version = raw.get('version', (0, 0, 0))
    # Upgrade pre-(0,7,0) propellants in memory so transport keys are
    # always present (sentinel-seeded) — mirrors openMotor's fileIO
    # migration. Real values come either from an already-migrated .ric or
    # from migrate_ric_transport folding in a .transport.yaml sidecar.
    _migrate_ric_propellant(data['propellant'], version)

    return {
        'grains': data['grains'],
        'nozzle': data['nozzle'],
        'propellant': data['propellant'],
        'config': data['config'],
        # v0.8.0 Phase 4: igniter block (None on un-migrated motors).
        'igniter': data.get('igniter'),
        'version': version,
    }


def load_transport(transport_path):
    """
    Load a srm_1d transport YAML sibling file (alongside a .ric).

    The transport YAML is a srm_1d-specific extension that supplies
    combustion gas transport properties not present in the openMotor
    schema (.ric files only carry combustion thermo: γ, T_flame, MW).

    Schema:
        mu: <Pa·s>
        k:  <W/(m·K)>
        Cp: <J/(kg·K)>

    Returns
    -------
    dict shaped like ``gas_props`` for ``convert_propellant``.
    """
    if not HAS_YAML:
        raise ImportError(
            "PyYAML is required to read transport YAML files."
        )
    with open(transport_path, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)
    required = {'mu', 'k', 'Cp'}
    missing = required - set(data.keys())
    if missing:
        raise ValueError(
            f"Transport YAML {transport_path} missing required keys: "
            f"{sorted(missing)}"
        )
    return {'mu': data['mu'], 'k': data['k'], 'Cp': data['Cp']}


# ================================================================
# v0.8.0 Phase 3 — per-tab gas transport in the .ric format
# ================================================================
# openMotor's v0.7.0 file format carries combustion-gas transport on each
# propellant tab: a single ``mu`` (viscosity is invariant of the
# frozen/effective equilibrium-chemistry shift) plus frozen + effective
# conductivity/specific-heat pairs, with a propellant-level
# ``transportVariant`` selector. The default 0.0 is the "not provided"
# sentinel (D7): srm_1d refuses to fabricate transport and hard-faults
# instead. This block reads that schema and migrates old files
# (mirroring openMotor's chained uilib/fileIO migration).

RIC_FORMAT_VERSION = (0, 7, 0)
TRANSPORT_UNSET = 0.0
_TRANSPORT_TAB_KEYS = (
    'mu', 'kThermalFrozen', 'cpFrozen', 'kThermalEffective', 'cpEffective',
)


def _seed_tab_transport_sentinels(propellant):
    """Mirror openMotor's migrateMotor_0_6_1_to_0_7_0: seed the transport
    keys (sentinel) + transportVariant on a propellant dict if absent. In
    place; returns the propellant."""
    propellant.setdefault('transportVariant', 'frozen')
    for tab in propellant.get('tabs', []):
        for key in _TRANSPORT_TAB_KEYS:
            tab.setdefault(key, TRANSPORT_UNSET)
    return propellant


def _migrate_ric_propellant(propellant, version):
    """Upgrade an in-memory .ric propellant dict to the current format.
    For pre-(0,7,0) files this seeds the transport sentinels so downstream
    reads always see the keys (and hard-fault on the sentinel rather than
    KeyError)."""
    if tuple(version) < RIC_FORMAT_VERSION:
        _seed_tab_transport_sentinels(propellant)
    return propellant


def _tab_transport(tab, variant):
    """Return ``{'mu','k','Cp'}`` for a .ric tab dict and variant
    ('effective'|'frozen'), or None if any value is the unset sentinel."""
    mu = float(tab.get('mu', TRANSPORT_UNSET))
    if variant == 'frozen':
        k = float(tab.get('kThermalFrozen', TRANSPORT_UNSET))
        cp = float(tab.get('cpFrozen', TRANSPORT_UNSET))
    else:
        k = float(tab.get('kThermalEffective', TRANSPORT_UNSET))
        cp = float(tab.get('cpEffective', TRANSPORT_UNSET))
    if mu <= 0.0 or k <= 0.0 or cp <= 0.0:
        return None
    return {'mu': mu, 'k': k, 'Cp': cp}


def _ric_raw_loader():
    """A SafeLoader that flattens openMotor's python tags (the version
    tuple and the fileTypes enum) to plain Python so the full
    {version, type, data} structure can be read for migration."""
    class _RicLoader(yaml.SafeLoader):
        pass
    _RicLoader.add_multi_constructor(
        'tag:yaml.org,2002:python/',
        lambda loader, suffix, node: loader.construct_mapping(node)
        if isinstance(node, yaml.MappingNode)
        else loader.construct_sequence(node)
        if isinstance(node, yaml.SequenceNode)
        else loader.construct_scalar(node),
    )
    return _RicLoader


class _RicFileType:
    """Stand-in for openMotor's ``fileTypes`` enum so we can re-emit the
    ``!!python/object/apply:uilib.fileIO.fileTypes [code]`` tag on write
    without importing uilib (which needs PyQt). MOTOR == 3."""
    def __init__(self, code):
        self.code = int(code)


def _ric_dumper():
    class _Dumper(yaml.Dumper):
        pass
    _Dumper.add_representer(
        _RicFileType,
        lambda dumper, data: dumper.represent_sequence(
            'tag:yaml.org,2002:python/object/apply:uilib.fileIO.fileTypes',
            [data.code]),
    )
    return _Dumper


def _load_sidecar(path):
    return load_transport(path) if (path and os.path.exists(path)) else None


def build_transport_library(motors_dir):
    """Group .ric files by propellant name and resolve a shared
    frozen+effective transport set per propellant.

    Propellants are shared across motors (e.g. several motors use the same
    'Chunc' propellant), so a sidecar found on ANY motor of a propellant
    supplies transport for ALL of them — filling motors that ship no
    sidecar of their own.

    Per propellant:
    - frozen  ← any ``<stem>.frozen.transport.yaml``; else the primary
      ``<stem>.transport.yaml`` (treat the lone primary as frozen, matching
      the 'frozen' default).
    - effective ← any primary ``<stem>.transport.yaml``; else the frozen.
    - mu (invariant) ← whichever sidecar is available.

    Returns ``{pname: {'mu','frozen':{'k','Cp'},'effective':{'k','Cp'}}}``,
    or ``{pname: None}`` when a propellant has no transport anywhere.
    """
    import glob
    from collections import defaultdict
    groups = defaultdict(list)
    for ric in sorted(glob.glob(os.path.join(motors_dir, '*.ric'))):
        try:
            motor = load_ric(ric)
        except Exception:
            continue
        groups[motor['propellant'].get('name', '?')].append(os.path.splitext(ric)[0])

    library = {}
    for pname, stems in groups.items():
        eff = froz = None
        for stem in stems:
            if eff is None:
                eff = _load_sidecar(stem + '.transport.yaml')
            if froz is None:
                froz = _load_sidecar(stem + '.frozen.transport.yaml')
        if eff is None and froz is None:
            library[pname] = None
            continue
        mu = (froz or eff)['mu']
        frozen_kcp = froz if froz is not None else eff
        effective_kcp = eff if eff is not None else froz
        library[pname] = {
            'mu': mu,
            'frozen': {'k': frozen_kcp['k'], 'Cp': frozen_kcp['Cp']},
            'effective': {'k': effective_kcp['k'], 'Cp': effective_kcp['Cp']},
        }
    return library


def migrate_ric_transport(ric_path, transport=None, variant='frozen',
                          write=True, out_path=None):
    """Write per-tab transport (v0.7.0 format) into a .ric file, set
    ``transportVariant``, bump the version, and (optionally) save —
    preserving openMotor's YAML tags so the file stays openMotor-loadable.

    Parameters
    ----------
    ric_path : str
    transport : dict or None
        Resolved transport ``{'mu', 'frozen':{'k','Cp'},
        'effective':{'k','Cp'}}`` (e.g. from :func:`build_transport_library`,
        which shares across same-propellant motors). If None, falls back to
        THIS motor's own ``<stem>.transport.yaml`` / ``.frozen.`` siblings
        (no sharing). If still nothing, tabs keep the unset sentinel and the
        motor will hard-fault on run (D7).
    variant : str
        The active ``transportVariant`` written into the file ('frozen' |
        'effective'). Default 'frozen'.
    write : bool
    out_path : str or None

    Returns
    -------
    dict : the migrated raw {version, type, data} structure.
    """
    if not HAS_YAML:
        raise ImportError("PyYAML is required to migrate .ric files.")

    if transport is None:
        stem, _ = os.path.splitext(ric_path)
        eff = _load_sidecar(stem + '.transport.yaml')
        froz = _load_sidecar(stem + '.frozen.transport.yaml')
        if eff is not None or froz is not None:
            f_kcp = froz if froz is not None else eff
            e_kcp = eff if eff is not None else froz
            transport = {
                'mu': (froz or eff)['mu'],
                'frozen': {'k': f_kcp['k'], 'Cp': f_kcp['Cp']},
                'effective': {'k': e_kcp['k'], 'Cp': e_kcp['Cp']},
            }

    with open(ric_path, 'r', encoding='utf-8') as f:
        raw = yaml.load(f, Loader=_ric_raw_loader())

    propellant = raw['data']['propellant']
    propellant['transportVariant'] = variant
    for tab in propellant.get('tabs', []):
        if transport is None:
            for key in _TRANSPORT_TAB_KEYS:
                tab.setdefault(key, TRANSPORT_UNSET)
        else:
            tab['mu'] = transport['mu']
            tab['kThermalFrozen'] = transport['frozen']['k']
            tab['cpFrozen'] = transport['frozen']['Cp']
            tab['kThermalEffective'] = transport['effective']['k']
            tab['cpEffective'] = transport['effective']['Cp']

    # v0.8.0 Phase 4: seed a default igniter block so the motor is
    # self-describing (bpnv + forward_plenum auto-sizing) if absent.
    _seed_motor_igniter(raw['data'])

    raw['version'] = tuple(RIC_FORMAT_VERSION)
    if 'type' in raw:
        code = raw['type'][0] if isinstance(raw['type'], (list, tuple)) else raw['type']
        raw['type'] = _RicFileType(code)
    else:
        raw['type'] = _RicFileType(3)  # fileTypes.MOTOR

    if write:
        target = out_path or ric_path
        with open(target, 'w', encoding='utf-8') as f:
            yaml.dump(raw, f, Dumper=_ric_dumper(), default_flow_style=False)

    return raw


def migrate_all_motors(motors_dir, variant='frozen', write=True, verbose=True):
    """Migrate every .ric in ``motors_dir`` to the v0.7.0 transport format,
    sharing transport across motors with the same propellant (see
    :func:`build_transport_library`). Returns a list of
    ``(filename, propellant, had_transport)`` tuples."""
    import glob
    library = build_transport_library(motors_dir)
    results = []
    for ric in sorted(glob.glob(os.path.join(motors_dir, '*.ric'))):
        try:
            pname = load_ric(ric)['propellant'].get('name', '?')
        except Exception as exc:
            results.append((os.path.basename(ric), f'ERR:{exc}', False))
            continue
        transport = library.get(pname)
        migrate_ric_transport(ric, transport=transport, variant=variant, write=write)
        results.append((os.path.basename(ric), pname, transport is not None))
        if verbose:
            tag = 'transport' if transport is not None else 'SENTINEL (no data)'
            print(f"  {os.path.basename(ric):<52} {pname:<24} {tag}")
    return results


def _builtin_pyrogen_path(name):
    return os.path.join(
        os.path.dirname(__file__), 'pyrogens', f'{name}.yaml'
    )


def load_pyrogen(path_or_name):
    """
    Load a pyrogen YAML datasheet.

    ``path_or_name`` may be a filesystem path or a built-in pyrogen name
    such as ``"bpnv"``.
    """
    if not HAS_YAML:
        raise ImportError(
            "PyYAML is required to read pyrogen YAML files."
        )

    candidate = path_or_name
    if isinstance(path_or_name, str) and not os.path.exists(candidate):
        builtin = _builtin_pyrogen_path(path_or_name.lower())
        if os.path.exists(builtin):
            candidate = builtin

    if not os.path.exists(candidate):
        raise ValueError(
            f"Unknown pyrogen '{path_or_name}'. Provide a YAML path or one "
            "of the built-in pyrogen names such as 'bpnv' or 'mtv'."
        )

    with open(candidate, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)

    required = {'name', 'a', 'n', 'rho', 'T_flame', 'M', 'gamma'}
    missing = required - set(data.keys())
    if missing:
        raise ValueError(
            f"Pyrogen YAML {candidate} missing required keys: "
            f"{sorted(missing)}"
        )

    return Pyrogen(
        name=str(data['name']),
        a=float(data['a']),
        n=float(data['n']),
        rho=float(data['rho']),
        T_flame=float(data['T_flame']),
        M=float(data['M']),
        gamma=float(data['gamma']),
        impetus_W=float(data.get('impetus_W', 0.0)),
        heat_flux_cal_cm2_s=(
            None if data.get('heat_flux_cal_cm2_s') is None
            else float(data['heat_flux_cal_cm2_s'])
        ),
        kappa_jet=float(data.get('kappa_jet', 8.0)),
        # v0.8.x condensed-phase gas fraction (ProPep/CEA gas_mass/charge_mass).
        gas_mass_fraction=float(data.get('gas_mass_fraction', 1.0)),
        # v0.7.3 Phase B.3 / B.4: optional YAML fields with defensible
        # defaults baked into the Pyrogen dataclass.
        form=str(data.get('form', 'pellets')),
        heat_delivery_mode=str(data.get('heat_delivery_mode', 'demar')),
        pellet_emissivity=float(data.get('pellet_emissivity', 0.7)),
        radiation_absorption_length_m=float(
            data.get('radiation_absorption_length_m', 1.0)
        ),
        # v0.7.4 Phase C.1: explicit particle dimensions (replaces
        # form-archetype A_burn multipliers).
        particle_diameter_m=float(data.get('particle_diameter_m', 5.0e-3)),
        particle_LD_ratio=float(data.get('particle_LD_ratio', 3.0)),
    )


def build_pyrogen_chamber(
    pyrogen, geo, nozzle,
    pyrogen_mass=None,
    pyrogen_throat_area=None,
    pyrogen_volume=None,
    pyrogen_burn_area=None,
    pyrogen_burn_law='0d',
    injection_topology='forward_plenum',
    cartridge_length_m=-1.0,
    basket_fill_fraction=0.5,
    pellet_packing_fraction=0.60,
):
    """
    Build a PyrogenChamber using v0.7.0 default sizing rules.

    v0.7.3 Phase A — ``injection_topology`` selects between the
    plenum-with-orifice (``'forward_plenum'``, v0.7.0+ default) and
    uncontained submerged-pyrogen models (``'head_basket'`` /
    ``'aft_basket'``). For uncontained topologies, ``A_throat`` and
    ``V_plenum`` are vestigial (PyrogenChamber still validates them at
    construction so existing motor configs don't break, but the time
    loop ignores them). ``cartridge_length_m=-1.0`` (default sentinel)
    derives the cartridge length from pyrogen mass via
    ``L_cart = m_pyrogen / (rho_p * A_port_avg)`` at sim init.
    """
    if pyrogen_mass is None:
        case_volume = np.pi / 4.0 * geo.D_outer ** 2 * geo.L_motor
        free_volume_m3 = max(case_volume - geo.total_propellant_volume(), 0.0)
        free_volume_in3 = free_volume_m3 / (0.0254 ** 3)
        pyrogen_mass = sutton_pyrogen_mass(free_volume_in3)

    if pyrogen_volume is None:
        pyrogen_volume = 1.5 * pyrogen_mass / pyrogen.rho

    if pyrogen_burn_area is None:
        # v0.7.4 Phase C.1: compute total burning surface area from
        # physical particle geometry (replaces v0.7.3 Phase B.3
        # ×1/×5/×20 form-archetype multipliers).
        d = float(pyrogen.particle_diameter_m)
        ld = float(pyrogen.particle_LD_ratio)
        if d <= 0.0:
            raise ValueError(
                f"Pyrogen '{pyrogen.name}': particle_diameter_m must be "
                f"positive; got {d}"
            )
        if ld <= 1.0:
            # Sphere of diameter d: A_total = 6·m / (ρ·d)
            pyrogen_burn_area = 6.0 * pyrogen_mass / (pyrogen.rho * d)
        else:
            # Cylinder D=d, L=ld·d: A_total = m·(4λ+2) / (ρ·λ·d)
            pyrogen_burn_area = (
                pyrogen_mass * (4.0 * ld + 2.0)
                / (pyrogen.rho * ld * d)
            )

    if pyrogen_throat_area is None:
        # v0.7.3.2 (2026-05-27): switch from A_throat = 0.01 * A_main
        # (which gave 9 mm² for Hasegawa A, far too small under
        # Phase B.0's cold-bore IC — drove plenum P_ig past 600 MPa
        # equilibrium and tripped numerical collapse within ~2 ms)
        # to a Kn-based sizing rule:
        #     A_throat = A_burn / Kn_design
        # where Kn_design = 100 is the Sutton 9e §14.5 mid-range for
        # pellet-form BKNO3 / MTV pyrogens at the 5-30 MPa working
        # pressure target. With A_burn ≈ 46 cm² for Hasegawa A
        # (Mizushima 3.2 mm pellets), this gives A_throat ≈ 46 mm²
        # — close to the calibrated test value (38.5 mm²).
        # Lower bound 1 mm² preserves choked-flow validity at the
        # very small pyrogens; upper bound 100 mm² prevents the
        # opposite degenerate (no throat confinement).
        # See srm_1d/docs/v0_7_4/references/frozen_collapse_investigation.md
        # for the regression timeline and physical motivation.
        Kn_design_pyrogen = 100.0
        pyrogen_throat_area = pyrogen_burn_area / Kn_design_pyrogen
        pyrogen_throat_area = min(max(pyrogen_throat_area, 1.0e-6), 1.0e-4)

    return PyrogenChamber(
        pyrogen=pyrogen,
        m_pyrogen_initial=pyrogen_mass,
        A_burn_initial=pyrogen_burn_area,
        A_throat=pyrogen_throat_area,
        V_plenum=pyrogen_volume,
        burn_law=pyrogen_burn_law,
        injection_topology=injection_topology,
        cartridge_length_m=cartridge_length_m,
        basket_fill_fraction=basket_fill_fraction,
        pellet_packing_fraction=pellet_packing_fraction,
    )


# ================================================================
# v0.8.0 Phase 4 — igniter as data (library + motor block, D3)
# ================================================================
# A motor's igniter is a ``data.igniter`` block mirroring how the motor
# embeds its propellant: an embedded pyrogen MATERIAL datasheet (the
# reusable library item, same field set as srm_1d/pyrogens/*.yaml)
# plus the per-motor chamber SIZING / topology. ``-1.0`` sizing values are
# the "auto" sentinel (Sutton sizing / derive at build). openMotor ignores
# the unknown ``data.igniter`` key on load (Phase 5 wires its GUI to it).

# Sizing defaults match build_pyrogen_chamber's auto behavior (so a
# migrated motor reproduces the current forward_plenum default run).
_IGNITER_SIZING_DEFAULTS = {
    'mass': -1.0,                 # -1 => Sutton auto
    'throat_area': -1.0,          # -1 => Kn-design auto
    'volume': -1.0,               # -1 => 1.5*m/rho auto
    'burn_area': -1.0,            # -1 => particle-geometry auto
    'burn_law': '0d',
    'injection_topology': 'forward_plenum',
    'cartridge_length_m': -1.0,   # -1 => derive at sim init
    'basket_fill_fraction': 0.5,
    'pellet_packing_fraction': 0.60,
}


def _pyrogen_to_block(pyrogen):
    """Serialize a Pyrogen dataclass to the ``igniter.pyrogen`` dict."""
    return {
        'name': pyrogen.name, 'a': pyrogen.a, 'n': pyrogen.n,
        'rho': pyrogen.rho, 'T_flame': pyrogen.T_flame, 'M': pyrogen.M,
        'gamma': pyrogen.gamma, 'impetus_W': pyrogen.impetus_W,
        'heat_flux_cal_cm2_s': pyrogen.heat_flux_cal_cm2_s,
        'Cp_gas': pyrogen.Cp_gas, 'kappa_jet': pyrogen.kappa_jet,
        'form': pyrogen.form,
        'particle_diameter_m': pyrogen.particle_diameter_m,
        'particle_LD_ratio': pyrogen.particle_LD_ratio,
        'heat_delivery_mode': pyrogen.heat_delivery_mode,
        'pellet_emissivity': pyrogen.pellet_emissivity,
        'radiation_absorption_length_m': pyrogen.radiation_absorption_length_m,
    }


def _block_to_pyrogen(pdict):
    """Build a Pyrogen dataclass from an ``igniter.pyrogen`` dict."""
    return Pyrogen(
        name=str(pdict['name']), a=float(pdict['a']), n=float(pdict['n']),
        rho=float(pdict['rho']), T_flame=float(pdict['T_flame']),
        M=float(pdict['M']), gamma=float(pdict['gamma']),
        impetus_W=float(pdict.get('impetus_W', 0.0)),
        heat_flux_cal_cm2_s=(None if pdict.get('heat_flux_cal_cm2_s') is None
                             else float(pdict['heat_flux_cal_cm2_s'])),
        Cp_gas=(None if pdict.get('Cp_gas') is None else float(pdict['Cp_gas'])),
        kappa_jet=float(pdict.get('kappa_jet', 8.0)),
        form=str(pdict.get('form', 'pellets')),
        particle_diameter_m=float(pdict.get('particle_diameter_m', 5.0e-3)),
        particle_LD_ratio=float(pdict.get('particle_LD_ratio', 3.0)),
        heat_delivery_mode=str(pdict.get('heat_delivery_mode', 'demar')),
        pellet_emissivity=float(pdict.get('pellet_emissivity', 0.7)),
        radiation_absorption_length_m=float(
            pdict.get('radiation_absorption_length_m', 1.0)),
    )


def default_igniter_block(pyrogen_name='bpnv'):
    """A default ``data.igniter`` block: the named built-in pyrogen material
    + auto chamber sizing (forward_plenum). Used by the migration to make
    motors self-describing."""
    block = {'pyrogen': _pyrogen_to_block(load_pyrogen(pyrogen_name))}
    block.update(dict(_IGNITER_SIZING_DEFAULTS))
    return block


def load_igniter(motor):
    """Return ``(Pyrogen, sizing_kwargs)`` from a motor dict's ``igniter``
    block, or None if the motor has no igniter block. ``sizing_kwargs`` maps
    directly onto :func:`build_pyrogen_chamber` (``-1.0`` mass/throat/volume/
    burn_area sentinels become None = auto)."""
    ig = motor.get('igniter')
    if ig is None:
        return None
    pyro = _block_to_pyrogen(ig['pyrogen'])

    def _auto(value):
        return None if (value is None or value == -1.0) else value

    sizing = {
        'pyrogen_mass': _auto(ig.get('mass', -1.0)),
        'pyrogen_throat_area': _auto(ig.get('throat_area', -1.0)),
        'pyrogen_volume': _auto(ig.get('volume', -1.0)),
        'pyrogen_burn_area': _auto(ig.get('burn_area', -1.0)),
        'pyrogen_burn_law': ig.get('burn_law', '0d'),
        'injection_topology': ig.get('injection_topology', 'forward_plenum'),
        'cartridge_length_m': ig.get('cartridge_length_m', -1.0),
        'basket_fill_fraction': ig.get('basket_fill_fraction', 0.5),
        'pellet_packing_fraction': ig.get('pellet_packing_fraction', 0.60),
    }
    return pyro, sizing


def _seed_motor_igniter(data, pyrogen_name='bpnv'):
    """Migration: add a default igniter block to a motor data dict if absent.
    Mirrors openMotor's migrateMotor seeding pattern. In place; returns data."""
    if 'igniter' not in data:
        data['igniter'] = default_igniter_block(pyrogen_name)
    return data


# ================================================================
# Conversion: openMotor → srm_1d data structures
# ================================================================

# Maps openMotor inhibitedEnds strings to (inhibit_fwd, inhibit_aft)
_INHIBIT_MAP = {
    'Neither': (False, False),
    'Top': (True, False),       # Top = forward (head end)
    'Bottom': (False, True),    # Bottom = aft (nozzle end)
    'Both': (True, True),
}


def _ric_tab_to_srm(ric_tab):
    """Convert one openMotor tab dict to srm_1d PropellantTab."""
    return PropellantTab(
        min_pressure=ric_tab.get('minPressure', 0.0),
        max_pressure=ric_tab.get('maxPressure', 2.0e7),
        a=ric_tab['a'],                 # m/s per Pa^n (same convention)
        n=ric_tab['n'],
        gamma=ric_tab['k'],             # openMotor: 'k' = γ
        T_flame=ric_tab['t'],           # openMotor: 't' = T_flame
        molecular_weight=ric_tab['m'] / 1000.0,  # g/mol → kg/mol
    )


def _representative_ric_tab(tabs):
    """Pick the tab spanning the widest pressure range — used to source
    transport-property estimation when gas_props isn't supplied."""
    if len(tabs) == 1:
        return tabs[0]
    return max(tabs,
               key=lambda t: t.get('maxPressure', 0) - t.get('minPressure', 0))


def _default_radiation_emissivity(propellant_name):
    """Return the adjacent-ignition radiation emissivity default.

    Always returns 0.0 -- adjacent-cell radiation is opt-in. The Phase 4
    radiation sweep on Hasegawa A showed the constant-T_flame -> T[neighbor]
    chain drives ignition spread at ~1 ms/cell, pushing interior flow
    supersonic faster than the signed-throat PISO boundary can vent.
    Sutton 9e Section 15.3 also documents pyrogen ignition as primarily
    convective rather than radiative. Set ``radiation_emissivity`` in the
    .ric file (or override on the Propellant) to opt back in once the
    spread-rate / numerical-stability interaction is understood.
    """
    return 0.0


def convert_propellant(ric_propellant, gas_props=None):
    """
    Convert an openMotor propellant dict to srm_1d Propellant. All
    `tabs` in the .ric file are preserved one-to-one as PropellantTab
    entries (full multi-tab support since v0.4.0).

    Parameters
    ----------
    ric_propellant : dict
        From the .ric file's 'propellant' key.
    gas_props : dict or None
        Explicit transport override {'mu': Pa·s, 'k': W/(m·K), 'Cp':
        J/(kg·K)}. If None, transport is read from the propellant's own
        per-tab schema (v0.7.0 .ric format) using the active
        ``transportVariant``. A missing/sentinel value hard-faults (D7) —
        srm_1d does not fabricate transport.

    Returns
    -------
    Propellant
    """
    ric_tabs = ric_propellant['tabs']
    tabs = [_ric_tab_to_srm(t) for t in ric_tabs]

    name = ric_propellant.get('name', 'openMotor propellant')
    density = ric_propellant['density']

    if gas_props is None:
        # Read transport from the .ric propellant (v0.7.0 per-tab schema),
        # collapsed to the representative tab for the scalar solver.
        variant = ric_propellant.get('transportVariant', 'frozen')
        rep = _representative_ric_tab(ric_tabs)
        gas_props = _tab_transport(rep, variant)
        if gas_props is None:
            raise ValueError(
                f"Propellant {name!r}: {variant} gas transport "
                f"(mu/k/Cp) is not provided in the .ric file (unset "
                f"sentinel). Supply it by migrating the .ric with "
                f"migrate_ric_transport (folds a .transport.yaml sibling), "
                f"or pass gas_props={{'mu','k','Cp'}} explicitly. srm_1d "
                f"does not fabricate transport properties."
            )
    mu = gas_props['mu']
    k_gas = gas_props['k']
    Cp = gas_props['Cp']

    return Propellant(
        name=name,
        tabs=tabs,
        rho_propellant=density,
        Cps=1500.0,        # Default solid Cp — not in .ric
        T_surface=1000.0,  # Default surface temp — not in .ric
        T_initial=293.0,
        mu_gas=mu,
        k_gas=k_gas,
        Cp_gas=Cp,
        radiation_emissivity=float(
            ric_propellant.get(
                'radiation_emissivity',
                _default_radiation_emissivity(name),
            )
        ),
    )


def convert_geometry(ric_grains, target_propellant_cells=100,
                     fmm_map_dim=1001):
    """
    Convert openMotor grain list to srm_1d MotorGeometry.
    The throat lives on a separate Nozzle object — see convert_nozzle().

    Supports BATES + Conical (analytic) and all 7 of openMotor's FMM
    grain types (Finocyl, Star, Moonburner, X, C, D, Custom). FMM grains
    have their regression maps built via openMotor (see srm_1d.fmm_grain)
    and attached as `GrainSegment.fmm_table`.

    Routes through ``build_snapped_geometry`` so cell boundaries align
    with segment edges and gap widths are guaranteed to be ≥1 cell.
    Inter-segment gaps default to ``max(3mm, 5%·D_outer)``.

    Parameters
    ----------
    ric_grains : list of dict
        From the .ric file's 'grains' key.
    target_propellant_cells : int
        Approximate number of cells to spend on propellant. Cell width
        is computed from this and the total propellant length, then
        used to integer-snap segment lengths and gaps.
    fmm_map_dim : int
        FMM regression-map resolution for FMM grain types. openMotor
        default is 1001. Higher = more accurate perimeter/port-area
        sampling but quadratically slower setup (skfmm.distance is
        O(mapDim²)). Sim-level config knob.

    Returns
    -------
    MotorGeometry
    """
    segments_spec = []
    D_outer = None

    for i, grain in enumerate(ric_grains):
        gtype = grain['type']
        props = grain['properties']
        inhibit_str = props.get('inhibitedEnds', 'Neither')
        if inhibit_str not in _INHIBIT_MAP:
            raise ValueError(
                f"Grain {i} has unknown inhibitedEnds '{inhibit_str}'. "
                f"Valid values: {list(_INHIBIT_MAP.keys())}"
            )
        inh_fwd, inh_aft = _INHIBIT_MAP[inhibit_str]
        seg_D_outer = props['diameter']

        taper_def = props.get('taper')
        bore_tapered = isinstance(taper_def, dict) and taper_def.get('enabled')
        # OD / end taper (independent of the bore taper): the casing diameter
        # shrinks over an end region (aft nozzle cone, fwd/aft dome). Applies
        # to FMM and analytic (BATES/Conical) grains alike.
        od_ends = []
        if isinstance(taper_def, dict):
            from .fmm_grain import od_ends_from_taper
            od_ends = od_ends_from_taper(taper_def)
        has_taper = bore_tapered or bool(od_ends)

        # A *bore* taper of a BATES/Conical is just a Conical — reject it (an
        # OD-only taper on them IS supported, analytically).
        if bore_tapered and gtype in ('BATES', 'Conical'):
            raise NotImplementedError(
                f"Grain {i}: the transient adapter supports BORE tapering only "
                f"for FMM cross-section grains (Finocyl, Star, ...), not "
                f"'{gtype}'. A linear bore taper of a BATES grain is a "
                f"Conical grain; use that instead. (An OD/end taper on "
                f"BATES/Conical is supported.)"
            )

        # A domed / coned end is bonded to the closure, so its end face does
        # not burn — force-inhibit it (mirrors the QS taper expander).
        eff_inh_fwd = inh_fwd or any(e.get('end') == 'fwd' for e in od_ends)
        eff_inh_aft = inh_aft or any(e.get('end') == 'aft' for e in od_ends)
        od_spec = (od_ends or None)

        if has_taper and gtype not in ('BATES', 'Conical'):
            # Axially-tapered FMM grain (bore and/or OD): build a TaperSpec and
            # route through build_snapped_geometry's per-cell taper path. The
            # transient solver resolves real per-station FMM tables (OD-clipped
            # where an OD taper is present); the discretization (station count)
            # is mesh-based and independent of the QS solver's L/D slicing.
            from .fmm_grain import taper_spec_from_props
            taper_spec = taper_spec_from_props(
                gtype, props, taper_def, map_dim=fmm_map_dim,
            )
            spec = {
                'D_bore_fwd': seg_D_outer,
                'D_bore_aft': seg_D_outer,
                'length': props['length'],
                'inhibit_fwd': eff_inh_fwd,
                'inhibit_aft': eff_inh_aft,
                'taper': taper_spec,
                'od_ends': od_spec,
            }
        elif gtype == 'BATES':
            spec = {
                'D_bore_fwd': props['coreDiameter'],
                'D_bore_aft': props['coreDiameter'],
                'length': props['length'],
                'inhibit_fwd': eff_inh_fwd,
                'inhibit_aft': eff_inh_aft,
                'od_ends': od_spec,
            }
        elif gtype == 'Conical':
            spec = {
                'D_bore_fwd': props['forwardCoreDiameter'],
                'D_bore_aft': props['aftCoreDiameter'],
                'length': props['length'],
                'inhibit_fwd': eff_inh_fwd,
                'inhibit_aft': eff_inh_aft,
                'od_ends': od_spec,
            }
        else:
            from .fmm_grain import from_ric_grain
            try:
                fmm_table = from_ric_grain(grain, map_dim=fmm_map_dim)
            except ValueError as e:
                raise ValueError(
                    f"Grain {i} has unsupported type '{gtype}'. "
                    f"BATES uses analytic; FMM types must be registered. "
                    f"Inner error: {e}"
                ) from e
            # FMM has no circular bore; D_bore_* are placeholders that
            # get overwritten by FmmTable in compile_geometry_arrays.
            spec = {
                'D_bore_fwd': seg_D_outer,
                'D_bore_aft': seg_D_outer,
                'length': props['length'],
                'inhibit_fwd': inh_fwd,
                'inhibit_aft': inh_aft,
                'fmm_table': fmm_table,
            }

        if D_outer is None:
            D_outer = seg_D_outer
        elif abs(seg_D_outer - D_outer) > 1e-9:
            warnings.warn(
                f"Grain {i} has diameter {seg_D_outer} != motor D_outer {D_outer}. "
                f"Using motor-level D_outer; per-segment outer-diameter is "
                f"not yet supported."
            )
        segments_spec.append(spec)

    if not segments_spec:
        raise ValueError("No supported grain segments found in .ric file.")

    # Default inter-segment gap: max(3mm, 5%·D_outer), but only when
    # at least one face at the interface is uninhibited. If both the
    # aft face of segment i and the forward face of segment i+1 are
    # inhibited, treat the interface as bonded/touching. This preserves
    # multi-slice .ric grains that use separate grain entries to describe
    # a continuous inhibited grain profile.
    inter_gap = max(0.003, D_outer * 0.05)
    for i, spec in enumerate(segments_spec[:-1]):
        next_spec = segments_spec[i + 1]
        interface_bonded = (
            spec.get('inhibit_aft', False)
            and next_spec.get('inhibit_fwd', False)
        )
        spec['gap_after'] = 0.0 if interface_bonded else inter_gap

    return build_snapped_geometry(
        segments_spec, D_outer,
        target_propellant_cells=target_propellant_cells,
    )


def convert_nozzle(ric_nozzle):
    """
    Convert an openMotor nozzle dict to srm_1d Nozzle.

    Field mapping (openMotor → srm_1d):
        throat        → D_throat              (m, same units)
        exit          → D_exit                (m, same units)
        efficiency    → efficiency            (—)
        divAngle      → div_angle             (deg, same)
        convAngle     → conv_angle            (deg, same)
        throatLength  → throat_length         (m, same)
        erosionCoeff  → erosion_coeff         m/(s·Pa) → μm/(s·MPa) (×1e12)
        slagCoeff     → slag_coeff            (m·Pa)/s → (m·MPa)/s (÷1e6)
    """
    # openMotor stores erosionCoeff in m/(s·Pa); our convention: μm/(s·MPa).
    # 1 m/(s·Pa) = 1e12 μm/(s·MPa).
    erosion_ours = ric_nozzle.get('erosionCoeff', 0.0) * 1e12
    # openMotor computes slag_rate = slagCoeff / pressure with pressure in Pa.
    # srm_1d uses pressure in MPa, so the coefficient must be divided by 1e6.
    slag_ours = ric_nozzle.get('slagCoeff', 0.0) / 1e6

    return Nozzle(
        D_throat=ric_nozzle['throat'],
        D_exit=ric_nozzle['exit'],
        efficiency=ric_nozzle.get('efficiency', 0.95),
        div_angle=ric_nozzle.get('divAngle', 15.0),
        conv_angle=ric_nozzle.get('convAngle', 30.0),
        throat_length=ric_nozzle.get('throatLength', 0.0),
        erosion_coeff=erosion_ours,
        slag_coeff=slag_ours,
    )


# ================================================================
# High-level convenience functions
# ================================================================

def ric_to_sim_args(motor, gas_props=None, target_propellant_cells=100,
                    **sim_overrides):
    """
    Convert a loaded .ric motor dict to run_simulation keyword arguments.

    Parameters
    ----------
    motor : dict
        Output from load_ric().
    gas_props : dict or None
        Transport properties: {'mu': Pa·s, 'k': W/(m·K), 'Cp': J/(kg·K)}.
    target_propellant_cells : int
        Approximate cell count to spend on propellant. Cell width is
        derived from this and used to integer-snap segment lengths
        and inter-segment gaps.
    **sim_overrides
        Additional keyword arguments passed to run_simulation
        (e.g. roughness, kappa, igniter params).

    Returns
    -------
    dict of keyword arguments for run_simulation. Includes 'geo',
    'propellant', 'nozzle', 'P_ambient', and 'P_cutoff' by default.
    """
    prop = convert_propellant(motor['propellant'], gas_props)
    geo = convert_geometry(motor['grains'],
                           target_propellant_cells=target_propellant_cells)
    nozzle = convert_nozzle(motor['nozzle'])

    P_amb = motor['config'].get('ambPressure', 101325.0)

    args = {
        'geo': geo,
        'propellant': prop,
        'nozzle': nozzle,
        'P_ambient': P_amb,
        'P_cutoff': P_amb * 5,
    }
    args.update(sim_overrides)

    return args


def run_from_ric(filepath, gas_props=None, transport_path=None,
                 pyrogen=None, pyrogen_mass=None,
                 pyrogen_throat_area=None, pyrogen_volume=None,
                 pyrogen_burn_area=None, pyrogen_burn_law='0d',
                 pyrogen_heat_flux_cal_cm2_s=None,
                 injection_topology='forward_plenum',
                 cartridge_length_m=-1.0,
                 basket_fill_fraction=0.5,
                 pellet_packing_fraction=0.60,
                 particle_diameter_m=None,
                 particle_LD_ratio=None,
                 T_ignition=756.0, k_solid=None,  # T_ign: v0.7.5 re-LHS (was 850)
                 radiation_emissivity=None,
                 flame_front_enabled=None, flame_front_velocity=None,
                 zn_enabled=None, kappa_zn=None,
                 verbose=True, **sim_overrides):
    """
    Load a .ric file, run the 1D simulation, compute performance.

    If ``gas_props`` is not given, looks for a sibling ``<stem>.transport.yaml``
    next to the .ric file. Falls back to estimated transport properties
    if neither is supplied.

    Parameters
    ----------
    filepath : str
        Path to the .ric file.
    gas_props : dict or None
        Explicit transport override: {'mu', 'k', 'Cp'}.
    transport_path : str or None
        Explicit transport YAML path. If None, sibling auto-resolution
        is attempted.
    pyrogen : Pyrogen, str, or None
        Explicit pyrogen object, built-in name, or YAML path. If None,
        a sibling ``<stem>.pyrogen.yaml`` must exist.
    verbose : bool
        If True, print simulation and performance summary blocks. Set
        False for large parameter sweeps.

    Returns
    -------
    result, perf, nozzle, geo, prop
    """
    motor = load_ric(filepath)

    stem, _ = os.path.splitext(filepath)

    # v0.8.0 Phase 3: gas transport now lives in the .ric propellant
    # (per-tab, selected by transportVariant). The .transport.yaml sidecar
    # is RETIRED — only an explicit transport_path / gas_props override is
    # honored; otherwise convert_propellant reads transport from the .ric
    # (and hard-faults on the unset sentinel rather than fabricating it).
    if gas_props is None and transport_path is not None:
        gas_props = load_transport(transport_path)

    args = ric_to_sim_args(
        motor, gas_props=gas_props,
        **sim_overrides,
    )

    geo = args.pop('geo')
    prop = args.pop('propellant')
    nozzle = args['nozzle']
    P_amb = args.get('P_ambient', 101325.0)

    # Propellant attribute overrides for LHS sweeps. None == use the
    # value from the .ric / propellant_overrides path; passing a number
    # mutates the propellant in place before the sim runs.
    if k_solid is not None:
        prop.k_solid = float(k_solid)
    if radiation_emissivity is not None:
        prop.radiation_emissivity = float(radiation_emissivity)
    # v0.7.4 Phase F: opt-in flame-spread front gate (default off in the
    # Propellant; a run script / sweep passes True to enable it).
    if flame_front_enabled is not None:
        prop.flame_front_enabled = bool(flame_front_enabled)
    if flame_front_velocity is not None:
        prop.flame_front_velocity = float(flame_front_velocity)
    # v0.7.4 Phase Z: opt-in Z-N dynamic burn-rate relaxation.
    if zn_enabled is not None:
        prop.zn_enabled = bool(zn_enabled)
    if kappa_zn is not None:
        prop.kappa_zn = float(kappa_zn)

    # v0.8.0 Phase 4: when no pyrogen kwarg is given, the motor file's
    # igniter block is authoritative (self-describing) — supplying both the
    # pyrogen material and the chamber sizing. An explicit pyrogen kwarg
    # takes the legacy path and uses the run_from_ric sizing kwargs.
    igniter_sizing = None
    if pyrogen is None:
        ig = load_igniter(motor)
        if ig is not None:
            pyrogen_obj, igniter_sizing = ig
        else:
            candidate = stem + '.pyrogen.yaml'
            if os.path.exists(candidate):
                pyrogen_obj = load_pyrogen(candidate)
            else:
                raise ValueError(
                    f"No pyrogen specified for {filepath} and the motor file "
                    "has no igniter block. Pass pyrogen='bpnv', "
                    "pyrogen=<Pyrogen>, migrate the motor "
                    "(migrate_all_motors), or add a sibling "
                    "<motor>.pyrogen.yaml."
                )
    elif isinstance(pyrogen, Pyrogen):
        pyrogen_obj = pyrogen
    else:
        pyrogen_obj = load_pyrogen(pyrogen)

    # Pyrogen attribute override for LHS sweeps (same pattern as
    # k_solid / radiation_emissivity above). The pyrogen-to-surface
    # sensible-power cap in Phase 3.5 changed effective heat delivery
    # noticeably; exposing this as a calibration knob lets the LHS
    # compensate explicitly instead of via degenerate routes (e.g.
    # oversized pyrogen volume).
    if pyrogen_heat_flux_cal_cm2_s is not None:
        pyrogen_obj.heat_flux_cal_cm2_s = float(pyrogen_heat_flux_cal_cm2_s)

    # v0.7.4 Phase C.1: per-run particle-geometry overrides. Lets a run
    # script tweak particle dimensions without editing the pyrogen YAML
    # (matches the user-flagged tunability requirement).
    if particle_diameter_m is not None:
        pyrogen_obj.particle_diameter_m = float(particle_diameter_m)
    if particle_LD_ratio is not None:
        pyrogen_obj.particle_LD_ratio = float(particle_LD_ratio)

    if igniter_sizing is not None:
        # Sizing from the motor file's igniter block (self-describing path).
        chamber_kwargs = igniter_sizing
    else:
        chamber_kwargs = dict(
            pyrogen_mass=pyrogen_mass,
            pyrogen_throat_area=pyrogen_throat_area,
            pyrogen_volume=pyrogen_volume,
            pyrogen_burn_area=pyrogen_burn_area,
            pyrogen_burn_law=pyrogen_burn_law,
            injection_topology=injection_topology,
            cartridge_length_m=cartridge_length_m,
            basket_fill_fraction=basket_fill_fraction,
            pellet_packing_fraction=pellet_packing_fraction,
        )
    args['pyrogen_chamber'] = build_pyrogen_chamber(
        pyrogen_obj, geo, nozzle, **chamber_kwargs)
    args['T_ignition'] = T_ignition
    args['verbose'] = verbose

    result = run_simulation(geo, prop, **args)

    perf = compute_motor_performance(result, nozzle, prop, P_ambient=P_amb)
    if verbose:
        print_performance_summary(perf, nozzle)

    return result, perf, nozzle, geo, prop


# ================================================================
# Per-grain metrics
# ================================================================

def compute_grain_metrics(result, geo, propellant):
    """
    Compute per-grain time histories from snapshot data.

    Returns regression depth, web remaining, and mass remaining for
    each grain at each snapshot time. Also computes Kn from the
    pressure trace using the equilibrium relation.

    Parameters
    ----------
    result : dict or SimulationChannels
        Output from run_simulation (must have snapshot 'D_port' field), or
        the channel-model equivalent. A dict is re-shaped via
        ``as_channels`` (pure, no recompute).
    geo : MotorGeometry
        The geometry used for the simulation.
    propellant : Propellant
        The propellant used.

    Returns
    -------
    dict with:
        'snap_times': ndarray (n_snaps,)
        'regression': ndarray (n_snaps, n_grains) — radial regression [m]
        'web': ndarray (n_snaps, n_grains) — web remaining [m]
        'grain_mass': ndarray (n_snaps, n_grains) — propellant mass [kg]
        'kn': ndarray (n_time,) — Kn from equilibrium at every time step
        'kn_times': ndarray (n_time,) — corresponding times
    """
    from .propellant import critical_flow_function, R_UNIVERSAL
    from .channels import as_channels

    sc = as_channels(result)
    dport_chan = sc.axial.get('D_port')
    ga = geo.compile_geometry_arrays()
    N_seg = ga['N_seg']
    dx = geo.dx

    n_snaps = dport_chan.n_frames if dport_chan is not None else 0
    snap_times = dport_chan.times if dport_chan is not None else np.zeros(0)
    dport_frames = dport_chan.getData() if dport_chan is not None else None
    regression = np.zeros((n_snaps, N_seg))
    web = np.zeros((n_snaps, N_seg))
    grain_mass = np.zeros((n_snaps, N_seg))

    for s in range(n_snaps):
        D_port = dport_frames[s]
        for k in range(N_seg):
            # Cells belonging to this grain
            mask = ga['cell_segment_id'] == k
            if not np.any(mask):
                continue

            D_bore_init = np.mean(ga['cell_D_bore_init'][mask])
            D_cells = D_port[mask]
            avg_D = np.mean(D_cells)

            regression[s, k] = (avg_D - D_bore_init) / 2.0
            web[s, k] = max(0.0, (geo.D_outer - avg_D) / 2.0)

            # Mass remaining: annular volume × density
            # Use per-cell D_port for accuracy
            cell_volume = np.sum(
                np.pi / 4.0 * (geo.D_outer**2 - D_cells**2) * dx
            )
            grain_mass[s, k] = cell_volume * propellant.rho_propellant

    # Kn from equilibrium: P^(1-n) = rho_p * a * Kn * c*
    # → Kn = P^(1-n) / (rho_p * a * c*)
    rep_tab = propellant.representative_tab()
    R_spec = R_UNIVERSAL / rep_tab.molecular_weight
    Gamma = critical_flow_function(rep_tab.gamma)
    c_star = np.sqrt(R_spec * rep_tab.T_flame) / Gamma

    t_arr = sc['time']
    P_arr = sc['P_head']
    # NOTE: Kn equilibrium uses representative-tab a/n. For multi-tab
    # propellants, pressure-dependent Kn would require per-sample tab
    # lookup — left as a follow-up since most exports use a single tab.
    kn_denom = propellant.rho_propellant * rep_tab.a * c_star
    kn = np.where(
        P_arr > 1e4,
        np.power(P_arr, 1.0 - rep_tab.n) / kn_denom,
        0.0,
    )

    snap_times = np.asarray(snap_times)

    return {
        'snap_times': snap_times,
        'regression': regression,
        'web': web,
        'grain_mass': grain_mass,
        'kn': kn,
        'kn_times': t_arr,
    }


# ================================================================
# CSV export (openMotor-compatible)
# ================================================================

def result_to_csv(result, perf=None, geo=None, propellant=None,
                  dt_sample=None, separator=','):
    """
    Export simulation results as openMotor-compatible CSV.

    Parameters
    ----------
    result : dict or SimulationChannels
        Output from run_simulation, or the channel-model equivalent.
    perf : dict or None
        Output from compute_motor_performance.
    geo : MotorGeometry or None
        If provided (with propellant), per-grain columns are included.
    propellant : Propellant or None
        Required alongside geo for grain metrics.
    dt_sample : float or None
        Resample to this timestep [s] for manageable file size.
        If None, uses 0.001s (1 kHz). openMotor default is 0.025s.
    separator : str
        Column separator.

    Returns
    -------
    str : CSV content.
    """
    from .channels import as_channels
    sc = as_channels(result)
    t = sc['time']
    P = sc['P_head']
    P_exit = sc['P_exit']
    D_throat = sc['D_throat']

    # Compute Kn and grain metrics if geometry is available
    grain_met = None
    if geo is not None and propellant is not None:
        grain_met = compute_grain_metrics(sc, geo, propellant)

    # Downsample
    if dt_sample is None:
        dt_sample = 0.001
    if dt_sample > 0 and len(t) > 1:
        dt_actual = t[1] - t[0]
        stride = max(1, int(dt_sample / dt_actual))
    else:
        stride = 1
    indices = np.arange(0, len(t), stride)

    # Build headers
    headers = ['Time (s)', 'Kn', 'Pressure (Pa)', 'Exit Pressure (Pa)',
               'dThroat (m)']
    if perf is not None:
        headers.extend(['Force (N)', 'Isp (s)'])
    if grain_met is not None:
        n_grains = grain_met['regression'].shape[1]
        for k in range(n_grains):
            headers.append(f'Regression G{k} (m)')
        for k in range(n_grains):
            headers.append(f'Web G{k} (m)')

    lines = [separator.join(headers)]

    # Interpolate per-grain data to time history if available
    reg_interp = None
    web_interp = None
    if grain_met is not None and len(grain_met['snap_times']) > 1:
        n_grains = grain_met['regression'].shape[1]
        reg_interp = np.zeros((len(t), n_grains))
        web_interp = np.zeros((len(t), n_grains))
        for k in range(n_grains):
            reg_interp[:, k] = np.interp(
                t, grain_met['snap_times'], grain_met['regression'][:, k]
            )
            web_interp[:, k] = np.interp(
                t, grain_met['snap_times'], grain_met['web'][:, k]
            )

    for i in indices:
        # Kn from equilibrium
        kn_val = grain_met['kn'][i] if grain_met is not None else 0.0

        row = [
            f"{t[i]:.6f}",
            f"{kn_val:.2f}",
            f"{P[i]:.2f}",
            f"{P_exit[i]:.2f}",
            f"{D_throat[i] - D_throat[0]:.8f}",
        ]
        if perf is not None:
            row.append(f"{perf['thrust'][i]:.2f}")
            row.append(f"{perf['Isp'][i]:.2f}")
        if reg_interp is not None:
            for k in range(n_grains):
                row.append(f"{reg_interp[i, k]:.8f}")
            for k in range(n_grains):
                row.append(f"{web_interp[i, k]:.8f}")
        lines.append(separator.join(row))

    return '\n'.join(lines)


def save_csv(filepath, result, perf=None, geo=None, propellant=None,
             dt_sample=None, separator=','):
    """Write result_to_csv output to a file."""
    csv = result_to_csv(result, perf, geo, propellant, dt_sample, separator)
    with open(filepath, 'w') as f:
        f.write(csv)
    n_lines = csv.count('\n')
    print(f"Saved {filepath} ({n_lines} rows, dt={dt_sample or 0.001:.3f}s)")


# ================================================================
# openMotor CSV loader (for comparison)
# ================================================================

def load_openmotor_csv(filepath, separator=','):
    """
    Load an openMotor CSV export for comparison plotting.

    Parameters
    ----------
    filepath : str
        Path to the CSV file exported from openMotor.

    Returns
    -------
    dict with numpy arrays keyed by column name (lowercase, no units).
    Always includes 'time' and 'pressure'. Other keys depend on what
    openMotor exported.
    """
    with open(filepath, 'r') as f:
        header_line = f.readline().strip()

    headers_raw = [h.strip() for h in header_line.split(separator)]

    # Normalize header names: "Time (s)" → "time", "Pressure (Pa)" → "pressure"
    def _normalize(h):
        h = h.split('(')[0].strip().lower()
        h = h.replace(' ', '_')
        return h

    col_names = [_normalize(h) for h in headers_raw]
    data = np.loadtxt(filepath, delimiter=separator, skiprows=1)

    result = {}
    for i, name in enumerate(col_names):
        if i < data.shape[1]:
            result[name] = data[:, i]

    return result


# ================================================================
# Summary / comparison helper
# ================================================================

def print_ric_summary(filepath):
    """Print a human-readable summary of a .ric file."""
    motor = load_ric(filepath)

    grains = motor['grains']
    noz = motor['nozzle']
    prop = motor['propellant']

    print(f"openMotor file: {filepath}")
    print(f"  Propellant: {prop.get('name', 'unnamed')}")
    print(f"    density = {prop['density']:.0f} kg/m³")
    for i, tab in enumerate(prop['tabs']):
        print(f"    tab {i}: a={tab['a']:.4e}  n={tab['n']:.2f}  "
              f"γ={tab['k']:.2f}  MW={tab['m']:.1f}g/mol  "
              f"T_flame={tab['t']:.0f}K  "
              f"P=[{tab.get('minPressure',0)/1e6:.1f}-{tab.get('maxPressure',0)/1e6:.1f}]MPa")

    print(f"  Grains: {len(grains)}")
    for i, g in enumerate(grains):
        p = g['properties']
        if g['type'] == 'Conical':
            bore_str = (f"D_bore_fwd={p['forwardCoreDiameter']*1e3:.1f}mm  "
                        f"D_bore_aft={p['aftCoreDiameter']*1e3:.1f}mm")
        elif 'coreDiameter' in p:
            bore_str = f"D_bore={p['coreDiameter']*1e3:.1f}mm"
        else:
            bore_str = "(FMM core)"
        print(f"    [{i}] {g['type']}: {bore_str}  "
              f"D_outer={p['diameter']*1e3:.1f}mm  L={p['length']*1e3:.1f}mm  "
              f"inhibited={p.get('inhibitedEnds', 'Neither')}")

    print(f"  Nozzle: throat={noz['throat']*1e3:.1f}mm  "
          f"exit={noz['exit']*1e3:.1f}mm  "
          f"div={noz.get('divAngle',15):.0f}°  "
          f"efficiency={noz.get('efficiency',0.95):.3f}")

    erosion = noz.get('erosionCoeff', 0)
    slag = noz.get('slagCoeff', 0)
    if erosion > 0 or slag > 0:
        print(f"    erosion={erosion:.2e} m/(s·Pa) "
              f"[= {erosion*1e12:.1f} μm/(s·MPa)]  "
              f"slag={slag:.2e} (m·MPa)/s")
    print()
