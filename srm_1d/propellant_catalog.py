"""Read an openMotor propellant catalog without changing solver physics.

The bundled catalog is reference data. Missing gas transport values are
reported, never estimated or replaced with defaults.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path

import yaml


DEFAULT_CATALOG = Path(__file__).resolve().parent / "catalogs" / "propellants.yaml"


class _CatalogLoader(yaml.SafeLoader):
    """Recognize only the two openMotor YAML metadata tags."""


def _sequence(loader, node):
    return loader.construct_sequence(node)


_CatalogLoader.add_constructor("tag:yaml.org,2002:python/tuple", _sequence)
_CatalogLoader.add_constructor(
    "tag:yaml.org,2002:python/object/apply:uilib.fileIO.fileTypes",
    _sequence,
)

_REQUIRED_TAB_FIELDS = frozenset(
    {"a", "n", "k", "m", "t", "minPressure", "maxPressure"}
)


def load_propellant_catalog(path=None):
    """Return an ordered dict mapping names to raw openMotor propellant data.

    Defaults to the bundled catalog. Duplicate names or malformed records
    fail loudly rather than silently overwriting entries.
    """
    source = Path(path) if path is not None else DEFAULT_CATALOG
    with source.open("r", encoding="utf-8") as stream:
        document = yaml.load(stream, Loader=_CatalogLoader)
    if not isinstance(document, dict) or not isinstance(document.get("data"), list):
        raise ValueError(f"Invalid propellant catalog {source}: expected a 'data' list")

    by_name = {}
    for i, record in enumerate(document["data"]):
        if not isinstance(record, dict):
            raise ValueError(f"Catalog entry {i} must be a mapping")
        name = record.get("name")
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"Catalog entry {i} is missing a name")
        if name in by_name:
            raise ValueError(f"Duplicate propellant name: {name!r}")
        if "density" not in record or not isinstance(record.get("tabs"), list) or not record["tabs"]:
            raise ValueError(f"Propellant {name!r} must provide density and nonempty tabs")
        for j, tab in enumerate(record["tabs"]):
            if not isinstance(tab, dict) or (_REQUIRED_TAB_FIELDS - tab.keys()):
                raise ValueError(f"Propellant {name!r} tab {j} lacks openMotor combustion fields")
        by_name[name] = record
    return by_name


def list_propellants(path=None):
    """List names in catalog order."""
    return list(load_propellant_catalog(path))


def get_propellant(name, path=None):
    """Return a detached copy of the matching raw openMotor definition."""
    catalog = load_propellant_catalog(path)
    try:
        return deepcopy(catalog[name])
    except KeyError:
        from difflib import get_close_matches
        close = get_close_matches(name, catalog.keys(), n=3)
        hint = f" Closest names: {', '.join(close)}" if close else ""
        raise KeyError(f"Propellant {name!r} not found.{hint}") from None


def missing_transport_fields(propellant):
    """Report zero/missing gas transport fields for the active variant.

    Returns {zero_based_tab_index: [missing_field_names]}.
    Values are never filled or inferred.
    """
    variant = propellant.get("transportVariant", "frozen")
    if variant not in ("frozen", "effective"):
        raise ValueError(f"Unrecognized transportVariant: {variant!r}")
    keys = ("mu", "kThermalFrozen", "cpFrozen") if variant == "frozen" else (
        "mu", "kThermalEffective", "cpEffective"
    )
    missing = {}
    for index, tab in enumerate(propellant["tabs"]):
        absent = []
        for key in keys:
            try:
                valid = float(tab.get(key, 0)) > 0.0
            except (TypeError, ValueError):
                valid = False
            if not valid:
                absent.append(key)
        if absent:
            missing[index] = absent
    return missing


def with_catalog_propellant(motor_dict, name, path=None):
    """Copy a flat load_ric motor dict and replace only its propellant.

    Compatible with the existing ric_to_sim_args adapter; that adapter
    still requires physically sourced nonzero transport before simulation.
    """
    if not isinstance(motor_dict, dict) or "propellant" not in motor_dict:
        raise ValueError("Expected a flat motor dict from load_ric")
    copied = deepcopy(motor_dict)
    copied["propellant"] = get_propellant(name, path)
    return copied


def main(argv=None):
    parser = argparse.ArgumentParser(description="Inspect an openMotor propellant catalog")
    parser.add_argument("--catalog", type=Path, help="Path to another library YAML")
    parser.add_argument("--show", metavar="NAME", help="Inspect an exact propellant name")
    args = parser.parse_args(argv)
    if args.show:
        propellant = get_propellant(args.show, args.catalog)
        missing = missing_transport_fields(propellant)
        print(f"Name: {propellant['name']}")
        print(f"Pressure tabs: {len(propellant['tabs'])}")
        print(f"Transport inputs complete: {not bool(missing)}")
        if missing:
            print(f"Missing transport fields by zero-based tab: {missing}")
    else:
        for name in list_propellants(args.catalog):
            print(name)


if __name__ == "__main__":
    main()
