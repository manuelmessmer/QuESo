"""Validate and normalize public multi-component QuESo settings."""

import re
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _require_mapping(value: Any, location: str) -> Mapping[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{location} must be a JSON object.")
    return value


def _reject_unknown_fields(
    value: Mapping[str, Any], allowed: set[str], location: str
) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValueError(f"{location} contains unknown fields: {unknown}.")


def _validate_identifier(value: Any, location: str) -> str:
    if not isinstance(value, str) or not IDENTIFIER_PATTERN.fullmatch(value):
        raise ValueError(f"{location} must match [A-Za-z_][A-Za-z0-9_]*.")
    return value


def _resolve_path(value: Any, settings_directory: Path, location: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{location} must be a non-empty path string.")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = settings_directory / path
    return str(path.resolve())


def _validate_global_settings(
    raw_settings: Mapping[str, Any], settings_directory: Path
) -> dict[str, Any]:
    global_settings = {
        "output_directory_name": "queso_output",
        "echo_level": 1,
        "write_output_to_file": True,
    }
    supplied = raw_settings.get("global_settings", {})
    supplied_mapping = _require_mapping(supplied, "global_settings")
    _reject_unknown_fields(supplied_mapping, set(global_settings), "global_settings")
    global_settings.update(supplied_mapping)

    if not isinstance(global_settings["echo_level"], int) or isinstance(
        global_settings["echo_level"], bool
    ):
        raise ValueError("global_settings.echo_level must be an integer.")
    if not isinstance(global_settings["write_output_to_file"], bool):
        raise ValueError("global_settings.write_output_to_file must be a boolean.")
    global_settings["output_directory_name"] = _resolve_path(
        global_settings["output_directory_name"],
        settings_directory,
        "global_settings.output_directory_name",
    )
    return global_settings


def _validate_components(raw_components: Any) -> dict[str, dict[str, Any]]:
    if not isinstance(raw_components, list) or not raw_components:
        raise ValueError("components must be a non-empty JSON array.")

    components: dict[str, dict[str, Any]] = {}
    for index, raw_component in enumerate(raw_components):
        location = f"components[{index}]"
        component = _require_mapping(raw_component, location)
        for field in ("component_name", "input_filename"):
            if field not in component:
                raise ValueError(f"{location} must define '{field}'.")

        component_name = _validate_identifier(
            component["component_name"], f"{location}.component_name"
        )
        if component_name in components:
            raise ValueError(f"Duplicate component name '{component_name}'.")
        components[component_name] = deepcopy(dict(component))

    for component_name, component in components.items():
        if "background_grid_settings" not in component:
            continue
        location = f"component '{component_name}'.background_grid_settings"
        grid = _require_mapping(component["background_grid_settings"], location)
        if "from_other_component" in grid:
            if set(grid) != {"from_other_component"}:
                raise ValueError(
                    f"{location} must contain only 'from_other_component' "
                    "when aliasing."
                )
            source_name = _validate_identifier(
                grid["from_other_component"], f"{location}.from_other_component"
            )
            if source_name == component_name:
                raise ValueError(f"Component '{component_name}' cannot alias itself.")
            if source_name not in components:
                raise ValueError(
                    f"Component '{component_name}' aliases unknown component "
                    f"'{source_name}'."
                )
            source_grid = components[source_name].get("background_grid_settings")
            if not isinstance(source_grid, dict):
                raise ValueError(
                    f"Component '{component_name}' aliases component "
                    f"'{source_name}' without concrete grid settings."
                )
            if "from_other_component" in source_grid:
                raise ValueError("Chained background-grid aliases are not allowed.")
    return components


def _validate_conditions(
    raw_conditions: Any, component_names: set[str], settings_directory: Path
) -> dict[str, list[dict[str, Any]]]:
    if raw_conditions is None:
        raw_conditions = []
    if not isinstance(raw_conditions, list):
        raise ValueError("conditions must be a JSON array.")

    conditions_by_component = {name: [] for name in component_names}
    condition_ids: set[int] = set()
    for index, raw_condition in enumerate(raw_conditions):
        location = f"conditions[{index}]"
        condition = _require_mapping(raw_condition, location)
        for field in ("condition_id", "component_name", "input_filename"):
            if field not in condition:
                raise ValueError(f"{location} must define '{field}'.")

        condition_id = condition["condition_id"]
        if (
            not isinstance(condition_id, int)
            or isinstance(condition_id, bool)
            or condition_id < 1
        ):
            raise ValueError(f"{location}.condition_id must be a positive integer.")
        if condition_id in condition_ids:
            raise ValueError(f"Duplicate condition_id {condition_id}.")
        condition_ids.add(condition_id)

        component_name = _validate_identifier(
            condition["component_name"], f"{location}.component_name"
        )
        if component_name not in component_names:
            raise ValueError(
                f"{location} references unknown component '{component_name}'."
            )

        normalized = deepcopy(dict(condition))
        normalized.pop("component_name")
        normalized["input_filename"] = _resolve_path(
            normalized["input_filename"],
            settings_directory,
            f"{location}.input_filename",
        )
        if "coupling_partner" in normalized:
            partner = _validate_identifier(
                normalized["coupling_partner"], f"{location}.coupling_partner"
            )
            if partner not in component_names:
                raise ValueError(
                    f"{location} references unknown coupling partner '{partner}'."
                )
            if partner == component_name:
                raise ValueError(f"{location} cannot couple a component to itself.")
            normalized.setdefault("slip", False)
        conditions_by_component[component_name].append(normalized)
    return conditions_by_component


def parse_settings_by_component(
    raw_settings: Mapping[str, Any], settings_path: Path
) -> dict[str, dict[str, Any]]:
    """Validate public settings and return flattened settings by component.

    Args:
        raw_settings: Parsed public QuESo settings.
        settings_path: Path of the JSON file that declared the settings.

    Returns:
        Normalized C++ settings dictionaries keyed by component name.

    Raises:
        ValueError: If the public schema is invalid.
    """
    root = _require_mapping(raw_settings, "QuESo settings")
    _reject_unknown_fields(
        root, {"global_settings", "components", "conditions"}, "QuESo settings"
    )
    if "components" not in root:
        raise ValueError(
            "QuESo settings must use the multi-component schema with 'components'."
        )

    settings_path = Path(settings_path).expanduser().resolve()
    settings_directory = settings_path.parent
    global_settings = _validate_global_settings(root, settings_directory)
    components = _validate_components(root["components"])
    conditions = _validate_conditions(
        root.get("conditions", []), set(components), settings_directory
    )

    normalized_components: dict[str, dict[str, Any]] = {}
    for component_name, source_component in components.items():
        component = deepcopy(source_component)
        component["input_filename"] = _resolve_path(
            component["input_filename"],
            settings_directory,
            f"component '{component_name}'.input_filename",
        )
        grid = component.get("background_grid_settings")
        if isinstance(grid, dict) and "from_other_component" in grid:
            source_name = grid["from_other_component"]
            component["background_grid_settings"] = deepcopy(
                components[source_name]["background_grid_settings"]
            )

        component["echo_level"] = global_settings["echo_level"]
        component["write_output_to_file"] = global_settings["write_output_to_file"]
        component["output_directory_name"] = str(
            Path(global_settings["output_directory_name"]) / component_name
        )
        if conditions[component_name]:
            component["conditions_settings_list"] = conditions[component_name]
        normalized_components[component_name] = component

    return normalized_components
