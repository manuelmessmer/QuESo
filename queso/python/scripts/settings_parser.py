"""Helpers for parsing and normalizing raw QuESo settings."""

from copy import deepcopy
from pathlib import Path
from typing import Any

REQUIRED_COMPONENT_FIELDS: dict[str, type] = {
    "component_name": str,
    "background_grid_settings": dict,
}

OPTIONAL_COMPONENT_FIELDS: dict[str, type] = {
    "input_filename": str,
    "trimmed_quadrature_rule_settings": dict,
    "non_trimmed_quadrature_rule_settings": dict,
    "conditions_settings_list": list,
}


def _resolve_path(filename: str, settings_dir: Path) -> str:
    """Resolve a settings path against the current environment and settings dir."""
    path = Path(filename)
    if path.is_absolute():
        return str(path)

    # Preserve repo-root-relative paths like "queso/tests/..." when they already
    # exist from the current working directory used by the test and example entrypoints.
    if path.exists():
        return str(path.resolve())

    return str((settings_dir / path).resolve())


def _resolve_input_filename(settings: dict[str, Any], settings_dir: Path) -> None:
    """Resolve ``input_filename`` in-place when the setting is present."""
    input_filename = settings.get("input_filename")
    if input_filename is not None:
        settings["input_filename"] = _resolve_path(str(input_filename), settings_dir)


def _set_component_output_directory(
    settings: dict[str, Any], component_name: str
) -> None:
    """Rewrite the output directory to a component-specific subdirectory."""
    output_directory_name = str(
        settings.get("output_directory_name", "queso_output")
    ).strip()
    settings["output_directory_name"] = str(
        Path(output_directory_name) / component_name
    )


def _validate_and_extract_global_settings(
    raw_settings: dict[str, Any],
) -> dict[str, Any]:
    """Validate the top-level settings wrapper and return normalized globals."""
    if "components" not in raw_settings:
        return _legacy_settings_to_global_settings(raw_settings)

    global_settings = raw_settings.get("global_settings", {})
    if not isinstance(global_settings, dict):
        raise ValueError("QuESoSettings field 'global_settings' must be a JSON object.")

    components = raw_settings.get("components")
    if not isinstance(components, list) or not components:
        raise ValueError("QuESoSettings must contain a non-empty 'components' list.")

    component_names: set[str] = set()
    component_map: dict[str, dict[str, Any]] = {}
    for component in components:
        if not isinstance(component, dict):
            raise ValueError("Each component must be a JSON object.")

        for field, expected_type in REQUIRED_COMPONENT_FIELDS.items():
            if field not in component:
                raise ValueError(f"Each component must define '{field}'.")
            if not isinstance(component[field], expected_type):
                raise ValueError(
                    f"Component field '{field}' must be {expected_type.__name__}."
                )

        for field, expected_type in OPTIONAL_COMPONENT_FIELDS.items():
            if field in component and not isinstance(component[field], expected_type):
                raise ValueError(
                    f"Component field '{field}' must be {expected_type.__name__}."
                )

        component_name = str(component["component_name"]).strip()
        if not component_name:
            raise ValueError(
                "Each component must define a non-empty string 'component_name'."
            )
        if component_name in component_names:
            raise ValueError(f"Duplicate component name '{component_name}'.")

        component_names.add(component_name)
        component_map[component_name] = component

    _validate_background_grid_aliases(component_map)
    return global_settings


def _legacy_settings_to_global_settings(raw_settings: dict[str, Any]) -> dict[str, Any]:
    """Extract global settings from the legacy single-component schema."""
    general_settings = raw_settings.get("general_settings")
    if not isinstance(general_settings, dict):
        raise ValueError(
            "QuESoSettings must define either 'components' or legacy 'general_settings'."
        )
    return {
        "output_directory_name": general_settings.get(
            "output_directory_name", "queso_output"
        ),
        "echo_level": general_settings.get("echo_level", 1),
        "write_output_to_file": general_settings.get("write_output_to_file", True),
    }


def _validate_background_grid_aliases(component_map: dict[str, dict[str, Any]]) -> None:
    """Validate cross-component background grid alias references."""
    for component_name, component in component_map.items():
        background_grid_settings = component["background_grid_settings"]
        alias_name = _get_grid_alias_name(background_grid_settings)
        if alias_name is None:
            continue

        if alias_name == component_name:
            raise ValueError(
                f"Component '{component_name}' cannot reference itself as background grid source."
            )
        if alias_name not in component_map:
            raise ValueError(
                f"Component '{component_name}' references unknown background grid component '{alias_name}'."
            )

        referenced_grid_settings = component_map[alias_name]["background_grid_settings"]
        if _get_grid_alias_name(referenced_grid_settings) is not None:
            raise ValueError(
                f"Component '{component_name}' references '{alias_name}', but chained background grid aliases are not allowed."
            )


def _get_grid_alias_name(background_grid_settings: dict[str, Any]) -> str | None:
    """Return the referenced component name for an alias-only grid definition."""
    alias_name = background_grid_settings.get("from_other_component")
    if alias_name is None:
        return None

    if set(background_grid_settings.keys()) != {"from_other_component"}:
        raise ValueError(
            "background_grid_settings must define either a full grid or exactly {'from_other_component': '<name>'}."
        )

    if not isinstance(alias_name, str) or not alias_name.strip():
        raise ValueError(
            "background_grid_settings.from_other_component must be a non-empty string."
        )

    return alias_name.strip()


def _resolve_background_grid_settings(
    component_name: str,
    component_map: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Return the concrete background grid settings for a component."""
    background_grid_settings = deepcopy(
        component_map[component_name]["background_grid_settings"]
    )
    alias_name = _get_grid_alias_name(background_grid_settings)
    if alias_name is None:
        return background_grid_settings
    return deepcopy(component_map[alias_name]["background_grid_settings"])


def parse_settings_by_component(
    raw_settings: dict[str, Any],
    settings_dir: Path,
) -> dict[str, dict[str, Any]]:
    """Normalize raw QuESo settings into a per-component settings map."""
    global_settings = deepcopy(_validate_and_extract_global_settings(raw_settings))
    components = raw_settings.get("components")
    if components is None:
        components = [
            {
                "component_name": "main",
                "background_grid_settings": raw_settings["background_grid_settings"],
                "trimmed_quadrature_rule_settings": raw_settings.get(
                    "trimmed_quadrature_rule_settings", {}
                ),
                "non_trimmed_quadrature_rule_settings": raw_settings.get(
                    "non_trimmed_quadrature_rule_settings", {}
                ),
                "conditions_settings_list": raw_settings.get(
                    "conditions_settings_list", []
                ),
            }
        ]
        input_filename = raw_settings["general_settings"].get("input_filename")
        if input_filename is not None:
            components[0]["input_filename"] = input_filename
    component_map = {
        str(component["component_name"]).strip(): component for component in components
    }

    component_settings: dict[str, dict[str, Any]] = {}
    for component_name, component in component_map.items():
        settings = deepcopy(component)
        settings["background_grid_settings"] = _resolve_background_grid_settings(
            component_name, component_map
        )
        settings.update(deepcopy(global_settings))
        _set_component_output_directory(settings, component_name)

        _resolve_input_filename(settings, settings_dir)

        conditions = settings.get("conditions_settings_list", [])
        if not isinstance(conditions, list):
            raise ValueError(
                f"Component '{component_name}' field 'conditions_settings_list' must be a list."
            )
        if not conditions:
            settings.pop("conditions_settings_list", None)
        else:
            for condition in conditions:
                if not isinstance(condition, dict):
                    raise ValueError(
                        f"Component '{component_name}' contains a condition that is not an object."
                    )
                if "condition_id" not in condition:
                    raise ValueError(
                        f"Component '{component_name}' has a condition without 'condition_id'."
                    )
                _resolve_input_filename(condition, settings_dir)

        component_settings[component_name] = settings

    return component_settings
