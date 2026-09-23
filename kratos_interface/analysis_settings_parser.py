"""Validation and normalization of QuESo-to-Kratos analysis settings."""

import re
from collections.abc import Mapping
from typing import Any

IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _identifier(value: Any, location: str) -> str:
    if not isinstance(value, str) or not IDENTIFIER_PATTERN.fullmatch(value):
        raise ValueError(f"{location} must be a valid Kratos identifier.")
    return value


def _model_part_path(value: Any, location: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{location} must be a non-empty model-part path.")
    for segment in value.split("."):
        _identifier(segment, location)
    return value


def parse_analysis_settings(
    raw_settings: Mapping[str, Any], queso_model: Any
) -> dict[str, Any]:
    """Validate analysis settings and resolve component-to-geometry ownership.

    Args:
        raw_settings: Parsed AnalysisSettings JSON object.
        queso_model: Configured and created pyqueso model.

    Returns:
        Normalized analysis settings.

    Raises:
        ValueError: If the settings violate the analysis contract.
    """
    if not isinstance(raw_settings, dict) or set(raw_settings) != {"model_parts"}:
        raise ValueError("AnalysisSettings must contain only 'model_parts'.")
    model_parts = raw_settings["model_parts"]
    if not isinstance(model_parts, list) or not model_parts:
        raise ValueError("AnalysisSettings.model_parts must be a non-empty array.")

    configured_components = set(queso_model.component_names)
    _validate_kratos_conditions(queso_model)
    assigned_components: set[str] = set()
    geometry_names: set[str] = set()
    normalized_geometries: dict[str, dict[str, Any]] = {}
    root_name = ""

    for model_part_index, raw_model_part in enumerate(model_parts):
        location = f"model_parts[{model_part_index}]"
        if not isinstance(raw_model_part, dict) or set(raw_model_part) != {
            "name",
            "geometries",
        }:
            raise ValueError(f"{location} must define only 'name' and 'geometries'.")
        model_part_name = _model_part_path(raw_model_part["name"], f"{location}.name")
        current_root = model_part_name.split(".")[0]
        if root_name and current_root != root_name:
            raise ValueError("All analysis model parts must share one root.")
        root_name = current_root

        geometries = raw_model_part["geometries"]
        if not isinstance(geometries, list) or not geometries:
            raise ValueError(f"{location}.geometries must be a non-empty array.")
        for geometry_index, raw_geometry in enumerate(geometries):
            geometry_location = f"{location}.geometries[{geometry_index}]"
            if not isinstance(raw_geometry, dict) or set(raw_geometry) != {
                "name",
                "queso_components",
            }:
                raise ValueError(
                    f"{geometry_location} must define only 'name' and "
                    "'queso_components'."
                )
            geometry_name = _identifier(
                raw_geometry["name"], f"{geometry_location}.name"
            )
            if geometry_name in geometry_names:
                raise ValueError(f"Duplicate analysis geometry '{geometry_name}'.")
            geometry_names.add(geometry_name)

            raw_components = raw_geometry["queso_components"]
            if not isinstance(raw_components, list) or not raw_components:
                raise ValueError(
                    f"{geometry_location}.queso_components must be a non-empty array."
                )
            components: list[dict[str, Any]] = []
            reference_grid: Any = None
            for component_index, raw_component in enumerate(raw_components):
                component_location = (
                    f"{geometry_location}.queso_components[{component_index}]"
                )
                if not isinstance(raw_component, dict) or set(raw_component) != {
                    "component_name",
                    "element_settings",
                }:
                    raise ValueError(
                        f"{component_location} must define only 'component_name' "
                        "and 'element_settings'."
                    )
                component_name = _identifier(
                    raw_component["component_name"],
                    f"{component_location}.component_name",
                )
                if component_name not in configured_components:
                    raise ValueError(
                        f"{component_location} references unknown component "
                        f"'{component_name}'."
                    )
                if component_name in assigned_components:
                    raise ValueError(
                        f"Component '{component_name}' is assigned to more than "
                        "one geometry."
                    )
                assigned_components.add(component_name)

                element_settings = raw_component["element_settings"]
                if not isinstance(element_settings, dict):
                    raise ValueError(
                        f"{component_location}.element_settings must be an object."
                    )
                if not set(element_settings).issubset({"property_id", "element_type"}):
                    raise ValueError(
                        f"{component_location}.element_settings contains "
                        "unknown fields."
                    )
                property_id = element_settings.get("property_id")
                if (
                    not isinstance(property_id, int)
                    or isinstance(property_id, bool)
                    or property_id < 1
                ):
                    raise ValueError(
                        f"{component_location}.element_settings.property_id "
                        "must be positive."
                    )
                element_type = element_settings.get(
                    "element_type", "UpdatedLagrangianElement3D8N"
                )
                if not isinstance(element_type, str) or not element_type:
                    raise ValueError(
                        f"{component_location}.element_settings.element_type "
                        "must be non-empty."
                    )

                grid = queso_model.settings(component_name)["background_grid_settings"]
                if reference_grid is None:
                    reference_grid = grid
                elif grid != reference_grid:
                    raise ValueError(
                        f"Geometry '{geometry_name}' contains components with "
                        "different grids."
                    )
                components.append(
                    {
                        "component_name": component_name,
                        "element_settings": {
                            "property_id": property_id,
                            "element_type": element_type,
                        },
                    }
                )

            normalized_geometries[geometry_name] = {
                "model_part_name": model_part_name,
                "background_grid_settings": reference_grid,
                "queso_components": components,
            }

    missing_components = configured_components - assigned_components
    if missing_components:
        raise ValueError(
            "Components are missing analysis geometry assignments: "
            f"{sorted(missing_components)}."
        )

    component_geometry = {
        component["component_name"]: geometry_name
        for geometry_name, geometry in normalized_geometries.items()
        for component in geometry["queso_components"]
    }
    for component_name in queso_model.component_names:
        for condition in queso_model.settings(component_name).get(
            "conditions_settings_list", []
        ):
            partner = condition.get("coupling_partner")
            if (
                partner is not None
                and component_geometry[component_name] == component_geometry[partner]
            ):
                raise ValueError(
                    f"Coupling condition {condition['condition_id']} must connect "
                    "distinct analysis geometries."
                )

    return {
        "root_model_part_name": root_name,
        "geometries": normalized_geometries,
        "component_geometry": component_geometry,
    }


def _validate_kratos_conditions(queso_model: Any) -> None:
    """Validate only the condition fields interpreted by the Kratos interface."""
    required_fields = {
        "PenaltySupportCondition": ("value", "penalty_factor"),
        "LagrangeSupportCondition": ("value",),
        "SurfaceLoadCondition": ("modulus", "direction"),
        "PressureLoadCondition": ("modulus",),
        "CouplingPenaltyCondition": ("coupling_partner", "penalty_factor"),
    }
    vector_fields = {"value", "direction"}
    scalar_fields = {"modulus", "penalty_factor"}

    for component_name in queso_model.component_names:
        conditions = queso_model.settings(component_name).get(
            "conditions_settings_list", []
        )
        for condition in conditions:
            condition_id = condition["condition_id"]
            location = f"condition {condition_id}"
            condition_type = condition.get("condition_type")
            if condition_type not in required_fields:
                raise ValueError(
                    f"{location} has unsupported Kratos condition type "
                    f"{condition_type!r}."
                )
            missing = [
                field
                for field in required_fields[condition_type]
                if field not in condition
            ]
            if missing:
                raise ValueError(f"{location} is missing Kratos fields: {missing}.")

            for field in vector_fields & condition.keys():
                value = condition[field]
                if (
                    not isinstance(value, (list, tuple))
                    or len(value) != 3
                    or any(
                        not isinstance(entry, (int, float)) or isinstance(entry, bool)
                        for entry in value
                    )
                ):
                    raise ValueError(f"{location}.{field} must contain three numbers.")
            for field in scalar_fields & condition.keys():
                value = condition[field]
                if not isinstance(value, (int, float)) or isinstance(value, bool):
                    raise ValueError(f"{location}.{field} must be numeric.")

            if condition_type == "CouplingPenaltyCondition":
                if not isinstance(condition.get("slip", False), bool):
                    raise ValueError(f"{location}.slip must be a boolean.")
            elif "coupling_partner" in condition or "slip" in condition:
                raise ValueError(
                    f"{location} uses coupling fields for a non-coupling condition."
                )
