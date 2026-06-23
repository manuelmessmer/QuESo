"""Helpers for parsing and validating the AnalysisSettings JSON structure."""

from typing import Any


REQUIRED_QUESO_COMPONENT_FIELDS: dict[str, type] = {
    "component_name": str,
    "element_settings": dict,
}

VALID_FIXED_DOF_NAMES = {
    "DISPLACEMENT_X",
    "DISPLACEMENT_Y",
    "DISPLACEMENT_Z",
}


def parse_geometry_analysis_settings(
    analysis_settings: dict[str, Any],
    queso_model: Any,
) -> dict[str, Any]:
    """Normalize AnalysisSettings into a fully validated analysis settings map."""
    root_model_part_name = _normalize_root_model_part_name(analysis_settings)
    model_parts = analysis_settings.get("model_parts")
    if not isinstance(model_parts, list) or not model_parts:
        raise ValueError("AnalysisSettings must contain a non-empty 'model_parts' list.")

    component_names = set(queso_model.component_names)
    geometries: dict[str, dict[str, Any]] = {}
    lagrange_dofs_required = False

    for model_part_entry in model_parts:
        if not isinstance(model_part_entry, dict):
            raise ValueError("Each model part entry must be a JSON object.")

        model_part_name = str(model_part_entry.get("name", "")).strip()
        if not model_part_name:
            raise ValueError("Each model part must define a non-empty 'name'.")

        kratos_model_part_name = f"{root_model_part_name}.{model_part_name}"

        model_part_geometries = model_part_entry.get("geometries")
        if not isinstance(model_part_geometries, list) or not model_part_geometries:
            raise ValueError(f"Model part '{model_part_name}' must define a non-empty 'geometries' list.")

        for geometry_entry in model_part_geometries:
            if not isinstance(geometry_entry, dict):
                raise ValueError(f"Each geometry in model part '{model_part_name}' must be a JSON object.")

            geometry_name = str(geometry_entry.get("name", "")).strip()
            if not geometry_name:
                raise ValueError(f"Each geometry in model part '{model_part_name}' must define a non-empty 'name'.")
            if geometry_name in geometries:
                raise ValueError(f"Duplicate geometry name '{geometry_name}'.")

            queso_components = geometry_entry.get("queso_components")
            if not isinstance(queso_components, list) or not queso_components:
                raise ValueError(f"Geometry '{geometry_name}' must define a non-empty 'queso_components' list.")

            normalized_components = _validate_and_normalize_queso_components(
                geometry_name,
                queso_components,
                component_names,
                queso_model,
            )
            background_grid_settings = _resolve_geometry_background_grid_settings(
                geometry_name,
                normalized_components,
                queso_model,
            )
            geometry_requires_lagrange_dofs = _geometry_requires_lagrange_dofs(
                normalized_components,
                queso_model,
            )
            lagrange_dofs_required = lagrange_dofs_required or geometry_requires_lagrange_dofs

            geometries[geometry_name] = {
                "model_part_name": kratos_model_part_name,
                "background_grid_settings": background_grid_settings,
                "queso_components": normalized_components,
                "fixed_model_parts": _validate_fixed_model_parts(
                    geometry_name,
                    geometry_entry.get("fixed_model_parts", []),
                    component_names,
                ),
            }

    return {
        "root_model_part_name": root_model_part_name,
        "lagrange_dofs_required": lagrange_dofs_required,
        "geometries": geometries,
    }


def _normalize_root_model_part_name(analysis_settings: dict[str, Any]) -> str:
    """Return the normalized root model part name."""
    root_model_part_name = str(analysis_settings.get("root_model_part_name", "Structure")).strip()
    return root_model_part_name or "Structure"


def _validate_and_normalize_queso_components(
    geometry_name: str,
    queso_components: list[dict[str, Any]],
    component_names: set[str],
    queso_model: Any,
) -> list[dict[str, Any]]:
    """Validate QuESo component entries for a geometry."""
    normalized_components: list[dict[str, Any]] = []
    for queso_component in queso_components:
        if not isinstance(queso_component, dict):
            raise ValueError(f"Each QuESo component in geometry '{geometry_name}' must be a JSON object.")

        for field, expected_type in REQUIRED_QUESO_COMPONENT_FIELDS.items():
            if field not in queso_component:
                raise ValueError(
                    f"QuESo component in geometry '{geometry_name}' is missing required field '{field}'."
                )
            if not isinstance(queso_component[field], expected_type):
                raise ValueError(
                    f"QuESo component in geometry '{geometry_name}' field '{field}' must be "
                    f"{expected_type.__name__}."
                )

        component_name = str(queso_component["component_name"]).strip()
        if component_name not in component_names:
            raise ValueError(
                f"Geometry '{geometry_name}' references unknown component '{component_name}'."
            )

        element_settings = _validate_element_settings(geometry_name, component_name, queso_component["element_settings"])
        condition_settings = _validate_condition_settings(
            geometry_name,
            component_name,
            queso_component.get("condition_settings", {}),
            queso_model,
        )

        normalized_components.append(
            {
                "component_name": component_name,
                "element_settings": element_settings,
                "condition_settings": condition_settings,
            }
        )

    return normalized_components


def _validate_element_settings(
    geometry_name: str,
    component_name: str,
    element_settings: dict[str, Any],
) -> dict[str, Any]:
    """Validate element settings for one geometry/component pair."""
    if not isinstance(element_settings, dict):
        raise ValueError(
            f"Geometry '{geometry_name}' component '{component_name}' field 'element_settings' must be an object."
        )

    property_id = element_settings.get("property_id")
    if not isinstance(property_id, int):
        raise ValueError(
            f"Geometry '{geometry_name}' component '{component_name}' must define integer element_settings.property_id."
        )

    normalized_settings = {"property_id": property_id}

    element_type = element_settings.get("element_type")
    if element_type is not None:
        if not isinstance(element_type, str):
            raise ValueError(
                f"Geometry '{geometry_name}' component '{component_name}' field 'element_settings.element_type' must be str."
            )
        normalized_settings["element_type"] = element_type

    polynomial_order = element_settings.get("polynomial_order")
    if polynomial_order is not None:
        normalized_settings["polynomial_order"] = _validate_polynomial_order(
            geometry_name,
            component_name,
            polynomial_order,
        )

    return normalized_settings


def _validate_polynomial_order(
    geometry_name: str,
    component_name: str,
    polynomial_order: Any,
) -> list[int]:
    """Validate an optional per-geometry polynomial-order override."""
    if not isinstance(polynomial_order, list) or len(polynomial_order) != 3:
        raise ValueError(
            f"Geometry '{geometry_name}' component '{component_name}' field "
            "'element_settings.polynomial_order' must be a list of three integers."
        )

    normalized_order: list[int] = []
    for value in polynomial_order:
        if not isinstance(value, int):
            raise ValueError(
                f"Geometry '{geometry_name}' component '{component_name}' field "
                "'element_settings.polynomial_order' must contain only integers."
            )
        if value < 1:
            raise ValueError(
                f"Geometry '{geometry_name}' component '{component_name}' field "
                "'element_settings.polynomial_order' must contain only integers >= 1."
            )
        normalized_order.append(value)

    return normalized_order


def _validate_condition_settings(
    geometry_name: str,
    component_name: str,
    condition_settings: dict[str, Any],
    queso_model: Any,
) -> dict[str, Any]:
    """Validate condition settings and referenced condition ids."""
    if not isinstance(condition_settings, dict):
        raise ValueError(
            f"Geometry '{geometry_name}' component '{component_name}' field 'condition_settings' must be an object."
        )

    active_condition_ids = condition_settings.get("active_condition_ids", [])
    if not isinstance(active_condition_ids, list):
        raise ValueError(
            f"Geometry '{geometry_name}' component '{component_name}' field 'condition_settings.active_condition_ids' must be a list."
        )

    available_ids = _get_conditions_by_id(component_name, queso_model)
    normalized_condition_ids: list[int] = []
    for condition_id in active_condition_ids:
        normalized_condition_id = int(condition_id)
        if normalized_condition_id not in available_ids:
            raise ValueError(
                f"Geometry '{geometry_name}' component '{component_name}' "
                f"references unknown condition_id {condition_id}."
            )
        normalized_condition_ids.append(normalized_condition_id)

    return {"active_condition_ids": normalized_condition_ids}


def _validate_fixed_model_parts(
    geometry_name: str,
    fixed_model_parts: Any,
    component_names: set[str],
) -> list[dict[str, Any]]:
    """Validate solver-side nodal fixing instructions for component submodelparts."""
    if fixed_model_parts in (None, []):
        return []
    if not isinstance(fixed_model_parts, list):
        raise ValueError(
            f"Geometry '{geometry_name}' field 'fixed_model_parts' must be a list."
        )

    normalized_entries: list[dict[str, Any]] = []
    for entry in fixed_model_parts:
        if not isinstance(entry, dict):
            raise ValueError(
                f"Each fixed model part in geometry '{geometry_name}' must be an object."
            )

        model_part_name = str(entry.get("model_part_name", "")).strip()
        if not model_part_name:
            raise ValueError(
                f"Each fixed model part in geometry '{geometry_name}' must define non-empty 'model_part_name'."
            )
        if model_part_name not in component_names:
            raise ValueError(
                f"Geometry '{geometry_name}' fixed model part references unknown component '{model_part_name}'."
            )

        fix = entry.get("fix")
        if isinstance(fix, str):
            fix_names = [fix]
        elif isinstance(fix, list) and all(isinstance(item, str) for item in fix):
            fix_names = list(fix)
        else:
            raise ValueError(
                f"Geometry '{geometry_name}' fixed model part '{model_part_name}' field 'fix' must be a string or list of strings."
            )

        normalized_fix_names = []
        for fix_name in fix_names:
            dof_name = str(fix_name).strip().upper()
            if dof_name not in VALID_FIXED_DOF_NAMES:
                raise ValueError(
                    f"Geometry '{geometry_name}' fixed model part '{model_part_name}' uses unsupported dof '{fix_name}'."
                )
            normalized_fix_names.append(dof_name)

        normalized_entries.append(
            {
                "model_part_name": model_part_name,
                "fix": normalized_fix_names,
            }
        )

    return normalized_entries


def _resolve_geometry_background_grid_settings(
    geometry_name: str,
    normalized_components: list[dict[str, Any]],
    queso_model: Any,
) -> dict[str, Any]:
    """Resolve and validate the shared background grid for a geometry."""
    component_names = [str(component["component_name"]) for component in normalized_components]
    if not component_names:
        raise ValueError(f"Geometry '{geometry_name}' does not reference any components.")

    required_grid_keys = [
        "lower_bound_xyz",
        "upper_bound_xyz",
        "lower_bound_uvw",
        "upper_bound_uvw",
        "number_of_elements",
        "polynomial_order",
    ]

    reference = dict(
        queso_model.settings(component_names[0]).get("background_grid_settings", {})
    )
    for key in required_grid_keys:
        if key not in reference:
            raise ValueError(
                f"Component '{component_names[0]}' is missing background_grid_settings.{key}."
            )

    overridden_polynomial_order = normalized_components[0]["element_settings"].get(
        "polynomial_order"
    )
    if overridden_polynomial_order is not None:
        reference["polynomial_order"] = overridden_polynomial_order

    mismatched: list[str] = []
    reference_component_name = component_names[0]
    for component_settings in normalized_components[1:]:
        component_name = str(component_settings["component_name"])
        current_grid_settings = queso_model.settings(component_name).get(
            "background_grid_settings", {}
        )
        current_polynomial_order = component_settings["element_settings"].get(
            "polynomial_order"
        )

        if current_polynomial_order is None:
            current_effective_polynomial_order = current_grid_settings.get("polynomial_order")
        else:
            current_effective_polynomial_order = current_polynomial_order

        current_reference = dict(current_grid_settings)
        current_reference["polynomial_order"] = current_effective_polynomial_order

        if current_reference != reference:
            mismatched.append(component_name)

    if mismatched:
        raise ValueError(
            f"Geometry '{geometry_name}' references components with different background grids or "
            f"element_settings.polynomial_order overrides: reference='{reference_component_name}', "
            f"mismatched={mismatched}."
        )

    return reference


def _geometry_requires_lagrange_dofs(
    normalized_components: list[dict[str, Any]],
    queso_model: Any,
) -> bool:
    """Return whether any active condition in the geometry requires Lagrange dofs."""
    for queso_component in normalized_components:
        component_name = str(queso_component["component_name"])
        conditions_by_id = _get_conditions_by_id(component_name, queso_model)
        active_condition_ids = queso_component["condition_settings"].get("active_condition_ids", [])
        for condition_id in active_condition_ids:
            condition_settings = conditions_by_id[int(condition_id)].GetSettings()
            if condition_settings.GetString("condition_type") == "LagrangeSupportCondition":
                return True
    return False


def _get_conditions_by_id(component_name: str, queso_model: Any) -> dict[int, Any]:
    """Return generated conditions for one component keyed by condition id."""
    return {
        condition.GetSettings().GetInt("condition_id"): condition
        for condition in queso_model.conditions(component_name)
    }

