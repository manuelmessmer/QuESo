"""Unit tests for public settings normalization."""

from pathlib import Path

import pytest
from pyqueso.scripts.settings_parser import parse_settings_by_component


def _grid() -> dict:
    return {
        "grid_type": "b_spline_grid",
        "lower_bound_xyz": [0.0, 0.0, 0.0],
        "upper_bound_xyz": [1.0, 1.0, 1.0],
        "lower_bound_uvw": [0.0, 0.0, 0.0],
        "upper_bound_uvw": [1.0, 1.0, 1.0],
        "polynomial_order": [2, 2, 2],
        "number_of_elements": [1, 1, 1],
    }


def test_defaults_aliases_and_master_coupling(tmp_path: Path) -> None:
    raw_settings = {
        "components": [
            {
                "component_name": "left",
                "input_filename": "left.stl",
                "background_grid_settings": _grid(),
            },
            {
                "component_name": "right",
                "input_filename": "right.stl",
                "background_grid_settings": {"from_other_component": "left"},
            },
        ],
        "conditions": [
            {
                "condition_id": 10,
                "condition_type": "CouplingPenaltyCondition",
                "component_name": "left",
                "coupling_partner": "right",
                "input_filename": "interface.stl",
                "penalty_factor": 1.0e10,
            }
        ],
    }
    settings_path = tmp_path / "case" / "QuESoSettings.json"
    normalized = parse_settings_by_component(raw_settings, settings_path)

    assert tuple(normalized) == ("left", "right")
    assert (
        normalized["left"]["background_grid_settings"]
        == normalized["right"]["background_grid_settings"]
    )
    assert (
        normalized["left"]["background_grid_settings"]
        is not normalized["right"]["background_grid_settings"]
    )
    assert normalized["left"]["conditions_settings_list"][0]["slip"] is False
    assert "conditions_settings_list" not in normalized["right"]
    assert normalized["left"]["echo_level"] == 1
    assert normalized["left"]["write_output_to_file"] is True
    assert (
        Path(normalized["left"]["input_filename"]) == settings_path.parent / "left.stl"
    )


def test_rejects_legacy_schema() -> None:
    with pytest.raises(ValueError, match="unknown fields"):
        parse_settings_by_component(
            {"general_settings": {}}, Path("QuESoSettings.json")
        )


def test_rejects_duplicate_condition_ids() -> None:
    raw_settings = {
        "components": [
            {
                "component_name": "main",
                "input_filename": "main.stl",
                "background_grid_settings": _grid(),
            }
        ],
        "conditions": [
            {
                "condition_id": 1,
                "condition_type": "PressureLoadCondition",
                "component_name": "main",
                "input_filename": "a.stl",
                "modulus": 1.0,
            },
            {
                "condition_id": 1,
                "condition_type": "PressureLoadCondition",
                "component_name": "main",
                "input_filename": "b.stl",
                "modulus": 1.0,
            },
        ],
    }
    with pytest.raises(ValueError, match="Duplicate condition_id"):
        parse_settings_by_component(raw_settings, Path("QuESoSettings.json"))


def test_passes_cpp_owned_fields_through() -> None:
    raw_settings = {
        "components": [
            {
                "component_name": "main",
                "input_filename": "main.stl",
                "future_cpp_setting": 42,
                "background_grid_settings": {
                    **_grid(),
                    "future_grid_setting": "kept",
                },
            }
        ],
        "conditions": [
            {
                "condition_id": 1,
                "condition_type": "CustomCondition",
                "component_name": "main",
                "input_filename": "custom.stl",
                "future_condition_setting": 3.5,
            }
        ],
    }
    normalized = parse_settings_by_component(raw_settings, Path("QuESoSettings.json"))[
        "main"
    ]
    assert normalized["future_cpp_setting"] == 42
    assert normalized["background_grid_settings"]["future_grid_setting"] == "kept"
    assert normalized["conditions_settings_list"][0]["future_condition_setting"] == 3.5
