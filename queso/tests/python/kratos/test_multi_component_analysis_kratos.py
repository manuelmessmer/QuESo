"""Kratos integration tests for multi-component analysis assembly."""

from collections.abc import Callable
from pathlib import Path

import pytest

pytest.importorskip("KratosMultiphysics")

from pyqueso.kratos_interface import Analysis
from pyqueso.kratos_interface.analysis_settings_parser import parse_analysis_settings

pytestmark = pytest.mark.kratos


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


def test_coupled_cantilever(
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("coupled_cantilever")
    analysis = Analysis(
        queso_settings_path=directory / "QuESoSettings.json",
        analysis_settings_path=directory / "AnalysisSettings.json",
        kratos_parameters_path=directory / "KratosParameters.json",
    )
    with pytest.raises(RuntimeError):
        _ = analysis.queso_model
    with pytest.raises(RuntimeError):
        _ = analysis.kratos_model

    analysis.run()
    root = analysis.kratos_model.GetModelPart("Structure")
    assert root.HasSubModelPart("left")
    assert root.HasSubModelPart("right")
    assert root.GetSubModelPart("left").NumberOfElements() > 0
    assert root.GetSubModelPart("right").NumberOfElements() > 0
    assert root.GetSubModelPart("left").NumberOfConditions() > 0
    assert root.GetSubModelPart("right").NumberOfConditions() > 0
    assert any("Coupling" in condition.Info() for condition in root.Conditions)
    assert (
        directory.parent / "kratos_output/left/EmbeddedModelPart_left_0_0.vtk"
    ).is_file()
    with pytest.raises(RuntimeError):
        analysis.run()


def test_rejects_unsupported_condition_semantics() -> None:
    class FakeModel:
        component_names = ("main",)

        @staticmethod
        def settings(component_name: str) -> dict:
            return {
                "background_grid_settings": _grid(),
                "conditions_settings_list": [
                    {
                        "condition_id": 1,
                        "condition_type": "CustomCondition",
                    }
                ],
            }

    analysis_settings = {
        "model_parts": [
            {
                "name": "Structure",
                "geometries": [
                    {
                        "name": "Volume",
                        "queso_components": [
                            {
                                "component_name": "main",
                                "element_settings": {"property_id": 1},
                            }
                        ],
                    }
                ],
            }
        ]
    }
    with pytest.raises(ValueError, match="unsupported Kratos condition"):
        parse_analysis_settings(analysis_settings, FakeModel())
