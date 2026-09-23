"""Tests for multi-component settings, model creation, and Kratos assembly."""

import tempfile
import unittest
from pathlib import Path

import pyqueso
from pyqueso.scripts.settings_parser import parse_settings_by_component

PROJECT_ROOT = Path(__file__).resolve().parents[3]


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


class TestSettingsParser(unittest.TestCase):
    """Validate the clean-break public settings contract."""

    def test_defaults_aliases_and_master_coupling(self) -> None:
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
        settings_path = Path(tempfile.gettempdir()) / "case" / "QuESoSettings.json"
        normalized = parse_settings_by_component(raw_settings, settings_path)

        self.assertEqual(tuple(normalized), ("left", "right"))
        self.assertEqual(
            normalized["left"]["background_grid_settings"],
            normalized["right"]["background_grid_settings"],
        )
        self.assertIsNot(
            normalized["left"]["background_grid_settings"],
            normalized["right"]["background_grid_settings"],
        )
        self.assertEqual(
            normalized["left"]["conditions_settings_list"][0]["slip"], False
        )
        self.assertNotIn("conditions_settings_list", normalized["right"])
        self.assertEqual(normalized["left"]["echo_level"], 1)
        self.assertTrue(normalized["left"]["write_output_to_file"])
        self.assertEqual(
            Path(normalized["left"]["input_filename"]),
            settings_path.parent / "left.stl",
        )

    def test_rejects_legacy_schema(self) -> None:
        with self.assertRaisesRegex(ValueError, "unknown fields"):
            parse_settings_by_component(
                {"general_settings": {}}, Path("QuESoSettings.json")
            )

    def test_rejects_duplicate_condition_ids(self) -> None:
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
        with self.assertRaisesRegex(ValueError, "Duplicate condition_id"):
            parse_settings_by_component(raw_settings, Path("QuESoSettings.json"))

    def test_passes_cpp_owned_fields_through(self) -> None:
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
        normalized = parse_settings_by_component(
            raw_settings, Path("QuESoSettings.json")
        )["main"]
        self.assertEqual(normalized["future_cpp_setting"], 42)
        self.assertEqual(
            normalized["background_grid_settings"]["future_grid_setting"], "kept"
        )
        self.assertEqual(
            normalized["conditions_settings_list"][0]["future_condition_setting"],
            3.5,
        )


class TestMultiComponentModel(unittest.TestCase):
    """Exercise named model access and lifecycle behavior."""

    def test_create_and_named_access(self) -> None:
        model = pyqueso.Model(
            PROJECT_ROOT
            / "examples/kratos_analysis/cantilever_coupled/QuESoSettings.json"
        )
        self.assertEqual(model.component_names, ("left", "right"))
        with self.assertRaises(RuntimeError):
            model.elements("left")

        model.create()
        self.assertGreater(len(model.elements("left")), 0)
        self.assertGreater(len(model.elements("right")), 0)
        self.assertEqual(
            [
                condition.settings.get_int("condition_id")
                for condition in model.conditions("left")
            ],
            [1, 10],
        )
        self.assertTrue(
            all(
                isinstance(segment.is_in_active_element, bool)
                for condition in model.conditions("left")
                for segment in condition.segments
            )
        )
        with self.assertRaises(RuntimeError):
            model.create()
        with self.assertRaises(KeyError):
            model.settings("missing")


try:
    from pyqueso.kratos_interface import Analysis
    from pyqueso.kratos_interface.analysis_settings_parser import (
        parse_analysis_settings,
    )

    KRATOS_AVAILABLE = True
except ImportError:
    KRATOS_AVAILABLE = False


@unittest.skipUnless(KRATOS_AVAILABLE, "KratosMultiphysics is unavailable")
class TestKratosAnalysis(unittest.TestCase):
    """Exercise automatic modeler injection and coupled entity assembly."""

    def test_coupled_cantilever(self) -> None:
        directory = PROJECT_ROOT / "examples/kratos_analysis/cantilever_coupled"
        analysis = Analysis(
            queso_settings_path=directory / "QuESoSettings.json",
            analysis_settings_path=directory / "AnalysisSettings.json",
            kratos_parameters_path=directory / "KratosParameters.json",
        )
        with self.assertRaises(RuntimeError):
            _ = analysis.queso_model
        with self.assertRaises(RuntimeError):
            _ = analysis.kratos_model

        analysis.run()
        root = analysis.kratos_model.GetModelPart("Structure")
        self.assertTrue(root.HasSubModelPart("left"))
        self.assertTrue(root.HasSubModelPart("right"))
        self.assertGreater(root.GetSubModelPart("left").NumberOfElements(), 0)
        self.assertGreater(root.GetSubModelPart("right").NumberOfElements(), 0)
        self.assertGreater(root.GetSubModelPart("left").NumberOfConditions(), 0)
        self.assertGreater(root.GetSubModelPart("right").NumberOfConditions(), 0)
        self.assertTrue(
            any("Coupling" in condition.Info() for condition in root.Conditions)
        )
        self.assertTrue(
            (
                PROJECT_ROOT / "kratos_output/left/EmbeddedModelPart_left_0_0.vtk"
            ).is_file()
        )
        with self.assertRaises(RuntimeError):
            analysis.run()

    def test_rejects_condition_semantics_it_does_not_support(self) -> None:
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
        with self.assertRaisesRegex(ValueError, "unsupported Kratos condition"):
            parse_analysis_settings(analysis_settings, FakeModel())


if __name__ == "__main__":
    unittest.main()
