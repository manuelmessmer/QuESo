# Project imports
# Unittest import
import unittest
from pathlib import Path

import pyqueso
from pyqueso.scripts.queso_unit_test import QuESoTestCase

PROJECT_ROOT = Path(__file__).resolve().parents[3]


class TestBoundaryConditions(QuESoTestCase):
    def check_values(self, model, active_areas):
        for condition in model.conditions("main"):
            condition_settings = condition.settings
            if (
                condition_settings.get_string("condition_type")
                == "SurfaceLoadCondition"
            ):
                input_filename = condition_settings.get_string("input_filename")
                self.assertEqual(
                    Path(input_filename),
                    Path("queso/tests/steering_knuckle_kratos/data/N1.stl").resolve(),
                )
                modulus = condition_settings.get_double("modulus")
                self.assertAlmostEqual(modulus, 5.0, 10)
                direction = condition_settings.get_double_vector("direction")
                self.assertListsAlmostEqual(direction, [-1.0, 2.0, 3.0], 10)
                area_segmented1 = 0
                for segment in condition.segments:
                    triangle_mesh_seg1 = segment.triangle_mesh
                    area_segmented1 += pyqueso.mesh.area(triangle_mesh_seg1)
                self.assertAlmostEqual(
                    area_segmented1, active_areas["SurfaceLoadCondition"], 5
                )
                area_segmented2 = 0
                for segment in condition.segments:
                    triangle_mesh_seg2 = segment.triangle_mesh
                    area_segmented2 += pyqueso.mesh.area(triangle_mesh_seg2)
                self.assertAlmostEqual(
                    area_segmented2, active_areas["SurfaceLoadCondition"], 5
                )
            elif (
                condition_settings.get_string("condition_type")
                == "PressureLoadCondition"
            ):
                input_filename = condition_settings.get_string("input_filename")
                self.assertEqual(
                    Path(input_filename),
                    Path("queso/tests/steering_knuckle_kratos/data/N2.stl").resolve(),
                )
                modulus = condition_settings.get_double("modulus")
                self.assertAlmostEqual(modulus, 2.0, 10)
                area_segmented1 = 0
                for segment in condition.segments:
                    triangle_mesh_seg1 = segment.triangle_mesh
                    area_segmented1 += pyqueso.mesh.area(triangle_mesh_seg1)
                self.assertAlmostEqual(
                    area_segmented1, active_areas["PressureLoadCondition"], 5
                )
                area_segmented2 = 0
                for segment in condition.segments:
                    triangle_mesh_seg2 = segment.triangle_mesh
                    area_segmented2 += pyqueso.mesh.area(triangle_mesh_seg2)
                self.assertAlmostEqual(
                    area_segmented2, active_areas["PressureLoadCondition"], 5
                )
            elif (
                condition_settings.get_string("condition_type")
                == "LagrangeSupportCondition"
            ):
                input_filename = condition_settings.get_string("input_filename")
                self.assertEqual(
                    Path(input_filename),
                    Path("queso/tests/steering_knuckle_kratos/data/N3.stl").resolve(),
                )
                value = condition_settings.get_double_vector("value")
                self.assertListsAlmostEqual(value, [0.0, 0.3, 0.0], 10)
                area_segmented1 = 0
                for segment in condition.segments:
                    triangle_mesh_seg1 = segment.triangle_mesh
                    area_segmented1 += pyqueso.mesh.area(triangle_mesh_seg1)
                self.assertAlmostEqual(
                    area_segmented1, active_areas["LagrangeSupportCondition"], 5
                )
                area_segmented2 = 0
                for segment in condition.segments:
                    triangle_mesh_seg2 = segment.triangle_mesh
                    area_segmented2 += pyqueso.mesh.area(triangle_mesh_seg2)
                self.assertAlmostEqual(
                    area_segmented2, active_areas["LagrangeSupportCondition"], 5
                )
            elif (
                condition_settings.get_string("condition_type")
                == "PenaltySupportCondition"
            ):
                input_filename = condition_settings.get_string("input_filename")
                self.assertEqual(
                    Path(input_filename),
                    Path("queso/tests/steering_knuckle_kratos/data/D1.stl").resolve(),
                )
                value = condition_settings.get_double_vector("value")
                self.assertListsAlmostEqual(value, [0.0, 0.0, 0.0], 10)
                penalty_factor = condition_settings.get_double("penalty_factor")
                self.assertAlmostEqual(penalty_factor, 1e10, 10)
                area_segmented1 = 0
                for segment in condition.segments:
                    triangle_mesh_seg1 = segment.triangle_mesh
                    area_segmented1 += pyqueso.mesh.area(triangle_mesh_seg1)
                self.assertAlmostEqual(
                    area_segmented1, active_areas["PenaltySupportCondition"], 5
                )
                area_segmented2 = 0
                for segment in condition.segments:
                    triangle_mesh_seg2 = segment.triangle_mesh
                    area_segmented2 += pyqueso.mesh.area(triangle_mesh_seg2)
                self.assertAlmostEqual(
                    area_segmented2, active_areas["PenaltySupportCondition"], 5
                )
            else:
                raise Exception(
                    "TestBoundaryConditions :: Given condition type does not exist."
                )

    def test_1(self):
        model = pyqueso.Model("queso/tests/boundary_conditions/QuESoSettings1.json")
        model.create()
        self.check_values(
            model,
            {
                "SurfaceLoadCondition": 332.37753991575306,
                "PressureLoadCondition": 577.978140756255,
                "LagrangeSupportCondition": 921.1636351534689,
                "PenaltySupportCondition": 1183.5430441558533,
            },
        )

    def test_2(self):
        model = pyqueso.Model("queso/tests/boundary_conditions/QuESoSettings2.json")
        model.create()
        self.check_values(
            model,
            {
                "SurfaceLoadCondition": 332.377539915753,
                "PressureLoadCondition": 577.978140756255,
                "LagrangeSupportCondition": 921.1636351534688,
                "PenaltySupportCondition": 1183.5430441558533,
            },
        )

    def test_coupling_penalty(self) -> None:
        settings_directory = (
            PROJECT_ROOT / "examples/kratos_analysis/cantilever_coupled"
        )
        model = pyqueso.Model(settings_directory / "QuESoSettings.json")
        model.create()

        left_conditions = {
            condition.settings.get_int("condition_id"): condition
            for condition in model.conditions("left")
        }
        right_condition_ids = {
            condition.settings.get_int("condition_id")
            for condition in model.conditions("right")
        }
        coupling = left_conditions[10]
        settings = coupling.settings

        self.assertNotIn(10, right_condition_ids)
        self.assertEqual(
            settings.get_string("condition_type"), "CouplingPenaltyCondition"
        )
        self.assertEqual(settings.get_string("coupling_partner"), "right")
        self.assertAlmostEqual(settings.get_double("penalty_factor"), 1.0e10)
        self.assertFalse(settings.get_bool("slip"))
        self.assertGreater(coupling.num_segments, 0)
        self.assertTrue(
            all(segment.is_in_active_element for segment in coupling.segments)
        )

        interface_mesh = pyqueso.mesh.TriangleMesh()
        pyqueso.io.read_mesh_from_stl(
            interface_mesh, str(settings_directory / "data/interface.stl")
        )
        segmented_area = sum(
            pyqueso.mesh.area(segment.triangle_mesh) for segment in coupling.segments
        )
        self.assertAlmostEqual(
            segmented_area, pyqueso.mesh.area(interface_mesh), places=8
        )


if __name__ == "__main__":
    unittest.main()
