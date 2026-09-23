# Project imports
import pyqueso
from pyqueso.scripts.helper import *
from pyqueso.scripts.queso_unit_test import QuESoTestCase
# Unittest import
import unittest

class TestBoundaryConditions(QuESoTestCase):
    def check_values(self, model, active_areas):
        for condition in model.conditions:
            condition_settings = condition.settings
            if condition_settings.get_string("condition_type") == "SurfaceLoadCondition":
                input_filename = condition_settings.get_string("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N1.stl")
                modulus = condition_settings.get_double("modulus")
                self.assertAlmostEqual(modulus, 5.0, 10)
                direction = condition_settings.get_double_vector("direction")
                self.assertListsAlmostEqual(direction, [-1.0, 2.0, 3.0], 10)
                area_segmented1 = 0
                for segment in condition.segments:
                    triangle_mesh_seg1 = segment.triangle_mesh
                    area_segmented1 += pyqueso.mesh.area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, active_areas["SurfaceLoadCondition"], 5)
                area_segmented2 = 0
                for segment in condition.segments:
                    triangle_mesh_seg2 = segment.triangle_mesh
                    area_segmented2 += pyqueso.mesh.area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, active_areas["SurfaceLoadCondition"], 5)
            elif condition_settings.get_string("condition_type") == "PressureLoadCondition":
                input_filename = condition_settings.get_string("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N2.stl")
                modulus = condition_settings.get_double("modulus")
                self.assertAlmostEqual(modulus, 2.0, 10)
                area_segmented1 = 0
                for segment in condition.segments:
                    triangle_mesh_seg1 = segment.triangle_mesh
                    area_segmented1 += pyqueso.mesh.area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, active_areas["PressureLoadCondition"], 5)
                area_segmented2 = 0
                for segment in condition.segments:
                    triangle_mesh_seg2 = segment.triangle_mesh
                    area_segmented2 += pyqueso.mesh.area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, active_areas["PressureLoadCondition"], 5)
            elif condition_settings.get_string("condition_type") == "LagrangeSupportCondition":
                input_filename = condition_settings.get_string("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N3.stl")
                value = condition_settings.get_double_vector("value")
                self.assertListsAlmostEqual(value, [0.0, 0.3, 0.0], 10)
                area_segmented1 = 0
                for segment in condition.segments:
                    triangle_mesh_seg1 = segment.triangle_mesh
                    area_segmented1 += pyqueso.mesh.area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, active_areas["LagrangeSupportCondition"], 5)
                area_segmented2 = 0
                for segment in condition.segments:
                    triangle_mesh_seg2 = segment.triangle_mesh
                    area_segmented2 += pyqueso.mesh.area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, active_areas["LagrangeSupportCondition"], 5)
            elif condition_settings.get_string("condition_type") == "PenaltySupportCondition":
                input_filename = condition_settings.get_string("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/D1.stl")
                value = condition_settings.get_double_vector("value")
                self.assertListsAlmostEqual(value, [0.0, 0.0, 0.0], 10)
                penalty_factor = condition_settings.get_double("penalty_factor")
                self.assertAlmostEqual(penalty_factor, 1e10, 10)
                area_segmented1 = 0
                for segment in condition.segments:
                    triangle_mesh_seg1 = segment.triangle_mesh
                    area_segmented1 += pyqueso.mesh.area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, active_areas["PenaltySupportCondition"], 5)
                area_segmented2 = 0
                for segment in condition.segments:
                    triangle_mesh_seg2 = segment.triangle_mesh
                    area_segmented2 += pyqueso.mesh.area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, active_areas["PenaltySupportCondition"], 5)
            else:
                raise Exception("TestBoundaryConditions :: Given condition type does not exist.")

    def test_1(self):
        model = pyqueso.Model(json_filename="queso/tests/boundary_conditions/QuESoSettings1.json")
        model.create()
        self.check_values(model, {
            "SurfaceLoadCondition": 332.37753991575306,
            "PressureLoadCondition": 577.978140756255,
            "LagrangeSupportCondition": 921.1636351534689,
            "PenaltySupportCondition": 1183.5430441558533,
        })

    def test_2(self):
        model = pyqueso.Model(json_filename="queso/tests/boundary_conditions/QuESoSettings2.json")
        model.create()
        self.check_values(model, {
            "SurfaceLoadCondition": 332.377539915753,
            "PressureLoadCondition": 577.978140756255,
            "LagrangeSupportCondition": 921.1636351534688,
            "PenaltySupportCondition": 1183.5430441558533,
        })

if __name__ == "__main__":
    unittest.main()
