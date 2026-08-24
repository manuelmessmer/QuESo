# Project imports
from QuESoPythonModule.PyQuESo import PyQuESo
from QuESoPythonModule import MeshUtil
from QuESoPythonModule.scripts.helper import *
from QuESoPythonModule.scripts.queso_unit_test import QuESoTestCase
# Unittest import
import unittest

class TestBoundaryConditions(QuESoTestCase):
    def check_values(self, pyqueso, active_areas):
        for condition in pyqueso.GetConditions():
            condition_settings = condition.GetSettings()
            if condition_settings.GetString("condition_type") == "SurfaceLoadCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N1.stl")
                modulus = condition_settings.GetDouble("modulus")
                self.assertAlmostEqual(modulus, 5.0, 10)
                direction = condition_settings.GetDoubleVector("direction")
                self.assertListsAlmostEqual(direction, [-1.0, 2.0, 3.0], 10)
                area_segmented1 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg1 = segment.GetTriangleMesh()
                    area_segmented1 += MeshUtil.Area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, active_areas["SurfaceLoadCondition"], 5)
                area_segmented2 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg2 = segment.GetTriangleMesh()
                    area_segmented2 += MeshUtil.Area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, active_areas["SurfaceLoadCondition"], 5)
            elif condition_settings.GetString("condition_type") == "PressureLoadCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N2.stl")
                modulus = condition_settings.GetDouble("modulus")
                self.assertAlmostEqual(modulus, 2.0, 10)
                area_segmented1 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg1 = segment.GetTriangleMesh()
                    area_segmented1 += MeshUtil.Area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, active_areas["PressureLoadCondition"], 5)
                area_segmented2 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg2 = segment.GetTriangleMesh()
                    area_segmented2 += MeshUtil.Area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, active_areas["PressureLoadCondition"], 5)
            elif condition_settings.GetString("condition_type") == "LagrangeSupportCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N3.stl")
                value = condition_settings.GetDoubleVector("value")
                self.assertListsAlmostEqual(value, [0.0, 0.3, 0.0], 10)
                area_segmented1 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg1 = segment.GetTriangleMesh()
                    area_segmented1 += MeshUtil.Area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, active_areas["LagrangeSupportCondition"], 5)
                area_segmented2 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg2 = segment.GetTriangleMesh()
                    area_segmented2 += MeshUtil.Area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, active_areas["LagrangeSupportCondition"], 5)
            elif condition_settings.GetString("condition_type") == "PenaltySupportCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/D1.stl")
                value = condition_settings.GetDoubleVector("value")
                self.assertListsAlmostEqual(value, [0.0, 0.0, 0.0], 10)
                penalty_factor = condition_settings.GetDouble("penalty_factor")
                self.assertAlmostEqual(penalty_factor, 1e10, 10)
                area_segmented1 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg1 = segment.GetTriangleMesh()
                    area_segmented1 += MeshUtil.Area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, active_areas["PenaltySupportCondition"], 5)
                area_segmented2 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg2 = segment.GetTriangleMesh()
                    area_segmented2 += MeshUtil.Area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, active_areas["PenaltySupportCondition"], 5)
            else:
                raise Exception("TestBoundaryConditions :: Given condition type does not exist.")

    def test_1(self):
        pyqueso = PyQuESo("queso/tests/boundary_conditions/QuESoSettings1.json")
        pyqueso.Run()
        self.check_values(pyqueso, {
            "SurfaceLoadCondition": 332.37753991575306,
            "PressureLoadCondition": 577.978140756255,
            "LagrangeSupportCondition": 921.1636351534689,
            "PenaltySupportCondition": 1183.5430441558533,
        })

    def test_2(self):
        pyqueso = PyQuESo("queso/tests/boundary_conditions/QuESoSettings2.json")
        pyqueso.Run()
        self.check_values(pyqueso, {
            "SurfaceLoadCondition": 332.377539915753,
            "PressureLoadCondition": 577.978140756255,
            "LagrangeSupportCondition": 921.1636351534688,
            "PenaltySupportCondition": 1183.5430441558533,
        })

if __name__ == "__main__":
    unittest.main()
