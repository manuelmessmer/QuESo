# Project imports
from QuESoPythonModule.PyQuESo import PyQuESo
from QuESoPythonModule import MeshUtilities as MeshUtil
from QuESoPythonModule.scripts.helper import *
from QuESoPythonModule.scripts.queso_unit_test import QuESoTestCase
# Unittest import
import unittest

class TestBoundaryConditions(QuESoTestCase):
    def check_values(self, pyqueso):
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
                for segment in condition:
                    triangle_mesh_seg1 = segment.GetTriangleMesh()
                    area_segmented1 += MeshUtil.Area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, 332.37754, 5)
                area_segmented2 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg2 = segment.GetTriangleMesh()
                    area_segmented2 += MeshUtil.Area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, 332.37754, 5)
            elif condition_settings.GetString("condition_type") == "PressureLoadCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N2.stl")
                modulus = condition_settings.GetDouble("modulus")
                self.assertAlmostEqual(modulus, 2.0, 10)
                area_segmented1 = 0
                for segment in condition:
                    triangle_mesh_seg1 = segment.GetTriangleMesh()
                    area_segmented1 += MeshUtil.Area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, 577.978141, 5)
                area_segmented2 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg2 = segment.GetTriangleMesh()
                    area_segmented2 += MeshUtil.Area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, 577.978141, 5)
            elif condition_settings.GetString("condition_type") == "LagrangeSupportCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N3.stl")
                value = condition_settings.GetDoubleVector("value")
                self.assertListsAlmostEqual(value, [0.0, 0.3, 0.0], 10)
                area_segmented1 = 0
                for segment in condition:
                    triangle_mesh_seg1 = segment.GetTriangleMesh()
                    area_segmented1 += MeshUtil.Area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, 921.163635, 5)
                area_segmented2 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg2 = segment.GetTriangleMesh()
                    area_segmented2 += MeshUtil.Area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, 921.163635, 5)
            elif condition_settings.GetString("condition_type") == "PenaltySupportCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/D1.stl")
                value = condition_settings.GetDoubleVector("value")
                self.assertListsAlmostEqual(value, [0.0, 0.0, 0.0], 10)
                penalty_factor = condition_settings.GetDouble("penalty_factor")
                self.assertAlmostEqual(penalty_factor, 1e10, 10)
                area_segmented1 = 0
                for segment in condition:
                    triangle_mesh_seg1 = segment.GetTriangleMesh()
                    area_segmented1 += MeshUtil.Area(triangle_mesh_seg1)
                self.assertAlmostEqual(area_segmented1, 1183.54304, 5)
                area_segmented2 = 0
                for segment in condition.GetSegments():
                    triangle_mesh_seg2 = segment.GetTriangleMesh()
                    area_segmented2 += MeshUtil.Area(triangle_mesh_seg2)
                self.assertAlmostEqual(area_segmented2, 1183.54304, 5)
            else:
                raise Exception("TestBoundaryConditions :: Given condition type does not exist.")

    def test_1(self):
        pyqueso = PyQuESo("queso/tests/boundary_conditions/QuESoSettings1.json")
        pyqueso.Run()
        self.check_values(pyqueso)

    def test_2(self):
        pyqueso = PyQuESo("queso/tests/boundary_conditions/QuESoSettings2.json")
        pyqueso.Run()
        self.check_values(pyqueso)

if __name__ == "__main__":
    unittest.main()
