# Project imports
from QuESo_PythonApplication.PyQuESo import PyQuESo
from QuESo_PythonApplication import ConditionSettings as ConSet
from queso.python_scripts.helper import *
from queso.python_scripts.QuESoUnittest import QuESoTestCase
# Unittest import
import unittest

class TestBoundaryConditions(QuESoTestCase):
    def check_triangle_mesh(self, triangle_mesh, ref_area):
        num_triangles = triangle_mesh.NumOfTriangles()
        area = 0.0
        for id in range(num_triangles):
            area += triangle_mesh.Area(id)
        self.assertAlmostEqual(area, ref_area, 5)

    def check_values(self, pyqueso):
        for condition in pyqueso.GetConditions():
            condition_settings = condition.GetSettings()
            if condition_settings.GetString(ConSet.condition_type) == "SurfaceLoadCondition":
                triangle_mesh = condition.GetTriangleMesh()
                input_filename = condition_settings.GetString(ConSet.input_filename)
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N1.stl")
                modulus = condition_settings.GetDouble(ConSet.modulus)
                self.assertAlmostEqual(modulus, 5.0, 10)
                direction = condition_settings.GetDoubleVector(ConSet.direction)
                self.assertListsAlmostEqual(direction, [-1.0, 2.0, 3.0], 10)

                triangle_mesh = condition.GetTriangleMesh()
                self.check_triangle_mesh(triangle_mesh,  332.37754)
            elif condition_settings.GetString(ConSet.condition_type) == "PressureLoadCondition":
                input_filename = condition_settings.GetString(ConSet.input_filename)
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N2.stl")
                modulus = condition_settings.GetDouble(ConSet.modulus)
                self.assertAlmostEqual(modulus, 2.0, 10)

                triangle_mesh = condition.GetTriangleMesh()
                self.check_triangle_mesh(triangle_mesh, 577.978141)
            elif condition_settings.GetString(ConSet.condition_type) == "LagrangeSupportCondition":
                input_filename = condition_settings.GetString(ConSet.input_filename)
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/N3.stl")
                value = condition_settings.GetDoubleVector(ConSet.value)
                self.assertListsAlmostEqual(value, [0.0, 0.3, 0.0], 10)

                triangle_mesh = condition.GetTriangleMesh()
                self.check_triangle_mesh(triangle_mesh, 921.163635)
            elif condition_settings.GetString(ConSet.condition_type) == "PenaltySupportCondition":
                input_filename = condition_settings.GetString(ConSet.input_filename)
                self.assertEqual(input_filename, "queso/tests/steering_knuckle_kratos/data/D1.stl")
                value = condition_settings.GetDoubleVector(ConSet.value)
                self.assertListsAlmostEqual(value, [0.0, 0.0, 0.0], 10)
                penalty_factor = condition_settings.GetDouble(ConSet.penalty_factor)
                self.assertAlmostEqual(penalty_factor, 1e10, 10)

                triangle_mesh = condition.GetTriangleMesh()
                self.check_triangle_mesh(triangle_mesh, 1183.54304)
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
