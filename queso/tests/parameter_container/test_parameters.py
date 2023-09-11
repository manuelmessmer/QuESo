# Project imports
import QuESo_PythonApplication
from queso.python_scripts.helper import *

import unittest

class TestParametersContainer(unittest.TestCase):
    def check_values(self, parameters):

        echo_level = parameters.EchoLevel()
        self.assertEqual(echo_level, 1)

        lower_bound_xyz = parameters.LowerBoundXYZ()
        self.assertAlmostEqual(lower_bound_xyz, [-130.0, -110.0, -110.0], 10)
        lower_bound_uvw = parameters.LowerBoundUVW()
        self.assertAlmostEqual(lower_bound_uvw, [-130.0, -110.0, -110.0], 10)

        upper_bound_xyz = parameters.UpperBoundXYZ()
        self.assertAlmostEqual(upper_bound_xyz, [20.0, 190.0, 190.0], 10)
        upper_bound_uvw = parameters.UpperBoundUVW()
        self.assertAlmostEqual(upper_bound_uvw, [20.0, 190.0, 190.0], 10)

        polynomial_order = parameters.Order()
        self.assertAlmostEqual(polynomial_order, [2, 2, 2], 10)

        number_of_elements = parameters.NumberOfElements()
        self.assertAlmostEqual(number_of_elements, [5, 5, 5], 10)

        integration_method = parameters.IntegrationMethod()
        self.assertEqual(integration_method, QuESo_PythonApplication.IntegrationMethod.Gauss)

        global_parameters = parameters.GetGlobalSettings()

        input_filename = global_parameters.GetString("input_filename")
        self.assertEqual(input_filename, "dummy.stl")

        moment_fitting_residual = global_parameters.GetDouble("moment_fitting_residual")
        self.assertAlmostEqual(moment_fitting_residual, 1e-8, 10)

        min_element_volume_ratio = global_parameters.GetDouble("min_element_volume_ratio")
        self.assertAlmostEqual(min_element_volume_ratio, 1e-3, 10)

        b_slpine_mesh = global_parameters.GetBool("b_spline_mesh")
        self.assertEqual(b_slpine_mesh, True)

        for condition_settings in parameters.GetConditionsSettingsVector():
            if condition_settings.GetString("type") == "SurfaceLoadCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "SurfaceLoadCondition.stl")
                modulus = condition_settings.GetDouble("modulus")
                self.assertAlmostEqual(modulus, 5.0, 10)
                direction = condition_settings.GetDoubleVector("direction")
                self.assertAlmostEqual(direction, [-1.0, 2.0, 3.0], 10)
            elif condition_settings.GetString("type") == "PressureLoadCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "PressureLoadCondition.stl")
                modulus = condition_settings.GetDouble("modulus")
                self.assertAlmostEqual(modulus, 2.0, 10)
            elif condition_settings.GetString("type") == "LagrangeSupportCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "LagrangeSupportCondition.stl")
                value = condition_settings.GetDoubleVector("value")
                self.assertAlmostEqual(value, [0.0, 0.3, 0.0], 10)
            elif condition_settings.GetString("type") == "PenaltySupportCondition":
                input_filename = condition_settings.GetString("input_filename")
                self.assertEqual(input_filename, "PenaltySupportCondition.stl")
                value = condition_settings.GetDoubleVector("value")
                self.assertAlmostEqual(value, [0.0, 0.0, 0.0], 10)
                penalty_factor = condition_settings.GetDouble("penalty_factor")
                self.assertAlmostEqual(penalty_factor, 1e10, 10)

    def test_1(self):
        parameters = ReadParameters("queso/tests/parameter_container/QuESoParameters1.json")
        self.check_values(parameters)

    def test_2(self):
        parameters = ReadParameters("queso/tests/parameter_container/QuESoParameters2.json")
        self.check_values(parameters)

if __name__ == "__main__":
    unittest.main()
