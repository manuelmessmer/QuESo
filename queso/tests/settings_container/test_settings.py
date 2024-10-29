# Project imports
import QuESo_PythonApplication as QuESo
from queso.python_scripts.json_import import JsonImport
from queso.python_scripts.QuESoUnittest import QuESoTestCase
# External imports
import unittest

class TestSettingsContainer(QuESoTestCase):
    def check_customized_values(self, settings):
        # Check general_settings
        general_settings = settings["general_settings"]

        self.assertTrue(general_settings.IsSet("input_filename"))
        input_filename = general_settings.GetString("input_filename")
        self.assertEqual(input_filename, "dummy.stl")

        self.assertTrue(general_settings.IsSet("output_directory_name"))
        output_directory_name = general_settings.GetString("output_directory_name")
        self.assertEqual(output_directory_name, "new_output")

        self.assertTrue(general_settings.IsSet("echo_level"))
        echo_level = general_settings.GetInt("echo_level")
        self.assertEqual(echo_level, 2)

        self.assertTrue(general_settings.IsSet("write_output_to_file"))
        write_output_to_file = general_settings.GetBool("write_output_to_file")
        self.assertFalse(write_output_to_file)

        # Check background_grid_settings
        background_grid_settings = settings["background_grid_settings"]

        self.assertTrue(background_grid_settings.IsSet("grid_type"))
        grid_type = background_grid_settings.GetGridType("grid_type")
        self.assertEqual(grid_type, QuESo.GridType.b_spline_grid)

        self.assertTrue(background_grid_settings.IsSet("lower_bound_xyz"))
        lower_bound_xyz = background_grid_settings.GetDoubleVector("lower_bound_xyz")
        self.assertListsAlmostEqual(lower_bound_xyz, [-130, -110, -110], 8 )

        self.assertTrue(background_grid_settings.IsSet("upper_bound_xyz"))
        upper_bound_xyz = background_grid_settings.GetDoubleVector("upper_bound_xyz")
        self.assertListsAlmostEqual(upper_bound_xyz, [20, 190, 190], 8 )

        self.assertTrue(background_grid_settings.IsSet("lower_bound_uvw"))
        lower_bound_uvw = background_grid_settings.GetDoubleVector("lower_bound_uvw")
        self.assertListsAlmostEqual(lower_bound_uvw, [1.23, 3.334, 5.66], 8 )

        self.assertTrue(background_grid_settings.IsSet("upper_bound_uvw"))
        upper_bound_uvw = background_grid_settings.GetDoubleVector("upper_bound_uvw")
        self.assertListsAlmostEqual(upper_bound_uvw, [4.4, 5.5, 2.2], 8 )

        self.assertTrue(background_grid_settings.IsSet("polynomial_order"))
        polynomial_order = background_grid_settings.GetIntVector("polynomial_order")
        self.assertListsEqual(polynomial_order, [2, 3, 2] )

        self.assertTrue(background_grid_settings.IsSet("number_of_elements"))
        number_of_elements = background_grid_settings.GetIntVector("number_of_elements")
        self.assertListsEqual(number_of_elements, [5, 2, 13] )

        # Check trimmed_quadrature_rule_settings
        trimmed_quadrature_rule_settings = settings["trimmed_quadrature_rule_settings"]

        self.assertTrue(trimmed_quadrature_rule_settings.IsSet("moment_fitting_residual"))
        moment_fitting_residual = trimmed_quadrature_rule_settings.GetDouble("moment_fitting_residual")
        self.assertAlmostEqual(moment_fitting_residual, 0.0023, 8)

        self.assertTrue(trimmed_quadrature_rule_settings.IsSet("min_element_volume_ratio"))
        min_element_volume_ratio = trimmed_quadrature_rule_settings.GetDouble("min_element_volume_ratio")
        self.assertAlmostEqual(min_element_volume_ratio, 0.012, 8)

        self.assertTrue(trimmed_quadrature_rule_settings.IsSet("min_num_boundary_triangles"))
        min_num_boundary_triangles = trimmed_quadrature_rule_settings.GetInt("min_num_boundary_triangles")
        self.assertEqual(min_num_boundary_triangles, 233)

        # Check non_trimmed_quadrature_rule_settings
        non_trimmed_quadrature_rule_settings = settings["non_trimmed_quadrature_rule_settings"]

        self.assertTrue(non_trimmed_quadrature_rule_settings.IsSet("integration_method"))
        integration_method = non_trimmed_quadrature_rule_settings.GetIntegrationMethod("integration_method")
        self.assertEqual(integration_method, QuESo.IntegrationMethod.GGQ_Optimal)

        conditions_settings_list = settings.GetList("conditions_settings_list")
        self.assertEqual(len(conditions_settings_list), 4)

        ref_condition_settings = {
        0:  {
                "condition_id" : 3,
                "condition_type" : "test_type1",
                "input_filename" : "test_filename_1.stl",
                "modulus" : 5.0,
                "direction" : [-1, 2, 3],
                "value" : [1, 2, 2],
                "penalty_factor" : 100000.0
            },
        1:  {
                "condition_id" : 1,
                "condition_type" : "test_type2",
                "input_filename" : "test_filename_2.stl",
                "modulus" : 2.0,
                "direction" : "Not set",
                "value" : "Not set",
                "penalty_factor" : "Not set"
            },
        2:  {
                "condition_id" : "Not set",
                "condition_type" : "Not set",
                "input_filename" : "test_filename_3.stl",
                "modulus" : "Not set",
                "direction" : "Not set",
                "value" : [0, 0.3, 0],
                "penalty_factor" : "Not set"
            },
        3:  {
                "condition_id" : "Not set",
                "condition_type" : "Not set",
                "input_filename" : "Not set",
                "modulus" : "Not set",
                "direction" : "Not set",
                "value" : [0, 0, 0],
                "penalty_factor" : 1e10
            }
        }

        for i, condition_setting in enumerate(conditions_settings_list):
            self.check_condition_customized_values(ref_condition_settings[i], condition_setting)
            self.check_condition_customized_values(ref_condition_settings[i], conditions_settings_list[i])

        condition_setting_1 = conditions_settings_list[0]

    def check_condition_customized_values(self, ref_settings, settings):
        if ref_settings["condition_id"] == "Not set":
            self.assertFalse(settings.IsSet("condition_id"))
        else:
            self.assertTrue(settings.IsSet("condition_id"))
            self.assertEqual(ref_settings["condition_id"], settings.GetInt("condition_id") )

        if ref_settings["condition_type"] == "Not set":
            self.assertFalse(settings.IsSet("condition_type"))
        else:
            self.assertTrue(settings.IsSet("condition_type"))
            self.assertEqual(ref_settings["condition_type"], settings.GetString("condition_type") )

        if ref_settings["input_filename"] == "Not set":
            self.assertFalse(settings.IsSet("input_filename"))
        else:
            self.assertTrue(settings.IsSet("input_filename"))
            self.assertEqual(ref_settings["input_filename"], settings.GetString("input_filename") )

        if ref_settings["modulus"] == "Not set":
                self.assertFalse(settings.IsSet("modulus"))
        else:
            self.assertTrue(settings.IsSet("modulus"))
            self.assertAlmostEqual(ref_settings["modulus"], settings.GetDouble("modulus"), 8 )

        if ref_settings["direction"] == "Not set":
            self.assertFalse(settings.IsSet("direction"))
        else:
            self.assertTrue(settings.IsSet("direction"))
            self.assertListsAlmostEqual(ref_settings["direction"], settings.GetDoubleVector("direction"), 8 )

        if ref_settings["value"] == "Not set":
            self.assertFalse(settings.IsSet("value"))
        else:
            self.assertTrue(settings.IsSet("value"))
            self.assertListsAlmostEqual(ref_settings["value"], settings.GetDoubleVector("value"), 8 )

        if ref_settings["penalty_factor"] == "Not set":
            self.assertFalse(settings.IsSet("penalty_factor"))
        else:
            self.assertTrue(settings.IsSet("penalty_factor"))
            self.assertAlmostEqual(ref_settings["penalty_factor"], settings.GetDouble("penalty_factor"), 8 )

    def check_default_values(self, settings):
        # Check general_settings
        general_settings = settings["general_settings"]

        self.assertFalse(general_settings.IsSet("input_filename"))

        self.assertTrue(general_settings.IsSet("output_directory_name"))
        output_directory_name = general_settings.GetString("output_directory_name")
        self.assertEqual(output_directory_name, "queso_output")

        self.assertTrue(general_settings.IsSet("echo_level"))
        echo_level = general_settings.GetInt("echo_level")
        self.assertEqual(echo_level, 1)

        self.assertTrue(general_settings.IsSet("write_output_to_file"))
        write_output_to_file = general_settings.GetBool("write_output_to_file")
        self.assertTrue(write_output_to_file)

        # Check background_grid_settings
        background_grid_settings = settings["background_grid_settings"]

        self.assertFalse(background_grid_settings.IsSet("grid_type"))

        self.assertFalse(background_grid_settings.IsSet("lower_bound_xyz"))

        self.assertFalse(background_grid_settings.IsSet("upper_bound_xyz"))

        self.assertFalse(background_grid_settings.IsSet("lower_bound_uvw"))

        self.assertFalse(background_grid_settings.IsSet("upper_bound_uvw"))

        self.assertFalse(background_grid_settings.IsSet("polynomial_order"))

        self.assertFalse(background_grid_settings.IsSet("number_of_elements"))

        # Check trimmed_quadrature_rule_settings
        trimmed_quadrature_rule_settings = settings["trimmed_quadrature_rule_settings"]

        self.assertTrue(trimmed_quadrature_rule_settings.IsSet("moment_fitting_residual"))
        moment_fitting_residual = trimmed_quadrature_rule_settings.GetDouble("moment_fitting_residual")
        self.assertAlmostEqual(moment_fitting_residual, 1e-10, 12)

        self.assertTrue(trimmed_quadrature_rule_settings.IsSet("min_element_volume_ratio"))
        min_element_volume_ratio = trimmed_quadrature_rule_settings.GetDouble("min_element_volume_ratio")
        self.assertAlmostEqual(min_element_volume_ratio, 0.001, 8)

        self.assertTrue(trimmed_quadrature_rule_settings.IsSet("min_num_boundary_triangles"))
        min_num_boundary_triangles = trimmed_quadrature_rule_settings.GetInt("min_num_boundary_triangles")
        self.assertEqual(min_num_boundary_triangles, 100)

        # Check non_trimmed_quadrature_rule_settings
        non_trimmed_quadrature_rule_settings = settings["non_trimmed_quadrature_rule_settings"]

        self.assertTrue(non_trimmed_quadrature_rule_settings.IsSet("integration_method"))
        integration_method = non_trimmed_quadrature_rule_settings.GetIntegrationMethod("integration_method")
        self.assertEqual(integration_method, QuESo.IntegrationMethod.Gauss)

        conditions_settings_list = settings.GetList("conditions_settings_list")
        self.assertEqual(len(conditions_settings_list), 4)

        for i, condition_setting in enumerate(conditions_settings_list):
            self.check_condition_default_values(condition_setting)
            self.check_condition_default_values(conditions_settings_list[i])

    def check_condition_default_values(self, settings):
        self.assertFalse(settings.IsSet("condition_id"))
        self.assertFalse(settings.IsSet("condition_type"))
        self.assertFalse(settings.IsSet("input_filename"))
        self.assertFalse(settings.IsSet("modulus"))
        self.assertFalse(settings.IsSet("direction"))
        self.assertFalse(settings.IsSet("value"))
        self.assertFalse(settings.IsSet("penalty_factor"))

    def test_customized_values(self):
        settings = JsonImport.ReadSettings("queso/tests/settings_container/QuESoSettings_custom_1.json")
        self.check_customized_values(settings)
        settings2 = JsonImport.ReadSettings("queso/tests/settings_container/QuESoSettings_custom_2.json")
        self.check_customized_values(settings2)

    def test_default_values(self):
        settings = JsonImport.ReadSettings("queso/tests/settings_container/QuESoSettings_default.json")
        self.check_default_values(settings)


if __name__ == "__main__":
    unittest.main()

