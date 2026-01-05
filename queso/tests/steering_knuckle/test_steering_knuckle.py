# Project imports
from QuESoPythonModule.PyQuESo import PyQuESo
from QuESoPythonModule.scripts.queso_unit_test import QuESoTestCase

# External imports
import unittest
import shutil
import os
import json

class TestStrainEnergySteeringKnuckleKratos(QuESoTestCase):
    def run_test(self, filename, tolerance):
        pyqueso = PyQuESo(filename)
        pyqueso.Run()

        json_filename = "queso/tests/steering_knuckle/output/model_info.json"
        # Note: Precision of json_dict is slightly lower. Double are written with std::setprecision n=10.
        with open(json_filename, 'r') as file:
            json_dict = json.load(file)

        volume = 0.0
        for element in pyqueso.GetElements():
            for point in element.GetIntegrationPoints():
                volume += point.Weight()
        model_info = pyqueso.GetModelInfo()
        settings = pyqueso.GetSettings()

        volume_info = model_info["quadrature_info"].GetDouble("represented_volume")
        self.assertAlmostEqual(volume, volume_info, places=7)

        ## Check model info
        # embedded_geometry_info
        self.assertEqual(model_info["embedded_geometry_info"].GetBool("is_closed"), True)
        self.assertAlmostEqual(model_info["embedded_geometry_info"].GetDouble("volume"), volume, places=5)

        self.assertEqual(json_dict["embedded_geometry_info"]["is_closed"], True)
        self.assertAlmostEqual(json_dict["embedded_geometry_info"]["volume"], volume, places=4)
        # quadrature_info
        self.assertAlmostEqual(model_info["quadrature_info"].GetDouble("percentage_of_geometry_volume"), 100.0, places=5)
        self.assertEqual(model_info["quadrature_info"].GetInt("tot_num_points"), 13251)
        self.assertAlmostEqual(model_info["quadrature_info"].GetDouble("num_of_points_per_full_element"), 23.25, places=5)
        self.assertGreater(model_info["quadrature_info"].GetDouble("num_of_points_per_trimmed_element"), 26)
        self.assertLess(model_info["quadrature_info"].GetDouble("num_of_points_per_trimmed_element"), 27)

        self.assertAlmostEqual(json_dict["quadrature_info"]["percentage_of_geometry_volume"], 100.0, places=5)
        self.assertEqual(json_dict["quadrature_info"]["tot_num_points"], 13251)
        self.assertAlmostEqual(json_dict["quadrature_info"]["num_of_points_per_full_element"], 23.25, places=5)
        self.assertGreater(json_dict["quadrature_info"]["num_of_points_per_trimmed_element"], 26)
        self.assertLess(json_dict["quadrature_info"]["num_of_points_per_trimmed_element"], 27)
        # background_grid_info
        self.assertEqual(model_info["background_grid_info"].GetInt("num_active_elements"), 498)
        self.assertEqual(model_info["background_grid_info"].GetInt("num_trimmed_elements"), 486)
        self.assertEqual(model_info["background_grid_info"].GetInt("num_full_elements"), 12)
        self.assertEqual(model_info["background_grid_info"].GetInt("num_inactive_elements"), 5752)

        self.assertEqual(json_dict["background_grid_info"]["num_active_elements"], 498)
        self.assertEqual(json_dict["background_grid_info"]["num_trimmed_elements"], 486)
        self.assertEqual(json_dict["background_grid_info"]["num_full_elements"], 12)
        self.assertEqual(json_dict["background_grid_info"]["num_inactive_elements"], 5752)
        # elapsed_time_info
        self.assertGreater(model_info["elapsed_time_info"].GetDouble("total"), 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["total"], 0.0)

        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].GetDouble("total"), 0.0)
        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].GetDouble("classification_of_elements"), 0.0)
        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].GetDouble("computation_of_intersections"), 0.0)
        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].GetDouble("solution_of_moment_fitting_eqs"), 0.0)
        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].GetDouble("construction_of_ggq_rules"), 0.0)

        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["total"], 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["classification_of_elements"], 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["computation_of_intersections"], 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["solution_of_moment_fitting_eqs"], 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["construction_of_ggq_rules"], 0.0)

        self.assertGreater(model_info["elapsed_time_info"]["conditions_time_info"].GetDouble("total"), 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["conditions_time_info"]["total"], 0.0)

        self.assertGreater(model_info["elapsed_time_info"]["write_files_time_info"].GetDouble("total"), 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["write_files_time_info"]["total"], 0.0)

        # conditions_infos_list
        conditions_info_list = model_info.GetList("conditions_infos_list")
        conditions_info_list_json = json_dict["conditions_infos_list"]

        surf_areas_ref = [332.3775399, 577.9781408, 921.1636352, 1183.543044]
        for i, (info, info_json) in enumerate(zip(conditions_info_list, conditions_info_list_json)):
            self.assertEqual(info.GetInt("condition_id"), (i+1) )
            self.assertEqual(info_json["condition_id"], (i+1) )
            self.assertAlmostEqual(info.GetDouble("surf_area"), surf_areas_ref[i], places=5 )
            self.assertAlmostEqual(info_json["surf_area"], surf_areas_ref[i], places=5 )
            self.assertAlmostEqual(info.GetDouble("perc_surf_area_in_active_domain"), 100.0, places=5 )
            self.assertAlmostEqual(info_json["perc_surf_area_in_active_domain"], 100.0, places=5 )

    def test_1(self):
        self.run_test("queso/tests/steering_knuckle/QuESoSettings1.json", 0.005)

    def tearDown(self):
        dir_name = "queso/tests/steering_knuckle/output"
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == "__main__":
    unittest.main()



