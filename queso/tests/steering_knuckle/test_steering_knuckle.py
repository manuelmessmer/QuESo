# Project imports
import json
import os
import shutil

# External imports
import unittest

import pyqueso
from pyqueso.scripts.queso_unit_test import QuESoTestCase


class TestStrainEnergySteeringKnuckleKratos(QuESoTestCase):
    def run_test(self, filename, tolerance):
        model = pyqueso.Model(filename)
        model.create()

        json_filename = "queso/tests/steering_knuckle/output/main/component_info.json"
        # Note: Precision of json_dict is slightly lower. Double are written with std::setprecision n=10.
        with open(json_filename, 'r') as file:
            json_dict = json.load(file)

        volume = 0.0
        for element in model.elements("main"):
            for point in element.integration_points:
                volume += point.weight
        model_info = model.component_info("main")

        volume_info = model_info["quadrature_info"].get_double("represented_volume")
        self.assertAlmostEqual(volume, volume_info, places=7)

        ## Check model info
        # embedded_geometry_info
        self.assertEqual(model_info["embedded_geometry_info"].get_bool("is_closed"), True)
        self.assertAlmostEqual(model_info["embedded_geometry_info"].get_double("volume"), volume, places=5)

        self.assertEqual(json_dict["embedded_geometry_info"]["is_closed"], True)
        self.assertAlmostEqual(json_dict["embedded_geometry_info"]["volume"], volume, places=4)
        # quadrature_info
        self.assertAlmostEqual(model_info["quadrature_info"].get_double("percentage_of_geometry_volume"), 100.0, places=5)
        self.assertEqual(model_info["quadrature_info"].get_int("tot_num_points"), 13252)
        self.assertAlmostEqual(model_info["quadrature_info"].get_double("num_of_points_per_full_element"), 23.25, places=5)
        self.assertGreater(model_info["quadrature_info"].get_double("num_of_points_per_trimmed_element"), 26)
        self.assertLess(model_info["quadrature_info"].get_double("num_of_points_per_trimmed_element"), 27)

        self.assertAlmostEqual(json_dict["quadrature_info"]["percentage_of_geometry_volume"], 100.0, places=5)
        self.assertEqual(json_dict["quadrature_info"]["tot_num_points"], 13252)
        self.assertAlmostEqual(json_dict["quadrature_info"]["num_of_points_per_full_element"], 23.25, places=5)
        self.assertGreater(json_dict["quadrature_info"]["num_of_points_per_trimmed_element"], 26)
        self.assertLess(json_dict["quadrature_info"]["num_of_points_per_trimmed_element"], 27)
        # background_grid_info
        self.assertEqual(model_info["background_grid_info"].get_int("num_active_elements"), 498)
        self.assertEqual(model_info["background_grid_info"].get_int("num_trimmed_elements"), 486)
        self.assertEqual(model_info["background_grid_info"].get_int("num_full_elements"), 12)
        self.assertEqual(model_info["background_grid_info"].get_int("num_inactive_elements"), 5752)

        self.assertEqual(json_dict["background_grid_info"]["num_active_elements"], 498)
        self.assertEqual(json_dict["background_grid_info"]["num_trimmed_elements"], 486)
        self.assertEqual(json_dict["background_grid_info"]["num_full_elements"], 12)
        self.assertEqual(json_dict["background_grid_info"]["num_inactive_elements"], 5752)
        # elapsed_time_info
        self.assertGreater(model_info["elapsed_time_info"].get_double("total"), 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["total"], 0.0)

        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].get_double("total"), 0.0)
        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].get_double("classification_of_elements"), 0.0)
        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].get_double("computation_of_intersections"), 0.0)
        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].get_double("solution_of_moment_fitting_eqs"), 0.0)
        self.assertGreater(model_info["elapsed_time_info"]["volume_time_info"].get_double("construction_of_ggq_rules"), 0.0)

        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["total"], 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["classification_of_elements"], 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["computation_of_intersections"], 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["solution_of_moment_fitting_eqs"], 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["volume_time_info"]["construction_of_ggq_rules"], 0.0)

        self.assertGreater(model_info["elapsed_time_info"]["conditions_time_info"].get_double("total"), 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["conditions_time_info"]["total"], 0.0)

        self.assertGreater(model_info["elapsed_time_info"]["write_files_time_info"].get_double("total"), 0.0)
        self.assertGreater(json_dict["elapsed_time_info"]["write_files_time_info"]["total"], 0.0)

        # conditions_infos_list
        conditions_info_list = model_info.get_list("conditions_infos_list")
        conditions_info_list_json = json_dict["conditions_infos_list"]

        surf_areas_ref = [332.3775399, 577.9781408, 921.1636352, 1183.543044]
        active_area_percentages_ref = [100.0, 100.0, 100.0, 100.0]
        for i, (info, info_json) in enumerate(zip(conditions_info_list, conditions_info_list_json)):
            self.assertEqual(info.get_int("condition_id"), (i+1) )
            self.assertEqual(info_json["condition_id"], (i+1) )
            self.assertAlmostEqual(info.get_double("surf_area"), surf_areas_ref[i], places=5 )
            self.assertAlmostEqual(info_json["surf_area"], surf_areas_ref[i], places=5 )
            self.assertAlmostEqual(info.get_double("perc_surf_area_in_active_domain"), active_area_percentages_ref[i], places=5 )
            self.assertAlmostEqual(info_json["perc_surf_area_in_active_domain"], active_area_percentages_ref[i], places=5 )

    def test_1(self):
        self.run_test("queso/tests/steering_knuckle/QuESoSettings1.json", 0.005)

    def tearDown(self):
        dir_name = "queso/tests/steering_knuckle/output"
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == "__main__":
    unittest.main()
