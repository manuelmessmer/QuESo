# Project imports
import QuESoPythonModule as QuESo_App
from QuESoPythonModule.kratos_interface.kratos_analysis import KratosAnalysis
from QuESoPythonModule.model import Model
from QuESoPythonModule.scripts.helper import *
from QuESoPythonModule.kratos_interface.model_part_utilities import add_conditions
from QuESoPythonModule.scripts.queso_unit_test import QuESoTestCase
# Kratos imports
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
# Unittest imports
import unittest

class TestBoundaryConditionsKratos(QuESoTestCase):
    def run_modelers(self, current_model, modelers_list):
        from KratosMultiphysics.modeler_factory import KratosModelerFactory
        factory = KratosModelerFactory()
        list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

        for modeler in list_of_modelers:
            modeler.SetupGeometryModel()

        for modeler in list_of_modelers:
            modeler.PrepareGeometryModel()

        for modeler in list_of_modelers:
            modeler.SetupModelPart()

    def check_surface_area(self, model_part, ref_value):
        area = 0.0
        for condition in model_part.Conditions:
            det_J = condition.GetGeometry().DeterminantOfJacobian()[0]
            weight = 0.5 # For triangles
            area += weight*det_J
        self.assertAlmostEqual(area,ref_value, 5)

    def construct_b_spline_volume(self, model, QuESoSettings):
        modeler_settings = KM.Parameters("""
            [{
                "modeler_name": "NurbsGeometryModeler",
                "Parameters": {
                    "model_part_name" : "NurbsMesh",
                    "geometry_name"   : "NurbsVolume"
                }
            }]
            """)

        grid_settings = QuESoSettings["background_grid_settings"]

        tmp_parameters = modeler_settings[0]["Parameters"]
        tmp_parameters.AddEmptyValue("lower_point_xyz")
        tmp_parameters["lower_point_xyz"].SetVector(grid_settings["lower_bound_xyz"])
        tmp_parameters.AddEmptyValue("upper_point_xyz")
        tmp_parameters["upper_point_xyz"].SetVector(grid_settings["upper_bound_xyz"])

        tmp_parameters.AddEmptyValue("lower_point_uvw")
        tmp_parameters["lower_point_uvw"].SetVector(grid_settings["lower_bound_uvw"])
        tmp_parameters.AddEmptyValue("upper_point_uvw")
        tmp_parameters["upper_point_uvw"].SetVector(grid_settings["upper_bound_uvw"])

        tmp_parameters.AddEmptyValue("polynomial_order")
        tmp_parameters["polynomial_order"].SetVector(grid_settings["polynomial_order"])
        tmp_parameters.AddEmptyValue("number_of_knot_spans")
        tmp_parameters["number_of_knot_spans"].SetVector(grid_settings["number_of_elements"])

        self.run_modelers(model, modeler_settings)

    def get_active_condition_ids(self, pyqueso, component_name):
        return [
            condition.GetSettings().GetInt("condition_id")
            for condition in pyqueso.conditions(component_name)
        ]

    def get_property_id_by_condition_id(self, pyqueso, model_part, component_name):
        property_id_by_condition_id = {}
        next_property_id = model_part.GetRootModelPart().NumberOfProperties() + 1
        for condition_id in self.get_active_condition_ids(pyqueso, component_name):
            while model_part.GetRootModelPart().HasProperties(next_property_id):
                next_property_id += 1
            model_part.CreateNewProperties(next_property_id)
            property_id_by_condition_id[int(condition_id)] = next_property_id
            next_property_id += 1
        return property_id_by_condition_id

    def get_total_force(self, model_part):
        force = [0.0, 0.0, 0.0]
        for condition in model_part.Conditions:
            force[0] += condition.GetValue(SMA.POINT_LOAD_X)
            force[1] += condition.GetValue(SMA.POINT_LOAD_Y)
            force[2] += condition.GetValue(SMA.POINT_LOAD_Z)
        return force

    def test_penalty_support(self):
        pyqueso = Model("queso/tests/boundary_conditions_kratos/QuESoSettings_Penalty.json")
        pyqueso.run()

        model = KM.Model()
        model_part = model.CreateModelPart("NurbsMesh")
        model_part.CreateSubModelPart("main")

        settings = pyqueso.settings("main")
        grid_settings = settings["background_grid_settings"]
        bounds_xyz = [grid_settings["lower_bound_xyz"], grid_settings["upper_bound_xyz"]]
        bounds_uvw = [grid_settings["lower_bound_uvw"], grid_settings["upper_bound_uvw"]]
        self.construct_b_spline_volume(model, settings)
        add_conditions(
            geometry_model_part=model_part,
            component_model_part=model_part.GetSubModelPart("main"),
            conditions=pyqueso.conditions("main"),
            active_condition_ids=self.get_active_condition_ids(pyqueso, "main"),
            property_id_by_condition_id=self.get_property_id_by_condition_id(pyqueso, model_part, "main"),
            bounds_xyz=bounds_xyz,
            bounds_uvw=bounds_uvw,
        )

        properties = model_part.GetProperties()[1]
        process_info = KM.ProcessInfo()
        penalty = properties.GetValue(IGA.PENALTY_FACTOR)
        self.assertAlmostEqual(penalty, 1e10, 10)

        for condition in model_part.Conditions:
            value = condition.GetValue(KM.DISPLACEMENT)
            self.assertListsAlmostEqual(value, [0.0, 0.0, 1.0], 10)

        self.check_surface_area(model_part, 1183.54304)

    def test_lagrange_support(self):
        pyqueso = Model("queso/tests/boundary_conditions_kratos/QuESoSettings_Lagrange.json")
        pyqueso.run()

        model = KM.Model()
        model_part = model.CreateModelPart("NurbsMesh")
        model_part.CreateSubModelPart("main")

        settings = pyqueso.settings("main")
        grid_settings = settings["background_grid_settings"]
        bounds_xyz = [grid_settings["lower_bound_xyz"], grid_settings["upper_bound_xyz"]]
        bounds_uvw = [grid_settings["lower_bound_uvw"], grid_settings["upper_bound_uvw"]]
        self.construct_b_spline_volume(model, settings)
        add_conditions(
            geometry_model_part=model_part,
            component_model_part=model_part.GetSubModelPart("main"),
            conditions=pyqueso.conditions("main"),
            active_condition_ids=self.get_active_condition_ids(pyqueso, "main"),
            property_id_by_condition_id=self.get_property_id_by_condition_id(pyqueso, model_part, "main"),
            bounds_xyz=bounds_xyz,
            bounds_uvw=bounds_uvw,
        )

        for condition in model_part.Conditions:
            value = condition.GetValue(KM.DISPLACEMENT)
            self.assertListsAlmostEqual(value, [0.0, 0.3, 0.0], 10)

        self.check_surface_area(model_part, 921.163635)

    def test_surface_load(self):
        pyqueso = Model("queso/tests/boundary_conditions_kratos/QuESoSettings_SurfaceLoad.json")
        pyqueso.run()

        model = KM.Model()
        model_part = model.CreateModelPart("NurbsMesh")
        model_part.CreateSubModelPart("main")

        settings = pyqueso.settings("main")
        grid_settings = settings["background_grid_settings"]
        bounds_xyz = [grid_settings["lower_bound_xyz"], grid_settings["upper_bound_xyz"]]
        bounds_uvw = [grid_settings["lower_bound_uvw"], grid_settings["upper_bound_uvw"]]
        self.construct_b_spline_volume(model, settings)
        add_conditions(
            geometry_model_part=model_part,
            component_model_part=model_part.GetSubModelPart("main"),
            conditions=pyqueso.conditions("main"),
            active_condition_ids=self.get_active_condition_ids(pyqueso, "main"),
            property_id_by_condition_id=self.get_property_id_by_condition_id(pyqueso, model_part, "main"),
            bounds_xyz=bounds_xyz,
            bounds_uvw=bounds_uvw,
        )

        force = self.get_total_force(model_part)
        ref_value = 9594.9131
        self.assertListsAlmostEqual(force, [ref_value] * 3, 4)

    def test_total_load(self):
        pyqueso = Model("queso/tests/boundary_conditions_kratos/QuESoSettings_TotalLoad.json")
        pyqueso.run()

        model = KM.Model()
        model_part = model.CreateModelPart("NurbsMesh")
        model_part.CreateSubModelPart("main")

        settings = pyqueso.settings("main")
        grid_settings = settings["background_grid_settings"]
        bounds_xyz = [grid_settings["lower_bound_xyz"], grid_settings["upper_bound_xyz"]]
        bounds_uvw = [grid_settings["lower_bound_uvw"], grid_settings["upper_bound_uvw"]]
        self.construct_b_spline_volume(model, settings)
        created_loads = add_conditions(
            geometry_model_part=model_part,
            component_model_part=model_part.GetSubModelPart("main"),
            conditions=pyqueso.conditions("main"),
            active_condition_ids=self.get_active_condition_ids(pyqueso, "main"),
            property_id_by_condition_id=self.get_property_id_by_condition_id(pyqueso, model_part, "main"),
            bounds_xyz=bounds_xyz,
            bounds_uvw=bounds_uvw,
        )

        self.assertGreater(len(created_loads), 0)
        force = self.get_total_force(model_part)
        ref_value = 50.0 / (3.0 ** 0.5)
        self.assertListsAlmostEqual(force, [ref_value] * 3, 5)

    def test_pressure_load(self):
        pyqueso = Model("queso/tests/boundary_conditions_kratos/QuESoSettings_Pressure.json")
        pyqueso.run()

        model = KM.Model()
        model_part = model.CreateModelPart("NurbsMesh")
        model_part.CreateSubModelPart("main")

        settings = pyqueso.settings("main")
        grid_settings = settings["background_grid_settings"]
        bounds_xyz = [grid_settings["lower_bound_xyz"], grid_settings["upper_bound_xyz"]]
        bounds_uvw = [grid_settings["lower_bound_uvw"], grid_settings["upper_bound_uvw"]]
        self.construct_b_spline_volume(model, settings)
        add_conditions(
            geometry_model_part=model_part,
            component_model_part=model_part.GetSubModelPart("main"),
            conditions=pyqueso.conditions("main"),
            active_condition_ids=self.get_active_condition_ids(pyqueso, "main"),
            property_id_by_condition_id=self.get_property_id_by_condition_id(pyqueso, model_part, "main"),
            bounds_xyz=bounds_xyz,
            bounds_uvw=bounds_uvw,
        )

        force = self.get_total_force(model_part)
        ref_value = -577.978141*2
        self.assertListsAlmostEqual(force, [0.0, 0.0, ref_value], 5)

    def test_total_load_ramp_in_kratos_analysis(self):
        analysis = KratosAnalysis(
            KM.Model(),
            queso_settings_filename="queso/tests/boundary_conditions_kratos/QuESoSettings_TotalLoad.json",
            analysis_settings_filename="queso/tests/boundary_conditions_kratos/AnalysisSettings.json",
            kratos_settings_filename="queso/tests/boundary_conditions_kratos/KratosParameters.json",
        )
        analysis.Initialize()

        model_part = analysis.GetModelPart()
        root_model_part = model_part.GetRootModelPart()
        root_model_part.ProcessInfo[KM.STEP] = 1
        analysis.ApplyBoundaryConditions()
        step_1_force = self.get_total_force(model_part)

        root_model_part.ProcessInfo[KM.STEP] = 2
        analysis.ApplyBoundaryConditions()
        step_2_force = self.get_total_force(model_part)

        root_model_part.ProcessInfo[KM.STEP] = 5
        analysis.ApplyBoundaryConditions()
        step_5_force = self.get_total_force(model_part)

        ref_value = 50.0 / (3.0 ** 0.5)
        self.assertListsAlmostEqual(step_1_force, [ref_value / 5.0] * 3, 5)
        self.assertListsAlmostEqual(step_2_force, [2.0 * ref_value / 5.0] * 3, 5)
        self.assertListsAlmostEqual(step_5_force, [ref_value] * 3, 5)

    def test_fixed_model_part_in_kratos_analysis(self):
        analysis = KratosAnalysis(
            KM.Model(),
            queso_settings_filename="queso/tests/boundary_conditions_kratos/QuESoSettings_TotalLoad.json",
            analysis_settings_filename="queso/tests/boundary_conditions_kratos/AnalysisSettings_FixedModelPart.json",
            kratos_settings_filename="queso/tests/boundary_conditions_kratos/KratosParameters.json",
        )
        analysis.Initialize()

        model_part = analysis.GetModelPart()
        component_model_part = model_part.GetSubModelPart("main")
        self.assertGreater(component_model_part.NumberOfNodes(), 0)

        for node in component_model_part.Nodes:
            self.assertTrue(node.IsFixed(KM.DISPLACEMENT_X))
            self.assertFalse(node.IsFixed(KM.DISPLACEMENT_Y))
            self.assertFalse(node.IsFixed(KM.DISPLACEMENT_Z))

if __name__ == "__main__":
    unittest.main()
