# Project imports
import QuESo_PythonApplication as QuESo_App
from QuESo_PythonApplication.PyQuESo import PyQuESo
from queso.python_scripts.helper import *
from kratos_interface.model_part_utilities import ModelPartUtilities
from queso.python_scripts.QuESoUnittest import QuESoTestCase
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
        tmp_parameters["lower_point_xyz"].SetVector(grid_settings.GetDoubleVector("lower_bound_xyz"))
        tmp_parameters.AddEmptyValue("upper_point_xyz")
        tmp_parameters["upper_point_xyz"].SetVector(grid_settings.GetDoubleVector("upper_bound_xyz"))

        tmp_parameters.AddEmptyValue("lower_point_uvw")
        tmp_parameters["lower_point_uvw"].SetVector(grid_settings.GetDoubleVector("lower_bound_uvw"))
        tmp_parameters.AddEmptyValue("upper_point_uvw")
        tmp_parameters["upper_point_uvw"].SetVector(grid_settings.GetDoubleVector("upper_bound_uvw"))

        tmp_parameters.AddEmptyValue("polynomial_order")
        tmp_parameters["polynomial_order"].SetVector(grid_settings.GetIntVector("polynomial_order"))
        tmp_parameters.AddEmptyValue("number_of_knot_spans")
        tmp_parameters["number_of_knot_spans"].SetVector(grid_settings.GetIntVector("number_of_elements"))

        self.run_modelers(model, modeler_settings)

    def test_penalty_support(self):
        pyqueso = PyQuESo("queso/tests/boundary_conditions_kratos/QuESoSettings_Penalty.json")
        pyqueso.Run()

        model = KM.Model()
        model_part = model.CreateModelPart("NurbsMesh")

        settings = pyqueso.GetSettings()
        grid_settings = settings["background_grid_settings"]
        bounds_xyz = [grid_settings.GetDoubleVector("lower_bound_xyz"),
                      grid_settings.GetDoubleVector("upper_bound_xyz") ]
        bounds_uvw = [grid_settings.GetDoubleVector("lower_bound_uvw"),
                      grid_settings.GetDoubleVector("upper_bound_uvw") ]
        self.construct_b_spline_volume(model, settings)
        ModelPartUtilities.AddConditionsToModelPart(model_part, pyqueso.GetConditions(), bounds_xyz, bounds_uvw)

        properties = model_part.GetProperties()[1]
        process_info = KM.ProcessInfo()
        penalty = properties.GetValue(IGA.PENALTY_FACTOR)
        self.assertAlmostEqual(penalty, 1e10, 10)

        for condition in model_part.Conditions:
            value = condition.GetValue(KM.DISPLACEMENT)
            self.assertListsAlmostEqual(value, [0.0, 0.0, 1.0], 10)

        self.check_surface_area(model_part, 1183.54304)

    def test_lagrange_support(self):
        pyqueso = PyQuESo("queso/tests/boundary_conditions_kratos/QuESoSettings_Lagrange.json")
        pyqueso.Run()

        model = KM.Model()
        model_part = model.CreateModelPart("NurbsMesh")

        settings = pyqueso.GetSettings()
        grid_settings = settings["background_grid_settings"]
        bounds_xyz = [grid_settings.GetDoubleVector("lower_bound_xyz"),
                      grid_settings.GetDoubleVector("upper_bound_xyz") ]
        bounds_uvw = [grid_settings.GetDoubleVector("lower_bound_uvw"),
                      grid_settings.GetDoubleVector("upper_bound_uvw") ]
        self.construct_b_spline_volume(model, settings)
        ModelPartUtilities.AddConditionsToModelPart(model_part, pyqueso.GetConditions(), bounds_xyz, bounds_uvw)

        for condition in model_part.Conditions:
            value = condition.GetValue(KM.DISPLACEMENT)
            self.assertListsAlmostEqual(value, [0.0, 0.3, 0.0], 10)

        self.check_surface_area(model_part, 921.163635)

    def test_surface_load(self):
        pyqueso = PyQuESo("queso/tests/boundary_conditions_kratos/QuESoSettings_SurfaceLoad.json")
        pyqueso.Run()

        model = KM.Model()
        model_part = model.CreateModelPart("NurbsMesh")

        settings = pyqueso.GetSettings()
        grid_settings = settings["background_grid_settings"]
        bounds_xyz = [grid_settings.GetDoubleVector("lower_bound_xyz"),
                      grid_settings.GetDoubleVector("upper_bound_xyz") ]
        bounds_uvw = [grid_settings.GetDoubleVector("lower_bound_uvw"),
                      grid_settings.GetDoubleVector("upper_bound_uvw") ]
        self.construct_b_spline_volume(model, settings)
        ModelPartUtilities.AddConditionsToModelPart(model_part, pyqueso.GetConditions(), bounds_xyz, bounds_uvw)

        force = [0.0, 0.0, 0.0]
        for condition in model_part.Conditions:
            value = condition.GetValue(SMA.POINT_LOAD)
            force[0] += value[0]
            force[1] += value[1]
            force[2] += value[2]
        ref_value = 959.49131
        self.assertListsAlmostEqual(force, [ref_value]*3, 5)

    def test_pressure_load(self):
        pyqueso = PyQuESo("queso/tests/boundary_conditions_kratos/QuESoSettings_Pressure.json")
        pyqueso.Run()

        model = KM.Model()
        model_part = model.CreateModelPart("NurbsMesh")

        settings = pyqueso.GetSettings()
        grid_settings = settings["background_grid_settings"]
        bounds_xyz = [grid_settings.GetDoubleVector("lower_bound_xyz"),
                      grid_settings.GetDoubleVector("upper_bound_xyz") ]
        bounds_uvw = [grid_settings.GetDoubleVector("lower_bound_uvw"),
                      grid_settings.GetDoubleVector("upper_bound_uvw") ]
        self.construct_b_spline_volume(model, settings)
        ModelPartUtilities.AddConditionsToModelPart(model_part, pyqueso.GetConditions(), bounds_xyz, bounds_uvw)

        force = [0.0, 0.0, 0.0]
        for condition in model_part.Conditions:
            value = condition.GetValue(SMA.POINT_LOAD)
            force[0] += value[0]
            force[1] += value[1]
            force[2] += value[2]
        ref_value = -577.978141*2
        self.assertListsAlmostEqual(force, [0.0, 0.0, ref_value], 5)

if __name__ == "__main__":
    unittest.main()