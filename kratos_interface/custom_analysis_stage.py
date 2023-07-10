# Import Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.IgaApplication.map_nurbs_volume_results_to_embedded_geometry_process import MapNurbsVolumeResultsToEmbeddedGeometryProcess
#from KratosMultiphysics.IgaApplication.assign_integration_points_to_background_elements_process import AssignIntegrationPointsToBackgroundElementsProcess
from kratos_interface.model_part_io import *
from tibra.python_scripts.helper import *

class CustomAnalysisStage(StructuralMechanicsAnalysis):
    """Customized Kratos Analysis Stage.

    Overrides the StructuralMechanicsAnalysis Stage from Kratos.
    """
    def __init__(self, model, tibra_parameters, kratos_settings_filename, elements, boundary_conditions, triangle_mesh):
        """The constructor."""
        # Read kratos settings
        with open(kratos_settings_filename,'r') as parameter_file:
            analysis_parameters = KM.Parameters(parameter_file.read())

        self.boundary_conditions = boundary_conditions
        self.elements = elements
        self.triangle_mesh = triangle_mesh
        self.lower_bound = tibra_parameters.LowerBound()
        self.upper_bound = tibra_parameters.UpperBound()
        #Override the NurbsGeometryModeler input parameters
        for modeler in analysis_parameters["modelers"].values():
            if modeler["modeler_name"].GetString() == "NurbsGeometryModeler":
                parameters = modeler["Parameters"]
                parameters.AddEmptyValue("lower_point")
                parameters["lower_point"].SetVector(tibra_parameters.LowerBound())
                parameters.AddEmptyValue("upper_point")
                parameters["upper_point"].SetVector(tibra_parameters.UpperBound())
                parameters.AddEmptyValue("polynomial_order")
                parameters["polynomial_order"].SetVector(tibra_parameters.Order())
                parameters.AddEmptyValue("number_of_knot_spans")
                parameters["number_of_knot_spans"].SetVector(tibra_parameters.NumberOfElements())

        self.Initialized = False
        super().__init__(model, analysis_parameters)


    def _ModelersSetupModelPart(self):
        """Override BaseClass to run NURBS modelers."""
        model_part = self.model.GetModelPart('NurbsMesh')
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)

        model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE,3)

        embedded_model_part = self.model.CreateModelPart("EmbeddedModelPart")

        embedded_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        embedded_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE,3)
        WriteKratosModelPart(self.triangle_mesh, embedded_model_part)

        embedded_model_part2 = self.model.CreateModelPart("EmbeddedModelPart2")
        embedded_model_part2.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        WriteKratosModelPart(self.triangle_mesh, embedded_model_part2)


        # Convert the geometry model or import analysis suitable models.
        for modeler in self._GetListOfModelers():
            print(modeler)
            if self.echo_level > 1:
                KM.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart started.")
            modeler.SetupModelPart()
            if self.echo_level > 1:
                KM.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart finished.")

        return super()._ModelersSetupModelPart()


    def ModifyInitialGeometry(self):
        """Override BaseClass to pass integration points to Kratos."""
        model_part = self.model.GetModelPart('NurbsMesh')
        nurbs_volume = model_part.GetGeometry("NurbsVolume")
        volume_properties = model_part.GetProperties()[1]

        el_count = 0
        for element in self.elements:
            integration_points = []
            if element.IsTrimmed():
                for point in element.GetIntegrationPoints():
                    weight = point.GetWeight()
                    if( weight > 0):
                        integration_points.append([point.GetX(), point.GetY(), point.GetZ(), point.GetWeight()])
            else:
                for point in element.GetIntegrationPoints():
                    integration_points.append([point.GetX(), point.GetY(), point.GetZ(), point.GetWeight()])

            if( len(integration_points) > 0 ):
                el_count += 1
                # Create quadrature_point_geometries
                # These are basically just integration points.
                quadrature_point_geometries = KM.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2, integration_points)
                el = model_part.CreateNewElement('SmallDisplacementElement3D8N', el_count, quadrature_point_geometries[0], volume_properties)

        print("Number of Elements/Integration Points (In BSpline Volume): ", model_part.NumberOfElements())

        # Add Dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

    # def OutputSolutionStep(self):

    #     execute_was_called = False
    #     for output_process in self._GetListOfOutputProcesses():
    #         if output_process.IsOutputStep():
    #             if not execute_was_called:
    #                 for process in self._GetListOfProcesses():
    #                     process.ExecuteBeforeOutputStep()
    #                 execute_was_called = True

    #     #         output_process.PrintOutput()

    #     # if execute_was_called:
    #     #     for process in self._GetListOfProcesses():
    #     #         process.ExecuteAfterOutputStep()


    #     embedded_model_part = self.model.GetModelPart("EmbeddedModelPart")
    #     for element in embedded_model_part.Elements:
    #         #print("0")
    #         print(element.GetGeometry())
    #         values = element.CalculateOnIntegrationPoints(KM.STRAIN_ENERGY, embedded_model_part.ProcessInfo)
    #         #print("1")

        #super().OutputSolutionStep()
        #return super().OutputSolutionStep()
    #     model_part = self.model.GetModelPart('NurbsMesh')
    #     embedded_model_part = self.model.GetModelPart("EmbeddedModelPart")
    #     embedded_model_part.SetProperties(model_part.GetProperties())
    #     return super().OutputSolutionStep()
    #     #  Map nodal values
    #     process_params = KM.Parameters(
    #     """ {
    #             "main_model_part_name"                    : "NurbsMesh",
    #             "nurbs_volume_name"                       : "NurbsVolume",
    #             "embedded_model_part_name"                : "EmbeddedModelPart",
    #             "nodal_results": ["DISPLACEMENT"]
    #     } """ )
    #     process = MapNurbsVolumeResultsToEmbeddedGeometryProcess(self.model, process_params)
    #     process.ExecuteBeforeOutputStep()

    #     # point_process_params = KM.Parameters(
    #     # """ {
    #     #         "main_model_part_name"                    : "NurbsMesh",
    #     #         "nurbs_volume_name"                       : "NurbsVolume",
    #     #         "embedded_model_part_name"                : "EmbeddedModelPart"
    #     # } """ )
    #     # point_process = KM.AssignIntegrationPointsToBackgroundElementsProcess(model, point_process_params)
    #     # point_process.ExecuteBeforeOutputStep()

    #     # Map integration point values
    #     # model_part = self.model.GetModelPart("NurbsMesh")
    #     # embedded_model_part = self.model.GetModelPart("EmbeddedModelPart")

    #     # for element in embedded_model_part.Elements:
    #     #     center = element.GetGeometry().Center()
    #     #     local_point = PointFromGlobalToParamSpace(center, self.lower_bound, self.upper_bound)
    #     #     local_point_kratos = KM.Vector(3)
    #     #     local_point_kratos[0] = local_point[0]
    #     #     local_point_kratos[1] = local_point[1]
    #     #     local_point_kratos[2] = local_point[2]

    #     #     values = element.CalculateOnIntegrationPoints(KM.CAUCHY_STRESS_VECTROR, model_part.ProcessInfo)

    #     return super().OutputSolutionStep()


    def InitializeSolutionStep(self):
        if not self.Initialized:
            model_part = self.model.GetModelPart("NurbsMesh")
            for bc in self.boundary_conditions:
                bc.apply(model_part)
            self.Initialized = True

        return super().InitializeSolutionStep()


