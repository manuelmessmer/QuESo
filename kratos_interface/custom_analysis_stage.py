# Import Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class CustomAnalysisStage(StructuralMechanicsAnalysis):
    """Customized Kratos Analysis Stage.

    Overrides the StructuralMechanicsAnalysis Stage from Kratos.
    """
    def __init__(self, model, general_settings, kratos_settings_filename, integration_points_embedder, boundary_conditions):
        """The constructor."""
        # Read kratos settings
        with open(kratos_settings_filename,'r') as parameter_file:
            analysis_parameters = KM.Parameters(parameter_file.read())

        self.boundary_conditions = boundary_conditions
        self.integration_points_embedder = integration_points_embedder

        #Override the NurbsGeometryModeler input parameters
        for modeler in analysis_parameters["modelers"]:
            if modeler["modeler_name"].GetString() == "NurbsGeometryModeler":
                parameters = modeler["Parameters"]
                parameters.AddEmptyValue("lower_point")
                parameters["lower_point"].SetVector(general_settings["lower_point"])
                parameters.AddEmptyValue("upper_point")
                parameters["upper_point"].SetVector(general_settings["upper_point"])
                parameters.AddEmptyValue("polynomial_order")
                parameters["polynomial_order"].SetVector(general_settings["polynomial_order"])
                parameters.AddEmptyValue("number_of_knot_spans")
                parameters["number_of_knot_spans"].SetVector(general_settings["number_of_knot_spans"])

        self.Initialized = False
        super().__init__(model, analysis_parameters)


    def _ModelersSetupModelPart(self):
        """Override BaseClass to run NURBS modelers."""
        model_part = self.model.GetModelPart('NurbsMesh')
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)

        model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE,3)

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

        integration_points = []
        # Loop over all given integration points
        for point in self.integration_points_embedder:
            integration_points.append([point.GetX(), point.GetY(), point.GetZ(), point.GetWeight()])

        # Create quadrature_point_geometries
        # These are basically just integration points.
        quadrature_point_geometries = KM.GeometriesVector()
        nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2,
           integration_points)

        print("Number quadrature point geometries: ", len(quadrature_point_geometries))
        # Assing element formulation and constitutive law to each integration point.
        for i in range(0, len(quadrature_point_geometries)):
            el = model_part.CreateNewElement('SmallDisplacementElement3D8N', i+1, quadrature_point_geometries[i], volume_properties)

        print("Number of Elements/Integration Points (In BSpline Volume): ", model_part.NumberOfElements())

        # Add Dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

    def InitializeSolutionStep(self):
        if not self.Initialized:
            model_part = self.model.GetModelPart("NurbsMesh")
            for bc in self.boundary_conditions:
                bc.apply(model_part)
            self.Initialized = True

        return super().InitializeSolutionStep()


