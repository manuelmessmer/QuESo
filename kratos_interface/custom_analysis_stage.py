# Import Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from kratos_interface.model_part_utilities import ModelPartUtilities

class CustomAnalysisStage(StructuralMechanicsAnalysis):
    """Customized Kratos Analysis Stage.

    Overrides the StructuralMechanicsAnalysis Stage from Kratos.
    """
    def __init__(self, model, tibra_parameters, kratos_settings_filename, elements, boundary_conditions):
        """The constructor."""
        # Read kratos settings
        with open(kratos_settings_filename,'r') as parameter_file:
            analysis_parameters = KM.Parameters(parameter_file.read())

        self.boundary_conditions = boundary_conditions
        self.elements = elements
        self.tibra_parameters = tibra_parameters
        #Override the NurbsGeometryModeler input parameters
        for modeler in analysis_parameters["modelers"].values():
            if modeler["modeler_name"].GetString() == "NurbsGeometryModeler":
                parameters = modeler["Parameters"]
                parameters.AddEmptyValue("lower_point")
                parameters["lower_point"].SetVector(self.tibra_parameters.LowerBound())
                parameters.AddEmptyValue("upper_point")
                parameters["upper_point"].SetVector(self.tibra_parameters.UpperBound())
                parameters.AddEmptyValue("polynomial_order")
                parameters["polynomial_order"].SetVector(self.tibra_parameters.Order())
                parameters.AddEmptyValue("number_of_knot_spans")
                parameters["number_of_knot_spans"].SetVector(self.tibra_parameters.NumberOfElements())

        self.Initialized = False
        super().__init__(model, analysis_parameters)


    def _ModelersSetupModelPart(self):
        """Override BaseClass to run NURBS modelers."""
        model_part = self.model.GetModelPart('NurbsMesh')
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)

        model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        # Convert the geometry model or import analysis suitable models.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KM.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart started.")
            modeler.SetupModelPart()
            if self.echo_level > 1:
                KM.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart finished.")

        return super()._ModelersSetupModelPart()


    def ModifyInitialGeometry(self):
        """Override BaseClass to pass integration points to Kratos."""
        model_part = self.model.GetModelPart('NurbsMesh')
        ModelPartUtilities.RemoveAllElements(model_part)
        ModelPartUtilities.RemoveAllConditions(model_part)
        ModelPartUtilities.AddElementsToModelPart(model_part, self.elements)
        ModelPartUtilities.AddConditionsToModelPart(model_part, self.boundary_conditions, self.tibra_parameters.LowerBound(), self.tibra_parameters.UpperBound())

        # Add Dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)


