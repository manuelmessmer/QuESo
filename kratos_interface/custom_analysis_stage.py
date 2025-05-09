# Import QuESo
import QuESo_PythonApplication as QuESo

# Import Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IgaApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from kratos_interface.model_part_utilities import ModelPartUtilities
import KratosMultiphysics.LinearSolversApplication

class CustomAnalysisStage(StructuralMechanicsAnalysis):
    """
    Customized Kratos Analysis Stage.

    This class overrides the `StructuralMechanicsAnalysis` stage from Kratos to include custom
    behavior specific to the needs of a QuESo-based analysis.

    The Kratos analysis uses a BSpline box as background grid.
    """
    def __init__(self,
            model: KM.Model,
            queso_settings: QuESo.Settings, # type: ignore (TODO: add .pyi)
            kratos_settings_filename: str,
            elements: QuESo.ElementVector, # type: ignore (TODO: add .pyi)
            boundary_conditions: QuESo.ConditionVector # type: ignore (TODO: add .pyi)
        ) -> None:
        """
        Constructor for CustomAnalysisStage.

        Initializes the class, reads the Kratos settings, configures the boundary conditions and
        elements, and sets up the model part.

        Args:
            model (KM.Model): The Kratos model to be used for the analysis.
            queso_settings (QuESo.Settings): The settings for the analysis, including conditions and grid settings.
            kratos_settings_filename (str): Path to the Kratos settings JSON file.
            elements (QuESo.ElementVector): List of elements to be added to the model part.
            boundary_conditions (QuESo.ConditionVector): List of boundary conditions to be applied to the model.
        """

        # Read kratos settings
        with open(kratos_settings_filename,'r') as parameter_file:
            analysis_parameters = KM.Parameters(parameter_file.read())

        # Initialize members
        self.boundary_conditions = boundary_conditions
        self.elements = elements
        self.queso_settings = queso_settings
        self.lagrange_dofs_required = False

        # Set up model part
        for condition_param in self.queso_settings.GetList("conditions_settings_list"):
            if( condition_param.GetString("condition_type") == "LagrangeSupportCondition" ):
                self.lagrange_dofs_required = True
        nurbs_model_part = model.CreateModelPart("NurbsMesh")
        nurbs_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        nurbs_model_part.AddNodalSolutionStepVariable(KM.REACTION)
        if self.lagrange_dofs_required:
            nurbs_model_part.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)
            nurbs_model_part.AddNodalSolutionStepVariable(IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION) # type: ignore

        # Override the NurbsGeometryModeler input parameters
        grid_settings = self.queso_settings["background_grid_settings"]
        for modeler in analysis_parameters["modelers"].values():
            if modeler["modeler_name"].GetString() == "NurbsGeometryModeler":
                parameters = modeler["Parameters"]
                parameters.AddEmptyValue("lower_point_xyz")
                parameters["lower_point_xyz"].SetVector(grid_settings.GetDoubleVector("lower_bound_xyz"))
                parameters.AddEmptyValue("upper_point_xyz")
                parameters["upper_point_xyz"].SetVector(grid_settings.GetDoubleVector("upper_bound_xyz"))

                parameters.AddEmptyValue("lower_point_uvw")
                parameters["lower_point_uvw"].SetVector(grid_settings.GetDoubleVector("lower_bound_uvw"))
                parameters.AddEmptyValue("upper_point_uvw")
                parameters["upper_point_uvw"].SetVector(grid_settings.GetDoubleVector("upper_bound_uvw"))

                parameters.AddEmptyValue("polynomial_order")
                parameters["polynomial_order"].SetVector(grid_settings.GetIntVector("polynomial_order"))
                parameters.AddEmptyValue("number_of_knot_spans")
                parameters["number_of_knot_spans"].SetVector(grid_settings.GetIntVector("number_of_elements"))

        self.Initialized = False
        super().__init__(model, analysis_parameters)


    def _ModelersSetupModelPart(self) -> None:
        """
        Override BaseClass to run NURBS modelers.

        This method creates a new model part for the embedded mesh, adds necessary variables,
        reads the mesh, and sets up modelers.
        """
        embedded_model_part = self.model.CreateModelPart('EmbeddedModelPart')
        embedded_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        embedded_model_part.AddNodalSolutionStepVariable(KM.REACTION)
        embedded_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        filename = self.queso_settings["general_settings"].GetString("input_filename")
        self.triangle_mesh = QuESo.TriangleMesh() # type: ignore (TODO: add .pyi)
        QuESo.IO.ReadMeshFromSTL(self.triangle_mesh, filename) # type: ignore (TODO: add .pyi)
        ModelPartUtilities.read_model_part_from_triangle_mesh(embedded_model_part, self.triangle_mesh)

        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KM.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart started.")
            modeler.SetupModelPart()
            if self.echo_level > 1:
                KM.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart finished.")
        super()._ModelersSetupModelPart()


    def ModifyInitialGeometry(self) -> None:
        """
        Override BaseClass to pass integration points to Kratos.

        This method modifies the initial geometry by clearing existing elements and conditions
        from the model part, adding new ones (with QuESo's integration points), and setting up DOFs.
        """
        model_part = self.model.GetModelPart('NurbsMesh')
        ModelPartUtilities.remove_all_elements(model_part)
        ModelPartUtilities.remove_all_conditions(model_part)
        ModelPartUtilities.add_elements_to_model_part(model_part, self.elements)
        grid_settings = self.queso_settings["background_grid_settings"]
        bounds_xyz = (grid_settings.GetDoubleVector("lower_bound_xyz"),
                      grid_settings.GetDoubleVector("upper_bound_xyz"))
        bounds_uvw = (grid_settings.GetDoubleVector("lower_bound_uvw"),
                      grid_settings.GetDoubleVector("upper_bound_uvw"))
        ModelPartUtilities.add_conditions_to_model_part(model_part, self.boundary_conditions, bounds_xyz, bounds_uvw)

        # Add Dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

        if self.lagrange_dofs_required:
            KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_X, IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION_X, model_part) # type: ignore
            KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Y, IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION_Y, model_part) # type: ignore
            KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Z, IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION_Z, model_part) # type: ignore


