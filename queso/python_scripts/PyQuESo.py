import os
import shutil
# import QuESo
import QuESo_PythonApplication as QuESo_App
from queso.python_scripts.b_spline_volume import BSplineVolume
from queso.python_scripts.helper import *
from queso.python_scripts.json_io import JsonIO

try:
    # TODO: Move to analysis
    import KratosMultiphysics as KM
    kratos_available = True
except:
    kratos_available = False

if kratos_available:
    from kratos_interface.kratos_analysis import Analysis
    from kratos_interface.model_part_utilities import ModelPartUtilities

class PyQuESo:
    """Main QuESo Python class.

    This class provides an interface to run and interact with the QuESo embedded model.
    """
    def __init__(self, json_filename: str) -> None:
        """Initializes PyQuESo with the provided JSON configuration.

        Args:
            json_filename (str): Path to the JSON configuration file.
        """
        self.settings = JsonIO.ReadSettings(json_filename)
        write_output_to_file = self.settings["general_settings"].GetBool("write_output_to_file")
        output_directory_name = self.settings["general_settings"].GetString("output_directory_name")
        if write_output_to_file:
            folder_path = "./" + output_directory_name + '/'
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
            os.mkdir(folder_path)

    def Run(self) -> None:
        """Runs the QuESo embedded model initialization and generation.
        """
        self.embedded_model = QuESo_App.EmbeddedModel(self.settings)
        self.embedded_model.CreateAllFromSettings()

    def GetElements(self) -> QuESo_App.ElementVector:
        """Returns a list of active elements from the embedded model.

        Returns:
            QuESo_App.ElementVector: List of active elements.
        """
        return self.embedded_model.GetElements()

    def GetConditions(self) -> QuESo_App.ConditionVector:
        """Returns a list of conditions from the embedded model.

        Returns:
            list: List of boundary or loading conditions.
        """
        return self.embedded_model.GetConditions()

    def GetSettings(self) -> QuESo_App.Settings:
        """Returns the QuESo settings used for initialization.

        Returns:
            QuESo_App.Settings: QuESo configuration settings.
        """
        return self.embedded_model.GetSettings()

    def GetModelInfo(self) -> QuESo_App.ModelInfo:
        """Returns information about the current embedded model.

        Returns:
            QuESo_App.ModelInfo: A dictionary containing model metadata.
        """
        return self.embedded_model.GetModelInfo()

    def GetBSplineVolume(self, knot_vector_type: str) -> BSplineVolume:
        """Generates a B-Spline volume from current settings.

        Args:
            knot_vector_type (str): Type of knot vector to use (e.g., "open_knot_vector").

        Returns:
            BSplineVolume: The generated B-Spline volume object.
        """
        return BSplineVolume(self.settings, knot_vector_type)

    def GetIntegrationPoints(self) -> QuESo_App.IntegrationPointVector:
        """Retrieves all integration points used for numerical computations.

        Returns:
            QuESo_App.IntegrationPointVector: Collection of integration points.
        """
        integration_points = QuESo_App.IntegrationPointVector()
        # Gather all poitnts (TODO: make this in C++)
        for element in self.GetElements():
            if element.IsTrimmed():
                for point_trimmed_reduced in element.GetIntegrationPoints():
                    weight = point_trimmed_reduced.Weight()
                    if( weight > 0):
                        integration_points.append(point_trimmed_reduced)
            else:
                for point_inside in element.GetIntegrationPoints():
                    integration_points.append(point_inside)
        return integration_points

    def GetAnalysis(self):
        """Returns the Kratos analysis object.

        Returns:
            Analysis: Kratos analysis object, if it exists.

        Raises:
            Exception: If Kratos analysis is not available.
        """
        if not self.analysis:
            raise Exception(f'Kratos analysis has not been created.')
        return self.analysis

    #########################################
    #### Kratos related member functions ####
    #########################################
    def UpdateKratosNurbsVolumeModelPart(self, kratos_model_part) -> None:
        """Updates a Kratos model part with QuESo geometry and conditions.

        Args:
            kratos_model_part (KM.ModelPart): Kratos model part to update.

        Raises:
            Exception: If Kratos is not available.
        """
        if kratos_available:
            ModelPartUtilities.RemoveAllElements(kratos_model_part)
            ModelPartUtilities.RemoveAllConditions(kratos_model_part)
            ModelPartUtilities.AddElementsToModelPart(kratos_model_part, self.GetElements())
            ModelPartUtilities.AddConditionsToModelPart(kratos_model_part, self.conditions, self.GetBoundsXYZ(), self.GetBoundsUVW())
        else:
            raise Exception("UpdateKratosNurbsVolumeModelPart :: Kratos is not available.")


    def RunKratosAnalysis(self, kratos_settings_filename: str="KratosParameters.json") -> None:
        """Runs a Kratos-based analysis using the QuESo model and Kratos parameters.

        Args:
            kratos_settings (str, optional): Path to Kratos settings JSON. Defaults to "KratosParameters.json".

        Raises:
            Exception: If Kratos is not available.
        """
        if kratos_available:
            self.analysis = Analysis( self.settings, kratos_settings_filename, self.GetElements(), self.GetConditions() )
        else:
            raise Exception("RunKratosAnalysis :: Kratos is not available.")
