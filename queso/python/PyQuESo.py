import os
import shutil
# import QuESo
import QuESoPythonModule as QuESo_App
from .scripts.b_spline_volume import BSplineVolume
from .scripts.helper import *
from .scripts.json_io import JsonIO

class PyQuESo:
    """Main QuESo Python class.

    This class provides an interface to run and interact with the QuESo embedded model.
    """
    def __init__(self, json_filename: str) -> None:
        """Initializes PyQuESo with the provided JSON configuration.

        Args:
            json_filename (str): Path to the JSON configuration file.
        """
        self._settings_holder = JsonIO.read_settings(json_filename)
        settings = self._settings_holder.GetObject()
        write_output_to_file = settings["general_settings"].GetBool("write_output_to_file")
        output_directory_name = settings["general_settings"].GetString("output_directory_name")
        if write_output_to_file:
            folder_path = "./" + output_directory_name + '/'
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
            os.mkdir(folder_path)

    def Run(self) -> None:
        """Runs the QuESo embedded model initialization and generation.
        """
        self._embedded_model = QuESo_App.EmbeddedModel(self._settings_holder) # type: ignore (TODO: add .pyi)
        self._settings_holder = None
        self._embedded_model.CreateAllFromSettings()

    def GetElements(self) -> QuESo_App.ElementVector: # type: ignore (TODO: add .pyi)
        """Returns a list of active elements from the embedded model.

        Returns:
            QuESo_App.ElementVector: List of active elements.
        """
        return self._embedded_model.GetElements()

    def GetConditions(self) -> QuESo_App.ConditionVector: # type: ignore (TODO: add .pyi)
        """Returns a list of conditions from the embedded model.

        Returns:
            list: List of boundary or loading conditions.
        """
        return self._embedded_model.GetConditions()

    def GetSettings(self) -> QuESo_App.Dictionary: # type: ignore (TODO: add .pyi)
        """Returns the QuESo settings used for initialization.

        Returns:
            QuESo_App.Settings: QuESo configuration settings.
        """
        if self._settings_holder and self._settings_holder.GetObject():
            return self._settings_holder.GetObject()
        else:
            return self._embedded_model.GetSettings()

    def GetModelInfo(self) -> QuESo_App.Dictionary: # type: ignore (TODO: add .pyi)
        """Returns information about the current embedded model.

        Returns:
            QuESo_App.ModelInfo: A dictionary containing model metadata.
        """
        return self._embedded_model.GetModelInfo()

    def GetBSplineVolume(self, knot_vector_type: str) -> BSplineVolume:
        """Generates a B-Spline volume from current settings.

        Args:
            knot_vector_type (str): Type of knot vector to use (e.g., "open_knot_vector").

        Returns:
            BSplineVolume: The generated B-Spline volume object.
        """
        return BSplineVolume(self.GetSettings(), knot_vector_type)

    def GetIntegrationPoints(self) -> QuESo_App.IntegrationPointVector: # type: ignore (TODO: add .pyi)
        """Retrieves all integration points used for numerical computations.

        Returns:
            QuESo_App.IntegrationPointVector: Collection of integration points.
        """
        integration_points = QuESo_App.IntegrationPointVector() # type: ignore (TODO: add .pyi)
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
    def RunKratosAnalysis(self, kratos_settings_filename: str="KratosParameters.json") -> None:
        """Runs a Kratos-based analysis using the QuESo model and Kratos parameters.

        Args:
            kratos_settings (str, optional): Path to Kratos settings JSON. Defaults to "KratosParameters.json".

        Raises:
            Exception: If Kratos is not available.
        """
        try:
            import KratosMultiphysics as KM
        except ImportError:
            raise ImportError("RunKratosAnalysis :: Kratos is not available.")

        from kratos_interface.kratos_analysis import Analysis
        self.analysis = Analysis( self.GetSettings(), kratos_settings_filename, self.GetElements(), self.GetConditions() )

