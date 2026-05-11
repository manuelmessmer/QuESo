import os
import json
from copy import deepcopy
from pathlib import Path
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
        self._settings_filename = str(json_filename)
        self._settings_dir = Path(json_filename).resolve().parent
        with open(json_filename, "r") as file:
            self._raw_settings = json.load(file)

        if "models" not in self._raw_settings or "shared_analysis_settings" not in self._raw_settings:
            raise RuntimeError("QuESoSettings must contain 'models' and 'shared_analysis_settings'.")

        self._embedded_models = {}
        self._settings_holder = None

        global_general = self._raw_settings.get("queso", {}).get("general_settings", {})
        output_directory_name = str(global_general.get("output_directory_name", "")).strip()
        if output_directory_name:
            os.makedirs(output_directory_name, exist_ok=True)
            for model in self._raw_settings.get("models", []):
                model_name = str(model.get("name", "")).strip()
                if model_name:
                    os.makedirs(os.path.join(output_directory_name, model_name), exist_ok=True)

    def _ResolvePath(self, filename: str) -> str:
        path = Path(filename)
        if path.is_absolute():
            return str(path)
        return str((self._settings_dir / path).resolve())

    def _BuildMergedModelSettings(self, model_entry: dict) -> dict:
        settings = dict(model_entry.get("queso", {}))
        if "general_settings" not in settings:
            settings["general_settings"] = {}
        general_settings = settings["general_settings"]

        global_general = self._raw_settings.get("queso", {}).get("general_settings", {})
        if "echo_level" in global_general and "echo_level" not in general_settings:
            general_settings["echo_level"] = global_general["echo_level"]

        global_output = global_general.get("output_directory_name")
        model_name = str(model_entry.get("name", "")).strip()
        if global_output and model_name:
            general_settings["output_directory_name"] = str(Path(global_output) / model_name)

        conditions = settings.get("conditions_settings_list", [])
        if isinstance(conditions, list) and len(conditions) == 0:
            settings.pop("conditions_settings_list", None)
            return settings

        for cond in conditions:
            if "condition_id" not in cond:
                raise RuntimeError("Each condition in queso.conditions_settings_list must define 'condition_id'.")
            input_filename = cond.get("input_filename")
            if input_filename is not None:
                cond["input_filename"] = self._ResolvePath(str(input_filename))
        return settings


    def Run(self) -> None:
        """Runs the QuESo embedded model initialization and generation.
        """
        self._embedded_models = {}
        for model_entry in self._raw_settings.get("models", []):
            model_name = str(model_entry.get("name", "")).strip()
            if not model_name:
                raise RuntimeError("Each models entry must contain a non-empty 'name'.")
            settings_dict = self._BuildMergedModelSettings(model_entry)
            settings_holder = JsonIO.read_settings(settings_dict)
            embedded_model = QuESo_App.EmbeddedModel(settings_holder) # type: ignore
            embedded_model.CreateAllFromSettings()
            self._embedded_models[model_name] = embedded_model

            analysis_settings = model_entry.get("analysis_settings", {})
            for element_group in analysis_settings.get("element_groups", []):
                source_name = str(element_group.get("source_name", "")).strip()
                if not source_name:
                    raise RuntimeError(
                        f"Model '{model_name}' has element_group without source_name. "
                        "source_name must be identical to model name."
                    )
                if source_name != model_name:
                    raise RuntimeError(
                        f"Model '{model_name}' has source_name '{source_name}', "
                        "but source_name must match the model name exactly."
                    )

    def _GetModelEntry(self, model_name: str) -> dict:
        model_key = str(model_name).strip()
        for model_entry in self._raw_settings.get("models", []):
            if str(model_entry.get("name", "")).strip() == model_key:
                return model_entry
        raise RuntimeError(f"Model '{model_key}' not available.")

    def GetModelNames(self) -> list[str]:
        return [str(model.get("name", "")).strip() for model in self._raw_settings.get("models", [])]

    def GetElements(self, model_name) -> QuESo_App.ElementVector: # type: ignore (TODO: add .pyi)
        """Returns a list of active elements from the embedded model.

        Returns:
            QuESo_App.ElementVector: List of active elements.
        """
        model_key = str(model_name)
        if model_key not in self._embedded_models:
            raise RuntimeError(f"Model '{model_key}' not available.")
        return self._embedded_models[model_key].GetElements()

    def GetSettings(self, model_name: str | None = None):
        """Returns the QuESo settings used for initialization.

        Returns merged root settings, or one model entry if model_name is provided.
        """
        if model_name is None:
            return deepcopy(self._raw_settings)
        return deepcopy(self._GetModelEntry(model_name))

    def GetConditions(self, model_name: str):
        model_entry = self._GetModelEntry(model_name)
        return deepcopy(model_entry.get("queso", {}).get("conditions_settings_list", []))

    def GetModelInfo(self, model_name: str): # type: ignore (TODO: add .pyi)
        """Returns information about the current embedded model.

        Returns model info for the requested model.
        """
        model_key = str(model_name).strip()
        if model_key not in self._embedded_models:
            raise RuntimeError(f"Model '{model_key}' not available.")
        return self._embedded_models[model_key].GetModelInfo()

    def GetBSplineVolume(self, model_name: str, knot_vector_type: str) -> BSplineVolume:
        """Generates a B-Spline volume from current settings.

        Args:
            knot_vector_type (str): Type of knot vector to use (e.g., "open_knot_vector").

        Returns model-specific B-Spline volume.
        """
        model_entry = self._GetModelEntry(model_name)
        model_queso = deepcopy(model_entry.get("queso", {}))
        settings_holder = JsonIO.read_settings(model_queso)
        return BSplineVolume(settings_holder.GetObject(), knot_vector_type)

    def GetIntegrationPoints(self, model_name: str) -> QuESo_App.IntegrationPointVector: # type: ignore (TODO: add .pyi)
        """Retrieves all integration points used for numerical computations.

        Returns integration points for one model.
        """
        integration_points = QuESo_App.IntegrationPointVector() # type: ignore
        for element in self.GetElements(model_name):
            if element.IsTrimmed():
                for point_trimmed_reduced in element.GetIntegrationPoints():
                    weight = point_trimmed_reduced.Weight()
                    if weight > 0:
                        integration_points.append(point_trimmed_reduced)
            else:
                for point_inside in element.GetIntegrationPoints():
                    integration_points.append(point_inside)
        return integration_points
