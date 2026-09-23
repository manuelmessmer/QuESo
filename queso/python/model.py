"""High-level QuESo model API."""

import os
import shutil
from typing import Optional

from . import IntegrationPointVector
from ._internal import _EmbeddedModel
from .scripts.b_spline_volume import BSplineVolume
from .scripts.json_io import JsonIO


class Model:
    """Load settings and create an embedded finite-element model."""

    def __init__(self, json_filename: str) -> None:
        """Load model settings from ``json_filename``."""
        self._settings_holder = JsonIO.read_settings(json_filename)
        self._component: Optional[_EmbeddedModel] = None
        self._analysis = None

        settings = self._settings_holder.dictionary
        general_settings = settings["general_settings"]
        write_output_to_file = general_settings.get_bool("write_output_to_file")
        output_directory_name = general_settings.get_string("output_directory_name")
        if write_output_to_file:
            folder_path = os.path.join(".", output_directory_name)
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
            os.mkdir(folder_path)

    def create(self) -> None:
        """Create all model data from the loaded settings."""
        if self._component is not None:
            raise RuntimeError("Model has already been created.")
        self._component = _EmbeddedModel(self._settings_holder)
        self._settings_holder = None
        self._component.create_all_from_settings()

    @property
    def elements(self):
        """Active elements in the created model."""
        return self._created_component.elements

    @property
    def conditions(self):
        """Conditions in the created model."""
        return self._created_component.conditions

    @property
    def settings(self):
        """Settings loaded for this model."""
        if self._settings_holder is not None:
            return self._settings_holder.dictionary
        return self._created_component.settings

    @property
    def model_info(self):
        """Information produced while creating the model."""
        return self._created_component.model_info

    def b_spline_volume(self, knot_vector_type: str) -> BSplineVolume:
        """Construct a B-spline volume from the model settings."""
        return BSplineVolume(self.settings, knot_vector_type)

    @property
    def integration_points(self) -> IntegrationPointVector:
        """All positive-weight integration points in the created model."""
        integration_points = IntegrationPointVector()
        for element in self.elements:
            if element.is_trimmed:
                for point in element.integration_points:
                    if point.weight > 0:
                        integration_points.append(point)
            else:
                for point in element.integration_points:
                    integration_points.append(point)
        return integration_points

    @property
    def analysis(self):
        """Kratos analysis created by :meth:`run_kratos_analysis`."""
        if self._analysis is None:
            raise RuntimeError("Kratos analysis has not been created.")
        return self._analysis

    def run_kratos_analysis(
        self, kratos_settings_filename: str = "KratosParameters.json"
    ) -> None:
        """Run a Kratos analysis using the created QuESo model."""
        try:
            import KratosMultiphysics  # noqa: F401
        except ImportError as exc:
            raise ImportError("KratosMultiphysics is not available.") from exc

        from .kratos_interface.kratos_analysis import Analysis

        self._analysis = Analysis(
            self.settings,
            kratos_settings_filename,
            self.elements,
            self.conditions,
        )

    @property
    def _created_component(self) -> _EmbeddedModel:
        if self._component is None:
            raise RuntimeError(
                "Model.create() must be called before accessing model results."
            )
        return self._component
