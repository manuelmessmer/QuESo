import json
from pathlib import Path
from typing import Any, TypeAlias

import QuESoPythonModule as queso

from .scripts.b_spline_volume import BSplineVolume
from .scripts.json_io import JsonIO
from .scripts.settings_parser import parse_settings_by_component

MainSettingsDict: TypeAlias = dict[str, Any]


class Model:
    """QuESo model containing one or more embedded components.

    The model reads a QuESo settings file, normalizes it to a per-component
    structure, and can then build all configured components via :meth:`run`.
    Accessors that expose generated components require ``run()`` to be called
    first.

    Typical usage:

        >>> model = Model("QuESoSettings.json")
        >>> model.run()
        >>> elements = model.elements("my_component")
    """

    def __init__(self, settings_path: Path) -> None:
        """Initialize the model from a QuESo settings file.

        Args:
            settings_path: Path to the QuESo JSON settings file.

        Raises:
            OSError: If the settings file cannot be read.
            json.JSONDecodeError: If the settings file is not valid JSON.
            ValueError: If the settings file content is invalid.
        """
        self._settings_path = Path(settings_path).resolve()
        self._settings_dir = self._settings_path.parent
        with self._settings_path.open("r") as file:
            self._raw_settings = json.load(file)

        self._component_settings = parse_settings_by_component(
            self._raw_settings, self._settings_dir
        )
        self._embedded_components: dict[str, Any] = {}

    def _verify_component_name(self, component_name: str):
        if component_name not in self._embedded_components:
            raise KeyError(f"Component '{component_name}' not available.")

    def _verify_configured_component_name(self, component_name: str) -> None:
        if component_name not in self._component_settings:
            raise KeyError(f"Component '{component_name}' not available.")

    @property
    def component_names(self) -> list[str]:
        """Return the configured component names in definition order."""
        return list(self._component_settings.keys())

    def run(self) -> None:
        """Build all configured components from the normalized settings.

        Re-running clears previously built components and creates them again
        from the current normalized settings.
        """
        self._embedded_components = {}
        for component_name, component_settings in self._component_settings.items():
            output_directory_name = str(
                component_settings.get("output_directory_name", "")
            ).strip()
            if output_directory_name:
                Path(output_directory_name).mkdir(parents=True, exist_ok=True)

            settings_holder = JsonIO.read_settings(component_settings)
            embedded_component = queso.EmbeddedComponent(settings_holder)  # type: ignore
            embedded_component.CreateAllFromSettings()
            self._embedded_components[component_name] = embedded_component

    def elements(self, component_name: str):
        """Return the generated elements for a component.

        Args:
            component_name: Name of the component to query.

        Returns:
            The QuESo element container for the requested component.

        Raises:
            KeyError: If the component was not configured or has not been built.
        """
        self._verify_component_name(component_name)
        return self._embedded_components[component_name].GetElements()

    def settings(self, component_name: str) -> MainSettingsDict:
        """Return the normalized settings dictionary for a component.

        Args:
            component_name: Name of the component to query.

        Returns:
            The normalized per-component settings dictionary.

        Raises:
            KeyError: If the component is not configured.
        """
        self._verify_configured_component_name(component_name)
        return self._component_settings[component_name]

    def conditions(self, component_name: str):
        """Return the generated conditions for a component.

        Args:
            component_name: Name of the component to query.

        Returns:
            The QuESo condition container for the requested component.

        Raises:
            KeyError: If the component was not configured or has not been built.
        """
        self._verify_component_name(component_name)
        return self._embedded_components[component_name].GetConditions()

    def component_info(self, component_name: str):
        """Return the generated component info dictionary for a component.

        Args:
            component_name: Name of the component to query.

        Returns:
            The QuESo component info dictionary.

        Raises:
            KeyError: If the component was not configured or has not been built.
        """
        self._verify_component_name(component_name)
        return self._embedded_components[component_name].GetComponentInfo()

    def b_spline_volume(
        self, component_name: str, knot_vector_type: str
    ) -> BSplineVolume:
        """Construct a B-spline volume helper for a component.

        Args:
            component_name: Name of the component to query.
            knot_vector_type: Requested knot vector representation.

        Returns:
            A Python ``BSplineVolume`` wrapper for the component.

        Raises:
            KeyError: If the component is not configured.
        """
        self._verify_configured_component_name(component_name)
        settings_holder = JsonIO.read_settings(self._component_settings[component_name])
        return BSplineVolume(settings_holder.GetObject(), knot_vector_type)
