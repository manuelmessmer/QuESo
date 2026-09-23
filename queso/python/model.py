"""High-level multi-component QuESo model API."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import IntegrationPointVector  # type: ignore
from ._internal import _EmbeddedComponent  # type: ignore
from .scripts.b_spline_volume import BSplineVolume
from .scripts.json_io import JsonIO
from .scripts.settings_parser import parse_settings_by_component


class Model:
    """Load, validate, and create one or more embedded components."""

    def __init__(self, settings_path: str | Path) -> None:
        """Load public QuESo settings from a JSON file.

        Args:
            settings_path: Path to the public multi-component settings file.

        Raises:
            OSError: If the settings file cannot be read.
            json.JSONDecodeError: If the file does not contain valid JSON.
            ValueError: If the settings schema is invalid.
        """
        self._settings_path = Path(settings_path).expanduser().resolve()
        with self._settings_path.open("r", encoding="utf-8") as settings_file:
            raw_settings = json.load(settings_file)
        self._component_settings = parse_settings_by_component(
            raw_settings, self._settings_path
        )
        self._components: dict[str, _EmbeddedComponent] | None = None

    @property
    def component_names(self) -> tuple[str, ...]:
        """Configured component names in definition order."""
        return tuple(self._component_settings)

    def create(self) -> None:
        """Create every configured embedded component exactly once."""
        if self._components is not None:
            raise RuntimeError("Model.create() may only be called once.")

        components: dict[str, _EmbeddedComponent] = {}
        for component_name, settings in self._component_settings.items():
            if settings["write_output_to_file"]:
                Path(settings["output_directory_name"]).mkdir(
                    parents=True, exist_ok=True
                )
            settings_holder = JsonIO.read_settings(settings)
            component = _EmbeddedComponent(settings_holder)
            component.create_all_from_settings()
            components[component_name] = component
        self._components = components

    def settings(self, component_name: str) -> dict[str, Any]:
        """Return normalized settings for a configured component."""
        self._verify_configured_component(component_name)
        return self._component_settings[component_name]

    def elements(self, component_name: str):
        """Return the created component's lazy element range."""
        return self._component(component_name).elements

    def conditions(self, component_name: str):
        """Return the created component's conditions."""
        return self._component(component_name).conditions

    def component_info(self, component_name: str):
        """Return information generated while creating a component."""
        return self._component(component_name).component_info

    def integration_points(self, component_name: str) -> IntegrationPointVector:
        """Return all positive-weight volume integration points for a component."""
        integration_points = IntegrationPointVector()
        for element in self.elements(component_name):
            for point in element.integration_points:
                if point.weight > 0.0:
                    integration_points.append(point)
        return integration_points

    def b_spline_volume(
        self, component_name: str, knot_vector_type: str
    ) -> BSplineVolume:
        """Construct a B-spline helper from one component's grid settings."""
        settings_holder = JsonIO.read_settings(self.settings(component_name))
        return BSplineVolume(settings_holder.dictionary, knot_vector_type)

    def _verify_configured_component(self, component_name: str) -> None:
        if component_name not in self._component_settings:
            raise KeyError(f"Unknown component '{component_name}'.")

    def _component(self, component_name: str) -> _EmbeddedComponent:
        self._verify_configured_component(component_name)
        if self._components is None:
            raise RuntimeError(
                "Model.create() must be called before accessing generated results."
            )
        return self._components[component_name]
