"""High-level Kratos analysis that owns a multi-component pyqueso model."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import KratosMultiphysics as KM  # type: ignore
import KratosMultiphysics.IgaApplication as IGA  # type: ignore
import KratosMultiphysics.LinearSolversApplication  # type: ignore
import pyqueso  # type: ignore
from KratosMultiphysics.StructuralMechanicsApplication import (  # type: ignore
    structural_mechanics_analysis,
)

from .analysis_settings_parser import parse_analysis_settings
from .model_part_utilities import (
    add_conditions,
    add_elements,
    embedded_model_part_name,
    ensure_model_part,
    read_triangle_mesh_to_model_part,
    remove_all_conditions,
    remove_all_elements,
)

PathLike = str | Path


class Analysis:
    """Run a Kratos analysis from QuESo, analysis, and Kratos JSON inputs."""

    def __init__(
        self,
        queso_settings_path: PathLike,
        analysis_settings_path: PathLike,
        kratos_parameters_path: PathLike,
    ) -> None:
        """Store required input paths without constructing either model."""
        self._queso_settings_path = Path(queso_settings_path).expanduser().resolve()
        self._analysis_settings_path = (
            Path(analysis_settings_path).expanduser().resolve()
        )
        self._kratos_parameters_path = (
            Path(kratos_parameters_path).expanduser().resolve()
        )
        self._queso_model: pyqueso.Model | None = None
        self._kratos_model: Any | None = None
        self._has_run = False

    @property
    def queso_model(self) -> pyqueso.Model:
        """Return the owned QuESo model after analysis construction."""
        if self._queso_model is None:
            raise RuntimeError("Analysis.run() has not created the QuESo model, yet.")
        return self._queso_model

    @property
    def kratos_model(self) -> Any:
        """Return the owned native Kratos model after analysis construction."""
        if self._kratos_model is None:
            raise RuntimeError("Analysis.run() has not created the Kratos model, yet.")
        return self._kratos_model

    @property
    def kratos_analysis_stage(self) -> Any:
        """Return the owned native Kratos analysis stage after analysis construction."""
        if self._stage is None:
            raise RuntimeError(
                "Analysis.run() has not created the Kratos analysis stage, yet."
            )
        return self._stage

    def run(self) -> None:
        """Create both models and delegate execution to Kratos' analysis stage."""
        if self._has_run:
            raise RuntimeError("Analysis.run() may only be called once.")
        self._has_run = True

        queso_model = pyqueso.Model(self._queso_settings_path)
        queso_model.create()
        self._queso_model = queso_model

        with self._analysis_settings_path.open("r", encoding="utf-8") as file:
            raw_analysis_settings = json.load(file)
        analysis_settings = parse_analysis_settings(raw_analysis_settings, queso_model)

        with self._kratos_parameters_path.open("r", encoding="utf-8") as file:
            raw_kratos_parameters = json.load(file)
        _resolve_kratos_paths(
            raw_kratos_parameters, self._kratos_parameters_path.parent
        )
        _inject_modelers(raw_kratos_parameters, analysis_settings)
        _inject_embedded_geometry_output_processes(
            raw_kratos_parameters, analysis_settings
        )
        raw_kratos_parameters["solver_settings"]["model_part_name"] = analysis_settings[
            "root_model_part_name"
        ]

        kratos_model = KM.Model()
        self._kratos_model = kratos_model
        self._stage = _AnalysisStage(
            kratos_model,
            KM.Parameters(json.dumps(raw_kratos_parameters)),
            queso_model,
            analysis_settings,
        )
        self._stage.Run()


class _AnalysisStage(structural_mechanics_analysis.StructuralMechanicsAnalysis):
    """Internal Kratos stage that assembles QuESo-generated entities."""

    def __init__(
        self,
        model: Any,
        project_parameters: Any,
        queso_model: pyqueso.Model,
        analysis_settings: dict[str, Any],
    ) -> None:
        self._queso_model = queso_model
        self._analysis_settings = analysis_settings
        self._condition_property_ids: dict[int, int] = {}
        self._setup_model_parts(model)
        super().__init__(model, project_parameters)

    def _setup_model_parts(self, model: Any) -> None:
        root_name = self._analysis_settings["root_model_part_name"]
        root = ensure_model_part(model, root_name)
        root.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        root.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
        root.AddNodalSolutionStepVariable(KM.REACTION)
        if self._requires_lagrange_dofs():
            root.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)
            root.AddNodalSolutionStepVariable(IGA.VECTOR_LAGRANGE_MULTIPLIER_REACTION)

        for geometry in self._analysis_settings["geometries"].values():
            model_part = ensure_model_part(model, geometry["model_part_name"])
            for component in geometry["queso_components"]:
                name = component["component_name"]
                if not model_part.HasSubModelPart(name):
                    model_part.CreateSubModelPart(name)

        for component_name in self._queso_model.component_names:
            embedded_model_part = ensure_model_part(
                model, embedded_model_part_name(component_name)
            )
            embedded_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            embedded_model_part.AddNodalSolutionStepVariable(KM.REACTION)
            embedded_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            read_triangle_mesh_to_model_part(
                embedded_model_part,
                self._queso_model.settings(component_name)["input_filename"],
            )

    def _requires_lagrange_dofs(self) -> bool:
        for component_name in self._queso_model.component_names:
            for condition in self._queso_model.settings(component_name).get(
                "conditions_settings_list", []
            ):
                if condition["condition_type"] == "LagrangeSupportCondition":
                    return True
        return False

    def ModifyInitialGeometry(self) -> None:
        """Replace modeler elements with QuESo quadrature entities and conditions."""
        processed_model_parts = set()
        for geometry in self._analysis_settings["geometries"].values():
            model_part_name = geometry["model_part_name"]
            if model_part_name not in processed_model_parts:
                model_part = self.model.GetModelPart(model_part_name)
                remove_all_elements(model_part)
                remove_all_conditions(model_part)
                processed_model_parts.add(model_part_name)

        for geometry_name, geometry in self._analysis_settings["geometries"].items():
            model_part = self.model.GetModelPart(geometry["model_part_name"])
            grid = geometry["background_grid_settings"]
            bounds_xyz = (grid["lower_bound_xyz"], grid["upper_bound_xyz"])
            bounds_uvw = (grid["lower_bound_uvw"], grid["upper_bound_uvw"])
            for component in geometry["queso_components"]:
                component_name = component["component_name"]
                component_model_part = model_part.GetSubModelPart(component_name)
                element_settings = component["element_settings"]
                add_elements(
                    model_part,
                    component_model_part,
                    self._queso_model.elements(component_name),
                    element_settings["property_id"],
                    geometry_name,
                    element_settings["element_type"],
                )
                add_conditions(
                    model_part,
                    component_model_part,
                    self._queso_model.conditions(component_name),
                    bounds_xyz,
                    bounds_uvw,
                    geometry_name,
                    self._condition_property_id,
                    self._coupling_target,
                )

            KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)
            if self._requires_lagrange_dofs():
                KM.VariableUtils().AddDof(
                    KM.VECTOR_LAGRANGE_MULTIPLIER_X,
                    IGA.VECTOR_LAGRANGE_MULTIPLIER_REACTION_X,
                    model_part,
                )
                KM.VariableUtils().AddDof(
                    KM.VECTOR_LAGRANGE_MULTIPLIER_Y,
                    IGA.VECTOR_LAGRANGE_MULTIPLIER_REACTION_Y,
                    model_part,
                )
                KM.VariableUtils().AddDof(
                    KM.VECTOR_LAGRANGE_MULTIPLIER_Z,
                    IGA.VECTOR_LAGRANGE_MULTIPLIER_REACTION_Z,
                    model_part,
                )

    def _condition_property_id(self, condition_id: int) -> int:
        if condition_id in self._condition_property_ids:
            return self._condition_property_ids[condition_id]
        root = self.model.GetModelPart(
            self._analysis_settings["root_model_part_name"]
        ).GetRootModelPart()
        property_id = 1
        while root.HasProperties(property_id):
            property_id += 1
        root.CreateNewProperties(property_id)
        self._condition_property_ids[condition_id] = property_id
        return property_id

    def _coupling_target(self, component_name: str) -> dict[str, Any]:
        geometry_name = self._analysis_settings["component_geometry"][component_name]
        geometry = self._analysis_settings["geometries"][geometry_name]
        grid = geometry["background_grid_settings"]
        return {
            "model_part": self.model.GetModelPart(geometry["model_part_name"]),
            "geometry_name": geometry_name,
            "bounds_xyz": (grid["lower_bound_xyz"], grid["upper_bound_xyz"]),
            "bounds_uvw": (grid["lower_bound_uvw"], grid["upper_bound_uvw"]),
        }


def _inject_modelers(
    kratos_parameters: dict[str, Any], analysis_settings: dict[str, Any]
) -> None:
    modelers = kratos_parameters.setdefault("modelers", [])
    if not isinstance(modelers, list):
        raise ValueError("KratosParameters.modelers must be an array.")
    if any(
        modeler.get("modeler_name") == "NurbsGeometryModeler" for modeler in modelers
    ):
        raise ValueError("NurbsGeometryModeler entries are injected automatically.")

    for geometry_name, geometry in analysis_settings["geometries"].items():
        grid = geometry["background_grid_settings"]
        modelers.append(
            {
                "modeler_name": "NurbsGeometryModeler",
                "Parameters": {
                    "model_part_name": geometry["model_part_name"],
                    "geometry_name": geometry_name,
                    "lower_point_xyz": grid["lower_bound_xyz"],
                    "upper_point_xyz": grid["upper_bound_xyz"],
                    "lower_point_uvw": grid["lower_bound_uvw"],
                    "upper_point_uvw": grid["upper_bound_uvw"],
                    "polynomial_order": grid["polynomial_order"],
                    "number_of_knot_spans": grid["number_of_elements"],
                },
            }
        )


def _inject_embedded_geometry_output_processes(
    kratos_parameters: dict[str, Any], analysis_settings: dict[str, Any]
) -> None:
    output_processes = kratos_parameters.setdefault("output_processes", {})
    if not isinstance(output_processes, dict):
        raise ValueError("KratosParameters.output_processes must be an object.")
    for processes in output_processes.values():
        if not isinstance(processes, list):
            raise ValueError("Each Kratos output process group must be an array.")
        if any(
            process.get("process_name") == "VtkEmbeddedGeometryOutputProcess"
            for process in processes
        ):
            raise ValueError(
                "VtkEmbeddedGeometryOutputProcess entries are injected automatically."
            )

    vtk_processes = output_processes.setdefault("vtk_output", [])
    for geometry_name, geometry in analysis_settings["geometries"].items():
        for component in geometry["queso_components"]:
            component_name = component["component_name"]
            embedded_name = embedded_model_part_name(component_name)
            vtk_processes.append(
                {
                    "python_module": "vtk_embedded_geometry_output_process",
                    "kratos_module": "KratosMultiphysics.IgaApplication",
                    "process_name": "VtkEmbeddedGeometryOutputProcess",
                    "help": (
                        "Write embedded-geometry VTK output for one QuESo component."
                    ),
                    "Parameters": {
                        "mapping_parameters": {
                            "main_model_part_name": geometry["model_part_name"],
                            "nurbs_volume_name": geometry_name,
                            "embedded_model_part_name": embedded_name,
                        },
                        "vtk_parameters": {
                            "model_part_name": embedded_name,
                            "output_control_type": "step",
                            "output_interval": 1,
                            "file_format": "ascii",
                            "output_precision": 7,
                            "output_sub_model_parts": False,
                            "output_path": f"kratos_output/{component_name}",
                            "save_output_files_in_folder": True,
                            "nodal_solution_step_data_variables": ["DISPLACEMENT"],
                            "nodal_data_value_variables": [
                                "CAUCHY_STRESS_VECTOR",
                                "VON_MISES_STRESS",
                            ],
                            "nodal_flags": [],
                            "element_data_value_variables": [],
                            "element_flags": [],
                            "condition_data_value_variables": [],
                            "condition_flags": [],
                            "gauss_point_variables_extrapolated_to_nodes": [],
                        },
                    },
                }
            )


def _resolve_kratos_paths(parameters: dict[str, Any], directory: Path) -> None:
    solver_settings = parameters.get("solver_settings", {})
    material_settings = solver_settings.get("material_import_settings", {})
    material_filename = material_settings.get("materials_filename")
    if isinstance(material_filename, str) and material_filename:
        path = Path(material_filename).expanduser()
        if not path.is_absolute():
            material_settings["materials_filename"] = str((directory / path).resolve())
