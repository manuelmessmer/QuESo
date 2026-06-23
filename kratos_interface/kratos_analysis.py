import json
import math
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import QuESoPythonModule as QuESo
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)

from QuESoPythonModule.model import Model
from QuESoPythonModule.kratos_interface.analysis_settings_parser import (
    parse_geometry_analysis_settings,
)
from QuESoPythonModule.kratos_interface.model_part_utilities import (
    remove_all_elements,
    remove_all_conditions,
    remove_all_nodes,
    add_conditions,
    add_elements,
    read_model_part_from_triangle_mesh,
)


class KratosAnalysis(StructuralMechanicsAnalysis):
    def __init__(
        self,
        model: KM.Model,
        queso_settings_filename: str = "QuESoSettings.json",
        analysis_settings_filename: str = "AnalysisSettings.json",
        kratos_settings_filename: str = "KratosParameters.json",
    ) -> None:
        self.model = model
        self.queso_settings_filename = str(Path(queso_settings_filename).resolve())
        self._settings_dir = Path(self.queso_settings_filename).parent
        self.analysis_settings_filename = self._resolve_path(analysis_settings_filename)
        self.kratos_settings_filename = self._resolve_path(kratos_settings_filename)

        self.queso_model = Model(self.queso_settings_filename)
        self.queso_model.run()

        with open(self.analysis_settings_filename, "r") as settings_file:
            raw_analysis_settings = json.load(settings_file)
        with open(self.kratos_settings_filename, "r") as parameter_file:
            analysis_parameters = KM.Parameters(parameter_file.read())

        self.analysis_settings = parse_geometry_analysis_settings(
            analysis_settings=raw_analysis_settings,
            queso_model=self.queso_model,
        )
        self._max_step = self._compute_max_step(analysis_parameters)
        self._ramped_condition_loads: dict[int, tuple[float, float, float]] = {}
        analysis_parameters["solver_settings"]["model_part_name"].SetString(
            self.analysis_settings["root_model_part_name"]
        )

        self._setup_model_parts()
        self._set_nurbs_volume_modeler_parameters(analysis_parameters)

        super().__init__(model, analysis_parameters)

    def _resolve_path(self, filename: str) -> str:
        path = Path(filename)
        if path.is_absolute():
            return str(path)

        if path.exists():
            return str(path.resolve())

        return str((self._settings_dir / path).resolve())

    def _compute_max_step(self, analysis_parameters: KM.Parameters) -> int:
        start_time = analysis_parameters["problem_data"]["start_time"].GetDouble()
        end_time = analysis_parameters["problem_data"]["end_time"].GetDouble()
        time_step = analysis_parameters["solver_settings"]["time_stepping"][
            "time_step"
        ].GetDouble()

        if time_step <= 0.0:
            raise ValueError("Kratos time_stepping.time_step must be > 0.")

        duration = end_time - start_time
        if duration < 0.0:
            raise ValueError("Kratos end_time must be >= start_time.")

        if duration == 0.0:
            return 1

        return max(1, int(round(duration / time_step)))

    def _create_condition_property_ids(
        self,
        model_part: KM.ModelPart,
        component_name: str,
        active_condition_ids: list[int],
    ) -> dict[int, int]:
        property_id_by_condition_id: dict[int, int] = {}
        next_property_id = model_part.GetRootModelPart().NumberOfProperties() + 1
        available_condition_ids = {
            condition.GetSettings().GetInt("condition_id")
            for condition in self.queso_model.conditions(component_name)
        }

        for condition_id in active_condition_ids:
            normalized_condition_id = int(condition_id)
            if normalized_condition_id not in available_condition_ids:
                continue

            while model_part.GetRootModelPart().HasProperties(next_property_id):
                next_property_id += 1
            model_part.CreateNewProperties(next_property_id)
            property_id_by_condition_id[normalized_condition_id] = next_property_id
            next_property_id += 1

        return property_id_by_condition_id

    def _setup_model_parts(self) -> None:
        root_model_part_name = self.analysis_settings["root_model_part_name"]
        if not self.model.HasModelPart(root_model_part_name):
            self.model.CreateModelPart(root_model_part_name)

        for geometry_settings in self.analysis_settings["geometries"].values():
            model_part_name = str(geometry_settings["model_part_name"])
            model_part = (
                self.model.GetModelPart(model_part_name)
                if self.model.HasModelPart(model_part_name)
                else self.model.CreateModelPart(model_part_name)
            )
            model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            model_part.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
            model_part.AddNodalSolutionStepVariable(KM.REACTION)
            if self.analysis_settings["lagrange_dofs_required"]:
                model_part.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)
                model_part.AddNodalSolutionStepVariable(IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION)  # type: ignore

            for queso_component in geometry_settings.get("queso_components", []):
                group_name = str(queso_component.get("component_name", "")).strip()
                if group_name and not model_part.HasSubModelPart(group_name):
                    model_part.CreateSubModelPart(group_name)

        for component_name in self.queso_model.component_names:
            component_settings = self.queso_model.settings(component_name)
            filename = str(component_settings.get("input_filename", "")).strip()
            if not filename:
                continue

            model_part_name = f"embedded_{component_name}"
            if self.model.HasModelPart(model_part_name):
                continue

            embedded_model_part = self.model.CreateModelPart(model_part_name)
            embedded_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            embedded_model_part.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
            embedded_model_part.AddNodalSolutionStepVariable(KM.REACTION)
            embedded_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    def _set_nurbs_volume_modeler_parameters(
        self, analysis_parameters: KM.Parameters
    ) -> None:
        for modeler in analysis_parameters["modelers"].values():
            if modeler["modeler_name"].GetString() != "NurbsGeometryModeler":
                continue
            parameters = modeler["Parameters"]
            model_part_name = parameters["model_part_name"].GetString()
            requested_geometry_name = (
                parameters["geometry_name"].GetString()
                if parameters.Has("geometry_name")
                else "NurbsVolume"
            )

            geometry_settings = self.analysis_settings["geometries"].get(
                requested_geometry_name
            )
            expected_model_part_name = (
                geometry_settings["model_part_name"] if geometry_settings is not None else None
            )
            matches_model_part = expected_model_part_name == model_part_name or (
                expected_model_part_name is not None
                and str(expected_model_part_name).split(".")[-1] == str(model_part_name)
            )
            if geometry_settings is None or not matches_model_part:
                continue
            grid_settings = geometry_settings["background_grid_settings"]

            for parameter_name, grid_name in (
                ("lower_point_xyz", "lower_bound_xyz"),
                ("upper_point_xyz", "upper_bound_xyz"),
                ("lower_point_uvw", "lower_bound_uvw"),
                ("upper_point_uvw", "upper_bound_uvw"),
            ):
                if parameters.Has(parameter_name):
                    parameters.RemoveValue(parameter_name)
                parameters.AddEmptyValue(parameter_name)
                parameters[parameter_name].SetVector(grid_settings[grid_name])

            parameters["model_part_name"].SetString(geometry_settings["model_part_name"])

            if parameters.Has("polynomial_order"):
                parameters.RemoveValue("polynomial_order")
            parameters.AddEmptyValue("polynomial_order")
            parameters["polynomial_order"].SetVector(
                grid_settings["polynomial_order"]
            )

            if parameters.Has("number_of_knot_spans"):
                parameters.RemoveValue("number_of_knot_spans")
            parameters.AddEmptyValue("number_of_knot_spans")
            parameters["number_of_knot_spans"].SetVector(
                grid_settings["number_of_elements"]
            )

            if parameters.Has("geometry_name"):
                parameters["geometry_name"].SetString(requested_geometry_name)
            else:
                parameters.AddEmptyValue("geometry_name")
                parameters["geometry_name"].SetString(requested_geometry_name)

    def ModifyInitialGeometry(self) -> None:
        self._ramped_condition_loads = {}

        for component_name in self.queso_model.component_names:
            component_settings = self.queso_model.settings(component_name)
            filename = str(component_settings.get("input_filename", "")).strip()
            if not filename:
                continue

            embedded_model_part = self.model.GetModelPart(f"embedded_{component_name}")
            remove_all_elements(embedded_model_part)
            remove_all_conditions(embedded_model_part)
            remove_all_nodes(embedded_model_part)

            triangle_mesh = QuESo.TriangleMesh()  # type: ignore
            QuESo.IO.ReadMeshFromSTL(triangle_mesh, filename)  # type: ignore
            read_model_part_from_triangle_mesh(embedded_model_part, triangle_mesh)

        seen_model_parts: set[str] = set()
        for geometry_settings in self.analysis_settings["geometries"].values():
            model_part_name = geometry_settings["model_part_name"]
            if model_part_name not in seen_model_parts:
                model_part = self.model.GetModelPart(model_part_name)
                remove_all_elements(model_part)
                remove_all_conditions(model_part)
                seen_model_parts.add(model_part_name)

        next_element_id = 1

        for geometry_name, geometry_settings in self.analysis_settings[
            "geometries"
        ].items():
            model_part = self.model.GetModelPart(geometry_settings["model_part_name"])
            grid_settings = geometry_settings["background_grid_settings"]
            bounds_xyz = (
                grid_settings["lower_bound_xyz"],
                grid_settings["upper_bound_xyz"],
            )
            bounds_uvw = (
                grid_settings["lower_bound_uvw"],
                grid_settings["upper_bound_uvw"],
            )

            for queso_component in geometry_settings["queso_components"]:
                component_name = str(queso_component["component_name"])
                component_model_part = model_part.GetSubModelPart(component_name)
                element_settings = queso_component["element_settings"]
                active_condition_ids = queso_component["condition_settings"].get(
                    "active_condition_ids", []
                )
                next_element_id = add_elements(
                    geometry_model_part=model_part,
                    component_model_part=component_model_part,
                    elements=self.queso_model.elements(component_name),
                    property_id=int(element_settings["property_id"]),
                    nurbs_volume_name=geometry_name,
                    first_element_id=next_element_id,
                    element_name=str(
                        element_settings.get(
                            "element_type",
                            "UpdatedLagrangianElement3D8N",
                        )
                    ),
                )
                property_id_by_condition_id = self._create_condition_property_ids(
                    model_part,
                    component_name,
                    active_condition_ids,
                )
                created_loads = add_conditions(
                    geometry_model_part=model_part,
                    component_model_part=component_model_part,
                    conditions=self.queso_model.conditions(component_name),
                    active_condition_ids=active_condition_ids,
                    property_id_by_condition_id=property_id_by_condition_id,
                    bounds_xyz=bounds_xyz,
                    bounds_uvw=bounds_uvw,
                    nurbs_volume_name=geometry_name,
                    get_geometry_grid_settings=lambda gname: self.analysis_settings[
                        "geometries"
                    ][gname]["background_grid_settings"],
                )
                self._ramped_condition_loads.update(dict(created_loads))

            KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

            if self.analysis_settings["lagrange_dofs_required"]:
                KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_X, IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION_X, model_part)  # type: ignore
                KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Y, IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION_Y, model_part)  # type: ignore
                KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Z, IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION_Z, model_part)  # type: ignore

    def ApplyBoundaryConditions(self) -> None:
        super().ApplyBoundaryConditions()
        model_part = self.model.GetModelPart(
            self.analysis_settings["root_model_part_name"]
        )
        step = model_part.ProcessInfo[KM.STEP]
        multiplier = min(max(float(step) / float(self._max_step), 0.0), 1.0)
        print("multiplier:" ,multiplier)
        total_force = 0.0
        root_model_part = model_part.GetRootModelPart()
        for condition_id, target_load in self._ramped_condition_loads.items():
            condition = root_model_part.GetCondition(condition_id)
            value_x = target_load[0] * multiplier
            value_y = target_load[1] * multiplier
            value_z = target_load[2] * multiplier

            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, value_x)
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, value_y)
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, value_z)

            total_force += math.sqrt(
                value_x * value_x + value_y * value_y + value_z * value_z
            )
        print(f"The total applied force: {total_force}N")

    def GetModelPart(self, name: str = "NurbsMesh") -> KM.ModelPart:
        geometries = self.analysis_settings["geometries"]
        if name in geometries:
            return self.model.GetModelPart(geometries[name]["model_part_name"])

        if self.model.HasModelPart(name):
            return self.model.GetModelPart(name)

        root_model_part_name = self.analysis_settings["root_model_part_name"]
        if root_model_part_name and self.model.HasModelPart(
            f"{root_model_part_name}.{name}"
        ):
            return self.model.GetModelPart(f"{root_model_part_name}.{name}")

        raise RuntimeError(
            f"ModelPart '{name}' not available. Known geometries: {list(geometries.keys())}"
        )

    def GetModelParts(self) -> dict[str, KM.ModelPart]:
        return {
            name: self.model.GetModelPart(settings["model_part_name"])
            for name, settings in self.analysis_settings["geometries"].items()
        }


Analysis = KratosAnalysis

