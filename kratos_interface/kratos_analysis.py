from pathlib import Path

import KratosMultiphysics as KM
import QuESoPythonModule as QuESo

from QuESoPythonModule.kratos_interface.custom_analysis_stage import CustomAnalysisStage


class _TriangleMeshSegment:
    def __init__(self, triangle_mesh) -> None:
        self._triangle_mesh = triangle_mesh

    def GetTriangleMesh(self):
        return self._triangle_mesh


class KratosAnalysis:
    def __init__(self, pyqueso, kratos_settings_filename: str = "KratosParameters.json") -> None:
        self.pyqueso = pyqueso
        self.kratos_settings_filename = kratos_settings_filename
        self.model = None
        self.root_model_part_name = None
        self.mesh_model_part_map = {}

    def _resolve_path(self, filename: str) -> str:
        path = Path(filename)
        if path.is_absolute():
            return str(path)
        return str((self.pyqueso._settings_dir / path).resolve())

    def _build_analysis_definitions(self):
        settings = self.pyqueso.GetSettings()
        models = settings.get("models", [])
        if not isinstance(models, list) or not models:
            raise RuntimeError("QuESoSettings must contain a non-empty 'models' list.")

        root_model_part_name = str(
            settings.get("shared_analysis_settings", {}).get("root_model_part_name", "Structure")
        )

        mesh_definitions = []
        boundary_condition_definitions = {}
        embedded_model_parts = []

        for model_entry in models:
            model_name = str(model_entry.get("name", "")).strip()
            if not model_name:
                raise RuntimeError("Each model must define a non-empty 'name'.")

            analysis_settings = model_entry.get("analysis_settings", {})
            if "kratos_model_part_name" not in analysis_settings:
                raise RuntimeError(f"Model '{model_name}' misses analysis_settings.kratos_model_part_name.")
            if "geometry_name" not in analysis_settings:
                raise RuntimeError(f"Model '{model_name}' misses analysis_settings.geometry_name.")
            if "element_groups" not in analysis_settings:
                raise RuntimeError(f"Model '{model_name}' misses analysis_settings.element_groups.")
            if "condition_ids" not in analysis_settings:
                raise RuntimeError(f"Model '{model_name}' misses analysis_settings.condition_ids.")

            queso_settings = model_entry.get("queso", {})
            grid_settings = queso_settings.get("background_grid_settings", {})
            required_grid_keys = [
                "lower_bound_xyz",
                "upper_bound_xyz",
                "lower_bound_uvw",
                "upper_bound_uvw",
                "number_of_elements",
                "polynomial_order",
            ]
            for key in required_grid_keys:
                if key not in grid_settings:
                    raise RuntimeError(f"Model '{model_name}' misses queso.background_grid_settings.{key}.")

            conditions = queso_settings.get("conditions_settings_list", [])
            by_index = {}
            for cond in conditions:
                if "condition_id" not in cond:
                    raise RuntimeError(f"Model '{model_name}' has condition without condition_id.")
                idx = int(cond["condition_id"])
                if idx in by_index:
                    raise RuntimeError(f"Model '{model_name}' has duplicate condition_id {idx}.")
                by_index[idx] = cond

            model_collection = f"model::{model_name}"
            boundary_condition_definitions[model_collection] = []
            for cond_idx in analysis_settings.get("condition_ids", []):
                idx = int(cond_idx)
                if idx not in by_index:
                    raise RuntimeError(f"Model '{model_name}' references unknown condition_id {idx}.")
                bc_def = dict(by_index[idx])
                if "input_filename" not in bc_def:
                    raise RuntimeError(f"Model '{model_name}' condition_id {idx} misses input_filename.")
                tri_mesh = QuESo.TriangleMesh()  # type: ignore
                QuESo.IO.ReadMeshFromSTL(tri_mesh, self._resolve_path(str(bc_def["input_filename"])))  # type: ignore
                bc_def["segments"] = [_TriangleMeshSegment(tri_mesh)]
                boundary_condition_definitions[model_collection].append(bc_def)

            element_groups = []
            for group in analysis_settings.get("element_groups", []):
                source_name = str(group.get("source_name", "")).strip()
                if not source_name:
                    raise RuntimeError(f"Model '{model_name}' has element_group without source_name.")
                if source_name != model_name:
                    raise RuntimeError(
                        f"Model '{model_name}' has source_name '{source_name}', "
                        "but source_name must match the model name exactly."
                    )
                elements = self.pyqueso.GetElements(model_name)
                element_groups.append(
                    {
                        "source_name": source_name,
                        "property_id": int(group["property_id"]),
                        "element_type": str(group.get("element_type", "UpdatedLagrangianElement3D8N")),
                        "elements": elements,
                    }
                )

            mesh_definitions.append(
                {
                    "name": model_name,
                    "kratos_model_part_name": str(analysis_settings["kratos_model_part_name"]),
                    "geometry_name": str(analysis_settings["geometry_name"]),
                    "grid_settings": grid_settings,
                    "element_groups": element_groups,
                    "boundary_condition_collections": [model_collection],
                }
            )

            input_filename = model_entry.get("queso", {}).get("general_settings", {}).get("input_filename")
            if input_filename is not None:
                embedded_model_parts.append(
                    {
                        "model_name": model_name,
                        "filename": self._resolve_path(str(input_filename)),
                    }
                )

        coupling_definitions = settings.get("shared_analysis_settings", {}).get("couplings", [])
        return (
            root_model_part_name,
            mesh_definitions,
            boundary_condition_definitions,
            coupling_definitions,
            embedded_model_parts,
        )

    def Run(self) -> None:
        if not self.pyqueso._embedded_models:
            self.pyqueso.Run()

        (
            root_model_part_name,
            mesh_definitions,
            boundary_condition_definitions,
            coupling_definitions,
            embedded_model_parts,
        ) = self._build_analysis_definitions()

        self.model = KM.Model()
        any_embedded = next(iter(self.pyqueso._embedded_models.values()), None)
        if any_embedded is None:
            raise RuntimeError("No embedded models available. Call PyQuESo.Run() first.")
        analysis_stage = CustomAnalysisStage(
            self.model,
            any_embedded.GetSettings(),
            self.kratos_settings_filename,
            root_model_part_name,
            mesh_definitions,
            boundary_condition_definitions,
            coupling_definitions,
            embedded_model_parts,
        )
        analysis_stage.Run()
        self.root_model_part_name = analysis_stage.root_model_part_name
        self.mesh_model_part_map = analysis_stage.mesh_model_part_map

    def GetModelPart(self, name: str = "NurbsMesh") -> KM.ModelPart:
        if self.model is None:
            raise RuntimeError("KratosAnalysis has not been run.")
        if name in self.mesh_model_part_map:
            return self.model.GetModelPart(self.mesh_model_part_map[name])

        if self.model.HasModelPart(name):
            return self.model.GetModelPart(name)

        if self.root_model_part_name and self.model.HasModelPart(f"{self.root_model_part_name}.{name}"):
            return self.model.GetModelPart(f"{self.root_model_part_name}.{name}")

        raise RuntimeError(
            f"ModelPart '{name}' not available. Known meshes: {list(self.mesh_model_part_map.keys())}"
        )

    def GetModelParts(self):
        if self.model is None:
            raise RuntimeError("KratosAnalysis has not been run.")
        return {name: self.model.GetModelPart(full_name) for name, full_name in self.mesh_model_part_map.items()}
