# Import QuESo
import QuESoPythonModule as QuESo

# Import Kratos
import math
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IgaApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from QuESoPythonModule.kratos_interface.model_part_utilities import ModelPartUtilities
from QuESoPythonModule.kratos_interface.weak_bcs import CouplingPenalty
import KratosMultiphysics.LinearSolversApplication


class CustomAnalysisStage(StructuralMechanicsAnalysis):
    def __init__(
        self,
        model: KM.Model,
        queso_settings: QuESo.Dictionary, # type: ignore
        kratos_settings_filename: str,
        root_model_part_name: str,
        mesh_definitions,
        boundary_condition_definitions,
        coupling_definitions,
        embedded_model_parts,
    ) -> None:
        with open(kratos_settings_filename, "r") as parameter_file:
            analysis_parameters = KM.Parameters(parameter_file.read())

        self.queso_settings = queso_settings
        self.mesh_definitions = mesh_definitions
        self.boundary_condition_definitions = boundary_condition_definitions
        self.coupling_definitions = coupling_definitions
        self.embedded_model_parts = embedded_model_parts
        self.lagrange_dofs_required = self._HasLagrangeSupportCondition()

        self.root_model_part_name = root_model_part_name
        if not self.root_model_part_name:
            self.root_model_part_name = "Structure"
        analysis_parameters["solver_settings"]["model_part_name"].SetString(self.root_model_part_name)

        if not model.HasModelPart(self.root_model_part_name):
            model.CreateModelPart(self.root_model_part_name)

        self.mesh_model_part_map = {}
        for mesh in self.mesh_definitions:
            mesh_name = mesh["name"]
            full_name = mesh.get("kratos_model_part_name", f"{self.root_model_part_name}.{mesh_name}")
            self.mesh_model_part_map[mesh_name] = full_name

            mesh_model_part = model.GetModelPart(full_name) if model.HasModelPart(full_name) else model.CreateModelPart(full_name)
            mesh_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            mesh_model_part.AddNodalSolutionStepVariable(KM.REACTION)
            if self.lagrange_dofs_required:
                mesh_model_part.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)
                mesh_model_part.AddNodalSolutionStepVariable(IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION) # type: ignore

            for element_group in mesh.get("element_groups", []):
                group_name = str(element_group.get("source_name", "")).strip()
                if group_name and not mesh_model_part.HasSubModelPart(group_name):
                    mesh_model_part.CreateSubModelPart(group_name)

        self.coupling_model_part_name = f"{self.root_model_part_name}.Coupling"
        if not model.HasModelPart(self.coupling_model_part_name):
            coupling_model_part = model.CreateModelPart(self.coupling_model_part_name)
            coupling_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            coupling_model_part.AddNodalSolutionStepVariable(KM.REACTION)

        self._OverrideNurbsModelerParameters(analysis_parameters)
        self._initialized = False

        self.Initialized = False
        super().__init__(model, analysis_parameters)

    def _HasLagrangeSupportCondition(self) -> bool:
        for bc_collection in self.boundary_condition_definitions.values():
            for bc in bc_collection:
                if bc.get("condition_type", "") == "LagrangeSupportCondition":
                    return True
        return False

    def _init_boundary_conditions(self):
        if not self._initialized: 
            self._condition_loads = []
            model_part = self.model.GetModelPart('Structure')
            for condition in model_part.Conditions:
                if condition.Has(KM.StructuralMechanicsApplication.POINT_LOAD_X):
                    value_x = condition.GetValue(KM.StructuralMechanicsApplication.POINT_LOAD_X)
                    value_y = condition.GetValue(KM.StructuralMechanicsApplication.POINT_LOAD_Y)
                    value_z = condition.GetValue(KM.StructuralMechanicsApplication.POINT_LOAD_Z)
                    self._condition_loads.append((value_x, value_y, value_z))
            self._initialized = True


    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()
        self._init_boundary_conditions()
        model_part = self.model.GetModelPart('Structure')
        step = model_part.ProcessInfo[KratosMultiphysics.STEP]
        mulitplier = step / 5
        total_force = 0.0
        for i, condition in enumerate(model_part.Conditions):

            if condition.Has(KM.StructuralMechanicsApplication.POINT_LOAD_X):
                value_x = self._condition_loads[i][0] * mulitplier
                value_y = self._condition_loads[i][1] * mulitplier
                value_z = self._condition_loads[i][2] * mulitplier

                condition.SetValue(KM.StructuralMechanicsApplication.POINT_LOAD_X, value_x) 
                condition.SetValue(KM.StructuralMechanicsApplication.POINT_LOAD_Y, value_y) 
                condition.SetValue(KM.StructuralMechanicsApplication.POINT_LOAD_Z, value_z) 

                total_force += math.sqrt( value_x*value_x + value_y*value_y + value_z*value_z )
        print(f"The total applied force: {total_force}N")

    def _OverrideNurbsModelerParameters(self, analysis_parameters: KM.Parameters) -> None:
        mesh_by_name = {mesh["name"]: mesh for mesh in self.mesh_definitions}
        mesh_by_full_name = {}
        mesh_by_full_and_geometry = {}
        for mesh in self.mesh_definitions:
            full_name = mesh.get("kratos_model_part_name", f"{self.root_model_part_name}.{mesh['name']}")
            geometry_name = mesh.get("geometry_name", "NurbsVolume")
            mesh_by_full_name.setdefault(full_name, []).append(mesh)
            mesh_by_full_and_geometry[(full_name, geometry_name)] = mesh

        for modeler in analysis_parameters["modelers"].values():
            if modeler["modeler_name"].GetString() != "NurbsGeometryModeler":
                continue

            parameters = modeler["Parameters"]
            model_part_name = parameters["model_part_name"].GetString()
            requested_geometry_name = parameters["geometry_name"].GetString() if parameters.Has("geometry_name") else "NurbsVolume"

            mesh = mesh_by_name.get(model_part_name)
            if mesh is None:
                mesh = mesh_by_full_and_geometry.get((model_part_name, requested_geometry_name))
            if mesh is None:
                candidates = mesh_by_full_name.get(model_part_name, [])
                if len(candidates) == 1:
                    mesh = candidates[0]
            if mesh is None:
                continue

            grid_settings = mesh["grid_settings"]

            if parameters.Has("lower_point_xyz"):
                parameters.RemoveValue("lower_point_xyz")
            parameters.AddEmptyValue("lower_point_xyz")
            parameters["lower_point_xyz"].SetVector(grid_settings["lower_bound_xyz"])

            if parameters.Has("upper_point_xyz"):
                parameters.RemoveValue("upper_point_xyz")
            parameters.AddEmptyValue("upper_point_xyz")
            parameters["upper_point_xyz"].SetVector(grid_settings["upper_bound_xyz"])

            if parameters.Has("lower_point_uvw"):
                parameters.RemoveValue("lower_point_uvw")
            parameters.AddEmptyValue("lower_point_uvw")
            parameters["lower_point_uvw"].SetVector(grid_settings["lower_bound_uvw"])

            if parameters.Has("upper_point_uvw"):
                parameters.RemoveValue("upper_point_uvw")
            parameters.AddEmptyValue("upper_point_uvw")
            parameters["upper_point_uvw"].SetVector(grid_settings["upper_bound_uvw"])

            if parameters.Has("polynomial_order"):
                parameters.RemoveValue("polynomial_order")
            parameters.AddEmptyValue("polynomial_order")
            # if "polynomial_order" in grid_settings:
            #     parameters["polynomial_order"].SetVector(grid_settings["polynomial_order"])
            # else:
            parameters["polynomial_order"].SetVector([1, 1, 1])

            if parameters.Has("number_of_knot_spans"):
                parameters.RemoveValue("number_of_knot_spans")
            parameters.AddEmptyValue("number_of_knot_spans")
            parameters["number_of_knot_spans"].SetVector(grid_settings["number_of_elements"])

            geometry_name = mesh.get("geometry_name", "NurbsVolume")
            if parameters.Has("geometry_name"):
                parameters["geometry_name"].SetString(geometry_name)
            else:
                parameters.AddEmptyValue("geometry_name")
                parameters["geometry_name"].SetString(geometry_name)


    def _ModelersSetupModelPart(self) -> None:
        for embedded_part in self.embedded_model_parts:
            raise RuntimeError(embedded_part)
            model_name = str(embedded_part.get("model_name", "")).strip()
            filename = str(embedded_part.get("filename", "")).strip()
            if not model_name or not filename:
                continue

            embedded_model_part_name = f"embedded_{model_name}"
            if self.model.HasModelPart(embedded_model_part_name):
                continue

            embedded_model_part = self.model.CreateModelPart(embedded_model_part_name)
            embedded_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            embedded_model_part.AddNodalSolutionStepVariable(KM.REACTION)
            embedded_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

            triangle_mesh = QuESo.TriangleMesh() # type: ignore
            QuESo.IO.ReadMeshFromSTL(triangle_mesh, filename) # type: ignore
            ModelPartUtilities.read_model_part_from_triangle_mesh(embedded_model_part, triangle_mesh)

        super()._ModelersSetupModelPart()

    def ModifyInitialGeometry(self) -> None:
        root_model_part = self.model.GetModelPart(self.root_model_part_name)
        next_element_id = root_model_part.NumberOfElements() + 1

        for mesh in self.mesh_definitions:
            model_part = self.model.GetModelPart(self.mesh_model_part_map[mesh["name"]])
            geometry_name = mesh.get("geometry_name", "NurbsVolume")
            ModelPartUtilities.remove_all_elements(model_part)
            ModelPartUtilities.remove_all_conditions(model_part)

            for element_group in mesh["element_groups"]:
                next_element_id = ModelPartUtilities.add_elements_to_model_part(
                    model_part,
                    element_group["elements"],
                    int(element_group["property_id"]),
                    geometry_name,
                    next_element_id,
                    str(element_group.get("element_type", "UpdatedLagrangianElement3D8N")),
                )
                next_element_id += 1

            bounds_xyz = (
                mesh["grid_settings"]["lower_bound_xyz"],
                mesh["grid_settings"]["upper_bound_xyz"],
            )
            bounds_uvw = (
                mesh["grid_settings"]["lower_bound_uvw"],
                mesh["grid_settings"]["upper_bound_uvw"],
            )

            boundary_conditions = []
            for collection_name in mesh.get("boundary_condition_collections", []):
                if collection_name not in self.boundary_condition_definitions:
                    raise RuntimeError(f"Boundary condition collection '{collection_name}' is not defined.")
                boundary_conditions.extend(self.boundary_condition_definitions[collection_name])

            ModelPartUtilities.add_conditions_to_model_part(
                model_part,
                boundary_conditions,
                bounds_xyz,
                bounds_uvw,
                geometry_name,
            )

            KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

            if self.lagrange_dofs_required:
                KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_X, IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION_X, model_part) # type: ignore
                KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Y, IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION_Y, model_part) # type: ignore
                KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Z, IgaApplication.VECTOR_LAGRANGE_MULTIPLIER_REACTION_Z, model_part) # type: ignore

        self._ValidateMeshIndependenceAndProperties()
        self._CreateCouplingConditions()

    def _ModelersSetupModelPart(self) -> None:
        for embedded_part in self.embedded_model_parts:
            filename = str(embedded_part.get("filename", "")).strip()
            model_name = str(embedded_part.get("model_name", "")).strip()
            if not filename or not model_name:
                continue
            mp_name = f"embedded_{model_name}"
            if self.model.HasModelPart(mp_name):
                continue

            embedded_model_part = self.model.CreateModelPart(mp_name)
            embedded_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            embedded_model_part.AddNodalSolutionStepVariable(KM.REACTION)
            embedded_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

            triangle_mesh = QuESo.TriangleMesh() # type: ignore
            QuESo.IO.ReadMeshFromSTL(triangle_mesh, filename) # type: ignore
            ModelPartUtilities.read_model_part_from_triangle_mesh(embedded_model_part, triangle_mesh)

        super()._ModelersSetupModelPart()

    def _ValidateMeshIndependenceAndProperties(self) -> None:
        mesh_names = [mesh["name"] for mesh in self.mesh_definitions]
        model_part_to_allowed_ids = {}
        for mesh in self.mesh_definitions:
            mp_name = self.mesh_model_part_map[mesh["name"]]
            allowed_property_ids = {int(group["property_id"]) for group in mesh.get("element_groups", [])}
            if mp_name not in model_part_to_allowed_ids:
                model_part_to_allowed_ids[mp_name] = set()
            model_part_to_allowed_ids[mp_name].update(allowed_property_ids)

        for mp_name, allowed_property_ids in model_part_to_allowed_ids.items():
            mp = self.model.GetModelPart(mp_name)
            violating_ids = set()
            for element in mp.Elements:
                prop_id = int(element.Properties.Id)
                if prop_id not in allowed_property_ids:
                    violating_ids.add(prop_id)
            if violating_ids:
                raise RuntimeError(
                    f"ModelPart '{mp_name}' contains elements with unexpected property ids {sorted(violating_ids)}. "
                    f"Allowed ids are {sorted(allowed_property_ids)}."
                )

        for i in range(len(mesh_names)):
            for j in range(i + 1, len(mesh_names)):
                mp_i = self.model.GetModelPart(self.mesh_model_part_map[mesh_names[i]])
                mp_j = self.model.GetModelPart(self.mesh_model_part_map[mesh_names[j]])
                if mp_i.FullName() == mp_j.FullName():
                    continue
                node_ids_i = {node.Id for node in mp_i.Nodes}
                node_ids_j = {node.Id for node in mp_j.Nodes}
                shared_ids = node_ids_i.intersection(node_ids_j)
                if shared_ids:
                    preview = sorted(shared_ids)[:20]
                    raise RuntimeError(
                        f"Meshes '{mesh_names[i]}' and '{mesh_names[j]}' share node IDs. "
                        f"Expected disjoint IDs. shared_count={len(shared_ids)}, preview={preview}"
                    )

    def _CreateCouplingConditions(self) -> None:
        if not self.coupling_definitions:
            return

        for coupling in self.coupling_definitions:
            master_mp = self.model.GetModelPart(self.mesh_model_part_map[coupling["master_mesh"]])
            slave_mp = self.model.GetModelPart(self.mesh_model_part_map[coupling["slave_mesh"]])

            master_mesh = next(mesh for mesh in self.mesh_definitions if mesh["name"] == coupling["master_mesh"])
            slave_mesh = next(mesh for mesh in self.mesh_definitions if mesh["name"] == coupling["slave_mesh"])
            master_geom_name = coupling.get("master_geometry_name", master_mesh.get("geometry_name", "NurbsVolume"))
            slave_geom_name = coupling.get("slave_geometry_name", slave_mesh.get("geometry_name", "NurbsVolume"))
            master_bounds_xyz = (
                master_mesh["grid_settings"]["lower_bound_xyz"],
                master_mesh["grid_settings"]["upper_bound_xyz"],
            )
            master_bounds_uvw = (
                master_mesh["grid_settings"]["lower_bound_uvw"],
                master_mesh["grid_settings"]["upper_bound_uvw"],
            )
            slave_bounds_xyz = (
                slave_mesh["grid_settings"]["lower_bound_xyz"],
                slave_mesh["grid_settings"]["upper_bound_xyz"],
            )
            slave_bounds_uvw = (
                slave_mesh["grid_settings"]["lower_bound_uvw"],
                slave_mesh["grid_settings"]["upper_bound_uvw"],
            )

            master_triangles = QuESo.TriangleMesh() # type: ignore
            slave_triangles = QuESo.TriangleMesh() # type: ignore
            QuESo.IO.ReadMeshFromSTL(master_triangles, coupling["master_interface_stl"]) # type: ignore
            QuESo.IO.ReadMeshFromSTL(slave_triangles, coupling["slave_interface_stl"]) # type: ignore
            coupling_bc = CouplingPenalty(
                master_triangles,
                slave_triangles,
                master_bounds_xyz,
                master_bounds_uvw,
                slave_bounds_xyz,
                slave_bounds_uvw,
                penalty_factor=float(coupling["penalty_factor"]),
                slip=bool(coupling.get("slip", False)),
                property_id=int(coupling.get("property_id", 1000)),
                master_nurbs_volume_name=master_geom_name,
                slave_nurbs_volume_name=slave_geom_name,
            )
            # Apply on master model part (and its root) where both geometries exist.
            coupling_bc.apply(master_mp)
