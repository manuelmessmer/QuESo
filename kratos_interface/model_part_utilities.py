from typing import Literal, Tuple

import KratosMultiphysics as KM
import QuESoPythonModule as QuESo

from QuESoPythonModule.kratos_interface.weak_bcs import CouplingPenalty
from QuESoPythonModule.kratos_interface.weak_bcs import LagrangeSupport
from QuESoPythonModule.kratos_interface.weak_bcs import PenaltySupport
from QuESoPythonModule.kratos_interface.weak_bcs import PressureLoad
from QuESoPythonModule.kratos_interface.weak_bcs import SurfaceLoad
from QuESoPythonModule.kratos_interface.weak_bcs import TotalLoad


Point3D = Tuple[float, float, float]
LoadVector = Tuple[float, float, float]


def _get_geometry_from_model_part(
    model_part: KM.ModelPart,
    geometry_name: str,
) -> KM.Geometry:
    """Resolve a geometry from the current model part, root, or submodelparts."""
    if model_part.HasGeometry(geometry_name):
        return model_part.GetGeometry(geometry_name)

    root_model_part = model_part.GetRootModelPart()
    if root_model_part.HasGeometry(geometry_name):
        return root_model_part.GetGeometry(geometry_name)

    for sub_model_part_name in model_part.GetSubModelPartNames():
        sub_model_part = model_part.GetSubModelPart(sub_model_part_name)
        if sub_model_part.HasGeometry(geometry_name):
            return sub_model_part.GetGeometry(geometry_name)

    raise RuntimeError(f"Geometry '{geometry_name}' not found in model part '{model_part.FullName()}'.")


def read_model_part_from_triangle_mesh(
    kratos_model_part: KM.ModelPart,
    triangle_mesh: QuESo.TriangleMesh,  # type: ignore[valid-type]
) -> None:
    """Populate a Kratos model part from a QuESo triangle mesh."""
    vertices = list(triangle_mesh.Vertices())
    id_map = {tuple(vertex): vertex_id for vertex_id, vertex in enumerate(vertices)}
    for vertex_id, vertex in enumerate(vertices):
        kratos_model_part.CreateNewNode(vertex_id + 1, vertex[0], vertex[1], vertex[2])

    for triangle_id, triangle in enumerate(triangle_mesh.Triangles()):
        node_ids = [
            id_map[tuple(triangle.p1)] + 1,
            id_map[tuple(triangle.p2)] + 1,
            id_map[tuple(triangle.p3)] + 1,
        ]
        kratos_model_part.CreateNewElement(
            "ShellThinElement3D3N",
            triangle_id + 1,
            node_ids,
            kratos_model_part.GetProperties()[1],
        )


def read_triangle_mesh_from_model_part(
    triangle_mesh: QuESo.TriangleMesh,  # type: ignore[valid-type]
    kratos_model_part: KM.ModelPart,
    entity_type: Literal["Elements", "Conditions"] = "Elements",
) -> None:
    """Populate a QuESo triangle mesh from triangle entities in a Kratos model part."""
    id_map = {node.Id: queso_id for queso_id, node in enumerate(kratos_model_part.Nodes)}
    for node in kratos_model_part.Nodes:  # type: ignore[attr-defined]
        triangle_mesh.AddVertex((node.X, node.Y, node.Z))

    if entity_type == "Elements":
        entity_list = kratos_model_part.Elements
    elif entity_type == "Conditions":
        entity_list = kratos_model_part.Conditions
    else:
        message = (
            f"Given entity type '{entity_type}' is not available. "
            "Available options are: 'Elements' and 'Conditions'."
        )
        raise ValueError(message)

    for entity in entity_list:  # type: ignore[attr-defined]
        geometry = entity.GetGeometry()
        if len(geometry) != 3:
            raise ValueError(
                "model_part_utilities.read_triangle_mesh_from_model_part only supports triangles."
            )

        node_ids = tuple(id_map[node.Id] for node in geometry)
        triangle_mesh.AddTriangle(node_ids)


def add_elements(
    geometry_model_part: KM.ModelPart,
    component_model_part: KM.ModelPart,
    elements: list[QuESo.Element],  # type: ignore[valid-type]
    property_id: int,
    nurbs_volume_name: str = "NurbsVolume",
    first_element_id: int | None = None,
    element_name: str = "UpdatedLagrangianElement3D8N",
) -> int:
    """Create QuESo elements on a geometry model part and assign them to a component submodelpart."""
    nurbs_volume = _get_geometry_from_model_part(geometry_model_part, nurbs_volume_name)
    volume_properties = geometry_model_part.GetProperties()[property_id]
    next_element_id = (
        geometry_model_part.GetRootModelPart().NumberOfElements() + 1
        if first_element_id is None
        else int(first_element_id)
    )
    created_element_ids: list[int] = []

    for element in elements:
        integration_points = [
            [point.x, point.y, point.z, point.weight]
            for point in element.GetIntegrationPoints()
            if point.weight > 0
        ]
        if not integration_points:
            continue

        quadrature_point_geometries = KM.GeometriesVector()
        nurbs_volume.CreateQuadraturePointGeometries(
            quadrature_point_geometries, 2, integration_points
        )
        geometry_model_part.CreateNewElement(
            element_name,
            next_element_id,
            quadrature_point_geometries[0],
            volume_properties,
        )
        created_element_ids.append(next_element_id)
        next_element_id += 1

    _assign_entities_to_sub_model_part(component_model_part, created_element_ids, "elements")
    return next_element_id


def add_conditions(
    geometry_model_part: KM.ModelPart,
    component_model_part: KM.ModelPart,
    conditions,
    active_condition_ids: list[int],
    property_id_by_condition_id: dict[int, int],
    bounds_xyz: Tuple[Point3D, Point3D],
    bounds_uvw: Tuple[Point3D, Point3D],
    nurbs_volume_name: str = "NurbsVolume",
    get_geometry_grid_settings=None,
) -> list[tuple[int, LoadVector]]:
    """Create active QuESo conditions on a geometry model part and assign them to a component submodelpart."""
    conditions_by_id = {
        condition.GetSettings().GetInt("condition_id"): condition
        for condition in conditions
    }
    active_conditions = _collect_active_conditions(conditions_by_id, active_condition_ids)
    boundary_conditions = []
    created_loads: list[tuple[int, LoadVector]] = []

    for boundary_condition in active_conditions:
        if isinstance(boundary_condition, dict):
            condition_id = int(boundary_condition["condition_id"])
            if "segments" not in boundary_condition:
                raise RuntimeError("Boundary condition definition must contain 'segments'.")
            if boundary_condition.get("condition_type") == "TotalLoadCondition":
                weak_bc = TotalLoad(
                    [segment.GetTriangleMesh() for segment in boundary_condition["segments"]],
                    bounds_xyz,
                    bounds_uvw,
                    float(boundary_condition["modulus"]),
                    boundary_condition["direction"],
                    property_id_by_condition_id[condition_id],
                    nurbs_volume_name,
                )
                boundary_conditions.append(weak_bc)
                continue

            for segment in boundary_condition["segments"]:
                weak_bc = _create_bc_from_settings(
                    boundary_condition,
                    segment,
                    property_id_by_condition_id[condition_id],
                    nurbs_volume_name,
                )
                weak_bc.bounds_xyz = bounds_xyz
                weak_bc.bounds_uvw = bounds_uvw
                boundary_conditions.append(weak_bc)
            continue

        if getattr(boundary_condition, "is_weak_condition", lambda: False)():
            condition_settings = boundary_condition.GetSettings()
            source_condition_id = condition_settings.GetInt("condition_id")
            type_name = condition_settings.GetString("condition_type")

            if type_name == "TotalLoadCondition":
                boundary_conditions.append(
                    TotalLoad(
                        [segment.GetTriangleMesh() for segment in boundary_condition.GetSegments()],
                        bounds_xyz,
                        bounds_uvw,
                        (
                            condition_settings.GetDouble("modulus")
                            if condition_settings.IsSet("modulus")
                            else 0.0
                        ),
                        (
                            condition_settings.GetDoubleVector("direction")
                            if condition_settings.IsSet("direction")
                            else [0.0, 0.0, 0.0]
                        ),
                        property_id_by_condition_id[source_condition_id],
                        nurbs_volume_name,
                    )
                )
                continue

            for segment in boundary_condition.GetSegments():
                if type_name == "CouplingPenaltyCondition":
                    if get_geometry_grid_settings is None:
                        raise RuntimeError(
                            "CouplingPenaltyCondition requires geometry grid settings lookup."
                        )
                    if not condition_settings.IsSet("coupling_partner"):
                        raise RuntimeError(
                            "CouplingPenaltyCondition must define 'coupling_partner'."
                        )

                    coupling_partner = condition_settings.GetString("coupling_partner")
                    partner_grid_settings = get_geometry_grid_settings(coupling_partner)
                    partner_bounds_xyz = (
                        partner_grid_settings["lower_bound_xyz"],
                        partner_grid_settings["upper_bound_xyz"],
                    )
                    partner_bounds_uvw = (
                        partner_grid_settings["lower_bound_uvw"],
                        partner_grid_settings["upper_bound_uvw"],
                    )

                    boundary_conditions.append(
                        CouplingPenalty(
                            segment.GetTriangleMesh(),
                            segment.GetTriangleMesh(),
                            bounds_xyz,
                            bounds_uvw,
                            partner_bounds_xyz,
                            partner_bounds_uvw,
                            penalty_factor=(
                                condition_settings.GetDouble("penalty_factor")
                                if condition_settings.IsSet("penalty_factor")
                                else 0.0
                            ),
                            property_id=property_id_by_condition_id[source_condition_id],
                            slip=(
                                condition_settings.GetBool("slip")
                                if condition_settings.IsSet("slip")
                                else False
                            ),
                            master_nurbs_volume_name=nurbs_volume_name,
                            slave_nurbs_volume_name=coupling_partner,
                        )
                    )
                    continue

                weak_bc = _create_bc_from_settings(
                    {
                        "condition_type": type_name,
                        "value": (
                            condition_settings.GetDoubleVector("value")
                            if condition_settings.IsSet("value")
                            else [0.0, 0.0, 0.0]
                        ),
                        "penalty_factor": (
                            condition_settings.GetDouble("penalty_factor")
                            if condition_settings.IsSet("penalty_factor")
                            else 0.0
                        ),
                        "modulus": (
                            condition_settings.GetDouble("modulus")
                            if condition_settings.IsSet("modulus")
                            else 0.0
                        ),
                        "direction": (
                            condition_settings.GetDoubleVector("direction")
                            if condition_settings.IsSet("direction")
                            else [0.0, 0.0, 0.0]
                        ),
                    },
                    segment,
                    property_id_by_condition_id[source_condition_id],
                    nurbs_volume_name,
                )
                weak_bc.bounds_xyz = bounds_xyz
                weak_bc.bounds_uvw = bounds_uvw
                boundary_conditions.append(weak_bc)
            continue

        boundary_conditions.append(boundary_condition)

    first_condition_id = geometry_model_part.GetRootModelPart().NumberOfConditions() + 1
    for boundary_condition in boundary_conditions:
        created_loads.extend(boundary_condition.apply(geometry_model_part))
    last_condition_id = geometry_model_part.GetRootModelPart().NumberOfConditions()
    created_condition_ids = list(range(first_condition_id, last_condition_id + 1))
    _assign_entities_to_sub_model_part(
        component_model_part,
        created_condition_ids,
        "conditions",
    )
    return created_loads


def remove_all_elements(kratos_model_part: KM.ModelPart) -> None:
    """Remove all elements from a Kratos model part."""
    for element in kratos_model_part.Elements:
        element.Set(KM.TO_ERASE, True)
    kratos_model_part.RemoveElements(KM.TO_ERASE)


def remove_all_conditions(kratos_model_part: KM.ModelPart) -> None:
    """Remove all conditions from a Kratos model part."""
    for condition in kratos_model_part.Conditions:
        condition.Set(KM.TO_ERASE, True)
    kratos_model_part.RemoveConditions(KM.TO_ERASE)


def remove_all_nodes(kratos_model_part: KM.ModelPart) -> None:
    """Remove all nodes from a Kratos model part."""
    for node in kratos_model_part.Nodes:
        node.Set(KM.TO_ERASE, True)
    kratos_model_part.RemoveNodes(KM.TO_ERASE)


def _create_bc_from_settings(
    condition_settings,
    segment,
    property_id: int,
    nurbs_volume_name: str = "NurbsVolume",
):
    type_name = condition_settings["condition_type"]
    if type_name == "PenaltySupportCondition":
        return PenaltySupport(
            segment.GetTriangleMesh(),
            None,
            None,
            condition_settings["value"],
            condition_settings["penalty_factor"],
            property_id,
            nurbs_volume_name,
        )
    if type_name == "LagrangeSupportCondition":
        return LagrangeSupport(
            segment.GetTriangleMesh(),
            None,
            None,
            condition_settings["value"],
            property_id,
            nurbs_volume_name,
        )
    if type_name == "SurfaceLoadCondition":
        return SurfaceLoad(
            segment.GetTriangleMesh(),
            None,
            None,
            condition_settings["modulus"],
            condition_settings["direction"],
            property_id,
            nurbs_volume_name,
        )
    if type_name == "TotalLoadCondition":
        return TotalLoad(
            segment.GetTriangleMesh(),
            None,
            None,
            condition_settings["modulus"],
            condition_settings["direction"],
            property_id,
            nurbs_volume_name,
        )
    if type_name == "PressureLoadCondition":
        return PressureLoad(
            segment.GetTriangleMesh(),
            None,
            None,
            condition_settings["modulus"],
            property_id,
            nurbs_volume_name,
        )

    available_conditions = (
        "PenaltySupportCondition",
        "LagrangeSupportCondition",
        "SurfaceLoadCondition",
        "TotalLoadCondition",
        "PressureLoadCondition",
    )
    options = ", ".join(f'"{condition}"' for condition in available_conditions)
    raise ValueError(
        f"Given condition type '{type_name}' is not available. Available options are: {options}."
    )


def _collect_active_conditions(conditions_by_id: dict[int, object], active_condition_ids: list[int]) -> list[object]:
    return [conditions_by_id[int(condition_id)] for condition_id in active_condition_ids]


def _assign_entities_to_sub_model_part(
    component_model_part: KM.ModelPart,
    entity_ids: list[int],
    entity_type: Literal["elements", "conditions"],
) -> None:
    if not entity_ids:
        return
    if entity_type == "elements":
        component_model_part.AddElements(entity_ids)
        root_model_part = component_model_part.GetRootModelPart()
        node_ids = {
            node.Id
            for entity_id in entity_ids
            for node in root_model_part.GetElement(entity_id).GetGeometry()
        }
        if node_ids:
            component_model_part.AddNodes(sorted(node_ids))
        return
    component_model_part.AddConditions(entity_ids)
    root_model_part = component_model_part.GetRootModelPart()
    node_ids = {
        node.Id
        for entity_id in entity_ids
        for node in root_model_part.GetCondition(entity_id).GetGeometry()
    }
    if node_ids:
        component_model_part.AddNodes(sorted(node_ids))
