"""Utilities for transferring QuESo entities into Kratos model parts."""

from collections.abc import Sequence
from typing import Any, Callable

import KratosMultiphysics as KM
import pyqueso

from .weak_bcs import (
    CouplingPenalty,
    LagrangeSupport,
    PenaltySupport,
    PressureLoad,
    SurfaceLoad,
)

Bounds = tuple[Sequence[float], Sequence[float]]


def embedded_model_part_name(component_name: str) -> str:
    """Return the Kratos model-part name used for one source volume mesh."""
    return f"EmbeddedModelPart_{component_name}"


def add_triangle_mesh_to_model_part(model_part: Any, triangle_mesh: Any) -> None:
    """Populate a Kratos model part with shell elements from a QuESo triangle mesh."""
    root = model_part.GetRootModelPart()
    if not root.HasProperties(1):
        root.CreateNewProperties(1)
    vertices = list(triangle_mesh.vertices())
    vertex_ids = {vertex: index + 1 for index, vertex in enumerate(vertices)}
    for vertex_id, vertex in enumerate(vertices, start=1):
        model_part.CreateNewNode(vertex_id, *vertex)
    for triangle_id, triangle in enumerate(triangle_mesh.triangles(), start=1):
        node_ids = [
            vertex_ids[triangle.p1],
            vertex_ids[triangle.p2],
            vertex_ids[triangle.p3],
        ]
        model_part.CreateNewElement(
            "ShellThinElement3D3N",
            triangle_id,
            node_ids,
            root.GetProperties()[1],
        )


def read_triangle_mesh_to_model_part(model_part: Any, filename: str) -> None:
    """Read an STL source volume and populate a Kratos model part from it."""
    triangle_mesh = pyqueso.mesh.TriangleMesh()
    pyqueso.io.read_mesh_from_stl(triangle_mesh, filename)
    add_triangle_mesh_to_model_part(model_part, triangle_mesh)


def ensure_model_part(model: Any, full_name: str) -> Any:
    """Create a full Kratos model-part path and return its final model part."""
    segments = full_name.split(".")
    root_name = segments[0]
    if model.HasModelPart(root_name):
        model_part = model.GetModelPart(root_name)
    else:
        model_part = model.CreateModelPart(root_name)
    current_name = root_name
    for segment in segments[1:]:
        current_name = f"{current_name}.{segment}"
        if model.HasModelPart(current_name):
            model_part = model.GetModelPart(current_name)
        else:
            model_part = model_part.CreateSubModelPart(segment)
    return model_part


def add_elements(
    geometry_model_part: Any,
    component_model_part: Any,
    elements: Any,
    property_id: int,
    geometry_name: str,
    element_type: str,
) -> None:
    """Create QuESo volume elements and assign them to a component submodelpart."""
    volume = geometry_model_part.GetGeometry(geometry_name)
    properties = geometry_model_part.GetProperties()[property_id]
    root = geometry_model_part.GetRootModelPart()
    created_ids: list[int] = []
    for element in elements:
        integration_points = [
            [point.x, point.y, point.z, point.weight]
            for point in element.integration_points
            if point.weight > 0.0
        ]
        if not integration_points:
            continue
        geometries = KM.GeometriesVector()
        volume.CreateQuadraturePointGeometries(geometries, 2, integration_points)
        element_id = root.NumberOfElements() + 1
        geometry_model_part.CreateNewElement(
            element_type, element_id, geometries[0], properties
        )
        created_ids.append(element_id)
    _assign_entities(component_model_part, created_ids, "elements")


def add_conditions(
    geometry_model_part: Any,
    component_model_part: Any,
    conditions: Any,
    bounds_xyz: Bounds,
    bounds_uvw: Bounds,
    geometry_name: str,
    property_id_for_condition: Callable[[int], int],
    coupling_target: Callable[[str], dict[str, Any]],
) -> None:
    """Create all active ordinary and coupling conditions for one component."""
    for condition in conditions:
        settings = condition.settings
        condition_id = settings.get_int("condition_id")
        condition_type = settings.get_string("condition_type")
        property_id = property_id_for_condition(condition_id)
        created_ids: list[int] = []
        for segment in condition.segments:
            if not segment.is_in_active_element:
                continue
            common = {
                "triangle_mesh": segment.triangle_mesh,
                "bounds_xyz": bounds_xyz,
                "bounds_uvw": bounds_uvw,
                "property_id": property_id,
                "geometry_name": geometry_name,
            }
            if condition_type == "PenaltySupportCondition":
                boundary_condition = PenaltySupport(
                    **common,
                    value=settings.get_double_vector("value"),
                    penalty_factor=settings.get_double("penalty_factor"),
                )
                created_ids.extend(boundary_condition.apply(geometry_model_part))
            elif condition_type == "LagrangeSupportCondition":
                boundary_condition = LagrangeSupport(
                    **common, value=settings.get_double_vector("value")
                )
                created_ids.extend(boundary_condition.apply(geometry_model_part))
            elif condition_type == "SurfaceLoadCondition":
                boundary_condition = SurfaceLoad(
                    **common,
                    modulus=settings.get_double("modulus"),
                    direction=settings.get_double_vector("direction"),
                )
                created_ids.extend(boundary_condition.apply(geometry_model_part))
            elif condition_type == "PressureLoadCondition":
                boundary_condition = PressureLoad(
                    **common, modulus=settings.get_double("modulus")
                )
                created_ids.extend(boundary_condition.apply(geometry_model_part))
            elif condition_type == "CouplingPenaltyCondition":
                partner_name = settings.get_string("coupling_partner")
                target = coupling_target(partner_name)
                boundary_condition = CouplingPenalty(
                    **common,
                    slave_bounds_xyz=target["bounds_xyz"],
                    slave_bounds_uvw=target["bounds_uvw"],
                    slave_geometry_name=target["geometry_name"],
                    penalty_factor=settings.get_double("penalty_factor"),
                    slip=settings.get_bool("slip"),
                )
                created_ids.extend(
                    boundary_condition.apply(geometry_model_part, target["model_part"])
                )
            else:
                raise ValueError(f"Unsupported condition type '{condition_type}'.")
        _assign_entities(component_model_part, created_ids, "conditions")


def remove_all_elements(model_part: Any) -> None:
    """Remove all elements from a Kratos model part."""
    for element in model_part.Elements:
        element.Set(KM.TO_ERASE, True)
    model_part.RemoveElements(KM.TO_ERASE)


def remove_all_conditions(model_part: Any) -> None:
    """Remove all conditions from a Kratos model part."""
    for condition in model_part.Conditions:
        condition.Set(KM.TO_ERASE, True)
    model_part.RemoveConditions(KM.TO_ERASE)


def _assign_entities(
    component_model_part: Any, entity_ids: list[int], entity_type: str
) -> None:
    if not entity_ids:
        return
    root = component_model_part.GetRootModelPart()
    if entity_type == "elements":
        component_model_part.AddElements(entity_ids)
        entities = [root.GetElement(entity_id) for entity_id in entity_ids]
    else:
        component_model_part.AddConditions(entity_ids)
        entities = [root.GetCondition(entity_id) for entity_id in entity_ids]
    node_ids = sorted({node.Id for entity in entities for node in entity.GetGeometry()})
    if node_ids:
        component_model_part.AddNodes(node_ids)
