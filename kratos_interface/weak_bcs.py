"""Create Kratos weak boundary and coupling conditions from QuESo surfaces."""

from collections.abc import Sequence
from typing import Any

import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np
from pyqueso.scripts.helper import point_from_global_to_param_space

Point3D = tuple[float, float, float]
Bounds = tuple[Sequence[float], Sequence[float]]


def _triangle_geometry(
    triangle: Any, bounds_xyz: Bounds, bounds_uvw: Bounds, reverse: bool = False
) -> Any:
    points = [triangle.p1, triangle.p2, triangle.p3]
    if reverse:
        points[1], points[2] = points[2], points[1]
    parameters = [
        point_from_global_to_param_space(point, bounds_xyz, bounds_uvw)
        for point in points
    ]
    nodes = [
        KM.Node(index + 1, *parameter) for index, parameter in enumerate(parameters)
    ]
    return KM.Triangle3D3(*nodes)


def _surface_geometry(volume: Any, triangle_geometry: Any) -> Any:
    return KM.SurfaceInNurbsVolumeGeometry(volume, triangle_geometry)


def _quadrature_geometry(surface: Any) -> Any:
    quadrature_geometries = KM.GeometriesVector()
    surface.CreateQuadraturePointGeometries(quadrature_geometries, 2)
    return quadrature_geometries[0]


class WeakCondition:
    """Base data shared by Kratos weak conditions."""

    def __init__(
        self,
        triangle_mesh: Any,
        bounds_xyz: Bounds,
        bounds_uvw: Bounds,
        property_id: int,
        geometry_name: str,
    ) -> None:
        self.triangle_mesh = triangle_mesh
        self.bounds_xyz = bounds_xyz
        self.bounds_uvw = bounds_uvw
        self.property_id = property_id
        self.geometry_name = geometry_name


class PenaltySupport(WeakCondition):
    """Create penalty support conditions from one condition segment."""

    def __init__(
        self, *args: Any, value: Sequence[float], penalty_factor: float, **kwargs: Any
    ) -> None:
        super().__init__(*args, **kwargs)
        self.value = value
        self.penalty_factor = penalty_factor

    def apply(self, model_part: Any) -> list[int]:
        """Create conditions and return their root-model-part IDs."""
        properties = model_part.GetProperties()[self.property_id]
        properties.SetValue(IGA.PENALTY_FACTOR, self.penalty_factor)
        volume = model_part.GetGeometry(self.geometry_name)
        condition_ids: list[int] = []
        for triangle in self.triangle_mesh.triangles():
            if triangle.aspect_ratio() >= 1.0e8:
                continue
            surface = _surface_geometry(
                volume, _triangle_geometry(triangle, self.bounds_xyz, self.bounds_uvw)
            )
            if surface.Area() <= 1.0e-14:
                continue
            condition_id = model_part.GetRootModelPart().NumberOfConditions() + 1
            condition = model_part.CreateNewCondition(
                "SupportPenaltyCondition",
                condition_id,
                _quadrature_geometry(surface),
                properties,
            )
            condition.SetValue(KM.DISPLACEMENT, KM.Vector(self.value))
            condition_ids.append(condition_id)
        return condition_ids


class LagrangeSupport(WeakCondition):
    """Create Lagrange-multiplier support conditions."""

    def __init__(self, *args: Any, value: Sequence[float], **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.value = value

    def apply(self, model_part: Any) -> list[int]:
        """Create conditions and return their root-model-part IDs."""
        properties = model_part.GetProperties()[self.property_id]
        volume = model_part.GetGeometry(self.geometry_name)
        condition_ids: list[int] = []
        for triangle in self.triangle_mesh.triangles():
            if triangle.aspect_ratio() >= 1.0e8:
                continue
            surface = _surface_geometry(
                volume, _triangle_geometry(triangle, self.bounds_xyz, self.bounds_uvw)
            )
            if surface.Area() <= 1.0e-14:
                continue
            condition_id = model_part.GetRootModelPart().NumberOfConditions() + 1
            condition = model_part.CreateNewCondition(
                "SupportLagrangeCondition",
                condition_id,
                _quadrature_geometry(surface),
                properties,
            )
            condition.SetValue(KM.DISPLACEMENT, KM.Vector(self.value))
            condition_ids.append(condition_id)
        return condition_ids


class SurfaceLoad(WeakCondition):
    """Create area-distributed load conditions."""

    def __init__(
        self, *args: Any, modulus: float, direction: Sequence[float], **kwargs: Any
    ) -> None:
        super().__init__(*args, **kwargs)
        direction_array = np.asarray(direction, dtype=float)
        norm = np.linalg.norm(direction_array)
        if norm <= 1.0e-14:
            raise ValueError("SurfaceLoadCondition direction must be non-zero.")
        self.force_density = float(modulus) * direction_array / norm

    def apply(self, model_part: Any) -> list[int]:
        """Create conditions and return their root-model-part IDs."""
        properties = model_part.GetProperties()[self.property_id]
        volume = model_part.GetGeometry(self.geometry_name)
        condition_ids: list[int] = []
        for triangle in self.triangle_mesh.triangles():
            for point in triangle.integration_points_global(1):
                if point.weight <= 1.0e-14:
                    continue
                local_point = point_from_global_to_param_space(
                    (point.x, point.y, point.z), self.bounds_xyz, self.bounds_uvw
                )
                geometries = KM.GeometriesVector()
                volume.CreateQuadraturePointGeometries(
                    geometries, 2, [[*local_point, point.weight]]
                )
                condition_id = model_part.GetRootModelPart().NumberOfConditions() + 1
                condition = model_part.CreateNewCondition(
                    "LoadCondition", condition_id, geometries[0], properties
                )
                force = point.weight * self.force_density
                condition.SetValue(SMA.POINT_LOAD_X, float(force[0]))
                condition.SetValue(SMA.POINT_LOAD_Y, float(force[1]))
                condition.SetValue(SMA.POINT_LOAD_Z, float(force[2]))
                condition_ids.append(condition_id)
        return condition_ids


class PressureLoad(WeakCondition):
    """Create pressure load conditions."""

    def __init__(self, *args: Any, modulus: float, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.modulus = float(modulus)

    def apply(self, model_part: Any) -> list[int]:
        """Create conditions and return their root-model-part IDs."""
        properties = model_part.GetProperties()[self.property_id]
        volume = model_part.GetGeometry(self.geometry_name)
        condition_ids: list[int] = []
        for triangle in self.triangle_mesh.triangles():
            for point in triangle.integration_points_global(1):
                if point.weight <= 1.0e-14:
                    continue
                local_point = point_from_global_to_param_space(
                    (point.x, point.y, point.z), self.bounds_xyz, self.bounds_uvw
                )
                geometries = KM.GeometriesVector()
                volume.CreateQuadraturePointGeometries(
                    geometries, 2, [[*local_point, point.weight]]
                )
                condition_id = model_part.GetRootModelPart().NumberOfConditions() + 1
                condition = model_part.CreateNewCondition(
                    "LoadCondition", condition_id, geometries[0], properties
                )
                force = -point.weight * self.modulus * np.asarray(point.normal)
                condition.SetValue(SMA.POINT_LOAD_X, float(force[0]))
                condition.SetValue(SMA.POINT_LOAD_Y, float(force[1]))
                condition.SetValue(SMA.POINT_LOAD_Z, float(force[2]))
                condition_ids.append(condition_id)
        return condition_ids


class CouplingPenalty(WeakCondition):
    """Create penalty coupling conditions from master-owned interface segments."""

    def __init__(
        self,
        *args: Any,
        slave_bounds_xyz: Bounds,
        slave_bounds_uvw: Bounds,
        slave_geometry_name: str,
        penalty_factor: float,
        slip: bool,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.slave_bounds_xyz = slave_bounds_xyz
        self.slave_bounds_uvw = slave_bounds_uvw
        self.slave_geometry_name = slave_geometry_name
        self.penalty_factor = float(penalty_factor)
        self.slip = bool(slip)

        #TODO: Add KM.CouplingGeometry to Kratos master.
        if not hasattr(KM, "CouplingGeometry"):
            raise RuntimeError(
                "CouplingPenaltyCondition requires "
                "KratosMultiphysics.CouplingGeometry, which is not yet part of "
                "Kratos master. Use https://github.com/manuelmessmer/Kratos instead."
            )

    def apply(self, master_model_part: Any, slave_model_part: Any) -> list[int]:
        """Create coupling conditions and return their root-model-part IDs."""
        properties = master_model_part.GetProperties()[self.property_id]
        properties.SetValue(IGA.PENALTY_FACTOR, self.penalty_factor)
        if hasattr(IGA, "COUPLING_SLIP"):
            properties.SetValue(IGA.COUPLING_SLIP, self.slip)
        master_volume = master_model_part.GetGeometry(self.geometry_name)
        slave_volume = slave_model_part.GetGeometry(self.slave_geometry_name)
        condition_ids: list[int] = []

        # TODO: Clip the interface against both grids and merge both clips into
        # one common interface partition before creating coupling triangles.
        for triangle in self.triangle_mesh.triangles():
            if triangle.aspect_ratio() >= 1.0e8:
                continue
            master_surface = _surface_geometry(
                master_volume,
                _triangle_geometry(triangle, self.bounds_xyz, self.bounds_uvw),
            )
            slave_surface = _surface_geometry(
                slave_volume,
                _triangle_geometry(
                    triangle,
                    self.slave_bounds_xyz,
                    self.slave_bounds_uvw,
                    reverse=True,
                ),
            )
            condition_id = master_model_part.GetRootModelPart().NumberOfConditions() + 1
            coupling_geometry = KM.CouplingGeometry(
                _quadrature_geometry(master_surface),
                _quadrature_geometry(slave_surface),
            )
            master_model_part.AddGeometry(coupling_geometry)
            condition = master_model_part.CreateNewCondition(
                "CouplingPenaltyCondition",
                condition_id,
                coupling_geometry,
                properties,
            )
            condition.Set(IGA.IgaFlags.FIX_DISPLACEMENT_X, True)
            condition.Set(IGA.IgaFlags.FIX_DISPLACEMENT_Y, True)
            condition.Set(IGA.IgaFlags.FIX_DISPLACEMENT_Z, True)

            condition_ids.append(condition_id)
        return condition_ids
