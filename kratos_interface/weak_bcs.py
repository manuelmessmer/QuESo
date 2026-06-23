import numpy as np
from typing import List, Tuple
# Import QuESo
import QuESoPythonModule as QuESo
from QuESoPythonModule.scripts.helper import *
# Import Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Type definition
Point3D = Tuple[float, float, float]
LoadVector = Tuple[float, float, float]

class WeakBcsBase():
    """Base class for applying weak boundary conditions.

    Derived classes must override `apply()`.
    """
    def __init__(self,
            bcs_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            bounds_xyz: Tuple[Point3D, Point3D],
            bounds_uvw: Tuple[Point3D, Point3D],
            property_id: int,
            nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        """Initializes WeakBcsBase.

        Args:
            bcs_triangles (QuESo.TriangleMesh): Triangle mesh of boundary.
            bounds_xyz (Tuple[Point3D, Point3D]): Global coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            bounds_uvw (Tuple[Point3D, Point3D]): Parametric coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
        """
        self.bcs_triangles = bcs_triangles
        self.bounds_xyz = bounds_xyz
        self.bounds_uvw = bounds_uvw
        self.property_id = property_id
        self.nurbs_volume_name = nurbs_volume_name

    def apply(self, model_part: KM.ModelPart) -> list[tuple[int, LoadVector]]:
        """Applies the boundary condition (to be overridden).

        Args:
            model_part (KM.ModelPart): Kratos model part.

        Raises:
            Exception: Always raised. This function must be overridden.
        """
        raise Exception("WeakBcsBase :: Function of base class is called!")

    def _collect_integration_points(self) -> list[tuple[Point3D, float, np.ndarray]]:
        integration_points: list[tuple[Point3D, float, np.ndarray]] = []
        triangle_meshes = self.bcs_triangles
        if not isinstance(triangle_meshes, list):
            triangle_meshes = [triangle_meshes]

        for triangle_mesh in triangle_meshes:
            for triangle in triangle_mesh.Triangles():
                for point in triangle.GetIPsGlobal(1):
                    weight = point.weight
                    if weight < 1e-14:
                        continue

                    global_point = (point.x, point.y, point.z)
                    local_point = point_from_global_to_param_space(
                        global_point,
                        self.bounds_xyz,
                        self.bounds_uvw,
                    )
                    integration_points.append(
                        (local_point, weight, np.asarray(point.normal, dtype=float))
                    )
        return integration_points

    @staticmethod
    def is_weak_condition() -> bool:
        """Checks if the condition is weak.

        Returns:
            bool: Always True.
        """
        return True

class PenaltySupport(WeakBcsBase):
    """PenaltySupport.

    Derived from WeakBcsBase.
    """
    def __init__(self,
            bcs_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            bounds_xyz: Tuple[Point3D, Point3D],
            bounds_uvw: Tuple[Point3D, Point3D],
            prescribed: Tuple[float, float, float],
            penalty: float,
            property_id: int,
            nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        """Initializes PenaltySupport.

        Args:
            bcs_triangles (QuESo.TriangleMesh): Triangle mesh for BC.
            bounds_xyz (Tuple[Point3D, Point3D]): Global coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            bounds_uvw (Tuple[Point3D, Point3D]): Parametric coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            prescribed (Tuple[float, float, float]): Prescribed displacement vector.
            penalty (float): Penalty factor.
        """
        super(PenaltySupport, self).__init__(
            bcs_triangles,
            bounds_xyz,
            bounds_uvw,
            property_id,
            nurbs_volume_name,
        )
        self.prescribed = prescribed
        self.penalty = penalty

    def apply(self, model_part: KM.ModelPart) -> list[tuple[int, LoadVector]]:
        """Applies the penalty support condition.

        Args:
            model_part (KM.ModelPart): Kratos model part.
        """
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1
        properties = model_part.GetProperties()[self.property_id]
        properties.SetValue(IgaApplication.PENALTY_FACTOR, self.penalty) # type: ignore
        nurbs_volume = model_part.GetGeometry(self.nurbs_volume_name)
        kratos_prescribed = KM.Vector([self.prescribed[0], self.prescribed[1], self.prescribed[2]])

        for tri in self.bcs_triangles.Triangles():
            # Map triangle points to parametric space
            params = [
                point_from_global_to_param_space(tri.p1, self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(tri.p2, self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(tri.p3, self.bounds_xyz, self.bounds_uvw),
            ]

            if tri.AspectRatio() >= 1e8:
                continue  # Skip badly shaped triangles

            # Create triangle geometry
            nodes = [KM.Node(i + 1, *param) for i, param in enumerate(params)]
            geom = KM.Triangle3D3(*nodes)
            quadrature_point_geometries = KM.GeometriesVector()

            # Create conditions
            surface_in_nurbs_volume = KM.SurfaceInNurbsVolumeGeometry(nurbs_volume, geom)
            surface_in_nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)

            if surface_in_nurbs_volume.Area() <= 1e-14:
                continue  # Skip negligible surfaces

            condition = model_part.CreateNewCondition(
                'SupportPenaltyCondition',
                id_counter,
                quadrature_point_geometries[0],
                properties
            )
            condition.SetValue(KM.DISPLACEMENT, kratos_prescribed)
            id_counter += 1
        return []

class LagrangeSupport(WeakBcsBase):
    """LagrangeMultiplierSupport.

    Derived from WeakBcsBase.
    """
    def __init__(self,
            bcs_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            bounds_xyz: Tuple[Point3D, Point3D],
            bounds_uvw: Tuple[Point3D, Point3D],
            prescribed: Tuple[float, float, float],
            property_id: int,
            nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        """Initializes LagrangeMultiplierSupport.

        Args:
            bcs_triangles (QuESo.TriangleMes): Triangle mesh.
            bounds_xyz (Tuple[Point3D, Point3D]): Global coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            bounds_uvw (Tuple[Point3D, Point3D]): Parametric coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            prescribed (Tuple[float, float, float]): Prescribed displacement.
        """
        super(LagrangeSupport, self).__init__(
            bcs_triangles,
            bounds_xyz,
            bounds_uvw,
            property_id,
            nurbs_volume_name,
        )
        self.prescribed = prescribed

    def apply(self, model_part: KM.ModelPart) -> list[tuple[int, LoadVector]]:
        """Applies the Lagrange support condition.

        Args:
            model_part (KM.ModelPart): Kratos model part.
        """
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1
        properties = model_part.GetProperties()[self.property_id]
        nurbs_volume = model_part.GetGeometry(self.nurbs_volume_name)
        kratos_prescribed = KM.Vector(self.prescribed)

        # Iterate over all triangles
        for tri in self.bcs_triangles.Triangles():
            # Map triangle vertices to parametric space
            params = [
                point_from_global_to_param_space(tri.p1, self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(tri.p2, self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(tri.p3, self.bounds_xyz, self.bounds_uvw),
            ]

            # Skip bad quality triangles
            if tri.AspectRatio() >= 1e8:
                continue

            # Create triangle geometry
            nodes = [KM.Node(i + 1, *param) for i, param in enumerate(params)]
            geom = KM.Triangle3D3(*nodes)
            quadrature_point_geometries = KM.GeometriesVector()

            # Generate quadrature points and area
            surface_in_volume = KM.SurfaceInNurbsVolumeGeometry(nurbs_volume, geom)
            surface_in_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)

            if surface_in_volume.Area() <= 1e-14:
                continue

            # Create condition
            condition = model_part.CreateNewCondition(
                'SupportLagrangeCondition',
                id_counter,
                quadrature_point_geometries[0],
                properties
            )
            condition.SetValue(KM.DISPLACEMENT, kratos_prescribed)
            id_counter += 1
        return []

class SurfaceLoad(WeakBcsBase):
    """SurfaceLoad.

    Derived from WeakBcsBase.
    """
    def __init__(self,
            bcs_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            bounds_xyz: Tuple[Point3D, Point3D],
            bounds_uvw: Tuple[Point3D, Point3D],
            modulus: float,
            direction: Tuple[float, float, float],
            property_id: int,
            nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        """Initializes SurfaceLoad.

        Args:
            bcs_triangles (QuESo.TriangleMesh): Triangle mesh.
            bounds_xyz (Tuple[Point3D, Point3D]): Global coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            bounds_uvw (Tuple[Point3D, Point3D]): Parametric coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            modulus (float): Load magnitude.
            direction (Tuple[float, float, float]): Load direction vector.
        """
        super(SurfaceLoad, self).__init__(
            bcs_triangles,
            bounds_xyz,
            bounds_uvw,
            property_id,
            nurbs_volume_name,
        )

        direction_array = np.array(direction)
        norm_direction = np.linalg.norm(direction_array)
        if norm_direction < 1e-10:
            self.force_density = np.zeros(3)
        else:
            normalized_direction = direction_array / norm_direction
            self.force_density = modulus * normalized_direction

    def apply(self, model_part: KM.ModelPart) -> list[tuple[int, LoadVector]]:
        """Applies the surface load.

        Args:
            model_part (KM.ModelPart): Kratos model part.
        """
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1
        properties = model_part.GetProperties()[self.property_id]
        nurbs_volume = model_part.GetGeometry(self.nurbs_volume_name)
        created_loads: list[tuple[int, LoadVector]] = []

        for local_point, weight, _ in self._collect_integration_points():
            integration_points = [[local_point[0], local_point[1], local_point[2], weight]]
            quadrature_point_geometries_boundary = KM.GeometriesVector()
            nurbs_volume.CreateQuadraturePointGeometries(
                quadrature_point_geometries_boundary,
                2,
                integration_points,
            )

            force = weight * self.force_density
            condition_id = id_counter
            condition = model_part.CreateNewCondition(
                "LoadCondition",
                condition_id,
                quadrature_point_geometries_boundary[0],
                properties
            )
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, force[0]) # type: ignore
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, force[1]) # type: ignore
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, force[2]) # type: ignore
            created_loads.append((condition_id, (float(force[0]), float(force[1]), float(force[2]))))
            id_counter += 1

        return created_loads


class TotalLoad(WeakBcsBase):
    """Total resultant force distributed over the active surface."""

    def __init__(self,
            bcs_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            bounds_xyz: Tuple[Point3D, Point3D],
            bounds_uvw: Tuple[Point3D, Point3D],
            modulus: float,
            direction: Tuple[float, float, float],
            property_id: int,
            nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        super(TotalLoad, self).__init__(
            bcs_triangles,
            bounds_xyz,
            bounds_uvw,
            property_id,
            nurbs_volume_name,
        )

        direction_array = np.array(direction)
        norm_direction = np.linalg.norm(direction_array)
        if norm_direction < 1e-10:
            self.total_force = np.zeros(3)
        else:
            normalized_direction = direction_array / norm_direction
            self.total_force = modulus * normalized_direction

    def apply(self, model_part: KM.ModelPart) -> list[tuple[int, LoadVector]]:
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1
        properties = model_part.GetProperties()[self.property_id]
        nurbs_volume = model_part.GetGeometry(self.nurbs_volume_name)
        created_loads: list[tuple[int, LoadVector]] = []
        collected_points = self._collect_integration_points()
        total_area = sum(weight for _, weight, _ in collected_points)

        if total_area <= 1e-14:
            raise RuntimeError("TotalLoadCondition requires a non-zero active surface area.")

        for local_point, weight, _ in collected_points:
            integration_points = [[local_point[0], local_point[1], local_point[2], weight]]
            quadrature_point_geometries_boundary = KM.GeometriesVector()
            nurbs_volume.CreateQuadraturePointGeometries(
                quadrature_point_geometries_boundary,
                2,
                integration_points,
            )

            force = (weight / total_area) * self.total_force
            condition_id = id_counter
            condition = model_part.CreateNewCondition(
                "LoadCondition",
                condition_id,
                quadrature_point_geometries_boundary[0],
                properties
            )
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, force[0]) # type: ignore
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, force[1]) # type: ignore
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, force[2]) # type: ignore
            created_loads.append((condition_id, (float(force[0]), float(force[1]), float(force[2]))))
            id_counter += 1

        return created_loads

class PressureLoad(WeakBcsBase):
    """PressureLoad.

    Derived from WeakBcsBase.
    """
    def __init__(self,
            bcs_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            bounds_xyz: Tuple[Point3D, Point3D],
            bounds_uvw: Tuple[Point3D, Point3D],
            modulus: float,
            property_id: int,
            nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        """Initializes PressureLoad.

        Args:
            bcs_triangles (QuESo.TriangleMesh): Triangle mesh.
            bounds_xyz (Tuple[Point3D, Point3D]): Global coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            bounds_uvw (Tuple[Point3D, Point3D]): Parametric coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            modulus (float): Pressure value (magnitude).
        """
        super(PressureLoad, self).__init__(
            bcs_triangles,
            bounds_xyz,
            bounds_uvw,
            property_id,
            nurbs_volume_name,
        )
        self.modulus = modulus

    def apply(self, model_part: KM.ModelPart) -> list[tuple[int, LoadVector]]:
        """Applies the PressureLoad.

        Args:
            model_part (KM.ModelPart): Kratos model part.
        """
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1
        properties = model_part.GetProperties()[self.property_id]
        nurbs_volume = model_part.GetGeometry(self.nurbs_volume_name)
        created_loads: list[tuple[int, LoadVector]] = []

        for local_point, weight, normal in self._collect_integration_points():
            integration_points = [[local_point[0], local_point[1], local_point[2], weight]]
            quadrature_point_geometries_boundary = KM.GeometriesVector()
            nurbs_volume.CreateQuadraturePointGeometries(
                quadrature_point_geometries_boundary,
                2,
                integration_points,
            )

            force = -1.0 * weight * self.modulus * normal
            condition_id = id_counter
            condition = model_part.CreateNewCondition(
                "LoadCondition",
                condition_id,
                quadrature_point_geometries_boundary[0],
                properties
            )

            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, force[0]) # type: ignore
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, force[1]) # type: ignore
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, force[2]) # type: ignore
            created_loads.append((condition_id, (float(force[0]), float(force[1]), float(force[2]))))
            id_counter += 1

        return created_loads


class CouplingPenalty(WeakBcsBase):
    def __init__(self,
            master_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            slave_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            master_bounds_xyz: Tuple[Point3D, Point3D],
            master_bounds_uvw: Tuple[Point3D, Point3D],
            slave_bounds_xyz: Tuple[Point3D, Point3D],
            slave_bounds_uvw: Tuple[Point3D, Point3D],
            penalty_factor: float,
            property_id: int,
            slip: bool = False,
            master_nurbs_volume_name: str = "NurbsVolume",
            slave_nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        super(CouplingPenalty, self).__init__(
            master_triangles,
            master_bounds_xyz,
            master_bounds_uvw,
            property_id,
            master_nurbs_volume_name,
        )
        self.master_triangles = master_triangles
        self.slave_triangles = slave_triangles
        self.slave_bounds_xyz = slave_bounds_xyz
        self.slave_bounds_uvw = slave_bounds_uvw
        self.penalty_factor = float(penalty_factor)
        self.slip = bool(slip)
        self.slave_nurbs_volume_name = slave_nurbs_volume_name
        self._cache = []

    @staticmethod
    def is_weak_condition() -> bool:
        return True


    def apply(self, model_part: KM.ModelPart) -> list[tuple[int, LoadVector]]:
        if self.master_triangles.NumOfTriangles() != self.slave_triangles.NumOfTriangles():
            raise RuntimeError(
                "CouplingPenalty :: Master/slave triangle counts differ: "
                f"{self.master_triangles.NumOfTriangles()} != {self.slave_triangles.NumOfTriangles()}"
            )
        properties = model_part.GetProperties()[self.property_id]
        properties.SetValue(IgaApplication.PENALTY_FACTOR, self.penalty_factor) # type: ignore
        properties.SetValue(IgaApplication.COUPLING_SLIP, self.slip) # type: ignore

        master_volume = model_part.GetGeometry(self.nurbs_volume_name)
        slave_volume = model_part.GetGeometry(self.slave_nurbs_volume_name)
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1

        for master_triangle, slave_triangle in zip(
            self.master_triangles.Triangles(),
            self.slave_triangles.Triangles(),
        ):
            master_params = [
                point_from_global_to_param_space(master_triangle.p1, self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(master_triangle.p2, self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(master_triangle.p3, self.bounds_xyz, self.bounds_uvw),
            ]
            slave_params = [
                point_from_global_to_param_space(slave_triangle.p1, self.slave_bounds_xyz, self.slave_bounds_uvw),
                point_from_global_to_param_space(slave_triangle.p2, self.slave_bounds_xyz, self.slave_bounds_uvw),
                point_from_global_to_param_space(slave_triangle.p3, self.slave_bounds_xyz, self.slave_bounds_uvw),
            ]

            if master_triangle.AspectRatio() >= 1e8: 
                continue
            if slave_triangle.AspectRatio() >= 1e8: 
                continue

            master_nodes = [KM.Node(i + 1, *param) for i, param in enumerate(master_params)]
            slave_nodes = [KM.Node(i + 1, *param) for i, param in enumerate(slave_params)]
            master_geom = KM.Triangle3D3(*master_nodes)
            slave_geom = KM.Triangle3D3(*slave_nodes)

            master_surface = KM.SurfaceInNurbsVolumeGeometry(master_volume, master_geom)
            slave_surface = KM.SurfaceInNurbsVolumeGeometry(slave_volume, slave_geom)

            master_qpg = KM.GeometriesVector()
            slave_qpg = KM.GeometriesVector()
            master_surface.CreateQuadraturePointGeometries(master_qpg, 2)
            slave_surface.CreateQuadraturePointGeometries(slave_qpg, 2)

            model_part.CreateNewCouplingCondition(
                id_counter,
                master_qpg[0],
                slave_qpg[0],
                properties,
            )
            id_counter += 1

        return []
