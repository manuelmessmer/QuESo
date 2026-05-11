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

class WeakBcsBase():
    """Base class for applying weak boundary conditions.

    Derived classes must override `apply()`.
    """
    def __init__(self,
            bcs_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            bounds_xyz: Tuple[Point3D, Point3D],
            bounds_uvw: Tuple[Point3D, Point3D],
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
        self.nurbs_volume_name = nurbs_volume_name

    def apply(self, model_part: KM.ModelPart) -> None:
        """Applies the boundary condition (to be overridden).

        Args:
            model_part (KM.ModelPart): Kratos model part.

        Raises:
            Exception: Always raised. This function must be overridden.
        """
        raise Exception("WeakBcsBase :: Function of base class is called!")

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
        super(PenaltySupport, self).__init__(bcs_triangles, bounds_xyz, bounds_uvw, nurbs_volume_name)
        self.prescribed = prescribed
        self.penalty = penalty

    def apply(self, model_part: KM.ModelPart) -> None:
        """Applies the penalty support condition.

        Args:
            model_part (KM.ModelPart): Kratos model part.
        """
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1
        properties = model_part.GetProperties()[1]
        properties.SetValue(IgaApplication.PENALTY_FACTOR, self.penalty) # type: ignore
        nurbs_volume = model_part.GetGeometry(self.nurbs_volume_name)
        kratos_prescribed = KM.Vector([self.prescribed[0], self.prescribed[1], self.prescribed[2]])

        for triangle_id in range(self.bcs_triangles.NumOfTriangles()):
            # Map triangle points to parametric space
            params = [
                point_from_global_to_param_space(self.bcs_triangles.P1(triangle_id), self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(self.bcs_triangles.P2(triangle_id), self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(self.bcs_triangles.P3(triangle_id), self.bounds_xyz, self.bounds_uvw),
            ]

            if QuESo.TriangleMesh.AspectRatioStatic(*params) >= 1e8: # type: ignore (TODO: add .pyi)
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

class LagrangeSupport(WeakBcsBase):
    """LagrangeMultiplierSupport.

    Derived from WeakBcsBase.
    """
    def __init__(self,
            bcs_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            bounds_xyz: Tuple[Point3D, Point3D],
            bounds_uvw: Tuple[Point3D, Point3D],
            prescribed: Tuple[float, float, float],
            nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        """Initializes LagrangeMultiplierSupport.

        Args:
            bcs_triangles (QuESo.TriangleMes): Triangle mesh.
            bounds_xyz (Tuple[Point3D, Point3D]): Global coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            bounds_uvw (Tuple[Point3D, Point3D]): Parametric coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            prescribed (Tuple[float, float, float]): Prescribed displacement.
        """
        super(LagrangeSupport, self).__init__(bcs_triangles, bounds_xyz, bounds_uvw, nurbs_volume_name)
        self.prescribed = prescribed

    def apply(self, model_part: KM.ModelPart) -> None:
        """Applies the Lagrange support condition.

        Args:
            model_part (KM.ModelPart): Kratos model part.
        """
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1
        properties = model_part.GetProperties()[1]
        nurbs_volume = model_part.GetGeometry(self.nurbs_volume_name)
        kratos_prescribed = KM.Vector(self.prescribed)

        # Iterate over all triangles
        for triangle_id in range(self.bcs_triangles.NumOfTriangles()):
            # Map triangle vertices to parametric space
            params = [
                point_from_global_to_param_space(self.bcs_triangles.P1(triangle_id), self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(self.bcs_triangles.P2(triangle_id), self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(self.bcs_triangles.P3(triangle_id), self.bounds_xyz, self.bounds_uvw),
            ]

            # Skip bad quality triangles
            if QuESo.TriangleMesh.AspectRatioStatic(*params) >= 1e8: # type: ignore (TODO: add .pyi)
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
        super(SurfaceLoad, self).__init__(bcs_triangles, bounds_xyz, bounds_uvw, nurbs_volume_name)

        direction_array = np.array(direction)
        norm_direction = np.linalg.norm(direction_array)
        if norm_direction < 1e-10:
            raise RuntimeError("SurfaceLoad :: Norm of 'direction' is close to zero.")
        normalized_direction = direction_array / norm_direction
        self.force = modulus * normalized_direction

    def apply(self, model_part: KM.ModelPart) -> None:
        """Applies the surface load.

        Args:
            model_part (KM.ModelPart): Kratos model part.
        """
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1
        properties = model_part.GetProperties()[1]
        nurbs_volume = model_part.GetGeometry(self.nurbs_volume_name)

        for triangle_id in range(self.bcs_triangles.NumOfTriangles()):
            #Get points in physical space.
            points = self.bcs_triangles.GetIntegrationPointsGlobal(triangle_id, 1)

            #Create kratos condition on each point.
            for point in points:
                global_point = (point.X(), point.Y(), point.Z())
                #Map points to local space of B-Spline box
                local_point = point_from_global_to_param_space(global_point, self.bounds_xyz, self.bounds_uvw)

                # Create quadrature points
                integration_points = [[local_point[0], local_point[1], local_point[2], point.Weight()]]
                quadrature_point_geometries_boundary = KM.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries_boundary, 2, integration_points)

                weight = point.Weight() # Weight contains all mapping terms.
                if weight < 1e-14:
                    continue # Skip insignificant weights

                condition = model_part.CreateNewCondition(
                    "LoadCondition",
                    id_counter,
                    quadrature_point_geometries_boundary[0],
                    properties
                )
                force = weight * np.array(self.force)
                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, force[0]) # type: ignore
                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, force[1]) # type: ignore
                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, force[2]) # type: ignore
                id_counter += 1
                geo = condition.GetGeometry()

class PressureLoad(WeakBcsBase):
    """PressureLoad.

    Derived from WeakBcsBase.
    """
    def __init__(self,
            bcs_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            bounds_xyz: Tuple[Point3D, Point3D],
            bounds_uvw: Tuple[Point3D, Point3D],
            modulus: float,
            nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        """Initializes PressureLoad.

        Args:
            bcs_triangles (QuESo.TriangleMesh): Triangle mesh.
            bounds_xyz (Tuple[Point3D, Point3D]): Global coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            bounds_uvw (Tuple[Point3D, Point3D]): Parametric coordinate bounds ([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            modulus (float): Pressure value (magnitude).
        """
        super(PressureLoad, self).__init__(bcs_triangles, bounds_xyz, bounds_uvw, nurbs_volume_name)
        self.modulus = modulus

    def apply(self, model_part: KM.ModelPart) -> None:
        """Applies the PressureLoad.

        Args:
            model_part (KM.ModelPart): Kratos model part.
        """
        id_counter = model_part.GetRootModelPart().NumberOfConditions() + 1
        properties = model_part.GetProperties()[1]
        nurbs_volume = model_part.GetGeometry(self.nurbs_volume_name)

        for triangle_id in range(self.bcs_triangles.NumOfTriangles()):
            #Get points in physical space.
            points = self.bcs_triangles.GetIntegrationPointsGlobal(triangle_id, 1)

            #Create kratos condition on each point.
            for point in points:
                global_point = (point.X(), point.Y(), point.Z())

                # Map points to local space of B-Spline box
                local_point = point_from_global_to_param_space(global_point, self.bounds_xyz, self.bounds_uvw)

                # Create quadrature points
                integration_points = [[local_point[0], local_point[1], local_point[2], point.Weight()]]
                quadrature_point_geometries_boundary = KM.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries_boundary, 2, integration_points)

                weight = point.Weight() # Weight contains all mapping terms.
                if weight < 1e-14:
                    continue  # Skip insignificant weights

                condition = model_part.CreateNewCondition(
                    "LoadCondition",
                    id_counter,
                    quadrature_point_geometries_boundary[0],
                    properties
                )

                normal = point.Normal()
                force = -1.0 * weight * self.modulus * np.array(normal)

                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, force[0]) # type: ignore
                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, force[1]) # type: ignore
                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, force[2]) # type: ignore
                id_counter += 1


class CouplingPenalty(WeakBcsBase):
    def __init__(self,
            master_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            slave_triangles: QuESo.TriangleMesh, # type: ignore (TODO: add .pyi)
            master_bounds_xyz: Tuple[Point3D, Point3D],
            master_bounds_uvw: Tuple[Point3D, Point3D],
            slave_bounds_xyz: Tuple[Point3D, Point3D],
            slave_bounds_uvw: Tuple[Point3D, Point3D],
            penalty_factor: float,
            slip: bool = False,
            property_id: int = 1000,
            master_nurbs_volume_name: str = "NurbsVolume",
            slave_nurbs_volume_name: str = "NurbsVolume"
        ) -> None:
        super(CouplingPenalty, self).__init__(master_triangles, master_bounds_xyz, master_bounds_uvw, master_nurbs_volume_name)
        self.master_triangles = master_triangles
        self.slave_triangles = slave_triangles
        self.slave_bounds_xyz = slave_bounds_xyz
        self.slave_bounds_uvw = slave_bounds_uvw
        self.penalty_factor = float(penalty_factor)
        self.slip = bool(slip)
        self.property_id = int(property_id)
        self.slave_nurbs_volume_name = slave_nurbs_volume_name
        self._cache = []

    @staticmethod
    def is_weak_condition() -> bool:
        return False

    def apply(self, model_part: KM.ModelPart) -> None:
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
        id_counter_0 = id_counter

        created_conditions = 0
        skipped_zero_area = 0
        skipped_qp_mismatch = 0
        skipped_empty_parts = 0

        for triangle_id in range(self.master_triangles.NumOfTriangles()):
            master_params = [
                point_from_global_to_param_space(self.master_triangles.P1(triangle_id), self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(self.master_triangles.P2(triangle_id), self.bounds_xyz, self.bounds_uvw),
                point_from_global_to_param_space(self.master_triangles.P3(triangle_id), self.bounds_xyz, self.bounds_uvw),
            ]
            slave_params = [
                point_from_global_to_param_space(self.slave_triangles.P1(triangle_id), self.slave_bounds_xyz, self.slave_bounds_uvw),
                point_from_global_to_param_space(self.slave_triangles.P2(triangle_id), self.slave_bounds_xyz, self.slave_bounds_uvw),
                point_from_global_to_param_space(self.slave_triangles.P3(triangle_id), self.slave_bounds_xyz, self.slave_bounds_uvw),
            ]

            if QuESo.TriangleMesh.AspectRatioStatic(*master_params) >= 1e8: # type: ignore
                continue
            if QuESo.TriangleMesh.AspectRatioStatic(*slave_params) >= 1e8: # type: ignore
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
            created_conditions += 1

        # KM.Logger.PrintInfo(
        #     "CouplingPenalty",
        #     (
        #         f"created={created_conditions}, skipped_zero_area={skipped_zero_area}, "
        #         f"skipped_qp_mismatch={skipped_qp_mismatch}, skipped_empty_parts={skipped_empty_parts}, "
        #         "skipped_empty_coupling_geom=0"
        #     )
        # )
        if created_conditions == 0:
            raise RuntimeError(
                "CouplingPenalty :: No coupling conditions were created. "
                f"Details: triangles={self.master_triangles.NumOfTriangles()}, "
                f"created={created_conditions}, skipped_zero_area={skipped_zero_area}, "
                f"skipped_qp_mismatch={skipped_qp_mismatch}, skipped_empty_parts={skipped_empty_parts}, "
                "skipped_empty_coupling_geom=0."
            )
