import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from queso.python_scripts.helper import *
from QuESo_PythonApplication import TriangleMesh as TriangleUtilities
import numpy as np

class WeakBcsBase():
    """Base Class to provide interface for the application of boundary conditions.

    Derived class must override 'apply()'.
    """
    def __init__(self, bcs_triangles, bounds_xyz, bounds_uvw):
        """The constructor."""
        self.bcs_triangles = bcs_triangles
        self.bounds_xyz = bounds_xyz
        self.bounds_uvw = bounds_uvw

    def IsWeakCondition():
        return True

    def apply(self, model_part):
        raise Exception("WeakBcsBase :: Function of base class is called!")

class PenaltySupport(WeakBcsBase):
    """PenaltySupport.

    Derived from WeakBcsBase.
    """
    def __init__(self, bcs_triangles, bounds_xyz, bounds_uvw, prescribed, penalty):
        """The constructor."""
        super(PenaltySupport, self).__init__(bcs_triangles, bounds_xyz, bounds_uvw)
        self.penalty = penalty
        self.prescribed = prescribed # Not used here

    def apply(self, model_part):
        """Overrides base class."""
        id_counter = model_part.NumberOfConditions() + 1
        properties = model_part.GetProperties()[1]
        properties.SetValue(IgaApplication.PENALTY_FACTOR, self.penalty)
        process_info = KM.ProcessInfo()
        nurbs_volume = model_part.GetGeometry("NurbsVolume")

        num_triangles = self.bcs_triangles.NumOfTriangles()
        for id in range(num_triangles):
            # Create kratos nodes on each vertex
            param1 = KM.Vector(3)
            param1 = PointFromGlobalToParamSpace( self.bcs_triangles.P1(id), self.bounds_xyz, self.bounds_uvw)
            node1 = KM.Node(1, param1[0], param1[1], param1[2])

            param2 = KM.Vector(3)
            param2 = PointFromGlobalToParamSpace( self.bcs_triangles.P2(id), self.bounds_xyz, self.bounds_uvw)
            node2 = KM.Node(2, param2[0], param2[1], param2[2])

            param3 = KM.Vector(3)
            param3 = PointFromGlobalToParamSpace( self.bcs_triangles.P3(id), self.bounds_xyz, self.bounds_uvw)
            node3 = KM.Node(3, param3[0], param3[1], param3[2])

            if TriangleUtilities.AspectRatioStatic(param1, param2, param3) < 1e8:
                # Create kratos triangles
                geom = KM.Triangle3D3(node1, node2, node3)
                quadrature_point_geometries = KM.GeometriesVector()

                # Create conditions
                surface_in_nurbs_volume = KM.SurfaceInNurbsVolumeGeometry(nurbs_volume, geom)
                surface_in_nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)

                surface_area = surface_in_nurbs_volume.Area()
                if( surface_area > 1e-14):
                    condition = model_part.CreateNewCondition('SupportPenaltyCondition', id_counter, quadrature_point_geometries[0], properties)
                    kratos_prescribed = KM.Vector([self.prescribed[0], self.prescribed[1], self.prescribed[2]])
                    condition.SetValue(KM.DISPLACEMENT, kratos_prescribed)
                    id_counter += 1

class LagrangeSupport(WeakBcsBase):
    """PenaltySupport.

    Derived from WeakBcsBase.
    """
    def __init__(self, bcs_triangles, bounds_xyz, bounds_uvw, prescribed):
        """The constructor."""
        super(LagrangeSupport, self).__init__(bcs_triangles, bounds_xyz, bounds_uvw)
        self.prescribed = prescribed # Not used here

    def apply(self, model_part):
        """Overrides base class."""
        id_counter = model_part.NumberOfConditions() + 1
        properties = model_part.GetProperties()[1]
        process_info = KM.ProcessInfo()
        nurbs_volume = model_part.GetGeometry("NurbsVolume")

        num_triangles = self.bcs_triangles.NumOfTriangles()
        for id in range(num_triangles):
            # Create kratos nodes on each vertex
            param1 = KM.Vector(3)
            param1 = PointFromGlobalToParamSpace( self.bcs_triangles.P1(id), self.bounds_xyz, self.bounds_uvw)
            node1 = KM.Node(1, param1[0], param1[1], param1[2])

            param2 = KM.Vector(3)
            param2 = PointFromGlobalToParamSpace( self.bcs_triangles.P2(id), self.bounds_xyz, self.bounds_uvw)
            node2 = KM.Node(2, param2[0], param2[1], param2[2])

            param3 = KM.Vector(3)
            param3 = PointFromGlobalToParamSpace( self.bcs_triangles.P3(id), self.bounds_xyz, self.bounds_uvw)
            node3 = KM.Node(3, param3[0], param3[1], param3[2])

            if TriangleUtilities.AspectRatioStatic(param1, param2, param3) < 1e8:
                # Create kratos triangles
                geom = KM.Triangle3D3(node1, node2, node3)
                quadrature_point_geometries = KM.GeometriesVector()

                # Create conditions
                surface_in_nurbs_volume = KM.SurfaceInNurbsVolumeGeometry(nurbs_volume, geom)
                surface_in_nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)

                surface_area = surface_in_nurbs_volume.Area()
                if( surface_area > 1e-14):
                    condition = model_part.CreateNewCondition('SupportLagrangeCondition', id_counter, quadrature_point_geometries[0], properties)
                    kratos_prescribed = KM.Vector([self.prescribed[0], self.prescribed[1], self.prescribed[2]])
                    condition.SetValue(KM.DISPLACEMENT, kratos_prescribed)
                    id_counter += 1

class SurfaceLoad(WeakBcsBase):
    """SurfaceLoad.

    Derived from WeakBcsBase.
    """
    def __init__(self, bcs_triangles, lower_point, upper_point, modulus, direction):
        """The constructor."""
        super(SurfaceLoad, self).__init__(bcs_triangles, lower_point, upper_point)

        norm_direction = np.sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2])
        if norm_direction < 1e-10:
            Exception("SurfaceLoad :: Norm of 'direction' is close to zero.")
        normalized_direction = direction / norm_direction

        self.force = modulus * normalized_direction
        self.conditions = []

    def apply(self, model_part ):
        """Overrides base class."""
        id_counter = model_part.NumberOfConditions() + 1
        properties = model_part.GetProperties()[1]
        nurbs_volume = model_part.GetGeometry("NurbsVolume")

        num_triangles = self.bcs_triangles.NumOfTriangles()
        for id in range(num_triangles):

            #Get points in physical space.
            points = self.bcs_triangles.GetIntegrationPointsGlobal(id, 1)
            #Create kratos condition on each point.
            for point in points:
                integration_points = []
                global_point = [point.X(), point.Y(), point.Z()]
                #Map points to local space of B-Spline box
                local_point = PointFromGlobalToParamSpace(global_point, self.bounds_xyz, self.bounds_uvw)

                integration_points.append([local_point[0], local_point[1], local_point[2], point.Weight()])
                quadrature_point_geometries_boundary = KM.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries_boundary, 2, integration_points)

                weight = point.Weight() # Weight contains all mapping terms.
                if weight > 1e-14:
                    condition = model_part.CreateNewCondition("LoadCondition", id_counter, quadrature_point_geometries_boundary[0], properties)

                    force_x = weight * self.force[0]
                    force_y = weight * self.force[1]
                    force_z = weight * self.force[2]
                    condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, force_x)
                    condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, force_y)
                    condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, force_z)
                    self.conditions.append(condition)
                    id_counter += 1

class PressureLoad(WeakBcsBase):
    """PressureLoad.

    Derived from WeakBcsBase.
    """
    def __init__(self, bcs_triangles, lower_point, upper_point, modulus):
        """The constructor."""
        super(PressureLoad, self).__init__(bcs_triangles, lower_point, upper_point)
        self.modulus = modulus
        self.conditions = []

    def apply(self, model_part ):
        """Overrides from base class."""
        id_counter = model_part.NumberOfConditions() + 1
        properties = model_part.GetProperties()[1]
        nurbs_volume = model_part.GetGeometry("NurbsVolume")

        num_triangles = self.bcs_triangles.NumOfTriangles()
        for id in range(num_triangles):

            #Get points in physical space.
            points = self.bcs_triangles.GetIntegrationPointsGlobal(id, 1)
            #Create kratos condition on each point.
            for point in points:
                integration_points = []
                global_point = [point.X(), point.Y(), point.Z()]
                #Map points to local space of B-Spline box
                local_point = PointFromGlobalToParamSpace(global_point, self.bounds_xyz, self.bounds_uvw)

                integration_points.append([local_point[0], local_point[1], local_point[2], point.Weight()])
                quadrature_point_geometries_boundary = KM.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries_boundary, 2, integration_points)

                weight = point.Weight() # Weight contains all mapping terms.
                if weight > 1e-14:
                    condition = model_part.CreateNewCondition("LoadCondition", id_counter, quadrature_point_geometries_boundary[0], properties)

                    normal = point.Normal()
                    force_x = -1.0*normal[0] * weight * self.modulus
                    force_y = -1.0*normal[1] * weight * self.modulus
                    force_z = -1.0*normal[2] * weight * self.modulus

                    condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, force_x)
                    condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, force_y)
                    condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, force_z)
                    self.conditions.append(condition)
                    id_counter += 1