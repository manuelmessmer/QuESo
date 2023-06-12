import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from tibra.python_scripts.helper import *

class WeakBcsBase():
    """Base Class to provide interface for the application of boundary conditions.

    Derived class must override 'apply()'.
    """
    def __init__(self, bcs_triangles, lower_point, upper_point):
        """The constructor."""
        self.bcs_triangles = bcs_triangles
        self.lower_point = lower_point
        self.upper_point = upper_point

    def apply(self, model_part):
        raise Exception("WeakBcsBase :: Function of base class is called!")

class PenaltySupport(WeakBcsBase):
    """PenaltySupport.

    Derived from WeakBcsBase.
    """
    def __init__(self, bcs_triangles, lower_point, upper_point, prescribed, penalty):
        """The constructor."""
        super(PenaltySupport, self).__init__(bcs_triangles, lower_point, upper_point)
        self.penalty = penalty
        self.prescribed = prescribed # Not used here

    def apply(self, model_part):
        """Overrides base class."""
        id_counter = model_part.NumberOfConditions() + 1
        properties = model_part.GetProperties()[0]
        properties.SetValue(IgaApplication.PENALTY_FACTOR, self.penalty)
        process_info = KM.ProcessInfo()
        nurbs_volume = model_part.GetGeometry("NurbsVolume")

        num_triangles = self.bcs_triangles.NumOfTriangles()
        for id in range(num_triangles):
            # Create kratos nodes on each vertex
            param1 = KM.Vector(3)
            param1 = PointFromGlobalToParamSpace( self.bcs_triangles.P1(id), self.lower_point, self.upper_point)
            node1 = KM.Node(1, param1[0], param1[1], param1[2])

            param2 = KM.Vector(3)
            param2 = PointFromGlobalToParamSpace( self.bcs_triangles.P2(id), self.lower_point, self.upper_point)
            node2 = KM.Node(2, param2[0], param2[1], param2[2])

            param3 = KM.Vector(3)
            param3 = PointFromGlobalToParamSpace( self.bcs_triangles.P3(id), self.lower_point, self.upper_point)
            node3 = KM.Node(3, param3[0], param3[1], param3[2])

            # Create kratos triangles
            geom = KM.Triangle3D3(node1, node2, node3)
            quadrature_point_geometries = KM.GeometriesVector()

            # Create conditions
            surface_in_nurbs_volume = KM.SurfaceInNurbsVolumeGeometry(nurbs_volume, geom)
            surface_in_nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)

            cond = model_part.CreateNewCondition('SupportPenaltyCondition', id_counter, quadrature_point_geometries[0], properties)
            displacement = [KM.Vector([0.0, 0.0, 0.0])]
            cond.SetValuesOnIntegrationPoints(KM.DISPLACEMENT,displacement,process_info)
            id_counter += 1

class SurfaceLoad(WeakBcsBase):
    """SurfaceLoad.

    Derived from WeakBcsBase.
    """
    def __init__(self, bcs_triangles, lower_point, upper_point, force, normal_flag):
        """The constructor."""
        super(SurfaceLoad, self).__init__(bcs_triangles, lower_point, upper_point)
        self.force = force
        self.normal_flag = normal_flag
        self.conditions = []

    def apply(self, model_part ):
        """Overrides base class."""
        id_counter = model_part.NumberOfConditions() + 1
        properties = model_part.GetProperties()[0]
        nurbs_volume = model_part.GetGeometry("NurbsVolume")

        num_triangles = self.bcs_triangles.NumOfTriangles()
        for id in range(num_triangles):

            #Get points in physical space.
            points = self.bcs_triangles.GetIntegrationPointsGlobal(id, 1)
            #Create kratos condition on each point.
            for point in points:
                integration_points = []
                global_point = [point.GetX(), point.GetY(), point.GetZ()]
                #Map points to local space of B-Spline box
                local_point = PointFromGlobalToParamSpace(global_point, self.lower_point, self.upper_point)

                integration_points.append([local_point[0], local_point[1], local_point[2], point.GetWeight()])
                quadrature_point_geometries_boundary = KM.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries_boundary, 2, integration_points)

                condition = model_part.CreateNewCondition("LoadCondition", id_counter, quadrature_point_geometries_boundary[0], properties)
                weight = point.GetWeight() # Weight contains all mapping terms.
                if self.normal_flag:
                    normal = point.Normal()
                    force_x = normal[0] * weight * self.force
                    force_y = normal[1] * weight * self.force
                    force_z = normal[2] * weight * self.force
                else:
                    force_x = weight * self.force[0]
                    force_y = weight * self.force[1]
                    force_z = weight * self.force[2]

                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, force_x)
                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, force_y)
                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, force_z)
                self.conditions.append(condition)
                id_counter += 1



