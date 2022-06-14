import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from src.python_scripts.helper import *
import numpy as np

class WeakBcsBase():
    def __init__(self, bcs_triangles, lower_point, upper_point):
        self.bcs_triangles = bcs_triangles
        self.lower_point = lower_point
        self.upper_point = upper_point

    def apply(self, model_part):
        raise Exception("WeakBcsBase :: Function of base class is called!")

class PenaltySupport(WeakBcsBase):
    def __init__(self, bcs_triangles, lower_point, upper_point, penalty):
        super(PenaltySupport, self).__init__(bcs_triangles, lower_point, upper_point)
        self.penalty = penalty

    def apply(self, model_part):
        id_counter = model_part.NumberOfConditions() + 1
        properties = model_part.GetProperties()[0]
        properties.SetValue(IgaApplication.PENALTY_FACTOR,self.penalty)
        process_info = KM.ProcessInfo()
        nurbs_volume = model_part.GetGeometry("NurbsVolume")

        for triangle in self.bcs_triangles:
            param1 = KM.Vector(3)
            param1 = FromGlobalToParamSpace( triangle.P1(), self.lower_point, self.upper_point)
            node1 = KM.Node(1, param1[0], param1[1], param1[2])

            param2 = KM.Vector(3)
            param2 = FromGlobalToParamSpace( triangle.P2(), self.lower_point, self.upper_point)
            node2 = KM.Node(2, param2[0], param2[1], param2[2])

            param3 = KM.Vector(3)
            param3 = FromGlobalToParamSpace( triangle.P3(), self.lower_point, self.upper_point)
            node3 = KM.Node(3, param3[0], param3[1], param3[2])

            # Create kratos triangles
            geom = KM.Triangle3D3(node1, node2, node3)
            quadrature_point_geometries = KM.GeometriesVector()

            surface_in_nurbs_volume = KM.SurfaceInNurbsVolumeGeometry(nurbs_volume, geom)
            surface_in_nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)

            cond = model_part.CreateNewCondition('SupportPenaltyCondition', id_counter, quadrature_point_geometries[0], properties)
            displacement = [KM.Vector([0.0,0.0,0])]
            cond.SetValuesOnIntegrationPoints(KM.DISPLACEMENT,displacement,process_info)
            id_counter += 1

    def change(self, model_part, time ):
        return 0

class SurfaceLoad(WeakBcsBase):
    def __init__(self, bcs_triangles, lower_point, upper_point, force, normal_flag):
        super(SurfaceLoad, self).__init__(bcs_triangles, lower_point, upper_point)
        self.force = force
        self.normal_flag = normal_flag
        self.p0 = []
        self.conditions = []

    def apply(self, model_part ):
        id_counter = model_part.NumberOfConditions() + 1
        properties = model_part.GetProperties()[0]
        process_info = KM.ProcessInfo()
        nurbs_volume = model_part.GetGeometry("NurbsVolume")
        force_x_total = 0.0
        force_y_total = 0.0
        force_z_total = 0.0
        force_total = 0.0

        for triangle in self.bcs_triangles:
            points = triangle.GetIntegrationPointsGlobal(1)
            for point in points:
                integration_points = []
                tmp_global = [point.GetX(), point.GetY(), point.GetZ()]
                center = FromGlobalToParamSpace(tmp_global, self.lower_point, self.upper_point)
                weight = triangle.Area()*point.GetWeight()*2
                normal = triangle.Normal()

                integration_points.append([center[0], center[1], center[2], point.GetWeight()])
                quadrature_point_geometries_boundary = KM.GeometriesVector()
                nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries_boundary, 2,
                integration_points)

                condition = model_part.CreateNewCondition("LoadCondition", id_counter, quadrature_point_geometries_boundary[0], properties)

                if self.normal_flag:
                    force_x = -normal[0] * weight * self.force
                    force_x_total += force_x
                    force_y = -normal[1] * weight * self.force
                    force_y_total += force_y
                    force_z = -normal[2] * weight * self.force
                    force_z_total += force_z
                    import numpy as np
                    force_total += np.sqrt(force_x*force_x + force_y*force_y + force_z*force_z)

                else:
                    force_x = -weight * self.force[0]
                    force_x_total += force_x
                    force_y = -weight * self.force[1]
                    force_y_total += force_y
                    force_z = -weight * self.force[2]
                    force_z_total += force_z
                    import numpy as np
                    force_total += np.sqrt(force_x*force_x + force_y*force_y + force_z*force_z)

                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, force_x)
                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, force_y)
                condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, force_z)
                self.p0.append(force_z)
                self.conditions.append(condition)
                id_counter += 1

        force_x_total = 0.0
        force_y_total = 0.0
        force_z_total = 0.0
        for condition in self.conditions:
            value_x = condition.GetValue(StructuralMechanicsApplication.POINT_LOAD_X)
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_X, value_x)
            force_x_total += value_x

            value_y = condition.GetValue(StructuralMechanicsApplication.POINT_LOAD_Y)
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Y, value_y)
            force_y_total += value_y

            value_z = condition.GetValue(StructuralMechanicsApplication.POINT_LOAD_Z)
            condition.SetValue(StructuralMechanicsApplication.POINT_LOAD_Z, value_z)
            force_z_total += value_z

        print("force_x_total: ", force_x_total)
        print("force_y_total: ", force_y_total)
        print("force_z_total: ", force_z_total)



