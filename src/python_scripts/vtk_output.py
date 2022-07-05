import pyevtk
from pyevtk.hl import pointsToVTK
import numpy as np
import meshio
from src.python_scripts.helper import *

class KnotspanWriter:
    def __init__(self, lower_point, upper_point):
        self.vtk_cells = []
        self.vtk_points = []
        self.lower_point = lower_point
        self.upper_point = upper_point

    def add_element(self,element):
        point_a = element.GetLocalLowerPoint()
        point_b = element.GetLocalUpperPoint()

        tmp_points = [ [point_a[0], point_a[1], point_a[2]],
                       [point_b[0], point_a[1], point_a[2]],
                       [point_b[0], point_b[1], point_a[2]],
                       [point_a[0], point_b[1], point_a[2]],
                       [point_a[0], point_a[1], point_b[2]],
                       [point_b[0], point_a[1], point_b[2]],
                       [point_b[0], point_b[1], point_b[2]],
                       [point_a[0], point_b[1], point_b[2]] ]

        current_id = len(self.vtk_points)
        for point in tmp_points:
            self.vtk_points.append(FromParamToGlobalSpace(point, self.lower_point, self.upper_point))
        self.vtk_cells.append(np.arange(current_id,current_id+8))

    def write_file(self, filename):
        cell_container = {"hexahedron" : np.array(self.vtk_cells)}
        meshio.write_points_cells( "./output/" + filename + ".vtu", np.array(self.vtk_points), cell_container)

class PointsWriter:
    def write_points_to_vtk_with_mapping(points, filename, lower_point, upper_point):
        point_x = []
        point_y = []
        point_z = []
        weights = []
        for point in points:
            if( point.GetWeight() > 0.0 ):
                tmp_point = [point.GetX(), point.GetY(), point.GetZ()]
                tmp_point = FromParamToGlobalSpace(tmp_point, lower_point, upper_point)
                point_x.append(tmp_point[0])
                point_y.append(tmp_point[1])
                point_z.append(tmp_point[2])
                weights.append(point.GetWeight())

        if( len(point_x) > 0):
            pointsToVTK("./output/" + (filename), np.array(point_x), np.array(point_y), np.array(point_z), {"weights" : np.array(weights)})

    def write_points_to_vtk(points, filename):
        point_x = []
        point_y = []
        point_z = []
        weights = []
        for point in points:
            if( point.GetWeight() > 0.0 ):
                point_x.append(point.GetX())
                point_y.append(point.GetY())
                point_z.append(point.GetZ())
                weights.append(point.GetWeight())

        if( len(point_x) > 0):
            pointsToVTK("./output/" + (filename), np.array(point_x), np.array(point_y), np.array(point_z), {"weights" : np.array(weights)})

    def write_bcs_to_vtk(triangles, filename):
        points_x = []
        points_y = []
        points_z = []
        weights = []
        for triangle in triangles:
            points_x.append(triangle.Center()[0])
            points_y.append(triangle.Center()[1])
            points_z.append(triangle.Center()[2])
            weights.append(triangle.Area())

        if( len(points_x) > 0):
            pointsToVTK("./output/" + (filename), np.array(points_x), np.array(points_y), np.array(points_z), {"weights" : np.array(weights)})