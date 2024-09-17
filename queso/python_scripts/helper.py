def PointFromGlobalToParamSpace(point, bound_xyz, bound_uvw):
    tmp_point = [0.0, 0.0, 0.0]
    lower_point_xyz = bound_xyz[0]
    upper_point_xyz = bound_xyz[1]
    lower_point_uvw = bound_uvw[0]
    upper_point_uvw = bound_uvw[1]
    tmp_point[0] = (point[0] - lower_point_xyz[0]) * abs(upper_point_uvw[0] - lower_point_uvw[0]) / abs(upper_point_xyz[0]-lower_point_xyz[0]) + lower_point_uvw[0]
    tmp_point[1] = (point[1] - lower_point_xyz[1]) * abs(upper_point_uvw[1] - lower_point_uvw[1]) / abs(upper_point_xyz[1]-lower_point_xyz[1]) + lower_point_uvw[1]
    tmp_point[2] = (point[2] - lower_point_xyz[2]) * abs(upper_point_uvw[2] - lower_point_uvw[2]) / abs(upper_point_xyz[2]-lower_point_xyz[2]) + lower_point_uvw[2]

    return tmp_point

def PointFromParamToGlobalSpace(point, bound_xyz, bound_uvw):
    tmp_point = [0.0, 0.0, 0.0]
    lower_point_xyz = bound_xyz[0]
    upper_point_xyz = bound_xyz[1]
    lower_point_uvw = bound_uvw[0]
    upper_point_uvw = bound_uvw[1]
    tmp_point[0] = (point[0] - lower_point_uvw[0]) * abs(upper_point_xyz[0] - lower_point_xyz[0]) / abs(upper_point_uvw[0]-lower_point_uvw[0]) + lower_point_xyz[0]
    tmp_point[1] = (point[1] - lower_point_uvw[1]) * abs(upper_point_xyz[1] - lower_point_xyz[1]) / abs(upper_point_uvw[1]-lower_point_uvw[1]) + lower_point_xyz[1]
    tmp_point[2] = (point[2] - lower_point_uvw[2]) * abs(upper_point_xyz[2] - lower_point_xyz[2]) / abs(upper_point_uvw[2]-lower_point_uvw[2]) + lower_point_xyz[2]

    return tmp_point
