from typing import Tuple, cast

# Type definitions
Point3D = Tuple[float, float, float]

def point_from_global_to_param_space(
        point: Point3D,
        bound_xyz: Tuple[Point3D, Point3D],
        bound_uvw: Tuple[Point3D, Point3D]
    ) -> Point3D:
    """
    Maps a point from global Cartesian space to parametric space.

    Performs a linear transformation of a point from the physical domain
    defined by `bound_xyz` to the parametric domain defined by `bound_uvw`.

    Args:
        point (List[float]): A list of 3 floats representing the Cartesian coordinates [x, y, z].
        bound_xyz (List[List[float]]): A list of two 3D points [[x_min, y_min, z_min], [x_max, y_max, z_max]]
            representing the bounds in physical space.
        bound_uvw (List[List[float]]): A list of two 3D points [[u_min, v_min, w_min], [u_max, v_max, w_max]]
            representing the bounds in parametric space.

    Returns:
        List[float]: A list of 3 floats representing the corresponding parametric coordinates [u, v, w].
    """
    lower_xyz, upper_xyz = bound_xyz
    lower_uvw, upper_uvw = bound_uvw

    result = tuple(
        (point[i] - lower_xyz[i]) * abs(upper_uvw[i] - lower_uvw[i]) / abs(upper_xyz[i] - lower_xyz[i]) + lower_uvw[i]
        for i in range(3)
    )

    return cast(Point3D, result)

def point_from_param_to_global_space(
        point: Point3D,
        bound_xyz: Tuple[Point3D, Point3D],
        bound_uvw: Tuple[Point3D, Point3D]
    ) -> Point3D:
    """
    Maps a point from parametric space to global Cartesian space.

    Performs a linear transformation of a point from the parametric domain
    defined by `bound_uvw` to the physical domain defined by `bound_xyz`.

    Args:
        point (List[float]): A list of 3 floats representing the parametric coordinates [u, v, w].
        bound_xyz (List[List[float]]): A list of two 3D points [[x_min, y_min, z_min], [x_max, y_max, z_max]]
            representing the bounds in physical space.
        bound_uvw (List[List[float]]): A list of two 3D points [[u_min, v_min, w_min], [u_max, v_max, w_max]]
            representing the bounds in parametric space.

    Returns:
        List[float]: A list of 3 floats representing the corresponding Cartesian coordinates [x, y, z].
    """
    lower_xyz, upper_xyz = bound_xyz
    lower_uvw, upper_uvw = bound_uvw

    result = tuple(
        (point[i] - lower_uvw[i]) * abs(upper_xyz[i] - lower_xyz[i]) / abs(upper_uvw[i] - lower_uvw[i]) + lower_xyz[i]
        for i in range(3)
    )
    return cast(Point3D, result)
