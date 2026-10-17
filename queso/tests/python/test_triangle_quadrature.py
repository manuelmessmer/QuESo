"""Python triangle quadrature bindings, defaults, and independent collections."""

from enum import Enum, IntEnum

import pytest
from pyqueso import BoundaryIntegrationPoint, TriangleQuadratureOrder, mesh


class OtherOrder(Enum):
    DEGREE_2 = 2


class OtherIntegerOrder(IntEnum):
    DEGREE_2 = 2


def _triangle_mesh():
    """Build a translated, tilted triangle with area 3.5 and normal (-3, -2, 6)/7."""
    owner = mesh.TriangleMesh()
    for point in ((1, 2, 3), (3, 2, 4), (1, 5, 4)):
        owner.add_vertex(point)
    owner.add_triangle((0, 1, 2))
    return owner


def test_native_quadrature_enum():
    assert issubclass(TriangleQuadratureOrder, Enum)
    assert not issubclass(TriangleQuadratureOrder, IntEnum)
    assert TriangleQuadratureOrder.__module__ == "pyqueso"
    assert not hasattr(mesh, "TriangleQuadratureOrder")
    assert [order.name for order in TriangleQuadratureOrder] == [
        "DEGREE_1",
        "DEGREE_2",
        "DEGREE_4",
        "DEGREE_6",
    ]


@pytest.mark.parametrize(
    "order,count",
    [
        (TriangleQuadratureOrder.DEGREE_1, 1),
        (TriangleQuadratureOrder.DEGREE_2, 3),
        (TriangleQuadratureOrder.DEGREE_4, 6),
        (TriangleQuadratureOrder.DEGREE_6, 12),
    ],
)
def test_global_points(order, count):
    owner = _triangle_mesh()
    triangle = owner.triangles()[0]
    global_points = triangle.integration_points_global(order=order)

    assert len(global_points) == count
    assert sum(point.weight for point in global_points) == pytest.approx(3.5)
    for physical in global_points:
        assert isinstance(physical, BoundaryIntegrationPoint)
        xi = (physical.x - 1) / 2
        eta = (physical.y - 2) / 3
        assert xi >= 0 and eta >= 0 and xi + eta <= 1
        assert physical.z == pytest.approx(3 + xi + eta)
        assert tuple(physical.normal) == pytest.approx((-3 / 7, -2 / 7, 6 / 7))
        assert physical.weight > 0
    for axis, centroid in enumerate((5 / 3, 3, 11 / 3)):
        moment = sum(point[axis] * point.weight for point in global_points)
        assert moment == pytest.approx(3.5 * centroid)


def test_defaults_and_independent_collections():
    owner = _triangle_mesh()
    triangle = owner.triangles()[0]
    get_points = triangle.integration_points_global
    points = get_points()
    explicit = get_points(order=TriangleQuadratureOrder.DEGREE_2)
    expected = [(tuple(point), point.weight) for point in explicit]
    assert len(points) == 3
    assert [(tuple(point), point.weight) for point in points] == expected
    expected_coordinates = ((4 / 3, 2.5, 10 / 3), (7 / 3, 2.5, 23 / 6), (4 / 3, 4, 23 / 6))
    for point, coordinates in zip(points, expected_coordinates):
        assert tuple(point) == pytest.approx(coordinates)
        assert point.weight == pytest.approx(7 / 6)
    # Each call owns a mutable collection, not the internal reference table.
    del points[0]
    fresh = get_points()
    assert len(points) == 2
    assert len(fresh) == 3
    assert [(tuple(point), point.weight) for point in fresh] == expected


def test_normals_follow_winding():
    owner = _triangle_mesh()
    owner.add_triangle((0, 2, 1))
    points = owner.triangles()[1].integration_points_global()
    assert sum(point.weight for point in points) == pytest.approx(3.5)
    for point in points:
        assert tuple(point.normal) == pytest.approx((3 / 7, 2 / 7, -6 / 7))


@pytest.mark.parametrize(
    "order",
    [-1, 0, 1, 2, 3, 4, 5, 6, "invalid", None, False, OtherOrder.DEGREE_2, OtherIntegerOrder.DEGREE_2],
)
def test_invalid_orders(order):
    owner = _triangle_mesh()
    with pytest.raises(TypeError):
        owner.triangles()[0].integration_points_global(order=order)


def test_only_global_quadrature_is_exposed():
    owner = _triangle_mesh()
    triangle = owner.triangles()[0]
    assert not hasattr(triangle, "integration_points_local")
    assert not hasattr(triangle, "integration_points_reference")
