"""Ownership policies for borrowed range items and independent tuple values."""

import gc
import weakref

import pytest
from pyqueso import mesh, Model


def _triangle_mesh():
    result = mesh.TriangleMesh()
    for point in ((0, 0, 0), (1, 0, 0), (0, 1, 0)):
        result.add_vertex(point)
    result.add_triangle((0, 1, 2))
    return result


@pytest.mark.parametrize("access", ["iterate", "index"])
def test_triangle_proxy_retains_mesh(access):
    owner = _triangle_mesh()
    owner_ref = weakref.ref(owner)
    triangle = (
        next(iter(owner.triangles())) if access == "iterate" else owner.triangles()[-1]
    )
    del owner
    gc.collect()

    assert owner_ref() is not None
    assert triangle.p1 == (0, 0, 0)
    assert triangle.p2 == (1, 0, 0)
    assert triangle.p3 == (0, 1, 0)
    assert triangle.normal == (0, 0, 1)
    assert triangle.area() == pytest.approx(0.5)

    del triangle
    gc.collect()
    assert owner_ref() is None


@pytest.mark.parametrize("access", ["iterate", "index"])
def test_vertex_tuple_is_independent_of_mesh(access):
    owner = _triangle_mesh()
    owner_ref = weakref.ref(owner)
    vertex = (
        next(iter(owner.vertices())) if access == "iterate" else owner.vertices()[0]
    )
    assert isinstance(vertex, tuple)
    del owner
    gc.collect()

    assert owner_ref() is None
    assert vertex == (0, 0, 0)


@pytest.mark.parametrize("access", ["iterate", "index"])
def test_condition_segment_retains_condition(copy_test_data, access):
    directory = copy_test_data("coupled_cantilever")
    model = Model(directory / "QuESoSettings.json")
    model.create()
    condition = model.conditions("left")[0]
    condition_ref = weakref.ref(condition)
    segment = (
        next(iter(condition.segments))
        if access == "iterate"
        else condition.segments[-1]
    )
    expected_count = segment.triangle_mesh.num_triangles
    expected_active = segment.is_in_active_element
    del condition, model
    gc.collect()

    assert condition_ref() is not None
    assert segment.triangle_mesh.num_triangles == expected_count
    assert segment.is_in_active_element == expected_active

    del segment
    gc.collect()
    assert condition_ref() is None
