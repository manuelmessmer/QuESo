"""Integration tests for the multi-component model API."""

from collections.abc import Callable
from pathlib import Path

import pytest
from pyqueso import Model


def test_create_and_named_access(
    copy_test_data: Callable[[str], Path],
) -> None:
    directory = copy_test_data("coupled_cantilever")
    model = Model(directory / "QuESoSettings.json")
    assert model.component_names == ("left", "right")
    with pytest.raises(RuntimeError):
        model.elements("left")

    model.create()
    assert len(model.elements("left")) > 0
    assert len(model.elements("right")) > 0
    assert [
        condition.settings.get_int("condition_id")
        for condition in model.conditions("left")
    ] == [1, 10]
    assert all(
        isinstance(segment.is_in_active_element, bool)
        for condition in model.conditions("left")
        for segment in condition.segments
    )
    with pytest.raises(RuntimeError):
        model.create()
    with pytest.raises(KeyError):
        model.settings("missing")
