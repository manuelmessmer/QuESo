"""Smoke tests for the explicitly supported published examples."""

import runpy
import shutil
from pathlib import Path

import pytest

pytestmark = pytest.mark.examples


EXAMPLES = [
    pytest.param("stanford_bunny", False, id="stanford-bunny"),
    pytest.param("steering_knuckle", False, id="steering-knuckle"),
    pytest.param(
        "kratos_analysis/cantilever",
        True,
        marks=pytest.mark.kratos,
        id="kratos-cantilever",
    ),
    pytest.param(
        "kratos_analysis/cantilever_coupled",
        True,
        marks=pytest.mark.kratos,
        id="kratos-coupled-cantilever",
    ),
]


def _copy_example(source: Path, destination: Path) -> None:
    destination.mkdir()
    for path in source.iterdir():
        if path.name in {"queso_output", "kratos_output", "__pycache__"}:
            continue
        target = destination / path.name
        if path.is_dir():
            shutil.copytree(path, target)
        else:
            shutil.copy2(path, target)


@pytest.mark.parametrize(("relative_path", "requires_kratos"), EXAMPLES)
def test_example_runs(
    relative_path: str,
    requires_kratos: bool,
    project_root: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    if requires_kratos:
        pytest.importorskip("KratosMultiphysics")

    source = project_root / "examples" / relative_path
    destination = tmp_path / source.name
    _copy_example(source, destination)
    monkeypatch.chdir(destination)
    runpy.run_path(str(destination / "QuESo_main.py"), run_name="__main__")
