"""Shared pytest configuration and test-data fixtures."""

from __future__ import annotations

import shutil
from collections.abc import Callable
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Return the repository root."""
    return Path(__file__).resolve().parents[3]


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """Return the canonical test-data directory."""
    return Path(__file__).resolve().parents[1] / "data"


@pytest.fixture
def copy_test_data(
    test_data_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> Callable[[str], Path]:
    """Copy a test-data file or directory into an isolated working directory."""
    monkeypatch.chdir(tmp_path)

    def copy(relative_path: str) -> Path:
        source = test_data_dir / relative_path
        destination = tmp_path / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source.is_dir():
            shutil.copytree(source, destination)
        else:
            shutil.copy2(source, destination)
        return destination

    return copy


def pytest_configure(config: pytest.Config) -> None:
    """Register the suite's execution-context markers."""
    config.addinivalue_line("markers", "kratos: requires KratosMultiphysics")
    config.addinivalue_line("markers", "examples: runs a published example")
    if not config.option.markexpr:
        config.option.markexpr = "not kratos"
