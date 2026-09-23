"""Run the single-component cantilever Kratos analysis."""

from pathlib import Path

from pyqueso.kratos_interface import Analysis


def main() -> None:
    """Construct and run the standalone Kratos analysis."""
    directory = Path(__file__).parent
    analysis = Analysis(
        queso_settings_path=directory / "QuESoSettings.json",
        analysis_settings_path=directory / "AnalysisSettings.json",
        kratos_parameters_path=directory / "KratosParameters.json",
    )
    analysis.run()


if __name__ == "__main__":
    main()
