"""Run the coupled two-component cantilever Kratos analysis."""

from pathlib import Path

from pyqueso.kratos_interface import Analysis

# INFO:
# To run this example, use the Kratos fork from
# https://github.com/manuelmessmer/Kratos. CouplingPenaltyCondition requires
# KratosMultiphysics.CouplingGeometry, which is not yet part of Kratos master.


def main() -> None:
    directory = Path(__file__).parent
    analysis = Analysis(
        queso_settings_path=directory / "QuESoSettings.json",
        analysis_settings_path=directory / "AnalysisSettings.json",
        kratos_parameters_path=directory / "KratosParameters.json",
    )
    analysis.run()


if __name__ == "__main__":
    main()
