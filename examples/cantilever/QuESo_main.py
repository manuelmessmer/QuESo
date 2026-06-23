# Project imports
import KratosMultiphysics as KM
from QuESoPythonModule.kratos_interface.kratos_analysis import KratosAnalysis


def main():
    analysis = KratosAnalysis(KM.Model())
    analysis.Run()


if __name__ == "__main__":
    main()
