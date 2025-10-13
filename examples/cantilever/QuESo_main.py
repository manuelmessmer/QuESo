# Project imports
from QuESoPythonModule.PyQuESo import PyQuESo

def main():
    pyqueso = PyQuESo("QuESoSettings.json")
    pyqueso.Run()

    # Direct Analysis with kratos
    pyqueso.RunKratosAnalysis()

if __name__ == "__main__":
    main()
