# Project imports
from TIBRA_PythonApplication.PyTIBRA import PyTIBRA

def main():
    pytibra = PyTIBRA("TIBRAParameters.json")
    pytibra.Run()

    # Direct Analysis with kratos
    pytibra.RunKratosAnalysis()

if __name__ == "__main__":
    main()
