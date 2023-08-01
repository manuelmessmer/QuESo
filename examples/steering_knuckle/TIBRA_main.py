# Project imports
from QuESo_PythonApplication.PyQuESo import PyQuESo

def main():
    pyqueso = PyQuESo("QuESoParameters.json")
    pyqueso.Run()

if __name__ == "__main__":
    main()
