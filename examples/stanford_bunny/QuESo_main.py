# Project imports
from QuESoPythonModule.PyQuESo import PyQuESo

def main():
    pyqueso = PyQuESo("QuESoSettings.json")
    pyqueso.Run()

if __name__ == "__main__":
    main()
