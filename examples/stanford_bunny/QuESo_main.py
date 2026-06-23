# Project imports
from QuESoPythonModule.model import Model


def main():
    pyqueso = Model("QuESoSettings.json")
    pyqueso.run()


if __name__ == "__main__":
    main()
