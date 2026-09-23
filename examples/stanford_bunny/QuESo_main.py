# Project imports
import pyqueso


def main():
    model = pyqueso.Model("QuESoSettings.json")
    model.create()

if __name__ == "__main__":
    main()
