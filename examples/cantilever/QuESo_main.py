# Project imports
import pyqueso

def main():
    model = pyqueso.Model(json_filename="QuESoSettings.json")
    model.create()

    # Direct Analysis with kratos
    model.run_kratos_analysis()

if __name__ == "__main__":
    main()
