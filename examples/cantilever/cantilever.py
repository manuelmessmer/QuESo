# Project imports
from TIBRA_PythonApplication.PyTIBRA import PyTIBRA
import matplotlib.pyplot as plt

def neumann_condition(x, y, z):
    return (z > (10 - 1e-6) )

def dirichlet_condition(x, y, z):
    return (z < (0.0 + 1e-6))

def main():
    pytibra = PyTIBRA("TIBRAParameters.json")
    pytibra.Run()

    # Direct Analysis with kratos
    surface_force = [0, -0.1, 0]
    neumann_boundaries = [[neumann_condition, surface_force]]
    penalty_factor = 1e10
    dirichlet_boundaries = [[dirichlet_condition, penalty_factor]]
    pytibra.RunKratosAnalysis(dirichlet_boundaries, neumann_boundaries)
    pytibra.PostProcess()



if __name__ == "__main__":
    main()
