# Project imports
from TrIGA_PythonApplication.PyTrIGA import PyTrIGA
import matplotlib.pyplot as plt

def neumann_condition(x, y, z):
    return (z > (10 - 1e-6) )

def dirichlet_condition(x, y, z):
    return (z < (0.0 + 1e-6))

def main():
    pytriga = PyTrIGA("TrIGAParameters.json")
    pytriga.Run()

    # Direct Analysis with kratos
    surface_force = [0, -0.1, 0]
    neumann_boundaries = [[neumann_condition, surface_force]]
    penalty_factor = 1e10
    dirichlet_boundaries = [[dirichlet_condition, penalty_factor]]
    pytriga.RunKratosAnalysis(dirichlet_boundaries, neumann_boundaries)
    pytriga.PostProcess()



if __name__ == "__main__":
    main()
