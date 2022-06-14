# Project imports
import TrIGA_PythonApplication as TrIGA_Application
from src.python_scripts.helper import *

try:
    import KratosMultiphysics as KM
    kratos_available = True
except:
    print("KratosMultiphysics is not available")
    kratos_available = False

if kratos_available:
    from kratos_interface.kratos_analysis import Analysis
    from kratos_interface.weak_bcs import PenaltySupport
    from kratos_interface.weak_bcs import SurfaceLoad
    from kratos_interface.bounding_box_bcs import DirichletCondition
    from kratos_interface.bounding_box_bcs import NeumannCondition

# External imports
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import shutil

def neumann_condition(x, y, z):
    return (z > (10 - 1e-6) )

def dirichlet_condition(x, y, z):
    return (z < (0.0 + 1e-6))

def GetParams():
    analysis_parameters = KM.Parameters("""{
            "problem_data"    : {
                "parallel_type" : "OpenMP",
                "echo_level"    : 1,
                "start_time"    : 0.0,
                "end_time"      : 1.0
            },
            "solver_settings" : {
                "solver_type"              : "Static",
                "analysis_type"            : "linear",
                "model_part_name"          : "NurbsMesh",
                "echo_level"    : 1,
                "domain_size"              : 3,
                "model_import_settings"    : {
                    "input_type"     : "use_input_model_part"
                },
                "time_stepping"            : {
                    "time_step" : 1.1
                },
                "linear_solver_settings":{
                    "preconditioner_type" : "additive_schwarz",
                    "solver_type": "bicgstab",
                    "max_iteration" : 5000,
                    "tolerance" : 1e-10
                },
                "rotation_dofs"            : false,
                "builder_and_solver_settings" : {
                    "use_block_builder" : true
                },
                "residual_relative_tolerance"        : 0.000001
            },
            "modelers": [{
                "modeler_name": "NurbsGeometryModeler",
                "Parameters": {
                    "model_part_name" : "NurbsMesh",
                    "geometry_name"   : "NurbsVolume",
                    "lower_point": [-2, -2,  0],
                    "upper_point": [2, 2, 5],
                    "polynomial_order" : [2, 2, 2],
                    "number_of_knot_spans" : [17,17,20] }
            }]
        }""")

    return analysis_parameters

def main():
    folder_path = "./output/"
    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)
    os.mkdir(folder_path)

    filename = "data/cylinder.stl"

    lower_point = [-1.5, -1.5, -1]
    upper_point = [1.5, 1.5, 11]
    number_of_elements = [5, 5, 20]
    order = [2, 2, 2]

    point_distribution_factor = 2
    initial_triangle_edge_length = 1
    minimum_number_of_triangles = 2000
    moment_fitting_residual = 1e-12
    integration_method = "Gauss"
    echo_level = 2
    E = 100
    nu = 0.0
    material_properties = [E, nu, 7.8e-6]
    embedding_flag = True
    triga = TrIGA_Application.TrIGA(filename, lower_point, upper_point, number_of_elements, order,
                                    initial_triangle_edge_length,
                                    minimum_number_of_triangles,
                                    moment_fitting_residual,
                                    point_distribution_factor,
                                    integration_method,
                                    echo_level,
                                    embedding_flag)

    points_all = TrIGA_Application.VectorOfIntegrationPoints()
    dirichlet_triangles = TrIGA_Application.VectorOfTriangles()
    neumann_triangles = TrIGA_Application.VectorOfTriangles()

    elements = triga.GetElements()
    #Loop over all elements
    for element in elements:
        if element.IsTrimmed():
            for point_trimmed_reduced in element.GetIntegrationPointsTrimmed():
                weight = point_trimmed_reduced.GetWeight()
                if( weight > 0):
                    points_all.append(point_trimmed_reduced)
                else:
                    print("Show: ", point_trimmed_reduced.GetWeight())
            for triangle in element.GetDirichletTriangles(dirichlet_condition):
                dirichlet_triangles.append(triangle)
            for triangle in element.GetNeumannTriangles(neumann_condition):
                neumann_triangles.append(triangle)
        else:
            for point_inside in element.GetIntegrationPointsInside():
                points_all.append(point_inside)


    # Define boundary conditions
    boundary_conditions = []
    boundary_conditions.append( PenaltySupport(dirichlet_triangles, lower_point, upper_point, 1e10) )
    boundary_conditions.append( SurfaceLoad(neumann_triangles, lower_point, upper_point, [0, -0.1, 0], False ))

    # Start Kratos Analysis
    if kratos_available:
        postprocess_flag = False
        surface_mesh_file = "dummy"
        gid_output_dest = "./output/GiD/"
        analysis_parameters = GetParams()
        for modeler in analysis_parameters["modelers"]:
            if modeler["modeler_name"].GetString() == "NurbsGeometryModeler":
                modeler["Parameters"]["lower_point"].SetVector(lower_point)
                modeler["Parameters"]["upper_point"].SetVector(upper_point)
                modeler["Parameters"]["number_of_knot_spans"].SetVector(number_of_elements)
                modeler["Parameters"]["polynomial_order"].SetVector(order)

        analysis = Analysis(surface_mesh_file, points_all, analysis_parameters, boundary_conditions, material_properties, postprocess_flag, gid_output_dest)

        # Compare results to analytical solution
        #########################################################
        model_part = analysis.GetModelPart()
        nurbs_volume = model_part.GetGeometry("NurbsVolume")

        #print(nurbs_volume)
        param = KM.Vector(3)
        param[0] = 0.5
        param[1] = (0 - lower_point[0])/ abs(lower_point[0] - upper_point[0])
        # Parameters
        I = math.pi / 4
        L = 10
        p = -0.1*math.pi
        G = E/ ( 2*(1+0.0))
        kappa_circle = (6 + 12*nu + 6*nu**2)/(7+12*nu+4*nu**2)
        #https://asmedigitalcollection.asme.org/appliedmechanics/article/68/1/87/461375/Shear-Coefficients-for-Timoshenko-Beam-Theory
        M = 100
        x_array = np.arange(0,10.001,0.1)
        x_simulation = []
        y_simulation = []
        x_timoshenko = []
        y_timoshenko = []
        xx = 10
        rel_error = []
        point_a_ref = p*xx*xx*(3*L-xx)/(6*E*I) + p*xx / (G*math.pi*kappa_circle)
        for x in x_array:
            param[2] = (x - lower_point[2])/ abs(lower_point[2] - upper_point[2])
            x_simulation.append( nurbs_volume.GlobalCoordinates(param)[2] )
            y_simulation.append( -nurbs_volume.GlobalCoordinates(param)[1] )

            ref = p*x*x*(3*L-x)/(6*E*I) + p*x / (G*math.pi*kappa_circle)
            y_timoshenko.append( ref )

        # Plot figure
        fig, ax = plt.subplots(1,1,figsize = (16,8))
        ax.plot(x_simulation,y_simulation,color='r',label="Simulation")
        ax.plot(x_array,y_timoshenko,"--",color='b',label="Analytical Solution")

        plt.ylabel("Displacement y")
        plt.xlabel("x direction")
        #plt.ylim([-1e-2, 1e-2])
        plt.xlim([0, 10])
        plt.legend()
        plt.grid(True)
        plt.savefig("output/cylinder_displacement.png")
        plt.show()


if __name__ == "__main__":
    main()
