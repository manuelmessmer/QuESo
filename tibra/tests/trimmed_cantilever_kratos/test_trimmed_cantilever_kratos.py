# Project imports
from platform import release
import re
from TIBRA_PythonApplication.PyTIBRA import PyTIBRA

try:
    import KratosMultiphysics as KM
    kratos_available = True
except:
    print("KratosMultiphysics is not available")
    kratos_available = False

import unittest
import math
import numpy as np

def neumann_condition(x, y, z):
    return (z > (10 - 1e-6) )

def dirichlet_condition(x, y, z):
    return (z < (0.0 + 1e-6))

class TestTrimmedCantileverKratos(unittest.TestCase):
    def test_1(self):
        #p=2
        #"number_of_knot_spans" : [2,2,10]
        #el=1000
        self.RunTest("tibra/tests/trimmed_cantilever_kratos/TIBRAParameters1.json", 0.002)

    def test_2(self):
        #p=2
        #"number_of_knot_spans" : [2,2,3]
        #el=3000
        self.RunTest("tibra/tests/trimmed_cantilever_kratos/TIBRAParameters2.json", 0.015)

    def test_3(self):
        #p=2
        #"number_of_knot_spans" : [8,8,10]
        #"integration_method" : "Gauss"
        #el=1000
        self.RunTest("tibra/tests/trimmed_cantilever_kratos/TIBRAParameters3.json", 0.0005)
        ips_inside = 0
        for element in self.pytibra.GetElements():
            if element.IsTrimmed():
                self.assertLessEqual(len(element.GetIntegrationPoints()), 27)
            else:
                ips_inside += len(element.GetIntegrationPoints())

        self.assertEqual(ips_inside, 2592)

    def test_4(self):
        #p=2
        #"number_of_knot_spans" : [8,8,10]
        #"integration_method : "ReducedExact"
        #el=1000
        self.RunTest("tibra/tests/trimmed_cantilever_kratos/TIBRAParameters4.json", 0.0005)

        ips_inside = 0
        for element in self.pytibra.GetElements():
            if element.IsTrimmed():
                self.assertLessEqual(len(element.GetIntegrationPoints()), 27)
            else:
                ips_inside += len(element.GetIntegrationPoints())
        self.assertEqual(ips_inside, 1275)

    def test_5(self):
        #p=2
        #"number_of_knot_spans" : [8,8,10]
        #"integration_method : "ReducedOrder1"
        #el=1000
        self.RunTest("tibra/tests/trimmed_cantilever_kratos/TIBRAParameters5.json", 0.0005)

        ips_inside = 0
        for element in self.pytibra.GetElements():
            if element.IsTrimmed():
                self.assertLessEqual(len(element.GetIntegrationPoints()), 27)
            else:
                ips_inside += len(element.GetIntegrationPoints())
        self.assertEqual(ips_inside, 572)

    def test_6(self):
        #p=2
        #"number_of_knot_spans" : [8,8,10]
        #"integration_method : "ReducedOrder2"
        #el=1000
        self.RunTest("tibra/tests/trimmed_cantilever_kratos/TIBRAParameters6.json", 0.0005)

        ips_inside = 0
        for element in self.pytibra.GetElements():
            if element.IsTrimmed():
                self.assertLessEqual(len(element.GetIntegrationPoints()), 27)
            else:
                ips_inside += len(element.GetIntegrationPoints())
        self.assertEqual(ips_inside, 243)

    def test_7(self):
        #p=3
        #"number_of_knot_spans" : [2,2,2]
        #"integration_method : "Gauss"
        self.RunTest("tibra/tests/trimmed_cantilever_kratos/TIBRAParameters7.json", 0.0008)
        for element in self.pytibra.GetElements():
            if element.IsTrimmed():
                self.assertLessEqual(len(element.GetIntegrationPoints()), 4*4*4)

    def test_8(self):
        #p=3
        #"number_of_knot_spans" : [2,2,2]
        #"integration_method : "Gauss"
        self.RunTest("tibra/tests/trimmed_cantilever_kratos/TIBRAParameters8.json", 0.0008)
        for element in self.pytibra.GetElements():
            if element.IsTrimmed():
                self.assertLessEqual(len(element.GetIntegrationPoints()), 5*5*5)

    def RunTest(self,filename, tolerance):
        if kratos_available:
            self.pytibra = PyTIBRA(filename)
            self.pytibra.Run()

            # Direct Analysis with kratos
            surface_force = [0, 0.1, 0]
            neumann_boundaries = [[neumann_condition, surface_force]]
            penalty_factor = 1e10
            dirichlet_boundaries = [[dirichlet_condition, penalty_factor]]
            self.pytibra.RunKratosAnalysis(dirichlet_boundaries, neumann_boundaries, "tibra/tests/trimmed_cantilever_kratos/KratosParameters.json")

            model_part = self.pytibra.GetAnalysis().GetModelPart()
            nurbs_volume = model_part.GetGeometry("NurbsVolume")

            self.CheckErrorInDisplacement(self.pytibra.GetLowerBound(), self.pytibra.GetUpperBound(), nurbs_volume, tolerance)

    def CheckErrorInDisplacement(self,lower_point, upper_point, nurbs_volume, tolerance):
        # Compare to analytical solution (Timoshenko Beam)
        I = math.pi / 4
        L = 10
        E = 100
        p = -0.1*math.pi
        nu = 0
        kappa_circle = (6 + 12*nu + 6*nu**2)/(7+12*nu+4*nu**2)
        G = E/ ( 2*(1+nu))

        x_test_array = np.arange(0.0,10.001,0.1)
        relative_errors = []
        xx = L
        u_ref_inf = -(p*xx*xx*(3*L-xx)/(6*E*I) + p*xx / (G*math.pi*kappa_circle))
        for xx in x_test_array:
            param = KM.Vector(3)
            param[0] = (0 - lower_point[0])/ abs(lower_point[0] - upper_point[0])
            param[1] = (0 - lower_point[1])/ abs(lower_point[1] - upper_point[1])
            param[2] = (xx- lower_point[2])/ abs(lower_point[2]-upper_point[2])
            u_y = nurbs_volume.GlobalCoordinates(param)[1]
            u_ref =  -(p*xx*xx*(3*L-xx)/(6*E*I) + p*xx / (G*math.pi*kappa_circle))
            relative_errors.append( abs(u_y-u_ref)/ abs(u_ref_inf))

        inf_norm = max(relative_errors)
        self.assertLess(inf_norm, tolerance)

if __name__ == "__main__":
    unittest.main()
