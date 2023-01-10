# Project imports
from TIBRA_PythonApplication.PyTIBRA import PyTIBRA
from tibra.python_scripts.b_spline_volume import BSplineVolume
from tibra.python_scripts.helper import *

try:
    import KratosMultiphysics as KM
    kratos_available = True
except:
    print("KratosMultiphysics is not available")
    kratos_available = False

import unittest
import math
import numpy as np
import json


class TestTrimmedCantileverKratos(unittest.TestCase):

    def test1(self):
        input_filename = "tibra/tests/b_spline_volume/TIBRAParameters1.json"
        results_filename = "tibra/tests/b_spline_volume/result_1.json"
        self.RunTest(input_filename, results_filename)

    def RunTest(self,input_filename, results_filename):
        if kratos_available:
            self.pytibra = PyTIBRA(input_filename)
            self.pytibra.Run()

            # Direct Analysis with kratos
            # surface_force = [0, 0.1, 0]
            neumann_boundaries = []
            # penalty_factor = 1e10
            dirichlet_boundaries = []
            self.pytibra.RunKratosAnalysis(dirichlet_boundaries, neumann_boundaries, "tibra/tests/trimmed_cantilever_kratos/KratosParameters.json")

            model_part = self.pytibra.GetAnalysis().GetModelPart()
            nurbs_volume = model_part.GetGeometry("NurbsVolume")


            parameters = ReadParameters(input_filename)
            volume = BSplineVolume(parameters)
            cps = volume.ControlPoints()
            knots_u = volume.KnotsU()
            knots_v = volume.KnotsV()
            knots_w = volume.KnotsW()

            for p1, p2 in zip(nurbs_volume, cps):
                self.assertAlmostEqual(p1.X, p2[0], 12)
                self.assertAlmostEqual(p1.Y, p2[1], 12)
                self.assertAlmostEqual(p1.Z, p2[2], 12)

            nurbs_u_ref = np.array([0])
            nurbs_u_ref = np.append(nurbs_u_ref, nurbs_volume.KnotsU())
            nurbs_u_ref = np.append(nurbs_u_ref, 1)
            for k1, k2 in zip(nurbs_u_ref, knots_u):
                self.assertAlmostEqual(k1, k2, 12)

            nurbs_v_ref = np.array([0])
            nurbs_v_ref = np.append(nurbs_v_ref, nurbs_volume.KnotsV())
            nurbs_v_ref = np.append(nurbs_v_ref, 1)
            for k1, k2 in zip(nurbs_v_ref, knots_v):
                self.assertAlmostEqual(k1, k2, 12)

            nurbs_w_ref = np.array([0])
            nurbs_w_ref = np.append(nurbs_w_ref, nurbs_volume.KnotsW())
            nurbs_w_ref = np.append(nurbs_w_ref, 1)
            print(nurbs_w_ref)
            print(knots_w)
            for k1, k2 in zip(nurbs_w_ref, knots_w):
                self.assertAlmostEqual(k1, k2, 12)

            # dict_volume = {}
            # dict_volume["cps"] = cps
            # dict_volume["knots_u"] = knots_u
            # dict_volume["knots_v"] = knots_v
            # dict_volume["knots_w"] = knots_w
            # with open("test.json", 'r') as file:
            #     json.dump(dict_volume, file)

            with open(results_filename, 'r') as file:
                res = json.load(file)

            # Check Control Points
            self.assertEqual(len(cps), 192)
            for cp1, cp2 in zip(res["cps"], cps):
                self.assertAlmostEqual(cp1[0], cp2[0], 12)
                self.assertAlmostEqual(cp1[1], cp2[1], 12)
                self.assertAlmostEqual(cp1[2], cp2[2], 12)

            # Check Knot Vectors
            self.assertEqual(len(knots_u), 7)
            for k1, k2 in zip(res["knots_u"], knots_u):
                self.assertAlmostEqual(k1, k2, 12)

            self.assertEqual(len(knots_v), 7)
            for k1, k2 in zip(res["knots_v"], knots_v):
                self.assertAlmostEqual(k1, k2, 12)

            self.assertEqual(len(knots_w), 15)
            for k1, k2 in zip(res["knots_w"], knots_w):
                self.assertAlmostEqual(k1, k2, 12)

if __name__ == "__main__":
    unittest.main()