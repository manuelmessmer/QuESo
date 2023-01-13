# Project imports
from TIBRA_PythonApplication.PyTIBRA import PyTIBRA
from tibra.python_scripts.helper import *

import unittest
import json

class TestBSplineVolume(unittest.TestCase):
    def test1(self):
        input_filename = "tibra/tests/b_spline_volume/TIBRAParameters1.json"
        results_filename = "tibra/tests/b_spline_volume/results_1.json"
        self.RunTest(input_filename, results_filename)

    def test2(self):
        input_filename = "tibra/tests/b_spline_volume/TIBRAParameters2.json"
        results_filename = "tibra/tests/b_spline_volume/results_2.json"
        self.RunTest(input_filename, results_filename)

    def test3(self):
        input_filename = "tibra/tests/b_spline_volume/TIBRAParameters3.json"
        results_filename = "tibra/tests/b_spline_volume/results_3.json"
        self.RunTest(input_filename, results_filename)

    def test4(self):
        input_filename = "tibra/tests/b_spline_volume/TIBRAParameters4.json"
        results_filename = "tibra/tests/b_spline_volume/results_4.json"
        self.RunTest(input_filename, results_filename)

    def RunTest(self,input_filename, results_filename):
        pytibra = PyTIBRA(input_filename)

        volume = pytibra.GetBSplineVolume()
        cps = volume.ControlPoints()
        knots_u = volume.KnotsU()
        knots_v = volume.KnotsV()
        knots_w = volume.KnotsW()
        polynomial_order = volume.PolynomialOrder()

        # Read results file
        with open(results_filename, 'r') as file:
            res = json.load(file)

        # Check number of control points
        n_cps_u = volume.NumberControlPointsInU()
        n_cps_v = volume.NumberControlPointsInV()
        n_cps_w = volume.NumberControlPointsInW()
        n_cps = n_cps_u*n_cps_v*n_cps_w
        self.assertEqual(n_cps, len(res["cps"]))
        # Check number of knot vectros
        self.assertAlmostEqual(len(knots_u), n_cps_u + polynomial_order[0] + 1)
        self.assertAlmostEqual(len(knots_v), n_cps_v + polynomial_order[1] + 1)
        self.assertAlmostEqual(len(knots_w), n_cps_w + polynomial_order[2] + 1)

        # Check control points list (linearized)
        self.assertEqual(len(cps), n_cps)
        for cp1, cp2 in zip(res["cps"], cps):
            self.assertAlmostEqual(cp1[0], cp2[0], 12)
            self.assertAlmostEqual(cp1[1], cp2[1], 12)
            self.assertAlmostEqual(cp1[2], cp2[2], 12)

        # Check control points matrix
        cps_matrix = volume.ControlPointsMatrix()
        self.assertEqual(n_cps_u, cps_matrix.shape[0])
        self.assertEqual(n_cps_v, cps_matrix.shape[1])
        self.assertEqual(n_cps_w, cps_matrix.shape[2])
        count = 0
        for i_w in range(n_cps_w):
            for i_v in range(n_cps_v):
                for i_u in range(n_cps_u):
                    test_cp = cps_matrix[i_u, i_v, i_w]
                    self.assertAlmostEqual(cps[count][0], test_cp[0], 12)
                    self.assertAlmostEqual(cps[count][1], test_cp[1], 12)
                    self.assertAlmostEqual(cps[count][2], test_cp[2], 12)
                    count += 1

        # Check knot vectors
        self.assertEqual(len(knots_u), len(res["knots_u"]))
        for k1, k2 in zip(res["knots_u"], knots_u):
            self.assertAlmostEqual(k1, k2, 12)

        self.assertEqual(len(knots_v), len(res["knots_v"]))
        for k1, k2 in zip(res["knots_v"], knots_v):
            self.assertAlmostEqual(k1, k2, 12)

        self.assertEqual(len(knots_w), len(res["knots_w"]))
        for k1, k2 in zip(res["knots_w"], knots_w):
            self.assertAlmostEqual(k1, k2, 12)

if __name__ == "__main__":
    unittest.main()