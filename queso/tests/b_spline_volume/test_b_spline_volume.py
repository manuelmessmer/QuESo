# Project imports
from QuESo_PythonApplication.PyQuESo import PyQuESo
from queso.python_scripts.helper import *
from queso.python_scripts.queso_unit_test import QuESoTestCase
# External imports
import unittest
import numpy as np
import json

class TestBSplineVolume(QuESoTestCase):
    def test1(self):
        input_filename = "queso/tests/b_spline_volume/QuESoSettings1.json"
        results_filename = "queso/tests/b_spline_volume/results_1.json"
        self.RunTest(input_filename, results_filename)

    def test2(self):
        input_filename = "queso/tests/b_spline_volume/QuESoSettings2.json"
        results_filename = "queso/tests/b_spline_volume/results_2.json"
        self.RunTest(input_filename, results_filename)

    def test3(self):
        input_filename = "queso/tests/b_spline_volume/QuESoSettings3.json"
        results_filename = "queso/tests/b_spline_volume/results_3.json"
        self.RunTest(input_filename, results_filename)

    def test4(self):
        input_filename = "queso/tests/b_spline_volume/QuESoSettings4.json"
        results_filename = "queso/tests/b_spline_volume/results_4.json"
        self.RunTest(input_filename, results_filename)

    def test5(self):
        input_filename = "queso/tests/b_spline_volume/QuESoSettings5.json"
        results_filename = "queso/tests/b_spline_volume/results_5.json"
        self.RunTest(input_filename, results_filename)

    def test6(self):
        input_filename = "queso/tests/b_spline_volume/QuESoSettings6.json"
        results_filename = "queso/tests/b_spline_volume/results_6.json"
        self.RunTest(input_filename, results_filename)

    def test7(self):
        input_filename = "queso/tests/b_spline_volume/QuESoSettings7.json"
        results_filename = "queso/tests/b_spline_volume/results_7.json"
        self.RunTest(input_filename, results_filename)

    def test8(self):
        input_filename = "queso/tests/b_spline_volume/QuESoSettings8.json"
        results_filename = "queso/tests/b_spline_volume/results_8.json"
        self.RunTest(input_filename, results_filename)

    def __CompareOpenVsNonOpen(self, open_spline, non_open_spline):
        knots_u_open = open_spline.t
        knots_u_non_open = non_open_spline.t
        delta_init = knots_u_non_open[1] - knots_u_non_open[0]
        for i in range(len(knots_u_non_open)-1):
            delta = knots_u_non_open[i+1] - knots_u_non_open[i]
            self.assertAlmostEqual(delta_init, delta, 12)

        # Check if two spline represent same line
        uu = np.arange(knots_u_open[0], knots_u_open[-1], 100)
        for u in uu:
            x_1 = non_open_spline(u)
            x_2 = open_spline(u)
            self.assertAlmostEqual(x_1, x_2, 12)

        # Check if every knot in open exists in non_open
        for u1 in knots_u_open:
            found = False
            for u2 in knots_u_non_open:
                if np.allclose(u1, u2):
                    found = True

            self.assertTrue(found)


    def RunTest(self, input_filename, results_filename):
        pyqueso = PyQuESo(input_filename)
        volume_open = pyqueso.GetBSplineVolume("open_knot_vector")
        cps = volume_open.ControlPoints()
        knots_u = volume_open.KnotsU()
        knots_v = volume_open.KnotsV()
        knots_w = volume_open.KnotsW()
        polynomial_order = volume_open.PolynomialOrder()

        # To write results:
        # dict_a = {}
        # dict_a["cps"] = cps
        # dict_a["knots_u"] = knots_u
        # dict_a["knots_v"] = knots_v
        # dict_a["knots_w"] = knots_w

        # # Read results file
        # with open("test.json", 'w') as file:
        #     json.dump(dict_a, file )

        # Read results file
        with open(results_filename, 'r') as file:
            res = json.load(file)

        # Check number of control points
        n_cps_u = volume_open.NumberControlPointsInU()
        n_cps_v = volume_open.NumberControlPointsInV()
        n_cps_w = volume_open.NumberControlPointsInW()
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
        cps_matrix = volume_open.ControlPointsMatrix()
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

        volume_non_open = pyqueso.GetBSplineVolume("non_open_knot_vector")

        self.__CompareOpenVsNonOpen(volume_open.GetSpline(0), volume_non_open.GetSpline(0))
        self.__CompareOpenVsNonOpen(volume_open.GetSpline(1), volume_non_open.GetSpline(1))
        self.__CompareOpenVsNonOpen(volume_open.GetSpline(2), volume_non_open.GetSpline(2))

if __name__ == "__main__":
    unittest.main()