# Project imports
from TIBRA_PythonApplication.PyTIBRA import PyTIBRA
from tibra.python_scripts.b_spline_volume import BSplineVolume
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
        # Construct B-Spline Volume from parameters
        parameters = ReadParameters(input_filename)
        volume = BSplineVolume(parameters)
        cps = volume.ControlPoints()
        knots_u = volume.KnotsU()
        knots_v = volume.KnotsV()
        knots_w = volume.KnotsW()

        # Read results file
        with open(results_filename, 'r') as file:
            res = json.load(file)

        # Check Control Points
        self.assertEqual(len(cps), len(res["cps"]))
        for cp1, cp2 in zip(res["cps"], cps):
            self.assertAlmostEqual(cp1[0], cp2[0], 12)
            self.assertAlmostEqual(cp1[1], cp2[1], 12)
            self.assertAlmostEqual(cp1[2], cp2[2], 12)

        # Check Knot Vectors
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