# Project imports
import TIBRA_PythonApplication

import unittest
from scipy.interpolate import BSpline
import numpy as np

class TestGGQ1d(unittest.TestCase):
    #p=2
    def test_1(self):
        p = 4
        r = 0
        e = 34
        a = 0.0
        b = 1.0
        points = TIBRA_PythonApplication.GGQRule.GetGGQ_Rule(p, r, e, a, b)
        cps = []
        weights = []
        for point in points:
            cps.append(point[0])
            weights.append(point[1])

        knots = a*np.ones(p)
        knots = np.concatenate( (knots, np.linspace(a, b, e+1)), axis=0)
        knots = np.concatenate( (knots, b*np.ones(p)), axis=0)

        spline = BSpline(knots, cps, p)
        print( spline.integrate(0,1.0, extrapolate=False) )

        test = 0.0
        for point in points:
            test += spline(point[0])*point[1]
        print(test)


if __name__ == "__main__":
    unittest.main()
