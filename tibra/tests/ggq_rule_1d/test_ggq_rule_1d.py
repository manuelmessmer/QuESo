# Project imports
import TIBRA_PythonApplication

import unittest
import numpy as np
from scipy import __version__ as __scipy_version__
from scipy._lib import _pep440
from scipy.interpolate import BSpline

#import scipy

class TestGGQ1d(unittest.TestCase):
    def check_ggq_rules(self, points, p, r, e, a, b):
        cps = []
        weights = []
        for point in points:
            cps.append(point[0])
            weights.append(point[1])

        knots = a*np.ones(p)
        knots = np.concatenate( (knots, np.linspace(a, b, e+1)), axis=0)
        knots = np.concatenate( (knots, b*np.ones(p)), axis=0)

        np_minversion = '1.8.0'
        if _pep440.parse(__scipy_version__) < _pep440.Version(np_minversion):
            print("Testing :: QQGRule1d is not executed. Requires scipy version: >=" + np_minversion)

            # This is a 'weaker' test. Does not check the positions, but only the weights.
            spline = BSpline(knots, cps, p)
            target = spline.integrate(0,1.0, extrapolate=False)
            area = 0.0
            for point in points:
                area += spline(point[0])*point[1]
            self.assertLess( np.abs(area - target), 1e-14 )
        else:
            target = []
            dim = len(knots) - p - 1
            for k in range(dim):
                target.append( (knots[p+k+1] - knots[k])/(p+1) )
            # Get B-Spline collocation matrix
            B = BSpline.design_matrix(cps, knots, p).toarray()
            error = np.linalg.norm( (target - np.array(weights).dot(np.array(B)) )/dim )
            self.assertLess( error, 1e-14 )
    #p=2
    #def test_1(self):
        # PolynomialDegree = 2
        # for e in range(1,101):
        #     points = TIBRA_PythonApplication.IntegrationPointFactory.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedExact)
        #     print(e)
        #     p = 4
        #     r = 0
        #     self.check_ggq_rules(points, p, r, e, 0.0, 1.0)

    def test_2(self):
        PolynomialDegree = 3
        for e in range(1,101):
            points = TIBRA_PythonApplication.IntegrationPointFactory.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedExact)
            for point in points:
                print(point)
            print(e)
            p = 6
            r = 1
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0)


if __name__ == "__main__":
    unittest.main()
