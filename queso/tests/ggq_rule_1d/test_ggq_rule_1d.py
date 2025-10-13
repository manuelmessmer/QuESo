# Project imports
import QuESoPythonModule

import unittest
import numpy as np
from scipy import __version__ as __scipy_version__
from scipy._lib import _pep440
from scipy.interpolate import BSpline

#import scipy

class TestGGQ1d(unittest.TestCase):
    def check_ggq_rules(self, points, p, r, e, a, b, must_pass):
        n = (p+1)*2 + (e-1)*(p-r) - p-1 # mumber dofs
        m = int(np.ceil(n/2))           # numer quadrature points

        nodes = []
        weights = []
        for point in points:
            nodes.append(point[0])
            weights.append(point[1])

        knots = a*np.ones(p+1)
        knots = np.concatenate( (knots, np.repeat(np.linspace(a, b, e+1)[1:-1], p-r)), axis=0)
        knots = np.concatenate( (knots, b*np.ones(p+1)), axis=0)

        np_minversion = '1.8.0'
        if _pep440.parse(__scipy_version__) >= _pep440.Version(np_minversion):
            target = []
            dim = len(knots) - p - 1
            for k in range(dim):
                target.append( (knots[p+k+1] - knots[k])/(p+1) )
            # Get B-Spline collocation matrix
            B = BSpline.design_matrix(nodes, knots, p).toarray()
            error = np.linalg.norm( (target - np.array(weights).dot(np.array(B)) )/dim )
            if must_pass:
                self.assertLess( error, 1e-15 )
                self.assertEqual(len(nodes), m)
            else:
                self.assertGreater( error, 1e-10)
        else:
            print("Testing :: QQGRule1d is skipped. Requires scipy version: >=" + np_minversion)

    def test_1(self):
        '''p=2: GGQ Optimal'''
        PolynomialDegree = 2
        for e in range(1,101):
            points = QuESoPythonModule.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, QuESoPythonModule.IntegrationMethod.GGQ_Optimal)
            p, r = 4, 0
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            if e > 1:
                p, r = 5, 0
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_2(self):
        '''p=3: GGQ Optimal'''
        PolynomialDegree = 3
        for e in range(1,101):
            points = QuESoPythonModule.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, QuESoPythonModule.IntegrationMethod.GGQ_Optimal)
            p, r = 6, 1
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            if e > 1:
                p, r = 7, 1
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)
                p, r = 7, 0
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_3(self):
        '''p=4: GGQ Optimal'''
        PolynomialDegree = 4
        for e in range(1,101):
            points = QuESoPythonModule.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, QuESoPythonModule.IntegrationMethod.GGQ_Optimal)
            p, r = 8, 2
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            if e > 1:
                p, r = 9, 2
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)
            if e > 2:
                p, r = 8, 1
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_4(self):
        '''p=2: GGQ Reduced1'''
        PolynomialDegree = 2
        for e in range(1,101):
            points = QuESoPythonModule.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, QuESoPythonModule.IntegrationMethod.GGQ_Reduced1)
            p, r = 3, 0
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            p, r = 4, 0 # Expected to fail.
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_5(self):
        '''p=3: GGQ Reduced1'''
        PolynomialDegree = 3
        for e in range(1,101):
            points = QuESoPythonModule.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, QuESoPythonModule.IntegrationMethod.GGQ_Reduced1)
            p, r = 5, 1
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            p, r = 6, 1 # Expected to fail.
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)
            if e > 1:
                p, r = 5, 0 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_6(self):
        '''p=4: GGQ Reduced1'''
        PolynomialDegree = 4
        for e in range(2,101):
            points = QuESoPythonModule.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, QuESoPythonModule.IntegrationMethod.GGQ_Reduced1)
            p, r = 7, 2
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            p, r = 8, 2 # Expected to fail.
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)
            if e > 2:
                p, r = 7, 1 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_7(self):
        '''p=2: GGQ Reduced2'''
        PolynomialDegree = 2
        for e in range(1,101):
            points = QuESoPythonModule.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, QuESoPythonModule.IntegrationMethod.GGQ_Reduced2)
            p, r = 2, 0
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            if e > 1:
                p, r = 3, 0 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_8(self):
        '''p=3: GGQ Reduced2'''
        PolynomialDegree = 3
        for e in range(1,101):
            points = QuESoPythonModule.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, QuESoPythonModule.IntegrationMethod.GGQ_Reduced2)
            p, r = 4, 1
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            if e > 1:
                p, r = 5, 1 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)
                p, r = 4, 0 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_9(self):
        '''p=4: GGQ Reduced2'''
        PolynomialDegree = 4
        for e in range(1,101):
            points = QuESoPythonModule.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, QuESoPythonModule.IntegrationMethod.GGQ_Reduced2)
            p, r = 6, 2
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            if e > 1:
                p, r = 7, 2 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)
            if e > 2:
                p, r = 6, 1 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

if __name__ == "__main__":
    unittest.main()
