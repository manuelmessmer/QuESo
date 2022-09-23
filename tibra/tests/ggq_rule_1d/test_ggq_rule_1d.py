# Project imports
import TIBRA_PythonApplication

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

        cps = []
        weights = []
        for point in points:
            cps.append(point[0])
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
            B = BSpline.design_matrix(cps, knots, p).toarray()
            error = np.linalg.norm( (target - np.array(weights).dot(np.array(B)) )/dim )
            if must_pass:
                self.assertLess( error, 1e-15 )
                self.assertEqual(len(cps), m)
            else:
                self.assertGreater( error, 1e-10)
        else:
            print("Testing :: QQGRule1d is skipped. Requires scipy version: >=" + np_minversion)

    def test_1(self):
        '''p=2: GGQ Optimal'''
        PolynomialDegree = 2
        for e in range(1,101):
            points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedExact)
            p, r = 4, 0
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            if e > 1:
                p, r = 5, 0
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_2(self):
        '''p=3: GGQ Optimal'''
        PolynomialDegree = 3
        for e in range(1,101):
            points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedExact)
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
            points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedExact)
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
            points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder1)
            p, r = 3, 0
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            p, r = 4, 0 # Expected to fail.
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_5(self):
        '''p=3: GGQ Reduced1'''
        PolynomialDegree = 3
        for e in range(1,101):
            points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder1)
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
            points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder1)
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
            points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder2)
            p, r = 2, 0
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            if e > 1:
                p, r = 3, 0 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    def test_8(self):
        '''p=3: GGQ Reduced2'''
        PolynomialDegree = 3
        for e in range(1,101):
            points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder2)
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
            points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder2)
            p, r = 6, 2
            self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)
            if e > 1:
                p, r = 7, 2 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)
            if e > 2:
                p, r = 6, 1 # Expected to fail.
                self.check_ggq_rules(points, p, r, e, 0.0, 1.0, False)

    # def test_6(self):
    #     '''p=2: GGQ Reduced2'''
    #     PolynomialDegree = 2
    #     for e in range(24,29):
    #         points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder2)

    #         p, r = 2, 0
    #         self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)

    # def test_6(self):
    #     '''p=3: GGQ Reduced2'''
    #     PolynomialDegree = 3
    #     for e in range(1,101):
    #         points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder2)
    #         print("e", e, len(points)%2)
    #         #print(len(points))
    #         # weight = 0.0
    #         # for i, p in enumerate(points):
    #         #     print(i+1, p[0], p[1])
    #         # #     weight += p[1]
    #         #print("Weight: ", weight)
    #         p, r = 4, 1
    #         self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)

    # def test_7(self):
    #     '''p=4: GGQ Reduced2'''
    #     PolynomialDegree = 4
    #     for e in range(1,101):
    #         points = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder2)
    #         print("e", e, len(points)%2)
    #         #print(len(points))
    #         # weight = 0.0
    #         # for i, p in enumerate(points):
    #         #     print(i+1, p[0], p[1])
    #         # #     weight += p[1]
    #         #print("Weight: ", weight)
    #         p, r = 6, 2
    #         self.check_ggq_rules(points, p, r, e, 0.0, 1.0, True)

    def f2(self):
        '''p=4: GGQ Reduced1'''
        PolynomialDegree = 2
        e1 = 37

        points1 = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e1, TIBRA_PythonApplication.IntegrationMethod.ReducedExact)
        #points2 = TIBRA_PythonApplication.IntegrationPointFactory1D.GetGGQ(PolynomialDegree, e2, TIBRA_PythonApplication.IntegrationMethod.ReducedOrder1)
        # i = 1
        # for points37 in points1:
        #     #print(i, "%.15f" % (points37[0]*37), "%.15f" % (points38[0]*38))
        #     print(i, "%.16f" % (points37[0]*37), "%.16f" % (points37[1]*37))
        #     i += 1

        # return 0

        # for points in points_S_4_1_38:
        #     #print(i, "%.15f" % (points37[0]*37), "%.15f" % (points38[0]*38))
        #     print(i, "%.16f" % (points[0]*38), "%.16f" % (points[1]*38))
        #     #print(i, "%.16f" % (points37[0]*37), "%.16f" % (points37[1]*37))
        #     i += 1
        #     #print("%.16f" % (points38[0]*38), "%.16f" % (points38[1]*38))
        # i = 1
        # for points in points_S_4_1_37:
        #     #print(i, "%.15f" % (points37[0]*37), "%.15f" % (points38[0]*38))
        #     print(i, "%.16f" % (points[0]*37), "%.16f" % (points[1]*37))
        #     i += 1

        i = 1
        print("len 31: ", len(points_S_6_2_20))
        for points in points_S_6_2_20:
            #print(i, "%.15f" % (points37[0]*37), "%.15f" % (points38[0]*38))
            print(i, "%.16f" % (points[0]*31), "%.16f" % (points[1]*31))
            i += 1


if __name__ == "__main__":
    unittest.main()

# for val1, val2 in zip(points1, points2):
        #     print("%.16f" % (val1[0]*e1), ", ", val2[0]*e2)
        #     #print("%.16f" % (val2[0]*e2), ", ", val2[1]*e2)

        # p21_4_1 = [
        #         [ 0.0080575753327455, 0.0195878352740383 ],
        #         [ 0.0340076706183558, 0.0288018925883000 ],
        #         [ 0.0638223881196816, 0.0314977459489792 ],
        #         [ 0.0952768913212719, 0.0308270251766042 ],
        #         [ 0.1269842461096242, 0.0321426846934533 ],
        #         [ 0.1587301757402807, 0.0321428226974250 ],
        #         [ 0.1904761924388986, 0.0309523745735811 ],
        #         [ 0.2222222222222225, 0.0321428571428569 ],
        #         [ 0.2539682539682541, 0.0321428571428571 ],
        #         [ 0.2857142857046493, 0.0309523809210623 ],
        #         [ 0.3174603173768016, 0.0321428569737375 ],
        #         [ 0.3492063486217188, 0.0321428562971816 ],
        #         [ 0.3809496707121625, 0.0309435762434679 ],
        #         [ 0.4126749532421955, 0.0320954744050789 ],
        #         [ 0.4442785595514147, 0.0318974900375299 ],
        #         [ 0.4745800954348970, 0.0268888571854336 ],
        #         [ 0.5000000000000000, 0.0253968253968260 ],
        #         [ 0.5254199045651029, 0.0268888571854336 ],
        #         [ 0.5557214404485854, 0.0318974900375299 ],
        #         [ 0.5873250467578045, 0.0320954744050789 ],
        #         [ 0.6190503292878377, 0.0309435762434679 ],
        #         [ 0.6507936513782812, 0.0321428562971816 ],
        #         [ 0.6825396826231985, 0.0321428569737375 ],
        #         [ 0.7142857142953507, 0.0309523809210623 ],
        #         [ 0.7460317460317459, 0.0321428571428571 ],
        #         [ 0.7777777777777776, 0.0321428571428569 ],
        #         [ 0.8095238075611015, 0.0309523745735811 ],
        #         [ 0.8412698242597194, 0.0321428226974250 ],
        #         [ 0.8730157538903759, 0.0321426846934533 ],
        #         [ 0.9047231086787281, 0.0308270251766042 ],
        #         [ 0.9361776118803184, 0.0314977459489792 ],
        #         [ 0.9659923293816443, 0.0288018925883000 ],
        #         [ 0.9919424246672546, 0.0195878352740383 ]]

        # for i, p in enumerate(p21_4_1):
        #     print(p[0]*21) #, "%.16f" % (p[1]*30))

        # p30 = [[ 0.0111111111111111, 0.0250000000000001 ],
        #         [ 0.0407407407407407, 0.0321428571428571 ],
        #         [ 0.0737373737373737, 0.0332417582417582 ],
        #         [ 0.1070460704607046, 0.0333267248215702 ],
        #         [ 0.1403776325344953, 0.0333328586888421 ],
        #         [ 0.1737108386845690, 0.0333332992544912 ],
        #         [ 0.2070441628864904, 0.0333333308865778 ],
        #         [ 0.2403774955642176, 0.0333333331576642 ],
        #         [ 0.2737108288504805, 0.0333333333207209 ],
        #         [ 0.3070441621804343, 0.0333333333324278 ],
        #         [ 0.3403774955135249, 0.0333333333332683 ],
        #         [ 0.3737108288468409, 0.0333333333333286 ],
        #         [ 0.4070441621801730, 0.0333333333333330 ],
        #         [ 0.4403774955135063, 0.0333333333333333 ],
        #         [ 0.4737108288468396, 0.0333333333333334 ],
        #         [ 0.5000000000000000, 0.0192450089729872 ],
        #         [ 0.5262891711531604, 0.0333333333333334 ],
        #         [ 0.5596225044864938, 0.0333333333333333 ],
        #         [ 0.5929558378198270, 0.0333333333333330 ],
        #         [ 0.6262891711531592, 0.0333333333333286 ],
        #         [ 0.6596225044864751, 0.0333333333332683 ],
        #         [ 0.6929558378195657, 0.0333333333324278 ],
        #         [ 0.7262891711495194, 0.0333333333207209 ],
        #         [ 0.7596225044357825, 0.0333333331576642 ],
        #         [ 0.7929558371135096, 0.0333333308865778 ],
        #         [ 0.8262891613154310, 0.0333332992544912 ],
        #         [ 0.8596223674655047, 0.0333328586888421 ],
        #         [ 0.8929539295392954, 0.0333267248215702 ],
        #         [ 0.9262626262626263, 0.0332417582417582 ],
        #         [ 0.9592592592592593, 0.0321428571428571 ],
        #         [ 0.9888888888888889, 0.0250000000000001 ]]

        # p31 = [
        #         [ 0.0107526881720430, 0.0241935483870970 ],
        #         [ 0.0394265232974910, 0.0311059907834101 ],
        #         [ 0.0713587487781036, 0.0321694434597661 ],
        #         [ 0.1035929714135851, 0.0322516691821647 ],
        #         [ 0.1358493218075761, 0.0322576051827504 ],
        #         [ 0.1681072632431312, 0.0322580315366044 ],
        #         [ 0.2003653189224101, 0.0322580621483012 ],
        #         [ 0.2326233828040815, 0.0322580643461267 ],
        #         [ 0.2648814472746585, 0.0322580645039234 ],
        #         [ 0.2971395117875170, 0.0322580645152527 ],
        #         [ 0.3293975763034114, 0.0322580645160662 ],
        #         [ 0.3616556408195234, 0.0322580645161245 ],
        #         [ 0.3939137053356511, 0.0322580645161287 ],
        #         [ 0.4261717698517802, 0.0322580645161290 ],
        #         [ 0.4584298343679094, 0.0322580645161291 ],
        #         [ 0.4877445857152969, 0.0254411333740261 ],
        #         [ 0.5122554142847031, 0.0254411333740261 ],
        #         [ 0.5415701656320907, 0.0322580645161291 ],
        #         [ 0.5738282301482197, 0.0322580645161290 ],
        #         [ 0.6060862946643488, 0.0322580645161287 ],
        #         [ 0.6383443591804766, 0.0322580645161245 ],
        #         [ 0.6706024236965886, 0.0322580645160662 ],
        #         [ 0.7028604882124830, 0.0322580645152527 ],
        #         [ 0.7351185527253415, 0.0322580645039234 ],
        #         [ 0.7673766171959185, 0.0322580643461267 ],
        #         [ 0.7996346810775899, 0.0322580621483012 ],
        #         [ 0.8318927367568687, 0.0322580315366044 ],
        #         [ 0.8641506781924240, 0.0322576051827504 ],
        #         [ 0.8964070285864149, 0.0322516691821647 ],
        #         [ 0.9286412512218963, 0.0321694434597661 ],
        #         [ 0.9605734767025089, 0.0311059907834101 ],
        #         [ 0.9892473118279570, 0.0241935483870970 ] ]

        # for p30_, p31_ in zip(p30, p31):
        #     #print(i, p_[0]*30, "%.16f" % (p_[1]*30))
        #     #print(i, "%.16f" % (p30_[1]*30))
        #     print(p30_[0]*30, p31_[0]*31)
        #print(cps)