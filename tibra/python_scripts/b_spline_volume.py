import numpy as np
import scipy.interpolate as si

class BSplineVolume:
    def __init__(self, Order, NumElements, LowerBound, UpperBound):
        self.spline_u = self.__construct_b_spline(Order[0], NumElements[0], LowerBound[0], UpperBound[0])
        self.spline_v = self.__construct_b_spline(Order[1], NumElements[1], LowerBound[1], UpperBound[1])
        self.spline_w = self.__construct_b_spline(Order[2], NumElements[2], LowerBound[2], UpperBound[2])

    def ControlPoints(self):
        cps = []
        for z in self.spline_w.c:
            for y in self.spline_v.c:
                for x in self.spline_u.c:
                    cps.append( [x, y, z] )
        return cps

    def KnotsU(self):
        return self.spline_u.t

    def KnotsV(self):
        return self.spline_u.t

    def KnotsW(self):
        return self.spline_u.t

    def __construct_b_spline(self, Order, NumElements, LowerBound, UpperBound):
        delta_u = 1.0/NumElements
        knots_u = np.array( (Order+1)*[0] )
        knots_u = np.append( knots_u, (Order+1)*[1] )

        delta_x = (UpperBound-LowerBound) / Order

        cps_x = np.arange(LowerBound, UpperBound+0.5*delta_x,  delta_x )
        spline_u = si.BSpline(knots_u, cps_x, Order)
        knots_u_to_insert = np.arange(delta_u, 1.0-0.5*delta_u, delta_u)

        for knot in knots_u_to_insert:
            spline_u = si.insert(knot, spline_u)

        num_cps = len(spline_u.t) - Order - 1
        spline_u.c = spline_u.c[:num_cps]

        return spline_u
