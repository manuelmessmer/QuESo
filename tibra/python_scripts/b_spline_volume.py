import numpy as np
import scipy.interpolate as si

class BSplineVolume:
    def __init__(self, Parameters):
        ''' Constructor of BSplineVolume \n
        Takes TIBRA::Parameters as input.
        '''
        self.Order = Parameters.Order()
        NumElements = Parameters.NumberOfElements()
        LowerBound = Parameters.LowerBound()
        UpperBound = Parameters.UpperBound()
        self.spline_u = self.__construct_b_spline(self.Order[0], NumElements[0], LowerBound[0], UpperBound[0])
        self.spline_v = self.__construct_b_spline(self.Order[1], NumElements[1], LowerBound[1], UpperBound[1])
        self.spline_w = self.__construct_b_spline(self.Order[2], NumElements[2], LowerBound[2], UpperBound[2])

    def ControlPoints(self):
        ''' Returns control points of B-Spline volume in a list.
            The point indices are linearized and can be accessed the following:

            cps = self.ControlPoints() \n
            count = 0 \n
            for i_u in range(NumberControlPointsInU()):
                for i_v in range(NumberControlPointsInV()):
                    for w in range(NumberControlPointsInW()):
                        point = cps[count] \n
                        count += 1
        '''
        cps = []
        for z in self.spline_w.c:
            for y in self.spline_v.c:
                for x in self.spline_u.c:
                    cps.append( [x, y, z] )
        return cps

    def ControlPointsMatrix(self):
        ''' Returns control points of B-Spline volume in a matrix.
            Points can be accessed as:

            cps = self.ControlPoints() \n
            for i_w in range(NumberControlPointsInW()):
                for i_v in range(NumberControlPointsInV()):
                    for i_u in range(NumberControlPointsInU()):
                        point = cps[i_u, i_v, i_w]
        '''
        n_cps_u = self.NumberControlPointsInU()
        n_cps_v = self.NumberControlPointsInV()
        n_cps_w = self.NumberControlPointsInW()
        n_cps = n_cps_u*n_cps_v*n_cps_w

        cps = np.zeros(n_cps*3).reshape(n_cps_u, n_cps_v, n_cps_w, 3)
        for i_z, z in enumerate(self.spline_w.c):
            for i_y, y in enumerate(self.spline_v.c):
                for i_x, x in enumerate(self.spline_u.c):
                    cps[i_x, i_y, i_z] = [x, y, z]

        return cps

    def KnotsU(self):
        ''' KnotsU \n
        Returns knot vector along u- (x) direction.
        '''
        return self.spline_u.t.tolist()

    def KnotsV(self):
        ''' KnotsV \n
        Returns knot vector along v- (y) direction.
        '''
        return self.spline_v.t.tolist()

    def KnotsW(self):
        ''' KnotsW \n
        Returns knot vector along w- (z) direction.
        '''
        return self.spline_w.t.tolist()

    def PolynomialOrder(self):
        ''' PolynomialOrder \n
        Returns polynomial order of b-spline volume: [p_x, p_y, p_z]
        '''
        return self.Order

    def NumberControlPointsInU(self):
        ''' NumberControlPointsInU \n
        Returns number of control points in u- (x) direction.
        '''
        return len(self.spline_u.c)

    def NumberControlPointsInV(self):
        ''' NumberControlPointsInV \n
        Returns number of control points in v- (y) direction.
        '''
        return len(self.spline_v.c)

    def NumberControlPointsInW(self):
        ''' NumberControlPointsInW \n
        Returns number of control points in w- (z) direction.
        '''
        return len(self.spline_w.c)

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

