import QuESo_PythonApplication as QuESo_App
import numpy as np
import scipy.interpolate as si

class BSplineVolume:
    def __init__(self, settings, knot_vector_type):
        ''' Constructor of BSplineVolume \n
        Takes QuESo::Parameters as input.
        '''
        grid_settings = settings[QuESo_App.MainSettings.background_grid_settings]
        self.Order = grid_settings.GetIntVector(QuESo_App.BackgroundGridSettings.polynomial_order)
        NumElements = grid_settings.GetIntVector(QuESo_App.BackgroundGridSettings.number_of_elements)
        LowerBoundXYZ = grid_settings.GetIntVector(QuESo_App.BackgroundGridSettings.lower_bound_xyz)
        UpperBoundXYZ = grid_settings.GetIntVector(QuESo_App.BackgroundGridSettings.upper_bound_xyz)
        LowerBoundUVW = grid_settings.GetIntVector(QuESo_App.BackgroundGridSettings.lower_bound_uvw)
        UpperBoundUVW = grid_settings.GetIntVector(QuESo_App.BackgroundGridSettings.upper_bound_uvw)
        if( knot_vector_type == "open_knot_vector" ):
            open_knot_vector = True
        elif( knot_vector_type == "non_open_knot_vector" ):
            open_knot_vector = False
        else:
            message = "BSplineVolume :: __init__ :: Given 'knot_vector_type': '" + knot_vector_type
            message += "' not valid. Available options are: 'open_knot_vector' and 'non_open_knot_vector'."
            raise Exception(message)
        self.spline_u = self.__construct_b_spline(
           self.Order[0], NumElements[0], LowerBoundXYZ[0], UpperBoundXYZ[0], LowerBoundUVW[0], UpperBoundUVW[0], open_knot_vector)
        self.spline_v = self.__construct_b_spline(
            self.Order[1], NumElements[1], LowerBoundXYZ[1], UpperBoundXYZ[1], LowerBoundUVW[1], UpperBoundUVW[1], open_knot_vector)
        self.spline_w = self.__construct_b_spline(
            self.Order[2], NumElements[2], LowerBoundXYZ[2], UpperBoundXYZ[2], LowerBoundUVW[2], UpperBoundUVW[2], open_knot_vector)

    def GetSpline(self, Index):
        if( Index == 0 ):
            return self.spline_u
        elif( Index == 1 ):
            return self.spline_v
        elif( Index == 2 ):
            return self.spline_w
        else:
            raise Exception("BSplineVolume :: GetSpline :: Index out of scope.")

    def ControlPoints(self):
        ''' Returns control points of B-Spline volume in a list.
            The point indices are linearized and can be accessed the following:

            cps = self.ControlPoints() \n
            count = 0 \n
            for i_w in range(NumberControlPointsInW()):
                for i_v in range(NumberControlPointsInV()):
                    for i_u in range(NumberControlPointsInU()):
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

    def __construct_b_spline(self, Order, NumElements, LowerBoundX, UpperBoundX, LowerBoundU, UpperBoundU, OpenKnotVector=False):
        delta_u = (UpperBoundU-LowerBoundU)/NumElements

        if OpenKnotVector:
            knots_u = np.array( (Order+1)*[LowerBoundU] )
            knots_u = np.append( knots_u, (Order+1)*[UpperBoundU] )
        else:
            knots_u = np.array( [LowerBoundU - (Order-i)*delta_u for i in range(Order+1)] )
            knots_u = np.append( knots_u, [(UpperBoundU + i*delta_u) for i in range(Order+1)] )

        delta_x = (UpperBoundX-LowerBoundX) / Order
        center = (UpperBoundX+LowerBoundX) / 2.0
        cps_x = np.arange(LowerBoundX, UpperBoundX+0.5*delta_x,  delta_x )

        if not OpenKnotVector:
            cps_x = [ val - (Order-1)*(center-val) / (NumElements) for val in cps_x  ]

        spline_u = si.BSpline(knots_u, cps_x, Order, extrapolate=False)
        knots_u_to_insert = np.arange(LowerBoundU+delta_u, UpperBoundU-0.5*delta_u, delta_u)

        for knot in knots_u_to_insert:
            spline_u = si.insert(knot, spline_u)

        num_cps = len(spline_u.t) - Order - 1
        spline_u.c = spline_u.c[:num_cps]

        return spline_u

