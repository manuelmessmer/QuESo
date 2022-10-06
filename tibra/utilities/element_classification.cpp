// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// External includes
#include <cmath>
#include <algorithm>
#include <iomanip>      // std::setprecision
#include <limits>

// Project includes
#include "utilities/element_classification.h"
#include <boost/multiprecision/cpp_dec_float.hpp>

# define My_PI          3.141592653589793238462643383279502884L /* pi */

float atan2_approx(float y, float x) {
	float abs_y = std::fabs(y) + 1e-10f;      // kludge to prevent 0/0 condition
	float r = (x - std::copysign(abs_y, x)) / (abs_y + std::fabs(x));
	float angle = M_PI/2.f - std::copysign(M_PI/4.f, x);

	angle += (0.1963f * r * r - 0.9817f) * r;
	return std::copysign(angle, y);
}

std::unique_ptr<std::vector<bool>> ElementClassification::PointsAreInside(const TriangleMesh& rTriangleMesh, const std::vector<PointType>& rPoints){

    typedef std::size_t IndexType;
    typedef boost::multiprecision::cpp_dec_float_50 Float;

    std::vector<long double> ret(rPoints.size(), 0.0);

    for( int triangle_id = 0; triangle_id < rTriangleMesh.NumOfTriangles(); ++triangle_id ){
        const auto& p1 = rTriangleMesh.P1(triangle_id);
        const auto& p2 = rTriangleMesh.P2(triangle_id);
        const auto& p3 = rTriangleMesh.P3(triangle_id);

        IndexType count = 0;
        for( const auto& point : rPoints){
        //const auto point = rPoints[0];
            PointType a{p1[0] - point[0], p1[1] - point[1], p1[2] - point[2]};
            PointType b{p2[0] - point[0], p2[1] - point[1], p2[2] - point[2]};
            PointType c{p3[0] - point[0], p3[1] - point[1], p3[2] - point[2]};

            // Compute determinant of:
            // [a0 b0 c0]
            // [a1 b1 c1]
            // [a2 b2 c2]
            const long double omega =  a[0] * (b[1]*c[2] - b[2]*c[1])
                                      -b[0] * (a[1]*c[2] - a[2]*c[1])
                                      +c[0] * (a[1]*b[2] - a[2]*b[1]);


            const long double anorm = std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
            const long double bnorm = std::sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
            const long double cnorm = std::sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

            long double k = anorm*bnorm*cnorm;

            k += cnorm * ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
            k += anorm * ( b[0]*c[0] + b[1]*c[1] + b[2]*c[2]);
            k += bnorm * ( a[0]*c[0] + a[1]*c[1] + a[2]*c[2]);

            ret[count] += std::atan2(omega,k);
            count++;
        }

    }

    std::unique_ptr<std::vector<bool>> p_result = std::make_unique<std::vector<bool>>(rPoints.size(), false);
    auto& r_result = *p_result;
    for( int i = 0; i < ret.size(); ++i ){
        if( ret[i] >= (2.0*My_PI-1e-10) ) // Due to limited precision of double small treshhold (1e-10) is required.
            r_result[i] = true;
    }

    return std::move(p_result);
}

bool ElementClassification::Anorm2(std::vector<double>& rResult, const std::vector<PointType>& rX){
    rResult.reserve(rX.size());
    std::for_each(rX.begin(), rX.end(), [&rResult](auto& point){
        rResult.push_back(std::sqrt( point[0]*point[0]+point[1]*point[1]+point[2]*point[2]) ); }
    );
    return true;
}
// def is_inside_turbo(triangles, X):
// 	# Compute euclidean norm along axis 1
// 	def anorm2(X):
// 		return numpy.sqrt(numpy.sum(X ** 2, axis = 1))



// 	# Compute 3x3 determinant along axis 1
// 	def adet(X, Y, Z):
// 		ret  = numpy.multiply(numpy.multiply(X[:,0], Y[:,1]), Z[:,2])
// 		ret += numpy.multiply(numpy.multiply(Y[:,0], Z[:,1]), X[:,2])
// 		ret += numpy.multiply(numpy.multiply(Z[:,0], X[:,1]), Y[:,2])
// 		ret -= numpy.multiply(numpy.multiply(Z[:,0], Y[:,1]), X[:,2])
// 		ret -= numpy.multiply(numpy.multiply(Y[:,0], X[:,1]), Z[:,2])
// 		ret -= numpy.multiply(numpy.multiply(X[:,0], Z[:,1]), Y[:,2])
// 		return ret



// 	# One generalized winding number per input vertex
// 	ret = numpy.zeros(X.shape[0], dtype = X.dtype)

// 	# Accumulate generalized winding number for each triangle
// 	for U, V, W in triangles:
// 		A, B, C = U - X, V - X, W - X
// 		omega = adet(A, B, C)

// 		a, b, c = anorm2(A), anorm2(B), anorm2(C)
// 		k  = a * b * c
// 		k += c * numpy.sum(numpy.multiply(A, B), axis = 1)
// 		k += a * numpy.sum(numpy.multiply(B, C), axis = 1)
// 		k += b * numpy.sum(numpy.multiply(C, A), axis = 1)

// 		ret += numpy.arctan2(omega, k)

// 	# Job done
// 	return ret >= 2 * numpy.pi