// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// External includes
#include <cmath>
#include <algorithm>
#include <iomanip>      // std::setprecision
#include <limits>
#include <chrono>

// Project includes
#include "utilities/geometrical_entity_classifier.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> drandon(0, 1);

bool GeometricalEntityClassifier::IsInside(const PointType& rPoint){

    if( mTree.IsWithinBoundingBox(rPoint)) {
        bool is_on_boundary = true;
        int intersection_count = 0;
        while( is_on_boundary ){

            // Get random direction
            TriangleMesh::Vector3d direction = {drandon(gen), drandon(gen), drandon(gen)};
            // Normalize
            double sum_direction = std::sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
            std::for_each(direction.begin(), direction.end(), [sum_direction](auto& rValue) { rValue /= sum_direction;} );

            // Construct ray
            Ray_AABB_primitive ray(rPoint, direction);
            // Get potential ray intersections from AABB tree.
            auto potential_intersections = mTree.Query(ray);

            // Test if potential intersections actually intersect.
            // If intersection lies on boundary cast a new ray.
            intersection_count = 0;
            is_on_boundary = false;
            for( auto r : potential_intersections){
                const auto& p1 = mTriangleMesh.P1(r);
                const auto& p2 = mTriangleMesh.P2(r);
                const auto& p3 = mTriangleMesh.P3(r);
                double t, u, v;
                if( ray.intersect(p1, p2, p3, t, u, v) ) {
                    intersection_count++;
                    double sum_u_v = u+v;
                    if( u < 0.0+1e-14 || v < 0.0+1e-14 || sum_u_v > 1.0-1e-14 ){
                        is_on_boundary = true;
                        break;
                    }
                }
            }
        }
        if( intersection_count % 2 == 1){ // If intersection count is odd, return true.
            return true;
        }
    }

    return false;

}

// Winding numbers algorithm: It actually works!!!

// std::cout << "Done" << std::endl;
// std::cout << "Tree timer: " << tree_timer << std::endl;
// std::cout << "ray timer: " << ray_timer << std::endl;
// return std::make_unique<std::vector<bool>>(rr);
// for( int triangle_id = 0; triangle_id < rTriangleMesh.NumOfTriangles(); ++triangle_id ){
//     const auto& p1 = rTriangleMesh.P1(triangle_id);
//     const auto& p2 = rTriangleMesh.P2(triangle_id);
//     const auto& p3 = rTriangleMesh.P3(triangle_id);

//     IndexType count = 0;
//     for( const auto& point : rPoints){
//     //const auto point = rPoints[0];
//         PointType a{p1[0] - point[0], p1[1] - point[1], p1[2] - point[2]};
//         PointType b{p2[0] - point[0], p2[1] - point[1], p2[2] - point[2]};
//         PointType c{p3[0] - point[0], p3[1] - point[1], p3[2] - point[2]};

//         // Compute determinant of:
//         // [a0 b0 c0]
//         // [a1 b1 c1]
//         // [a2 b2 c2]
//         const long double omega =  a[0] * (b[1]*c[2] - b[2]*c[1])
//                                   -b[0] * (a[1]*c[2] - a[2]*c[1])
//                                   +c[0] * (a[1]*b[2] - a[2]*b[1]);


//         const long double anorm = std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
//         const long double bnorm = std::sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
//         const long double cnorm = std::sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

//         long double k = anorm*bnorm*cnorm;

//         k += cnorm * ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
//         k += anorm * ( b[0]*c[0] + b[1]*c[1] + b[2]*c[2]);
//         k += bnorm * ( a[0]*c[0] + a[1]*c[1] + a[2]*c[2]);

//         ret[count] += std::atan2(omega,k);
//         count++;
//     }

// }

// std::unique_ptr<std::vector<bool>> p_result = std::make_unique<std::vector<bool>>(rPoints.size(), false);
// auto& r_result = *p_result;
// for( int i = 0; i < ret.size(); ++i ){
//     if( ret[i] >= (2.0*My_PI-1e-10) ) // Due to limited precision of double small treshhold (1e-10) is required.
//         r_result[i] = true;
// }

// return std::move(p_result);

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