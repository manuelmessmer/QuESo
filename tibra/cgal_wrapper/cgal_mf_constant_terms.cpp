// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// CGAL includes
#include <CGAL/Polygon_mesh_processing/measure.h>

// Project includes
#include "cgal_wrapper/cgal_mf_constant_terms.h"
#include "utilities/polynomial_utilities.h"

namespace cgal {

void ConstantTerms::Compute(const Element& rElement, VectorType& rConstantTerms, const Parameters& rParam){

    typedef Element::K K;
    typedef K::Point_3 Point_3;
    typedef K::Vector_3 Vector;
    typedef Element::vertex_descriptor vertex_descriptor;
    typedef Element::PositionType PositionType;

    const double jacobian_x = std::abs(rParam.PointB()[0] - rParam.PointA()[0]);
    const double jacobian_y = std::abs(rParam.PointB()[1] - rParam.PointA()[1]);
    const double jacobian_z = std::abs(rParam.PointB()[2] - rParam.PointA()[2]);

    auto& r_surface_mesh = rElement.GetSurfaceMesh();
    PointType a = rElement.GetLocalLowerPoint();
    PointType b = rElement.GetLocalUpperPoint();

    const int ffactor = 1;
    const int order_u = rParam.Order()[0];
    const int order_v = rParam.Order()[1];
    const int order_w = rParam.Order()[2];

    const std::size_t number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor+1) * (order_w*ffactor + 1);

    rConstantTerms.resize(number_of_functions, false);
    std::fill( rConstantTerms.begin(),rConstantTerms.end(), 0.0);

    std::size_t row_index = 0;
    for( int i_x = 0; i_x <= order_u*ffactor; ++i_x){
        for( int i_y = 0; i_y <= order_v*ffactor; ++i_y ){
            for( int i_z = 0; i_z <= order_w*ffactor; ++i_z){
                // Loop over all faces/triangels in r_surface_mesg
                for(auto fd : faces(r_surface_mesh)) {
                    // Compute normal for each face/triangle.
                    Vector normal = CGAL::Polygon_mesh_processing::compute_face_normal(fd, r_surface_mesh);
                    const double area = CGAL::Polygon_mesh_processing::face_area(fd, r_surface_mesh);

                    std::array<PointType, 3> coordinates;
                    PositionType positions = r_surface_mesh.points();
                    int index = 0;
                    for( auto vh : vertices_around_face(halfedge(fd, r_surface_mesh), r_surface_mesh) ){
                        coordinates[index] = {positions[vh].x(), positions[vh].y(), positions[vh].z()};
                        index++;
                    }

                    // Construct triangle to get integration points
                    Triangle3D3N triangle(coordinates[0], coordinates[1], coordinates[2]);
                    auto p_global_integration_points = triangle.GetIntegrationPointsGlobal(1);

                    for( auto global_point : (*p_global_integration_points)){
                        PointType local_point = MappingUtilities::FromGlobalToLocalSpace(global_point, rParam.PointA(), rParam.PointB());
                        PointType value;

                        const double f_x_x = Polynomial::f_x(local_point[0], i_x, a[0], b[0]);
                        const double f_x_y = Polynomial::f_x(local_point[1], i_y, a[1], b[1]);
                        const double f_x_z = Polynomial::f_x(local_point[2], i_z, a[2], b[2]);

                        value[0] = Polynomial::f_x_int(local_point[0], i_x, a[0], b[0])*f_x_y*f_x_z;
                        value[1] = f_x_x*Polynomial::f_x_int(local_point[1], i_y, a[1], b[1])*f_x_z;
                        value[2] = f_x_x*f_x_y*Polynomial::f_x_int(local_point[2], i_z, a[2], b[2]);

                        double integrand = normal[0]*value[0]*jacobian_x + normal[1]*value[1]*jacobian_y + normal[2]*value[2]*jacobian_z;

                        // Normalize weights to 1
                        const double factor = 2.0;
                        // Add contribution to constant terms
                        rConstantTerms[row_index] += 1.0/3.0*integrand * area * global_point.GetWeight() * factor;
                    }
                }
                row_index++;
            }
        }
    }
}

} // End namesapce cgal