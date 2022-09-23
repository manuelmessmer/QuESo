// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/centroid.h>

// External includes
#include <stdexcept>
#include <cmath>

// Project includes
#include "quadrature/moment_fitting_utilities.h"
#include "utilities/mapping_utilities.h"
#include "utilities/polynomial_utilities.h"
#include "geometries/element.h"
#include "geometries/triangle_3d_3n.h"
#include "solvers/nnls.h"
#include "io/io_utilities.h"

typedef std::size_t SizeType;
typedef std::array<double, 3> PointType;
typedef std::array<int, 3> IntArrayType;
typedef boost::numeric::ublas::matrix<double> MatrixType;
typedef boost::numeric::ublas::vector<double> VectorType;

typedef Element::K K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector;
typedef Element::vertex_descriptor vertex_descriptor;
typedef Element::PositionType PositionType;


void MomentFitting::DistributeInitialIntegrationPoints(const Element& rElement, IntegrationPointVectorType& rIntegrationPoint, const int PointDistributionFactor, const Parameters& rParam){

    const double factor = PointDistributionFactor;
    rIntegrationPoint.reserve( (int) factor*(rParam.Order()[0]+1)*factor*(rParam.Order()[1]+1)*factor*(rParam.Order()[2]+1) );

    auto bounding_box = rElement.ComputeTrimmedBoundingBox();

    const double delta_x = std::abs(bounding_box[1][0] - bounding_box[0][0]) / ( (double)factor*(rParam.Order()[0]+1) );
    const double delta_y = std::abs(bounding_box[1][1] - bounding_box[0][1]) / ( (double)factor*(rParam.Order()[1]+1) );
    const double delta_z = std::abs(bounding_box[1][2] - bounding_box[0][2]) / ( (double)factor*(rParam.Order()[2]+1) );

    double xx = bounding_box[0][0] + 0.5*delta_x;
    while( xx < bounding_box[1][0] ){
        double yy = bounding_box[0][1] + 0.5*delta_y;
        while( yy < bounding_box[1][1] ){
            double zz = bounding_box[0][2] + 0.5*delta_z;
            while( zz < bounding_box[1][2]){
                PointType tmp_point = {xx,yy,zz};
                if( rElement.IsPointInTrimmedDomain(tmp_point) ){
                    rIntegrationPoint.push_back(IntegrationPoint((xx - rParam.PointA()[0]) / std::abs(rParam.PointB()[0] - rParam.PointA()[0]),
                                                    (yy - rParam.PointA()[1]) / std::abs(rParam.PointB()[1] - rParam.PointA()[1]),
                                                    (zz - rParam.PointA()[2]) / std::abs(rParam.PointB()[2] - rParam.PointA()[2]),
                                                     0.0 ));
                }
                zz += delta_z;
            }
            yy += delta_y;
        }
        xx += delta_x;
    }
}

void MomentFitting::ComputeConstantTerms(const Element& rElement, VectorType& rConstantTerms, const Parameters& rParam){

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

void MomentFitting::CreateIntegrationPointsTrimmed(Element& rElement, const Parameters& rParam){

    double residual = 1e10;
    int point_distribution_factor = rParam.GetPointDistributionFactor();
    VectorType constant_terms{};
    ComputeConstantTerms(rElement, constant_terms, rParam);
    const int max_iteration = 3;
    int iteration = 1;
    while( residual > rParam.MomentFittingResidual() && iteration < max_iteration){
        auto& reduced_points = rElement.GetIntegrationPointsTrimmed();
        if( !rParam.UseCustomizedTrimmedPositions() ){
            // This is only used for test_moment_fitting.cpp
            reduced_points.clear();
        }
        residual = CreateIntegrationPointsTrimmed(rElement, constant_terms, point_distribution_factor, rParam);
        point_distribution_factor *= 2;
        if( residual > 1e-5 ) {
            //std::cout << "size: " << reduced_points.size() << std::endl;
            reduced_points.clear();
        }
        iteration++;
    }

    if( residual > rParam.MomentFittingResidual() && rParam.EchoLevel() > 2){
        std::cout << "size: " << rElement.GetIntegrationPointsTrimmed().size() << std::endl;
        std::cout << "Moment Fitting :: Targeted residual can not be achieved!: " << residual << std::endl;
        // IO::WriteMeshToVTK(rElement.GetSurfaceMesh(), "fail.vtk", true);
    }
}

double MomentFitting::CreateIntegrationPointsTrimmed(Element& rElement, const VectorType& rConstantTerms, const int PointDistributionFactor, const Parameters& rParam) {

    IntegrationPointVectorType new_integration_points;
    int maximum_iteration;
    if( !rParam.UseCustomizedTrimmedPositions() ){
        DistributeInitialIntegrationPoints( rElement, new_integration_points, PointDistributionFactor, rParam);
        maximum_iteration = 1000;
    }
    else { // This is only used for test_moment_fitting.cpp
        new_integration_points = rElement.GetIntegrationPointsTrimmed();
        rElement.GetIntegrationPointsTrimmed().clear();
        maximum_iteration = 1;
    }

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

    double global_residual = 1e-15;
    const double allowed_residual = rParam.MomentFittingResidual();
    int number_iterations = 0;
    bool change = false;
    IntegrationPointVectorType prev_solution{};
    double prev_residual = 0.0;
    while( change || (global_residual < allowed_residual && number_iterations < maximum_iteration) ){

        const std::size_t number_reduced_points = new_integration_points.size();

        MatrixType fitting_matrix(number_of_functions, number_reduced_points);

        std::size_t row_index = 0;
        for( int i_x = 0; i_x <= order_u*ffactor; ++i_x){
            for( int i_y = 0; i_y <= order_v*ffactor; ++i_y ){
                for( int i_z = 0; i_z <= order_w*ffactor; ++i_z){
                    // Loop over all faces/triangels in r_surface_mesg
                    const auto points_it_begin = new_integration_points.begin();
                    for( int column_index = 0; column_index < number_reduced_points; ++column_index ){
                        auto point_it = points_it_begin + column_index;

                        const double value = Polynomial::f_x(point_it->X(), i_x, a[0], b[0])
                                                * Polynomial::f_x(point_it->Y(), i_y, a[1], b[1])
                                                * Polynomial::f_x(point_it->Z(), i_z, a[2], b[2]);

                        fitting_matrix(row_index,column_index) = value;
                    }
                    row_index++;
                }
            }
        }
        VectorType weights(number_reduced_points);
        // Solve non-negative Least-Square-Error problem.
        global_residual = NNLS::nnls(fitting_matrix, rConstantTerms, weights)/number_of_functions;

        //Write computed weights onto reduced integration points
        for( int i = 0; i < number_reduced_points; ++i){
            // Divide by det_jacobian to account for the corresponding multiplication during the element integration within the used external solver.
            double new_weight = weights[i]/(jacobian_x*jacobian_y*jacobian_z);
            new_integration_points[i].SetWeight(new_weight);
        }
        change = false;
        if( !rParam.UseCustomizedTrimmedPositions() ){
            if( number_iterations == 0){
                if( global_residual > allowed_residual ){
                    //std::cout << "Moment Fitting :: Targeted residual can not be achieved!: " << global_residual << std::endl;
                }
                // Sort integration points according to weight
                sort(new_integration_points.begin(), new_integration_points.end(), [](const IntegrationPoint& point_a, const IntegrationPoint& point_b) -> bool {
                        return point_a.GetWeightConst() > point_b.GetWeightConst();
                    });
                new_integration_points.erase(new_integration_points.begin()+number_of_functions, new_integration_points.end());
                change = true;
            }
            else if( global_residual < allowed_residual ){
                prev_solution.clear();
                prev_solution.insert(prev_solution.begin(), new_integration_points.begin(), new_integration_points.end());
                prev_residual = global_residual;
                auto min_value_it = new_integration_points.begin();
                double min_value = 1e10;
                double max_value = -1e10;
                auto begin_it = new_integration_points.begin();
                for(int i = 0; i < new_integration_points.size(); i++){
                    auto it = begin_it + i;
                    if( it->GetWeight() < min_value ) {
                        min_value_it = it;
                        min_value = it->GetWeight();
                    }
                    if( it->GetWeight() > max_value ) {
                        max_value = it->GetWeight();
                    }
                }
                begin_it = new_integration_points.begin();
                int counter = 0;
                for(int i = 0; i < new_integration_points.size(); i++){
                    auto it = begin_it + i;
                    // TODO: Fix this > 2..4
                    if( it->GetWeight() < 1e-8*max_value && new_integration_points.size() > 4){
                        new_integration_points.erase(it);
                        change = true;
                        counter++;
                    }
                }
                if( counter == 0 && new_integration_points.size() > 4){
                    new_integration_points.erase(min_value_it);
                    change = true;
                }
                if( new_integration_points.size() == 4){
                    number_iterations = maximum_iteration + 1;
                }
            }
        }
        number_iterations++;
    }
    auto& reduced_points = rElement.GetIntegrationPointsTrimmed();

    if( (global_residual >= allowed_residual && prev_solution.size() > 0 && number_iterations < maximum_iteration) ) {
        reduced_points.insert(reduced_points.begin(), prev_solution.begin(), prev_solution.end());
        reduced_points.erase(std::remove_if(reduced_points.begin(), reduced_points.end(), [](const IntegrationPoint& point) {
            return point.GetWeightConst() < 1e-14; }), reduced_points.end());
        return prev_residual;
    }
    else{
        reduced_points.insert(reduced_points.begin(), new_integration_points.begin(), new_integration_points.end());
        reduced_points.erase(std::remove_if(reduced_points.begin(), reduced_points.end(), [](const IntegrationPoint& point) {
            return point.GetWeightConst() < 1e-14; }), reduced_points.end());
        return global_residual;
    }

}


// double MomentFitting::f_x(double x, int order, double a, double b){
//     double tmp_x = (2*x - a - b)/(b-a);
//     return p_n(tmp_x,order);
// }



// double MomentFitting::f_x(double x, int order){
//     if( order == 0){double
//     }
//     else {
//         return std::pow(x, order);
//     }
// }

// double MomentFitting::f_x_integral(double x, int order) {
//     if(  order == 0 ){
//         return x;
//     }
//     else {
//         return 1.0/ ( (double)order + 1.0) * std::pow(x, order+1);
//     }
// }

// double MomentFitting::f_x_integral(double x, int order, double a, double b){
//     switch(order)
//     {
//         case 0:
//             return x;
//         case 1:
//             return -std::pow((a + b - 2.0*x),2)/(4.0*(a - b));
//         case 2:
//             return - x/2.0 - std::pow((a + b - 2.0*x),3)/(4.0*std::pow((a - b),2));
//         case 3:
//             return (3.0*std::pow( (a + b - 2.0*x),2) )/(8.0*(a - b)) - (5*std::pow((a + b - 2.0*x),4))/(16*std::pow((a - b),3));
//         case 4:
//             return (3.0*x)/8.0 + (5.0*std::pow((a + b - 2.0*x),3))/(8*std::pow((a - b),2)) - (7.0*std::pow((a + b - 2*x),5))/(16*std::pow((a - b),4));
//         case 5:
//             return (35*std::pow((a + b - 2*x),4))/(32*std::pow((a - b),3)) - (15*std::pow((a + b - 2*x),2))/(32*(a - b)) - (21*std::pow((a + b - 2*x),6))/(32*std::pow((a - b),5));
//         case 6:
//             return (63*std::pow((a + b - 2*x),5))/(32*std::pow((a - b),4)) - (35*std::pow((a + b - 2*x),3))/(32*std::pow((a - b),2)) - (5*x)/16 - (33*std::pow((a + b - 2*x),7))/(32*std::pow((a - b),6));
//         case 7:
//             return (35.0*std::pow( (a + b - 2*x),2))/(64*(a - b)) - (315*std::pow((a + b - 2*x),4))/(128.0*std::pow((a - b),3)) + (231.0*std::pow((a + b - 2*x),6))/(64.0*std::pow((a - b),5)) - (429.0*std::pow((a + b - 2*x),8))/(256.0*std::pow((a - b),7));
//         case 8:
//             return (35.0*x)/128.0 + (105.0*std::pow( (a + b - 2*x),3) )/(64.0*std::pow( (a - b),2) ) - (693.0*std::pow( (a + b - 2*x),5) )/(128.0*std::pow((a - b),4) ) + (429.0*std::pow((a + b - 2*x),7))/(64.0*std::pow((a - b),6)) - (715.0*std::pow( (a + b - 2*x),9) )/(256.0*std::pow((a - b),8));
//     }

//     throw  std::invalid_argument("MomentFitting :: Order out of range!\n");
// }

// double MomentFitting::p_n(double x, int order) {
//     switch(order)
//     {
//         case 0:
//             return 1;
//         case 1:
//             return x;
//         case 2:
//             return 1.0/2.0*(3.0*std::pow(x,2)-1.0);
//         case 3:
//             return 1.0/2.0*(5.0*std::pow(x,3) - 3.0*x);
//         case 4:
//             return 1.0/8.0*(35.0*std::pow(x,4)-30.0*std::pow(x,2) +3.0);
//         case 5:
//             return 1.0/8.0*(63.0*std::pow(x,5)-70.0*std::pow(x,3)+15.0*x);
//         case 6:
//             return 1.0/16.0*(231.0*std::pow(x,6)-315.0*std::pow(x,4)+105.0*std::pow(x,2)-5.0);
//         case 7:
//             return 1.0/16.0*(429.0*std::pow(x,7)-693.0*std::pow(x,5)+315.0*std::pow(x,3)-35.0*x);
//         case 8:
//             return 1.0/128.0*(6435.0*std::pow(x,8) - 12012.0*std::pow(x,6)+6930.0*std::pow(x,4)-1260.0*std::pow(x,2)+35.0);
//     }

//     throw  std::invalid_argument("MomentFitting :: Order out of range!\n");
// }


