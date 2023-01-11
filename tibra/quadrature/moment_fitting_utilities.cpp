// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <stdexcept>
#include <cmath>
//// Project includes
#include "define.hpp"
#include "quadrature/moment_fitting_utilities.h"
#include "utilities/mapping_utilities.h"
#include "utilities/polynomial_utilities.h"
#include "containers/element.hpp"
#include "solvers/nnls.h"
#include "io/io_utilities.h"

namespace tibra {

typedef boost::numeric::ublas::matrix<double> MatrixType;
typedef boost::numeric::ublas::vector<double> VectorType;

void MomentFitting::DistributeInitialIntegrationPoints(const Element& rElement, IntegrationPointVectorType& rIntegrationPoint, SizeType PointDistributionFactor, const Parameters& rParam){

    const double factor = PointDistributionFactor;
    rIntegrationPoint.reserve( (int) factor*(rParam.Order()[0]+1)*factor*(rParam.Order()[1]+1)*factor*(rParam.Order()[2]+1) );

    const auto bounding_box = rElement.pGetTrimmedDomain()->GetBoundingBoxOfTrimmedDomain();

    const double delta_x = std::abs(bounding_box.second[0] - bounding_box.first[0]) / ( (double)factor*(rParam.Order()[0]+1) );
    const double delta_y = std::abs(bounding_box.second[1] - bounding_box.first[1]) / ( (double)factor*(rParam.Order()[1]+1) );
    const double delta_z = std::abs(bounding_box.second[2] - bounding_box.first[2]) / ( (double)factor*(rParam.Order()[2]+1) );

    const auto lower_bound = rParam.LowerBound();
    const auto upper_bound = rParam.UpperBound();

    double xx = bounding_box.first[0] + 0.5*delta_x;
    while( xx < bounding_box.second[0] ){
        double yy = bounding_box.first[1] + 0.5*delta_y;
        while( yy < bounding_box.second[1] ){
            double zz = bounding_box.first[2] + 0.5*delta_z;
            while( zz < bounding_box.second[2]){
                const PointType tmp_point = {xx,yy,zz};
                if( rElement.pGetTrimmedDomain()->IsInsideTrimmedDomain(tmp_point) ){
                    auto tmp_point_param = Mapping::GlobalToParam(tmp_point, lower_bound, upper_bound );
                    rIntegrationPoint.push_back(IntegrationPoint(tmp_point_param, 0.0 ));
                }
                zz += delta_z;
            }
            yy += delta_y;
        }
        xx += delta_x;
    }
}

void MomentFitting::ComputeConstantTerms(const Element& rElement, const BoundaryIPsVectorPtrType& pBoundaryIps,
                                         VectorType& rConstantTerms, const Parameters& rParam){

    const auto lower_bound = rParam.LowerBound();
    const auto upper_bound = rParam.UpperBound();

    const double jacobian_x = std::abs(lower_bound[0] - upper_bound[0]);
    const double jacobian_y = std::abs(lower_bound[1] - upper_bound[1]);
    const double jacobian_z = std::abs(lower_bound[2] - upper_bound[2]);

    const PointType& a = rElement.GetLowerBoundParam();
    const PointType& b = rElement.GetUpperBoundParam();

    const int ffactor = 1;
    const int order_u = rParam.Order()[0];
    const int order_v = rParam.Order()[1];
    const int order_w = rParam.Order()[2];

    const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor+1) * (order_w*ffactor + 1);

    rConstantTerms.resize(number_of_functions, false);
    std::fill( rConstantTerms.begin(),rConstantTerms.end(), 0.0);

    IndexType row_index = 0;
    // Loop over all points.
    const auto begin_points_it_ptr = pBoundaryIps->begin();

    // X-direction
    std::vector<double> f_x_x(order_u*ffactor+1);
    std::vector<double> f_x_int_x(order_u*ffactor+1);
    // Y-direction
    std::vector<double> f_x_y(order_v*ffactor+1);
    std::vector<double> f_x_int_y(order_v*ffactor+1);
    // Z-direction
    std::vector<double> f_x_z(order_w*ffactor+1);
    std::vector<double> f_x_int_z(order_w*ffactor+1);
    for( int i = 0; i < pBoundaryIps->size(); ++i ){
        // Note: The evaluation of polynomials is expensive. Therefore, precompute and store values
        // for f_x_x and f_x_int at each point.
        auto point_it = (begin_points_it_ptr + i);
        const auto& normal = point_it->Normal();
        PointType local_point = Mapping::GlobalToParam(*point_it, lower_bound, upper_bound);

        // X-Direction
        for( IndexType i_x = 0; i_x <= order_u*ffactor; ++i_x){
            f_x_x[i_x] = Polynomial::f_x(local_point[0], i_x, a[0], b[0]);
            f_x_int_x[i_x] = Polynomial::f_x_int(local_point[0], i_x, a[0], b[0]);
        }
        // Y-Direction
        for( IndexType i_y = 0; i_y <= order_v*ffactor; ++i_y){
            f_x_y[i_y] = Polynomial::f_x(local_point[1], i_y, a[1], b[1]);
            f_x_int_y[i_y] = Polynomial::f_x_int(local_point[1], i_y, a[1], b[1]);
        }
        // Z-Direction
        for( IndexType i_z = 0; i_z <= order_w*ffactor; ++i_z){
            f_x_z[i_z] = Polynomial::f_x(local_point[2], i_z, a[2], b[2]);
            f_x_int_z[i_z] = Polynomial::f_x_int(local_point[2], i_z, a[2], b[2]);
        }

        const double weight = 1.0/3.0*point_it->GetWeight();
        row_index = 0;
        // Assembly RHS
        for( int i_x = 0; i_x <= order_u*ffactor; ++i_x){
            for( int i_y = 0; i_y <= order_v*ffactor; ++i_y ){
                for( int i_z = 0; i_z <= order_w*ffactor; ++i_z){
                    // Compute normal for each face/triangle.
                    PointType value;
                    value[0] = f_x_int_x[i_x]*f_x_y[i_y]*f_x_z[i_z];
                    value[1] = f_x_x[i_x]*f_x_int_y[i_y]*f_x_z[i_z];
                    value[2] = f_x_x[i_x]*f_x_y[i_y]*f_x_int_z[i_z];

                    double integrand = normal[0]*value[0]*jacobian_x + normal[1]*value[1]*jacobian_y + normal[2]*value[2]*jacobian_z;
                    rConstantTerms[row_index] += integrand * weight;
                    row_index++;
                }

            }
        }
    }
}


void MomentFitting::CreateIntegrationPointsTrimmed(Element& rElement, const Parameters& rParam){

    double residual = 1e10;
    SizeType point_distribution_factor = rParam.GetPointDistributionFactor();
    VectorType constant_terms{};

    const auto p_trimmed_domain = rElement.pGetTrimmedDomain();
    const auto p_boundary_ips = p_trimmed_domain->pGetBoundaryIps();

    ComputeConstantTerms(rElement, p_boundary_ips, constant_terms, rParam);

    const int max_iteration = 3;
    int iteration = 1;
    while( residual > rParam.MomentFittingResidual() && iteration < max_iteration){
        auto& reduced_points = rElement.GetIntegrationPoints();
        if( !rParam.UseCustomizedTrimmedPositions() ){
            // This is only used for test_moment_fitting.cpp
            reduced_points.clear();
        }
        residual = CreateIntegrationPointsTrimmed(rElement, constant_terms, point_distribution_factor, rParam);

        point_distribution_factor *= 2;
        if( residual > 1e-2 ) {
            reduced_points.clear();
        }
        iteration++;
    }

    if( residual > rParam.MomentFittingResidual() && rParam.EchoLevel() > 2){
        //std::cout << "size: " << rElement.GetIntegrationPoints().size() << std::endl;
        std::cout << "Moment Fitting :: Targeted residual can not be achieved!: " << residual << std::endl;
    }
}

double MomentFitting::CreateIntegrationPointsTrimmed(Element& rElement, const VectorType& rConstantTerms, SizeType PointDistributionFactor, const Parameters& rParam) {

    IntegrationPointVectorType new_integration_points{};
    new_integration_points.resize(0UL);
    int maximum_iteration;
    if( !rParam.UseCustomizedTrimmedPositions() ){
        SizeType point_distribution_factor = PointDistributionFactor;
        const SizeType min_num_points = (rParam.Order()[0]+1)*(rParam.Order()[1]+1)*(rParam.Order()[2]+1)*(point_distribution_factor);
        while( new_integration_points.size() < min_num_points ){
            DistributeInitialIntegrationPoints( rElement, new_integration_points, point_distribution_factor, rParam);
            point_distribution_factor *= 2;
        }
        maximum_iteration = 1000;
    }
    else { // This is only used for test_moment_fitting.cpp
        new_integration_points = rElement.GetIntegrationPoints();
        rElement.GetIntegrationPoints().clear();
        maximum_iteration = 1;
    }

    const double jacobian_x = std::abs(rParam.UpperBound()[0] - rParam.LowerBound()[0]);
    const double jacobian_y = std::abs(rParam.UpperBound()[1] - rParam.LowerBound()[1]);
    const double jacobian_z = std::abs(rParam.UpperBound()[2] - rParam.LowerBound()[2]);

    PointType a = rElement.GetLowerBoundParam();
    PointType b = rElement.GetUpperBoundParam();

    const int ffactor = 1;
    const int order_u = rParam.Order()[0];
    const int order_v = rParam.Order()[1];
    const int order_w = rParam.Order()[2];

    const IndexType number_of_functions = (order_u*ffactor + 1) * (order_v*ffactor+1) * (order_w*ffactor + 1);

    double global_residual = 1e-15;
    const double allowed_residual = rParam.MomentFittingResidual();
    int number_iterations = 0;
    bool change = false;
    IntegrationPointVectorType prev_solution{};
    double prev_residual = 0.0;

    /// Enter point elimination algorithm
    while( change || (global_residual < allowed_residual && number_iterations < maximum_iteration) ){

        const IndexType number_reduced_points = new_integration_points.size();

        /// Assemble moment fitting matrix.
        MatrixType fitting_matrix(number_of_functions, number_reduced_points);
        IndexType row_index = 0;
        for( int i_x = 0; i_x <= order_u*ffactor; ++i_x){
            for( int i_y = 0; i_y <= order_v*ffactor; ++i_y ){
                for( int i_z = 0; i_z <= order_w*ffactor; ++i_z){
                    // Loop over all points
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
        global_residual = nnls::nnls(fitting_matrix, rConstantTerms, weights)/number_of_functions;

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
                std::sort(new_integration_points.begin(), new_integration_points.end(), [](const IntegrationPoint& point_a, const IntegrationPoint& point_b) -> bool {
                        return point_a.GetWeight() > point_b.GetWeight();
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
                    if( it->GetWeight() < EPS1*max_value && new_integration_points.size() > 4){
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
    auto& reduced_points = rElement.GetIntegrationPoints();

    if( (global_residual >= allowed_residual && prev_solution.size() > 0 && number_iterations < maximum_iteration) ) {
        reduced_points.insert(reduced_points.begin(), prev_solution.begin(), prev_solution.end());
        reduced_points.erase(std::remove_if(reduced_points.begin(), reduced_points.end(), [](const IntegrationPoint& point) {
            return point.GetWeight() < EPS4; }), reduced_points.end());
        return prev_residual;
    }
    else{
        reduced_points.insert(reduced_points.begin(), new_integration_points.begin(), new_integration_points.end());
        reduced_points.erase(std::remove_if(reduced_points.begin(), reduced_points.end(), [](const IntegrationPoint& point) {
            return point.GetWeight() < EPS4; }), reduced_points.end());
        return global_residual;
    }

}

} // End namespace tibra

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


