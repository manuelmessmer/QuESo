// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// External includes
#include <stdexcept>

// Project includes
#include "utilities/integration_points/integration_points_factory.h"
#include "utilities/integration_point_utilities.h"

void IntegrationPointUtilities::IntegrationPoints3D(
        IntersectionTest& rInsideTest,
        IntegrationPointType& rIntegrationPoints,
        SizeType PointsInU, SizeType PointsInV, SizeType PointsInW,
        double U0, double U1, double V0, double V1, double W0, double W1)
{

   if(PointsInU < 1 || PointsInV < 1 || PointsInW < 1){
        throw  std::invalid_argument("PointsInU, -V and -W need to be bigger than 0 - PointsInU:" + std::to_string(PointsInU)
            + ", PointsInV:" + std::to_string(PointsInV) + " and PointsInW:" + std::to_string(PointsInW) + "\n");
   }

    const double distance_u = U1 - U0;
    const double length_u = std::abs(U1 - U0);
    const double distance_v = V1 - V0;
    const double length_v = std::abs(V1 - V0);
    const double distance_w = W1 - W0;
    const double length_w = std::abs(W1 - W0);

    const std::vector<std::array<double, 2>>& integration_point_list_u = IntegrationPointFactory::GetIntegrationPoints(PointsInU, IntegrationPointFactory::IntegrationMethod::Gauss);
    const std::vector<std::array<double, 2>>& integration_point_list_v = IntegrationPointFactory::GetIntegrationPoints(PointsInV, IntegrationPointFactory::IntegrationMethod::Gauss);
    const std::vector<std::array<double, 2>>& integration_point_list_w = IntegrationPointFactory::GetIntegrationPoints(PointsInW, IntegrationPointFactory::IntegrationMethod::Gauss);

    for (SizeType u = 0; u < PointsInU; ++u) {
        for (SizeType v = 0; v < PointsInV; ++v) {
            for( SizeType w = 0; w < PointsInW; ++w) {
                std::array<double,3> tmp_point = {  U0 + distance_u * integration_point_list_u[u][0],
                                                    V0 + distance_v * integration_point_list_v[v][0],
                                                    W0 + distance_w * integration_point_list_w[w][0] };

                if( !rInsideTest.IsInsideLocalCoordinates(tmp_point) ){
                    rIntegrationPoints.push_back( IntegrationPoint( tmp_point[0],
                                                                    tmp_point[1],
                                                                    tmp_point[2],
                                                                    integration_point_list_u[u][1] * length_u *
                                                                    integration_point_list_v[v][1] * length_v *
                                                                    integration_point_list_w[w][1] * length_w ) );
                }
            }
        }
    }
}

void IntegrationPointUtilities::IntegrationPoints3D(
        IntegrationPointType& rIntegrationPoints,
        SizeType PointsInU, SizeType PointsInV, SizeType PointsInW,
        double U0, double U1, double V0, double V1, double W0, double W1)
{

   if(PointsInU < 1 || PointsInV < 1 || PointsInW < 1){
        throw  std::invalid_argument("PointsInU, -V and -W need to be bigger than 0 - PointsInU:" + std::to_string(PointsInU)
            + ", PointsInV:" + std::to_string(PointsInV) + " and PointsInW:" + std::to_string(PointsInW) + "\n");
   }

    const double distance_u = U1 - U0;
    const double length_u = std::abs(U1 - U0);
    const double distance_v = V1 - V0;
    const double length_v = std::abs(V1 - V0);
    const double distance_w = W1 - W0;
    const double length_w = std::abs(W1 - W0);

    const std::vector<std::array<double, 2>>& integration_point_list_u = IntegrationPointFactory::GetIntegrationPoints(PointsInU, IntegrationPointFactory::IntegrationMethod::Gauss);
    const std::vector<std::array<double, 2>>& integration_point_list_v = IntegrationPointFactory::GetIntegrationPoints(PointsInV, IntegrationPointFactory::IntegrationMethod::Gauss);
    const std::vector<std::array<double, 2>>& integration_point_list_w = IntegrationPointFactory::GetIntegrationPoints(PointsInW, IntegrationPointFactory::IntegrationMethod::Gauss);

    for (SizeType u = 0; u < PointsInU; ++u) {
        for (SizeType v = 0; v < PointsInV; ++v) {
            for( SizeType w = 0; w < PointsInW; ++w) {
                rIntegrationPoints.push_back( IntegrationPoint( U0 + distance_u * integration_point_list_u[u][0],
                                                                V0 + distance_v * integration_point_list_v[v][0],
                                                                W0 + distance_w * integration_point_list_w[w][0],
                                                                integration_point_list_u[u][1] * length_u *
                                                                integration_point_list_v[v][1] * length_v *
                                                                integration_point_list_w[w][1] * length_w ) );
            }
        }
    }
}

void IntegrationPointUtilities::CreateGaussLegendrePoints(
    IntersectionTest& rInsideTest,
    IntegrationPointType& rIntegrationPoints,
    std::array<double,3> point_A,
    std::array<double,3> point_B,
    SizeType OrderU,
    SizeType OrderV,
    SizeType OrderW)
{

    const SizeType number_of_integration_points = (OrderU + 1)*(OrderV + 1)*(OrderW + 1);

    IntegrationPointType tmp_integration_points;
    tmp_integration_points.reserve(number_of_integration_points);
    IntegrationPoints3D(
        rInsideTest,
        tmp_integration_points,
        OrderU+1, OrderV+1, OrderW+1,
        point_A[0], point_B[0],
        point_A[1], point_B[1],
        point_A[2], point_B[2]);

    //#pragma omp critical
    rIntegrationPoints.insert(rIntegrationPoints.end(), tmp_integration_points.begin(), tmp_integration_points.end());
}

void IntegrationPointUtilities::CreateGaussLegendrePoints(
    IntegrationPointType& rIntegrationPoints,
    std::array<double,3> point_A,
    std::array<double,3> point_B,
    SizeType OrderU,
    SizeType OrderV,
    SizeType OrderW)
{

    const SizeType number_of_integration_points = (OrderU +1)*(OrderV +1)*(OrderW +1);

    IntegrationPointType tmp_integration_points;
    tmp_integration_points.reserve(number_of_integration_points);
    IntegrationPointUtilities::IntegrationPoints3D(
        tmp_integration_points,
        OrderU+1, OrderV+1, OrderW+1,
        point_A[0], point_B[0],
        point_A[1], point_B[1],
        point_A[2], point_B[2]);

    //#pragma omp critical
    rIntegrationPoints.insert(rIntegrationPoints.end(), tmp_integration_points.begin(), tmp_integration_points.end());
}
