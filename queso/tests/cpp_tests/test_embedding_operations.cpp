// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "includes/checks.hpp"
#include "utilities/mapping_utilities.h"
#include "io/io_utilities.h"
#include "QuESo_main.h"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( EmbeddingOperationsTestSuite )

BOOST_AUTO_TEST_CASE(Intersection) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Intersected Knot Span" << std::endl;

    std::string filename = "queso/tests/cpp_tests/data/cylinder.stl";
    Parameters parameters( {Component("input_filename", filename),
                            Component("echo_level", 0UL),
                            Component("lower_bound_xyz", PointType(0.0, 0.0, 0.0)),
                            Component("upper_bound_xyz", PointType(2.0, 2.0, 1.0)),
                            Component("lower_bound_uvw", PointType(0.0, 0.0, 0.0)),
                            Component("upper_bound_uvw", PointType(4.0, 5.0, 3.0)),
                            Component("number_of_elements", Vector3i(1, 1, 1)),
                            Component("polynomial_order", Vector3i(2, 2, 2)),
                            Component("integration_method", IntegrationMethod::Gauss),
                            Component("init_point_distribution_factor", 3UL),
                            Component("min_num_boundary_triangles", 5000UL),
                            Component("moment_fitting_residual", 1e-8) });
    QuESo queso(parameters);
    queso.Run();

    const auto& elements = queso.GetElements();

    QuESo_CHECK_EQUAL(elements.size(), 1);

    const auto& points_reduced = (*elements.begin())->GetIntegrationPoints();
    QuESo_CHECK_LT(points_reduced.size(), 28);

    const auto& r_triangle_mesh = (*elements.begin())->pGetTrimmedDomain()->GetTriangleMesh();
    const IndexType num_triangles = r_triangle_mesh.NumOfTriangles();

    QuESo_CHECK_GT(num_triangles, parameters.MinimumNumberOfTriangles());

    double area = 0.0;
    for( IndexType triangle_id = 0; triangle_id < num_triangles; ++triangle_id){
        const PointType coordinates = r_triangle_mesh.Center(triangle_id);

        QuESo_CHECK_GT(coordinates[2], -1e-6);
        QuESo_CHECK_LT(coordinates[2], 1.0+1e-6);

        QuESo_CHECK_GT(coordinates[0], -1e-6);
        QuESo_CHECK_GT(coordinates[1], -1e-6);

        if( coordinates[0] > 1e-6 && coordinates[1] > 1e-6
            && coordinates[2] > 1e-6 && coordinates[2] < 1-1e-6){ // check in z richtung fehlt noch
            double radius = std::sqrt( std::pow(coordinates[0],2) + std::pow(coordinates[1], 2) );
            QuESo_CHECK_GT(radius, 0.998);
        }
        area += r_triangle_mesh.Area(triangle_id);
    }
    QuESo_CHECK_LT(area, 5.141592654);
    QuESo_CHECK_GT(area, 5.135);
}

void TestElephant( IntegrationMethodType IntegrationMethod, const Vector3i&  rOrder, IndexType NumPointsInside, double Tolerance,
                        bool BSplineMesh, const BoundingBoxType& rBoundsUVW, bool Large){

    Vector3i num_elements = (Large) ? Vector3i(14, 22, 12) : Vector3i(7, 11, 6);

    std::string filename = "queso/tests/cpp_tests/data/elephant.stl";
    Parameters parameters( {Component("input_filename", filename),
                            Component("echo_level", 0UL),
                            Component("lower_bound_xyz", PointType(-0.37, -0.55, -0.31)),
                            Component("upper_bound_xyz", PointType(0.37, 0.55, 0.31)),
                            Component("b_spline_mesh", BSplineMesh),
                            Component("lower_bound_uvw", rBoundsUVW.first),
                            Component("upper_bound_uvw", rBoundsUVW.second),
                            Component("number_of_elements", num_elements),
                            Component("polynomial_order", rOrder),
                            Component("integration_method", IntegrationMethod) });
    QuESo queso(parameters);
    queso.Run();

    const auto& elements = queso.GetElements();

    // Compute total volume
    double volume_trimmed = 0.0;
    double volume_inside = 0.0;
    const auto el_it_begin = elements.begin();
    int num_elements_inside = 0;
    int num_elements_trimmed = 0;
    int num_points_inside = 0;
    for( IndexType i = 0; i < elements.size(); ++i){
        auto el_it = *(el_it_begin+i);
        const double det_j = el_it->DetJ();
        if( el_it->IsTrimmed() ){
            const auto& points_trimmed = el_it->GetIntegrationPoints();
            QuESo_CHECK_GT(points_trimmed.size(), 0);
            QuESo_CHECK_LT(points_trimmed.size(), (rOrder[0]+1)*(rOrder[1]+1)*(rOrder[2]+1)+1);
            for( const auto& point : points_trimmed ){
                volume_trimmed += point.GetWeight()*det_j;
            }
            num_elements_trimmed++;
        } else {
            const auto& points_inside = el_it->GetIntegrationPoints();
            for( const auto& point : points_inside ){
                volume_inside += point.GetWeight()*det_j;
                num_points_inside++;
            }
            num_elements_inside++;
        }
    }
    if( Large ) {
        QuESo_CHECK_EQUAL(num_elements_inside, 108);
    } else {
        QuESo_CHECK_EQUAL(num_elements_inside, 5);
    }
    QuESo_CHECK_EQUAL(num_points_inside, static_cast<int>(NumPointsInside));

    // Check volume inside
    const BoundingBoxType el_bounding_box = (*el_it_begin)->GetBoundsXYZ();
    const PointType el_delta = el_bounding_box.second - el_bounding_box.first;
    const double ref_volume_inside = (el_delta[0]*el_delta[1]*el_delta[2])*num_elements_inside;
    QuESo_CHECK_RELATIVE_NEAR(volume_inside, ref_volume_inside, 1e-13);

    // Check volume inside + trimmed
    const double volume_tot = (volume_trimmed+volume_inside);
    const double ref_volume_tot = 0.0462012;
    QuESo_CHECK_RELATIVE_NEAR(volume_tot, ref_volume_tot, Tolerance);
}

// p=2
BOOST_AUTO_TEST_CASE(VolumeElephant1) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant Gauss (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::Gauss, Vector3i(2, 2, 2), 2916, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::Gauss, Vector3i(2, 2, 2), 2916, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
    TestElephant(IntegrationMethod::Gauss, Vector3i(2, 2, 2), 2916, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant2) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant Optimal (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::GGQ_Optimal, Vector3i(2, 2, 2), 1786, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::GGQ_Optimal, Vector3i(2, 2, 2), 1786, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant3) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant GGQ_Reduced1 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::GGQ_Reduced1, Vector3i(2, 2, 2), 673, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::GGQ_Reduced1, Vector3i(2, 2, 2), 673, 0.0002, use_b_spline_mesh, MakeBox({-6.0, -7.5, -1.2}, {22.0, 1.82, 2.8}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant4) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant GGQ_Reduced2 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::GGQ_Reduced2, Vector3i(2, 2, 2), 406, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::GGQ_Reduced2, Vector3i(2, 2, 2), 406, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), true);
}

// p=3
BOOST_AUTO_TEST_CASE(VolumeElephant5) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant Gauss (p=3)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::Gauss, Vector3i(3, 3, 3), 320, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), false);
    TestElephant(IntegrationMethod::Gauss, Vector3i(3, 3, 3), 320, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), false);
    TestElephant(IntegrationMethod::Gauss, Vector3i(3, 3, 3), 320, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), false);
}

// p=3
BOOST_AUTO_TEST_CASE(VolumeElephant6) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant GGQ_Optimal (p=3)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::GGQ_Optimal, Vector3i(3, 3, 3), 256, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), false);
    TestElephant(IntegrationMethod::GGQ_Optimal, Vector3i(3, 3, 3), 256, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), false);
}

// p=4
BOOST_AUTO_TEST_CASE(VolumeElephant7) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant Gauss (p=4)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::Gauss, Vector3i(4, 4, 4), 625, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), false);
    TestElephant(IntegrationMethod::Gauss, Vector3i(4, 4, 4), 625, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), false);
    TestElephant(IntegrationMethod::Gauss, Vector3i(4, 4, 4), 625, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), false);
}

// p=4
BOOST_AUTO_TEST_CASE(VolumeElephant8) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant GGQ_Optimal (p=4)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::GGQ_Optimal, Vector3i(4, 4, 4), 525, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), false);
}

// p=Mix
BOOST_AUTO_TEST_CASE(VolumeElephant9) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant Gauss (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::Gauss, Vector3i(2, 3, 4), 6480, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::Gauss, Vector3i(3, 4, 2), 6480, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
    TestElephant(IntegrationMethod::Gauss, Vector3i(4, 2, 3), 6480, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant10) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant Gauss_Reduced1 (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::Gauss_Reduced1, Vector3i(2, 3, 4), 2592, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::Gauss_Reduced1, Vector3i(3, 4, 2), 2592, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
    TestElephant(IntegrationMethod::Gauss_Reduced1, Vector3i(4, 2, 3), 2592, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant11) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant Gauss_Reduced2 (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::Gauss_Reduced2, Vector3i(2, 3, 4), 648, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::Gauss_Reduced2, Vector3i(3, 4, 2), 648, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
    TestElephant(IntegrationMethod::Gauss_Reduced2, Vector3i(4, 2, 3), 648, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant12) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant Optimal (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::GGQ_Optimal, Vector3i(2, 3, 4), 3657, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::GGQ_Optimal, Vector3i(3, 4, 2), 3682, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant13) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant GGQ_Reduced1 (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::GGQ_Reduced1, Vector3i(2, 3, 4), 1710, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::GGQ_Reduced1, Vector3i(3, 4, 2), 1713, 0.0002, use_b_spline_mesh, MakeBox({-6.0, -7.5, -1.2}, {22.0, 1.82, 2.8}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant14) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant GGQ_Reduced2 (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::GGQ_Reduced2, Vector3i(2, 3, 4), 1166, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::GGQ_Reduced2, Vector3i(3, 4, 2), 1192, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), true);
}

void TestSteeringKnuckle( IntegrationMethodType IntegrationMethod, IndexType p, IndexType NumPointsInside, double Tolerance,
                        bool BSplineMesh, const BoundingBoxType& rBoundsUVW, bool Large){

    Vector3i num_elements = (Large) ? Vector3i(20, 40, 40) : Vector3i(5, 10, 10);

    std::string filename = "queso/tests/cpp_tests/data/steering_knuckle.stl";
    Parameters parameters( {Component("input_filename", filename),
                            Component("echo_level", 0UL),
                            Component("lower_bound_xyz", PointType(-130.0, -110.0, -110.0)),
                            Component("upper_bound_xyz", PointType(20.0, 190.0, 190.0)),
                            Component("b_spline_mesh", BSplineMesh),
                            Component("lower_bound_uvw", rBoundsUVW.first),
                            Component("upper_bound_uvw", rBoundsUVW.second),
                            Component("number_of_elements", num_elements),
                            Component("polynomial_order", Vector3i(p, p, p)),
                            Component("integration_method", IntegrationMethod) });
    QuESo queso(parameters);
    queso.Run();

    const auto& elements = queso.GetElements();

    const PointType delta_uvw = rBoundsUVW.second - rBoundsUVW.first;
    // Compute total volume
    double volume_trimmed = 0.0;
    double volume_inside = 0.0;
    const auto el_it_begin = elements.begin();
    int num_elements_inside = 0;
    int num_elements_trimmed = 0;
    int num_points_inside = 0;
    for( IndexType i = 0; i < elements.size(); ++i){
        auto el_it = *(el_it_begin+i);
        const double det_j = el_it->DetJ();
        if( el_it->IsTrimmed() ){
            const auto& points_trimmed = el_it->GetIntegrationPoints();
            QuESo_CHECK_GT(points_trimmed.size(), 0);
            QuESo_CHECK_LT(points_trimmed.size(), (p+1)*(p+1)*(p+1)+1);
            for( const auto& point : points_trimmed ){
                volume_trimmed += point.GetWeight()*det_j;
            }
            num_elements_trimmed++;
        } else {
            const auto& points_inside = el_it->GetIntegrationPoints();
            for( const auto& point : points_inside ){
                volume_inside += point.GetWeight()*det_j;
                num_points_inside++;
            }
            num_elements_inside++;
        }
    }
    if( Large ) {
        QuESo_CHECK_EQUAL(num_elements_inside, 201);
    } else {
        QuESo_CHECK_EQUAL(num_elements_inside, 5);
    }
    QuESo_CHECK_EQUAL(num_points_inside, static_cast<int>(NumPointsInside));

    // Check volume inside
    const BoundingBoxType el_bounding_box = (*el_it_begin)->GetBoundsXYZ();
    const PointType delta = el_bounding_box.second - el_bounding_box.first;
    const double ref_volume_inside = (delta[0]*delta[1]*delta[2])*num_elements_inside;
    QuESo_CHECK_RELATIVE_NEAR(volume_inside, ref_volume_inside, 1e-13);

    // Check volume inside + trimmed
    const double volume_tot = (volume_trimmed+volume_inside);
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, filename.c_str());
    const double ref_volume_tot = MeshUtilities::Volume(triangle_mesh);
    QuESo_CHECK_RELATIVE_NEAR(volume_tot, ref_volume_tot, Tolerance);
}

// p=2
BOOST_AUTO_TEST_CASE(SteeringKnuckle1) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Steering Knuckle Gauss (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::Gauss, 2, 5427, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss, 2, 5427, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss, 2, 5427, 0.0002, use_b_spline_mesh, MakeBox({-12.0, -1110.0, -2110.0}, {220.0, 3190.0, 3190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss, 2, 5427, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss, 2, 5427, 0.0002, false, MakeBox({-2, -5.0, -10.0}, {3.0, 19.0, 10.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle2) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Steering Knuckle Gauss_Reduced1 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced1, 2, 1608, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced1, 2, 1608, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced1, 2, 1608, 0.0002, use_b_spline_mesh, MakeBox({-12.0, -1110.0, -2110.0}, {220.0, 3190.0, 3190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced1, 2, 1608, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced1, 2, 1608, 0.0002, false, MakeBox({-2, -5.0, -10.0}, {3.0, 19.0, 10.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle3) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Steering Knuckle Gauss_Reduced2 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced2, 2, 201, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced2, 2, 201, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced2, 2, 201, 0.0002, use_b_spline_mesh, MakeBox({-12.0, -1110.0, -2110.0}, {220.0, 3190.0, 3190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced2, 2, 201, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::Gauss_Reduced2, 2, 201, 0.0002, false, MakeBox({-2, -5.0, -10.0}, {3.0, 19.0, 10.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle4) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Steering Knuckle Optimal (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::GGQ_Optimal, 2, 3483, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::GGQ_Optimal, 2, 3483, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::GGQ_Optimal, 2, 3483, 0.0002, use_b_spline_mesh, MakeBox({-2.0, -5.0, -10.0}, {3.0, 19.0, 10.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle5) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Steering Knuckle GGQ_Reduced1 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::GGQ_Reduced1, 2, 1320, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::GGQ_Reduced1, 2, 1320, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle6) {
    QuESo_INFO << "Testing :: Test Embedding Operations :: Volume Elephant GGQ_Reduced2 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::GGQ_Reduced2, 2, 804, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::GGQ_Reduced2, 2, 804, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso