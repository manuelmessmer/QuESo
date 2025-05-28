//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#define BOOST_TEST_DYN_LINK

//// External includes
#include <boost/test/unit_test.hpp>
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/containers/background_grid.hpp"
#include "queso/quadrature/trimmed_element.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/brep_operator.h"
#include "queso/io/io_utilities.h"
#include "queso/tests/cpp_tests/class_testers/trimmed_element_tester.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( PointEliminationTestSuite )

void RunCylinder(const Vector3i& rOrder, double Residual){
    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-1.5, -1.5, -1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{1.5, 1.5, 12.0});
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{0.0, 0.0, 0.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{6, 6, 13});

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/cylinder.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 1e-3;
    const IndexType min_num_triangles = 500;

    GridIndexer grid_indexer(r_settings);
    IndexType number_trimmed_elements = 0;
    for( IndexType i = 0; i < grid_indexer.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
        Vector3d lower_bound_xyz = bounding_box.first;
        Vector3d upper_bound_xyz = bounding_box.second;

        const BoundingBoxType bounding_box_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(i);
        Vector3d lower_bound_uvw = bounding_box_uvw.first;
        Vector3d upper_bound_uvw = bounding_box_uvw.second;

        // Construct element
        ElementType element(1, MakeBox(lower_bound_xyz, upper_bound_xyz),
                               MakeBox(lower_bound_uvw, upper_bound_uvw));

        if( brep_operator.GetIntersectionState(lower_bound_xyz, upper_bound_xyz) == IntersectionState::trimmed){
            // Get trimmed domain
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound_xyz, upper_bound_xyz, min_vol_ratio, min_num_triangles);
            if( p_trimmed_domain ){
                ++number_trimmed_elements;

                // Add trimmed domain to element.
                element.SetIsTrimmed(true);
                element.pSetTrimmedDomain(p_trimmed_domain);

                // Run point elimination
                const auto residual = QuadratureTrimmedElementTester<ElementType>::AssembleIPs(element, rOrder, Residual);

                // Check if residual is smaller than targeted.
                QuESo_CHECK_LT(residual, 1e-6);

                // Must be more points than p*p*p.
                auto& r_points = element.GetIntegrationPoints();

                QuESo_CHECK_LT(r_points.size(), (rOrder[0]+1)*(rOrder[1]+1)*(rOrder[2]+1)+1);
                QuESo_CHECK_GT(r_points.size(), rOrder[0]*rOrder[1]*rOrder[2]);

                // Get copy of points.
                ElementType::IntegrationPointVectorType copy_points(r_points);

                // Compute constant terms.
                std::vector<double> constant_terms{};
                auto p_boundary_ips = element.pGetTrimmedDomain()->pGetBoundaryIps<BoundaryIntegrationPoint>();
                QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, p_boundary_ips, element, rOrder);

                // Run moment fitting again.
                const auto residual_2 = QuadratureTrimmedElementTester<ElementType>::MomentFitting(constant_terms, r_points, element, rOrder);

                // Check if residual and weights are the same.
                QuESo_CHECK_NEAR( residual, residual_2, EPS4 );
                double volume = 0.0;
                for( IndexType i = 0; i < r_points.size(); ++i){
                    const double weight1 = r_points[i].Weight();
                    QuESo_CHECK_GT(weight1, EPS4);
                    const double weight2 = copy_points[i].Weight();
                    QuESo_CHECK_RELATIVE_NEAR( weight1, weight2 , EPS2 );
                    volume += weight1*element.DetJ(); // Multiplied with det(J).
                }
                // Check if integration points contain correct volume;
                const auto& r_mesh = element.pGetTrimmedDomain()->GetTriangleMesh();
                const double ref_volume = MeshUtilities::Volume(r_mesh);
                QuESo_CHECK_RELATIVE_NEAR(volume, ref_volume, 100.0*Residual);
            }
        }
    }
    QuESo_CHECK_EQUAL(number_trimmed_elements, 120);
}

BOOST_AUTO_TEST_CASE(PointEliminationCylinder1Test) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Cylinder Quadratic" << std::endl;
    RunCylinder({2, 2, 2}, 1e-8);
}

BOOST_AUTO_TEST_CASE(PointEliminationCylinder2Test) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Cylinder Cubic" << std::endl;
    RunCylinder({3, 3, 3}, 1e-8);
}

BOOST_AUTO_TEST_CASE(PointEliminationCylinder4Test) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Cylinder Mixed" << std::endl;
    RunCylinder({2, 3, 4}, 1e-7);
}


BOOST_AUTO_TEST_CASE(PointEliminationKnuckleTest) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Knuckle" << std::endl;

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-130.0, -110.0, -110.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{-40, 10.0, 10.0});
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-130.0, -110.0, -110.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{-40, 10.0, 10.0});
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{9, 12, 12});

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/steering_knuckle.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 1e-3;
    const IndexType min_num_triangles = 500;

    GridIndexer grid_indexer(r_settings);
    IndexType number_trimmed_elements = 0;
    for( IndexType i = 0; i < grid_indexer.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
        Vector3d lower_bound_xyz = bounding_box.first;
        Vector3d upper_bound_xyz = bounding_box.second;

        const BoundingBoxType bounding_box_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(i);
        Vector3d lower_bound_uvw = bounding_box_uvw.first;
        Vector3d upper_bound_uvw = bounding_box_uvw.second;

        // Construct element
        ElementType element(1, MakeBox(lower_bound_xyz, upper_bound_xyz),
                               MakeBox(lower_bound_uvw, upper_bound_uvw));

        if( brep_operator.GetIntersectionState(lower_bound_xyz, upper_bound_xyz) == IntersectionState::trimmed){
            // Get trimmed domain
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound_xyz, upper_bound_xyz, min_vol_ratio, min_num_triangles);
            if( p_trimmed_domain ){
                ++number_trimmed_elements;

                // Add trimmed domain to element.
                element.SetIsTrimmed(true);
                element.pSetTrimmedDomain(p_trimmed_domain);

                // Run point elimination
                const auto residual = QuadratureTrimmedElementTester<ElementType>::AssembleIPs(element, {2, 2, 2}, 1e-8);

                // Check if residual is smaller than targeted.
                QuESo_CHECK_LT(residual, 1e-8);

                // Must be more points than p*p*p.
                auto& r_points = element.GetIntegrationPoints();
                QuESo_CHECK_LT(r_points.size(), 28);
                QuESo_CHECK_GT(r_points.size(), 7);

                // Get copy of points.
                ElementType::IntegrationPointVectorType copy_points(r_points);

                // Compute constant terms.
                std::vector<double> constant_terms{};
                auto p_boundary_ips = element.pGetTrimmedDomain()->pGetBoundaryIps<BoundaryIntegrationPoint>();
                QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, p_boundary_ips, element, {2, 2, 2});

                // Run moment fitting again.
                const auto residual_2 = QuadratureTrimmedElementTester<ElementType>::MomentFitting(constant_terms, r_points, element, {2, 2, 2});

                // Check if residual and weights are the same.
                QuESo_CHECK_LT( residual, residual_2+EPS4 );
                QuESo_CHECK_LT( residual_2, residual+EPS4 );
                double volume = 0.0;
                for( IndexType i = 0; i < r_points.size(); ++i){
                    const double weight1 = r_points[i].Weight();
                    QuESo_CHECK_GT(weight1, EPS4);
                    const double weight2 = copy_points[i].Weight();
                    const double error = std::abs(weight1 - weight2)/ weight1;
                    QuESo_CHECK_LT( error , EPS2 );
                    volume += weight1*element.DetJ(); // Multiplied with det(J).
                }
                // Check if integration points contain correct volume;
                const auto& r_mesh = element.pGetTrimmedDomain()->GetTriangleMesh();
                const double ref_volume = MeshUtilities::Volume(r_mesh);
                const double volume_error = std::abs(volume - ref_volume)/ ref_volume;
                QuESo_CHECK_LT(volume_error, 1e-6); // Note can not be better as moment fitting residual.
            }
        }
    }
    QuESo_CHECK_EQUAL(number_trimmed_elements, 80);
}

BOOST_AUTO_TEST_CASE(PointEliminationElephantTest) {
    QuESo_INFO << "Testing :: Test Point Elimination :: Elephant" << std::endl;

    typedef IntegrationPoint IntegrationPointType;
    typedef BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
    r_grid_settings.SetValue(BackgroundGridSettings::grid_type, GridType::hexahedral_fe_grid);
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-0.4, -0.6, -0.35});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{0.4, 0.6, 0.35});
    r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{-1.0, -1.0, -1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0,  1.0,  1.0});
    r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{8, 12, 7});

    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, "queso/tests/cpp_tests/data/elephant.stl");

    // Build brep_operator
    BRepOperator brep_operator(triangle_mesh);

    const double min_vol_ratio = 1e-3;
    const IndexType min_num_triangles = 500;

    GridIndexer grid_indexer(r_settings);
    IndexType number_trimmed_elements = 0;
    for( IndexType i = 0; i < grid_indexer.NumberOfElements(); ++i){
        const BoundingBoxType bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
        Vector3d lower_bound_xyz = bounding_box.first;
        Vector3d upper_bound_xyz = bounding_box.second;

        const BoundingBoxType bounding_box_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(i);
        Vector3d lower_bound_uvw = bounding_box_uvw.first;
        Vector3d upper_bound_uvw = bounding_box_uvw.second;

        // Construct element
        ElementType element(1, MakeBox(lower_bound_xyz, upper_bound_xyz),
                               MakeBox(lower_bound_uvw, upper_bound_uvw));

        if( brep_operator.GetIntersectionState(lower_bound_xyz, upper_bound_xyz) == IntersectionState::trimmed){
            // Get trimmed domain
            auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(lower_bound_xyz, upper_bound_xyz, min_vol_ratio, min_num_triangles);
            if( p_trimmed_domain ){
                ++number_trimmed_elements;

                // Add trimmed domain to element.
                element.SetIsTrimmed(true);
                element.pSetTrimmedDomain(p_trimmed_domain);

                // Run point elimination
                const auto residual = QuadratureTrimmedElementTester<ElementType>::AssembleIPs(element, {2, 2, 2}, 1e-8);

                // Check if residual is smaller than targeted.
                QuESo_CHECK_LT(residual, 1e-8);

                // Must be more points than p*p*p.
                auto& r_points = element.GetIntegrationPoints();
                QuESo_CHECK_LT(r_points.size(), 28);
                QuESo_CHECK_GT(r_points.size(), 7);

                // Get copy of points.
                ElementType::IntegrationPointVectorType copy_points(r_points);

                // Compute constant terms.
                std::vector<double> constant_terms{};
                auto p_boundary_ips = element.pGetTrimmedDomain()->pGetBoundaryIps<BoundaryIntegrationPointType>();
                QuadratureTrimmedElementTester<ElementType>::ComputeConstantTerms(constant_terms, p_boundary_ips, element, {2, 2, 2});

                // Run moment fitting again.
                const auto residual_2 = QuadratureTrimmedElementTester<ElementType>::MomentFitting(constant_terms, r_points, element, {2, 2, 2});

                // Check if residual and weights are the same.
                QuESo_CHECK_LT( residual, residual_2+EPS4 );
                QuESo_CHECK_LT( residual_2, residual+EPS4 );
                double volume = 0.0;
                for( IndexType i = 0; i < r_points.size(); ++i){
                    const double weight1 = r_points[i].Weight();
                    QuESo_CHECK_GT(weight1, EPS4);
                    const double weight2 = copy_points[i].Weight();
                    const double error = std::abs(weight1 - weight2)/ weight1;
                    QuESo_CHECK_LT( error , EPS2 );
                    volume += weight1*element.DetJ(); // Multiplied with det(J).
                }

                // Check if integration points contain correct volume;
                const auto& r_mesh = element.pGetTrimmedDomain()->GetTriangleMesh();
                const double ref_volume = MeshUtilities::Volume(r_mesh);
                const double volume_error = std::abs(volume - ref_volume);
                QuESo_CHECK_LT(volume_error, 1e-7);
            }
        }
    }
    QuESo_CHECK_EQUAL(number_trimmed_elements, 153);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso