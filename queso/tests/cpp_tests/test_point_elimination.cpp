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

//// External includes
#include <boost/test/unit_test.hpp>

//// STL includes
#include <string_view>

//// Project includes
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/containers/trimmed_element.hpp"
#include "queso/embedding/domain_mesh_embedder.h"
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/io/io_utilities.h"
#include "queso/quadrature/moment_fitting.hpp"
#include "queso/utilities/mesh_utilities.h"

#include "queso/tests/cpp_tests/global_config.hpp"

// This suite tests point elimination on real trimmed geometries and verifies that returned reduced rules have stable,
// current weights. Related coverage: test_moment_fitting_assembly.cpp checks assembly details, test_moment_fitting.cpp
// checks unreduced NNLS solves, and test_element_trimmed.cpp checks public Compute API behavior.

namespace queso::Testing {

namespace {

    using IntegrationPointType = IntegrationPoint;
    using BoundaryIntegrationPointType = BoundaryIntegrationPoint;
    using ElementType = TrimmedElement<IntegrationPointType, BoundaryIntegrationPointType>;

    struct PointEliminationCase
    {
        std::string_view stl_filename;
        GridType grid_type;
        PointType lower_bound_xyz;
        PointType upper_bound_xyz;
        PointType lower_bound_uvw;
        PointType upper_bound_uvw;
        Vector3i number_of_elements;
        Vector3i integration_order;
        double target_residual;
        double residual_tolerance;
        double volume_tolerance;
        IndexType max_number_of_points;
        IndexType expected_trimmed_elements;
    };

    Dictionary<queso::key::MainValuesTypeTag> MakeSettings(const PointEliminationCase& rCase)
    {
        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_settings = *p_settings;

        auto& r_grid_settings = r_settings[MainSettings::background_grid_settings];
        r_grid_settings.SetValue(BackgroundGridSettings::grid_type, rCase.grid_type);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, rCase.lower_bound_xyz);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, rCase.upper_bound_xyz);
        r_grid_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, rCase.lower_bound_uvw);
        r_grid_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, rCase.upper_bound_uvw);
        r_grid_settings.SetValue(BackgroundGridSettings::polynomial_order, rCase.integration_order);
        r_grid_settings.SetValue(BackgroundGridSettings::number_of_elements, rCase.number_of_elements);
        r_grid_settings.CheckRequired();

        return std::move(r_settings);
    }

    double RefitReturnedRuleAndCheckWeightsAreUnchanged(
        ElementType& rElement,
        ElementType::IntegrationPointVectorType& rPoints,
        const ElementType::IntegrationPointVectorType& rReferencePoints,
        const Vector3i& rIntegrationOrder
    )
    {
        auto boundary_ips =
            rElement.GetActiveDomainBoundaryIps<BoundaryIntegrationPointType, CoordinateSpace::global>();
        const auto order_info = quadrature::moment_fitting::detail::MakeIntegrationOrderInfo(rIntegrationOrder);
        auto constant_terms = quadrature::moment_fitting::detail::ComputeConstantTerms(
            boundary_ips, rElement.GetCellBounds<CoordinateSpace::global>(), order_info
        );

        const quadrature::moment_fitting::detail::MomentFittingProblem problem(constant_terms);
        quadrature::moment_fitting::detail::MomentFittingScratch scratch{};
        const quadrature::moment_fitting::detail::IntegrationGeometry geometry{
            rElement.GetCellBounds<CoordinateSpace::parametric>(), rElement.DetJ()
        };
        const auto fitting_result = quadrature::moment_fitting::detail::MomentFitting<ElementType>(
            rPoints, geometry, order_info, problem, scratch
        );

        for (IndexType i = 0; i < rPoints.size(); ++i) {
            QuESo_CHECK_GT(rPoints[i].Weight(), EPS4);
            QuESo_CHECK_RELATIVE_NEAR(rPoints[i].Weight(), rReferencePoints[i].Weight(), EPS2);
        }

        return fitting_result.residual;
    }

    double IntegratedVolume(const ElementType& rElement)
    {
        double volume = 0.0;
        for (const auto& rPoint : rElement.GetIntegrationPoints()) { volume += rPoint.Weight() * rElement.DetJ(); }
        return volume;
    }

    void CheckVolume(const ElementType& rElement, double Tolerance)
    {
        const double volume = IntegratedVolume(rElement);
        const double reference_volume = MeshUtilities::Volume(rElement.GetActiveDomainBoundaryMesh());
        QuESo_CHECK_RELATIVE_NEAR(volume, reference_volume, Tolerance);
    }

    void CheckReturnedRule(ElementType& rElement, const PointEliminationCase& rCase)
    {
        const auto residual =
            quadrature::moment_fitting::Compute(
                rElement, { .integration_order = rCase.integration_order, .residual = rCase.target_residual }
            )
                .value();
        QuESo_CHECK_NEAR(residual, 0.0, rCase.residual_tolerance);

        auto& r_points = rElement.GetIntegrationPoints();
        QuESo_CHECK_LT(r_points.size(), rCase.max_number_of_points);
        QuESo_CHECK_IS_FALSE(r_points.empty());

        const ElementType::IntegrationPointVectorType reference_points(r_points);
        // Refit the returned point set to verify point elimination did not return stale weights.
        const double residual_after_refitting_returned_rule =
            RefitReturnedRuleAndCheckWeightsAreUnchanged(rElement, r_points, reference_points, rCase.integration_order);
        QuESo_CHECK_NEAR(residual, residual_after_refitting_returned_rule, EPS4);

        CheckVolume(rElement, rCase.volume_tolerance);
    }

    void RunPointEliminationCase(const PointEliminationCase& rCase)
    {
        const auto settings = MakeSettings(rCase);
        TriangleMesh triangle_mesh{};
        const std::string stl_path = GlobalConfig::GetInstance().BaseDir + "/data/" + std::string(rCase.stl_filename);
        IO::ReadMeshFromSTL(triangle_mesh, stl_path);

        constexpr double min_vol_ratio = 1e-3;
        constexpr IndexType min_num_triangles = 500;

        GridIndexer grid_indexer(settings);
        embedding::DomainMeshEmbedder domain_embedder(triangle_mesh.View(), grid_indexer);
        const auto states = domain_embedder.Classify();
        IndexType number_trimmed_elements = 0;
        for (IndexType i = 0; i < grid_indexer.NumberOfElements(); ++i) {
            const BoundingBoxType bounding_box = grid_indexer.GetBoundingBoxXYZFromIndex(i);
            if (states[i] != IntersectionState::trimmed) { continue; }

            auto trimmed_domain = domain_embedder.MakeTrimmedDomain(i, min_num_triangles);
            const PointType delta = bounding_box.upper - bounding_box.lower;
            const double cell_volume = delta[0] * delta[1] * delta[2];
            if (MeshUtilities::Volume(trimmed_domain.GetBoundaryMesh()) / cell_volume <= min_vol_ratio) { continue; }

            ++number_trimmed_elements;
            const BoundingBoxType bounding_box_uvw = grid_indexer.GetBoundingBoxUVWFromIndex(i);
            ElementType element(1, ElementBounds{ bounding_box, bounding_box_uvw }, std::move(trimmed_domain));
            CheckReturnedRule(element, rCase);
        }
        QuESo_CHECK_EQUAL(number_trimmed_elements, rCase.expected_trimmed_elements);
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(PointEliminationTestSuite)

BOOST_AUTO_TEST_CASE(PointEliminationCylinder1Test)
{
    // Verifies reduced quadratic rules on trimmed cylinder cells and re-solve stability of returned weights.
    QuESo_INFO << "Testing :: Test Point Elimination :: Cylinder Quadratic" << std::endl;
    constexpr PointEliminationCase test_case{ .stl_filename = "cylinder.stl",
                                              .grid_type = GridType::b_spline_grid,
                                              .lower_bound_xyz = { -1.5, -1.5, -1.0 },
                                              .upper_bound_xyz = { 1.5, 1.5, 12.0 },
                                              .lower_bound_uvw = { 0.0, 0.0, 0.0 },
                                              .upper_bound_uvw = { 1.0, 1.0, 1.0 },
                                              .number_of_elements = { 6, 6, 13 },
                                              .integration_order = { 2, 2, 2 },
                                              .target_residual = 1e-8,
                                              .residual_tolerance = 1e-6,
                                              .volume_tolerance = 1e-6,
                                              .max_number_of_points = 28,
                                              .expected_trimmed_elements = 120 };
    RunPointEliminationCase(test_case);
}

BOOST_AUTO_TEST_CASE(PointEliminationCylinder2Test)
{
    // Verifies reduced cubic rules on trimmed cylinder cells and re-solve stability of returned weights.
    QuESo_INFO << "Testing :: Test Point Elimination :: Cylinder Cubic" << std::endl;
    constexpr PointEliminationCase test_case{ .stl_filename = "cylinder.stl",
                                              .grid_type = GridType::b_spline_grid,
                                              .lower_bound_xyz = { -1.5, -1.5, -1.0 },
                                              .upper_bound_xyz = { 1.5, 1.5, 12.0 },
                                              .lower_bound_uvw = { 0.0, 0.0, 0.0 },
                                              .upper_bound_uvw = { 1.0, 1.0, 1.0 },
                                              .number_of_elements = { 6, 6, 13 },
                                              .integration_order = { 3, 3, 3 },
                                              .target_residual = 1e-8,
                                              .residual_tolerance = 1e-6,
                                              .volume_tolerance = 1e-6,
                                              .max_number_of_points = 65,
                                              .expected_trimmed_elements = 120 };
    RunPointEliminationCase(test_case);
}

BOOST_AUTO_TEST_CASE(PointEliminationCylinder4Test)
{
    // Verifies reduced mixed-order rules on trimmed cylinder cells and re-solve stability of returned weights.
    QuESo_INFO << "Testing :: Test Point Elimination :: Cylinder Mixed" << std::endl;
    constexpr PointEliminationCase test_case{ .stl_filename = "cylinder.stl",
                                              .grid_type = GridType::b_spline_grid,
                                              .lower_bound_xyz = { -1.5, -1.5, -1.0 },
                                              .upper_bound_xyz = { 1.5, 1.5, 12.0 },
                                              .lower_bound_uvw = { 0.0, 0.0, 0.0 },
                                              .upper_bound_uvw = { 1.0, 1.0, 1.0 },
                                              .number_of_elements = { 6, 6, 13 },
                                              .integration_order = { 2, 3, 4 },
                                              .target_residual = 1e-7,
                                              .residual_tolerance = 1e-6,
                                              .volume_tolerance = 1e-5,
                                              .max_number_of_points = 61,
                                              .expected_trimmed_elements = 120 };
    RunPointEliminationCase(test_case);
}

BOOST_AUTO_TEST_CASE(PointEliminationKnuckleTest)
{
    // Verifies reduced rules on steering-knuckle geometry, including volume recovery and stable returned weights.
    QuESo_INFO << "Testing :: Test Point Elimination :: Knuckle" << std::endl;
    constexpr PointEliminationCase test_case{ .stl_filename = "steering_knuckle.stl",
                                              .grid_type = GridType::b_spline_grid,
                                              .lower_bound_xyz = { -130.0, -110.0, -110.0 },
                                              .upper_bound_xyz = { -40.0, 10.0, 10.0 },
                                              .lower_bound_uvw = { -130.0, -110.0, -110.0 },
                                              .upper_bound_uvw = { -40.0, 10.0, 10.0 },
                                              .number_of_elements = { 9, 12, 12 },
                                              .integration_order = { 2, 2, 2 },
                                              .target_residual = 1e-8,
                                              .residual_tolerance = 1e-8,
                                              .volume_tolerance = 1e-6,
                                              .max_number_of_points = 28,
                                              .expected_trimmed_elements = 80 };
    RunPointEliminationCase(test_case);
}

BOOST_AUTO_TEST_CASE(PointEliminationElephantTest)
{
    // Verifies reduced rules on elephant geometry, including volume recovery and stable returned weights.
    QuESo_INFO << "Testing :: Test Point Elimination :: Elephant" << std::endl;
    constexpr PointEliminationCase test_case{ .stl_filename = "elephant.stl",
                                              .grid_type = GridType::hexahedral_fe_grid,
                                              .lower_bound_xyz = { -0.4, -0.6, -0.35 },
                                              .upper_bound_xyz = { 0.4, 0.6, 0.35 },
                                              .lower_bound_uvw = { -1.0, -1.0, -1.0 },
                                              .upper_bound_uvw = { 1.0, 1.0, 1.0 },
                                              .number_of_elements = { 8, 12, 7 },
                                              .integration_order = { 2, 2, 2 },
                                              .target_residual = 1e-8,
                                              .residual_tolerance = 1e-8,
                                              .volume_tolerance = 1e-7,
                                              .max_number_of_points = 28,
                                              .expected_trimmed_elements = 153 };
    RunPointEliminationCase(test_case);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing
