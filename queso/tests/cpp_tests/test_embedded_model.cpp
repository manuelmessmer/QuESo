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
#include <numeric>      // std::accumulate
#include <memory>       //std::addressof
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"
#include "queso/containers/grid_indexer.hpp"
#include "queso/io/io_utilities.h"
#include "queso/embedded_model.h"

#include "queso/tests/cpp_tests/global_config.hpp"

namespace queso {
namespace Testing {

BOOST_AUTO_TEST_SUITE( EmbeddedModelTestSuite )

BOOST_AUTO_TEST_CASE(IntersectedElementTest) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Intersected Element" << std::endl;

    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    std::string filename = base_dir + "/data/cylinder.stl";

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::input_filename, filename);
    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::echo_level, 0u);
    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::write_output_to_file, false);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{0.0, 0.0, 0.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{2.0, 2.0, 1.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{0.0, 0.0, 0.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{4.0, 5.0, 3.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, Vector3i{1, 1, 1});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::polynomial_order, Vector3i{2, 2, 2});
    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::moment_fitting_residual, 1e-8);
    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, 5000u);
    r_settings[MainSettings::non_trimmed_quadrature_rule_settings].SetValue(NonTrimmedQuadratureRuleSettings::integration_method, IntegrationMethod::gauss);

    EmbeddedModel embedded_model = EmbeddedModel::Create(std::move(p_settings));
    embedded_model.CreateAllFromSettings();

    const auto& elements = embedded_model.GetElements();

    QuESo_CHECK_EQUAL(elements.size(), 1);

    const auto& points_reduced = (*elements.begin())->GetIntegrationPoints();
    QuESo_CHECK_LT(points_reduced.size(), 28);

    const auto& r_triangle_mesh = (*elements.begin())->pGetTrimmedDomain()->GetTriangleMesh();
    const IndexType num_triangles = r_triangle_mesh.NumOfTriangles();

    QuESo_CHECK_GT(num_triangles, 5000.0);

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

    Vector3i num_elements = (Large) ? Vector3i{14, 22, 12} : Vector3i{7, 11, 6};
    GridTypeType grid_type = (BSplineMesh) ? GridType::b_spline_grid : GridType::hexahedral_fe_grid;
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    std::string filename = base_dir + "/data/elephant.stl";

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::input_filename, filename);
    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::echo_level, 0u);
    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::write_output_to_file, false);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, grid_type);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-0.37, -0.55, -0.31});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{0.37, 0.55, 0.31});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, rBoundsUVW.first);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, rBoundsUVW.second);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, num_elements);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::polynomial_order, rOrder);
    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::moment_fitting_residual, 1e-6);
    r_settings[MainSettings::non_trimmed_quadrature_rule_settings].SetValue(NonTrimmedQuadratureRuleSettings::integration_method, IntegrationMethod);

    EmbeddedModel embedded_model = EmbeddedModel::Create(std::move(p_settings));
    embedded_model.CreateAllFromSettings();

    const auto& elements = embedded_model.GetElements();

    // Compute total volume
    double volume_trimmed = 0.0;
    double volume_inside = 0.0;
    const auto el_it_begin = elements.begin();
    int num_elements_inside = 0;
    int num_elements_trimmed = 0;
    int num_points_inside = 0;
    for( IndexType i = 0; i < elements.size(); ++i){
        const auto& el_it = *(el_it_begin+i);
        const double det_j = el_it->DetJ();
        if( el_it->IsTrimmed() ){
            const auto& points_trimmed = el_it->GetIntegrationPoints();
            QuESo_CHECK_GT(points_trimmed.size(), 0);
            QuESo_CHECK_LT(points_trimmed.size(), (rOrder[0]+1)*(rOrder[1]+1)*(rOrder[2]+1)+1);
            for( const auto& point : points_trimmed ){
                volume_trimmed += point.Weight()*det_j;
            }
            num_elements_trimmed++;
        } else {
            const auto& points_inside = el_it->GetIntegrationPoints();
            for( const auto& point : points_inside ){
                volume_inside += point.Weight()*det_j;
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
    const PointType el_delta = Math::Subtract( el_bounding_box.second, el_bounding_box.first );
    const double ref_volume_inside = (el_delta[0]*el_delta[1]*el_delta[2])*num_elements_inside;
    QuESo_CHECK_RELATIVE_NEAR(volume_inside, ref_volume_inside, 1e-13);

    // Check volume inside + trimmed
    const double volume_tot = (volume_trimmed+volume_inside);
    const double ref_volume_tot = 0.0462012;
    QuESo_CHECK_RELATIVE_NEAR(volume_tot, ref_volume_tot, Tolerance);
}

// p=2
BOOST_AUTO_TEST_CASE(VolumeElephant1Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant Gauss (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::gauss, Vector3i{2, 2, 2}, 2916, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::gauss, Vector3i{2, 2, 2}, 2916, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
    TestElephant(IntegrationMethod::gauss, Vector3i{2, 2, 2}, 2916, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant2Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant Optimal (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::ggq_optimal, Vector3i{2, 2, 2}, 1786, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::ggq_optimal, Vector3i{2, 2, 2}, 1786, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant3Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant GGQ_Reduced1 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::ggq_reduced_1, Vector3i{2, 2, 2}, 673, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::ggq_reduced_1, Vector3i{2, 2, 2}, 673, 0.0002, use_b_spline_mesh, MakeBox({-6.0, -7.5, -1.2}, {22.0, 1.82, 2.8}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant4Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant GGQ_Reduced2 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::ggq_reduced_2, Vector3i{2, 2, 2}, 406, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::ggq_reduced_2, Vector3i{2, 2, 2}, 406, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), true);
}

// p=3
BOOST_AUTO_TEST_CASE(VolumeElephant5Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant Gauss (p=3)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::gauss, Vector3i{3, 3, 3}, 320, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), false);
    TestElephant(IntegrationMethod::gauss, Vector3i{3, 3, 3}, 320, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), false);
    TestElephant(IntegrationMethod::gauss, Vector3i{3, 3, 3}, 320, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), false);
}

// p=3
BOOST_AUTO_TEST_CASE(VolumeElephant6Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant GGQ_Optimal (p=3)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::ggq_optimal, Vector3i{3, 3, 3}, 256, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), false);
    TestElephant(IntegrationMethod::ggq_optimal, Vector3i{3, 3, 3}, 256, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), false);
}

// p=4
BOOST_AUTO_TEST_CASE(VolumeElephant7Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant Gauss (p=4)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::gauss, Vector3i{4, 4, 4}, 625, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), false);
    TestElephant(IntegrationMethod::gauss, Vector3i{4, 4, 4}, 625, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), false);
    TestElephant(IntegrationMethod::gauss, Vector3i{4, 4, 4}, 625, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), false);
}

// p=4
BOOST_AUTO_TEST_CASE(VolumeElephant8Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant GGQ_Optimal (p=4)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::ggq_optimal, Vector3i{4, 4, 4}, 525, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), false);
}

// p=Mix
BOOST_AUTO_TEST_CASE(VolumeElephant9Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant Gauss (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::gauss, Vector3i{2, 3, 4}, 6480, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::gauss, Vector3i{3, 4, 2}, 6480, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
    TestElephant(IntegrationMethod::gauss, Vector3i{4, 2, 3}, 6480, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant10Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant Gauss_Reduced1 (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::gauss_reduced_1, Vector3i{2, 3, 4}, 2592, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::gauss_reduced_1, Vector3i{3, 4, 2}, 2592, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
    TestElephant(IntegrationMethod::gauss_reduced_1, Vector3i{4, 2, 3}, 2592, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant11Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant Gauss_Reduced2 (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::gauss_reduced_2, Vector3i{2, 3, 4}, 648, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::gauss_reduced_2, Vector3i{3, 4, 2}, 648, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
    TestElephant(IntegrationMethod::gauss_reduced_2, Vector3i{4, 2, 3}, 648, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant12Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant Optimal (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::ggq_optimal, Vector3i{2, 3, 4}, 3657, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::ggq_optimal, Vector3i{3, 4, 2}, 3682, 0.0002, use_b_spline_mesh, MakeBox({-1.0, -5.5, -2.2}, {44.0, 1.12, 2.0}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant13Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant GGQ_Reduced1 (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::ggq_reduced_1, Vector3i{2, 3, 4}, 1710, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::ggq_reduced_1, Vector3i{3, 4, 2}, 1713, 0.0002, use_b_spline_mesh, MakeBox({-6.0, -7.5, -1.2}, {22.0, 1.82, 2.8}), true);
}

BOOST_AUTO_TEST_CASE(VolumeElephant14Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume :: Elephant GGQ_Reduced2 (p=Mix)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestElephant(IntegrationMethod::ggq_reduced_2, Vector3i{2, 3, 4}, 1166, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestElephant(IntegrationMethod::ggq_reduced_2, Vector3i{3, 4, 2}, 1192, 0.0002, use_b_spline_mesh, MakeBox({-0.37, -0.55, -0.31}, {0.37, 0.55, 0.31}), true);
}

void TestSteeringKnuckle( IntegrationMethodType IntegrationMethod, IndexType p, IndexType NumPointsInside, double Tolerance,
                        bool BSplineMesh, const BoundingBoxType& rBoundsUVW, bool Large){

    Vector3i num_elements = (Large) ? Vector3i{20, 40, 40} : Vector3i{5, 10, 10};
    GridTypeType grid_type = (BSplineMesh) ? GridType::b_spline_grid : GridType::hexahedral_fe_grid;
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    std::string filename = base_dir + "/data/steering_knuckle.stl";

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::input_filename, filename);
    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::echo_level, 0u);
    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::write_output_to_file, false);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, grid_type);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-130.0, -110.0, -110.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{20.0, 190.0, 190.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, rBoundsUVW.first);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, rBoundsUVW.second);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, num_elements);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::polynomial_order, Vector3i{p, p, p});
    r_settings[MainSettings::non_trimmed_quadrature_rule_settings].SetValue(NonTrimmedQuadratureRuleSettings::integration_method, IntegrationMethod);

    auto& r_cond_settings_list = r_settings.GetList(MainSettings::conditions_settings_list);

    auto p_cond_settings_1 = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
    p_cond_settings_1->SetValue(ConditionSettings::condition_id, 1u);
    p_cond_settings_1->SetValue(ConditionSettings::input_filename, base_dir + "/data/steering_knuckle_D1.stl");
    p_cond_settings_1->SetValue(ConditionSettings::condition_type, std::string("PenaltySupportCondition"));
    r_cond_settings_list.push_back(std::move(p_cond_settings_1));

    auto p_cond_settings_2 = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
    p_cond_settings_2->SetValue(ConditionSettings::condition_id, 2u);
    p_cond_settings_2->SetValue(ConditionSettings::input_filename, base_dir + "/data/steering_knuckle_N1.stl");
    p_cond_settings_2->SetValue(ConditionSettings::condition_type, std::string("SurfaceLoadCondition"));
    r_cond_settings_list.push_back(std::move(p_cond_settings_2));

    auto p_cond_settings_3 = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
    p_cond_settings_3->SetValue(ConditionSettings::condition_id, 3u);
    p_cond_settings_3->SetValue(ConditionSettings::input_filename, base_dir + "/data/steering_knuckle_N2.stl");
    p_cond_settings_3->SetValue(ConditionSettings::condition_type, std::string("SurfaceLoadCondition"));
    r_cond_settings_list.push_back(std::move(p_cond_settings_3));

    auto p_cond_settings_4 = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
    p_cond_settings_4->SetValue(ConditionSettings::condition_id, 4u);
    p_cond_settings_4->SetValue(ConditionSettings::input_filename, base_dir + "/data/steering_knuckle_N3.stl");
    p_cond_settings_4->SetValue(ConditionSettings::condition_type, std::string("SurfaceLoadCondition"));
    r_cond_settings_list.push_back(std::move(p_cond_settings_4));

    EmbeddedModel embedded_model = EmbeddedModel::Create(std::move(p_settings));
    embedded_model.CreateAllFromSettings();

    /// Test volume
    const auto& elements = embedded_model.GetElements();

    // Compute total volume
    double volume_trimmed = 0.0;
    double volume_inside = 0.0;
    const auto el_it_begin = elements.begin();
    int num_elements_inside = 0;
    int num_elements_trimmed = 0;
    int num_points_inside = 0;
    for( IndexType i = 0; i < elements.size(); ++i){
        const auto& el_it = *(el_it_begin+i);
        const double det_j = el_it->DetJ();
        if( el_it->IsTrimmed() ){
            const auto& points_trimmed = el_it->GetIntegrationPoints();
            QuESo_CHECK_GT(points_trimmed.size(), 0);
            QuESo_CHECK_LT(points_trimmed.size(), (p+1)*(p+1)*(p+1)+1);
            for( const auto& point : points_trimmed ){
                volume_trimmed += point.Weight()*det_j;
            }
            num_elements_trimmed++;
        } else {
            const auto& points_inside = el_it->GetIntegrationPoints();
            for( const auto& point : points_inside ){
                volume_inside += point.Weight()*det_j;
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
    const PointType delta = Math::Subtract( el_bounding_box.second, el_bounding_box.first );
    const double ref_volume_inside = (delta[0]*delta[1]*delta[2])*num_elements_inside;
    QuESo_CHECK_RELATIVE_NEAR(volume_inside, ref_volume_inside, 1e-13);

    // Check volume inside + trimmed
    const double volume_tot = (volume_trimmed+volume_inside);
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, filename.c_str());
    const double ref_volume_tot = MeshUtilities::Volume(triangle_mesh);
    QuESo_CHECK_RELATIVE_NEAR(volume_tot, ref_volume_tot, Tolerance);

    /// Test conditions
    const auto& conditions = embedded_model.GetConditions();
    QuESo_CHECK_EQUAL(conditions.size(), 4);
    for( const auto& p_condition : conditions ){
        TriangleMesh triangle_mesh{};
        const auto& r_cond_setting = p_condition->GetSettings();
        std::string filename = r_cond_setting.GetValue<std::string>(ConditionSettings::input_filename);
        IO::ReadMeshFromSTL(triangle_mesh, filename);
        const double ref_area = MeshUtilities::Area(triangle_mesh);
        auto lambda = [](double result, const auto& r_segment){return result + MeshUtilities::Area(r_segment.GetTriangleMesh()); };
        const double area = std::accumulate(p_condition->SegmentsBegin(), p_condition->SegmentsEnd(), 0.0, lambda);
        QuESo_CHECK_NEAR(ref_area, area, 1e-10);
    }
}

// p=2
BOOST_AUTO_TEST_CASE(SteeringKnuckle1Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume And Conditions :: Steering Knuckle Gauss (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::gauss, 2, 5427, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss, 2, 5427, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss, 2, 5427, 0.0002, use_b_spline_mesh, MakeBox({-12.0, -1110.0, -2110.0}, {220.0, 3190.0, 3190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss, 2, 5427, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss, 2, 5427, 0.0002, false, MakeBox({-2, -5.0, -10.0}, {3.0, 19.0, 10.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle2Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume And Conditions :: Steering Knuckle Gauss_Reduced1 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_1, 2, 1608, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_1, 2, 1608, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_1, 2, 1608, 0.0002, use_b_spline_mesh, MakeBox({-12.0, -1110.0, -2110.0}, {220.0, 3190.0, 3190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_1, 2, 1608, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_1, 2, 1608, 0.0002, false, MakeBox({-2, -5.0, -10.0}, {3.0, 19.0, 10.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle3Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume And Conditions :: Steering Knuckle Gauss_Reduced2 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_2, 2, 201, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_2, 2, 201, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_2, 2, 201, 0.0002, use_b_spline_mesh, MakeBox({-12.0, -1110.0, -2110.0}, {220.0, 3190.0, 3190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_2, 2, 201, 0.0002, false, MakeBox({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::gauss_reduced_2, 2, 201, 0.0002, false, MakeBox({-2, -5.0, -10.0}, {3.0, 19.0, 10.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle4Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume And Conditions :: Steering Knuckle Optimal (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::ggq_optimal, 2, 3483, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::ggq_optimal, 2, 3483, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
    TestSteeringKnuckle(IntegrationMethod::ggq_optimal, 2, 3483, 0.0002, use_b_spline_mesh, MakeBox({-2.0, -5.0, -10.0}, {3.0, 19.0, 10.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle5Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume And Conditions :: Steering Knuckle GGQ_Reduced1 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::ggq_reduced_1, 2, 1320, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::ggq_reduced_1, 2, 1320, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
}

BOOST_AUTO_TEST_CASE(SteeringKnuckle6Test) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Create Volume And Conditions :: Volume Elephant GGQ_Reduced2 (p=2)" << std::endl;
    const bool use_b_spline_mesh = true;
    TestSteeringKnuckle(IntegrationMethod::ggq_reduced_2, 2, 804, 0.0002, use_b_spline_mesh, MakeBox({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), true);
    TestSteeringKnuckle(IntegrationMethod::ggq_reduced_2, 2, 804, 0.0002, use_b_spline_mesh, MakeBox({-130.0, -110.0, -110.0}, {20.0, 190.0, 190.0}), true);
}


BOOST_AUTO_TEST_CASE(SteeringKnuckleModelInfoTest) {
    QuESo_INFO << "Testing :: Test Embedded Model :: Model Info :: Steering Knuckle" << std::endl;

    Vector3i num_elements = Vector3i{10, 20, 20};
    GridTypeType grid_type = GridType::b_spline_grid;
    std::string base_dir = GlobalConfig::GetInstance().BaseDir;
    std::string filename = base_dir + "/data/steering_knuckle.stl";

    auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& r_settings = *p_settings;

    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::input_filename, filename);
    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::echo_level, 0u);
    r_settings[MainSettings::general_settings].SetValue(GeneralSettings::write_output_to_file, false);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, grid_type);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{-130.0, -110.0, -110.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{20.0, 190.0, 190.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{0.0, 0.0, 0.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 1.0, 1.0});
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, num_elements);
    r_settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::polynomial_order, Vector3i{2, 2, 2});
    r_settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio,  0.0);
    r_settings[MainSettings::non_trimmed_quadrature_rule_settings].SetValue(NonTrimmedQuadratureRuleSettings::integration_method,  IntegrationMethod::ggq_optimal);

    auto& r_cond_setttings_list = r_settings.GetList(MainSettings::conditions_settings_list);

    auto p_cond_settings_1 = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
    p_cond_settings_1->SetValue(ConditionSettings::condition_id, 1u);
    p_cond_settings_1->SetValue(ConditionSettings::input_filename, base_dir + "/data/steering_knuckle_D1.stl");
    p_cond_settings_1->SetValue(ConditionSettings::condition_type, std::string("PenaltySupportCondition"));
    r_cond_setttings_list.push_back(std::move(p_cond_settings_1));

    EmbeddedModel embedded_model = EmbeddedModel::Create(std::move(p_settings));
    embedded_model.CreateAllFromSettings();

    /// Check model info
    const auto& r_model_info = embedded_model.GetModelInfo();
    // embedded_geometry_info
    const auto& r_geo_info = r_model_info[MainInfo::embedded_geometry_info];
    QuESo_CHECK(r_geo_info.GetValue<bool>(EmbeddedGeometryInfo::is_closed));
    TriangleMesh triangle_mesh{};
    IO::ReadMeshFromSTL(triangle_mesh, filename);
    double volume_ref = MeshUtilities::VolumeOMP(triangle_mesh);
    QuESo_CHECK_NEAR(r_geo_info.GetValue<double>(EmbeddedGeometryInfo::volume), volume_ref, 1e-10);
    // quadrature_info
    const auto& r_quad_info = r_model_info[MainInfo::quadrature_info];
    QuESo_CHECK_RELATIVE_NEAR( r_quad_info.GetValue<double>(QuadratureInfo::represented_volume), volume_ref, 1e-5)
    QuESo_CHECK_RELATIVE_NEAR( r_quad_info.GetValue<double>(QuadratureInfo::percentage_of_geometry_volume), 100.0, 1e-5)
    QuESo_CHECK_EQUAL(r_quad_info.GetValue<IndexType>(QuadratureInfo::tot_num_points), 9508);
    QuESo_CHECK_RELATIVE_NEAR( r_quad_info.GetValue<double>(QuadratureInfo::num_of_points_per_full_element), 25.2, 1e-5)
    const double num_of_points_per_trimmed_element = r_quad_info.GetValue<double>(QuadratureInfo::num_of_points_per_trimmed_element);
    QuESo_CHECK_GT(num_of_points_per_trimmed_element, 26);
    QuESo_CHECK_LT(num_of_points_per_trimmed_element, 27);
    // background_grid_info
    const auto& r_grid_info = r_model_info[MainInfo::background_grid_info];
    QuESo_CHECK_EQUAL(r_grid_info.GetValue<IndexType>(BackgroundGridInfo::num_active_elements), 354);
    QuESo_CHECK_EQUAL(r_grid_info.GetValue<IndexType>(BackgroundGridInfo::num_trimmed_elements), 349);
    QuESo_CHECK_EQUAL(r_grid_info.GetValue<IndexType>(BackgroundGridInfo::num_full_elements), 5);
    QuESo_CHECK_EQUAL(r_grid_info.GetValue<IndexType>(BackgroundGridInfo::num_inactive_elements), 3646);
    // elapsed_time_info
    const auto& r_elapsed_time_info = r_model_info[MainInfo::elapsed_time_info];
    const double total_time = r_elapsed_time_info.GetValue<double>(ElapsedTimeInfo::total);

    const auto& r_volume_time_info = r_elapsed_time_info[ElapsedTimeInfo::volume_time_info];
    const double et_volume_total = r_volume_time_info.GetValue<double>(VolumeTimeInfo::total);
    const double et_coe = r_volume_time_info.GetValue<double>(VolumeTimeInfo::classification_of_elements);
    QuESo_CHECK_GT(et_coe, EPS1);
    const double et_coi = r_volume_time_info.GetValue<double>(VolumeTimeInfo::computation_of_intersections);
    QuESo_CHECK_GT(et_coi, EPS1);
    const double et_somf = r_volume_time_info.GetValue<double>(VolumeTimeInfo::solution_of_moment_fitting_eqs);
    QuESo_CHECK_GT(et_somf, EPS1);
    const double et_cogr = r_volume_time_info.GetValue<double>(VolumeTimeInfo::construction_of_ggq_rules);
    QuESo_CHECK_GT(et_cogr, EPS1);
    QuESo_CHECK_GT( et_volume_total, (et_coe+et_coi+et_somf+et_cogr) );

    const auto& r_condition_time_info = r_elapsed_time_info[ElapsedTimeInfo::conditions_time_info];
    const double et_condition_total = r_condition_time_info.GetValue<double>(ConditionsTimeInfo::total);
    QuESo_CHECK_GT(et_condition_total, EPS1);

    const auto& r_write_file_time_info = r_elapsed_time_info[ElapsedTimeInfo::write_files_time_info];
    const double et_write_file_total = r_write_file_time_info.GetValue<double>(WriteFilesTimeInfo::total);
    QuESo_CHECK_NEAR(et_write_file_total, 0.0, 1e-10);

    QuESo_CHECK_RELATIVE_NEAR(total_time, (et_volume_total+et_condition_total+et_write_file_total), 1e-5);

    const auto& r_condition_info_1 = *r_model_info.GetList(MainInfo::conditions_infos_list)[0];
    const auto& r_condition_1 = embedded_model.GetConditions()[0];
    const auto& r_condition_info_1_other = r_condition_1->GetInfo();
    QuESo_CHECK_EQUAL(std::addressof(r_condition_info_1), std::addressof(r_condition_info_1_other));
    QuESo_CHECK_EQUAL( r_condition_info_1.GetValue<IndexType>(ConditionInfo::condition_id), 1);

    const auto& r_settings_obtained = embedded_model.GetSettings();
    const auto& r_cond_settings_1 = *r_settings_obtained.GetList(MainSettings::conditions_settings_list)[0];
    const std::string& r_filename = r_cond_settings_1.GetValue<std::string>(ConditionSettings::input_filename);
    TriangleMesh triangle_mesh_cond{};
    IO::ReadMeshFromSTL(triangle_mesh_cond, r_filename);
    const double area_ref = MeshUtilities::Area(triangle_mesh_cond);
    QuESo_CHECK_NEAR( r_condition_info_1.GetValue<double>(ConditionInfo::surf_area), area_ref, 1e-10);
    QuESo_CHECK_NEAR( r_condition_info_1.GetValue<double>(ConditionInfo::perc_surf_area_in_active_domain), 100.0, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso
