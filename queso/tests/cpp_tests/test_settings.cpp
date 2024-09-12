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
#include <boost/test/execution_monitor.hpp>
//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/includes/settings.hpp"

namespace queso {


namespace Testing {

BOOST_AUTO_TEST_SUITE( SettingsTestSuite )

BOOST_AUTO_TEST_CASE(SettingsWrongTypeTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Wrong Types" << std::endl;
    Settings setting;

    BOOST_REQUIRE_THROW(setting[GeneralSettings::echo_level], std::exception); // Wrong Key type

    /// General settings
    auto& r_general_settings = setting[MainSettings::general_settings];
    BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception); // Wrong Key type

    BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::input_filename), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::output_directory_name), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::echo_level), std::exception); // Wrong Value type

    /// Mesh settings
    auto& r_mesh_settings = setting[MainSettings::background_grid_settings];
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<IndexType>(GeneralSettings::echo_level), std::exception); // Wrong Key type

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(GeneralSettings::echo_level, 2u), std::exception); // Wrong Key type
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{1.0, 2.0, 3.0});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(GeneralSettings::echo_level, IndexType(2)), std::exception); // Wrong Key type
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>(BackgroundGridSettings::lower_bound_xyz), std::exception); // Wrong Value type
    const PointType& r_lower_bound_xyz = r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz);
    QuESo_CHECK_POINT_NEAR(r_lower_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{1.0, 2.0, 3.0});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>(BackgroundGridSettings::upper_bound_xyz), std::exception); // Wrong Value type
    const PointType& r_upper_bound_xyz = r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz);
    QuESo_CHECK_POINT_NEAR(r_upper_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{1.0, 2.0, 3.0});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>(BackgroundGridSettings::lower_bound_uvw), std::exception); // Wrong Value type
    const PointType& r_lower_bound_uvw= r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw);
    QuESo_CHECK_POINT_NEAR(r_lower_bound_uvw, PointType({1.0, 2.0, 3.0}), 1e-10);

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 2.0, 1.0});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>(BackgroundGridSettings::upper_bound_uvw), std::exception); // Wrong Value type
    const PointType& r_upper_bound_uvw = r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw);
    QuESo_CHECK_Vector3i_EQUAL(r_upper_bound_uvw, Vector3i({1, 2, 3}) );

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::polynomial_order), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::polynomial_order, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{1, 2, 3});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::polynomial_order, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::polynomial_order), std::exception); // Wrong Value type
    const Vector3i& r_polynomial_order = r_mesh_settings.GetValue<Vector3i>(BackgroundGridSettings::polynomial_order);
    QuESo_CHECK_Vector3i_EQUAL(r_polynomial_order, Vector3i({1, 2, 3}) );

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::number_of_elements), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::number_of_elements, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{1, 2, 3});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::number_of_elements, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::number_of_elements), std::exception); // Wrong Value type
    const Vector3i& number_of_elements = r_mesh_settings.GetValue<Vector3i>(BackgroundGridSettings::number_of_elements);
    QuESo_CHECK_Vector3i_EQUAL(number_of_elements, Vector3i({1, 2, 3}) );

    /// Trimmed quadrature rule settings
    auto& r_trimmed_quad_rule_settings = setting[MainSettings::trimmed_quadrature_rule_settings];
    BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{1, 2, 3}), std::exception); // Wrong Key type
    BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>(TrimmedQuadratureRuleSettings::moment_fitting_residual), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), std::exception); // Wrong Value type

    /// Non trimmed quadrature rule settings
    auto& r_non_trimmed_quad_rule_settings = setting[MainSettings::non_trimmed_quadrature_rule_settings];
    BOOST_REQUIRE_THROW(r_non_trimmed_quad_rule_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{1, 2, 3}), std::exception); // Wrong Key type
    BOOST_REQUIRE_THROW(r_non_trimmed_quad_rule_settings.GetValue<IndexType>(NonTrimmedQuadratureRuleSettings::integration_method), std::exception); // Wrong Value type
}


BOOST_AUTO_TEST_CASE(SettingsDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Default Values" << std::endl;

    Settings settings;

    /// General settings
    QuESo_CHECK( !settings[MainSettings::general_settings].IsSet(GeneralSettings::input_filename) );
    BOOST_REQUIRE_THROW( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::input_filename), std::exception );

    QuESo_CHECK( settings[MainSettings::general_settings].IsSet(GeneralSettings::output_directory_name) );
    QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::output_directory_name), std::string("queso_output") );

    QuESo_CHECK( settings[MainSettings::general_settings].IsSet(GeneralSettings::echo_level) );
    QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level), 1u);

    /// Mesh settings
    QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::lower_bound_xyz) );
    BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz), std::exception );

    QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::upper_bound_xyz) );
    BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz), std::exception );

    QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::lower_bound_uvw) );
    BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception );

    QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::upper_bound_uvw) );
    BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw), std::exception );

    QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::polynomial_order) );
    BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::polynomial_order), std::exception );

    QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::number_of_elements) );
    BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::number_of_elements), std::exception );

    /// TrimmedQuadratureRuleSettings settings
    QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::moment_fitting_residual) );
    QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual), 1e-10,1e-10 );

    QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::min_element_volume_ratio) );
    QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 1e-3,1e-10 );

    QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::min_num_boundary_triangles) );
    QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 100 );

    // NonTrimmedQuadratureRuleSettings settings
    QuESo_CHECK( settings[MainSettings::non_trimmed_quadrature_rule_settings].IsSet(NonTrimmedQuadratureRuleSettings::integration_method) );
    QuESo_CHECK_EQUAL( settings[MainSettings::non_trimmed_quadrature_rule_settings].GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::Gauss );
}

BOOST_AUTO_TEST_CASE(SettingsCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Customized Values" << std::endl;

    Settings settings;

    /// General settings
    settings[MainSettings::general_settings].SetValue(GeneralSettings::input_filename, std::string("test_filename.stl"));
    settings[MainSettings::general_settings].SetValue(GeneralSettings::output_directory_name, std::string("new_output/"));
    settings[MainSettings::general_settings].SetValue(GeneralSettings::echo_level, 2u);


    QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::input_filename), std::string("test_filename.stl") );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::general_settings).GetValueNoCheck<std::string>(GeneralSettings::input_filename),
        std::string("test_filename.stl") );

    QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::output_directory_name), std::string("new_output/") );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::general_settings).GetValueNoCheck<std::string>(GeneralSettings::output_directory_name),
        std::string("new_output/") );

    QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level), 2u );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::general_settings).GetValueNoCheck<IndexType>(GeneralSettings::echo_level), 2UL );

    /// Mesh settings
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, PointType({1.0, 1.0, 2.0}));
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, PointType({0.0, 3.0, 2.0}));
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, PointType({0.1, -1.0, 2.0}));
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, PointType({0.66, 1.0, 2.2}));
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::polynomial_order, Vector3i({5, 6, 7}));
    settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, Vector3i({8, 9, 2}));

    QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz), PointType({1.0, 1.0, 2.0}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::background_grid_settings).GetValueNoCheck<PointType>(BackgroundGridSettings::lower_bound_xyz),
        PointType({1.0, 1.0, 2.0}) );

    QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz), PointType({0.0, 3.0, 2.0}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::background_grid_settings).GetValueNoCheck<PointType>(BackgroundGridSettings::upper_bound_xyz),
        PointType({0.0, 3.0, 2.0}) );

    QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), PointType({0.1, -1.0, 2.0}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::background_grid_settings).GetValueNoCheck<PointType>(BackgroundGridSettings::lower_bound_uvw),
        PointType({0.1, -1.0, 2.0}) );

    QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw), PointType({0.66, 1.0, 2.2}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::background_grid_settings).GetValueNoCheck<PointType>(BackgroundGridSettings::upper_bound_uvw),
        PointType({0.66, 1.0, 2.2}) );

    QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::polynomial_order), Vector3i({5, 6, 7}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::background_grid_settings).GetValueNoCheck<Vector3i>(BackgroundGridSettings::polynomial_order),
        Vector3i({5, 6, 7}) );

    QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::number_of_elements), Vector3i({8, 9, 2}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::background_grid_settings).GetValueNoCheck<Vector3i>(BackgroundGridSettings::number_of_elements),
        Vector3i({8, 9, 2}) );

    /// Trimmed quadrature rule settings
    settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::moment_fitting_residual, 5.6e-5);
    settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 0.45);
    settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, IndexType(234));

    QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual), 5.6e-5, 1e-10 );
    QuESo_CHECK_RELATIVE_NEAR( settings.GetSubDictionaryNoCheck(MainSettings::trimmed_quadrature_rule_settings).GetValueNoCheck<double>(
        TrimmedQuadratureRuleSettings::moment_fitting_residual), 5.6e-5, 1e-10 );

    QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 0.45, 1e-10 );
    QuESo_CHECK_RELATIVE_NEAR( settings.GetSubDictionaryNoCheck(MainSettings::trimmed_quadrature_rule_settings).GetValueNoCheck<double>(
        TrimmedQuadratureRuleSettings::min_element_volume_ratio), 0.45, 1e-10 );

    QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 234u );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::trimmed_quadrature_rule_settings).GetValueNoCheck<IndexType>(
        TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 234u);

    /// Non trimmed quadrature rule settings
    settings[MainSettings::non_trimmed_quadrature_rule_settings].SetValue(NonTrimmedQuadratureRuleSettings::integration_method, IntegrationMethod::GGQ_Optimal);

    QuESo_CHECK_EQUAL( settings[MainSettings::non_trimmed_quadrature_rule_settings].GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::GGQ_Optimal );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(MainSettings::non_trimmed_quadrature_rule_settings).GetValueNoCheck<IntegrationMethod>(
        NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::GGQ_Optimal);
}


BOOST_AUTO_TEST_CASE(SettingsConditionSettingsWrongTypeTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Condition Settings Wrong Types " << std::endl;

    Settings settings;

    auto& r_cond_settings_list = settings.CreateNewConditionSettings();
    BOOST_REQUIRE_THROW( r_cond_settings_list.IsSet(BackgroundGridSettings::lower_bound_uvw), std::exception ); // Wrong Key

    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(ConditionSettings::condition_id), std::exception ); // Not set
    BOOST_REQUIRE_THROW( r_cond_settings_list.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType({2.0, 3.0, 4.0})), std::exception ); // Wrong Key type
    BOOST_REQUIRE_THROW( r_cond_settings_list.SetValue(ConditionSettings::condition_id, 2.0), std::exception ); // Wrong Key type

    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(ConditionSettings::condition_id), std::exception ); // Not set
    r_cond_settings_list.SetValue(ConditionSettings::condition_id, 100u);
    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(GeneralSettings::echo_level), std::exception ); // Wrong Key type
    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<double>(ConditionSettings::condition_id), std::exception ); // Wrong Type

    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<std::string>(ConditionSettings::condition_type), std::exception ); // Not set
    r_cond_settings_list.SetValue(ConditionSettings::condition_type, std::string("dummy"));
    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(GeneralSettings::echo_level), std::exception ); // Wrong Key type
    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(ConditionSettings::condition_type), std::exception ); // Wrong Type

    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<std::string>(ConditionSettings::input_filename), std::exception ); // Not set
    r_cond_settings_list.SetValue(ConditionSettings::input_filename, std::string("dummy_2"));
    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(ConditionSettings::input_filename), std::exception ); // Wrong Type

    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<double>(ConditionSettings::modulus), std::exception ); // Not set
    r_cond_settings_list.SetValue(ConditionSettings::modulus, 2.0);
    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(ConditionSettings::modulus), std::exception ); // Wrong Type

    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<PointType>(ConditionSettings::direction), std::exception ); // Not set
    r_cond_settings_list.SetValue(ConditionSettings::direction, PointType({2.2, 3.0, 4.4}));
    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(ConditionSettings::direction), std::exception ); // Wrong Type

    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<PointType>(ConditionSettings::value), std::exception ); // Not set
    r_cond_settings_list.SetValue(ConditionSettings::value, PointType({2.1, 1.0, 2.4}));
    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(ConditionSettings::value), std::exception ); // Wrong Type

    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<double>(ConditionSettings::penalty_factor), std::exception ); // Not set
    r_cond_settings_list.SetValue(ConditionSettings::penalty_factor, 1e5);
    BOOST_REQUIRE_THROW( r_cond_settings_list.GetValue<IndexType>(ConditionSettings::penalty_factor), std::exception ); // Wrong Type
}


BOOST_AUTO_TEST_CASE(SettingsConditionSettingsDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Condition Settings Defaut Values" << std::endl;

    Settings settings;
    auto& r_cond_settings_list = settings.CreateNewConditionSettings();

    QuESo_CHECK( !r_cond_settings_list.IsSet(ConditionSettings::condition_id) );
    QuESo_CHECK( !r_cond_settings_list.IsSet(ConditionSettings::condition_type) );
    QuESo_CHECK( !r_cond_settings_list.IsSet(ConditionSettings::input_filename) );
    QuESo_CHECK( !r_cond_settings_list.IsSet(ConditionSettings::modulus) );
    QuESo_CHECK( !r_cond_settings_list.IsSet(ConditionSettings::direction) );
    QuESo_CHECK( !r_cond_settings_list.IsSet(ConditionSettings::value) );
    QuESo_CHECK( !r_cond_settings_list.IsSet(ConditionSettings::penalty_factor) );
};

BOOST_AUTO_TEST_CASE(SettingsConditionSettingsCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Condition Settings Customized Values" << std::endl;

    Settings settings;
    {
        auto& r_cond_settings_list = settings.CreateNewConditionSettings();

        r_cond_settings_list.SetValue(ConditionSettings::condition_id, 120UL);
        r_cond_settings_list.SetValue(ConditionSettings::condition_type, std::string("Hello"));
        r_cond_settings_list.SetValue(ConditionSettings::input_filename, std::string("hallo"));
        r_cond_settings_list.SetValue(ConditionSettings::modulus, 200.0);
        r_cond_settings_list.SetValue(ConditionSettings::penalty_factor, 300.0);
    }
    {
        auto& r_cond_settings_list = settings.CreateNewConditionSettings();

        r_cond_settings_list.SetValue(ConditionSettings::condition_id, 150UL);
        r_cond_settings_list.SetValue(ConditionSettings::condition_type, std::string("Hello2"));
        r_cond_settings_list.SetValue(ConditionSettings::input_filename, std::string("hallo2"));
        r_cond_settings_list.SetValue(ConditionSettings::direction, PointType{2.1, 3.2, 4.7});
        r_cond_settings_list.SetValue(ConditionSettings::value, PointType{2.2, 2.3, 1.4} );
    }
    {
        auto& r_cond_settings_list = settings[MainSettings::conditions_settings_list];
        QuESo_CHECK_EQUAL( r_cond_settings_list.NumberOfSubDictionaries(), 2 );

        QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::condition_id) );
        QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<IndexType>(ConditionSettings::condition_id), 120u );

        QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::condition_type) );
        QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<std::string>(ConditionSettings::condition_type), std::string("Hello") );

        QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::input_filename) );
        QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<std::string>(ConditionSettings::input_filename), std::string("hallo") );

        QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::modulus) );
        QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<double>(ConditionSettings::modulus), 200.0 );

        QuESo_CHECK( !r_cond_settings_list[0u].IsSet(ConditionSettings::direction) );
        BOOST_REQUIRE_THROW( r_cond_settings_list[0u].GetValue<PointType>(ConditionSettings::direction), std::exception );

        QuESo_CHECK( !r_cond_settings_list[0u].IsSet(ConditionSettings::value) );
        BOOST_REQUIRE_THROW( r_cond_settings_list[0u].GetValue<PointType>(ConditionSettings::value), std::exception );

        QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::penalty_factor) );
        QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<double>(ConditionSettings::penalty_factor), 300.0 );

    }
    {
        auto& r_cond_settings_list = settings[MainSettings::conditions_settings_list];
        QuESo_CHECK_EQUAL( r_cond_settings_list.NumberOfSubDictionaries(), 2 );

        QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::condition_id) );
        QuESo_CHECK_EQUAL( r_cond_settings_list[1u].GetValue<IndexType>(ConditionSettings::condition_id), 150u );

        QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::condition_type) );
        QuESo_CHECK_EQUAL( r_cond_settings_list[1u].GetValue<std::string>(ConditionSettings::condition_type), std::string("Hello2") );

        QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::input_filename) );
        QuESo_CHECK_EQUAL( r_cond_settings_list[1u].GetValue<std::string>(ConditionSettings::input_filename), std::string("hallo2") );

        QuESo_CHECK( !r_cond_settings_list[1u].IsSet(ConditionSettings::modulus) );
        BOOST_REQUIRE_THROW( r_cond_settings_list[1u].GetValue<double>(ConditionSettings::modulus), std::exception );

        QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::direction) );
        QuESo_CHECK_POINT_NEAR( r_cond_settings_list[1u].GetValue<PointType>(ConditionSettings::direction), PointType({2.1, 3.2, 4.7}), 1e-10 );

        QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::value) );
        QuESo_CHECK_POINT_NEAR( r_cond_settings_list[1u].GetValue<PointType>(ConditionSettings::value), PointType({2.2, 2.3, 1.4}), 1e-10 );

        QuESo_CHECK( !r_cond_settings_list[1u].IsSet(ConditionSettings::penalty_factor) );
        BOOST_REQUIRE_THROW( r_cond_settings_list[1u].GetValue<double>(ConditionSettings::penalty_factor), std::exception );
    }
};


BOOST_AUTO_TEST_CASE(SettingsCastAmbiguousTypesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Cast Ambiguous Types" << std::endl;

    Settings settings;

    // Cast IndexType to double
    BOOST_REQUIRE_THROW( settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::moment_fitting_residual, 5u), std::exception ); // Wrong Type
    // Lets set it with wrong type
    settings[MainSettings::trimmed_quadrature_rule_settings].SetValueWithAmbiguousType(TrimmedQuadratureRuleSettings::moment_fitting_residual, 5u);
    {
        const double r_double = settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual);
        QuESo_CHECK_RELATIVE_NEAR(r_double, 5.0, 1e-10);
    }
    // Lets set it with correct type
    settings[MainSettings::trimmed_quadrature_rule_settings].SetValueWithAmbiguousType(TrimmedQuadratureRuleSettings::moment_fitting_residual, 7.7);
    {
        const double r_double = settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual);
        QuESo_CHECK_RELATIVE_NEAR(r_double, 7.7, 1e-10);
    }

    // Cast Vector3i to Vector3d
    BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz), std::exception ); // Not set
    BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, Vector3i({0, 0, 0})), std::exception ); // Wrong Type
    // Lets set it with wrong type
    settings[MainSettings::background_grid_settings].SetValueWithAmbiguousType(BackgroundGridSettings::lower_bound_xyz, Vector3i({5, 3, 2}));
    {
        const PointType& r_point = settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz);
        QuESo_CHECK_POINT_RELATIVE_NEAR(r_point, PointType({5.0, 3.0, 2.0}), 1e-10);
    }
    // Lets set it with correct type
    settings[MainSettings::background_grid_settings].SetValueWithAmbiguousType(BackgroundGridSettings::lower_bound_xyz, Vector3d({5.5, 2.2, 1.1}));
    {
        const PointType& r_point = settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz);
        QuESo_CHECK_POINT_RELATIVE_NEAR(r_point, PointType({5.5, 2.2, 1.1}), 1e-10);
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

