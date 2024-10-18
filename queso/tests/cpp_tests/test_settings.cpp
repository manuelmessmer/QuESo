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
#include "queso/includes/settings.hpp"

namespace queso {


namespace Testing {

BOOST_AUTO_TEST_SUITE( SettingsTestSuite )

BOOST_AUTO_TEST_CASE(SettingsWrongTypeTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Wrong Types" << std::endl;

    {   /// Enum access
        Settings setting;
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW(setting[GeneralSettings::echo_level], std::exception); // Wrong Key type

            /// General settings
            auto& r_general_settings = setting[MainSettings::general_settings];
            BOOST_REQUIRE_THROW( r_general_settings.IsSet(MainSettings::general_settings), std::exception); // Wrong Key type

            BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception); // Wrong Key type

            BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::input_filename), std::exception); // Wrong Value type
            BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::output_directory_name), std::exception); // Wrong Value type
            BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::echo_level), std::exception); // Wrong Value type
            BOOST_REQUIRE_THROW(r_general_settings.SetValue(GeneralSettings::echo_level, -1), std::exception); // Negative value
            BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::write_output_to_file), std::exception); // Wrong Value type

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
            r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 2.0, 3.0});
            BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, IndexType(2)), std::exception); // Wrong Value type
            BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>(BackgroundGridSettings::upper_bound_uvw), std::exception); // Wrong Value type
            const PointType& r_upper_bound_uvw = r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw);
            QuESo_CHECK_POINT_NEAR(r_upper_bound_uvw, PointType({1.0, 2.0, 3.0}), 1e-10 );

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
    }
    {   /// String access
        Settings setting;

        BOOST_REQUIRE_THROW(setting["echo_level"], std::exception); // Wrong Key type

        /// General settings
        auto& r_general_settings = setting["general_settings"];
        BOOST_REQUIRE_THROW( r_general_settings.IsSet("general_settings"), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>("lower_bound_uvw"), std::exception); // Wrong Key type

        BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>("input_filename"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>("output_directory_name"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>("echo_level"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_general_settings.SetValue("echo_level", -1), std::exception); // Negative value
        BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>("write_output_to_file"), std::exception); // Wrong Value type

        /// Mesh settings
        auto& r_mesh_settings = setting["background_grid_settings"];
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<IndexType>("echo_level"), std::exception); // Wrong Key type

        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>("lower_bound_xyz"), std::exception); // Not set
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("echo_level", 2u), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("lower_bound_xyz", IndexType(2)), std::exception); // Wrong Value type
        r_mesh_settings.SetValue("lower_bound_xyz", PointType{1.0, 2.0, 3.0});
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("echo_level", IndexType(2)), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("lower_bound_xyz", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>("lower_bound_xyz"), std::exception); // Wrong Value type
        const PointType& r_lower_bound_xyz = r_mesh_settings.GetValue<PointType>("lower_bound_xyz");
        QuESo_CHECK_POINT_NEAR(r_lower_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>("upper_bound_xyz"), std::exception); // Not set
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("upper_bound_xyz", IndexType(2)), std::exception); // Wrong Value type
        r_mesh_settings.SetValue("upper_bound_xyz", PointType{1.0, 2.0, 3.0});
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("upper_bound_xyz", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>("upper_bound_xyz"), std::exception); // Wrong Value type
        const PointType& r_upper_bound_xyz = r_mesh_settings.GetValue<PointType>("upper_bound_xyz");
        QuESo_CHECK_POINT_NEAR(r_upper_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>("lower_bound_uvw"), std::exception); // Not set
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("lower_bound_uvw", IndexType(2)), std::exception); // Wrong Value type
        r_mesh_settings.SetValue("lower_bound_uvw", PointType{1.0, 2.0, 3.0});
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("lower_bound_uvw", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>("lower_bound_uvw"), std::exception); // Wrong Value type
        const PointType& r_lower_bound_uvw= r_mesh_settings.GetValue<PointType>("lower_bound_uvw");
        QuESo_CHECK_POINT_NEAR(r_lower_bound_uvw, PointType({1.0, 2.0, 3.0}), 1e-10);

        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>("upper_bound_uvw"), std::exception); // Not set
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("upper_bound_uvw", IndexType(2)), std::exception); // Wrong Value type
        r_mesh_settings.SetValue("upper_bound_uvw", PointType{1.0, 2.0, 3.0});
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("upper_bound_uvw", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>("upper_bound_uvw"), std::exception); // Wrong Value type
        const PointType& r_upper_bound_uvw = r_mesh_settings.GetValue<PointType>("upper_bound_uvw");
        QuESo_CHECK_POINT_NEAR(r_upper_bound_uvw, PointType({1.0, 2.0, 3.0}), 1e-10 );

        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>("polynomial_order"), std::exception); // Not set
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("polynomial_order", IndexType(2)), std::exception); // Wrong Value type
        r_mesh_settings.SetValue("polynomial_order", Vector3i{1, 2, 3});
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("polynomial_order", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>("polynomial_order"), std::exception); // Wrong Value type
        const Vector3i& r_polynomial_order = r_mesh_settings.GetValue<Vector3i>("polynomial_order");
        QuESo_CHECK_Vector3i_EQUAL(r_polynomial_order, Vector3i({1, 2, 3}) );

        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>("number_of_elements"), std::exception); // Not set
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("number_of_elements", IndexType(2)), std::exception); // Wrong Value type
        r_mesh_settings.SetValue("number_of_elements", Vector3i{1, 2, 3});
        BOOST_REQUIRE_THROW(r_mesh_settings.SetValue("number_of_elements", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>("number_of_elements"), std::exception); // Wrong Value type
        const Vector3i& number_of_elements = r_mesh_settings.GetValue<Vector3i>("number_of_elements");
        QuESo_CHECK_Vector3i_EQUAL(number_of_elements, Vector3i({1, 2, 3}) );

        /// Trimmed quadrature rule settings
        auto& r_trimmed_quad_rule_settings = setting["trimmed_quadrature_rule_settings"];
        BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.SetValue("number_of_elements", Vector3i{1, 2, 3}), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>("moment_fitting_residual"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>("min_element_volume_ratio"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>("min_num_boundary_triangles"), std::exception); // Wrong Value type

        /// Non trimmed quadrature rule settings
        auto& r_non_trimmed_quad_rule_settings = setting["non_trimmed_quadrature_rule_settings"];
        BOOST_REQUIRE_THROW(r_non_trimmed_quad_rule_settings.SetValue("number_of_elements", Vector3i{1, 2, 3}), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(r_non_trimmed_quad_rule_settings.GetValue<IndexType>("integration_method"), std::exception); // Wrong Value type
    }
}


BOOST_AUTO_TEST_CASE(SettingsDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Default Values" << std::endl;

    {   /// Enum access
        Settings settings;

        /// General settings
        QuESo_CHECK( !settings[MainSettings::general_settings].IsSet(GeneralSettings::input_filename) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::input_filename), std::exception );
        }
        QuESo_CHECK( settings[MainSettings::general_settings].IsSet(GeneralSettings::output_directory_name) );
        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::output_directory_name), std::string("queso_output") );

        QuESo_CHECK( settings[MainSettings::general_settings].IsSet(GeneralSettings::echo_level) );
        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level), 1u);

        QuESo_CHECK( settings[MainSettings::general_settings].IsSet(GeneralSettings::write_output_to_file) );
        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<bool>(GeneralSettings::write_output_to_file), true);

        /// Mesh settings
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::grid_type) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::grid_type), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::lower_bound_xyz) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::upper_bound_xyz) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::lower_bound_uvw) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::upper_bound_uvw) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::polynomial_order) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::polynomial_order), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::number_of_elements) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::number_of_elements), std::exception );
        }
        /// TrimmedQuadratureRuleSettings settings
        QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::moment_fitting_residual) );
        QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual), 1e-10,1e-10 );

        QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::min_element_volume_ratio) );
        QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 1e-3,1e-10 );

        QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::min_num_boundary_triangles) );
        QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 100 );

        QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed) );
        QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<bool>(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed), true );

        // NonTrimmedQuadratureRuleSettings settings
        QuESo_CHECK( settings[MainSettings::non_trimmed_quadrature_rule_settings].IsSet(NonTrimmedQuadratureRuleSettings::integration_method) );
        QuESo_CHECK_EQUAL( settings[MainSettings::non_trimmed_quadrature_rule_settings].GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::gauss );
    }
    {   /// String access
        Settings settings;

        /// General settings
        QuESo_CHECK( !settings["general_settings"].IsSet("input_filename") );
        BOOST_REQUIRE_THROW( settings["general_settings"].GetValue<std::string>("input_filename"), std::exception );

        QuESo_CHECK( settings["general_settings"].IsSet("output_directory_name") );
        QuESo_CHECK_EQUAL( settings["general_settings"].GetValue<std::string>("output_directory_name"), std::string("queso_output") );

        QuESo_CHECK( settings["general_settings"].IsSet("echo_level") );
        QuESo_CHECK_EQUAL( settings["general_settings"].GetValue<IndexType>("echo_level"), 1u);

        QuESo_CHECK( settings["general_settings"].IsSet("write_output_to_file") );
        QuESo_CHECK_EQUAL( settings["general_settings"].GetValue<bool>("write_output_to_file"), true);

        /// Mesh settings
        QuESo_CHECK( !settings["background_grid_settings"].IsSet("grid_type") );
        BOOST_REQUIRE_THROW( settings["background_grid_settings"].GetValue<PointType>("grid_type"), std::exception );

        QuESo_CHECK( !settings["background_grid_settings"].IsSet("lower_bound_xyz") );
        BOOST_REQUIRE_THROW( settings["background_grid_settings"].GetValue<PointType>("lower_bound_xyz"), std::exception );

        QuESo_CHECK( !settings["background_grid_settings"].IsSet("upper_bound_xyz") );
        BOOST_REQUIRE_THROW( settings["background_grid_settings"].GetValue<PointType>("upper_bound_xyz"), std::exception );

        QuESo_CHECK( !settings["background_grid_settings"].IsSet("lower_bound_uvw") );
        BOOST_REQUIRE_THROW( settings["background_grid_settings"].GetValue<PointType>("lower_bound_uvw"), std::exception );

        QuESo_CHECK( !settings["background_grid_settings"].IsSet(BackgroundGridSettings::upper_bound_uvw) );
        BOOST_REQUIRE_THROW( settings["background_grid_settings"].GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw), std::exception );

        QuESo_CHECK( !settings["background_grid_settings"].IsSet("polynomial_order") );
        BOOST_REQUIRE_THROW( settings["background_grid_settings"].GetValue<Vector3i>("polynomial_order"), std::exception );

        QuESo_CHECK( !settings["background_grid_settings"].IsSet("number_of_elements") );
        BOOST_REQUIRE_THROW( settings["background_grid_settings"].GetValue<Vector3i>("number_of_elements"), std::exception );

        /// TrimmedQuadratureRuleSettings settings
        QuESo_CHECK( settings["trimmed_quadrature_rule_settings"].IsSet("moment_fitting_residual") );
        QuESo_CHECK_RELATIVE_NEAR( settings["trimmed_quadrature_rule_settings"].GetValue<double>("moment_fitting_residual"), 1e-10,1e-10 );

        QuESo_CHECK( settings["trimmed_quadrature_rule_settings"].IsSet("min_element_volume_ratio") );
        QuESo_CHECK_RELATIVE_NEAR( settings["trimmed_quadrature_rule_settings"].GetValue<double>("min_element_volume_ratio"), 1e-3,1e-10 );

        QuESo_CHECK( settings["trimmed_quadrature_rule_settings"].IsSet("min_num_boundary_triangles") );
        QuESo_CHECK_EQUAL( settings["trimmed_quadrature_rule_settings"].GetValue<IndexType>("min_num_boundary_triangles"), 100 );

        QuESo_CHECK( settings["trimmed_quadrature_rule_settings"].IsSet("neglect_elements_if_stl_is_flawed") );
        QuESo_CHECK_EQUAL( settings["trimmed_quadrature_rule_settings"].GetValue<bool>("neglect_elements_if_stl_is_flawed"), true );

        // NonTrimmedQuadratureRuleSettings settings
        QuESo_CHECK( settings["non_trimmed_quadrature_rule_settings"].IsSet("integration_method") );
        QuESo_CHECK_EQUAL( settings["non_trimmed_quadrature_rule_settings"].GetValue<IntegrationMethod>("integration_method"), IntegrationMethod::gauss );
    }
}


BOOST_AUTO_TEST_CASE(SettingsCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Customized Values" << std::endl;

    {   /// Enum access
        Settings settings;

        /// General settings
        settings[MainSettings::general_settings].SetValue(GeneralSettings::input_filename, std::string("test_filename.stl"));
        settings[MainSettings::general_settings].SetValue(GeneralSettings::output_directory_name, std::string("new_output/"));
        settings[MainSettings::general_settings].SetValue(GeneralSettings::echo_level, 2u);
        settings[MainSettings::general_settings].SetValue(GeneralSettings::write_output_to_file, false);


        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::input_filename), std::string("test_filename.stl") );

        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::output_directory_name), std::string("new_output/") );

        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level), 2u );

        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<bool>(GeneralSettings::write_output_to_file), false );

        /// Mesh settings
        settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
        settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_xyz, PointType({1.0, 1.0, 2.0}));
        settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_xyz, PointType({0.0, 3.0, 2.0}));
        settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::lower_bound_uvw, PointType({0.1, -1.0, 2.0}));
        settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::upper_bound_uvw, PointType({0.66, 1.0, 2.2}));
        settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::polynomial_order, Vector3i({5, 6, 7}));
        settings[MainSettings::background_grid_settings].SetValue(BackgroundGridSettings::number_of_elements, Vector3i({8, 9, 2}));

        QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<GridType>(BackgroundGridSettings::grid_type), GridType::b_spline_grid );
        QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz), PointType({1.0, 1.0, 2.0}) );
        QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz), PointType({0.0, 3.0, 2.0}) );
        QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), PointType({0.1, -1.0, 2.0}) );
        QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw), PointType({0.66, 1.0, 2.2}) );
        QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::polynomial_order), Vector3i({5, 6, 7}) );
        QuESo_CHECK_EQUAL( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::number_of_elements), Vector3i({8, 9, 2}) );


        /// Trimmed quadrature rule settings
        settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::moment_fitting_residual, 5.6e-5);
        settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 0.45);
        settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, IndexType(234));
        settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed, false);

        QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual), 5.6e-5, 1e-10 );
        QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 0.45, 1e-10 );
        QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 234u );
        QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<bool>(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed), false );

        /// Non trimmed quadrature rule settings
        settings[MainSettings::non_trimmed_quadrature_rule_settings].SetValue(NonTrimmedQuadratureRuleSettings::integration_method, IntegrationMethod::ggq_optimal);

        QuESo_CHECK_EQUAL( settings[MainSettings::non_trimmed_quadrature_rule_settings].GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::ggq_optimal );
    }
    {   /// String acccess
        Settings settings;

        /// General settings
        settings["general_settings"].SetValue("input_filename", std::string("test_filename.stl"));
        settings["general_settings"].SetValue("output_directory_name", std::string("new_output/"));
        settings["general_settings"].SetValue("echo_level", 2u);
        settings["general_settings"].SetValue("write_output_to_file", false);


        QuESo_CHECK_EQUAL( settings["general_settings"].GetValue<std::string>("input_filename"), std::string("test_filename.stl") );

        QuESo_CHECK_EQUAL( settings["general_settings"].GetValue<std::string>("output_directory_name"), std::string("new_output/") );
        QuESo_CHECK_EQUAL( settings["general_settings"].GetValue<IndexType>("echo_level"), 2u );
        QuESo_CHECK_EQUAL( settings["general_settings"].GetValue<bool>("write_output_to_file"), false );

        /// Mesh settings
        settings["background_grid_settings"].SetValue("grid_type", GridType::b_spline_grid);
        settings["background_grid_settings"].SetValue("lower_bound_xyz", PointType({1.0, 1.0, 2.0}));
        settings["background_grid_settings"].SetValue("upper_bound_xyz", PointType({0.0, 3.0, 2.0}));
        settings["background_grid_settings"].SetValue("lower_bound_uvw", PointType({0.1, -1.0, 2.0}));
        settings["background_grid_settings"].SetValue("upper_bound_uvw", PointType({0.66, 1.0, 2.2}));
        settings["background_grid_settings"].SetValue("polynomial_order", Vector3i({5, 6, 7}));
        settings["background_grid_settings"].SetValue("number_of_elements", Vector3i({8, 9, 2}));

        QuESo_CHECK_EQUAL( settings["background_grid_settings"].GetValue<GridType>("grid_type"), GridType::b_spline_grid );
        QuESo_CHECK_EQUAL( settings["background_grid_settings"].GetValue<PointType>("lower_bound_xyz"), PointType({1.0, 1.0, 2.0}) );
        QuESo_CHECK_EQUAL( settings["background_grid_settings"].GetValue<PointType>("upper_bound_xyz"), PointType({0.0, 3.0, 2.0}) );
        QuESo_CHECK_EQUAL( settings["background_grid_settings"].GetValue<PointType>("lower_bound_uvw"), PointType({0.1, -1.0, 2.0}) );
        QuESo_CHECK_EQUAL( settings["background_grid_settings"].GetValue<PointType>("upper_bound_uvw"), PointType({0.66, 1.0, 2.2}) );
        QuESo_CHECK_EQUAL( settings["background_grid_settings"].GetValue<Vector3i>("polynomial_order"), Vector3i({5, 6, 7}) );
        QuESo_CHECK_EQUAL( settings["background_grid_settings"].GetValue<Vector3i>("number_of_elements"), Vector3i({8, 9, 2}) );


        /// Trimmed quadrature rule settings
        settings["trimmed_quadrature_rule_settings"].SetValue("moment_fitting_residual", 5.6e-5);
        settings["trimmed_quadrature_rule_settings"].SetValue("min_element_volume_ratio", 0.45);
        settings["trimmed_quadrature_rule_settings"].SetValue("min_num_boundary_triangles", IndexType(234));
        settings["trimmed_quadrature_rule_settings"].SetValue("neglect_elements_if_stl_is_flawed", false);

        QuESo_CHECK_RELATIVE_NEAR( settings["trimmed_quadrature_rule_settings"].GetValue<double>("moment_fitting_residual"), 5.6e-5, 1e-10 );
        QuESo_CHECK_RELATIVE_NEAR( settings["trimmed_quadrature_rule_settings"].GetValue<double>("min_element_volume_ratio"), 0.45, 1e-10 );
        QuESo_CHECK_EQUAL( settings["trimmed_quadrature_rule_settings"].GetValue<IndexType>("min_num_boundary_triangles"), 234u );
        QuESo_CHECK_EQUAL( settings["trimmed_quadrature_rule_settings"].GetValue<bool>("neglect_elements_if_stl_is_flawed"), false );

        /// Non trimmed quadrature rule settings
        settings["non_trimmed_quadrature_rule_settings"].SetValue("integration_method", IntegrationMethod::ggq_optimal);
        QuESo_CHECK_EQUAL( settings["non_trimmed_quadrature_rule_settings"].GetValue<IntegrationMethod>("integration_method"), IntegrationMethod::ggq_optimal );
    }
}

BOOST_AUTO_TEST_CASE(SettingsConditionSettingsWrongTypeTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Condition Settings Wrong Types " << std::endl;

    {   /// Enum access
        Settings settings;

        if( !NOTDEBUG ) {
            auto& r_cond_settings = settings.CreateNewConditionSettings();
            BOOST_REQUIRE_THROW( r_cond_settings.IsSet(BackgroundGridSettings::lower_bound_uvw), std::exception ); // Wrong Key

            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(ConditionSettings::condition_id), std::exception ); // Not set
            BOOST_REQUIRE_THROW( r_cond_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType({2.0, 3.0, 4.0})), std::exception ); // Wrong Key type
            BOOST_REQUIRE_THROW( r_cond_settings.SetValue(ConditionSettings::condition_id, 2.0), std::exception ); // Wrong Key type

            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(ConditionSettings::condition_id), std::exception ); // Not set
            r_cond_settings.SetValue(ConditionSettings::condition_id, 100u);
            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(GeneralSettings::echo_level), std::exception ); // Wrong Key type
            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<double>(ConditionSettings::condition_id), std::exception ); // Wrong Type

            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<std::string>(ConditionSettings::condition_type), std::exception ); // Not set
            r_cond_settings.SetValue(ConditionSettings::condition_type, std::string("dummy"));
            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(GeneralSettings::echo_level), std::exception ); // Wrong Key type
            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(ConditionSettings::condition_type), std::exception ); // Wrong Type

            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<std::string>(ConditionSettings::input_filename), std::exception ); // Not set
            r_cond_settings.SetValue(ConditionSettings::input_filename, std::string("dummy_2"));
            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(ConditionSettings::input_filename), std::exception ); // Wrong Type

            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<double>(ConditionSettings::modulus), std::exception ); // Not set
            r_cond_settings.SetValue(ConditionSettings::modulus, 2.0);
            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(ConditionSettings::modulus), std::exception ); // Wrong Type

            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<PointType>(ConditionSettings::direction), std::exception ); // Not set
            r_cond_settings.SetValue(ConditionSettings::direction, PointType({2.2, 3.0, 4.4}));
            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(ConditionSettings::direction), std::exception ); // Wrong Type

            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<PointType>(ConditionSettings::value), std::exception ); // Not set
            r_cond_settings.SetValue(ConditionSettings::value, PointType({2.1, 1.0, 2.4}));
            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(ConditionSettings::value), std::exception ); // Wrong Type

            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<double>(ConditionSettings::penalty_factor), std::exception ); // Not set
            r_cond_settings.SetValue(ConditionSettings::penalty_factor, 1e5);
            BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>(ConditionSettings::penalty_factor), std::exception ); // Wrong Type
        }
    }
    {   /// String access
        Settings settings;

        auto& r_cond_settings = settings.CreateNewConditionSettings();
        BOOST_REQUIRE_THROW( r_cond_settings.IsSet("lower_bound_uvw"), std::exception ); // Wrong Key

        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("condition_id"), std::exception ); // Not set
        BOOST_REQUIRE_THROW( r_cond_settings.SetValue("lower_bound_uvw", PointType({2.0, 3.0, 4.0})), std::exception ); // Wrong Key type
        BOOST_REQUIRE_THROW( r_cond_settings.SetValue("condition_id", 2.0), std::exception ); // Wrong Key type

        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("condition_id"), std::exception ); // Not set
        r_cond_settings.SetValue("condition_id", 100u);
        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("echo_level"), std::exception ); // Wrong Key type
        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<double>("condition_id"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<std::string>("condition_type"), std::exception ); // Not set
        r_cond_settings.SetValue("condition_type", std::string("dummy"));
        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("echo_level"), std::exception ); // Wrong Key type
        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("condition_type"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<std::string>("input_filename"), std::exception ); // Not set
        r_cond_settings.SetValue("input_filename", std::string("dummy_2"));
        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("input_filename"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<double>("modulus"), std::exception ); // Not set
        r_cond_settings.SetValue("modulus", 2.0);
        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("modulus"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<PointType>("direction"), std::exception ); // Not set
        r_cond_settings.SetValue("direction", PointType({2.2, 3.0, 4.4}));
        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("direction"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<PointType>("value"), std::exception ); // Not set
        r_cond_settings.SetValue("value", PointType({2.1, 1.0, 2.4}));
        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("value"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<double>("penalty_factor"), std::exception ); // Not set
        r_cond_settings.SetValue("penalty_factor", 1e5);
        BOOST_REQUIRE_THROW( r_cond_settings.GetValue<IndexType>("penalty_factor"), std::exception ); // Wrong Type
    }
}


BOOST_AUTO_TEST_CASE(SettingsConditionSettingsDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Condition Settings Defaut Values" << std::endl;

    {   /// Enum access
        Settings settings;
        auto& r_cond_settings = settings.CreateNewConditionSettings();

        QuESo_CHECK( !r_cond_settings.IsSet(ConditionSettings::condition_id) );
        QuESo_CHECK( !r_cond_settings.IsSet(ConditionSettings::condition_type) );
        QuESo_CHECK( !r_cond_settings.IsSet(ConditionSettings::input_filename) );
        QuESo_CHECK( !r_cond_settings.IsSet(ConditionSettings::modulus) );
        QuESo_CHECK( !r_cond_settings.IsSet(ConditionSettings::direction) );
        QuESo_CHECK( !r_cond_settings.IsSet(ConditionSettings::value) );
        QuESo_CHECK( !r_cond_settings.IsSet(ConditionSettings::penalty_factor) );
    }
    {   /// String access
        Settings settings;
        auto& r_cond_settings = settings.CreateNewConditionSettings();

        QuESo_CHECK( !r_cond_settings.IsSet("condition_id") );
        QuESo_CHECK( !r_cond_settings.IsSet("condition_type") );
        QuESo_CHECK( !r_cond_settings.IsSet("input_filename") );
        QuESo_CHECK( !r_cond_settings.IsSet("modulus") );
        QuESo_CHECK( !r_cond_settings.IsSet("direction") );
        QuESo_CHECK( !r_cond_settings.IsSet("value") );
        QuESo_CHECK( !r_cond_settings.IsSet("penalty_factor") );
    }
};

BOOST_AUTO_TEST_CASE(SettingsConditionSettingsCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Condition Settings Customized Values" << std::endl;

    {   /// Enum acess
        Settings settings;
        {
            auto& r_cond_settings = settings.CreateNewConditionSettings();

            r_cond_settings.SetValue(ConditionSettings::condition_id, 120UL);
            r_cond_settings.SetValue(ConditionSettings::condition_type, std::string("Hello"));
            r_cond_settings.SetValue(ConditionSettings::input_filename, std::string("hallo"));
            r_cond_settings.SetValue(ConditionSettings::modulus, 200.0);
            r_cond_settings.SetValue(ConditionSettings::penalty_factor, 300.0);
        }
        {
            auto& r_cond_settings = settings.CreateNewConditionSettings();

            r_cond_settings.SetValue(ConditionSettings::condition_id, 150UL);
            r_cond_settings.SetValue(ConditionSettings::condition_type, std::string("Hello2"));
            r_cond_settings.SetValue(ConditionSettings::input_filename, std::string("hallo2"));
            r_cond_settings.SetValue(ConditionSettings::direction, PointType{2.1, 3.2, 4.7});
            r_cond_settings.SetValue(ConditionSettings::value, PointType{2.2, 2.3, 1.4} );
        }
        {
            auto& r_cond_settings_list = settings.GetList(MainSettings::conditions_settings_list);
            QuESo_CHECK_EQUAL( r_cond_settings_list.size(), 2 );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::condition_id) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<IndexType>(ConditionSettings::condition_id), 120u );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::condition_type) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<std::string>(ConditionSettings::condition_type), std::string("Hello") );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::input_filename) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<std::string>(ConditionSettings::input_filename), std::string("hallo") );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::modulus) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<double>(ConditionSettings::modulus), 200.0 );

            QuESo_CHECK( !r_cond_settings_list[0u].IsSet(ConditionSettings::direction) );
            if( !NOTDEBUG )
                BOOST_REQUIRE_THROW( r_cond_settings_list[0u].GetValue<PointType>(ConditionSettings::direction), std::exception );

            QuESo_CHECK( !r_cond_settings_list[0u].IsSet(ConditionSettings::value) );
            if( !NOTDEBUG )
                BOOST_REQUIRE_THROW( r_cond_settings_list[0u].GetValue<PointType>(ConditionSettings::value), std::exception );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet(ConditionSettings::penalty_factor) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<double>(ConditionSettings::penalty_factor), 300.0 );
        }
        {
            auto& r_cond_settings_list = settings.GetList(MainSettings::conditions_settings_list);
            QuESo_CHECK_EQUAL( r_cond_settings_list.size(), 2 );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::condition_id) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[1u].GetValue<IndexType>(ConditionSettings::condition_id), 150u );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::condition_type) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[1u].GetValue<std::string>(ConditionSettings::condition_type), std::string("Hello2") );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::input_filename) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[1u].GetValue<std::string>(ConditionSettings::input_filename), std::string("hallo2") );

            QuESo_CHECK( !r_cond_settings_list[1u].IsSet(ConditionSettings::modulus) );
            if( !NOTDEBUG )
                BOOST_REQUIRE_THROW( r_cond_settings_list[1u].GetValue<double>(ConditionSettings::modulus), std::exception );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::direction) );
            QuESo_CHECK_POINT_NEAR( r_cond_settings_list[1u].GetValue<PointType>(ConditionSettings::direction), PointType({2.1, 3.2, 4.7}), 1e-10 );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet(ConditionSettings::value) );
            QuESo_CHECK_POINT_NEAR( r_cond_settings_list[1u].GetValue<PointType>(ConditionSettings::value), PointType({2.2, 2.3, 1.4}), 1e-10 );

            QuESo_CHECK( !r_cond_settings_list[1u].IsSet(ConditionSettings::penalty_factor) );
            if( !NOTDEBUG )
                BOOST_REQUIRE_THROW( r_cond_settings_list[1u].GetValue<double>(ConditionSettings::penalty_factor), std::exception );
        }
    }
    {   /// String access
        Settings settings;
        {
            auto& r_cond_settings = settings.CreateNewConditionSettings();

            r_cond_settings.SetValue("condition_id", 120UL);
            r_cond_settings.SetValue("condition_type", std::string("Hello"));
            r_cond_settings.SetValue("input_filename", std::string("hallo"));
            r_cond_settings.SetValue("modulus", 200.0);
            r_cond_settings.SetValue("penalty_factor", 300.0);
        }
        {
            auto& r_cond_settings = settings.CreateNewConditionSettings();

            r_cond_settings.SetValue("condition_id", 150UL);
            r_cond_settings.SetValue("condition_type", std::string("Hello2"));
            r_cond_settings.SetValue("input_filename", std::string("hallo2"));
            r_cond_settings.SetValue("direction", PointType{2.1, 3.2, 4.7});
            r_cond_settings.SetValue("value", PointType{2.2, 2.3, 1.4} );
        }
        {
            auto& r_cond_settings_list = settings.GetList("conditions_settings_list");
            QuESo_CHECK_EQUAL( r_cond_settings_list.size(), 2 );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet("condition_id") );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<IndexType>("condition_id"), 120u );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet("condition_type") );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<std::string>("condition_type"), std::string("Hello") );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet("input_filename") );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<std::string>("input_filename"), std::string("hallo") );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet("modulus") );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<double>("modulus"), 200.0 );

            QuESo_CHECK( !r_cond_settings_list[0u].IsSet("direction") );
            BOOST_REQUIRE_THROW( r_cond_settings_list[0u].GetValue<PointType>("direction"), std::exception );

            QuESo_CHECK( !r_cond_settings_list[0u].IsSet("value") );
            BOOST_REQUIRE_THROW( r_cond_settings_list[0u].GetValue<PointType>("value"), std::exception );

            QuESo_CHECK( r_cond_settings_list[0u].IsSet("penalty_factor") );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u].GetValue<double>("penalty_factor"), 300.0 );
        }
        {
            auto& r_cond_settings_list = settings.GetList("conditions_settings_list");
            QuESo_CHECK_EQUAL( r_cond_settings_list.size(), 2 );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet("condition_id") );
            QuESo_CHECK_EQUAL( r_cond_settings_list[1u].GetValue<IndexType>("condition_id"), 150u );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet("condition_type") );
            QuESo_CHECK_EQUAL( r_cond_settings_list[1u].GetValue<std::string>("condition_type"), std::string("Hello2") );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet("input_filename") );
            QuESo_CHECK_EQUAL( r_cond_settings_list[1u].GetValue<std::string>("input_filename"), std::string("hallo2") );

            QuESo_CHECK( !r_cond_settings_list[1u].IsSet("modulus") );
            BOOST_REQUIRE_THROW( r_cond_settings_list[1u].GetValue<double>("modulus"), std::exception );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet("direction") );
            QuESo_CHECK_POINT_NEAR( r_cond_settings_list[1u].GetValue<PointType>("direction"), PointType({2.1, 3.2, 4.7}), 1e-10 );

            QuESo_CHECK( r_cond_settings_list[1u].IsSet("value") );
            QuESo_CHECK_POINT_NEAR( r_cond_settings_list[1u].GetValue<PointType>("value"), PointType({2.2, 2.3, 1.4}), 1e-10 );

            QuESo_CHECK( !r_cond_settings_list[1u].IsSet("penalty_factor") );
            BOOST_REQUIRE_THROW( r_cond_settings_list[1u].GetValue<double>("penalty_factor"), std::exception );
        }
    }
};


BOOST_AUTO_TEST_CASE(SettingsCastAmbiguousTypesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Cast Ambiguous Types" << std::endl;

    Settings settings;

    /// Cast IndexType to double
    BOOST_REQUIRE_THROW( settings["trimmed_quadrature_rule_settings"].SetValue("moment_fitting_residual", 5u), std::exception ); // Wrong Type
    BOOST_REQUIRE_THROW( settings["trimmed_quadrature_rule_settings"].SetValueWithAmbiguousType("moment_fitting_residual", -1), std::exception ); // Negative value
    // Lets set it with wrong type
    settings["trimmed_quadrature_rule_settings"].SetValueWithAmbiguousType("moment_fitting_residual", 5u);
    {
        const double r_double = settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual);
        QuESo_CHECK_RELATIVE_NEAR(r_double, 5.0, 1e-10);
    }
    // Lets set it with correct type
    settings["trimmed_quadrature_rule_settings"].SetValueWithAmbiguousType("moment_fitting_residual", 7.7);
    {
        const double r_double = settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual);
        QuESo_CHECK_RELATIVE_NEAR(r_double, 7.7, 1e-10);
    }

    // Cast Vector3i to Vector3d
    BOOST_REQUIRE_THROW( settings["background_grid_settings"].GetValue<PointType>("lower_bound_xyz"), std::exception ); // Not set
    BOOST_REQUIRE_THROW( settings["background_grid_settings"].SetValue("lower_bound_xyz", Vector3i({0, 0, 0})), std::exception ); // Wrong Type
    // Lets set it with wrong type
    settings["background_grid_settings"].SetValueWithAmbiguousType("lower_bound_xyz", Vector3i({5, 3, 2}));
    {
        const PointType& r_point = settings["background_grid_settings"].GetValue<PointType>("lower_bound_xyz");
        QuESo_CHECK_POINT_RELATIVE_NEAR(r_point, PointType({5.0, 3.0, 2.0}), 1e-10);
    }
    // Lets set it with correct type
    settings["background_grid_settings"].SetValueWithAmbiguousType("lower_bound_xyz", Vector3d({5.5, 2.2, 1.1}));
    {
        const PointType& r_point = settings["background_grid_settings"].GetValue<PointType>("lower_bound_xyz");
        QuESo_CHECK_POINT_RELATIVE_NEAR(r_point, PointType({5.5, 2.2, 1.1}), 1e-10);
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

