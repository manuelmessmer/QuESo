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

//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/includes/dictionary_factory.hpp"

namespace queso {

namespace Testing {

BOOST_AUTO_TEST_SUITE( SettingsTestSuite )

BOOST_AUTO_TEST_CASE(SettingsWrongTypeTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Wrong Types" << std::endl;

    {   /// Key access
        auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& settings = *p_settings;

        auto& r_mesh_settings = settings[MainSettings::background_grid_settings];
        // GetValue must throw also in release mode, when value is not set.
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz), std::exception); // Not set
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz), std::exception); // Not set
        BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception); // Not set

        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW(settings[MainInfo::background_grid_info], std::exception); // Wrong Key type

            // General settings
            auto& r_general_settings = settings[MainSettings::general_settings];
            BOOST_REQUIRE_THROW(r_general_settings.IsSet(ConditionSettings::condition_id), std::exception); // Wrong Key type
            BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception); // Wrong Key type
            BOOST_REQUIRE_THROW(r_general_settings.GetValueFast<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception); // Wrong Key type

            /// Mesh settings
            BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<IndexType>(GeneralSettings::echo_level), std::exception); // Wrong Key type

            BOOST_REQUIRE_THROW(r_mesh_settings.GetValueFast<PointType>(BackgroundGridSettings::lower_bound_xyz), std::exception); // Not set
            BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(GeneralSettings::echo_level, 2u), std::exception); // Wrong Key type

            r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{1.0, 2.0, 3.0});
            BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(GeneralSettings::echo_level, IndexType(2)), std::exception); // Wrong Key type
            const PointType& r_lower_bound_xyz = r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz);
            QuESo_CHECK_POINT_NEAR(r_lower_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

            BOOST_REQUIRE_THROW(r_mesh_settings.GetValueFast<PointType>(BackgroundGridSettings::upper_bound_xyz), std::exception); // Not set
            r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{1.0, 2.0, 3.0});
            const PointType& r_upper_bound_xyz = r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz);
            QuESo_CHECK_POINT_NEAR(r_upper_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

            BOOST_REQUIRE_THROW(r_mesh_settings.GetValueFast<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception); // Not set
            r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{1.0, 2.0, 3.0});
            const PointType& r_lower_bound_uvw= r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw);
            QuESo_CHECK_POINT_NEAR(r_lower_bound_uvw, PointType({1.0, 2.0, 3.0}), 1e-10);

            BOOST_REQUIRE_THROW(r_mesh_settings.GetValueFast<PointType>(BackgroundGridSettings::upper_bound_uvw), std::exception); // Not set
            r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{1.0, 2.0, 3.0});
            const PointType& r_upper_bound_uvw = r_mesh_settings.GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw);
            QuESo_CHECK_POINT_NEAR(r_upper_bound_uvw, PointType({1.0, 2.0, 3.0}), 1e-10 );

            BOOST_REQUIRE_THROW(r_mesh_settings.GetValueFast<Vector3i>(BackgroundGridSettings::polynomial_order), std::exception); // Not set
            r_mesh_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{1, 2, 3});
            const Vector3i& r_polynomial_order = r_mesh_settings.GetValue<Vector3i>(BackgroundGridSettings::polynomial_order);
            QuESo_CHECK_Vector3i_EQUAL(r_polynomial_order, Vector3i({1, 2, 3}) );

            BOOST_REQUIRE_THROW(r_mesh_settings.GetValueFast<Vector3i>(BackgroundGridSettings::number_of_elements), std::exception); // Not set
            r_mesh_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{1, 2, 3});
            const Vector3i& number_of_elements = r_mesh_settings.GetValue<Vector3i>(BackgroundGridSettings::number_of_elements);
            QuESo_CHECK_Vector3i_EQUAL(number_of_elements, Vector3i({1, 2, 3}) );
        }
    }
    {   /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& settings = *p_settingss;

        BOOST_REQUIRE_THROW(StringAccess::GetSubDictionary(settings, "echo_level"), std::exception); // Wrong Key type

        /// General settings
        auto& r_general_settings = StringAccess::GetSubDictionary(settings, "general_settings");
        BOOST_REQUIRE_THROW(StringAccess::IsSet(r_general_settings, "general_settings"), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_general_settings, "lower_bound_uvw"), std::exception); // Wrong Key type

        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_general_settings, "input_filename"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_general_settings, "output_directory_name"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_general_settings, "echo_level"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_general_settings, "echo_level", -1), std::exception); // Negative value
        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_general_settings, "write_output_to_file"), std::exception); // Wrong Value type

        /// Mesh settings
        auto& r_mesh_settings = StringAccess::GetSubDictionary(settings, "background_grid_settings");
        BOOST_REQUIRE_THROW(StringAccess::GetValue<IndexType>(r_mesh_settings, "echo_level"), std::exception); // Wrong Key type

        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_mesh_settings, "lower_bound_xyz"), std::exception); // Not set
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "echo_level", 2u), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "lower_bound_xyz", IndexType(2)), std::exception); // Wrong Value type
        StringAccess::SetValue(r_mesh_settings, "lower_bound_xyz", PointType{1.0, 2.0, 3.0});
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "echo_level", IndexType(2)), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "lower_bound_xyz", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<Vector3i>(r_mesh_settings, "lower_bound_xyz"), std::exception); // Wrong Value type
        const PointType& r_lower_bound_xyz = StringAccess::GetValue<PointType>(r_mesh_settings, "lower_bound_xyz");
        QuESo_CHECK_POINT_NEAR(r_lower_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_mesh_settings, "upper_bound_xyz"), std::exception); // Not set
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "upper_bound_xyz", IndexType(2)), std::exception); // Wrong Value type
        StringAccess::SetValue(r_mesh_settings, "upper_bound_xyz", PointType{1.0, 2.0, 3.0});
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "upper_bound_xyz", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<Vector3i>(r_mesh_settings, "upper_bound_xyz"), std::exception); // Wrong Value type
        const PointType& r_upper_bound_xyz = StringAccess::GetValue<PointType>(r_mesh_settings, "upper_bound_xyz");
        QuESo_CHECK_POINT_NEAR(r_upper_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_mesh_settings, "lower_bound_uvw"), std::exception); // Not set
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "lower_bound_uvw", IndexType(2)), std::exception); // Wrong Value type
        StringAccess::SetValue(r_mesh_settings, "lower_bound_uvw", PointType{1.0, 2.0, 3.0});
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "lower_bound_uvw", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<Vector3i>(r_mesh_settings, "lower_bound_uvw"), std::exception); // Wrong Value type
        const PointType& r_lower_bound_uvw = StringAccess::GetValue<PointType>(r_mesh_settings, "lower_bound_uvw");
        QuESo_CHECK_POINT_NEAR(r_lower_bound_uvw, PointType({1.0, 2.0, 3.0}), 1e-10);

        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_mesh_settings, "upper_bound_uvw"), std::exception); // Not set
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "upper_bound_uvw", IndexType(2)), std::exception); // Wrong Value type
        StringAccess::SetValue(r_mesh_settings, "upper_bound_uvw", PointType{1.0, 2.0, 3.0});
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "upper_bound_uvw", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<Vector3i>(r_mesh_settings, "upper_bound_uvw"), std::exception); // Wrong Value type
        const PointType& r_upper_bound_uvw = StringAccess::GetValue<PointType>(r_mesh_settings, "upper_bound_uvw");
        QuESo_CHECK_POINT_NEAR(r_upper_bound_uvw, PointType({1.0, 2.0, 3.0}), 1e-10 );

        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_mesh_settings, "polynomial_order"), std::exception); // Not set
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "polynomial_order", IndexType(2)), std::exception); // Wrong Value type
        StringAccess::SetValue(r_mesh_settings, "polynomial_order", Vector3i{1, 2, 3});
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "polynomial_order", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_mesh_settings, "polynomial_order"), std::exception); // Wrong Value type
        const Vector3i& r_polynomial_order = StringAccess::GetValue<Vector3i>(r_mesh_settings, "polynomial_order");
        QuESo_CHECK_Vector3i_EQUAL(r_polynomial_order, Vector3i({1, 2, 3}) );

        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_mesh_settings, "number_of_elements"), std::exception); // Not set
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "number_of_elements", IndexType(2)), std::exception); // Wrong Value type
        StringAccess::SetValue(r_mesh_settings, "number_of_elements", Vector3i{1, 2, 3});
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_mesh_settings, "number_of_elements", IndexType(2)), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<PointType>(r_mesh_settings, "number_of_elements"), std::exception); // Wrong Value type
        const Vector3i& number_of_elements = StringAccess::GetValue<Vector3i>(r_mesh_settings, "number_of_elements");
        QuESo_CHECK_Vector3i_EQUAL(number_of_elements, Vector3i({1, 2, 3}) );

        /// Trimmed quadrature rule settings
        auto& r_trimmed_quad_rule_settings = StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings");
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_trimmed_quad_rule_settings, "number_of_elements", Vector3i{1, 2, 3}), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<Vector3i>(r_trimmed_quad_rule_settings, "moment_fitting_residual"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<Vector3i>(r_trimmed_quad_rule_settings, "min_element_volume_ratio"), std::exception); // Wrong Value type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<Vector3i>(r_trimmed_quad_rule_settings, "min_num_boundary_triangles"), std::exception); // Wrong Value type

        /// Non trimmed quadrature rule settings
        auto& r_non_trimmed_quad_rule_settings = StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings");
        BOOST_REQUIRE_THROW(StringAccess::SetValue(r_non_trimmed_quad_rule_settings, "number_of_elements", Vector3i{1, 2, 3}), std::exception); // Wrong Key type
        BOOST_REQUIRE_THROW(StringAccess::GetValue<IndexType>(r_non_trimmed_quad_rule_settings, "integration_method"), std::exception); // Wrong Value type
    }
}


BOOST_AUTO_TEST_CASE(SettingsDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Default Values" << std::endl;

    {   /// Key access
        auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& settings = *p_settingss;

        /// General settings
        QuESo_CHECK( !settings[MainSettings::general_settings].IsSet(GeneralSettings::input_filename) );
        BOOST_REQUIRE_THROW( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::input_filename), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::general_settings].GetValueFast<std::string>(GeneralSettings::input_filename), std::exception );
        }
        QuESo_CHECK( settings[MainSettings::general_settings].IsSet(GeneralSettings::output_directory_name) );
        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<std::string>(GeneralSettings::output_directory_name), std::string("queso_output") );
        if( !NOTDEBUG ) {
            QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValueFast<std::string>(GeneralSettings::output_directory_name), std::string("queso_output") );
        }
        QuESo_CHECK( settings[MainSettings::general_settings].IsSet(GeneralSettings::echo_level) );
        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<IndexType>(GeneralSettings::echo_level), 1u);
        if( !NOTDEBUG ) {
            QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValueFast<IndexType>(GeneralSettings::echo_level), 1u);
        }
        QuESo_CHECK( settings[MainSettings::general_settings].IsSet(GeneralSettings::write_output_to_file) );
        QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValue<bool>(GeneralSettings::write_output_to_file), true);
        if( !NOTDEBUG ) {
            QuESo_CHECK_EQUAL( settings[MainSettings::general_settings].GetValueFast<bool>(GeneralSettings::write_output_to_file), true);
        }
        /// Mesh settings
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::grid_type) );
        BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<GridType>(BackgroundGridSettings::grid_type), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValueFast<GridType>(BackgroundGridSettings::grid_type), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::lower_bound_xyz) );
        BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValueFast<PointType>(BackgroundGridSettings::lower_bound_xyz), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::upper_bound_xyz) );
        BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValueFast<PointType>(BackgroundGridSettings::upper_bound_xyz), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::lower_bound_uvw) );
        BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValueFast<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::upper_bound_uvw) );
        BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValueFast<PointType>(BackgroundGridSettings::upper_bound_uvw), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::polynomial_order) );
        BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::polynomial_order), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValueFast<Vector3i>(BackgroundGridSettings::polynomial_order), std::exception );
        }
        QuESo_CHECK( !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::number_of_elements) );
        BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::number_of_elements), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( settings[MainSettings::background_grid_settings].GetValueFast<Vector3i>(BackgroundGridSettings::number_of_elements), std::exception );
        }

        /// TrimmedQuadratureRuleSettings settings
        QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::moment_fitting_residual) );
        QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual), 1e-10,1e-10 );
        QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValueFast<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual), 1e-10,1e-10 );

        QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::min_element_volume_ratio) );
        QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 1e-3,1e-10 );
        QuESo_CHECK_RELATIVE_NEAR( settings[MainSettings::trimmed_quadrature_rule_settings].GetValueFast<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 1e-3,1e-10 );

        QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::min_num_boundary_triangles) );
        QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 100 );
        QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValueFast<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 100 );

        QuESo_CHECK( settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed) );
        QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValue<bool>(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed), true );
        QuESo_CHECK_EQUAL( settings[MainSettings::trimmed_quadrature_rule_settings].GetValueFast<bool>(TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed), true );

        // NonTrimmedQuadratureRuleSettings settings
        QuESo_CHECK( settings[MainSettings::non_trimmed_quadrature_rule_settings].IsSet(NonTrimmedQuadratureRuleSettings::integration_method) );
        QuESo_CHECK_EQUAL( settings[MainSettings::non_trimmed_quadrature_rule_settings].GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::gauss );
        QuESo_CHECK_EQUAL( settings[MainSettings::non_trimmed_quadrature_rule_settings].GetValueFast<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::gauss );
    }
    {   /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& settings = *p_settingss;

        /// General settings
        QuESo_CHECK( !StringAccess::IsSet( StringAccess::GetSubDictionary(settings, "general_settings"), "input_filename" ));
        BOOST_REQUIRE_THROW( StringAccess::GetValue<std::string>(StringAccess::GetSubDictionary(settings, "general_settings"), "input_filename"), std::exception );

        QuESo_CHECK( StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "general_settings"), "output_directory_name") );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<std::string>(StringAccess::GetSubDictionary(settings, "general_settings"), "output_directory_name"), std::string("queso_output") );

        QuESo_CHECK( StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "general_settings"), "echo_level") );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(StringAccess::GetSubDictionary(settings, "general_settings"), "echo_level"), 1u);

        QuESo_CHECK( StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "general_settings"), "write_output_to_file") );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<bool>(StringAccess::GetSubDictionary(settings, "general_settings"), "write_output_to_file"), true);

        /// Mesh settings
        QuESo_CHECK( !StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "grid_type") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "grid_type"), std::exception );

        QuESo_CHECK( !StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz"), std::exception );

        QuESo_CHECK( !StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_xyz") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_xyz"), std::exception );

        QuESo_CHECK( !StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_uvw") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_uvw"), std::exception );

        QuESo_CHECK( !StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_uvw") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_uvw"), std::exception );

        QuESo_CHECK( !StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "polynomial_order") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<Vector3i>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "polynomial_order"), std::exception );

        QuESo_CHECK( !StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "number_of_elements") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<Vector3i>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "number_of_elements"), std::exception );

        /// TrimmedQuadratureRuleSettings settings
        QuESo_CHECK( StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual") );
        QuESo_CHECK_RELATIVE_NEAR( StringAccess::GetValue<double>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual"), 1e-10,1e-10 );

        QuESo_CHECK( StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "min_element_volume_ratio") );
        QuESo_CHECK_RELATIVE_NEAR( StringAccess::GetValue<double>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "min_element_volume_ratio"), 1e-3,1e-10 );

        QuESo_CHECK( StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "min_num_boundary_triangles") );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "min_num_boundary_triangles"), 100 );

        QuESo_CHECK( StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "neglect_elements_if_stl_is_flawed") );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<bool>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "neglect_elements_if_stl_is_flawed"), true );

        // NonTrimmedQuadratureRuleSettings settings
        QuESo_CHECK( StringAccess::IsSet(StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings"), "integration_method") );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<IntegrationMethod>(StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings"), "integration_method"), IntegrationMethod::gauss );
    }
}


BOOST_AUTO_TEST_CASE(SettingsCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Customized Values" << std::endl;

    {   /// Key access
        auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& settings = *p_settingss;

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
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& settings = *p_settingss;

        /// General settings
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "general_settings"), "input_filename", std::string("test_filename.stl"));

        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "general_settings"), "output_directory_name", std::string("new_output/"));
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "general_settings"), "echo_level", 2u);
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "general_settings"), "write_output_to_file", false);

        QuESo_CHECK_EQUAL( StringAccess::GetValue<std::string>(StringAccess::GetSubDictionary(settings, "general_settings"), "input_filename"), std::string("test_filename.stl") );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<std::string>(StringAccess::GetSubDictionary(settings, "general_settings"), "output_directory_name"), std::string("new_output/") );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(StringAccess::GetSubDictionary(settings, "general_settings"), "echo_level"), 2u );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<bool>(StringAccess::GetSubDictionary(settings, "general_settings"), "write_output_to_file"), false );

        /// Mesh settings
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "background_grid_settings"), "grid_type", GridType::b_spline_grid);
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz", PointType({1.0, 1.0, 2.0}));
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_xyz", PointType({0.0, 3.0, 2.0}));
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_uvw", PointType({0.1, -1.0, 2.0}));
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_uvw", PointType({0.66, 1.0, 2.2}));
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "background_grid_settings"), "polynomial_order", Vector3i({5, 6, 7}));
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "background_grid_settings"), "number_of_elements", Vector3i({8, 9, 2}));

        QuESo_CHECK_EQUAL( StringAccess::GetValue<GridType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "grid_type"), GridType::b_spline_grid );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz"), PointType({1.0, 1.0, 2.0}) );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_xyz"), PointType({0.0, 3.0, 2.0}) );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_uvw"), PointType({0.1, -1.0, 2.0}) );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_uvw"), PointType({0.66, 1.0, 2.2}) );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<Vector3i>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "polynomial_order"), Vector3i({5, 6, 7}) );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<Vector3i>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "number_of_elements"), Vector3i({8, 9, 2}) );

        /// Trimmed quadrature rule settings
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual", 5.6e-5);
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "min_element_volume_ratio", 0.45);
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "min_num_boundary_triangles", IndexType(234));
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "neglect_elements_if_stl_is_flawed", false);

        QuESo_CHECK_RELATIVE_NEAR( StringAccess::GetValue<double>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual"), 5.6e-5, 1e-10 );
        QuESo_CHECK_RELATIVE_NEAR( StringAccess::GetValue<double>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "min_element_volume_ratio"), 0.45, 1e-10 );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "min_num_boundary_triangles"), 234u );
        QuESo_CHECK_EQUAL( StringAccess::GetValue<bool>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "neglect_elements_if_stl_is_flawed"), false );

        /// Non trimmed quadrature rule settings
        StringAccess::SetValue( StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings"), "integration_method", IntegrationMethod::ggq_optimal);
        QuESo_CHECK_EQUAL( StringAccess::GetValue<IntegrationMethod>(StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings"), "integration_method"), IntegrationMethod::ggq_optimal );
    }
}

BOOST_AUTO_TEST_CASE(SettingsConditionSettingsWrongTypeTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Condition Settings Wrong Types " << std::endl;

    {   /// Key access
        auto p_condition_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
        auto& cond_settings = *p_condition_settings;

        BOOST_REQUIRE_THROW( cond_settings.GetValue<IndexType>(ConditionSettings::condition_id), std::exception ); // Not set
        BOOST_REQUIRE_THROW( cond_settings.GetValue<IndexType>(ConditionSettings::condition_id), std::exception ); // Not set
        BOOST_REQUIRE_THROW( cond_settings.GetValue<std::string>(ConditionSettings::condition_type), std::exception ); // Not set
        BOOST_REQUIRE_THROW( cond_settings.GetValue<std::string>(ConditionSettings::input_filename), std::exception ); // Not set
        BOOST_REQUIRE_THROW( cond_settings.GetValue<double>(ConditionSettings::modulus), std::exception ); // Not set
        BOOST_REQUIRE_THROW( cond_settings.GetValue<PointType>(ConditionSettings::direction), std::exception ); // Not set
        BOOST_REQUIRE_THROW( cond_settings.GetValue<PointType>(ConditionSettings::value), std::exception ); // Not set
        BOOST_REQUIRE_THROW( cond_settings.GetValue<double>(ConditionSettings::penalty_factor), std::exception ); // Not set

        if( !NOTDEBUG ) {

            BOOST_REQUIRE_THROW( cond_settings.IsSet(BackgroundGridSettings::lower_bound_uvw), std::exception ); // Wrong Key

            BOOST_REQUIRE_THROW( cond_settings.GetValueFast<IndexType>(ConditionSettings::condition_id), std::exception ); // Not set
            BOOST_REQUIRE_THROW( cond_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType({2.0, 3.0, 4.0})), std::exception ); // Wrong Key type

            BOOST_REQUIRE_THROW( cond_settings.GetValueFast<IndexType>(ConditionSettings::condition_id), std::exception ); // Not set
            cond_settings.SetValue(ConditionSettings::condition_id, 100u);
            BOOST_REQUIRE_THROW( cond_settings.GetValue<IndexType>(GeneralSettings::echo_level), std::exception ); // Wrong Key type
            BOOST_REQUIRE_THROW( cond_settings.GetValueFast<IndexType>(GeneralSettings::echo_level), std::exception ); // Wrong Key type

            BOOST_REQUIRE_THROW( cond_settings.GetValue<std::string>(ConditionSettings::condition_type), std::exception ); // Not set
            cond_settings.SetValue(ConditionSettings::condition_type, std::string("dummy"));
            BOOST_REQUIRE_THROW( cond_settings.GetValue<IndexType>(GeneralSettings::echo_level), std::exception ); // Wrong Key type
            BOOST_REQUIRE_THROW( cond_settings.GetValueFast<IndexType>(GeneralSettings::echo_level), std::exception ); // Wrong Key type

            BOOST_REQUIRE_THROW( cond_settings.GetValueFast<std::string>(ConditionSettings::input_filename), std::exception ); // Not set
            BOOST_REQUIRE_THROW( cond_settings.GetValueFast<double>(ConditionSettings::modulus), std::exception ); // Not set
            BOOST_REQUIRE_THROW( cond_settings.GetValueFast<PointType>(ConditionSettings::direction), std::exception ); // Not set
            BOOST_REQUIRE_THROW( cond_settings.GetValueFast<PointType>(ConditionSettings::value), std::exception ); // Not set
            BOOST_REQUIRE_THROW( cond_settings.GetValueFast<double>(ConditionSettings::penalty_factor), std::exception ); // Not set
        }
    }
    {   /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_condition_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
        auto& cond_settings = *p_condition_settings;

        BOOST_REQUIRE_THROW( StringAccess::IsSet(cond_settings, "lower_bound_uvw"), std::exception ); // Wrong Key

        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "condition_id"), std::exception ); // Not set
        BOOST_REQUIRE_THROW( StringAccess::SetValue(cond_settings, "lower_bound_uvw", PointType({2.0, 3.0, 4.0})), std::exception ); // Wrong Key type
        BOOST_REQUIRE_THROW( StringAccess::SetValue(cond_settings, "condition_id", 2.0), std::exception ); // Wrong Value type

        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "condition_id"), std::exception ); // Not set
        StringAccess::SetValue(cond_settings, "condition_id", 100u);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "echo_level"), std::exception ); // Wrong Key type
        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(cond_settings, "condition_id"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( StringAccess::GetValue<std::string>(cond_settings, "condition_type"), std::exception ); // Not set
        StringAccess::SetValue(cond_settings, "condition_type", std::string("dummy"));
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "echo_level"), std::exception ); // Wrong Key type
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "condition_type"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( StringAccess::GetValue<std::string>(cond_settings, "input_filename"), std::exception ); // Not set
        StringAccess::SetValue(cond_settings, "input_filename", std::string("dummy_2"));
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "input_filename"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(cond_settings, "modulus"), std::exception ); // Not set
        StringAccess::SetValue(cond_settings, "modulus", 2.0);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "modulus"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(cond_settings, "direction"), std::exception ); // Not set
        StringAccess::SetValue(cond_settings, "direction", PointType({2.2, 3.0, 4.4}));
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "direction"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(cond_settings, "value"), std::exception ); // Not set
        StringAccess::SetValue(cond_settings, "value", PointType({2.1, 1.0, 2.4}));
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "value"), std::exception ); // Wrong Type

        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(cond_settings, "penalty_factor"), std::exception ); // Not set
        StringAccess::SetValue(cond_settings, "penalty_factor", 1e5);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(cond_settings, "penalty_factor"), std::exception ); // Wrong Type
    }
}


BOOST_AUTO_TEST_CASE(SettingsConditionSettingsDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Condition Settings Defaut Values" << std::endl;

    {   /// Key access
        auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
        auto& cond_settings = *p_cond_settings;

        QuESo_CHECK( !cond_settings.IsSet(ConditionSettings::condition_id) );
        QuESo_CHECK( !cond_settings.IsSet(ConditionSettings::condition_type) );
        QuESo_CHECK( !cond_settings.IsSet(ConditionSettings::input_filename) );
        QuESo_CHECK( !cond_settings.IsSet(ConditionSettings::modulus) );
        QuESo_CHECK( !cond_settings.IsSet(ConditionSettings::direction) );
        QuESo_CHECK( !cond_settings.IsSet(ConditionSettings::value) );
        QuESo_CHECK( !cond_settings.IsSet(ConditionSettings::penalty_factor) );
    }
    {   /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
        auto& cond_settings = *p_cond_settings;

        QuESo_CHECK( !StringAccess::IsSet(cond_settings, "condition_id") );
        QuESo_CHECK( !StringAccess::IsSet(cond_settings, "condition_type") );
        QuESo_CHECK( !StringAccess::IsSet(cond_settings, "input_filename") );
        QuESo_CHECK( !StringAccess::IsSet(cond_settings, "modulus") );
        QuESo_CHECK( !StringAccess::IsSet(cond_settings, "direction") );
        QuESo_CHECK( !StringAccess::IsSet(cond_settings, "value") );
        QuESo_CHECK( !StringAccess::IsSet(cond_settings, "penalty_factor") );
    }
};

BOOST_AUTO_TEST_CASE(SettingsConditionSettingsCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Condition Settings Customized Values" << std::endl;

    {   /// Key access
        auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& r_cond_setting_list = p_settingss->GetList(MainSettings::conditions_settings_list);

        {
            auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");

            p_cond_settings->SetValue(ConditionSettings::condition_id, 120UL);
            p_cond_settings->SetValue(ConditionSettings::condition_type, std::string("Hello"));
            p_cond_settings->SetValue(ConditionSettings::input_filename, std::string("hallo"));
            p_cond_settings->SetValue(ConditionSettings::modulus, 200.0);
            p_cond_settings->SetValue(ConditionSettings::penalty_factor, 300.0);

            r_cond_setting_list.push_back(std::move(p_cond_settings));
        }
        {
            auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");

            p_cond_settings->SetValue(ConditionSettings::condition_id, 150UL);
            p_cond_settings->SetValue(ConditionSettings::condition_type, std::string("Hello2"));
            p_cond_settings->SetValue(ConditionSettings::input_filename, std::string("hallo2"));
            p_cond_settings->SetValue(ConditionSettings::direction, PointType{2.1, 3.2, 4.7});
            p_cond_settings->SetValue(ConditionSettings::value, PointType{2.2, 2.3, 1.4} );

            r_cond_setting_list.push_back(std::move(p_cond_settings));
        }
        {
            auto& r_cond_settings_list = p_settingss->GetList(MainSettings::conditions_settings_list);
            QuESo_CHECK_EQUAL( r_cond_settings_list.size(), 2 );

            QuESo_CHECK( r_cond_settings_list[0u]->IsSet(ConditionSettings::condition_id) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u]->GetValue<IndexType>(ConditionSettings::condition_id), 120u );

            QuESo_CHECK( r_cond_settings_list[0u]->IsSet(ConditionSettings::condition_type) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u]->GetValue<std::string>(ConditionSettings::condition_type), std::string("Hello") );

            QuESo_CHECK( r_cond_settings_list[0u]->IsSet(ConditionSettings::input_filename) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u]->GetValue<std::string>(ConditionSettings::input_filename), std::string("hallo") );

            QuESo_CHECK( r_cond_settings_list[0u]->IsSet(ConditionSettings::modulus) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u]->GetValue<double>(ConditionSettings::modulus), 200.0 );

            QuESo_CHECK( !r_cond_settings_list[0u]->IsSet(ConditionSettings::direction) );
            BOOST_REQUIRE_THROW( r_cond_settings_list[0u]->GetValue<PointType>(ConditionSettings::direction), std::exception );
            if( !NOTDEBUG )
                BOOST_REQUIRE_THROW( r_cond_settings_list[0u]->GetValueFast<PointType>(ConditionSettings::direction), std::exception );

            QuESo_CHECK( !r_cond_settings_list[0u]->IsSet(ConditionSettings::value) );
            BOOST_REQUIRE_THROW( r_cond_settings_list[0u]->GetValue<PointType>(ConditionSettings::value), std::exception );
            if( !NOTDEBUG )
                BOOST_REQUIRE_THROW( r_cond_settings_list[0u]->GetValueFast<PointType>(ConditionSettings::value), std::exception );

            QuESo_CHECK( r_cond_settings_list[0u]->IsSet(ConditionSettings::penalty_factor) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[0u]->GetValue<double>(ConditionSettings::penalty_factor), 300.0 );
        }
        {
            auto& r_cond_settings_list = p_settingss->GetList(MainSettings::conditions_settings_list);
            QuESo_CHECK_EQUAL( r_cond_settings_list.size(), 2 );

            QuESo_CHECK( r_cond_settings_list[1u]->IsSet(ConditionSettings::condition_id) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[1u]->GetValue<IndexType>(ConditionSettings::condition_id), 150u );

            QuESo_CHECK( r_cond_settings_list[1u]->IsSet(ConditionSettings::condition_type) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[1u]->GetValue<std::string>(ConditionSettings::condition_type), std::string("Hello2") );

            QuESo_CHECK( r_cond_settings_list[1u]->IsSet(ConditionSettings::input_filename) );
            QuESo_CHECK_EQUAL( r_cond_settings_list[1u]->GetValue<std::string>(ConditionSettings::input_filename), std::string("hallo2") );

            QuESo_CHECK( !r_cond_settings_list[1u]->IsSet(ConditionSettings::modulus) );
            BOOST_REQUIRE_THROW( r_cond_settings_list[1u]->GetValue<double>(ConditionSettings::modulus), std::exception );
            if( !NOTDEBUG )
                BOOST_REQUIRE_THROW( r_cond_settings_list[1u]->GetValueFast<double>(ConditionSettings::modulus), std::exception );

            QuESo_CHECK( r_cond_settings_list[1u]->IsSet(ConditionSettings::direction) );
            QuESo_CHECK_POINT_NEAR( r_cond_settings_list[1u]->GetValue<PointType>(ConditionSettings::direction), PointType({2.1, 3.2, 4.7}), 1e-10 );

            QuESo_CHECK( r_cond_settings_list[1u]->IsSet(ConditionSettings::value) );
            QuESo_CHECK_POINT_NEAR( r_cond_settings_list[1u]->GetValue<PointType>(ConditionSettings::value), PointType({2.2, 2.3, 1.4}), 1e-10 );

            QuESo_CHECK( !r_cond_settings_list[1u]->IsSet(ConditionSettings::penalty_factor) );
            BOOST_REQUIRE_THROW( r_cond_settings_list[1u]->GetValue<double>(ConditionSettings::penalty_factor), std::exception );
            if( !NOTDEBUG )
                BOOST_REQUIRE_THROW( r_cond_settings_list[1u]->GetValueFast<double>(ConditionSettings::penalty_factor), std::exception );
        }
    }
    {   /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
        auto& settings = *p_settingss;
        auto& r_cond_setting_list = StringAccess::GetList(settings, "conditions_settings_list");

        {
            auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
            StringAccess::SetValue(*p_cond_settings, "condition_id", 120UL);
            StringAccess::SetValue(*p_cond_settings, "condition_type", std::string("Hello"));
            StringAccess::SetValue(*p_cond_settings, "input_filename", std::string("hallo"));
            StringAccess::SetValue(*p_cond_settings, "modulus", 200.0);
            StringAccess::SetValue(*p_cond_settings, "penalty_factor", 300.0);

            r_cond_setting_list.push_back(std::move(p_cond_settings));
        }
        {
            auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
            StringAccess::SetValue(*p_cond_settings, "condition_id", 150UL);
            StringAccess::SetValue(*p_cond_settings, "condition_type", std::string("Hello2"));
            StringAccess::SetValue(*p_cond_settings, "input_filename", std::string("hallo2"));
            StringAccess::SetValue(*p_cond_settings, "direction", PointType{2.1, 3.2, 4.7});
            StringAccess::SetValue(*p_cond_settings, "value", PointType{2.2, 2.3, 1.4});

            r_cond_setting_list.push_back(std::move(p_cond_settings));
        }
        {
            auto& r_cond_settings_list = StringAccess::GetList(settings, "conditions_settings_list");
            QuESo_CHECK_EQUAL( r_cond_settings_list.size(), 2 );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[0u], "condition_id") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(*r_cond_settings_list[0u], "condition_id"), 120u );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[0u], "condition_type") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<std::string>(*r_cond_settings_list[0u], "condition_type"), std::string("Hello") );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[0u], "input_filename") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<std::string>(*r_cond_settings_list[0u], "input_filename"), std::string("hallo") );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[0u], "modulus") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<double>(*r_cond_settings_list[0u], "modulus"), 200.0 );

            QuESo_CHECK( !StringAccess::IsSet(*r_cond_settings_list[0u], "direction") );
            BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(*r_cond_settings_list[0u], "direction"), std::exception );

            QuESo_CHECK( !StringAccess::IsSet(*r_cond_settings_list[0u], "value") );
            BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(*r_cond_settings_list[0u], "value"), std::exception );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[0u], "penalty_factor") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<double>(*r_cond_settings_list[0u], "penalty_factor"), 300.0 );
        }
        {
            auto& r_cond_settings_list = StringAccess::GetList(settings, "conditions_settings_list");
            QuESo_CHECK_EQUAL( r_cond_settings_list.size(), 2 );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[1u], "condition_id") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(*r_cond_settings_list[1u], "condition_id"), 150u );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[1u], "condition_type") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<std::string>(*r_cond_settings_list[1u], "condition_type"), std::string("Hello2") );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[1u], "input_filename") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<std::string>(*r_cond_settings_list[1u], "input_filename"), std::string("hallo2") );

            QuESo_CHECK( !StringAccess::IsSet(*r_cond_settings_list[1u], "modulus") );
            BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(*r_cond_settings_list[1u], "modulus"), std::exception );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[1u], "direction") );
            QuESo_CHECK_POINT_NEAR( StringAccess::GetValue<PointType>(*r_cond_settings_list[1u], "direction"), PointType({2.1, 3.2, 4.7}), 1e-10 );

            QuESo_CHECK( StringAccess::IsSet(*r_cond_settings_list[1u], "value") );
            QuESo_CHECK_POINT_NEAR( StringAccess::GetValue<PointType>(*r_cond_settings_list[1u], "value"), PointType({2.2, 2.3, 1.4}), 1e-10 );

            QuESo_CHECK( !StringAccess::IsSet(*r_cond_settings_list[1u], "penalty_factor") );
            BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(*r_cond_settings_list[1u], "penalty_factor"), std::exception );
        }
    }
};


BOOST_AUTO_TEST_CASE(SettingsCastAmbiguousTypesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Cast Ambiguous Types" << std::endl;

    using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
    using StringAccess = DictionaryStringAccess<DictionaryType>;

    auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& settings = *p_settingss;

    // Lets set it with wrong type
    StringAccess::SetValue(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual", 5u);
    {
        const double r_double = StringAccess::GetValue<double>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual");
        QuESo_CHECK_RELATIVE_NEAR(r_double, 5.0, 1e-10);
    }
    // Lets set it with correct type
    StringAccess::SetValue(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual", 7.7);
    {
        const double r_double = StringAccess::GetValue<double>(StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual");
        QuESo_CHECK_RELATIVE_NEAR(r_double, 7.7, 1e-10);
    }

    // Cast Vector3i to Vector3d
    BOOST_REQUIRE_THROW( StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz"), std::exception ); // Not set

    // Lets set it with wrong type
    StringAccess::SetValue(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz", Vector3i({5, 3, 2}));
    {
        const PointType& r_point = StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz");
        QuESo_CHECK_POINT_RELATIVE_NEAR(r_point, PointType({5.0, 3.0, 2.0}), 1e-10);
    }
    // Lets set it with correct type
    StringAccess::SetValue(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz", Vector3d({5.5, 2.2, 1.1}));
    {
        const PointType& r_point = StringAccess::GetValue<PointType>(StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz");
        QuESo_CHECK_POINT_RELATIVE_NEAR(r_point, PointType({5.5, 2.2, 1.1}), 1e-10);
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

