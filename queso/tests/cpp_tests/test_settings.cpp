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

/**
 * @class  Container to store all global parameters.
 * @author Manuel Messmer
**/

namespace Testing {


BOOST_AUTO_TEST_SUITE( SettingsTestSuite )

BOOST_AUTO_TEST_CASE(ParameterCustomConstructorTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Wrong Types" << std::endl;
    Settings setting;

    BOOST_REQUIRE_THROW(setting[GeneralSettings::echo_level], std::exception); // Wrong Key type

    /// General settings
    auto& r_general_settings = setting[Main::general_settings];
    BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(MeshSettings::lower_bound_uvw), std::exception); // Wrong Key type

    BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::input_filename), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::output_directory_name), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::echo_level), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_general_settings.GetValue<PointType>(GeneralSettings::embedding_flag), std::exception); // Wrong Value type

    /// Mesh settings
    auto& r_mesh_settings = setting[Main::mesh_settings];
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<IndexType>(GeneralSettings::echo_level), std::exception); // Wrong Key type

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(MeshSettings::lower_bound_xyz), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(GeneralSettings::echo_level, IndexType(2)), std::exception); // Wrong Key type
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::lower_bound_xyz, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(MeshSettings::lower_bound_xyz, PointType{1.0, 2.0, 3.0});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(GeneralSettings::echo_level, IndexType(2)), std::exception); // Wrong Key type
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::lower_bound_xyz, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>(MeshSettings::lower_bound_xyz), std::exception); // Wrong Value type
    const PointType& r_lower_bound_xyz = r_mesh_settings.GetValue<PointType>(MeshSettings::lower_bound_xyz);
    QuESo_CHECK_POINT_NEAR(r_lower_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(MeshSettings::upper_bound_xyz), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::upper_bound_xyz, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(MeshSettings::upper_bound_xyz, PointType{1.0, 2.0, 3.0});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::upper_bound_xyz, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>(MeshSettings::upper_bound_xyz), std::exception); // Wrong Value type
    const PointType& r_upper_bound_xyz = r_mesh_settings.GetValue<PointType>(MeshSettings::upper_bound_xyz);
    QuESo_CHECK_POINT_NEAR(r_upper_bound_xyz, PointType({1.0, 2.0, 3.0}), 1e-10);

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(MeshSettings::lower_bound_uvw), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::lower_bound_uvw, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(MeshSettings::lower_bound_uvw, PointType{1.0, 2.0, 3.0});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::lower_bound_uvw, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>(MeshSettings::lower_bound_uvw), std::exception); // Wrong Value type
    const PointType& r_lower_bound_uvw= r_mesh_settings.GetValue<PointType>(MeshSettings::lower_bound_uvw);
    QuESo_CHECK_POINT_NEAR(r_lower_bound_uvw, PointType({1.0, 2.0, 3.0}), 1e-10);

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(MeshSettings::upper_bound_uvw), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::upper_bound_uvw, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(MeshSettings::upper_bound_uvw, PointType{1.0, 2.0, 3.0});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::upper_bound_uvw, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<Vector3i>(MeshSettings::upper_bound_uvw), std::exception); // Wrong Value type
    const PointType& r_upper_bound_uvw = r_mesh_settings.GetValue<PointType>(MeshSettings::upper_bound_uvw);
    QuESo_CHECK_Vector3i_EQUAL(r_upper_bound_uvw, Vector3i({1, 2, 3}) );

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(MeshSettings::polynomial_order), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::polynomial_order, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(MeshSettings::polynomial_order, Vector3i{1, 2, 3});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::polynomial_order, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(MeshSettings::polynomial_order), std::exception); // Wrong Value type
    const Vector3i& r_polynomial_order = r_mesh_settings.GetValue<Vector3i>(MeshSettings::polynomial_order);
    QuESo_CHECK_Vector3i_EQUAL(r_polynomial_order, Vector3i({1, 2, 3}) );

    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(MeshSettings::number_of_elements), std::exception); // Not set
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::number_of_elements, IndexType(2)), std::exception); // Wrong Value type
    r_mesh_settings.SetValue(MeshSettings::number_of_elements, Vector3i{1, 2, 3});
    BOOST_REQUIRE_THROW(r_mesh_settings.SetValue(MeshSettings::number_of_elements, IndexType(2)), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_mesh_settings.GetValue<PointType>(MeshSettings::number_of_elements), std::exception); // Wrong Value type
    const Vector3i& number_of_elements = r_mesh_settings.GetValue<Vector3i>(MeshSettings::number_of_elements);
    QuESo_CHECK_Vector3i_EQUAL(number_of_elements, Vector3i({1, 2, 3}) );

    /// Trimmed quadrature rule settings
    auto& r_trimmed_quad_rule_settings = setting[Main::trimmed_quadrature_rule_settings];
    BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.SetValue(MeshSettings::number_of_elements, Vector3i{1, 2, 3}), std::exception); // Wrong Key type
    BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>(TrimmedQuadratureRuleSettings::moment_fitting_residual), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), std::exception); // Wrong Value type
    BOOST_REQUIRE_THROW(r_trimmed_quad_rule_settings.GetValue<Vector3i>(TrimmedQuadratureRuleSettings::use_customized_trimmed_points), std::exception); // Wrong Value type

    /// Non trimmed quadrature rule settings
    auto& r_non_trimmed_quad_rule_settings = setting[Main::non_trimmed_quadrature_rule_settings];
    BOOST_REQUIRE_THROW(r_non_trimmed_quad_rule_settings.SetValue(MeshSettings::number_of_elements, Vector3i{1, 2, 3}), std::exception); // Wrong Key type
    BOOST_REQUIRE_THROW(r_non_trimmed_quad_rule_settings.GetValue<IndexType>(NonTrimmedQuadratureRuleSettings::integration_method), std::exception); // Wrong Value type
}


BOOST_AUTO_TEST_CASE(SettingsDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Default Values" << std::endl;

    Settings settings;

    /// General settings
    QuESo_CHECK( !settings[Main::general_settings].IsSet(GeneralSettings::input_filename) );
    BOOST_REQUIRE_THROW( settings[Main::general_settings].GetValue<std::string>(GeneralSettings::input_filename), std::exception );

    QuESo_CHECK( settings[Main::general_settings].IsSet(GeneralSettings::output_directory_name) );
    QuESo_CHECK_EQUAL( settings[Main::general_settings].GetValue<std::string>(GeneralSettings::output_directory_name), std::string("queso_output") );

    QuESo_CHECK( settings[Main::general_settings].IsSet(GeneralSettings::echo_level) );
    QuESo_CHECK_EQUAL( settings[Main::general_settings].GetValue<IndexType>(GeneralSettings::echo_level), 1UL);

    QuESo_CHECK( settings[Main::general_settings].IsSet(GeneralSettings::embedding_flag) );
    QuESo_CHECK_EQUAL( settings[Main::general_settings].GetValue<bool>(GeneralSettings::embedding_flag), true);

    /// Mesh settings
    QuESo_CHECK( !settings[Main::mesh_settings].IsSet(MeshSettings::lower_bound_xyz) );
    BOOST_REQUIRE_THROW( settings[Main::mesh_settings].GetValue<PointType>(MeshSettings::lower_bound_xyz), std::exception );

    QuESo_CHECK( !settings[Main::mesh_settings].IsSet(MeshSettings::upper_bound_xyz) );
    BOOST_REQUIRE_THROW( settings[Main::mesh_settings].GetValue<PointType>(MeshSettings::upper_bound_xyz), std::exception );

    QuESo_CHECK( !settings[Main::mesh_settings].IsSet(MeshSettings::lower_bound_uvw) );
    BOOST_REQUIRE_THROW( settings[Main::mesh_settings].GetValue<PointType>(MeshSettings::lower_bound_uvw), std::exception );

    QuESo_CHECK( !settings[Main::mesh_settings].IsSet(MeshSettings::upper_bound_uvw) );
    BOOST_REQUIRE_THROW( settings[Main::mesh_settings].GetValue<PointType>(MeshSettings::upper_bound_uvw), std::exception );

    QuESo_CHECK( !settings[Main::mesh_settings].IsSet(MeshSettings::polynomial_order) );
    BOOST_REQUIRE_THROW( settings[Main::mesh_settings].GetValue<Vector3i>(MeshSettings::polynomial_order), std::exception );

    QuESo_CHECK( !settings[Main::mesh_settings].IsSet(MeshSettings::number_of_elements) );
    BOOST_REQUIRE_THROW( settings[Main::mesh_settings].GetValue<Vector3i>(MeshSettings::number_of_elements), std::exception );

    /// TrimmedQuadratureRuleSettings settings
    QuESo_CHECK( settings[Main::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::moment_fitting_residual) );
    QuESo_CHECK_RELATIVE_NEAR( settings[Main::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual), 1e-10,1e-10 );

    QuESo_CHECK( settings[Main::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::min_element_volume_ratio) );
    QuESo_CHECK_RELATIVE_NEAR( settings[Main::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 1e-3,1e-10 );

    QuESo_CHECK( settings[Main::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::min_num_boundary_triangles) );
    QuESo_CHECK_EQUAL( settings[Main::trimmed_quadrature_rule_settings].GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 100 );

    QuESo_CHECK( settings[Main::trimmed_quadrature_rule_settings].IsSet(TrimmedQuadratureRuleSettings::use_customized_trimmed_points) );
    QuESo_CHECK_EQUAL( settings[Main::trimmed_quadrature_rule_settings].GetValue<bool>(TrimmedQuadratureRuleSettings::use_customized_trimmed_points), false );

    // NonTrimmedQuadratureRuleSettings settings
    QuESo_CHECK( settings[Main::non_trimmed_quadrature_rule_settings].IsSet(NonTrimmedQuadratureRuleSettings::integration_method) );
    QuESo_CHECK_EQUAL( settings[Main::non_trimmed_quadrature_rule_settings].GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::Gauss );
}

BOOST_AUTO_TEST_CASE(SettingsCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test Settings :: Test Customized Values" << std::endl;

    Settings settings;

    /// General settings
    settings[Main::general_settings].SetValue(GeneralSettings::input_filename, std::string("test_filename.stl"));
    settings[Main::general_settings].SetValue(GeneralSettings::output_directory_name, std::string("new_output/"));
    settings[Main::general_settings].SetValue(GeneralSettings::echo_level, 2UL);
    settings[Main::general_settings].SetValue(GeneralSettings::embedding_flag, false);


    QuESo_CHECK_EQUAL( settings[Main::general_settings].GetValue<std::string>(GeneralSettings::input_filename), std::string("test_filename.stl") );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::general_settings).GetValueNoCheck<std::string>(GeneralSettings::input_filename),
        std::string("test_filename.stl") );

    QuESo_CHECK_EQUAL( settings[Main::general_settings].GetValue<std::string>(GeneralSettings::output_directory_name), std::string("new_output/") );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::general_settings).GetValueNoCheck<std::string>(GeneralSettings::output_directory_name),
        std::string("new_output/") );

    QuESo_CHECK_EQUAL( settings[Main::general_settings].GetValue<IndexType>(GeneralSettings::echo_level), 2Ul );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::general_settings).GetValueNoCheck<IndexType>(GeneralSettings::echo_level), 2UL );

    QuESo_CHECK_EQUAL( settings[Main::general_settings].GetValue<bool>(GeneralSettings::embedding_flag), false);
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::general_settings).GetValueNoCheck<bool>(GeneralSettings::embedding_flag), false );

    /// Mesh settings
    settings[Main::mesh_settings].SetValue(MeshSettings::lower_bound_xyz, PointType({1.0, 1.0, 2.0}));
    settings[Main::mesh_settings].SetValue(MeshSettings::upper_bound_xyz, PointType({0.0, 3.0, 2.0}));
    settings[Main::mesh_settings].SetValue(MeshSettings::lower_bound_uvw, PointType({0.1, -1.0, 2.0}));
    settings[Main::mesh_settings].SetValue(MeshSettings::upper_bound_uvw, PointType({0.66, 1.0, 2.2}));
    settings[Main::mesh_settings].SetValue(MeshSettings::polynomial_order, Vector3i({5, 6, 7}));
    settings[Main::mesh_settings].SetValue(MeshSettings::number_of_elements, Vector3i({8, 9, 2}));

    QuESo_CHECK_EQUAL( settings[Main::mesh_settings].GetValue<PointType>(MeshSettings::lower_bound_xyz), PointType({1.0, 1.0, 2.0}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::mesh_settings).GetValueNoCheck<PointType>(MeshSettings::lower_bound_xyz),
        PointType({1.0, 1.0, 2.0}) );

    QuESo_CHECK_EQUAL( settings[Main::mesh_settings].GetValue<PointType>(MeshSettings::upper_bound_xyz), PointType({0.0, 3.0, 2.0}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::mesh_settings).GetValueNoCheck<PointType>(MeshSettings::upper_bound_xyz),
        PointType({0.0, 3.0, 2.0}) );

    QuESo_CHECK_EQUAL( settings[Main::mesh_settings].GetValue<PointType>(MeshSettings::lower_bound_uvw), PointType({0.1, -1.0, 2.0}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::mesh_settings).GetValueNoCheck<PointType>(MeshSettings::lower_bound_uvw),
        PointType({0.1, -1.0, 2.0}) );

    QuESo_CHECK_EQUAL( settings[Main::mesh_settings].GetValue<PointType>(MeshSettings::upper_bound_uvw), PointType({0.66, 1.0, 2.2}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::mesh_settings).GetValueNoCheck<PointType>(MeshSettings::upper_bound_uvw),
        PointType({0.66, 1.0, 2.2}) );

    QuESo_CHECK_EQUAL( settings[Main::mesh_settings].GetValue<Vector3i>(MeshSettings::polynomial_order), Vector3i({5, 6, 7}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::mesh_settings).GetValueNoCheck<Vector3i>(MeshSettings::polynomial_order),
        Vector3i({5, 6, 7}) );

    QuESo_CHECK_EQUAL( settings[Main::mesh_settings].GetValue<Vector3i>(MeshSettings::number_of_elements), Vector3i({8, 9, 2}) );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::mesh_settings).GetValueNoCheck<Vector3i>(MeshSettings::number_of_elements),
        Vector3i({8, 9, 2}) );

    /// Trimmed quadrature rule settings
    settings[Main::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::moment_fitting_residual, 5.6e-5);
    settings[Main::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_element_volume_ratio, 0.45);
    settings[Main::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::min_num_boundary_triangles, IndexType(234));
    settings[Main::trimmed_quadrature_rule_settings].SetValue(TrimmedQuadratureRuleSettings::use_customized_trimmed_points, true);

    QuESo_CHECK_RELATIVE_NEAR( settings[Main::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::moment_fitting_residual), 5.6e-5, 1e-10 );
    QuESo_CHECK_RELATIVE_NEAR( settings.GetSubDictionaryNoCheck(Main::trimmed_quadrature_rule_settings).GetValueNoCheck<double>(
        TrimmedQuadratureRuleSettings::moment_fitting_residual), 5.6e-5, 1e-10 );

    QuESo_CHECK_RELATIVE_NEAR( settings[Main::trimmed_quadrature_rule_settings].GetValue<double>(TrimmedQuadratureRuleSettings::min_element_volume_ratio), 0.45, 1e-10 );
    QuESo_CHECK_RELATIVE_NEAR( settings.GetSubDictionaryNoCheck(Main::trimmed_quadrature_rule_settings).GetValueNoCheck<double>(
        TrimmedQuadratureRuleSettings::min_element_volume_ratio), 0.45, 1e-10 );

    QuESo_CHECK_EQUAL( settings[Main::trimmed_quadrature_rule_settings].GetValue<IndexType>(TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 234UL );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::trimmed_quadrature_rule_settings).GetValueNoCheck<IndexType>(
        TrimmedQuadratureRuleSettings::min_num_boundary_triangles), 234UL);

    QuESo_CHECK_EQUAL( settings[Main::trimmed_quadrature_rule_settings].GetValue<bool>(TrimmedQuadratureRuleSettings::use_customized_trimmed_points), true );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::trimmed_quadrature_rule_settings).GetValueNoCheck<bool>(
        TrimmedQuadratureRuleSettings::use_customized_trimmed_points), true );

    /// Non trimmed quadrature rule settings
    settings[Main::non_trimmed_quadrature_rule_settings].SetValue(NonTrimmedQuadratureRuleSettings::integration_method, IntegrationMethod::GGQ_Optimal);

    QuESo_CHECK_EQUAL( settings[Main::non_trimmed_quadrature_rule_settings].GetValue<IntegrationMethod>(NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::GGQ_Optimal );
    QuESo_CHECK_EQUAL( settings.GetSubDictionaryNoCheck(Main::non_trimmed_quadrature_rule_settings).GetValueNoCheck<IntegrationMethod>(
        NonTrimmedQuadratureRuleSettings::integration_method), IntegrationMethod::GGQ_Optimal);
}
// BOOST_AUTO_TEST_CASE(ParameterCastWrongTypesTest) {
//     QuESo_INFO << "Testing :: Test Parameter :: Test Cast Wrong Value types" << std::endl;

//     Parameters parameters{};
//     // Check if IndexTypes are properly converted to doubles
//     parameters.Set<double>("moment_fitting_residual", 3.0);

//     parameters.Set<unsigned long>("moment_fitting_residual", 4UL);
//     BOOST_CHECK_THROW( parameters.Get<unsigned long>("moment_fitting_residual"), std::exception );
//     double value = parameters.Get<double>("moment_fitting_residual");
//     QuESo_CHECK_NEAR(value, 4.0, 1e-14);

//     parameters.Set<unsigned long>("moment_fitting_residual", 5UL);
//     BOOST_CHECK_THROW( parameters.Get<unsigned long>("moment_fitting_residual"), std::exception );
//     value = parameters.Get<double>("moment_fitting_residual");
//     QuESo_CHECK_NEAR(value, 5.0, 1e-14);

//     parameters.Set<double>("moment_fitting_residual", 5UL);
//     BOOST_CHECK_THROW( parameters.Get<unsigned long>("moment_fitting_residual"), std::exception );
//     value = parameters.Get<double>("moment_fitting_residual");
//     QuESo_CHECK_NEAR(value, 5.0, 1e-14);

//     // Check if Vector3i are properly converted to Vector3d
//     parameters.Set<Vector3d>("lower_bound_xyz", Vector3d{1.0, 1.0, 1.0});

//     parameters.Set<Vector3i>("lower_bound_xyz", Vector3i{2, 2, 2});
//     BOOST_CHECK_THROW( parameters.Get<Vector3i>("lower_bound_xyz"), std::exception );
//     auto vector = parameters.Get<Vector3d>("lower_bound_xyz");
//     QuESo_CHECK_POINT_RELATIVE_NEAR(Vector3d({2.0, 2.0, 2.0}), vector, 1e-14);

//     parameters.Set<Vector3i>("lower_bound_xyz", Vector3i{54, 5, 5});
//     BOOST_CHECK_THROW( parameters.Get<Vector3i>("lower_bound_xyz"), std::exception );
//     vector = parameters.Get<Vector3d>("lower_bound_xyz");
//     QuESo_CHECK_POINT_RELATIVE_NEAR(Vector3d({54, 5, 5}), vector, 1e-14);

//     parameters.Set<Vector3d>("lower_bound_xyz", Vector3d{3.0, 4.0, 5.0});
//     BOOST_CHECK_THROW( parameters.Get<Vector3i>("lower_bound_xyz"), std::exception );
//     vector = parameters.Get<Vector3d>("lower_bound_xyz");
//     QuESo_CHECK_POINT_RELATIVE_NEAR(Vector3d({3.0, 4.0, 5.0}), vector, 1e-14);

//     parameters.Set<Vector3d>("lower_bound_xyz", Vector3d{3.0, 4.0, 5.0});
//     BOOST_CHECK_THROW( parameters.Get<Vector3i>("lower_bound_xyz"), std::exception );
//     vector = parameters.Get<Vector3d>("lower_bound_xyz");
//     QuESo_CHECK_POINT_RELATIVE_NEAR(Vector3d({3.0, 4.0, 5.0}), vector, 1e-14);
// }


// void CheckValuesConditionParameteters(const Parameters& rParameters){

//     for( const auto& r_cond : rParameters.GetConditionsSettingsVector() ){
//         if( r_cond.Get<std::string>("type") == "PenaltySupportCondition" ){
//             QuESo_CHECK_EQUAL(r_cond.Get<std::string>("input_filename") , std::string("test.stl"));
//             QuESo_CHECK_RELATIVE_NEAR(r_cond.Get<double>("penalty_factor"), 5.0, 1e-14 );
//             QuESo_CHECK_POINT_RELATIVE_NEAR(r_cond.Get<PointType>("value"), PointType({5.0, 3.0, 4.0}), 1e-14 );
//         }
//         if( r_cond.Get<std::string>("type") == "LagrangeSupportCondition" ){
//             QuESo_CHECK_EQUAL(r_cond.Get<std::string>("input_filename") , std::string("test.stl"));
//             QuESo_CHECK_POINT_RELATIVE_NEAR(r_cond.Get<PointType>("value"), PointType({5.0, 3.0, 4.0}), 1e-14 );
//         }
//         if( r_cond.Get<std::string>("type") == "SurfaceLoadCondition" ){
//             QuESo_CHECK_EQUAL(r_cond.Get<std::string>("input_filename") , std::string("test.stl"));
//             QuESo_CHECK_RELATIVE_NEAR(r_cond.Get<double>("modulus"), 10.0, 1e-14 );
//             QuESo_CHECK_POINT_RELATIVE_NEAR(r_cond.Get<PointType>("direction"), PointType({5.0, 3.0, 4.0}), 1e-14 );
//         }
//         if( r_cond.Get<std::string>("type") == "PressureLoadCondition" ){
//             QuESo_CHECK_EQUAL(r_cond.Get<std::string>("input_filename") , std::string("test.stl"));
//             QuESo_CHECK_RELATIVE_NEAR(r_cond.Get<double>("modulus"), 15.0, 1e-14 );
//         }
//     }
// }

// BOOST_AUTO_TEST_CASE(ConditionParametetersCustomSetTest) {
//     QuESo_INFO << "Testing :: Test Parameter :: Test Condition Parameters Custom Set" << std::endl;

//     Parameters parameters{};
//     ConditionParameters penalty_support_params("PenaltySupportCondition");
//     penalty_support_params.Set("input_filename", std::string("test.stl"));
//     penalty_support_params.Set("penalty_factor", 5.0);
//     penalty_support_params.Set("value", PointType{5.0, 3.0, 4.0});
//     parameters.AddConditionSettings(penalty_support_params);

//     ConditionParameters langrange_support_params("LagrangeSupportCondition");
//     langrange_support_params.Set("input_filename", std::string("test.stl"));
//     langrange_support_params.Set("value", PointType{5.0, 3.0, 4.0});
//     parameters.AddConditionSettings(langrange_support_params);

//     ConditionParameters surface_load_params("SurfaceLoadCondition");
//     surface_load_params.Set("input_filename", std::string("test.stl"));
//     surface_load_params.Set("modulus", 10.0);
//     surface_load_params.Set("direction", PointType{5.0, 3.0, 4.0});
//     parameters.AddConditionSettings(surface_load_params);

//     ConditionParameters surface_pressure_params("PressureLoadCondition");
//     surface_pressure_params.Set("input_filename", std::string("test.stl"));
//     surface_pressure_params.Set("modulus", 15.0);
//     parameters.AddConditionSettings(surface_pressure_params);

//     CheckValuesConditionParameteters(parameters);
//     Parameters parameters_copy_const(parameters);
//     CheckValuesConditionParameteters(parameters_copy_const);
//     Parameters parameters_copy_assign = parameters;
//     CheckValuesConditionParameteters(parameters_copy_assign);
// }

// BOOST_AUTO_TEST_CASE(ConditionParametetersCustomConstructorTest) {
//     QuESo_INFO << "Testing :: Test Parameter :: Test Condition Parameters Custom Constructor" << std::endl;

//     Parameters parameters{};
//     ConditionParameters penalty_support_params({Component("type", std::string("PenaltySupportCondition")),
//                                                 Component("input_filename", std::string("test.stl")),
//                                                 Component("penalty_factor", 5.0),
//                                                 Component("value", PointType{5.0, 3.0, 4.0}) });
//     parameters.AddConditionSettings(penalty_support_params);

//     ConditionParameters langrange_support_params({Component("type", std::string("LagrangeSupportCondition")),
//                                                   Component("input_filename", std::string("test.stl")),
//                                                   Component("value", PointType{5.0, 3.0, 4.0}) });
//     parameters.AddConditionSettings(langrange_support_params);

//     ConditionParameters surface_load_params({Component("type", std::string("SurfaceLoadCondition")),
//                                              Component("input_filename", std::string("test.stl")),
//                                              Component("direction", PointType{5.0, 3.0, 4.0}),
//                                              Component("modulus", 10.0) });
//     parameters.AddConditionSettings(surface_load_params);

//     ConditionParameters surface_pressure_params({Component("type", std::string("PressureLoadCondition")),
//                                                  Component("input_filename", std::string("test.stl")),
//                                                  Component("modulus", 15.0) });
//     parameters.AddConditionSettings(surface_pressure_params);

//     CheckValuesConditionParameteters(parameters);
//     Parameters parameters_copy_const(parameters);
//     CheckValuesConditionParameteters(parameters_copy_const);
//     Parameters parameters_copy_assign = parameters;
//     CheckValuesConditionParameteters(parameters_copy_assign);
// }

// BOOST_AUTO_TEST_CASE(ConditionParametetersWrongTypesTest) {
//     QuESo_INFO << "Testing :: Test Parameter :: Test Condition Parameters Wrong Value types" << std::endl;

//     ConditionParameters penalty_support_params("PenaltySupportCondition");
//     BOOST_CHECK_THROW( penalty_support_params.Set("penalty_factor", PointType{5.0, 3.0, 4.0}), std::exception );

//     ConditionParameters penalty_support_params_2("PenaltySupportCondition");
//     BOOST_CHECK_THROW( penalty_support_params_2.Set("modulus", 5.0), std::exception );

//     ConditionParameters penalty_support_params_3("PenaltySupportCondition");
//     BOOST_CHECK_THROW( penalty_support_params_3.Set("penalty_facto", 5.0), std::exception );
// }

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

