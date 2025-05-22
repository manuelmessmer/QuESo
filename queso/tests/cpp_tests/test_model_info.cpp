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

namespace queso {


namespace Testing {

BOOST_AUTO_TEST_SUITE( ModelInfoTestSuite )

BOOST_AUTO_TEST_CASE(ModelInfoDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test ModelInfo :: Test Default Values" << std::endl;

    {   /// Key access
        auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ModelInfo");
        auto& r_model_info = *p_model_info;

        /// embedded_geometry_info
        QuESo_CHECK( !r_model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::is_closed) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::embedded_geometry_info].GetValue<bool>(EmbeddedGeometryInfo::is_closed), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::embedded_geometry_info].GetValueFast<bool>(EmbeddedGeometryInfo::is_closed), std::exception );
        }
        QuESo_CHECK( !r_model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::volume) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::embedded_geometry_info].GetValue<double>(EmbeddedGeometryInfo::volume), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::embedded_geometry_info].GetValueFast<double>(EmbeddedGeometryInfo::volume), std::exception );
        }
        /// quadrature_info
        QuESo_CHECK( !r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::represented_volume) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::represented_volume), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValueFast<double>(QuadratureInfo::represented_volume), std::exception );
        }
        QuESo_CHECK( !r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::percentage_of_geometry_volume) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::percentage_of_geometry_volume), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValueFast<double>(QuadratureInfo::percentage_of_geometry_volume), std::exception );
        }
        QuESo_CHECK( !r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::tot_num_points) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValue<IndexType>(QuadratureInfo::tot_num_points), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValueFast<IndexType>(QuadratureInfo::tot_num_points), std::exception );
        }
        QuESo_CHECK( !r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_full_element) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::num_of_points_per_full_element), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValueFast<double>(QuadratureInfo::num_of_points_per_full_element), std::exception );
        }
        QuESo_CHECK( !r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_trimmed_element) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::num_of_points_per_trimmed_element), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::quadrature_info].GetValueFast<double>(QuadratureInfo::num_of_points_per_trimmed_element), std::exception );
        }
        /// background_grid_info
        QuESo_CHECK( !r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_active_elements) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_active_elements), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::background_grid_info].GetValueFast<IndexType>(BackgroundGridInfo::num_active_elements), std::exception );
        }
        QuESo_CHECK( !r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_trimmed_elements) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_trimmed_elements), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::background_grid_info].GetValueFast<IndexType>(BackgroundGridInfo::num_trimmed_elements), std::exception );
        }
        QuESo_CHECK( !r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_full_elements) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_full_elements), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::background_grid_info].GetValueFast<IndexType>(BackgroundGridInfo::num_full_elements), std::exception );
        }
        QuESo_CHECK( !r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_inactive_elements) );
        BOOST_REQUIRE_THROW( r_model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_inactive_elements), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_model_info[MainInfo::background_grid_info].GetValueFast<IndexType>(BackgroundGridInfo::num_inactive_elements), std::exception );
        }

        /// elapsed_time_info
        const auto& r_elpased_time_info = r_model_info[MainInfo::elapsed_time_info];
        QuESo_CHECK( r_elpased_time_info.IsSet(ElapsedTimeInfo::total) );
        QuESo_CHECK_NEAR(r_elpased_time_info.GetValue<double>(ElapsedTimeInfo::total), 0.0, EPS0);
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_elpased_time_info.GetValue<IndexType>(ConditionInfo::condition_id), std::exception ); // Wrong type
            BOOST_REQUIRE_THROW( r_elpased_time_info.GetValue<double>(ConditionInfo::surf_area), std::exception ); // Wrong key
            BOOST_REQUIRE_THROW( r_elpased_time_info.GetValueFast<IndexType>(ConditionInfo::condition_id), std::exception ); // Wrong type
            BOOST_REQUIRE_THROW( r_elpased_time_info.GetValueFast<double>(ConditionInfo::surf_area), std::exception ); // Wrong key
        }

        QuESo_CHECK( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::total) );
        QuESo_CHECK_NEAR(r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::total), 0.0, EPS0);
        if( !NOTDEBUG ) { // Wrong key type
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<IndexType>(ConditionInfo::condition_id), std::exception );
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValueFast<IndexType>(ConditionInfo::condition_id), std::exception );
        }
        QuESo_CHECK( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::classification_of_elements) );
        QuESo_CHECK_NEAR(r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::classification_of_elements), 0.0, EPS0);
        if( !NOTDEBUG ) { // Wrong key type
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(ConditionsTimeInfo::total), std::exception );
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValueFast<double>(ConditionsTimeInfo::total), std::exception );
        }
        QuESo_CHECK( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::computation_of_intersections) );
        QuESo_CHECK_NEAR(r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::computation_of_intersections), 0.0, EPS0);
        if( !NOTDEBUG ) { // Wrong type
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(ConditionsTimeInfo::total), std::exception );
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValueFast<double>(ConditionsTimeInfo::total), std::exception );
        }
        QuESo_CHECK( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::solution_of_moment_fitting_eqs) );
        QuESo_CHECK_NEAR(r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::solution_of_moment_fitting_eqs), 0.0, EPS0);

        QuESo_CHECK( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::construction_of_ggq_rules) );
        QuESo_CHECK_NEAR(r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::construction_of_ggq_rules), 0.0, EPS0);

        QuESo_CHECK( r_elpased_time_info[ElapsedTimeInfo::conditions_time_info].IsSet(ConditionsTimeInfo::total) );
        QuESo_CHECK_NEAR(r_elpased_time_info[ElapsedTimeInfo::conditions_time_info].GetValue<double>(ConditionsTimeInfo::total), 0.0, EPS0);

        QuESo_CHECK( r_elpased_time_info[ElapsedTimeInfo::write_files_time_info].IsSet(WriteFilesTimeInfo::total) );
        QuESo_CHECK_NEAR(r_elpased_time_info[ElapsedTimeInfo::write_files_time_info].GetValue<double>(WriteFilesTimeInfo::total), 0.0, EPS0);
    }
    {   /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ModelInfo");
        auto& r_model_info = *p_model_info;

        /// embedded_geometry_info
        auto& r_embedded_geometry_info = StringAccess::GetSubDictionary(r_model_info, "embedded_geometry_info");
        QuESo_CHECK( !StringAccess::IsSet(r_embedded_geometry_info, "is_closed") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<bool>(r_embedded_geometry_info, "is_closed"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_embedded_geometry_info, "volume") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(r_embedded_geometry_info, "volume"), std::exception );

        /// quadrature_info
        auto& r_quadrature_info = StringAccess::GetSubDictionary(r_model_info, "quadrature_info");
        QuESo_CHECK( !StringAccess::IsSet(r_quadrature_info, "represented_volume") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(r_quadrature_info, "represented_volume"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_quadrature_info, "percentage_of_geometry_volume") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(r_quadrature_info, "percentage_of_geometry_volume"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_quadrature_info, "tot_num_points") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_quadrature_info, "tot_num_points"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_quadrature_info, "num_of_points_per_full_element") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(r_quadrature_info, "num_of_points_per_full_element"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_quadrature_info, "num_of_points_per_trimmed_element") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(r_quadrature_info, "num_of_points_per_trimmed_element"), std::exception );

        /// background_grid_info
        auto& r_background_grid_info = StringAccess::GetSubDictionary(r_model_info, "background_grid_info");
        QuESo_CHECK( !StringAccess::IsSet(r_background_grid_info, "num_active_elements") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_background_grid_info, "num_active_elements"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_background_grid_info, "num_trimmed_elements") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_background_grid_info, "num_trimmed_elements"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_background_grid_info, "num_full_elements") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_background_grid_info, "num_full_elements"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_background_grid_info, "num_inactive_elements") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_background_grid_info, "num_inactive_elements"), std::exception );


        /// elapsed_time_info
        auto& r_elpased_time_info = StringAccess::GetSubDictionary(r_model_info, "elapsed_time_info");

        QuESo_CHECK( StringAccess::IsSet(r_elpased_time_info, "total") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_elpased_time_info, "total"), 0.0, EPS0);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_elpased_time_info, "total"), std::exception ); // Wrong type
        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(r_elpased_time_info, "total2"), std::exception ); // Wrong key

        /// volume_time_info
        auto& r_volume_time_info = StringAccess::GetSubDictionary(r_elpased_time_info, "volume_time_info");
        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "total") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_volume_time_info, "total"), 0.0, EPS0);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_volume_time_info, "total"), std::exception ); // Wrong type

        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "classification_of_elements") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_volume_time_info, "classification_of_elements"), 0.0, EPS0);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_volume_time_info, "classification_of_elements"), std::exception ); // Wrong type

        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "computation_of_intersections") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_volume_time_info, "computation_of_intersections"), 0.0, EPS0);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_volume_time_info, "computation_of_intersections"), std::exception ); // Wrong type

        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "solution_of_moment_fitting_eqs") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_volume_time_info, "solution_of_moment_fitting_eqs"), 0.0, EPS0);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_volume_time_info, "solution_of_moment_fitting_eqs"), std::exception );  // Wrong type

        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "construction_of_ggq_rules") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_volume_time_info, "construction_of_ggq_rules"), 0.0, EPS0);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_volume_time_info, "construction_of_ggq_rules"), std::exception );  // Wrong type

        /// conditions_time_info
        auto& r_conditions_time_info = StringAccess::GetSubDictionary(r_elpased_time_info, "conditions_time_info");
        QuESo_CHECK( StringAccess::IsSet(r_conditions_time_info, "total") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_conditions_time_info, "total"), 0.0, EPS0);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_conditions_time_info, "total"), std::exception );  // Wrong type

        /// write_files_time_info
        auto& r_write_files_time_info = StringAccess::GetSubDictionary(r_elpased_time_info, "write_files_time_info");
        QuESo_CHECK( StringAccess::IsSet(r_write_files_time_info, "total") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_write_files_time_info, "total"), 0.0, EPS0);
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_write_files_time_info, "total"), std::exception ); // Wrong type
    }
}


BOOST_AUTO_TEST_CASE(ModelInfoCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test ModelInfo :: Test Customized Values" << std::endl;

    {   /// Key access
        auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ModelInfo");
        auto& r_model_info = *p_model_info;

        /// General r_model_info
        r_model_info[MainInfo::embedded_geometry_info].SetValue(EmbeddedGeometryInfo::is_closed, true);
        QuESo_CHECK( r_model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::is_closed) );
        QuESo_CHECK_EQUAL(r_model_info[MainInfo::embedded_geometry_info].GetValue<bool>(EmbeddedGeometryInfo::is_closed), true);

        r_model_info[MainInfo::embedded_geometry_info].SetValue(EmbeddedGeometryInfo::volume, 2.0);
        QuESo_CHECK( r_model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::volume) );
        QuESo_CHECK_NEAR(r_model_info[MainInfo::embedded_geometry_info].GetValue<double>(EmbeddedGeometryInfo::volume), 2.0, EPS2);

        /// General quadrature_info
        r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::represented_volume, 2.0);
        QuESo_CHECK( r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::represented_volume) );
        QuESo_CHECK_NEAR(r_model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::represented_volume), 2.0, EPS2);

        r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::percentage_of_geometry_volume, 5.0);
        QuESo_CHECK( r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::percentage_of_geometry_volume) );
        QuESo_CHECK_NEAR(r_model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::percentage_of_geometry_volume), 5.0, EPS2);

        r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::tot_num_points, 3u);
        QuESo_CHECK( r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::tot_num_points) );
        QuESo_CHECK_EQUAL(r_model_info[MainInfo::quadrature_info].GetValue<IndexType>(QuadratureInfo::tot_num_points), 3);

        r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::num_of_points_per_full_element, 5.5);
        QuESo_CHECK( r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_full_element) );
        QuESo_CHECK_NEAR(r_model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::num_of_points_per_full_element), 5.5, EPS2);

        r_model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::num_of_points_per_trimmed_element, 3.3);
        QuESo_CHECK( r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_trimmed_element) );
        QuESo_CHECK_NEAR(r_model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::num_of_points_per_trimmed_element), 3.3, EPS2);

        /// General background_grid_info
        r_model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_active_elements, 5u);
        QuESo_CHECK( r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_active_elements) );
        QuESo_CHECK_EQUAL(r_model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_active_elements), 5);

        r_model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_trimmed_elements, 7u);
        QuESo_CHECK( r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_trimmed_elements) );
        QuESo_CHECK_EQUAL(r_model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_trimmed_elements), 7);

        r_model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_full_elements, 101u);
        QuESo_CHECK( r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_full_elements) );
        QuESo_CHECK_EQUAL(r_model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_full_elements), 101);

        r_model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_inactive_elements, 52u);
        QuESo_CHECK( r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_inactive_elements) );
        QuESo_CHECK_EQUAL(r_model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_inactive_elements), 52);

        /// elapsed_time_info
        auto& r_elapsed_time_info = r_model_info[MainInfo::elapsed_time_info];
        r_elapsed_time_info.SetValue(ElapsedTimeInfo::total, 101.0);
        QuESo_CHECK( r_elapsed_time_info.IsSet(ElapsedTimeInfo::total) );
        QuESo_CHECK_NEAR( r_elapsed_time_info.GetValue<double>(ElapsedTimeInfo::total), 101.0, EPS2);

        /// volume_time_info
        r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(VolumeTimeInfo::total, 12.2);
        QuESo_CHECK( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::total) );
        QuESo_CHECK_NEAR( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::total), 12.2, EPS2);

        r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(VolumeTimeInfo::classification_of_elements, 122.2);
        QuESo_CHECK( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::classification_of_elements) );
        QuESo_CHECK_NEAR( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::classification_of_elements), 122.2, EPS2);

        r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(VolumeTimeInfo::computation_of_intersections, 44.5);
        QuESo_CHECK( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::computation_of_intersections) );
        QuESo_CHECK_NEAR( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::computation_of_intersections), 44.5, EPS2);

        r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(VolumeTimeInfo::solution_of_moment_fitting_eqs, 33.2);
        QuESo_CHECK( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::solution_of_moment_fitting_eqs) );
        QuESo_CHECK_NEAR( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::solution_of_moment_fitting_eqs), 33.2, EPS2);

        r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(VolumeTimeInfo::construction_of_ggq_rules, 2.23);
        QuESo_CHECK( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::construction_of_ggq_rules) );
        QuESo_CHECK_NEAR( r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::construction_of_ggq_rules), 2.23, EPS2);

        /// conditions_time_info
        r_elapsed_time_info[ElapsedTimeInfo::conditions_time_info].SetValue(ConditionsTimeInfo::total, 1.123);
        QuESo_CHECK( r_elapsed_time_info[ElapsedTimeInfo::conditions_time_info].IsSet(ConditionsTimeInfo::total) );
        QuESo_CHECK_NEAR( r_elapsed_time_info[ElapsedTimeInfo::conditions_time_info].GetValue<double>(ConditionsTimeInfo::total), 1.123, EPS2);

        /// write_files_time_info
        r_elapsed_time_info[ElapsedTimeInfo::write_files_time_info].SetValue(WriteFilesTimeInfo::total, 5.123);
        QuESo_CHECK( r_elapsed_time_info[ElapsedTimeInfo::write_files_time_info].IsSet(WriteFilesTimeInfo::total) );
        QuESo_CHECK_NEAR( r_elapsed_time_info[ElapsedTimeInfo::write_files_time_info].GetValue<double>(WriteFilesTimeInfo::total), 5.123, EPS2);
    }

    {   /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ModelInfo");
        auto& r_model_info = *p_model_info;

        /// General embedded geometry info
        auto& r_embedded_geometry_info = StringAccess::GetSubDictionary(r_model_info, "embedded_geometry_info");
        StringAccess::SetValue(r_embedded_geometry_info, "is_closed", true);
        QuESo_CHECK( StringAccess::IsSet(r_embedded_geometry_info, "is_closed") );
        QuESo_CHECK_EQUAL(StringAccess::GetValue<bool>(r_embedded_geometry_info, "is_closed"), true);

        StringAccess::SetValue(r_embedded_geometry_info, "volume", 2.0);
        QuESo_CHECK( StringAccess::IsSet(r_embedded_geometry_info, "volume") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_embedded_geometry_info, "volume"), 2.0, EPS2);

        /// General quadrature_info
        auto& r_quadrature_info = StringAccess::GetSubDictionary(r_model_info, "quadrature_info");
        StringAccess::SetValue(r_quadrature_info, "represented_volume", 2.0);
        QuESo_CHECK( StringAccess::IsSet(r_quadrature_info, "represented_volume") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_quadrature_info, "represented_volume"), 2.0, EPS2);

        StringAccess::SetValue(r_quadrature_info, "percentage_of_geometry_volume", 5.0);
        QuESo_CHECK( StringAccess::IsSet(r_quadrature_info, "percentage_of_geometry_volume") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_quadrature_info, "percentage_of_geometry_volume"), 5.0, EPS2);

        StringAccess::SetValue(r_quadrature_info, "tot_num_points", 3u);
        QuESo_CHECK( StringAccess::IsSet(r_quadrature_info, "tot_num_points") );
        QuESo_CHECK_EQUAL(StringAccess::GetValue<IndexType>(r_quadrature_info, "tot_num_points"), 3);

        StringAccess::SetValue(r_quadrature_info, "num_of_points_per_full_element", 5.5);
        QuESo_CHECK( StringAccess::IsSet(r_quadrature_info, "num_of_points_per_full_element") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_quadrature_info, "num_of_points_per_full_element"), 5.5, EPS2);

        StringAccess::SetValue(r_quadrature_info, "num_of_points_per_trimmed_element", 3.3);
        QuESo_CHECK( StringAccess::IsSet(r_quadrature_info, "num_of_points_per_trimmed_element") );
        QuESo_CHECK_NEAR(StringAccess::GetValue<double>(r_quadrature_info, "num_of_points_per_trimmed_element"), 3.3, EPS2);

        /// General background_grid_info
        auto& r_background_grid_info = StringAccess::GetSubDictionary(r_model_info, "background_grid_info");
        StringAccess::SetValue(r_background_grid_info, "num_active_elements", 5u);
        QuESo_CHECK( StringAccess::IsSet(r_background_grid_info, "num_active_elements") );
        QuESo_CHECK_EQUAL(StringAccess::GetValue<IndexType>(r_background_grid_info, "num_active_elements"), 5);

        StringAccess::SetValue(r_background_grid_info, "num_trimmed_elements", 7u);
        QuESo_CHECK( StringAccess::IsSet(r_background_grid_info, "num_trimmed_elements") );
        QuESo_CHECK_EQUAL(StringAccess::GetValue<IndexType>(r_background_grid_info, "num_trimmed_elements"), 7);

        StringAccess::SetValue(r_background_grid_info, "num_full_elements", 101u);
        QuESo_CHECK( StringAccess::IsSet(r_background_grid_info, "num_full_elements") );
        QuESo_CHECK_EQUAL(StringAccess::GetValue<IndexType>(r_background_grid_info, "num_full_elements"), 101);

        StringAccess::SetValue(r_background_grid_info, "num_inactive_elements", 52u);
        QuESo_CHECK( StringAccess::IsSet(r_background_grid_info, "num_inactive_elements") );
        QuESo_CHECK_EQUAL(StringAccess::GetValue<IndexType>(r_background_grid_info, "num_inactive_elements"), 52);

        /// elapsed_time_info
        auto& r_elapsed_time_info = StringAccess::GetSubDictionary(r_model_info, "elapsed_time_info");
        StringAccess::SetValue(r_elapsed_time_info, "total", 101.0);
        QuESo_CHECK( StringAccess::IsSet(r_elapsed_time_info, "total") );
        QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_elapsed_time_info, "total"), 101.0, EPS2);

        /// volume_time_info
        auto& r_volume_time_info = StringAccess::GetSubDictionary(r_elapsed_time_info, "volume_time_info");
        StringAccess::SetValue(r_volume_time_info, "total", 12.2);
        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "total") );
        QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_volume_time_info, "total"), 12.2, EPS2);

        StringAccess::SetValue(r_volume_time_info, "classification_of_elements", 122.2);
        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "classification_of_elements") );
        QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_volume_time_info, "classification_of_elements"), 122.2, EPS2);

        StringAccess::SetValue(r_volume_time_info, "computation_of_intersections", 44.5);
        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "computation_of_intersections") );
        QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_volume_time_info, "computation_of_intersections"), 44.5, EPS2);

        StringAccess::SetValue(r_volume_time_info, "solution_of_moment_fitting_eqs", 33.2);
        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "solution_of_moment_fitting_eqs") );
        QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_volume_time_info, "solution_of_moment_fitting_eqs"), 33.2, EPS2);

        StringAccess::SetValue(r_volume_time_info, "construction_of_ggq_rules", 2.23);
        QuESo_CHECK( StringAccess::IsSet(r_volume_time_info, "construction_of_ggq_rules") );
        QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_volume_time_info, "construction_of_ggq_rules"), 2.23, EPS2);

        /// conditions_time_info
        auto& r_conditions_time_info = StringAccess::GetSubDictionary(r_elapsed_time_info, "conditions_time_info");
        StringAccess::SetValue(r_conditions_time_info, "total", 1.123);
        QuESo_CHECK( StringAccess::IsSet(r_conditions_time_info, "total") );
        QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_conditions_time_info, "total"), 1.123, EPS2);

        /// write_files_time_info
        auto& r_write_files_time_info = StringAccess::GetSubDictionary(r_elapsed_time_info, "write_files_time_info");
        StringAccess::SetValue(r_write_files_time_info, "total", 5.123);
        QuESo_CHECK( StringAccess::IsSet(r_write_files_time_info, "total") );
        QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_write_files_time_info, "total"), 5.123, EPS2);
    }
}

BOOST_AUTO_TEST_CASE(ModelInfoConditionsInfoDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test ModelInfo :: Test Condition Info Defaut Values" << std::endl;

    {   /// Key access
        auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
        auto& r_cond_info = *p_cond_info;

        QuESo_CHECK( !r_cond_info.IsSet(ConditionInfo::condition_id) );
        BOOST_REQUIRE_THROW( r_cond_info.GetValue<IndexType>(ConditionInfo::condition_id), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_cond_info.GetValueFast<IndexType>(ConditionInfo::condition_id), std::exception );
        }
        QuESo_CHECK( !r_cond_info.IsSet(ConditionInfo::surf_area) );
        BOOST_REQUIRE_THROW( r_cond_info.GetValue<double>(ConditionInfo::surf_area), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_cond_info.GetValueFast<double>(ConditionInfo::surf_area), std::exception );
        }
        QuESo_CHECK( !r_cond_info.IsSet(ConditionInfo::perc_surf_area_in_active_domain) );
        BOOST_REQUIRE_THROW( r_cond_info.GetValue<double>(ConditionInfo::perc_surf_area_in_active_domain), std::exception );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_cond_info.GetValueFast<double>(ConditionInfo::perc_surf_area_in_active_domain), std::exception );
        }
    }
    {   /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
        auto& r_cond_info = *p_cond_info;

        QuESo_CHECK( !StringAccess::IsSet(r_cond_info, "condition_id") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<IndexType>(r_cond_info, "condition_id"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_cond_info, "surf_area") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(r_cond_info, "surf_area"), std::exception );
        QuESo_CHECK( !StringAccess::IsSet(r_cond_info, "perc_surf_area_in_active_domain") );
        BOOST_REQUIRE_THROW( StringAccess::GetValue<double>(r_cond_info, "perc_surf_area_in_active_domain"), std::exception );
    }
};

BOOST_AUTO_TEST_CASE(ModelInfoConditionsInfoCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test ModelInfo :: Test Condition Info Customized Values" << std::endl;

    {   /// Key access
        auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ModelInfo");

        {
            auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
            auto& r_cond_info = *p_cond_info;

            r_cond_info.SetValue(ConditionInfo::condition_id, 101u);
            QuESo_CHECK( r_cond_info.IsSet(ConditionInfo::condition_id) );
            QuESo_CHECK_EQUAL( r_cond_info.GetValue<IndexType>(ConditionInfo::condition_id), 101);

            r_cond_info.SetValue(ConditionInfo::surf_area, 22.2);
            QuESo_CHECK( r_cond_info.IsSet(ConditionInfo::surf_area) );
            QuESo_CHECK_NEAR( r_cond_info.GetValue<double>(ConditionInfo::surf_area), 22.2, EPS2 );

            r_cond_info.SetValue(ConditionInfo::perc_surf_area_in_active_domain, 223.1);
            QuESo_CHECK( r_cond_info.IsSet(ConditionInfo::perc_surf_area_in_active_domain) );
            QuESo_CHECK_NEAR( r_cond_info.GetValue<double>(ConditionInfo::perc_surf_area_in_active_domain), 223.1, EPS2 );

            auto& r_cond_list = p_model_info->GetList(MainInfo::conditions_infos_list);
            r_cond_list.push_back( std::move(p_cond_info) );
        }
        {
            auto& p_cond_info = p_model_info->GetList(MainInfo::conditions_infos_list)[0];
            auto& r_cond_info = *p_cond_info;

            QuESo_CHECK( r_cond_info.IsSet(ConditionInfo::condition_id) );
            QuESo_CHECK_EQUAL( r_cond_info.GetValue<IndexType>(ConditionInfo::condition_id), 101 );

            QuESo_CHECK( r_cond_info.IsSet(ConditionInfo::surf_area) );
            QuESo_CHECK_NEAR( r_cond_info.GetValue<double>(ConditionInfo::surf_area), 22.2, EPS2 );

            QuESo_CHECK( r_cond_info.IsSet(ConditionInfo::perc_surf_area_in_active_domain) );
            QuESo_CHECK_NEAR( r_cond_info.GetValue<double>(ConditionInfo::perc_surf_area_in_active_domain), 223.1, EPS2 );
        }
    }
    {   /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ModelInfo");
        {
            auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
            auto& r_cond_info = *p_cond_info;

            StringAccess::SetValue(r_cond_info, "condition_id", 101);
            QuESo_CHECK( StringAccess::IsSet(r_cond_info, "condition_id") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(r_cond_info, "condition_id"), 101 );

            StringAccess::SetValue(r_cond_info, "surf_area", 22.2);
            QuESo_CHECK( StringAccess::IsSet(r_cond_info, "surf_area") );
            QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_cond_info, "surf_area"), 22.2, EPS2 );

            StringAccess::SetValue(r_cond_info, "perc_surf_area_in_active_domain", 223.1);
            QuESo_CHECK( StringAccess::IsSet(r_cond_info, "perc_surf_area_in_active_domain") );
            QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_cond_info, "perc_surf_area_in_active_domain"), 223.1, EPS2 );

            auto& r_cond_list = StringAccess::GetList(*p_model_info, "conditions_infos_list");
            r_cond_list.push_back( std::move(p_cond_info) );
        }
        {
            const auto& p_cond_list = StringAccess::GetList(*p_model_info, "conditions_infos_list")[0];
            const auto& r_cond_info = *p_cond_list;

            QuESo_CHECK( StringAccess::IsSet(r_cond_info, "condition_id") );
            QuESo_CHECK_EQUAL( StringAccess::GetValue<IndexType>(r_cond_info, "condition_id"), 101 );

            QuESo_CHECK( StringAccess::IsSet(r_cond_info, "surf_area") );
            QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_cond_info, "surf_area"), 22.2, EPS2 );

            QuESo_CHECK( StringAccess::IsSet(r_cond_info, "perc_surf_area_in_active_domain") );
            QuESo_CHECK_NEAR( StringAccess::GetValue<double>(r_cond_info, "perc_surf_area_in_active_domain"), 223.1, EPS2 );
        }
    }
};

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

