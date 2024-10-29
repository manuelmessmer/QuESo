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
#include "queso/includes/model_info.hpp"

namespace queso {


namespace Testing {

BOOST_AUTO_TEST_SUITE( ModelInfoTestSuite )

BOOST_AUTO_TEST_CASE(ModelInfoDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test ModelInfo :: Test Default Values" << std::endl;

    {   /// Enum access
        ModelInfo model_info;

        /// embedded_geometry_info
        QuESo_CHECK( !model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::is_closed) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::embedded_geometry_info].GetValue<bool>(EmbeddedGeometryInfo::is_closed), std::exception );
        }
        QuESo_CHECK( !model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::volume) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::embedded_geometry_info].GetValue<double>(EmbeddedGeometryInfo::volume), std::exception );
        }
        /// quadrature_info
        QuESo_CHECK( !model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::represented_volume) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::represented_volume), std::exception );
        }
        QuESo_CHECK( !model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::percentage_of_geometry_volume) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::percentage_of_geometry_volume), std::exception );
        }
        QuESo_CHECK( !model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::tot_num_points) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::quadrature_info].GetValue<IndexType>(QuadratureInfo::tot_num_points), std::exception );
        }
        QuESo_CHECK( !model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_full_element) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::quadrature_info].GetValue<IndexType>(QuadratureInfo::num_of_points_per_full_element), std::exception );
        }
        QuESo_CHECK( !model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_trimmed_element) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::quadrature_info].GetValue<IndexType>(QuadratureInfo::num_of_points_per_trimmed_element), std::exception );
        }
        /// background_grid_info
        QuESo_CHECK( !model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_active_elements) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_active_elements), std::exception );
        }
        QuESo_CHECK( !model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_trimmed_elements) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_trimmed_elements), std::exception );
        }
        QuESo_CHECK( !model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_full_elements) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_full_elements), std::exception );
        }
        QuESo_CHECK( !model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_inactive_elements) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_inactive_elements), std::exception );
        }

        /// elapsed_time_info
        const auto& r_elpased_time_info = model_info[MainInfo::elapsed_time_info];
        QuESo_CHECK( !r_elpased_time_info.IsSet(ElapsedTimeInfo::total) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_elpased_time_info.GetValue<double>(ElapsedTimeInfo::total), std::exception );
        }

        QuESo_CHECK( !r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::total) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::total), std::exception );
        }
        QuESo_CHECK( !r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::classification_of_elements) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::classification_of_elements), std::exception );
        }
        QuESo_CHECK( !r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::computation_of_intersections) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::computation_of_intersections), std::exception );
        }
        QuESo_CHECK( !r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::solution_of_moment_fitting_eqs) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::solution_of_moment_fitting_eqs), std::exception );
        }
        QuESo_CHECK( !r_elpased_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::construction_of_ggq_rules) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::volume_time_info].GetValue<double>(VolumeTimeInfo::construction_of_ggq_rules), std::exception );
        }

        QuESo_CHECK( !r_elpased_time_info[ElapsedTimeInfo::conditions_time_info].IsSet(ConditionsTimeInfo::total) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::conditions_time_info].GetValue<double>(ConditionsTimeInfo::total), std::exception );
        }
        QuESo_CHECK( !r_elpased_time_info[ElapsedTimeInfo::write_files_time_info].IsSet(WriteFilesTimeInfo::total) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_elpased_time_info[ElapsedTimeInfo::write_files_time_info].GetValue<double>(WriteFilesTimeInfo::total), std::exception );
        }
    }
    {   /// String access
        ModelInfo model_info;

        /// embedded_geometry_info
        QuESo_CHECK( !model_info["embedded_geometry_info"].IsSet("is_closed") );
        BOOST_REQUIRE_THROW( model_info["embedded_geometry_info"].GetValue<bool>("is_closed"), std::exception );
        QuESo_CHECK( !model_info["embedded_geometry_info"].IsSet("volume") );
        BOOST_REQUIRE_THROW( model_info["embedded_geometry_info"].GetValue<double>("volume"), std::exception );

        /// quadrature_info
        QuESo_CHECK( !model_info["quadrature_info"].IsSet("represented_volume") );
        BOOST_REQUIRE_THROW( model_info["quadrature_info"].GetValue<double>("represented_volume"), std::exception );
        QuESo_CHECK( !model_info["quadrature_info"].IsSet("percentage_of_geometry_volume") );
        BOOST_REQUIRE_THROW( model_info["quadrature_info"].GetValue<double>("percentage_of_geometry_volume"), std::exception );
        QuESo_CHECK( !model_info["quadrature_info"].IsSet("tot_num_points") );
        BOOST_REQUIRE_THROW( model_info["quadrature_info"].GetValue<IndexType>("tot_num_points"), std::exception );
        QuESo_CHECK( !model_info["quadrature_info"].IsSet("num_of_points_per_full_element") );
        BOOST_REQUIRE_THROW( model_info["quadrature_info"].GetValue<double>("num_of_points_per_full_element"), std::exception );
        QuESo_CHECK( !model_info["quadrature_info"].IsSet("num_of_points_per_trimmed_element") );
        BOOST_REQUIRE_THROW( model_info["quadrature_info"].GetValue<double>("num_of_points_per_trimmed_element"), std::exception );

        /// background_grid_info
        QuESo_CHECK( !model_info["background_grid_info"].IsSet("num_active_elements") );
        BOOST_REQUIRE_THROW( model_info["background_grid_info"].GetValue<IndexType>("num_active_elements"), std::exception );
        QuESo_CHECK( !model_info["background_grid_info"].IsSet("num_trimmed_elements") );
        BOOST_REQUIRE_THROW( model_info["background_grid_info"].GetValue<IndexType>("num_trimmed_elements"), std::exception );
        QuESo_CHECK( !model_info["background_grid_info"].IsSet("num_full_elements") );
        BOOST_REQUIRE_THROW( model_info["background_grid_info"].GetValue<IndexType>("num_full_elements"), std::exception );
        QuESo_CHECK( !model_info["background_grid_info"].IsSet("num_inactive_elements") );
        BOOST_REQUIRE_THROW( model_info["background_grid_info"].GetValue<IndexType>("num_inactive_elements"), std::exception );


        /// elapsed_time_info
        auto& r_elpased_time_info = model_info["elapsed_time_info"];

        QuESo_CHECK( !r_elpased_time_info.IsSet("total") );
        BOOST_REQUIRE_THROW( r_elpased_time_info.GetValue<double>("total"), std::exception );
        QuESo_CHECK( !r_elpased_time_info["volume_time_info"].IsSet("total") );
        BOOST_REQUIRE_THROW( r_elpased_time_info["volume_time_info"].GetValue<double>("total"), std::exception );
        QuESo_CHECK( !r_elpased_time_info["volume_time_info"].IsSet("classification_of_elements") );
        BOOST_REQUIRE_THROW( r_elpased_time_info["volume_time_info"].GetValue<double>("classification_of_elements"), std::exception );
        QuESo_CHECK( !r_elpased_time_info["volume_time_info"].IsSet("computation_of_intersections") );
        BOOST_REQUIRE_THROW( r_elpased_time_info["volume_time_info"].GetValue<double>("computation_of_intersections"), std::exception );
        QuESo_CHECK( !r_elpased_time_info["volume_time_info"].IsSet("solution_of_moment_fitting_eqs") );
        BOOST_REQUIRE_THROW( r_elpased_time_info["volume_time_info"].GetValue<double>("solution_of_moment_fitting_eqs"), std::exception );
        QuESo_CHECK( !r_elpased_time_info["volume_time_info"].IsSet("construction_of_ggq_rules") );
        BOOST_REQUIRE_THROW( r_elpased_time_info["volume_time_info"].GetValue<double>("construction_of_ggq_rules"), std::exception );

        QuESo_CHECK( !r_elpased_time_info["conditions_time_info"].IsSet("total") );
        BOOST_REQUIRE_THROW( r_elpased_time_info["conditions_time_info"].GetValue<double>("total"), std::exception );
        QuESo_CHECK( !r_elpased_time_info["write_files_time_info"].IsSet("total") );
        BOOST_REQUIRE_THROW( r_elpased_time_info["write_files_time_info"].GetValue<double>("total"), std::exception );
    }
}


BOOST_AUTO_TEST_CASE(ModelInfoCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test ModelInfo :: Test Customized Values" << std::endl;

    {   /// Enum access
        ModelInfo model_info;

        /// General model_info
        model_info[MainInfo::embedded_geometry_info].SetValue(EmbeddedGeometryInfo::is_closed, true);
        QuESo_CHECK( model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::is_closed) );
        QuESo_CHECK_EQUAL(model_info[MainInfo::embedded_geometry_info].GetValue<bool>(EmbeddedGeometryInfo::is_closed), true);

        model_info[MainInfo::embedded_geometry_info].SetValue(EmbeddedGeometryInfo::volume, 2.0);
        QuESo_CHECK( model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::volume) );
        QuESo_CHECK_NEAR(model_info[MainInfo::embedded_geometry_info].GetValue<double>(EmbeddedGeometryInfo::volume), 2.0, EPS2);

        /// General quadrature_info
        model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::represented_volume, 2.0);
        QuESo_CHECK( model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::represented_volume) );
        QuESo_CHECK_NEAR(model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::represented_volume), 2.0, EPS2);

        model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::percentage_of_geometry_volume, 5.0);
        QuESo_CHECK( model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::percentage_of_geometry_volume) );
        QuESo_CHECK_NEAR(model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::percentage_of_geometry_volume), 5.0, EPS2);

        model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::tot_num_points, 3u);
        QuESo_CHECK( model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::tot_num_points) );
        QuESo_CHECK_EQUAL(model_info[MainInfo::quadrature_info].GetValue<IndexType>(QuadratureInfo::tot_num_points), 3);

        model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::num_of_points_per_full_element, 5.5);
        QuESo_CHECK( model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_full_element) );
        QuESo_CHECK_NEAR(model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::num_of_points_per_full_element), 5.5, EPS2);

        model_info[MainInfo::quadrature_info].SetValue(QuadratureInfo::num_of_points_per_trimmed_element, 3.3);
        QuESo_CHECK( model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_trimmed_element) );
        QuESo_CHECK_NEAR(model_info[MainInfo::quadrature_info].GetValue<double>(QuadratureInfo::num_of_points_per_trimmed_element), 3.3, EPS2);

        /// General background_grid_info
        model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_active_elements, 5);
        QuESo_CHECK( model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_active_elements) );
        QuESo_CHECK_EQUAL(model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_active_elements), 5);

        model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_trimmed_elements, 7);
        QuESo_CHECK( model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_trimmed_elements) );
        QuESo_CHECK_EQUAL(model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_trimmed_elements), 7);

        model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_full_elements, 101);
        QuESo_CHECK( model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_full_elements) );
        QuESo_CHECK_EQUAL(model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_full_elements), 101);

        model_info[MainInfo::background_grid_info].SetValue(BackgroundGridInfo::num_inactive_elements, 52);
        QuESo_CHECK( model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_inactive_elements) );
        QuESo_CHECK_EQUAL(model_info[MainInfo::background_grid_info].GetValue<IndexType>(BackgroundGridInfo::num_inactive_elements), 52);

        /// elapsed_time_info
        auto& r_elapsed_time_info = model_info[MainInfo::elapsed_time_info];
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
        ModelInfo model_info;

        /// General model_info
        model_info["embedded_geometry_info"].SetValue("is_closed", true);
        QuESo_CHECK( model_info["embedded_geometry_info"].IsSet("is_closed") );
        QuESo_CHECK_EQUAL(model_info["embedded_geometry_info"].GetValue<bool>("is_closed"), true);

        model_info["embedded_geometry_info"].SetValue("volume", 2.0);
        QuESo_CHECK( model_info["embedded_geometry_info"].IsSet("volume") );
        QuESo_CHECK_NEAR(model_info["embedded_geometry_info"].GetValue<double>("volume"), 2.0, EPS2);

        /// General quadrature_info
        model_info["quadrature_info"].SetValue("represented_volume", 2.0);
        QuESo_CHECK( model_info["quadrature_info"].IsSet("represented_volume") );
        QuESo_CHECK_NEAR(model_info["quadrature_info"].GetValue<double>("represented_volume"), 2.0, EPS2);

        model_info["quadrature_info"].SetValue("percentage_of_geometry_volume", 5.0);
        QuESo_CHECK( model_info["quadrature_info"].IsSet("percentage_of_geometry_volume") );
        QuESo_CHECK_NEAR(model_info["quadrature_info"].GetValue<double>("percentage_of_geometry_volume"), 5.0, EPS2);

        model_info["quadrature_info"].SetValue("tot_num_points", 3u);
        QuESo_CHECK( model_info["quadrature_info"].IsSet("tot_num_points") );
        QuESo_CHECK_EQUAL(model_info["quadrature_info"].GetValue<IndexType>("tot_num_points"), 3);

        model_info["quadrature_info"].SetValue("num_of_points_per_full_element", 5.5);
        QuESo_CHECK( model_info["quadrature_info"].IsSet("num_of_points_per_full_element") );
        QuESo_CHECK_NEAR(model_info["quadrature_info"].GetValue<double>("num_of_points_per_full_element"), 5.5, EPS2);

        model_info["quadrature_info"].SetValue("num_of_points_per_trimmed_element", 3.3);
        QuESo_CHECK( model_info["quadrature_info"].IsSet("num_of_points_per_trimmed_element") );
        QuESo_CHECK_NEAR(model_info["quadrature_info"].GetValue<double>("num_of_points_per_trimmed_element"), 3.3, EPS2);

        /// General background_grid_info
        model_info["background_grid_info"].SetValue("num_active_elements", 5);
        QuESo_CHECK( model_info["background_grid_info"].IsSet("num_active_elements") );
        QuESo_CHECK_EQUAL(model_info["background_grid_info"].GetValue<IndexType>("num_active_elements"), 5);

        model_info["background_grid_info"].SetValue("num_trimmed_elements", 7);
        QuESo_CHECK( model_info["background_grid_info"].IsSet("num_trimmed_elements") );
        QuESo_CHECK_EQUAL(model_info["background_grid_info"].GetValue<IndexType>("num_trimmed_elements"), 7);

        model_info["background_grid_info"].SetValue("num_full_elements", 101);
        QuESo_CHECK( model_info["background_grid_info"].IsSet("num_full_elements") );
        QuESo_CHECK_EQUAL(model_info["background_grid_info"].GetValue<IndexType>("num_full_elements"), 101);

        model_info["background_grid_info"].SetValue("num_inactive_elements", 52);
        QuESo_CHECK( model_info["background_grid_info"].IsSet("num_inactive_elements") );
        QuESo_CHECK_EQUAL(model_info["background_grid_info"].GetValue<IndexType>("num_inactive_elements"), 52);

        /// elapsed_time_info
        auto& r_elapsed_time_info = model_info["elapsed_time_info"];
        r_elapsed_time_info.SetValue("total", 101.0);
        QuESo_CHECK( r_elapsed_time_info.IsSet("total") );
        QuESo_CHECK_NEAR( r_elapsed_time_info.GetValue<double>("total"), 101.0, EPS2);

        /// volume_time_info
        r_elapsed_time_info["volume_time_info"].SetValue("total", 12.2);
        QuESo_CHECK( r_elapsed_time_info["volume_time_info"].IsSet("total") );
        QuESo_CHECK_NEAR( r_elapsed_time_info["volume_time_info"].GetValue<double>("total"), 12.2, EPS2);

        r_elapsed_time_info["volume_time_info"].SetValue("classification_of_elements", 122.2);
        QuESo_CHECK( r_elapsed_time_info["volume_time_info"].IsSet("classification_of_elements") );
        QuESo_CHECK_NEAR( r_elapsed_time_info["volume_time_info"].GetValue<double>("classification_of_elements"), 122.2, EPS2);

        r_elapsed_time_info["volume_time_info"].SetValue("computation_of_intersections", 44.5);
        QuESo_CHECK( r_elapsed_time_info["volume_time_info"].IsSet("computation_of_intersections") );
        QuESo_CHECK_NEAR( r_elapsed_time_info["volume_time_info"].GetValue<double>("computation_of_intersections"), 44.5, EPS2);

        r_elapsed_time_info["volume_time_info"].SetValue("solution_of_moment_fitting_eqs", 33.2);
        QuESo_CHECK( r_elapsed_time_info["volume_time_info"].IsSet("solution_of_moment_fitting_eqs") );
        QuESo_CHECK_NEAR( r_elapsed_time_info["volume_time_info"].GetValue<double>("solution_of_moment_fitting_eqs"), 33.2, EPS2);

        r_elapsed_time_info["volume_time_info"].SetValue("construction_of_ggq_rules", 2.23);
        QuESo_CHECK( r_elapsed_time_info["volume_time_info"].IsSet("construction_of_ggq_rules") );
        QuESo_CHECK_NEAR( r_elapsed_time_info["volume_time_info"].GetValue<double>("construction_of_ggq_rules"), 2.23, EPS2);

        /// conditions_time_info
        r_elapsed_time_info["conditions_time_info"].SetValue("total", 1.123);
        QuESo_CHECK( r_elapsed_time_info["conditions_time_info"].IsSet("total") );
        QuESo_CHECK_NEAR( r_elapsed_time_info["conditions_time_info"].GetValue<double>("total"), 1.123, EPS2);

        /// write_files_time_info
        r_elapsed_time_info["write_files_time_info"].SetValue("total", 5.123);
        QuESo_CHECK( r_elapsed_time_info["write_files_time_info"].IsSet("total") );
        QuESo_CHECK_NEAR( r_elapsed_time_info["write_files_time_info"].GetValue<double>("total"), 5.123, EPS2);
    }
}

BOOST_AUTO_TEST_CASE(ModelInfoConditionsInfoDefaultValuesTest) {
    QuESo_INFO << "Testing :: Test ModelInfo :: Test Condition Info Defaut Values" << std::endl;

    {   /// Enum access
        ModelInfo model_info;
        auto& r_cond_info = model_info.CreateNewConditionInfo();

        QuESo_CHECK( !r_cond_info.IsSet(ConditionInfo::condition_id) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_cond_info.GetValue<IndexType>(ConditionInfo::condition_id), std::exception );
        }
        QuESo_CHECK( !r_cond_info.IsSet(ConditionInfo::surf_area) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_cond_info.GetValue<double>(ConditionInfo::surf_area), std::exception );
        }
        QuESo_CHECK( !r_cond_info.IsSet(ConditionInfo::perc_surf_area_in_active_domain) );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_cond_info.GetValue<double>(ConditionInfo::perc_surf_area_in_active_domain), std::exception );
        }
    }
    {   /// String access
        ModelInfo model_info;
        auto& r_cond_info = model_info.CreateNewConditionInfo();

        QuESo_CHECK( !r_cond_info.IsSet("condition_id") );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_cond_info.GetValue<IndexType>("condition_id"), std::exception );
        }
        QuESo_CHECK( !r_cond_info.IsSet("surf_area") );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_cond_info.GetValue<double>("surf_area"), std::exception );
        }
        QuESo_CHECK( !r_cond_info.IsSet("perc_surf_area_in_active_domain") );
        if( !NOTDEBUG ) {
            BOOST_REQUIRE_THROW( r_cond_info.GetValue<double>("perc_surf_area_in_active_domain"), std::exception );
        }
    }
};

BOOST_AUTO_TEST_CASE(ModelInfoConditionsInfoCustomizedValuesTest) {
    QuESo_INFO << "Testing :: Test ModelInfo :: Test Condition Info Customized Values" << std::endl;

    {   /// Enum access
        ModelInfo model_info;

        auto& r_cond_info = model_info.CreateNewConditionInfo();

        r_cond_info.SetValue(ConditionInfo::condition_id, 101);
        QuESo_CHECK( r_cond_info.IsSet(ConditionInfo::condition_id) );
        QuESo_CHECK_EQUAL( r_cond_info.GetValue<IndexType>(ConditionInfo::condition_id), 101 );

        r_cond_info.SetValue(ConditionInfo::surf_area, 22.2);
        QuESo_CHECK( r_cond_info.IsSet(ConditionInfo::surf_area) );
        QuESo_CHECK_NEAR( r_cond_info.GetValue<double>(ConditionInfo::surf_area), 22.2, EPS2 );

        r_cond_info.SetValue(ConditionInfo::perc_surf_area_in_active_domain, 223.1);
        QuESo_CHECK( r_cond_info.IsSet(ConditionInfo::perc_surf_area_in_active_domain) );
        QuESo_CHECK_NEAR( r_cond_info.GetValue<double>(ConditionInfo::perc_surf_area_in_active_domain), 223.1, EPS2 );

        auto& r_cond_info_2 = model_info.GetList(MainInfo::conditions_infos_list)[0];
        QuESo_CHECK( r_cond_info_2.IsSet(ConditionInfo::condition_id) );
        QuESo_CHECK_EQUAL( r_cond_info_2.GetValue<IndexType>(ConditionInfo::condition_id), 101 );

        QuESo_CHECK( r_cond_info_2.IsSet(ConditionInfo::surf_area) );
        QuESo_CHECK_NEAR( r_cond_info_2.GetValue<double>(ConditionInfo::surf_area), 22.2, EPS2 );

        QuESo_CHECK( r_cond_info_2.IsSet(ConditionInfo::perc_surf_area_in_active_domain) );
        QuESo_CHECK_NEAR( r_cond_info_2.GetValue<double>(ConditionInfo::perc_surf_area_in_active_domain), 223.1, EPS2 );
    }
    {   /// String access
        ModelInfo model_info;
        auto& r_cond_info = model_info.CreateNewConditionInfo();

        r_cond_info.SetValue("condition_id", 101);
        QuESo_CHECK( r_cond_info.IsSet("condition_id") );
        QuESo_CHECK_EQUAL( r_cond_info.GetValue<IndexType>("condition_id"), 101 );

        r_cond_info.SetValue("surf_area", 22.2);
        QuESo_CHECK( r_cond_info.IsSet("surf_area") );
        QuESo_CHECK_NEAR( r_cond_info.GetValue<double>("surf_area"), 22.2, EPS2 );

        r_cond_info.SetValue("perc_surf_area_in_active_domain", 223.1);
        QuESo_CHECK( r_cond_info.IsSet("perc_surf_area_in_active_domain") );
        QuESo_CHECK_NEAR( r_cond_info.GetValue<double>("perc_surf_area_in_active_domain"), 223.1, EPS2 );

        auto& r_cond_info_2 = model_info.GetList("conditions_infos_list")[0];
        QuESo_CHECK( r_cond_info_2.IsSet("condition_id") );
        QuESo_CHECK_EQUAL( r_cond_info_2.GetValue<IndexType>("condition_id"), 101 );

        QuESo_CHECK( r_cond_info_2.IsSet("surf_area") );
        QuESo_CHECK_NEAR( r_cond_info_2.GetValue<double>("surf_area"), 22.2, EPS2 );

        QuESo_CHECK( r_cond_info_2.IsSet("perc_surf_area_in_active_domain") );
        QuESo_CHECK_NEAR( r_cond_info_2.GetValue<double>("perc_surf_area_in_active_domain"), 223.1, EPS2 );
    }
};

BOOST_AUTO_TEST_SUITE_END()

} // End namespace Testing
} // End namespace queso

