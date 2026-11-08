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
#include <doctest/doctest.h>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/includes/dictionary_factory.hpp"

namespace queso::testing {

    TEST_SUITE("ModelInfoTestSuite")
    {

        TEST_CASE("ModelInfoTestSuite - Model Info Default Values Test")
        {


            {  /// Key access
                auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ComponentInfo");
                auto& r_model_info = *p_model_info;

                REQUIRE_THROWS_AS(r_model_info.CheckRequired(), std::exception);

                /// embedded_geometry_info
                CHECK(!r_model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::is_closed));
                CHECK(!r_model_info[MainInfo::embedded_geometry_info].IsSet(EmbeddedGeometryInfo::volume));
                /// quadrature_info
                CHECK(!r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::represented_volume));
                CHECK(!r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::percentage_of_geometry_volume));
                CHECK(!r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::tot_num_points));
                CHECK(!r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_full_element));
                CHECK(
                    !r_model_info[MainInfo::quadrature_info].IsSet(QuadratureInfo::num_of_points_per_trimmed_element)
                );
                /// background_grid_info
                CHECK(!r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_active_elements));
                CHECK(!r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_trimmed_elements));
                CHECK(!r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_full_elements));
                CHECK(!r_model_info[MainInfo::background_grid_info].IsSet(BackgroundGridInfo::num_inactive_elements));
                if (!NOTDEBUG) {
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::embedded_geometry_info].GetRequiredValue<bool>(
                            EmbeddedGeometryInfo::is_closed
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::embedded_geometry_info].GetRequiredValue<double>(
                            EmbeddedGeometryInfo::volume
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::quadrature_info].GetRequiredValue<double>(
                            QuadratureInfo::represented_volume
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::quadrature_info].GetRequiredValue<double>(
                            QuadratureInfo::percentage_of_geometry_volume
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::quadrature_info].GetRequiredValue<IndexType>(
                            QuadratureInfo::tot_num_points
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::quadrature_info].GetRequiredValue<double>(
                            QuadratureInfo::num_of_points_per_full_element
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::quadrature_info].GetRequiredValue<double>(
                            QuadratureInfo::num_of_points_per_trimmed_element
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::background_grid_info].GetRequiredValue<IndexType>(
                            BackgroundGridInfo::num_active_elements
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::background_grid_info].GetRequiredValue<IndexType>(
                            BackgroundGridInfo::num_trimmed_elements
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::background_grid_info].GetRequiredValue<IndexType>(
                            BackgroundGridInfo::num_full_elements
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        r_model_info[MainInfo::background_grid_info].GetRequiredValue<IndexType>(
                            BackgroundGridInfo::num_inactive_elements
                        ),
                        std::exception
                    );
                }

                /// elapsed_time_info
                const auto& r_elpased_time_info = r_model_info[MainInfo::elapsed_time_info];
                r_elpased_time_info.CheckRequired();
                CHECK(r_elpased_time_info.IsSet(ElapsedTimeInfo::total));
                CheckNear(r_elpased_time_info.GetRequiredValue<double>(ElapsedTimeInfo::total), 0.0, EPS0);
                if (!NOTDEBUG) {
                    REQUIRE_THROWS_AS(
                        r_elpased_time_info.GetRequiredValue<IndexType>(ConditionInfo::condition_id), std::exception
                    );  // Wrong key
                    REQUIRE_THROWS_AS(
                        r_elpased_time_info.GetRequiredValue<double>(ConditionInfo::surf_area), std::exception
                    );  // Wrong key
                }

                const auto& r_volume_time_info = r_elpased_time_info[ElapsedTimeInfo::volume_time_info];
                r_volume_time_info.CheckRequired();
                CHECK(r_volume_time_info.IsSet(VolumeTimeInfo::total));
                CheckNear(r_volume_time_info.GetRequiredValue<double>(VolumeTimeInfo::total), 0.0, EPS0);
                if (!NOTDEBUG) {  // Wrong key type
                    REQUIRE_THROWS_AS(
                        r_volume_time_info.GetRequiredValue<IndexType>(ConditionInfo::condition_id), std::exception
                    );
                }
                CHECK(r_volume_time_info.IsSet(VolumeTimeInfo::classification_of_elements));
                CheckNear(
                    r_volume_time_info.GetRequiredValue<double>(VolumeTimeInfo::classification_of_elements), 0.0, EPS0
                );
                if (!NOTDEBUG) {  // Wrong key type
                    REQUIRE_THROWS_AS(
                        r_volume_time_info.GetRequiredValue<double>(ConditionsTimeInfo::total), std::exception
                    );
                }
                CHECK(r_volume_time_info.IsSet(VolumeTimeInfo::computation_of_intersections));
                CheckNear(
                    r_volume_time_info.GetRequiredValue<double>(VolumeTimeInfo::computation_of_intersections), 0.0, EPS0
                );
                CHECK(r_volume_time_info.IsSet(VolumeTimeInfo::solution_of_moment_fitting_eqs));
                CheckNear(
                    r_volume_time_info.GetRequiredValue<double>(VolumeTimeInfo::solution_of_moment_fitting_eqs),
                    0.0,
                    EPS0
                );

                CHECK(r_volume_time_info.IsSet(VolumeTimeInfo::construction_of_ggq_rules));
                CheckNear(
                    r_volume_time_info.GetRequiredValue<double>(VolumeTimeInfo::construction_of_ggq_rules), 0.0, EPS0
                );

                const auto& r_conditions_time_info = r_elpased_time_info[ElapsedTimeInfo::conditions_time_info];
                r_conditions_time_info.CheckRequired();
                CHECK(r_conditions_time_info.IsSet(ConditionsTimeInfo::total));
                CheckNear(r_conditions_time_info.GetRequiredValue<double>(ConditionsTimeInfo::total), 0.0, EPS0);

                const auto& r_write_files_time_info = r_elpased_time_info[ElapsedTimeInfo::write_files_time_info];
                r_write_files_time_info.CheckRequired();
                CHECK(r_write_files_time_info.IsSet(WriteFilesTimeInfo::total));
                CheckNear(r_write_files_time_info.GetRequiredValue<double>(WriteFilesTimeInfo::total), 0.0, EPS0);
            }
            {  /// String access
                using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
                using StringAccess = DictionaryStringAccess<DictionaryType>;

                auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ComponentInfo");
                auto& r_model_info = *p_model_info;

                REQUIRE_THROWS_AS(r_model_info.CheckRequired(), std::exception);

                /// embedded_geometry_info
                auto& r_embedded_geometry_info = StringAccess::GetSubDictionary(r_model_info, "embedded_geometry_info");
                CHECK(!StringAccess::IsSet(r_embedded_geometry_info, "is_closed"));
                REQUIRE_THROWS_AS(StringAccess::GetValue<bool>(r_embedded_geometry_info, "is_closed"), std::exception);
                CHECK(!StringAccess::IsSet(r_embedded_geometry_info, "volume"));
                REQUIRE_THROWS_AS(StringAccess::GetValue<double>(r_embedded_geometry_info, "volume"), std::exception);

                /// quadrature_info
                auto& r_quadrature_info = StringAccess::GetSubDictionary(r_model_info, "quadrature_info");
                CHECK(!StringAccess::IsSet(r_quadrature_info, "represented_volume"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<double>(r_quadrature_info, "represented_volume"), std::exception
                );
                CHECK(!StringAccess::IsSet(r_quadrature_info, "percentage_of_geometry_volume"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<double>(r_quadrature_info, "percentage_of_geometry_volume"), std::exception
                );
                CHECK(!StringAccess::IsSet(r_quadrature_info, "tot_num_points"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_quadrature_info, "tot_num_points"), std::exception
                );
                CHECK(!StringAccess::IsSet(r_quadrature_info, "num_of_points_per_full_element"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<double>(r_quadrature_info, "num_of_points_per_full_element"), std::exception
                );
                CHECK(!StringAccess::IsSet(r_quadrature_info, "num_of_points_per_trimmed_element"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<double>(r_quadrature_info, "num_of_points_per_trimmed_element"),
                    std::exception
                );

                /// background_grid_info
                auto& r_background_grid_info = StringAccess::GetSubDictionary(r_model_info, "background_grid_info");
                CHECK(!StringAccess::IsSet(r_background_grid_info, "num_active_elements"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_background_grid_info, "num_active_elements"), std::exception
                );
                CHECK(!StringAccess::IsSet(r_background_grid_info, "num_trimmed_elements"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_background_grid_info, "num_trimmed_elements"), std::exception
                );
                CHECK(!StringAccess::IsSet(r_background_grid_info, "num_full_elements"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_background_grid_info, "num_full_elements"), std::exception
                );
                CHECK(!StringAccess::IsSet(r_background_grid_info, "num_inactive_elements"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_background_grid_info, "num_inactive_elements"), std::exception
                );


                /// elapsed_time_info
                auto& r_elpased_time_info = StringAccess::GetSubDictionary(r_model_info, "elapsed_time_info");

                CHECK(StringAccess::IsSet(r_elpased_time_info, "total"));
                CheckNear(StringAccess::GetValue<double>(r_elpased_time_info, "total"), 0.0, EPS0);
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_elpased_time_info, "total"), std::exception
                );  // Wrong type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<double>(r_elpased_time_info, "total2"), std::exception
                );  // Wrong key

                /// volume_time_info
                auto& r_volume_time_info = StringAccess::GetSubDictionary(r_elpased_time_info, "volume_time_info");
                CHECK(StringAccess::IsSet(r_volume_time_info, "total"));
                CheckNear(StringAccess::GetValue<double>(r_volume_time_info, "total"), 0.0, EPS0);
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_volume_time_info, "total"), std::exception
                );  // Wrong type

                CHECK(StringAccess::IsSet(r_volume_time_info, "classification_of_elements"));
                CheckNear(StringAccess::GetValue<double>(r_volume_time_info, "classification_of_elements"), 0.0, EPS0);
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_volume_time_info, "classification_of_elements"), std::exception
                );  // Wrong type

                CHECK(StringAccess::IsSet(r_volume_time_info, "computation_of_intersections"));
                CheckNear(
                    StringAccess::GetValue<double>(r_volume_time_info, "computation_of_intersections"), 0.0, EPS0
                );
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_volume_time_info, "computation_of_intersections"),
                    std::exception
                );  // Wrong type

                CHECK(StringAccess::IsSet(r_volume_time_info, "solution_of_moment_fitting_eqs"));
                CheckNear(
                    StringAccess::GetValue<double>(r_volume_time_info, "solution_of_moment_fitting_eqs"), 0.0, EPS0
                );
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_volume_time_info, "solution_of_moment_fitting_eqs"),
                    std::exception
                );  // Wrong type

                CHECK(StringAccess::IsSet(r_volume_time_info, "construction_of_ggq_rules"));
                CheckNear(StringAccess::GetValue<double>(r_volume_time_info, "construction_of_ggq_rules"), 0.0, EPS0);
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_volume_time_info, "construction_of_ggq_rules"), std::exception
                );  // Wrong type

                /// conditions_time_info
                auto& r_conditions_time_info =
                    StringAccess::GetSubDictionary(r_elpased_time_info, "conditions_time_info");
                CHECK(StringAccess::IsSet(r_conditions_time_info, "total"));
                CheckNear(StringAccess::GetValue<double>(r_conditions_time_info, "total"), 0.0, EPS0);
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_conditions_time_info, "total"), std::exception
                );  // Wrong type

                /// write_files_time_info
                auto& r_write_files_time_info =
                    StringAccess::GetSubDictionary(r_elpased_time_info, "write_files_time_info");
                CHECK(StringAccess::IsSet(r_write_files_time_info, "total"));
                CheckNear(StringAccess::GetValue<double>(r_write_files_time_info, "total"), 0.0, EPS0);
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_write_files_time_info, "total"), std::exception
                );  // Wrong type
            }
        }


        TEST_CASE("ModelInfoTestSuite - Model Info Customized Values Test")
        {


            {  /// Key access
                auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ComponentInfo");
                auto& r_model_info = *p_model_info;

                /// General r_model_info
                auto& r_embedded_geometry_info = r_model_info[MainInfo::embedded_geometry_info];
                r_embedded_geometry_info.SetValue(EmbeddedGeometryInfo::is_closed, true);
                r_embedded_geometry_info.SetValue(EmbeddedGeometryInfo::volume, 2.0);
                r_embedded_geometry_info.CheckRequired();
                CHECK(r_embedded_geometry_info.IsSet(EmbeddedGeometryInfo::is_closed));
                CHECK((r_embedded_geometry_info.GetRequiredValue<bool>(EmbeddedGeometryInfo::is_closed)) == (true));
                CHECK(r_embedded_geometry_info.IsSet(EmbeddedGeometryInfo::volume));
                CheckNear(r_embedded_geometry_info.GetRequiredValue<double>(EmbeddedGeometryInfo::volume), 2.0, EPS2);

                /// General quadrature_info
                auto& r_quadrature_info = r_model_info[MainInfo::quadrature_info];
                r_quadrature_info.SetValue(QuadratureInfo::represented_volume, 2.0);
                r_quadrature_info.SetValue(QuadratureInfo::percentage_of_geometry_volume, 5.0);
                r_quadrature_info.SetValue(QuadratureInfo::tot_num_points, 3u);
                r_quadrature_info.SetValue(QuadratureInfo::num_of_points_per_full_element, 5.5);
                r_quadrature_info.SetValue(QuadratureInfo::num_of_points_per_trimmed_element, 3.3);
                r_quadrature_info.CheckRequired();
                CHECK(r_quadrature_info.IsSet(QuadratureInfo::represented_volume));
                CheckNear(r_quadrature_info.GetRequiredValue<double>(QuadratureInfo::represented_volume), 2.0, EPS2);
                CHECK(r_quadrature_info.IsSet(QuadratureInfo::percentage_of_geometry_volume));
                CheckNear(
                    r_quadrature_info.GetRequiredValue<double>(QuadratureInfo::percentage_of_geometry_volume), 5.0, EPS2
                );
                CHECK(r_quadrature_info.IsSet(QuadratureInfo::tot_num_points));
                CHECK((r_quadrature_info.GetRequiredValue<IndexType>(QuadratureInfo::tot_num_points)) == (3));
                CHECK(r_quadrature_info.IsSet(QuadratureInfo::num_of_points_per_full_element));
                CheckNear(
                    r_quadrature_info.GetRequiredValue<double>(QuadratureInfo::num_of_points_per_full_element),
                    5.5,
                    EPS2
                );
                CHECK(r_quadrature_info.IsSet(QuadratureInfo::num_of_points_per_trimmed_element));
                CheckNear(
                    r_quadrature_info.GetRequiredValue<double>(QuadratureInfo::num_of_points_per_trimmed_element),
                    3.3,
                    EPS2
                );

                /// General background_grid_info
                auto& r_background_grid_info = r_model_info[MainInfo::background_grid_info];
                r_background_grid_info.SetValue(BackgroundGridInfo::num_active_elements, 5u);
                r_background_grid_info.SetValue(BackgroundGridInfo::num_trimmed_elements, 7u);
                r_background_grid_info.SetValue(BackgroundGridInfo::num_full_elements, 101u);
                r_background_grid_info.SetValue(BackgroundGridInfo::num_inactive_elements, 52u);
                r_background_grid_info.CheckRequired();
                CHECK(r_background_grid_info.IsSet(BackgroundGridInfo::num_active_elements));
                CHECK(
                    (r_background_grid_info.GetRequiredValue<IndexType>(BackgroundGridInfo::num_active_elements)) == (5)
                );
                CHECK(r_background_grid_info.IsSet(BackgroundGridInfo::num_trimmed_elements));
                CHECK(
                    (r_background_grid_info.GetRequiredValue<IndexType>(BackgroundGridInfo::num_trimmed_elements))
                    == (7)
                );
                CHECK(r_background_grid_info.IsSet(BackgroundGridInfo::num_full_elements));
                CHECK(
                    (r_background_grid_info.GetRequiredValue<IndexType>(BackgroundGridInfo::num_full_elements)) == (101)
                );
                CHECK(r_background_grid_info.IsSet(BackgroundGridInfo::num_inactive_elements));
                CHECK(
                    (r_background_grid_info.GetRequiredValue<IndexType>(BackgroundGridInfo::num_inactive_elements))
                    == (52)
                );

                /// elapsed_time_info
                auto& r_elapsed_time_info = r_model_info[MainInfo::elapsed_time_info];
                r_elapsed_time_info.SetValue(ElapsedTimeInfo::total, 101.0);
                r_elapsed_time_info.CheckRequired();
                CHECK(r_elapsed_time_info.IsSet(ElapsedTimeInfo::total));
                CheckNear(r_elapsed_time_info.GetRequiredValue<double>(ElapsedTimeInfo::total), 101.0, EPS2);

                /// volume_time_info
                r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(VolumeTimeInfo::total, 12.2);
                r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(
                    VolumeTimeInfo::classification_of_elements, 122.2
                );
                r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(
                    VolumeTimeInfo::computation_of_intersections, 44.5
                );
                r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(
                    VolumeTimeInfo::solution_of_moment_fitting_eqs, 33.2
                );
                r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].SetValue(
                    VolumeTimeInfo::construction_of_ggq_rules, 2.23
                );
                r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].CheckRequired();
                CHECK(r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(VolumeTimeInfo::total));
                CheckNear(
                    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetRequiredValue<double>(
                        VolumeTimeInfo::total
                    ),
                    12.2,
                    EPS2
                );
                CHECK(
                    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(
                        VolumeTimeInfo::classification_of_elements
                    )
                );
                CheckNear(
                    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetRequiredValue<double>(
                        VolumeTimeInfo::classification_of_elements
                    ),
                    122.2,
                    EPS2
                );
                CHECK(
                    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(
                        VolumeTimeInfo::computation_of_intersections
                    )
                );
                CheckNear(
                    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetRequiredValue<double>(
                        VolumeTimeInfo::computation_of_intersections
                    ),
                    44.5,
                    EPS2
                );
                CHECK(
                    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(
                        VolumeTimeInfo::solution_of_moment_fitting_eqs
                    )
                );
                CheckNear(
                    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetRequiredValue<double>(
                        VolumeTimeInfo::solution_of_moment_fitting_eqs
                    ),
                    33.2,
                    EPS2
                );
                CHECK(
                    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].IsSet(
                        VolumeTimeInfo::construction_of_ggq_rules
                    )
                );
                CheckNear(
                    r_elapsed_time_info[ElapsedTimeInfo::volume_time_info].GetRequiredValue<double>(
                        VolumeTimeInfo::construction_of_ggq_rules
                    ),
                    2.23,
                    EPS2
                );

                /// conditions_time_info
                r_elapsed_time_info[ElapsedTimeInfo::conditions_time_info].SetValue(ConditionsTimeInfo::total, 1.123);
                r_elapsed_time_info[ElapsedTimeInfo::conditions_time_info].CheckRequired();
                CHECK(r_elapsed_time_info[ElapsedTimeInfo::conditions_time_info].IsSet(ConditionsTimeInfo::total));
                CheckNear(
                    r_elapsed_time_info[ElapsedTimeInfo::conditions_time_info].GetRequiredValue<double>(
                        ConditionsTimeInfo::total
                    ),
                    1.123,
                    EPS2
                );

                /// write_files_time_info
                r_elapsed_time_info[ElapsedTimeInfo::write_files_time_info].SetValue(WriteFilesTimeInfo::total, 5.123);
                r_elapsed_time_info[ElapsedTimeInfo::write_files_time_info].CheckRequired();
                CHECK(r_elapsed_time_info[ElapsedTimeInfo::write_files_time_info].IsSet(WriteFilesTimeInfo::total));
                CheckNear(
                    r_elapsed_time_info[ElapsedTimeInfo::write_files_time_info].GetRequiredValue<double>(
                        WriteFilesTimeInfo::total
                    ),
                    5.123,
                    EPS2
                );

                r_elapsed_time_info.CheckRequired();
            }

            {  /// String access
                using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
                using StringAccess = DictionaryStringAccess<DictionaryType>;

                auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ComponentInfo");
                auto& r_model_info = *p_model_info;

                /// General embedded geometry info
                auto& r_embedded_geometry_info = StringAccess::GetSubDictionary(r_model_info, "embedded_geometry_info");
                StringAccess::SetValue(r_embedded_geometry_info, "is_closed", true);
                CHECK(StringAccess::IsSet(r_embedded_geometry_info, "is_closed"));
                CHECK((StringAccess::GetValue<bool>(r_embedded_geometry_info, "is_closed")) == (true));

                StringAccess::SetValue(r_embedded_geometry_info, "volume", 2.0);
                CHECK(StringAccess::IsSet(r_embedded_geometry_info, "volume"));
                CheckNear(StringAccess::GetValue<double>(r_embedded_geometry_info, "volume"), 2.0, EPS2);

                /// General quadrature_info
                auto& r_quadrature_info = StringAccess::GetSubDictionary(r_model_info, "quadrature_info");
                StringAccess::SetValue(r_quadrature_info, "represented_volume", 2.0);
                CHECK(StringAccess::IsSet(r_quadrature_info, "represented_volume"));
                CheckNear(StringAccess::GetValue<double>(r_quadrature_info, "represented_volume"), 2.0, EPS2);

                StringAccess::SetValue(r_quadrature_info, "percentage_of_geometry_volume", 5.0);
                CHECK(StringAccess::IsSet(r_quadrature_info, "percentage_of_geometry_volume"));
                CheckNear(
                    StringAccess::GetValue<double>(r_quadrature_info, "percentage_of_geometry_volume"), 5.0, EPS2
                );

                StringAccess::SetValue(r_quadrature_info, "tot_num_points", 3u);
                CHECK(StringAccess::IsSet(r_quadrature_info, "tot_num_points"));
                CHECK((StringAccess::GetValue<IndexType>(r_quadrature_info, "tot_num_points")) == (3));

                StringAccess::SetValue(r_quadrature_info, "num_of_points_per_full_element", 5.5);
                CHECK(StringAccess::IsSet(r_quadrature_info, "num_of_points_per_full_element"));
                CheckNear(
                    StringAccess::GetValue<double>(r_quadrature_info, "num_of_points_per_full_element"), 5.5, EPS2
                );

                StringAccess::SetValue(r_quadrature_info, "num_of_points_per_trimmed_element", 3.3);
                CHECK(StringAccess::IsSet(r_quadrature_info, "num_of_points_per_trimmed_element"));
                CheckNear(
                    StringAccess::GetValue<double>(r_quadrature_info, "num_of_points_per_trimmed_element"), 3.3, EPS2
                );

                /// General background_grid_info
                auto& r_background_grid_info = StringAccess::GetSubDictionary(r_model_info, "background_grid_info");
                StringAccess::SetValue(r_background_grid_info, "num_active_elements", 5u);
                CHECK(StringAccess::IsSet(r_background_grid_info, "num_active_elements"));
                CHECK((StringAccess::GetValue<IndexType>(r_background_grid_info, "num_active_elements")) == (5));

                StringAccess::SetValue(r_background_grid_info, "num_trimmed_elements", 7u);
                CHECK(StringAccess::IsSet(r_background_grid_info, "num_trimmed_elements"));
                CHECK((StringAccess::GetValue<IndexType>(r_background_grid_info, "num_trimmed_elements")) == (7));

                StringAccess::SetValue(r_background_grid_info, "num_full_elements", 101u);
                CHECK(StringAccess::IsSet(r_background_grid_info, "num_full_elements"));
                CHECK((StringAccess::GetValue<IndexType>(r_background_grid_info, "num_full_elements")) == (101));

                StringAccess::SetValue(r_background_grid_info, "num_inactive_elements", 52u);
                CHECK(StringAccess::IsSet(r_background_grid_info, "num_inactive_elements"));
                CHECK((StringAccess::GetValue<IndexType>(r_background_grid_info, "num_inactive_elements")) == (52));

                /// elapsed_time_info
                auto& r_elapsed_time_info = StringAccess::GetSubDictionary(r_model_info, "elapsed_time_info");
                StringAccess::SetValue(r_elapsed_time_info, "total", 101.0);
                CHECK(StringAccess::IsSet(r_elapsed_time_info, "total"));
                CheckNear(StringAccess::GetValue<double>(r_elapsed_time_info, "total"), 101.0, EPS2);

                /// volume_time_info
                auto& r_volume_time_info = StringAccess::GetSubDictionary(r_elapsed_time_info, "volume_time_info");
                StringAccess::SetValue(r_volume_time_info, "total", 12.2);
                CHECK(StringAccess::IsSet(r_volume_time_info, "total"));
                CheckNear(StringAccess::GetValue<double>(r_volume_time_info, "total"), 12.2, EPS2);

                StringAccess::SetValue(r_volume_time_info, "classification_of_elements", 122.2);
                CHECK(StringAccess::IsSet(r_volume_time_info, "classification_of_elements"));
                CheckNear(
                    StringAccess::GetValue<double>(r_volume_time_info, "classification_of_elements"), 122.2, EPS2
                );

                StringAccess::SetValue(r_volume_time_info, "computation_of_intersections", 44.5);
                CHECK(StringAccess::IsSet(r_volume_time_info, "computation_of_intersections"));
                CheckNear(
                    StringAccess::GetValue<double>(r_volume_time_info, "computation_of_intersections"), 44.5, EPS2
                );

                StringAccess::SetValue(r_volume_time_info, "solution_of_moment_fitting_eqs", 33.2);
                CHECK(StringAccess::IsSet(r_volume_time_info, "solution_of_moment_fitting_eqs"));
                CheckNear(
                    StringAccess::GetValue<double>(r_volume_time_info, "solution_of_moment_fitting_eqs"), 33.2, EPS2
                );

                StringAccess::SetValue(r_volume_time_info, "construction_of_ggq_rules", 2.23);
                CHECK(StringAccess::IsSet(r_volume_time_info, "construction_of_ggq_rules"));
                CheckNear(StringAccess::GetValue<double>(r_volume_time_info, "construction_of_ggq_rules"), 2.23, EPS2);

                /// conditions_time_info
                auto& r_conditions_time_info =
                    StringAccess::GetSubDictionary(r_elapsed_time_info, "conditions_time_info");
                StringAccess::SetValue(r_conditions_time_info, "total", 1.123);
                CHECK(StringAccess::IsSet(r_conditions_time_info, "total"));
                CheckNear(StringAccess::GetValue<double>(r_conditions_time_info, "total"), 1.123, EPS2);

                /// write_files_time_info
                auto& r_write_files_time_info =
                    StringAccess::GetSubDictionary(r_elapsed_time_info, "write_files_time_info");
                StringAccess::SetValue(r_write_files_time_info, "total", 5.123);
                CHECK(StringAccess::IsSet(r_write_files_time_info, "total"));
                CheckNear(StringAccess::GetValue<double>(r_write_files_time_info, "total"), 5.123, EPS2);

                r_elapsed_time_info.CheckRequired();
            }
        }

        TEST_CASE("ModelInfoTestSuite - Model Info Conditions Info Default Values Test"){


            { /// Key access
              auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
        auto& r_cond_info = *p_cond_info;

        REQUIRE_THROWS_AS(r_cond_info.CheckRequired(), std::exception);

        CHECK(!r_cond_info.IsSet(ConditionInfo::condition_id));
        CHECK(!r_cond_info.IsSet(ConditionInfo::surf_area));
        CHECK(!r_cond_info.IsSet(ConditionInfo::perc_surf_area_in_active_domain));
        if (!NOTDEBUG) {
            REQUIRE_THROWS_AS(r_cond_info.GetRequiredValue<IndexType>(ConditionInfo::condition_id), std::exception);
            REQUIRE_THROWS_AS(r_cond_info.GetRequiredValue<double>(ConditionInfo::surf_area), std::exception);
            REQUIRE_THROWS_AS(
                r_cond_info.GetRequiredValue<double>(ConditionInfo::perc_surf_area_in_active_domain), std::exception
            );
        }
    }
    {  /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
        auto& r_cond_info = *p_cond_info;

        REQUIRE_THROWS_AS(r_cond_info.CheckRequired(), std::exception);

        CHECK(!StringAccess::IsSet(r_cond_info, "condition_id"));
        REQUIRE_THROWS_AS(StringAccess::GetValue<IndexType>(r_cond_info, "condition_id"), std::exception);
        CHECK(!StringAccess::IsSet(r_cond_info, "surf_area"));
        REQUIRE_THROWS_AS(StringAccess::GetValue<double>(r_cond_info, "surf_area"), std::exception);
        CHECK(!StringAccess::IsSet(r_cond_info, "perc_surf_area_in_active_domain"));
        REQUIRE_THROWS_AS(
            StringAccess::GetValue<double>(r_cond_info, "perc_surf_area_in_active_domain"), std::exception
        );
    }
};

TEST_CASE("ModelInfoTestSuite - Model Info Conditions Info Customized Values Test"){


    { /// Key access
      auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ComponentInfo");

{
    auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
    auto& r_cond_info = *p_cond_info;

    r_cond_info.SetValue(ConditionInfo::condition_id, 101u);
    r_cond_info.SetValue(ConditionInfo::surf_area, 22.2);
    r_cond_info.SetValue(ConditionInfo::perc_surf_area_in_active_domain, 223.1);
    r_cond_info.CheckRequired();
    CHECK(r_cond_info.IsSet(ConditionInfo::condition_id));
    CHECK((r_cond_info.GetRequiredValue<IndexType>(ConditionInfo::condition_id)) == (101));

    CHECK(r_cond_info.IsSet(ConditionInfo::surf_area));
    CheckNear(r_cond_info.GetRequiredValue<double>(ConditionInfo::surf_area), 22.2, EPS2);

    CHECK(r_cond_info.IsSet(ConditionInfo::perc_surf_area_in_active_domain));
    CheckNear(r_cond_info.GetRequiredValue<double>(ConditionInfo::perc_surf_area_in_active_domain), 223.1, EPS2);

    auto& r_cond_list = p_model_info->GetList(MainInfo::conditions_infos_list);
    r_cond_list.push_back(std::move(p_cond_info));
}
{
    auto& p_cond_info = p_model_info->GetList(MainInfo::conditions_infos_list)[0];
    auto& r_cond_info = *p_cond_info;
    r_cond_info.CheckRequired();

    CHECK(r_cond_info.IsSet(ConditionInfo::condition_id));
    CHECK((r_cond_info.GetRequiredValue<IndexType>(ConditionInfo::condition_id)) == (101));

    CHECK(r_cond_info.IsSet(ConditionInfo::surf_area));
    CheckNear(r_cond_info.GetRequiredValue<double>(ConditionInfo::surf_area), 22.2, EPS2);

    CHECK(r_cond_info.IsSet(ConditionInfo::perc_surf_area_in_active_domain));
    CheckNear(r_cond_info.GetRequiredValue<double>(ConditionInfo::perc_surf_area_in_active_domain), 223.1, EPS2);

    r_cond_info.CheckRequired();
}
}  // namespace queso
{  /// String access
    using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
    using StringAccess = DictionaryStringAccess<DictionaryType>;

    auto p_model_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ComponentInfo");
    {
        auto p_cond_info = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionInfo");
        auto& r_cond_info = *p_cond_info;

        StringAccess::SetValue(r_cond_info, "condition_id", 101);
        CHECK(StringAccess::IsSet(r_cond_info, "condition_id"));
        CHECK((StringAccess::GetValue<IndexType>(r_cond_info, "condition_id")) == (101));

        StringAccess::SetValue(r_cond_info, "surf_area", 22.2);
        CHECK(StringAccess::IsSet(r_cond_info, "surf_area"));
        CheckNear(StringAccess::GetValue<double>(r_cond_info, "surf_area"), 22.2, EPS2);

        StringAccess::SetValue(r_cond_info, "perc_surf_area_in_active_domain", 223.1);
        CHECK(StringAccess::IsSet(r_cond_info, "perc_surf_area_in_active_domain"));
        CheckNear(StringAccess::GetValue<double>(r_cond_info, "perc_surf_area_in_active_domain"), 223.1, EPS2);

        r_cond_info.CheckRequired();

        auto& r_cond_list = StringAccess::GetList(*p_model_info, "conditions_infos_list");
        r_cond_list.push_back(std::move(p_cond_info));
    }
    {
        const auto& p_cond_list = StringAccess::GetList(*p_model_info, "conditions_infos_list")[0];
        const auto& r_cond_info = *p_cond_list;

        CHECK(StringAccess::IsSet(r_cond_info, "condition_id"));
        CHECK((StringAccess::GetValue<IndexType>(r_cond_info, "condition_id")) == (101));

        CHECK(StringAccess::IsSet(r_cond_info, "surf_area"));
        CheckNear(StringAccess::GetValue<double>(r_cond_info, "surf_area"), 22.2, EPS2);

        CHECK(StringAccess::IsSet(r_cond_info, "perc_surf_area_in_active_domain"));
        CheckNear(StringAccess::GetValue<double>(r_cond_info, "perc_surf_area_in_active_domain"), 223.1, EPS2);

        r_cond_info.CheckRequired();
    }
}
}
;

}  // TEST_SUITE

} // namespace queso::testing
