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

    TEST_SUITE("SettingsTestSuite")
    {

        TEST_CASE("SettingsTestSuite - Settings Wrong Type Test")
        {


            {  /// Key access
                auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
                auto& settings = *p_settings;

                auto& r_mesh_settings = settings[MainSettings::background_grid_settings];
                REQUIRE_THROWS_AS(r_mesh_settings.CheckRequired(), std::exception);

                if (!NOTDEBUG) {
                    REQUIRE_THROWS_AS(settings[MainInfo::background_grid_info], std::exception);  // Wrong Key type

                    // Main settings
                    REQUIRE_THROWS_AS(
                        settings.IsSet(ConditionSettings::condition_id), std::exception
                    );  // Wrong Key type
                    REQUIRE_THROWS_AS(
                        settings.GetRequiredValue<PointType>(BackgroundGridSettings::lower_bound_uvw), std::exception
                    );  // Wrong Key type

                    /// Mesh settings
                    REQUIRE_THROWS_AS(
                        r_mesh_settings.GetRequiredValue<IndexType>(MainSettings::echo_level), std::exception
                    );  // Wrong Key type

                    REQUIRE_THROWS_AS(
                        r_mesh_settings.SetValue(MainSettings::echo_level, 2u), std::exception
                    );  // Wrong Key type

                    r_mesh_settings.SetValue(BackgroundGridSettings::grid_type, GridType::b_spline_grid);
                    r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_xyz, PointType{ 1.0, 2.0, 3.0 });
                    REQUIRE_THROWS_AS(
                        r_mesh_settings.SetValue(MainSettings::echo_level, IndexType(2)), std::exception
                    );  // Wrong Key type
                    CHECK(r_mesh_settings.IsSet(BackgroundGridSettings::lower_bound_xyz));
                    REQUIRE_THROWS_AS(
                        r_mesh_settings.CheckRequired(), std::exception
                    );  // Not all required values are set.

                    r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_xyz, PointType{ 1.0, 2.0, 3.0 });
                    CHECK(r_mesh_settings.IsSet(BackgroundGridSettings::upper_bound_xyz));
                    REQUIRE_THROWS_AS(
                        r_mesh_settings.CheckRequired(), std::exception
                    );  // Not all required values are set.

                    r_mesh_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType{ 1.0, 2.0, 3.0 });
                    CHECK(r_mesh_settings.IsSet(BackgroundGridSettings::lower_bound_uvw));
                    REQUIRE_THROWS_AS(
                        r_mesh_settings.CheckRequired(), std::exception
                    );  // Not all required values are set.

                    r_mesh_settings.SetValue(BackgroundGridSettings::upper_bound_uvw, PointType{ 1.0, 2.0, 3.0 });
                    CHECK(r_mesh_settings.IsSet(BackgroundGridSettings::upper_bound_uvw));
                    REQUIRE_THROWS_AS(
                        r_mesh_settings.CheckRequired(), std::exception
                    );  // Not all required values are set.

                    r_mesh_settings.SetValue(BackgroundGridSettings::polynomial_order, Vector3i{ 1, 2, 3 });
                    CHECK(r_mesh_settings.IsSet(BackgroundGridSettings::polynomial_order));
                    REQUIRE_THROWS_AS(
                        r_mesh_settings.CheckRequired(), std::exception
                    );  // Not all required values are set.

                    r_mesh_settings.SetValue(BackgroundGridSettings::number_of_elements, Vector3i{ 1, 2, 3 });
                    CHECK(r_mesh_settings.IsSet(BackgroundGridSettings::number_of_elements));
                    r_mesh_settings.CheckRequired();

                    CHECK(
                        (r_mesh_settings.GetRequiredValue<GridType>(BackgroundGridSettings::grid_type))
                        == (GridType::b_spline_grid)
                    );
                    CheckPointNear(
                        r_mesh_settings.GetRequiredValue<PointType>(BackgroundGridSettings::lower_bound_xyz),
                        PointType({ 1.0, 2.0, 3.0 }),
                        1e-10
                    );
                    CheckPointNear(
                        r_mesh_settings.GetRequiredValue<PointType>(BackgroundGridSettings::upper_bound_xyz),
                        PointType({ 1.0, 2.0, 3.0 }),
                        1e-10
                    );
                    CheckPointNear(
                        r_mesh_settings.GetRequiredValue<PointType>(BackgroundGridSettings::lower_bound_uvw),
                        PointType({ 1.0, 2.0, 3.0 }),
                        1e-10
                    );
                    CheckPointNear(
                        r_mesh_settings.GetRequiredValue<PointType>(BackgroundGridSettings::upper_bound_uvw),
                        PointType({ 1.0, 2.0, 3.0 }),
                        1e-10
                    );
                    CHECK(
                        r_mesh_settings.GetRequiredValue<Vector3i>(BackgroundGridSettings::polynomial_order)
                        == Vector3i({ 1, 2, 3 })
                    );
                    CHECK(
                        r_mesh_settings.GetRequiredValue<Vector3i>(BackgroundGridSettings::number_of_elements)
                        == Vector3i({ 1, 2, 3 })
                    );
                }
            }
            {  /// String access
                using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
                using StringAccess = DictionaryStringAccess<DictionaryType>;

                auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
                auto& settings = *p_settingss;
                REQUIRE_THROWS_AS(settings.CheckRequired(), std::exception);

                REQUIRE_THROWS_AS(
                    StringAccess::GetSubDictionary(settings, "echo_level"), std::exception
                );  // Wrong Key type

                /// Main settings
                REQUIRE_THROWS_AS(StringAccess::IsSet(settings, "general_settings"), std::exception);  // Wrong key
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(settings, "lower_bound_uvw"), std::exception
                );  // Wrong Key type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(settings, "input_filename"), std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(settings, "output_directory_name"), std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(settings, "echo_level"), std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(settings, "echo_level", -1), std::exception
                );  // Negative value
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(settings, "write_output_to_file"), std::exception
                );  // Wrong Value type

                /// Mesh settings
                auto& r_mesh_settings = StringAccess::GetSubDictionary(settings, "background_grid_settings");
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_mesh_settings, "echo_level"), std::exception
                );  // Wrong Key type

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(r_mesh_settings, "lower_bound_xyz"), std::exception
                );  // Not set
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "echo_level", 2u), std::exception
                );  // Wrong Key type
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "lower_bound_xyz", IndexType(2)), std::exception
                );  // Wrong Value type
                StringAccess::SetValue(r_mesh_settings, "lower_bound_xyz", PointType{ 1.0, 2.0, 3.0 });
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "echo_level", IndexType(2)), std::exception
                );  // Wrong Key type
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "lower_bound_xyz", IndexType(2)), std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<Vector3i>(r_mesh_settings, "lower_bound_xyz"), std::exception
                );  // Wrong Value type
                const PointType& r_lower_bound_xyz =
                    StringAccess::GetValue<PointType>(r_mesh_settings, "lower_bound_xyz");
                CheckPointNear(r_lower_bound_xyz, PointType({ 1.0, 2.0, 3.0 }), 1e-10);

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(r_mesh_settings, "upper_bound_xyz"), std::exception
                );  // Not set
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "upper_bound_xyz", IndexType(2)), std::exception
                );  // Wrong Value type
                StringAccess::SetValue(r_mesh_settings, "upper_bound_xyz", PointType{ 1.0, 2.0, 3.0 });
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "upper_bound_xyz", IndexType(2)), std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<Vector3i>(r_mesh_settings, "upper_bound_xyz"), std::exception
                );  // Wrong Value type
                const PointType& r_upper_bound_xyz =
                    StringAccess::GetValue<PointType>(r_mesh_settings, "upper_bound_xyz");
                CheckPointNear(r_upper_bound_xyz, PointType({ 1.0, 2.0, 3.0 }), 1e-10);

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(r_mesh_settings, "lower_bound_uvw"), std::exception
                );  // Not set
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "lower_bound_uvw", IndexType(2)), std::exception
                );  // Wrong Value type
                StringAccess::SetValue(r_mesh_settings, "lower_bound_uvw", PointType{ 1.0, 2.0, 3.0 });
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "lower_bound_uvw", IndexType(2)), std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<Vector3i>(r_mesh_settings, "lower_bound_uvw"), std::exception
                );  // Wrong Value type
                const PointType& r_lower_bound_uvw =
                    StringAccess::GetValue<PointType>(r_mesh_settings, "lower_bound_uvw");
                CheckPointNear(r_lower_bound_uvw, PointType({ 1.0, 2.0, 3.0 }), 1e-10);

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(r_mesh_settings, "upper_bound_uvw"), std::exception
                );  // Not set
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "upper_bound_uvw", IndexType(2)), std::exception
                );  // Wrong Value type
                StringAccess::SetValue(r_mesh_settings, "upper_bound_uvw", PointType{ 1.0, 2.0, 3.0 });
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "upper_bound_uvw", IndexType(2)), std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<Vector3i>(r_mesh_settings, "upper_bound_uvw"), std::exception
                );  // Wrong Value type
                const PointType& r_upper_bound_uvw =
                    StringAccess::GetValue<PointType>(r_mesh_settings, "upper_bound_uvw");
                CheckPointNear(r_upper_bound_uvw, PointType({ 1.0, 2.0, 3.0 }), 1e-10);

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(r_mesh_settings, "polynomial_order"), std::exception
                );  // Not set
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "polynomial_order", IndexType(2)), std::exception
                );  // Wrong Value type
                StringAccess::SetValue(r_mesh_settings, "polynomial_order", Vector3i{ 1, 2, 3 });
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "polynomial_order", IndexType(2)), std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(r_mesh_settings, "polynomial_order"), std::exception
                );  // Wrong Value type
                const Vector3i& r_polynomial_order =
                    StringAccess::GetValue<Vector3i>(r_mesh_settings, "polynomial_order");
                CHECK(r_polynomial_order == Vector3i({ 1, 2, 3 }));

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(r_mesh_settings, "number_of_elements"), std::exception
                );  // Not set
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "number_of_elements", IndexType(2)), std::exception
                );  // Wrong Value type
                StringAccess::SetValue(r_mesh_settings, "number_of_elements", Vector3i{ 1, 2, 3 });
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_mesh_settings, "number_of_elements", IndexType(2)), std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(r_mesh_settings, "number_of_elements"), std::exception
                );  // Wrong Value type
                const Vector3i& number_of_elements =
                    StringAccess::GetValue<Vector3i>(r_mesh_settings, "number_of_elements");
                CHECK(number_of_elements == Vector3i({ 1, 2, 3 }));

                /// Trimmed quadrature rule settings
                auto& r_trimmed_quad_rule_settings =
                    StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings");
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_trimmed_quad_rule_settings, "number_of_elements", Vector3i{ 1, 2, 3 }),
                    std::exception
                );  // Wrong Key type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<Vector3i>(r_trimmed_quad_rule_settings, "moment_fitting_residual"),
                    std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<Vector3i>(r_trimmed_quad_rule_settings, "min_element_volume_ratio"),
                    std::exception
                );  // Wrong Value type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<Vector3i>(r_trimmed_quad_rule_settings, "min_num_boundary_triangles"),
                    std::exception
                );  // Wrong Value type

                /// Non trimmed quadrature rule settings
                auto& r_non_trimmed_quad_rule_settings =
                    StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings");
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(r_non_trimmed_quad_rule_settings, "number_of_elements", Vector3i{ 1, 2, 3 }),
                    std::exception
                );  // Wrong Key type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(r_non_trimmed_quad_rule_settings, "integration_method"),
                    std::exception
                );  // Wrong Value type
            }
        }

        TEST_CASE("SettingsTestSuite - Settings Check Required Test")
        {


            auto p_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
            auto& settings = *p_settings;

            REQUIRE_THROWS_AS(settings.CheckRequired(), std::exception);

            settings.SetValue(MainSettings::component_name, std::string("main"));
            settings.SetValue(MainSettings::input_filename, std::string("input.stl"));
            settings[MainSettings::background_grid_settings].SetValue(
                BackgroundGridSettings::grid_type, GridType::b_spline_grid
            );
            settings[MainSettings::background_grid_settings].SetValue(
                BackgroundGridSettings::lower_bound_xyz, PointType({ 0.0, 0.0, 0.0 })
            );
            settings[MainSettings::background_grid_settings].SetValue(
                BackgroundGridSettings::upper_bound_xyz, PointType({ 1.0, 1.0, 1.0 })
            );
            settings[MainSettings::background_grid_settings].SetValue(
                BackgroundGridSettings::lower_bound_uvw, PointType({ 0.0, 0.0, 0.0 })
            );
            settings[MainSettings::background_grid_settings].SetValue(
                BackgroundGridSettings::upper_bound_uvw, PointType({ 1.0, 1.0, 1.0 })
            );
            settings[MainSettings::background_grid_settings].SetValue(
                BackgroundGridSettings::polynomial_order, Vector3i({ 2, 2, 2 })
            );
            settings[MainSettings::background_grid_settings].SetValue(
                BackgroundGridSettings::number_of_elements, Vector3i({ 4, 4, 4 })
            );

            settings.CheckRequired();

            auto p_condition_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
            p_condition_settings->SetValue(ConditionSettings::condition_id, 5u);
            settings.GetList(MainSettings::conditions_settings_list).push_back(std::move(p_condition_settings));

            REQUIRE_THROWS_AS(settings.CheckRequired(), std::exception);

            settings.GetList(MainSettings::conditions_settings_list)[0]->SetValue(
                ConditionSettings::input_filename, std::string("condition.stl")
            );
            settings.CheckRequired();
        }

        TEST_CASE("SettingsTestSuite - Settings Default Values Test")
        {


            {  /// Key access
                auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
                auto& settings = *p_settingss;

                /// Main settings
                CHECK(!settings.IsSet(MainSettings::component_name));
                CHECK(!settings.IsSet(MainSettings::input_filename));
                REQUIRE_THROWS_AS(settings.CheckRequired(), std::exception);
                settings.SetValue(MainSettings::component_name, std::string("main"));
                settings.SetValue(MainSettings::input_filename, std::string("input.stl"));
                CHECK(settings.IsSet(MainSettings::output_directory_name));
                CHECK(settings.IsSet(MainSettings::echo_level));
                CHECK(settings.IsSet(MainSettings::write_output_to_file));
                /// Mesh settings
                CHECK(!settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::grid_type));
                REQUIRE_THROWS_AS(settings[MainSettings::background_grid_settings].CheckRequired(), std::exception);
                CHECK(!settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::lower_bound_xyz));
                CHECK(!settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::upper_bound_xyz));
                CHECK(!settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::lower_bound_uvw));
                CHECK(!settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::upper_bound_uvw));
                CHECK(
                    !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::polynomial_order)
                );
                CHECK(
                    !settings[MainSettings::background_grid_settings].IsSet(BackgroundGridSettings::number_of_elements)
                );
                if (!NOTDEBUG) {
                    REQUIRE_THROWS_AS(
                        settings[MainSettings::background_grid_settings].GetRequiredValue<GridType>(
                            BackgroundGridSettings::grid_type
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        settings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                            BackgroundGridSettings::lower_bound_xyz
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        settings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                            BackgroundGridSettings::upper_bound_xyz
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        settings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                            BackgroundGridSettings::lower_bound_uvw
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        settings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                            BackgroundGridSettings::upper_bound_uvw
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        settings[MainSettings::background_grid_settings].GetRequiredValue<Vector3i>(
                            BackgroundGridSettings::polynomial_order
                        ),
                        std::exception
                    );
                    REQUIRE_THROWS_AS(
                        settings[MainSettings::background_grid_settings].GetRequiredValue<Vector3i>(
                            BackgroundGridSettings::number_of_elements
                        ),
                        std::exception
                    );
                }
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::grid_type, GridType::b_spline_grid
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::lower_bound_xyz, PointType({ 0.0, 0.0, 0.0 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::upper_bound_xyz, PointType({ 1.0, 1.0, 1.0 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::lower_bound_uvw, PointType({ 0.0, 0.0, 0.0 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::upper_bound_uvw, PointType({ 1.0, 1.0, 1.0 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::polynomial_order, Vector3i({ 2, 2, 2 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::number_of_elements, Vector3i({ 4, 4, 4 })
                );
                settings.CheckRequired();

                CHECK(
                    (settings.GetRequiredValue<std::string>(MainSettings::output_directory_name))
                    == (std::string("queso_output"))
                );
                CHECK((settings.GetRequiredValue<IndexType>(MainSettings::echo_level)) == (1u));
                CHECK((settings.GetRequiredValue<bool>(MainSettings::write_output_to_file)) == (true));

                /// TrimmedQuadratureRuleSettings settings
                settings[MainSettings::trimmed_quadrature_rule_settings].CheckRequired();
                CHECK(
                    settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(
                        TrimmedQuadratureRuleSettings::moment_fitting_residual
                    )
                );
                CheckRelativeNear(
                    settings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<double>(
                        TrimmedQuadratureRuleSettings::moment_fitting_residual
                    ),
                    1e-10,
                    1e-10
                );

                CHECK(
                    settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(
                        TrimmedQuadratureRuleSettings::min_element_volume_ratio
                    )
                );
                CheckRelativeNear(
                    settings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<double>(
                        TrimmedQuadratureRuleSettings::min_element_volume_ratio
                    ),
                    1e-3,
                    1e-10
                );

                CHECK(
                    settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(
                        TrimmedQuadratureRuleSettings::min_num_boundary_triangles
                    )
                );
                CHECK(
                    (settings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<IndexType>(
                        TrimmedQuadratureRuleSettings::min_num_boundary_triangles
                    ))
                    == (100)
                );

                CHECK(
                    settings[MainSettings::trimmed_quadrature_rule_settings].IsSet(
                        TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed
                    )
                );
                CHECK(
                    (settings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<bool>(
                        TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed
                    ))
                    == (true)
                );

                // NonTrimmedQuadratureRuleSettings settings
                settings[MainSettings::non_trimmed_quadrature_rule_settings].CheckRequired();
                CHECK(
                    settings[MainSettings::non_trimmed_quadrature_rule_settings].IsSet(
                        NonTrimmedQuadratureRuleSettings::integration_method
                    )
                );
                CHECK(
                    (settings[MainSettings::non_trimmed_quadrature_rule_settings].GetRequiredValue<IntegrationMethod>(
                        NonTrimmedQuadratureRuleSettings::integration_method
                    ))
                    == (IntegrationMethod::gauss)
                );
            }
            {  /// String access
                using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
                using StringAccess = DictionaryStringAccess<DictionaryType>;

                auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
                auto& settings = *p_settingss;

                /// General settings
                CHECK(!StringAccess::IsSet(settings, "component_name"));
                CHECK(!StringAccess::IsSet(settings, "input_filename"));
                REQUIRE_THROWS_AS(StringAccess::GetValue<std::string>(settings, "input_filename"), std::exception);

                CHECK(StringAccess::IsSet(settings, "output_directory_name"));
                CHECK(
                    (StringAccess::GetValue<std::string>(settings, "output_directory_name"))
                    == (std::string("queso_output"))
                );

                CHECK(StringAccess::IsSet(settings, "echo_level"));
                CHECK((StringAccess::GetValue<IndexType>(settings, "echo_level")) == (1u));

                CHECK(StringAccess::IsSet(settings, "write_output_to_file"));
                CHECK((StringAccess::GetValue<bool>(settings, "write_output_to_file")) == (true));

                /// Mesh settings
                CHECK(!StringAccess::IsSet(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"), "grid_type"
                ));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "grid_type"
                    ),
                    std::exception
                );

                CHECK(!StringAccess::IsSet(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz"
                ));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz"
                    ),
                    std::exception
                );

                CHECK(!StringAccess::IsSet(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_xyz"
                ));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_xyz"
                    ),
                    std::exception
                );

                CHECK(!StringAccess::IsSet(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_uvw"
                ));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_uvw"
                    ),
                    std::exception
                );

                CHECK(!StringAccess::IsSet(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_uvw"
                ));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_uvw"
                    ),
                    std::exception
                );

                CHECK(!StringAccess::IsSet(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"), "polynomial_order"
                ));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<Vector3i>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "polynomial_order"
                    ),
                    std::exception
                );

                CHECK(!StringAccess::IsSet(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"), "number_of_elements"
                ));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<Vector3i>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "number_of_elements"
                    ),
                    std::exception
                );

                /// TrimmedQuadratureRuleSettings settings
                CHECK(
                    StringAccess::IsSet(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "moment_fitting_residual"
                    )
                );
                CheckRelativeNear(
                    StringAccess::GetValue<double>(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "moment_fitting_residual"
                    ),
                    1e-10,
                    1e-10
                );

                CHECK(
                    StringAccess::IsSet(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "min_element_volume_ratio"
                    )
                );
                CheckRelativeNear(
                    StringAccess::GetValue<double>(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "min_element_volume_ratio"
                    ),
                    1e-3,
                    1e-10
                );

                CHECK(
                    StringAccess::IsSet(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "min_num_boundary_triangles"
                    )
                );
                CHECK(
                    (StringAccess::GetValue<IndexType>(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "min_num_boundary_triangles"
                    ))
                    == (100)
                );

                CHECK(
                    StringAccess::IsSet(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "neglect_elements_if_stl_is_flawed"
                    )
                );
                CHECK(
                    (StringAccess::GetValue<bool>(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "neglect_elements_if_stl_is_flawed"
                    ))
                    == (true)
                );

                // NonTrimmedQuadratureRuleSettings settings
                CHECK(
                    StringAccess::IsSet(
                        StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings"),
                        "integration_method"
                    )
                );
                CHECK(
                    (StringAccess::GetValue<IntegrationMethod>(
                        StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings"),
                        "integration_method"
                    ))
                    == (IntegrationMethod::gauss)
                );
            }
        }


        TEST_CASE("SettingsTestSuite - Settings Customized Values Test")
        {


            {  /// Key access
                auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
                auto& settings = *p_settingss;

                /// Main settings
                settings.SetValue(MainSettings::component_name, std::string("main"));
                settings.SetValue(MainSettings::input_filename, std::string("test_filename.stl"));
                settings.SetValue(MainSettings::output_directory_name, std::string("new_output/"));
                settings.SetValue(MainSettings::echo_level, 2u);
                settings.SetValue(MainSettings::write_output_to_file, false);

                /// Mesh settings
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::grid_type, GridType::b_spline_grid
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::lower_bound_xyz, PointType({ 1.0, 1.0, 2.0 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::upper_bound_xyz, PointType({ 0.0, 3.0, 2.0 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::lower_bound_uvw, PointType({ 0.1, -1.0, 2.0 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::upper_bound_uvw, PointType({ 0.66, 1.0, 2.2 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::polynomial_order, Vector3i({ 5, 6, 7 })
                );
                settings[MainSettings::background_grid_settings].SetValue(
                    BackgroundGridSettings::number_of_elements, Vector3i({ 8, 9, 2 })
                );
                settings[MainSettings::background_grid_settings].CheckRequired();

                CHECK(
                    (settings[MainSettings::background_grid_settings].GetRequiredValue<GridType>(
                        BackgroundGridSettings::grid_type
                    ))
                    == (GridType::b_spline_grid)
                );
                CHECK(
                    (settings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                        BackgroundGridSettings::lower_bound_xyz
                    ))
                    == (PointType({ 1.0, 1.0, 2.0 }))
                );
                CHECK(
                    (settings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                        BackgroundGridSettings::upper_bound_xyz
                    ))
                    == (PointType({ 0.0, 3.0, 2.0 }))
                );
                CHECK(
                    (settings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                        BackgroundGridSettings::lower_bound_uvw
                    ))
                    == (PointType({ 0.1, -1.0, 2.0 }))
                );
                CHECK(
                    (settings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                        BackgroundGridSettings::upper_bound_uvw
                    ))
                    == (PointType({ 0.66, 1.0, 2.2 }))
                );
                CHECK(
                    (settings[MainSettings::background_grid_settings].GetRequiredValue<Vector3i>(
                        BackgroundGridSettings::polynomial_order
                    ))
                    == (Vector3i({ 5, 6, 7 }))
                );
                CHECK(
                    (settings[MainSettings::background_grid_settings].GetRequiredValue<Vector3i>(
                        BackgroundGridSettings::number_of_elements
                    ))
                    == (Vector3i({ 8, 9, 2 }))
                );


                /// Trimmed quadrature rule settings
                settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(
                    TrimmedQuadratureRuleSettings::moment_fitting_residual, 5.6e-5
                );
                settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(
                    TrimmedQuadratureRuleSettings::min_element_volume_ratio, 0.45
                );
                settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(
                    TrimmedQuadratureRuleSettings::min_num_boundary_triangles, IndexType(234)
                );
                settings[MainSettings::trimmed_quadrature_rule_settings].SetValue(
                    TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed, false
                );
                settings[MainSettings::trimmed_quadrature_rule_settings].CheckRequired();

                CheckRelativeNear(
                    settings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<double>(
                        TrimmedQuadratureRuleSettings::moment_fitting_residual
                    ),
                    5.6e-5,
                    1e-10
                );
                CheckRelativeNear(
                    settings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<double>(
                        TrimmedQuadratureRuleSettings::min_element_volume_ratio
                    ),
                    0.45,
                    1e-10
                );
                CHECK(
                    (settings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<IndexType>(
                        TrimmedQuadratureRuleSettings::min_num_boundary_triangles
                    ))
                    == (234u)
                );
                CHECK(
                    (settings[MainSettings::trimmed_quadrature_rule_settings].GetRequiredValue<bool>(
                        TrimmedQuadratureRuleSettings::neglect_elements_if_stl_is_flawed
                    ))
                    == (false)
                );

                /// Non trimmed quadrature rule settings
                settings[MainSettings::non_trimmed_quadrature_rule_settings].SetValue(
                    NonTrimmedQuadratureRuleSettings::integration_method, IntegrationMethod::ggq_optimal
                );
                settings[MainSettings::non_trimmed_quadrature_rule_settings].CheckRequired();

                settings.CheckRequired();

                CHECK(
                    (settings.GetRequiredValue<std::string>(MainSettings::input_filename))
                    == (std::string("test_filename.stl"))
                );
                CHECK(
                    (settings.GetRequiredValue<std::string>(MainSettings::output_directory_name))
                    == (std::string("new_output/"))
                );
                CHECK((settings.GetRequiredValue<IndexType>(MainSettings::echo_level)) == (2u));
                CHECK((settings.GetRequiredValue<bool>(MainSettings::write_output_to_file)) == (false));

                CHECK(
                    (settings[MainSettings::non_trimmed_quadrature_rule_settings].GetRequiredValue<IntegrationMethod>(
                        NonTrimmedQuadratureRuleSettings::integration_method
                    ))
                    == (IntegrationMethod::ggq_optimal)
                );
            }
            {  /// String acccess
                using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
                using StringAccess = DictionaryStringAccess<DictionaryType>;

                auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
                auto& settings = *p_settingss;

                /// General settings
                StringAccess::SetValue(settings, "component_name", std::string("main"));
                StringAccess::SetValue(settings, "input_filename", std::string("test_filename.stl"));

                StringAccess::SetValue(settings, "output_directory_name", std::string("new_output/"));
                StringAccess::SetValue(settings, "echo_level", 2u);
                StringAccess::SetValue(settings, "write_output_to_file", false);

                CHECK(
                    (StringAccess::GetValue<std::string>(settings, "input_filename"))
                    == (std::string("test_filename.stl"))
                );
                CHECK(
                    (StringAccess::GetValue<std::string>(settings, "output_directory_name"))
                    == (std::string("new_output/"))
                );
                CHECK((StringAccess::GetValue<IndexType>(settings, "echo_level")) == (2u));
                CHECK((StringAccess::GetValue<bool>(settings, "write_output_to_file")) == (false));

                /// Mesh settings
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"),
                    "grid_type",
                    GridType::b_spline_grid
                );
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"),
                    "lower_bound_xyz",
                    PointType({ 1.0, 1.0, 2.0 })
                );
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"),
                    "upper_bound_xyz",
                    PointType({ 0.0, 3.0, 2.0 })
                );
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"),
                    "lower_bound_uvw",
                    PointType({ 0.1, -1.0, 2.0 })
                );
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"),
                    "upper_bound_uvw",
                    PointType({ 0.66, 1.0, 2.2 })
                );
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"),
                    "polynomial_order",
                    Vector3i({ 5, 6, 7 })
                );
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "background_grid_settings"),
                    "number_of_elements",
                    Vector3i({ 8, 9, 2 })
                );

                CHECK(
                    (StringAccess::GetValue<GridType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "grid_type"
                    ))
                    == (GridType::b_spline_grid)
                );
                CHECK(
                    (StringAccess::GetValue<PointType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz"
                    ))
                    == (PointType({ 1.0, 1.0, 2.0 }))
                );
                CHECK(
                    (StringAccess::GetValue<PointType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_xyz"
                    ))
                    == (PointType({ 0.0, 3.0, 2.0 }))
                );
                CHECK(
                    (StringAccess::GetValue<PointType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_uvw"
                    ))
                    == (PointType({ 0.1, -1.0, 2.0 }))
                );
                CHECK(
                    (StringAccess::GetValue<PointType>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "upper_bound_uvw"
                    ))
                    == (PointType({ 0.66, 1.0, 2.2 }))
                );
                CHECK(
                    (StringAccess::GetValue<Vector3i>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "polynomial_order"
                    ))
                    == (Vector3i({ 5, 6, 7 }))
                );
                CHECK(
                    (StringAccess::GetValue<Vector3i>(
                        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "number_of_elements"
                    ))
                    == (Vector3i({ 8, 9, 2 }))
                );

                /// Trimmed quadrature rule settings
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                    "moment_fitting_residual",
                    5.6e-5
                );
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                    "min_element_volume_ratio",
                    0.45
                );
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                    "min_num_boundary_triangles",
                    IndexType(234)
                );
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                    "neglect_elements_if_stl_is_flawed",
                    false
                );

                CheckRelativeNear(
                    StringAccess::GetValue<double>(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "moment_fitting_residual"
                    ),
                    5.6e-5,
                    1e-10
                );
                CheckRelativeNear(
                    StringAccess::GetValue<double>(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "min_element_volume_ratio"
                    ),
                    0.45,
                    1e-10
                );
                CHECK(
                    (StringAccess::GetValue<IndexType>(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "min_num_boundary_triangles"
                    ))
                    == (234u)
                );
                CHECK(
                    (StringAccess::GetValue<bool>(
                        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"),
                        "neglect_elements_if_stl_is_flawed"
                    ))
                    == (false)
                );

                /// Non trimmed quadrature rule settings
                StringAccess::SetValue(
                    StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings"),
                    "integration_method",
                    IntegrationMethod::ggq_optimal
                );
                CHECK(
                    (StringAccess::GetValue<IntegrationMethod>(
                        StringAccess::GetSubDictionary(settings, "non_trimmed_quadrature_rule_settings"),
                        "integration_method"
                    ))
                    == (IntegrationMethod::ggq_optimal)
                );
            }
        }

        TEST_CASE("SettingsTestSuite - Settings Condition Settings Wrong Type Test")
        {


            {  /// Key access
                auto p_condition_settings =
                    DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
                auto& cond_settings = *p_condition_settings;

                REQUIRE_THROWS_AS(cond_settings.CheckRequired(), std::exception);
                CHECK(!cond_settings.IsSet(ConditionSettings::condition_id));
                CHECK(!cond_settings.IsSet(ConditionSettings::input_filename));
                REQUIRE(!cond_settings.TryGetValue<std::string>(ConditionSettings::condition_type).has_value());
                REQUIRE(!cond_settings.TryGetValue<double>(ConditionSettings::modulus).has_value());
                REQUIRE(!cond_settings.TryGetValue<PointType>(ConditionSettings::direction).has_value());
                REQUIRE(!cond_settings.TryGetValue<PointType>(ConditionSettings::value).has_value());
                REQUIRE(!cond_settings.TryGetValue<double>(ConditionSettings::penalty_factor).has_value());

                if (!NOTDEBUG) {

                    REQUIRE_THROWS_AS(
                        cond_settings.IsSet(BackgroundGridSettings::lower_bound_uvw), std::exception
                    );  // Wrong Key

                    REQUIRE_THROWS_AS(
                        cond_settings.SetValue(BackgroundGridSettings::lower_bound_uvw, PointType({ 2.0, 3.0, 4.0 })),
                        std::exception
                    );  // Wrong Key type

                    cond_settings.SetValue(ConditionSettings::condition_id, 100u);
                    CHECK(cond_settings.IsSet(ConditionSettings::condition_id));
                    REQUIRE_THROWS_AS(
                        cond_settings.GetRequiredValue<IndexType>(MainSettings::echo_level), std::exception
                    );  // Wrong Key type

                    REQUIRE(!cond_settings.TryGetValue<std::string>(ConditionSettings::condition_type).has_value());
                    cond_settings.SetValue(ConditionSettings::condition_type, std::string("dummy"));
                    CHECK(
                        (cond_settings.TryGetValue<std::string>(ConditionSettings::condition_type)->get())
                        == (std::string("dummy"))
                    );
                    REQUIRE_THROWS_AS(
                        cond_settings.GetRequiredValue<IndexType>(MainSettings::echo_level), std::exception
                    );  // Wrong Key type

                    REQUIRE_THROWS_AS(cond_settings.CheckRequired(), std::exception);
                    REQUIRE(!cond_settings.TryGetValue<double>(ConditionSettings::modulus).has_value());
                    REQUIRE(!cond_settings.TryGetValue<PointType>(ConditionSettings::direction).has_value());
                    REQUIRE(!cond_settings.TryGetValue<PointType>(ConditionSettings::value).has_value());
                    REQUIRE(!cond_settings.TryGetValue<double>(ConditionSettings::penalty_factor).has_value());
                }
            }
            {  /// String access
                using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
                using StringAccess = DictionaryStringAccess<DictionaryType>;

                auto p_condition_settings =
                    DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
                auto& cond_settings = *p_condition_settings;

                REQUIRE_THROWS_AS(StringAccess::IsSet(cond_settings, "lower_bound_uvw"), std::exception);  // Wrong Key

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "condition_id"), std::exception
                );  // Not set
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(cond_settings, "lower_bound_uvw", PointType({ 2.0, 3.0, 4.0 })),
                    std::exception
                );  // Wrong Key type
                REQUIRE_THROWS_AS(
                    StringAccess::SetValue(cond_settings, "condition_id", 2.0), std::exception
                );  // Wrong Value type

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "condition_id"), std::exception
                );  // Not set
                StringAccess::SetValue(cond_settings, "condition_id", 100u);
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "echo_level"), std::exception
                );  // Wrong Key type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<double>(cond_settings, "condition_id"), std::exception
                );  // Wrong Type

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<std::string>(cond_settings, "condition_type"), std::exception
                );  // Not set
                StringAccess::SetValue(cond_settings, "condition_type", std::string("dummy"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "echo_level"), std::exception
                );  // Wrong Key type
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "condition_type"), std::exception
                );  // Wrong Type

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<std::string>(cond_settings, "input_filename"), std::exception
                );  // Not set
                StringAccess::SetValue(cond_settings, "input_filename", std::string("dummy_2"));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "input_filename"), std::exception
                );  // Wrong Type

                REQUIRE_THROWS_AS(StringAccess::GetValue<double>(cond_settings, "modulus"), std::exception);  // Not set
                StringAccess::SetValue(cond_settings, "modulus", 2.0);
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "modulus"), std::exception
                );  // Wrong Type

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(cond_settings, "direction"), std::exception
                );  // Not set
                StringAccess::SetValue(cond_settings, "direction", PointType({ 2.2, 3.0, 4.4 }));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "direction"), std::exception
                );  // Wrong Type

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<PointType>(cond_settings, "value"), std::exception
                );  // Not set
                StringAccess::SetValue(cond_settings, "value", PointType({ 2.1, 1.0, 2.4 }));
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "value"), std::exception
                );  // Wrong Type

                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<double>(cond_settings, "penalty_factor"), std::exception
                );  // Not set
                StringAccess::SetValue(cond_settings, "penalty_factor", 1e5);
                REQUIRE_THROWS_AS(
                    StringAccess::GetValue<IndexType>(cond_settings, "penalty_factor"), std::exception
                );  // Wrong Type
            }
        }


        TEST_CASE("SettingsTestSuite - Settings Condition Settings Default Values Test"){


            { /// Key access
              auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
        auto& cond_settings = *p_cond_settings;

        CHECK(!cond_settings.IsSet(ConditionSettings::condition_id));
        CHECK(!cond_settings.IsSet(ConditionSettings::condition_type));
        CHECK(!cond_settings.IsSet(ConditionSettings::input_filename));
        CHECK(!cond_settings.IsSet(ConditionSettings::modulus));
        CHECK(!cond_settings.IsSet(ConditionSettings::direction));
        CHECK(!cond_settings.IsSet(ConditionSettings::value));
        CHECK(!cond_settings.IsSet(ConditionSettings::penalty_factor));
    }
    {  /// String access
        using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
        using StringAccess = DictionaryStringAccess<DictionaryType>;

        auto p_cond_settings = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("ConditionSettings");
        auto& cond_settings = *p_cond_settings;

        CHECK(!StringAccess::IsSet(cond_settings, "condition_id"));
        CHECK(!StringAccess::IsSet(cond_settings, "condition_type"));
        CHECK(!StringAccess::IsSet(cond_settings, "input_filename"));
        CHECK(!StringAccess::IsSet(cond_settings, "modulus"));
        CHECK(!StringAccess::IsSet(cond_settings, "direction"));
        CHECK(!StringAccess::IsSet(cond_settings, "value"));
        CHECK(!StringAccess::IsSet(cond_settings, "penalty_factor"));
    }
};

TEST_CASE("SettingsTestSuite - Settings Condition Settings Customized Values Test"){


    { /// Key access
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
    p_cond_settings->SetValue(ConditionSettings::direction, PointType{ 2.1, 3.2, 4.7 });
    p_cond_settings->SetValue(ConditionSettings::value, PointType{ 2.2, 2.3, 1.4 });

    r_cond_setting_list.push_back(std::move(p_cond_settings));
}
{
    auto& r_cond_settings_list = p_settingss->GetList(MainSettings::conditions_settings_list);
    REQUIRE((r_cond_settings_list.size()) == (2));
    r_cond_settings_list[0u]->CheckRequired();

    CHECK(r_cond_settings_list[0u]->IsSet(ConditionSettings::condition_id));
    CHECK((r_cond_settings_list[0u]->GetRequiredValue<IndexType>(ConditionSettings::condition_id)) == (120u));

    CHECK(r_cond_settings_list[0u]->IsSet(ConditionSettings::condition_type));
    CHECK(
        (r_cond_settings_list[0u]->TryGetValue<std::string>(ConditionSettings::condition_type)->get())
        == (std::string("Hello"))
    );

    CHECK(r_cond_settings_list[0u]->IsSet(ConditionSettings::input_filename));
    CHECK(
        (r_cond_settings_list[0u]->GetRequiredValue<std::string>(ConditionSettings::input_filename))
        == (std::string("hallo"))
    );

    CHECK(r_cond_settings_list[0u]->IsSet(ConditionSettings::modulus));
    CHECK((r_cond_settings_list[0u]->TryGetValue<double>(ConditionSettings::modulus)->get()) == (200.0));

    CHECK(!r_cond_settings_list[0u]->IsSet(ConditionSettings::direction));
    REQUIRE(!r_cond_settings_list[0u]->TryGetValue<PointType>(ConditionSettings::direction).has_value());

    CHECK(!r_cond_settings_list[0u]->IsSet(ConditionSettings::value));
    REQUIRE(!r_cond_settings_list[0u]->TryGetValue<PointType>(ConditionSettings::value).has_value());

    CHECK(r_cond_settings_list[0u]->IsSet(ConditionSettings::penalty_factor));
    CHECK((r_cond_settings_list[0u]->TryGetValue<double>(ConditionSettings::penalty_factor)->get()) == (300.0));
}
{
    auto& r_cond_settings_list = p_settingss->GetList(MainSettings::conditions_settings_list);
    REQUIRE((r_cond_settings_list.size()) == (2));
    r_cond_settings_list[1u]->CheckRequired();

    CHECK(r_cond_settings_list[1u]->IsSet(ConditionSettings::condition_id));
    CHECK((r_cond_settings_list[1u]->GetRequiredValue<IndexType>(ConditionSettings::condition_id)) == (150u));

    CHECK(r_cond_settings_list[1u]->IsSet(ConditionSettings::condition_type));
    CHECK(
        (r_cond_settings_list[1u]->TryGetValue<std::string>(ConditionSettings::condition_type)->get())
        == (std::string("Hello2"))
    );

    CHECK(r_cond_settings_list[1u]->IsSet(ConditionSettings::input_filename));
    CHECK(
        (r_cond_settings_list[1u]->GetRequiredValue<std::string>(ConditionSettings::input_filename))
        == (std::string("hallo2"))
    );

    CHECK(!r_cond_settings_list[1u]->IsSet(ConditionSettings::modulus));
    REQUIRE(!r_cond_settings_list[1u]->TryGetValue<double>(ConditionSettings::modulus).has_value());

    CHECK(r_cond_settings_list[1u]->IsSet(ConditionSettings::direction));
    CheckPointNear(
        r_cond_settings_list[1u]->TryGetValue<PointType>(ConditionSettings::direction)->get(),
        PointType({ 2.1, 3.2, 4.7 }),
        1e-10
    );

    CHECK(r_cond_settings_list[1u]->IsSet(ConditionSettings::value));
    CheckPointNear(
        r_cond_settings_list[1u]->TryGetValue<PointType>(ConditionSettings::value)->get(),
        PointType({ 2.2, 2.3, 1.4 }),
        1e-10
    );

    CHECK(!r_cond_settings_list[1u]->IsSet(ConditionSettings::penalty_factor));
    REQUIRE(!r_cond_settings_list[1u]->TryGetValue<double>(ConditionSettings::penalty_factor).has_value());
}
}  // namespace queso
{  /// String access
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
        StringAccess::SetValue(*p_cond_settings, "direction", PointType{ 2.1, 3.2, 4.7 });
        StringAccess::SetValue(*p_cond_settings, "value", PointType{ 2.2, 2.3, 1.4 });

        r_cond_setting_list.push_back(std::move(p_cond_settings));
    }
    {
        auto& r_cond_settings_list = StringAccess::GetList(settings, "conditions_settings_list");
        REQUIRE((r_cond_settings_list.size()) == (2));

        CHECK(StringAccess::IsSet(*r_cond_settings_list[0u], "condition_id"));
        CHECK((StringAccess::GetValue<IndexType>(*r_cond_settings_list[0u], "condition_id")) == (120u));

        CHECK(StringAccess::IsSet(*r_cond_settings_list[0u], "condition_type"));
        CHECK(
            (StringAccess::GetValue<std::string>(*r_cond_settings_list[0u], "condition_type")) == (std::string("Hello"))
        );

        CHECK(StringAccess::IsSet(*r_cond_settings_list[0u], "input_filename"));
        CHECK(
            (StringAccess::GetValue<std::string>(*r_cond_settings_list[0u], "input_filename")) == (std::string("hallo"))
        );

        CHECK(StringAccess::IsSet(*r_cond_settings_list[0u], "modulus"));
        CHECK((StringAccess::GetValue<double>(*r_cond_settings_list[0u], "modulus")) == (200.0));

        CHECK(!StringAccess::IsSet(*r_cond_settings_list[0u], "direction"));
        REQUIRE_THROWS_AS(StringAccess::GetValue<PointType>(*r_cond_settings_list[0u], "direction"), std::exception);

        CHECK(!StringAccess::IsSet(*r_cond_settings_list[0u], "value"));
        REQUIRE_THROWS_AS(StringAccess::GetValue<PointType>(*r_cond_settings_list[0u], "value"), std::exception);

        CHECK(StringAccess::IsSet(*r_cond_settings_list[0u], "penalty_factor"));
        CHECK((StringAccess::GetValue<double>(*r_cond_settings_list[0u], "penalty_factor")) == (300.0));
    }
    {
        auto& r_cond_settings_list = StringAccess::GetList(settings, "conditions_settings_list");
        REQUIRE((r_cond_settings_list.size()) == (2));

        CHECK(StringAccess::IsSet(*r_cond_settings_list[1u], "condition_id"));
        CHECK((StringAccess::GetValue<IndexType>(*r_cond_settings_list[1u], "condition_id")) == (150u));

        CHECK(StringAccess::IsSet(*r_cond_settings_list[1u], "condition_type"));
        CHECK(
            (StringAccess::GetValue<std::string>(*r_cond_settings_list[1u], "condition_type"))
            == (std::string("Hello2"))
        );

        CHECK(StringAccess::IsSet(*r_cond_settings_list[1u], "input_filename"));
        CHECK(
            (StringAccess::GetValue<std::string>(*r_cond_settings_list[1u], "input_filename"))
            == (std::string("hallo2"))
        );

        CHECK(!StringAccess::IsSet(*r_cond_settings_list[1u], "modulus"));
        REQUIRE_THROWS_AS(StringAccess::GetValue<double>(*r_cond_settings_list[1u], "modulus"), std::exception);

        CHECK(StringAccess::IsSet(*r_cond_settings_list[1u], "direction"));
        CheckPointNear(
            StringAccess::GetValue<PointType>(*r_cond_settings_list[1u], "direction"),
            PointType({ 2.1, 3.2, 4.7 }),
            1e-10
        );

        CHECK(StringAccess::IsSet(*r_cond_settings_list[1u], "value"));
        CheckPointNear(
            StringAccess::GetValue<PointType>(*r_cond_settings_list[1u], "value"), PointType({ 2.2, 2.3, 1.4 }), 1e-10
        );

        CHECK(!StringAccess::IsSet(*r_cond_settings_list[1u], "penalty_factor"));
        REQUIRE_THROWS_AS(StringAccess::GetValue<double>(*r_cond_settings_list[1u], "penalty_factor"), std::exception);
    }
}
}
;


TEST_CASE("SettingsTestSuite - Settings Cast Ambiguous Types Test")
{


    using DictionaryType = Dictionary<queso::key::MainValuesTypeTag>;
    using StringAccess = DictionaryStringAccess<DictionaryType>;

    auto p_settingss = DictionaryFactory<queso::key::MainValuesTypeTag>::Create("Settings");
    auto& settings = *p_settingss;

    // Lets set it with wrong type
    StringAccess::SetValue(
        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual", 5u
    );
    {
        const double r_double = StringAccess::GetValue<double>(
            StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual"
        );
        CheckRelativeNear(r_double, 5.0, 1e-10);
    }
    // Lets set it with correct type
    StringAccess::SetValue(
        StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual", 7.7
    );
    {
        const double r_double = StringAccess::GetValue<double>(
            StringAccess::GetSubDictionary(settings, "trimmed_quadrature_rule_settings"), "moment_fitting_residual"
        );
        CheckRelativeNear(r_double, 7.7, 1e-10);
    }

    // Cast Vector3i to Vector3d
    REQUIRE_THROWS_AS(
        StringAccess::GetValue<PointType>(
            StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz"
        ),
        std::exception
    );  // Not set

    // Lets set it with wrong type
    StringAccess::SetValue(
        StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz", Vector3i({ 5, 3, 2 })
    );
    {
        const PointType& r_point = StringAccess::GetValue<PointType>(
            StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz"
        );
        CheckPointRelativeNear(r_point, PointType({ 5.0, 3.0, 2.0 }), 1e-10);
    }
    // Lets set it with correct type
    StringAccess::SetValue(
        StringAccess::GetSubDictionary(settings, "background_grid_settings"),
        "lower_bound_xyz",
        Vector3d({ 5.5, 2.2, 1.1 })
    );
    {
        const PointType& r_point = StringAccess::GetValue<PointType>(
            StringAccess::GetSubDictionary(settings, "background_grid_settings"), "lower_bound_xyz"
        );
        CheckPointRelativeNear(r_point, PointType({ 5.5, 2.2, 1.1 }), 1e-10);
    }
}

}  // TEST_SUITE

} // namespace queso::testing
