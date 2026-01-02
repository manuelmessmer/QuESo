cmake_minimum_required(VERSION 3.21)

# Set the CPP standard.
if (NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
endif()

# Set the project name and language.
project(
    QuESo
    VERSION 0.0.1
    DESCRIPTION "Qadrature for Embedded Solids"
    LANGUAGES CXX C
)

# Create target interfaces.
add_library(QuESoWarnings INTERFACE)
add_library(QuESoOptions INTERFACE)

include(cmake/ProjectOptions.cmake)
queso_setup_options()

include(cmake/CompilerWarnings.cmake)
queso_set_project_warnings(
    QuESoWarnings
    ${QUESO_WARNINGS_AS_ERRORS}
    ""
    ""
    ""
    ""
)

target_compile_features(QuESoOptions INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_library(QuESo::QuESoOptions ALIAS QuESoOptions)
add_library(QuESo::QuESoWarnings ALIAS QuESoWarnings)

include(cmake/Common.cmake)

# Create externals target.
add_library(Externals STATIC)
target_sources(Externals
    PRIVATE
        ${CMAKE_SOURCE_DIR}/external_libraries/aabb_tree/AABB_base.cc
)
target_include_directories(Externals
    PUBLIC
        ${CMAKE_SOURCE_DIR}/external_libraries
)
set_target_properties(Externals PROPERTIES POSITION_INDEPENDENT_CODE ON)
queso_set_common_target_properties(Externals)
add_library(QuESo::Externals ALIAS Externals)

# Add main QuESo source.
add_subdirectory(queso)
