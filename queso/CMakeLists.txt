# This is the CMake script for compiling QuESo Python Module.
cmake_minimum_required(VERSION 3.16)
project(QuESo_Application)

if(NOT DEFINED ENV{QUESO_BUILD_TYPE})
    message(WARNING "QUESO_BUILD_TYPE environment variable is not set! Defaulting to Release.")
    set(CMAKE_BUILD_TYPE Release)
else()
    set(CMAKE_BUILD_TYPE $ENV{QUESO_BUILD_TYPE})
endif()
message("-- CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include local directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../external_libraries)

# Compiler flags
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    add_compile_options(/EHsc /W1 /std:c++17)
else()
    add_compile_options(-funroll-loops -Wall -Wsuggest-override -std=c++17)
endif()

# OMP support
find_package(OpenMP)
if(OPENMP_FOUND)
    message(STATUS "Using OpenMP for shared memory parallelization")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(STATUS "OpenMP not found, disabling parallelization.")
endif()

# Print compiler flags
message("-- CXX FLAGS: ${CMAKE_CXX_FLAGS}")
message("-- C FLAGS: ${CMAKE_C_FLAGS}")

# Use Python executable from environment if available
if(DEFINED ENV{PYTHON_EXECUTABLE})
    set(PYTHON_EXECUTABLE $ENV{PYTHON_EXECUTABLE})
endif(DEFINED ENV{PYTHON_EXECUTABLE})

# Include source directories
file(GLOB QuESo_ApplicationSource
    ${CMAKE_CURRENT_SOURCE_DIR}/embedded_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/quadrature/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/embedding/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/quadrature/integration_points_1d/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/containers/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/solvers/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/modeler/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/includes/*.cpp

    # Add external libraries
    ${CMAKE_CURRENT_SOURCE_DIR}/../external_libraries/aabb_tree/*.cc
)

# QuESo python interface sources
file(GLOB QuESo_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/python/*.cpp)

# Create Core Library
add_library(QuESo_ApplicationCore SHARED ${QuESo_ApplicationSource})
target_compile_definitions(QuESo_ApplicationCore PRIVATE QUESO_BUILD_TESTING=0)
target_link_libraries(QuESo_ApplicationCore PUBLIC ${PYTHON_LIBRARIES})
set_target_properties(QuESo_ApplicationCore PROPERTIES COMPILE_DEFINITIONS "QuESo_Application=EXPORT,API")

# Find pybind11
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
include(pybind11Tools)

# Create Python Module
pybind11_add_module(QuESo_Application MODULE THIN_LTO ${QuESo_PYTHON_INTERFACE})
target_link_libraries(QuESo_Application PRIVATE QuESo_ApplicationCore)
set_target_properties(QuESo_Application PROPERTIES PREFIX "")

# === TESTING ===
option(QUESO_BUILD_TESTING "Enable unit testing framework" OFF)
if(QUESO_BUILD_TESTING)
    message("-- QuESo C++ unit testing framework: enabled.")
else()
    message("-- QuESo C++ unit testing framework: disabled.")
endif()

if(QUESO_BUILD_TESTING)
    # Find Boost
    if(DEFINED ENV{BOOST_ROOT})
        set(BOOST_ROOT $ENV{BOOST_ROOT})
    endif(DEFINED ENV{BOOST_ROOT})
    option(Boost_USE_STATIC_LIBS "Use static versions of Boost libraries" OFF)
    find_package(Boost REQUIRED COMPONENTS unit_test_framework)

    # Add test sources
    file(GLOB TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    file(GLOB TEST_THINGI_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/tests/thingi10k_tests/*.cpp)

    # Extract test names
    get_filename_component(testName ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_queso.cpp NAME_WE)
    get_filename_component(testThingiName ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_queso_thingi10k.cpp NAME_WE)

    # Add test executables
    add_executable(${testName} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_queso.cpp ${TEST_SRCS})
    add_executable(${testThingiName} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_queso_thingi10k.cpp ${TEST_THINGI_SRCS})

    # Target compile definitions
    target_compile_definitions(${testName} PRIVATE QUESO_BUILD_TESTING=1)
    target_compile_definitions(${testThingiName} PRIVATE QUESO_BUILD_TESTING=1)

    # Link to Boost libraries and QuESo
    target_link_libraries(${testName} Boost::unit_test_framework QuESo_ApplicationCore)
    target_link_libraries(${testThingiName} Boost::unit_test_framework QuESo_ApplicationCore)

    # Move to TestExecutables
    set(TEST_EXECUTABLE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables)
    set_target_properties(${testName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXECUTABLE_DIR})
    set_target_properties(${testThingiName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_EXECUTABLE_DIR})

    # Register tests
    add_test(NAME ${testName} COMMAND ${TEST_EXECUTABLE_DIR}/${testName})
    add_test(NAME ${testThingiName} COMMAND ${TEST_EXECUTABLE_DIR}/${testThingiName})
endif()

# Setting the libs folder for the shared objects built in kratos
install(TARGETS QuESo_Application DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../libs)
install(TARGETS QuESo_ApplicationCore DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../libs)

# Install Python module
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/QuESo_Application.py DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../QuESo_PythonApplication RENAME __init__.py)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/PyQuESo.py DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../QuESo_PythonApplication)
