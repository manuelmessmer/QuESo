# This is the CMake script for compiling QuESo Python Module.
cmake_minimum_required(VERSION 3.1)
project(QuESo_Application)

# Set build type to Release
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")

# Include local directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../external_libraries)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /EHsc -W2 /std:c++17 /D _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -funroll-loops -Wsuggest-override -std=c++17 ")
endif()

# OMP support
find_package(OpenMP)

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")
    endif()

    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
    message(STATUS "OMP was not found")
    return()
endif()

message("-- CXX FLAGS: ${CMAKE_CXX_FLAGS}")

# Try to use python executable from env variable
if(DEFINED ENV{PYTHON_EXECUTABLE})
    set(PYTHON_EXECUTABLE $ENV{PYTHON_EXECUTABLE})
endif(DEFINED ENV{PYTHON_EXECUTABLE})

# ################## PYBIND11
# Set QuESo specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
include(pybind11Tools)

# include source directories
file(GLOB QuESo_ApplicationSource
    ${CMAKE_CURRENT_SOURCE_DIR}/QuESo_main.cpp
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

# # QuESo python interface sources
file(GLOB_RECURSE QuESo_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/python/*.cpp)

add_library(QuESo_ApplicationCore SHARED ${QuESo_ApplicationSource})
target_link_libraries(QuESo_ApplicationCore PUBLIC ${PYTHON_LIBRARIES})
set_target_properties(QuESo_ApplicationCore PROPERTIES COMPILE_DEFINITIONS "QuESo_Application=EXPORT,API")

pybind11_add_module(QuESo_Application MODULE THIN_LTO ${QuESo_PYTHON_INTERFACE})
target_link_libraries(QuESo_Application PRIVATE QuESo_ApplicationCore)
set_target_properties(QuESo_Application PROPERTIES PREFIX "")

if(${QUESO_BUILD_TESTING} MATCHES ON)
    if(DEFINED ENV{BOOST_ROOT})
        set(BOOST_ROOT $ENV{BOOST_ROOT})
    endif(DEFINED ENV{BOOST_ROOT})

    option(Boost_USE_STATIC_LIBS "Use static versions of Boost libraries" OFF)

    # include testing directories
    enable_testing()
    find_package(Boost COMPONENTS unit_test_framework REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
    file(GLOB TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    file(GLOB TEST_THINGI_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/tests/thingi10k_tests/*.cpp)

    # file(GLOB TEST_SRCS2 ${CMAKE_CURRENT_SOURCE_DIR}/tests/thingi10k_tests/*.cpp)
    # Run compile cpp tests
    # Extract the filename without an extension (NAME_WE)
    get_filename_component(testName ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_queso.cpp NAME_WE)
    get_filename_component(testThingiName ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_queso_thingi10k.cpp NAME_WE)

    # Add compile target
    add_executable(${testName} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_queso.cpp ${TEST_SRCS})
    add_executable(${testThingiName} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_queso_thingi10k.cpp ${TEST_THINGI_SRCS})

    # Link to Boost libraries and QuESo
    target_link_libraries(${testName} ${Boost_LIBRARIES} QuESo_ApplicationCore)
    target_link_libraries(${testThingiName} ${Boost_LIBRARIES} QuESo_ApplicationCore)

    # Move to TestExecutables
    set_target_properties(${testName} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables)
    set_target_properties(${testThingiName} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables)

    add_test(NAME ${testName}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables/${testName})

    add_test(NAME ${testThingiName}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables/${testThingiName})
endif(${QUESO_BUILD_TESTING} MATCHES ON)

# Setting the libs folder for the shared objects built in kratos
install(TARGETS QuESo_Application DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../libs)
install(TARGETS QuESo_ApplicationCore DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../libs)

# Install Python module
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/QuESo_Application.py DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../QuESo_PythonApplication RENAME __init__.py)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/PyQuESo.py DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../QuESo_PythonApplication)
