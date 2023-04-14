# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling CGAL Python Module.
cmake_minimum_required(VERSION 3.1...3.15)
project( TIBRA_Application )

# Set build type to Release
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")

# Include local directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../external_libraries)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Import CGAL and its components
if( USE_CGAL MATCHES ON )
  find_package( CGAL QUIET COMPONENTS  )

  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  if ( NOT CGAL_FOUND )
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
  endif()
  add_definitions( -DUSE_CGAL )
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /EHsc -W2 /std:c++17 /D _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -funroll-loops -std=c++17 ")
endif()

# OMP support
find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
  else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")
  endif()
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message(STATUS "OMP was not found")
  return()
endif()

message("-- CXX FLAGS: ${CMAKE_CXX_FLAGS}")

# Try to use python executable from env variable
if(DEFINED ENV{PYTHON_EXECUTABLE})
set(PYTHON_EXECUTABLE $ENV{PYTHON_EXECUTABLE})
endif(DEFINED ENV{PYTHON_EXECUTABLE})

################### PYBIND11
# Set TIBRA specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
include(pybind11Tools)

#include source directories
file(GLOB TIBRA_ApplicationSource
  ${CMAKE_CURRENT_SOURCE_DIR}/TIBRA_main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/quadrature/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/embedding/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/quadrature/integration_points_1d/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/containers/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/io/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/modeler/*.cpp
  # Add external libraries
  ${CMAKE_CURRENT_SOURCE_DIR}/../external_libraries/aabb_tree/*.cc
)

## If cgal support is enables
if( USE_CGAL MATCHES ON )
  file(GLOB TIBRA_CGALWrapperSource
    ${CMAKE_CURRENT_SOURCE_DIR}/cgal_wrapper/*.cpp
  )
endif()

## CGAL python interface sources
file(GLOB_RECURSE TIBRA_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/python/*.cpp)

if( USE_CGAL MATCHES ON )
  add_library(TIBRA_ApplicationCore SHARED ${TIBRA_ApplicationSource} ${TIBRA_CGALWrapperSource})
  target_link_libraries(TIBRA_ApplicationCore PUBLIC CGAL ${PYTHON_LIBRARIES})
else()
  add_library(TIBRA_ApplicationCore SHARED ${TIBRA_ApplicationSource} )
  target_link_libraries(TIBRA_ApplicationCore PUBLIC ${PYTHON_LIBRARIES})
endif()

set_target_properties(TIBRA_ApplicationCore PROPERTIES COMPILE_DEFINITIONS "TIBRA_Application=EXPORT,API")

pybind11_add_module(TIBRA_Application MODULE THIN_LTO ${TIBRA_PYTHON_INTERFACE})
target_link_libraries( TIBRA_Application PRIVATE TIBRA_ApplicationCore )
set_target_properties( TIBRA_Application PROPERTIES PREFIX "")

if(DEFINED ENV{BOOST_ROOT})
  set(BOOST_ROOT $ENV{BOOST_ROOT})
endif(DEFINED ENV{BOOST_ROOT})

option(Boost_USE_STATIC_LIBS "Use static versions of Boost libraries" OFF)
#include testing directories
enable_testing()
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories (${Boost_INCLUDE_DIRS})
file(GLOB TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
file(GLOB TEST_SRCS_2 ${CMAKE_CURRENT_SOURCE_DIR}/tests/thingi10k_tests/*.cpp)

#file(GLOB TEST_SRCS2 ${CMAKE_CURRENT_SOURCE_DIR}/tests/thingi10k_tests/*.cpp)
# Run compile cpp tests
# Extract the filename without an extension (NAME_WE)
get_filename_component(testName ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_tibra.cpp NAME_WE)
get_filename_component(testName_2 ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_tibra_big.cpp NAME_WE)
# Add compile target
add_executable(${testName} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_tibra.cpp ${TEST_SRCS})
add_executable(${testName_2} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_tibra_big.cpp ${TEST_SRCS_2})

# Link to Boost libraries and TIBRA
target_link_libraries(${testName} ${Boost_LIBRARIES} TIBRA_ApplicationCore)
target_link_libraries(${testName_2} ${Boost_LIBRARIES} TIBRA_ApplicationCore)

# Move to TestExecutables
set_target_properties(${testName} PROPERTIES
             RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables)
set_target_properties(${testName_2} PROPERTIES
             RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables)
add_test(NAME ${testName}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables/${testName} )

add_test(NAME ${testName_2}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables/${testName_2} )

#Setting the libs folder for the shared objects built in kratos
install(TARGETS TIBRA_Application DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../libs)
install(TARGETS TIBRA_ApplicationCore DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../libs)

# Install Python module
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/TIBRA_Application.py DESTINATION  ${CMAKE_CURRENT_SOURCE_DIR}/../TIBRA_PythonApplication RENAME __init__.py)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/PyTIBRA.py DESTINATION  ${CMAKE_CURRENT_SOURCE_DIR}/../TIBRA_PythonApplication)
