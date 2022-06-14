
# This is the CMake script for compiling TrIGA Python Module.
cmake_minimum_required(VERSION 3.1...3.15)
project( TrIGA_Application )

# Set build type to Release
set(CMAKE_BUILD_TYPE RELEASE)
# Include local directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/includes)

# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

# include vtk
find_package(VTK REQUIRED)
if(VTK_FOUND)
  if(VTK_USE_FILE)
    include(${VTK_USE_FILE})
  endif()
  if ("${VTK_VERSION_MAJOR}" GREATER "5" OR VTK_VERSION VERSION_GREATER 5)
    message("found VTK. Version:" ${VTK_VERSION}. VTK_DIR: ${VTK_DIR})
  else()
    message(STATUS "VTK version 6.0 or greater is required")
  endif()
else()
  message(STATUS "VTK was not found")
endif()

# OMP support
find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message(STATUS "OMP was not found")
  return()
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/includes/vtkbool)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/includes/vtkbool)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/includes/vtkbool/libs/vp)

################### PYBIND11
# Set kratos specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
include(pybind11Tools)

#include source directories
file(GLOB TrIGA_ApplicationSource
  ${CMAKE_CURRENT_SOURCE_DIR}/includes/isotropic_remeshing/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/geometries/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/TrIGA_main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/io/*.cpp
)
## TrIGA python interface sources
file(GLOB_RECURSE TrIGA_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/python/*.cpp)

add_library(TrIGA_ApplicationCore SHARED ${TrIGA_ApplicationSource})
target_link_libraries(TrIGA_ApplicationCore PUBLIC ${VTK_LIBRARIES} vtkbool)
set_target_properties(TrIGA_ApplicationCore PROPERTIES COMPILE_DEFINITIONS "TrIGA_Application=EXPORT,API")

pybind11_add_module(TrIGA_Application MODULE THIN_LTO ${TrIGA_PYTHON_INTERFACE})
target_link_libraries( TrIGA_Application PRIVATE TrIGA_ApplicationCore )
set_target_properties( TrIGA_Application PROPERTIES PREFIX "")

# Include testing directories
enable_testing()
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories (${Boost_INCLUDE_DIRS})
file(GLOB TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)

#Run through each source
foreach(testSrc ${TEST_SRCS})
        #Extract the filename without an extension (NAME_WE)
        get_filename_component(testName ${testSrc} NAME_WE)
        #Add compile target
        add_executable(${testName} ${testSrc})
        #link to Boost libraries AND your targets and dependencies
        target_link_libraries(${testName} ${Boost_LIBRARIES} TrIGA_ApplicationCore)
        #I like to move testing binaries into a testBin directory
        set_target_properties(${testName} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables)
        #Finally add it to test execution -
        #Notice the WORKING_DIRECTORY and COMMAND
        add_test(NAME ${testName}
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../TestExecutables/${testName} )
endforeach(testSrc)

# Setting the libs folder for the shared objects built in kratos
install(TARGETS TrIGA_Application DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../libs)
install(TARGETS TrIGA_ApplicationCore DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../libs)
install(TARGETS vtkbool DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../libs)

# Install Python module
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/TrIGA_Application.py DESTINATION  ${CMAKE_CURRENT_SOURCE_DIR}/../TrIGA_PythonApplication RENAME __init__.py)

