cmake_minimum_required(VERSION 3.21)

# Use NEW policy for case-sensitive package root variables
if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

# Use NEW policy for BoostConfig.cmake
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

# Find Boost
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

# Add test sources
file(GLOB TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/cpp/test_*.cpp)
set(TEST_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/cpp/test_queso.cpp)
set(TEST_THINGI_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/cpp/test_queso_thingi10k.cpp)
set(TEST_THINGI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/test_thingi10k.cpp)
list(REMOVE_ITEM TEST_SRCS ${TEST_MAIN} ${TEST_THINGI_MAIN} ${TEST_THINGI_SOURCE})

# Extract test names
get_filename_component(testName ${TEST_MAIN} NAME_WE)
get_filename_component(testThingiName ${TEST_THINGI_MAIN} NAME_WE)

# Add test executables
add_executable(${testName} ${TEST_MAIN} ${TEST_SRCS})
add_executable(${testThingiName} ${TEST_THINGI_MAIN} ${TEST_THINGI_SOURCE})

# Target compile definitions
target_compile_definitions(${testName} PRIVATE QUESO_BUILD_TESTING=1)
target_compile_definitions(${testThingiName} PRIVATE QUESO_BUILD_TESTING=1)

# Link to Boost libraries and QuESo
target_link_libraries(${testName} PRIVATE Boost::unit_test_framework QuESo::QuESo)
target_link_libraries(${testThingiName} PRIVATE Boost::unit_test_framework QuESo::QuESo OpenMP::OpenMP_CXX)
target_compile_definitions(${testName} PRIVATE BOOST_TEST_DYN_LINK)
target_compile_definitions(${testThingiName} PRIVATE BOOST_TEST_DYN_LINK)

# Move to TestExecutableDir
set(TestExecutableDir ${CMAKE_CURRENT_BINARY_DIR}/bin)
set_target_properties(${testName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TestExecutableDir})
set_target_properties(${testThingiName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TestExecutableDir})

# Configure shell script
set(PYTHON_TEST_SOURCE_DIR_IN ${CMAKE_CURRENT_SOURCE_DIR}/python)
set(TEST_DATA_DIR_IN ${CMAKE_CURRENT_SOURCE_DIR}/data)
set(CPP_TEST_EXE_DIR_IN ${TestExecutableDir})
set(THINGI10K_ID_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data/thingi10k)

if( CMAKE_HOST_WIN32 )
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cmd.in ${CMAKE_BINARY_DIR}/bin/run_tests.cmd @ONLY)
else()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_tests.in ${CMAKE_BINARY_DIR}/bin/run_tests @ONLY)
endif()
