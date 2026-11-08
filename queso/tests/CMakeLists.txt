cmake_minimum_required(VERSION 3.21)

# Add doctest
set(DoctestDirectory "${PROJECT_SOURCE_DIR}/external_libraries/doctest")
add_library(QuESoDoctest INTERFACE)
add_library(doctest::doctest ALIAS QuESoDoctest)
target_include_directories(QuESoDoctest SYSTEM INTERFACE "${DoctestDirectory}")
include("${DoctestDirectory}/scripts/cmake/doctest.cmake")

# Add test sources
file(GLOB TEST_SRCS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cpp/test_*.cpp)
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
target_sources(${testName} PRIVATE cpp/helpers/test_data.cpp)
target_sources(${testThingiName} PRIVATE cpp/helpers/thingi10k_options.cpp)

# Target compile definitions
target_compile_definitions(${testName} PRIVATE QUESO_BUILD_TESTING=1)
target_compile_definitions(${testThingiName} PRIVATE QUESO_BUILD_TESTING=1)
foreach(TestTarget IN ITEMS ${testName} ${testThingiName})
    target_compile_definitions(${TestTarget} PRIVATE QUESO_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
    target_link_libraries(${TestTarget} PRIVATE doctest::doctest QuESo::QuESo)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${TestTarget} PRIVATE OpenMP::OpenMP_CXX)
    endif()
endforeach()


# Move to TestExecutableDir
set(TestExecutableDir ${CMAKE_CURRENT_BINARY_DIR}/bin)
set_target_properties(${testName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TestExecutableDir})
set_target_properties(${testThingiName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TestExecutableDir})

# Discovery launches the executable during the build, before run_tests.cmd sets PATH.
if(WIN32)
    foreach(TestTarget IN ITEMS ${testName} ${testThingiName})
        add_custom_command(TARGET ${TestTarget} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:QuESo>" "$<TARGET_FILE_DIR:${TestTarget}>"
            VERBATIM)
    endforeach()
endif()

# Discovery only lists cases; test execution inherits the caller's QUESO_DEBUG setting.
# Thingi10K requires explicit corpus arguments and is intentionally not auto-registered.
doctest_discover_tests(${testName}
    TEST_PREFIX "queso/"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    EXTRA_ARGS "--data-dir=${CMAKE_CURRENT_SOURCE_DIR}/data")

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
