cmake_minimum_required(VERSION 3.21)

# QuESo python interface sources.
set(PythonSource
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_containers_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_dictionary_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_globals_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_io_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/QuESoPythonModule.cpp
)

if(QUESO_BUILD_TESTING)
    list(APPEND PythonSource ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_test_helpers_to_python.cpp)
endif()

# Find Python and the bundled pybind11 CMake package.
find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
set(pybind11_DIR "${CMAKE_SOURCE_DIR}/external_libraries/pybind11/share/cmake/pybind11")
find_package(pybind11 CONFIG REQUIRED)

# Create Python Module.
pybind11_add_module(_core ${PythonSource})
target_link_libraries(_core PRIVATE QuESo::QuESo)
if(QUESO_BUILD_TESTING)
    target_compile_definitions(_core PRIVATE QUESO_PYTHON_BUILD_TEST_HELPERS)
endif()
add_library(QuESo::PythonModule ALIAS _core)

set_target_properties(_core PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyqueso"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyqueso"
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyqueso"
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${CMAKE_BINARY_DIR}/pyqueso/__init__.py COPYONLY)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/model.py DESTINATION ${CMAKE_BINARY_DIR}/pyqueso)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts DESTINATION ${CMAKE_BINARY_DIR}/pyqueso)
file(COPY ${CMAKE_SOURCE_DIR}/kratos_interface DESTINATION ${CMAKE_BINARY_DIR}/pyqueso)

if(UNIX)
	# Set RPATH.
	set_target_properties(_core PROPERTIES
		BUILD_WITH_INSTALL_RPATH TRUE
		INSTALL_RPATH "$ORIGIN"  # relative to the Python module.
	)
endif()

# Copy libQuESo.so / QuESo.dll to the Python module directory.
if( UNIX )
	set(DllSource "${CMAKE_BINARY_DIR}/lib/libQuESo.so")
else()
	set(DllSource "${CMAKE_BINARY_DIR}/bin/QuESo.dll")
endif()
add_custom_command(TARGET _core POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
	"${DllSource}"
	"${CMAKE_BINARY_DIR}/pyqueso/"
) 

# Install bindings target.
install(TARGETS _core DESTINATION ${CMAKE_INSTALL_PREFIX}/pyqueso)

# Install Python module.
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py DESTINATION ${CMAKE_INSTALL_PREFIX}/pyqueso)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/model.py DESTINATION ${CMAKE_INSTALL_PREFIX}/pyqueso)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts DESTINATION ${CMAKE_INSTALL_PREFIX}/pyqueso)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/kratos_interface DESTINATION ${CMAKE_INSTALL_PREFIX}/pyqueso)
