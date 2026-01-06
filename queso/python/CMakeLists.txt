cmake_minimum_required(VERSION 3.21)

# QuESo python interface sources.
set(PythonSource
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_containers_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_dictionary_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_globals_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_io_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/add_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/QuESoPythonModule.cpp
)

# Find pybind11.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/pybind11Tools.cmake)

# Create Python Module.
pybind11_add_module(_QuESoPythonModule ${PythonSource})
target_link_libraries(_QuESoPythonModule PRIVATE QuESo::QuESo)
add_library(QuESo::QuESoPythonModule ALIAS _QuESoPythonModule)

set_target_properties(_QuESoPythonModule PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/QuESoPythonModule"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/QuESoPythonModule"
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/QuESoPythonModule"
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/QuESoPythonModule.py ${CMAKE_BINARY_DIR}/QuESoPythonModule/__init__.py COPYONLY)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/PyQuESo.py DESTINATION ${CMAKE_BINARY_DIR}/QuESoPythonModule)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts DESTINATION ${CMAKE_BINARY_DIR}/QuESoPythonModule)
file(COPY ${CMAKE_SOURCE_DIR}/kratos_interface DESTINATION ${CMAKE_BINARY_DIR}/QuESoPythonModule)

if(UNIX)
	# Set RPATH.
	set_target_properties(_QuESoPythonModule PROPERTIES
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
add_custom_command(TARGET _QuESoPythonModule POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
	"${DllSource}"
	"${CMAKE_BINARY_DIR}/QuESoPythonModule/"
) 

# Install bindings target.
install(TARGETS _QuESoPythonModule DESTINATION ${CMAKE_INSTALL_PREFIX}/QuESoPythonModule)

# Install Python module.
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/QuESoPythonModule.py DESTINATION ${CMAKE_INSTALL_PREFIX}/QuESoPythonModule RENAME __init__.py)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/PyQuESo.py DESTINATION ${CMAKE_INSTALL_PREFIX}/QuESoPythonModule)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts DESTINATION ${CMAKE_INSTALL_PREFIX}/QuESoPythonModule)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/kratos_interface DESTINATION ${CMAKE_INSTALL_PREFIX}/QuESoPythonModule)

