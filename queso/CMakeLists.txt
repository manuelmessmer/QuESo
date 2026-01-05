cmake_minimum_required(VERSION 3.21)

set(EMBEDDING_SOURCES
    embedding/aabb_primitive.cpp
    embedding/aabb_tree.cpp
    embedding/brep_operator.cpp
    embedding/clipper.cpp
    embedding/flood_fill.cpp
    embedding/geometry_query.cpp
    embedding/octree.cpp
    embedding/polygon.cpp
    embedding/ray_aabb_primitive.cpp
    embedding/trimmed_domain_on_plane.cpp
    embedding/trimmed_domain.cpp
)
set(UTILITIES_SOURCES
    utilities/mesh_utilities.cpp
)
set(QUADRATURE_SOURCES
    quadrature/integration_points_1d/gauss_legendre_integration_points.cpp
    quadrature/integration_points_1d/precomputed_ggq_rules_optimal.cpp
    quadrature/integration_points_1d/precomputed_ggq_rules_reduced1.cpp
    quadrature/integration_points_1d/precomputed_ggq_rules_reduced2.cpp
    quadrature/integration_points_1d/integration_points_factory_1d.cpp
)
set(SOLVERS_SOURCES
    solvers/nnls.cpp
)
set(IO_SOURCES
    io/io_utilities.cpp
)

if( WIN32 )
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

add_library(QuESo SHARED
	${EMBEDDING_SOURCES}
	${UTILITIES_SOURCES}
	${QUADRATURE_SOURCES}
	${SOLVERS_SOURCES}
	${IO_SOURCES}
	embedded_model.cpp
)
add_library(QuESo::QuESo ALIAS QuESo)

target_compile_definitions(QuESo PRIVATE QUESO_BUILD_TESTING=0)

queso_set_common_target_properties(QuESo)

# Add include directories.
target_include_directories(QuESo
    PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/queso
)

# Link necessary libraries.
target_link_libraries(QuESo
    PRIVATE
        QuESo::QuESoWarnings
    PUBLIC
        QuESo::Externals
        QuESo::QuESoOptions
)
# Add OMP support.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "Using OpenMP for shared memory parallelization.")
    target_link_libraries(QuESo PRIVATE OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP not found, disabling parallelization.")
endif()

# Install target.
install(TARGETS QuESo DESTINATION ${CMAKE_INSTALL_PREFIX}/QuESoPythonModule)

# Create Python module.
add_subdirectory(python)

# Create test module.
if( QUESO_BUILD_TESTING )
    add_subdirectory(tests)
endif()
