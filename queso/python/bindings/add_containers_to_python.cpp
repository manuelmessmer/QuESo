//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

// STL includes
#include <iterator>
#include <ranges>
#include <span>
#include <utility>

// Project includes
#include "queso/embedded_model.h"
#include "queso/python/bindings/add_containers_to_python.h"
#include "queso/python/bindings/define_python.hpp"
#include "queso/utilities/triangle_utilities.hpp"

// Note: PYBIND11_MAKE_OPAQUE must live at file scope.

using IntegrationPointType = queso::IntegrationPoint;
using BoundaryIntegrationPointType = queso::BoundaryIntegrationPoint;
using BackgroundGridType = queso::BackgroundGrid<IntegrationPointType, BoundaryIntegrationPointType>;

// Bound element types
using ElementViewType = BackgroundGridType::ElementViewType;

// Bound condition types
using ConditionType = queso::Condition<ElementViewType>;
using ConditionPtrType = queso::Unique<ConditionType>;

// Bound condition segment types
using ConditionSegmentType = queso::ConditionSegment<ElementViewType>;

// Bound integration point types
using IntegrationPointVectorType = ElementViewType::IntegrationPointVectorType;
PYBIND11_MAKE_OPAQUE(IntegrationPointVectorType)

using BoundaryIpVectorType = std::vector<BoundaryIntegrationPointType>;
PYBIND11_MAKE_OPAQUE(BoundaryIpVectorType);

namespace queso::python {

using MainDictionaryHolderType = UniqueHolder<MainDictionaryType>;

namespace py = pybind11;

namespace {

    using PythonTriangleProxy = TriangleProxy<WithNormals>;

    template<class TRange>
    struct PythonRange
    {
        TRange range;
    };

    template<class TRange, py::return_value_policy TPolicy = py::return_value_policy::reference_internal>
    void BindRange(py::module& rModule, const char* pName)
    {
        auto range_binder = py::class_<PythonRange<TRange>>(rModule, pName);

        if constexpr (std::ranges::sized_range<TRange>) {
            range_binder.def("__len__", [](const PythonRange<TRange>& self) { return std::ranges::size(self.range); });
        }

        if constexpr (std::ranges::random_access_range<TRange> && std::ranges::sized_range<TRange>) {
            range_binder.def(
                "__getitem__",
                [](PythonRange<TRange>& self, py::ssize_t Index) -> decltype(auto) {
                    const auto size = static_cast<py::ssize_t>(std::ranges::size(self.range));
                    if (Index < 0) { Index += size; }
                    if (Index < 0 || Index >= size) { throw py::index_error(); }

                    using DifferenceType = std::iter_difference_t<decltype(std::ranges::begin(self.range))>;
                    return *std::ranges::next(std::ranges::begin(self.range), static_cast<DifferenceType>(Index));
                },
                py::return_value_policy::reference_internal
            );
        }

        range_binder.def(
            "__iter__",
            [](PythonRange<TRange>& self) {
                return py::make_iterator<TPolicy>(std::ranges::begin(self.range), std::ranges::end(self.range));
            },
            py::keep_alive<0, 1>()
        );
    }

    inline std::tuple<double, double, double> ToTuple3(const Vector3d& rV)
    { return { rV[0], rV[1], rV[2] }; }

    auto VerticesAsTuples(const TriangleMesh& rMesh)
    {
        return rMesh.Vertices() | std::views::transform([](const Vector3d& rVertex) { return ToTuple3(rVertex); });
    }

}  // namespace

using namespace pybind11::literals;

void AddContainersToPython(pybind11::module& m, pybind11::module& rMeshModule)
{
    auto InternalModule = m.def_submodule("_internal");
    auto MeshInternalModule = rMeshModule.def_submodule("_internal");

    using TriangleRangeType = decltype(std::declval<const TriangleMesh&>().Triangles<WithNormals>());
    using PythonTriangleRange = PythonRange<TriangleRangeType>;
    using VertexRangeType = decltype(VerticesAsTuples(std::declval<const TriangleMesh&>()));
    using PythonVertexRange = PythonRange<VertexRangeType>;
    BindRange<TriangleRangeType>(MeshInternalModule, "_TriangleRange");
    BindRange<VertexRangeType>(MeshInternalModule, "_VertexRange");

    py::class_<PythonTriangleProxy>(rMeshModule, "Triangle", "Read-only view of a triangle mesh face.")
        .def_property_readonly(
            "p1",
            [](const PythonTriangleProxy& self) { return ToTuple3(Vector3d{ self.P1[0], self.P1[1], self.P1[2] }); },
            "First vertex coordinates."
        )
        .def_property_readonly(
            "p2",
            [](const PythonTriangleProxy& self) { return ToTuple3(Vector3d{ self.P2[0], self.P2[1], self.P2[2] }); },
            "Second vertex coordinates."
        )
        .def_property_readonly(
            "p3",
            [](const PythonTriangleProxy& self) { return ToTuple3(Vector3d{ self.P3[0], self.P3[1], self.P3[2] }); },
            "Third vertex coordinates."
        )
        .def_property_readonly(
            "normal",
            [](const PythonTriangleProxy& self) {
                return ToTuple3(Vector3d{ self.Normal[0], self.Normal[1], self.Normal[2] });
            },
            "Unit normal vector."
        )
        .def(
            "area",
            [](const PythonTriangleProxy& self) { return TriangleUtilities::Area(self); },
            "Return the triangle area."
        )
        .def(
            "center",
            [](const PythonTriangleProxy& self) { return ToTuple3(TriangleUtilities::Center(self)); },
            "Return the triangle center."
        )
        .def(
            "aspect_ratio",
            [](const PythonTriangleProxy& self) { return TriangleUtilities::AspectRatio(self); },
            "Return the triangle aspect ratio."
        )
        .def(
            "integration_points_global",
            [](const PythonTriangleProxy& self, IndexType Method) {
                return TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPointType>(self, Method);
            },
            "method"_a,
            py::return_value_policy::move,
            "Construct global boundary integration points for this triangle."
        );

    py::class_<IntegrationPointType>(m, "IntegrationPoint", "Weighted point in three-dimensional space.")
        .def(py::init<double, double, double, double>(), "x"_a, "y"_a, "z"_a, "weight"_a, "Create a weighted point.")
        .def_property_readonly(
            "x", [](const IntegrationPointType& self) { return self[0]; }, "X coordinate."
        )
        .def_property_readonly(
            "y", [](const IntegrationPointType& self) { return self[1]; }, "Y coordinate."
        )
        .def_property_readonly(
            "z", [](const IntegrationPointType& self) { return self[2]; }, "Z coordinate."
        )
        .def_property_readonly("weight", &IntegrationPointType::Weight, "Integration weight.")
        .def(
            "__getitem__",
            [](const IntegrationPointType& self, int i) {
                if (i < 0) i += 3;
                if (i < 0 || i >= 3) throw py::index_error();
                return self[static_cast<IndexType>(i)];
            }
        )
        .def("__len__", [](const IntegrationPointType&) { return 3; })
        .def(
            "__iter__",
            [](const IntegrationPointType& self) {
                return py::make_iterator(self.Point().begin(), self.Point().end());
            },
            py::keep_alive<0, 1>()
        )
        .def("__repr__", [](const IntegrationPointType& self) {
            std::ostringstream oss;
            oss << "IntegrationPoint(" << self[0] << ", " << self[1] << ", " << self[2] << ", "
                << "weight=" << self.Weight() << ")";
            return oss.str();
        });

    py::bind_vector<IntegrationPointVectorType>(
        m, "IntegrationPointVector", "Mutable collection of integration points."
    );

    py::class_<BoundaryIntegrationPointType>(m, "BoundaryIntegrationPoint", "Weighted point with an outward normal.")
        .def(
            py::init<double, double, double, double, const std::array<double, 3>&>(),
            "x"_a,
            "y"_a,
            "z"_a,
            "weight"_a,
            "normal"_a
        )
        .def_property_readonly(
            "x", [](const BoundaryIntegrationPointType& self) { return self[0]; }, "X coordinate."
        )
        .def_property_readonly(
            "y", [](const BoundaryIntegrationPointType& self) { return self[1]; }, "Y coordinate."
        )
        .def_property_readonly(
            "z", [](const BoundaryIntegrationPointType& self) { return self[2]; }, "Z coordinate."
        )
        .def_property_readonly("weight", &BoundaryIntegrationPointType::Weight, "Integration weight.")
        .def_property_readonly(
            "normal",
            &BoundaryIntegrationPointType::Normal,
            py::return_value_policy::reference_internal,
            "Outward normal vector."
        )
        .def(
            "__getitem__",
            [](const BoundaryIntegrationPointType& self, int i) {
                if (i < 0) i += 3;
                if (i < 0 || i >= 3) throw py::index_error();
                return self[static_cast<IndexType>(i)];
            }
        )
        .def("__len__", [](const BoundaryIntegrationPointType&) { return 3; })
        .def(
            "__iter__",
            [](const BoundaryIntegrationPointType& self) {
                return py::make_iterator(self.Point().begin(), self.Point().end());
            },
            py::keep_alive<0, 1>()
        )
        .def("__repr__", [](const BoundaryIntegrationPointType& self) {
            std::ostringstream oss;
            oss << "BoundaryIntegrationPoint(" << self[0] << ", " << self[1] << ", " << self[2]
                << ", weight=" << self.Weight() << ", normal=" << self.Normal() << ")";
            return oss.str();
        });

    py::bind_vector<BoundaryIpVectorType>(
        m, "BoundaryIntegrationPointVector", "Mutable collection of boundary integration points."
    );

    py::class_<TriangleMesh, Unique<TriangleMesh>>(rMeshModule, "TriangleMesh", "Mutable triangle mesh.")
        .def(py::init<>())
        .def_property_readonly("num_triangles", &TriangleMesh::NumOfTriangles, "Number of triangles in the mesh.")
        .def_property_readonly("num_vertices", &TriangleMesh::NumOfVertices, "Number of vertices in the mesh.")
        .def(
            "vertices",
            [](const TriangleMesh& rMesh) { return PythonVertexRange{ VerticesAsTuples(rMesh) }; },
            py::keep_alive<0, 1>(),
            "Iterate over mesh vertices."
        )
        .def(
            "triangles",
            [](const TriangleMesh& rMesh) { return PythonTriangleRange{ rMesh.Triangles<WithNormals>() }; },
            py::keep_alive<0, 1>(),
            "Iterate over mesh triangles."
        )
        .def("reserve", &TriangleMesh::Reserve, "num_triangles"_a, "Reserve storage for triangles.")
        .def(
            "add_vertex",
            [](TriangleMesh& self, const std::array<double, 3>& rVertex) { return self.AddVertex(PointType(rVertex)); },
            "vertex"_a,
            "Append a vertex and return its index."
        )
        .def(
            "add_triangle",
            [](TriangleMesh& self, const std::array<IndexType, 3>& rTriangle) {
                const Vector3i triangle_indices(rTriangle);
                const auto& vertices = self.Vertices();
                const TriangleProxy<WithoutNormals> triangle{
                    std::span<const double, 3>(vertices[triangle_indices[0]].data(), 3),
                    std::span<const double, 3>(vertices[triangle_indices[1]].data(), 3),
                    std::span<const double, 3>(vertices[triangle_indices[2]].data(), 3)
                };
                self.AddTriangle(triangle_indices, TriangleUtilities::Normal(triangle));
            },
            "vertex_indices"_a,
            "Append a triangle from three vertex indices."
        );

    py::class_<ElementViewType>(m, "Element", "Read-only view of an active background-grid element.")
        .def_property_readonly(
            "integration_points",
            [](const ElementViewType& rElement) -> const IntegrationPointVectorType& {
                return rElement.GetIntegrationPoints<CoordinateSpace::parametric>();
            },
            py::return_value_policy::reference_internal,
            "Parametric integration points."
        )
        .def_property_readonly(
            "lower_bound_xyz",
            [](const ElementViewType& rElement) { return rElement.GetCellBounds<CoordinateSpace::global>().lower; },
            "Lower global-coordinate bound."
        )
        .def_property_readonly(
            "upper_bound_xyz",
            [](const ElementViewType& rElement) { return rElement.GetCellBounds<CoordinateSpace::global>().upper; },
            "Upper global-coordinate bound."
        )
        .def_property_readonly(
            "lower_bound_uvw",
            [](const ElementViewType& rElement) { return rElement.GetCellBounds<CoordinateSpace::parametric>().lower; },
            "Lower parametric-coordinate bound."
        )
        .def_property_readonly(
            "upper_bound_uvw",
            [](const ElementViewType& rElement) { return rElement.GetCellBounds<CoordinateSpace::parametric>().upper; },
            "Upper parametric-coordinate bound."
        )
        .def_property_readonly("id", &ElementViewType::GetId, "Element identifier.")
        .def_property_readonly(
            "is_trimmed", &ElementViewType::IsTrimmed, "Whether the element intersects the boundary."
        );

    using GridType = BackgroundGridType;
    using ElementRangeType = decltype(std::declval<const GridType&>().GetElementViews());
    using PythonElementRange = PythonRange<ElementRangeType>;
    BindRange<ElementRangeType>(InternalModule, "_ElementRange");

    py::class_<GridType>(m, "BackgroundGrid", "Background grid containing active elements and conditions.")
        .def(py::init<const GridType::MainDictionaryType&>())
        .def_property_readonly(
            "elements",
            py::cpp_function(
                [](const GridType& rGrid) { return PythonElementRange{ rGrid.GetElementViews() }; },
                py::keep_alive<0, 1>()
            )
        )
        .def_property_readonly("num_active_elements", &GridType::NumberOfActiveElements, "Number of active elements.")
        .def_property_readonly(
            "conditions", &GridType::GetConditions, py::return_value_policy::reference_internal, "Read-only conditions."
        )
        .def_property_readonly("num_conditions", &GridType::NumberOfConditions, "Number of conditions.");

    py::class_<ConditionSegmentType>(m, "ConditionSegment", "Read-only clipped condition surface segment.")
        .def_property_readonly(
            "triangle_mesh",
            &ConditionSegmentType::GetTriangleMesh,
            py::return_value_policy::reference_internal,
            "Boundary triangle mesh for this segment."
        );

    using ConditionSegmentRangeType = std::span<const ConditionSegmentType>;
    using PythonConditionSegmentRange = PythonRange<ConditionSegmentRangeType>;
    BindRange<ConditionSegmentRangeType>(InternalModule, "_ConditionSegmentRange");

    py::class_<ConditionType, ConditionPtrType>(m, "Condition", "Boundary or loading condition.")
        .def_static("is_weak_condition", []() -> bool { return true; })
        .def_property_readonly(
            "settings",
            &ConditionType::GetSettings,
            py::return_value_policy::reference_internal,
            "Read-only condition settings."
        )
        .def_property_readonly(
            "segments",
            py::cpp_function(
                [](const ConditionType& rCondition) { return PythonConditionSegmentRange{ rCondition.GetSegments() }; },
                py::keep_alive<0, 1>()
            )
        )
        .def_property_readonly("num_segments", &ConditionType::NumberOfSegments, "Number of clipped surface segments.");

    py::class_<EmbeddedModel>(InternalModule, "_EmbeddedModel", "Private C++ model implementation.")
        .def(py::init([](MainDictionaryHolderType& rSettings) {
            return MakeUnique<EmbeddedModel>(EmbeddedModel::Create(rSettings.Release()));
        }))
        .def("create_all_from_settings", &EmbeddedModel::CreateAllFromSettings)
        .def_property_readonly(
            "elements",
            py::cpp_function(
                [](const EmbeddedModel& rEmbeddedModel) {
                    return PythonElementRange{ rEmbeddedModel.GetElementViews() };
                },
                py::keep_alive<0, 1>()
            )
        )
        .def_property_readonly("conditions", &EmbeddedModel::GetConditions, py::return_value_policy::reference_internal)
        .def_property_readonly("settings", &EmbeddedModel::GetSettings, py::return_value_policy::reference_internal)
        .def_property_readonly(
            "model_info",
            static_cast<const EmbeddedModel::MainDictionaryType& (EmbeddedModel::*)() const>(
                &EmbeddedModel::GetModelInfo
            ),
            py::return_value_policy::reference_internal
        );
}

}  // namespace queso::python
