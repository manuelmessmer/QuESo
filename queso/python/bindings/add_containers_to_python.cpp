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

/// Project inlcudes
#include "queso/python/bindings/define_python.hpp"
#include "queso/python/bindings/add_containers_to_python.h"
// To export
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "queso/embedded_model.h"
#include "queso/utilities/triangle_utilities.hpp"

// Note: PYBIND11_MAKE_OPAQUE must live at file scope.

using IntegrationPointType = queso::IntegrationPoint;
using BoundaryIntegrationPointType = queso::BoundaryIntegrationPoint;

// Element related types
using ElementType = queso::Element<IntegrationPointType, BoundaryIntegrationPointType>;
using ElementPtrType = queso::Unique<ElementType>;
using ElementVectorPtrType = std::vector<ElementPtrType>;
PYBIND11_MAKE_OPAQUE(ElementVectorPtrType);

// Condition related types
using ConditionType = queso::Condition<ElementType>;
using ConditionPtrType = queso::Unique<ConditionType>;
using ConditionVectorPtrType = std::vector<ConditionPtrType>;
PYBIND11_MAKE_OPAQUE(ConditionVectorPtrType);

// Condition segments related types
using ConditionSegmentType = queso::ConditionSegment<ElementType>;
using ConditionSegmentPtrType = queso::Unique<ConditionSegmentType>;
using ConditionSegmentVectorPtrType = std::vector<ConditionSegmentPtrType>;
PYBIND11_MAKE_OPAQUE(ConditionSegmentVectorPtrType)

// Point / integration point vector types
using IntegrationPointVectorType = ElementType::IntegrationPointVectorType;
PYBIND11_MAKE_OPAQUE(IntegrationPointVectorType)

using BoundaryIpVectorType = ElementType::BoundaryIntegrationPointVectorType;
PYBIND11_MAKE_OPAQUE(BoundaryIpVectorType);

using IntegrationPoint1DVectorType = std::vector<std::array<double,2>>;
PYBIND11_MAKE_OPAQUE(IntegrationPoint1DVectorType);

namespace queso {
namespace Python {

using MainDictionaryHolderType = UniqueHolder<MainDictionaryType>;

namespace py = pybind11;

namespace {

using PythonTriangleProxy = TriangleProxy<WithNormals>;

struct TriangleIteratorTag {};

struct VertexIteratorTag {};

template<class TTag>
struct PythonMeshIterator {
    const TriangleMesh *mpMesh;
    IndexType mIndex;
};

using PythonTriangleIterator = PythonMeshIterator<TriangleIteratorTag>;
using PythonVertexIterator = PythonMeshIterator<VertexIteratorTag>;

inline std::tuple<double, double, double> ToTuple3(const Vector3d &rV)
{
    return { rV[0], rV[1], rV[2] };
}

} // namespace

void AddContainersToPython(pybind11::module& m) {
    py::class_<PythonTriangleIterator>(m, "TriangleIterator")
        .def("__iter__", [](PythonTriangleIterator &self) -> PythonTriangleIterator& { return self; }, py::return_value_policy::reference_internal)
        .def("__next__", [](PythonTriangleIterator &self) {
            if (self.mIndex >= self.mpMesh->NumOfTriangles()) {
                throw py::stop_iteration();
            }
            return self.mpMesh->Triangle<WithNormals>(self.mIndex++);
        })
    ;

    py::class_<PythonVertexIterator>(m, "VertexIterator")
        .def("__iter__", [](PythonVertexIterator &self) -> PythonVertexIterator& { return self; }, py::return_value_policy::reference_internal)
        .def("__next__", [](PythonVertexIterator &self) {
            if (self.mIndex >= self.mpMesh->NumOfVertices()) {
                throw py::stop_iteration();
            }
            return ToTuple3(self.mpMesh->Vertex(self.mIndex++));
        })
    ;

    py::class_<PythonTriangleProxy>(m, "Triangle")
        .def_property_readonly("p1", [](const PythonTriangleProxy &self){ return ToTuple3(Vector3d{self.P1[0], self.P1[1], self.P1[2]}); })
        .def_property_readonly("p2", [](const PythonTriangleProxy &self){ return ToTuple3(Vector3d{self.P2[0], self.P2[1], self.P2[2]}); })
        .def_property_readonly("p3", [](const PythonTriangleProxy &self){ return ToTuple3(Vector3d{self.P3[0], self.P3[1], self.P3[2]}); })
        .def_property_readonly("normal", [](const PythonTriangleProxy &self){ return ToTuple3(Vector3d{self.Normal[0], self.Normal[1], self.Normal[2]}); })
        .def("Area", [](const PythonTriangleProxy &self){
            return TriangleUtilities::Area(self);
        })
        .def("Center", [](const PythonTriangleProxy &self){
            return ToTuple3(TriangleUtilities::Center(self));
        })
        .def("AspectRatio", [](const PythonTriangleProxy &self){
            return TriangleUtilities::AspectRatio(self);
        })
        .def("GetIPsGlobal", [](const PythonTriangleProxy &self, IndexType Method){
            return TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPointType>(self, Method);
        }, py::return_value_policy::move)
    ;

    /// Export Integration Points
    py::class_<IntegrationPointType>(m, "IntegrationPoint")
        .def(py::init<double, double, double, double>())
        .def("X", [](const IntegrationPointType& self){return self[0];} )
        .def("Y", [](const IntegrationPointType& self){return self[1];} )
        .def("Z", [](const IntegrationPointType& self){return self[2];} )
        .def("Weight", &IntegrationPointType::Weight)
        .def("SetWeight", &IntegrationPointType::SetWeight)
        .def("__getitem__",  [](const IntegrationPointType &self, IndexType i){return self[i];} )
        .def("__len__", [](const IntegrationPointType&) { return 3; })
        .def("__iter__", [](const IntegrationPointType &self) {
            return py::make_iterator(self.data().begin(), self.data().end());
        }, py::keep_alive<0, 1>())
        .def("__repr__", [](const IntegrationPointType& self) {
            std::ostringstream oss;
            oss << "IntegrationPoint(" << self[0] << ", " << self[1] << ", " << self[2] << ", "
                << "Weight=" << self.Weight() << ")";
            return oss.str();
        });
    ;

    /// Export IntegrationPoint Vector
    py::bind_vector<IntegrationPointVectorType>(m, "IntegrationPointVector");

    /// Export BoundaryIntegrationPoint
    py::class_<BoundaryIntegrationPointType, IntegrationPointType>(m, "BoundaryIntegrationPoint")
        .def(py::init<double, double, double, double, const std::array<double,3>& >())
        .def("Normal", &BoundaryIntegrationPointType::Normal )
        .def("__repr__", [](const BoundaryIntegrationPointType& self) {
            std::ostringstream oss;
            oss << "IntegrationPoint(" << self[0] << ", " << self[1] << ", " << self[2]
                << ", Weight=" << self.Weight()
                << ", Normal=" << self.Normal() << ")";
            return oss.str();
        });
    ;

    /// Export BoundaryIntegrationPoint Vector
    py::bind_vector<BoundaryIpVectorType>(m, "BoundaryIPVector");

    /// Export TriangleMesh
    py::class_<TriangleMesh, Unique<TriangleMesh>>(m,"TriangleMesh")
        .def(py::init<>())
        .def("NumOfTriangles", &TriangleMesh::NumOfTriangles)
        .def("NumOfVertices", &TriangleMesh::NumOfVertices)
        .def("Vertices", [](const TriangleMesh &self){
            return PythonVertexIterator{ &self, 0 };
        }, py::keep_alive<0, 1>())
        .def("Triangles", [](const TriangleMesh &self){
            return PythonTriangleIterator{ &self, 0 };
        }, py::keep_alive<0, 1>())
        .def("Reserve", &TriangleMesh::Reserve)
        .def("AddVertex", [](TriangleMesh& self, const std::array<double,3>& rVertex){
            return self.AddVertex( PointType(rVertex) ); })
        .def("AddTriangle", [](TriangleMesh& self, const std::array<IndexType,3>& rTriangle){
            const Vector3i triangle_indices(rTriangle);
            const auto &vertices = self.Vertices();
            const TriangleProxy<WithoutNormals> triangle {
                std::span<const double, 3>(vertices[triangle_indices[0]].data(), 3),
                std::span<const double, 3>(vertices[triangle_indices[1]].data(), 3),
                std::span<const double, 3>(vertices[triangle_indices[2]].data(), 3)
            };
            self.AddTriangle(triangle_indices, TriangleUtilities::Normal(triangle));
        })
    ;

    /// Export Element
    py::class_<ElementType, ElementPtrType>(m,"Element")
        .def("GetIntegrationPoints",  static_cast< const IntegrationPointVectorType& (ElementType::*)() const>(&ElementType::GetIntegrationPoints)
            ,py::return_value_policy::reference_internal ) // Export const version
        .def("LowerBoundXYZ", [](const ElementType& rElement ){ return rElement.GetBoundsXYZ().first; })
        .def("UpperBoundXYZ", [](const ElementType& rElement ){ return rElement.GetBoundsXYZ().second; })
        .def("LowerBoundUVW", [](const ElementType& rElement ){ return rElement.GetBoundsUVW().first; })
        .def("UpperBoundUVW", [](const ElementType& rElement ){ return rElement.GetBoundsUVW().second; })
        .def("GetNumberBoundaryTriangles", [](const ElementType& rElement ){
            return rElement.pGetTrimmedDomain()->GetTriangleMesh().NumOfTriangles();
        })
        .def("ID", &ElementType::GetId)
        .def("IsTrimmed", &ElementType::IsTrimmed)
    ;

    // Export Element Vector
    py::class_<ElementVectorPtrType>(m, "ElementVector")
        .def("__getitem__", [](const ElementVectorPtrType &self, const IndexType i)
            { return &(*self[i]); }, py::return_value_policy::reference_internal)
        .def("__len__", [](const ElementVectorPtrType &self) { return self.size(); })
        .def("__iter__", [](ElementVectorPtrType &self) {
            return py::make_iterator( dereference_iterator(self.begin()), dereference_iterator(self.end()) );
        }, py::keep_alive<0, 1>() )
    ;

    /// Export BackgroundGrid
    py::class_<BackgroundGrid<ElementType>>(m, "BackgroundGrid")
        .def(py::init<const BackgroundGrid<ElementType>::MainDictionaryType&>())
        .def("GetElements", &BackgroundGrid<ElementType>::GetElements, py::return_value_policy::reference_internal)
        .def("NumberOfActiveElements", &BackgroundGrid<ElementType>::NumberOfActiveElements)
        .def("GetConditions", &BackgroundGrid<ElementType>::GetConditions, py::return_value_policy::reference_internal)
        .def("NumberOfConditions", &BackgroundGrid<ElementType>::NumberOfConditions)
    ;

    /// Export Condition Segment
    py::class_<ConditionSegmentType, ConditionSegmentPtrType>(m,"ConditionSegment")
        .def("GetTriangleMesh", &ConditionSegmentType::GetTriangleMesh , py::return_value_policy::reference_internal )
    ;

    // Export ConditionSegment Vector
    py::class_<ConditionSegmentVectorPtrType>(m, "ConditionSegmentVector")
        .def("__getitem__", [](const ConditionSegmentVectorPtrType &self, const IndexType i)
            { return &(*self[i]); }, py::return_value_policy::reference_internal)
        .def("__len__", [](const ConditionSegmentVectorPtrType &self) { return self.size(); })
        .def("__iter__", [](ConditionSegmentVectorPtrType &self) {
            return py::make_iterator( dereference_iterator(self.begin()), dereference_iterator(self.end()) );
        }, py::keep_alive<0, 1>() )
    ;

    /// Export Condition
    py::class_<ConditionType, ConditionPtrType>(m,"Condition")
        .def_static("is_weak_condition", []()->bool { return true; } )
        .def("GetSettings", &ConditionType::GetSettings, py::return_value_policy::reference_internal)
        .def("GetSegments", &ConditionType::GetSegments, py::return_value_policy::reference_internal)
        .def("NumberOfSegments", &ConditionType::NumberOfSegments )
        .def("__getitem__", [](const ConditionType &self, const IndexType i)
            { return &(*self.GetSegments()[i]); }, py::return_value_policy::reference_internal)
        .def("__len__", [](const ConditionType &self) { return self.NumberOfSegments(); })
        .def("__iter__", [](ConditionType &self) {
            return py::make_iterator( self.SegmentsBegin(), self.SegmentsEnd() );
        }, py::keep_alive<0, 1>())
    ;

    // Export Condition Vector
    py::class_<ConditionVectorPtrType>(m, "ConditionVector")
        .def("__getitem__", [](const ConditionVectorPtrType &self, const IndexType i)
            { return &(*self[i]); }, py::return_value_policy::reference_internal)
        .def("__len__", [](const ConditionVectorPtrType &self) { return self.size(); })
        .def("__iter__", [](ConditionVectorPtrType &self) {
            return py::make_iterator( dereference_iterator(self.begin()), dereference_iterator(self.end()) );
        }, py::keep_alive<0, 1>())
    ;

    /// Export Integration Points 1D vector. Just a: (std::vector<std::array<double,2>>)
    py::bind_vector<IntegrationPoint1DVectorType>(m, "IntegrationPoint1DVector");
    /// Export IntegrationPointFactory1D (mainly for Testing in py)
    py::class_<IntegrationPointFactory1D>(m,"IntegrationPointFactory1D")
        .def_static("GetGGQ", &IntegrationPointFactory1D::GetGGQ, py::return_value_policy::move)
    ;

    /// Export EmbeddedModel
    py::class_<EmbeddedModel>(m,"EmbeddedModel")
        .def(py::init([](MainDictionaryHolderType& rSettings) {
            return MakeUnique<EmbeddedModel>(EmbeddedModel::Create(rSettings.Release()));
        }))
        .def("CreateAllFromSettings", &EmbeddedModel::CreateAllFromSettings)
        .def("GetElements", &EmbeddedModel::GetElements, py::return_value_policy::reference_internal)
        .def("GetConditions", &EmbeddedModel::GetConditions, py::return_value_policy::reference_internal)
        .def("GetSettings", &EmbeddedModel::GetSettings, py::return_value_policy::reference_internal)
        .def("GetModelInfo", static_cast<const EmbeddedModel::MainDictionaryType& (EmbeddedModel::*)() const>(&EmbeddedModel::GetModelInfo),
            py::return_value_policy::reference_internal)
    ;

} // End AddContainersToPython

} // End namespace Python
} // End namespace queso

