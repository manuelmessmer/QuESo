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
#include "queso/containers/triangle_mesh.hpp"
#include "queso/containers/background_grid.hpp"
#include "queso/containers/condition.hpp"
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "queso/embedded_model.h"

// Note: PYBIND11_MAKE_OPAQUE must live at file scope.

// Point types
using PointVectorType = std::vector<queso::PointType>;
PYBIND11_MAKE_OPAQUE(PointVectorType)
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

void AddContainersToPython(pybind11::module& m) {

    /// Export PointType
    py::class_<PointType>(m,"Point")
        .def(py::init<std::array<double,3>>())
        .def(py::init<double, double, double>())
        .def("X", [](const PointType& self){return self[0];} )
        .def("Y", [](const PointType& self){return self[1];} )
        .def("Z", [](const PointType& self){return self[2];} )
        .def("__getitem__",  [](const PointType &self, IndexType i){return self[i];} )
        .def("__len__", [](const PointType&) { return 3; })
        .def("__iter__", [](const PointType &self) {
                return py::make_iterator(self.begin(), self.end());
            }, py::keep_alive<0, 1>())  // Keep PointType alive while iterator is used
        .def("__repr__", [](const PointType &self) {
                return "Point(" + std::to_string(self[0])
                    + ", " + std::to_string(self[1])
                    + ", " + std::to_string(self[2]) + ")";
            })
        ;

    /// Export PointVector
    py::bind_vector<PointVectorType>(m, "PointVector");

    /// Export Vector3i
    py::class_<Vector3i>(m,"Vector3i")
        .def(py::init<std::array<IndexType,3>>())
        .def(py::init<IndexType, IndexType, IndexType>())
        .def("X", [](const Vector3i& self){return self[0];} )
        .def("Y", [](const Vector3i& self){return self[1];} )
        .def("Z", [](const Vector3i& self){return self[2];} )
        .def("__getitem__",  [](const Vector3i &self, IndexType i){return self[i];} )
        .def("__len__", [](const Vector3i&) { return 3; })
        .def("__iter__", [](const Vector3i &self) {
                return py::make_iterator(self.begin(), self.end());
            }, py::keep_alive<0, 1>())  // Keep Vector3i alive while iterator is used
        .def("__repr__", [](const Vector3i &self) {
                std::ostringstream oss;
                oss << "Vector3i(" << self[0] << ", " << self[1] << ", " << self[2] << ")";
                return oss.str();
            })
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

    /// Export TriangleMeshInterface
    py::class_<TriangleMeshInterface, Unique<TriangleMeshInterface>>(m,"TriangleMeshInterface")
    ;

    /// Export TriangleMesh
    py::class_<TriangleMesh, Unique<TriangleMesh>, TriangleMeshInterface>(m,"TriangleMesh")
        .def(py::init<>())
        .def("Center", &TriangleMesh::Center)
        .def("Normal", [](TriangleMesh& self, IndexType Id){
            return self.Normal(Id); })
        .def("Area", [](TriangleMesh& self, IndexType Id){
            return self.Area(Id); })
        .def("GetIntegrationPointsGlobal", [](TriangleMesh& self, IndexType Id, IndexType Method){
            return self.pGetIPsGlobal<BoundaryIntegrationPointType>(Id, Method); }, py::return_value_policy::move)
        .def("NumOfTriangles", &TriangleMesh::NumOfTriangles)
        .def("P1", &TriangleMesh::P1)
        .def("P2", &TriangleMesh::P2)
        .def("P3", &TriangleMesh::P3)
        .def("Reserve", &TriangleMesh::Reserve)
        .def("AddVertex", [](TriangleMesh& self, const std::array<double,3>& rVertex){
            return self.AddVertex( PointType(rVertex) ); })
        .def("AddTriangle", [](TriangleMesh& self, const std::array<IndexType,3>& rTriangle){
            return self.AddTriangle( Vector3i(rTriangle) ); })
        .def("AddNormal", [](TriangleMesh& self, const PointType& rNormal){
            return self.AddNormal( rNormal ); })
        .def("GetVertices",  static_cast< const std::vector<PointType>& (TriangleMesh::*)() const>(&TriangleMesh::GetVertices)
            , py::return_value_policy::reference_internal ) // Export const version
        .def("GetTriangles",  static_cast< const std::vector<Vector3i>& (TriangleMesh::*)() const>(&TriangleMesh::GetTriangles)
            , py::return_value_policy::reference_internal ) // Export const version
        .def_static("AspectRatioStatic", [](const std::array<double,3>& rV1, const std::array<double,3>& rV2, const std::array<double,3>& rV3){
            return TriangleMesh::AspectRatio( PointType(rV1), PointType(rV2), PointType(rV3) ); })
        .def_static("NormalStatic", [](const std::array<double,3>& rV1, const std::array<double,3>& rV2, const std::array<double,3>& rV3){
            return TriangleMesh::Normal( PointType(rV1), PointType(rV2), PointType(rV3)); })
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

