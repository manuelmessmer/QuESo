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
#include "queso/python/define_python.hpp"
#include "queso/python/add_containers_to_python.h"
// To export
#include "queso/containers/triangle_mesh.hpp"
#include "queso/containers/background_grid.hpp"
#include "queso/containers/condition.hpp"
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "queso/embedded_model.h"

// Note: PYBIND11_MAKE_OPAQUE can not be captured within namespace
typedef std::vector<queso::PointType> PointVectorType;
PYBIND11_MAKE_OPAQUE(PointVectorType);

typedef queso::IntegrationPoint IntegrationPointType;
typedef queso::BoundaryIntegrationPoint BoundaryIntegrationPointType;
typedef queso::Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;
typedef queso::Unique<ElementType> ElementPtrType;
PYBIND11_MAKE_OPAQUE(ElementPtrType);
typedef std::vector<ElementPtrType> ElementVectorPtrType;
PYBIND11_MAKE_OPAQUE(ElementVectorPtrType);

typedef queso::Condition<ElementType> ConditionType;
typedef queso::Unique<ConditionType> ConditionPtrType;
typedef std::vector<ConditionPtrType> ConditionVectorPtrType;
PYBIND11_MAKE_OPAQUE(ConditionVectorPtrType);

typedef queso::ConditionSegment<ElementType> ConditionSegmentType;
typedef queso::Unique<ConditionSegmentType> ConditionSegmentPtrType;
PYBIND11_MAKE_OPAQUE(ConditionSegmentPtrType)
typedef std::vector<ConditionSegmentPtrType> ConditionSegmentVectorPtrType;
PYBIND11_MAKE_OPAQUE(ConditionSegmentVectorPtrType)

typedef ElementType::IntegrationPointVectorType IntegrationPointVectorType;
PYBIND11_MAKE_OPAQUE(IntegrationPointVectorType);

typedef ElementType::BoundaryIntegrationPointVectorType BoundaryIpVectorType;
PYBIND11_MAKE_OPAQUE(BoundaryIpVectorType);

typedef std::vector<std::array<double,2>> IntegrationPoint1DVectorType;
PYBIND11_MAKE_OPAQUE(IntegrationPoint1DVectorType);

namespace queso {
namespace Python {

namespace py = pybind11;

void AddContainersToPython(pybind11::module& m) {

    /// Export PointType
    py::class_<PointType, Unique<PointType>>(m,"Point")
        .def(py::init<std::array<double,3>>())
        .def(py::init<double, double, double>())
        .def("X", [](const PointType& self){return self[0];} )
        .def("Y", [](const PointType& self){return self[1];} )
        .def("Z", [](const PointType& self){return self[2];} )
        .def("__getitem__",  [](const PointType &v, IndexType i){return v[i];} )
        ;

    /// Export PointVector
    py::bind_vector<PointVectorType, Unique<PointVectorType>>
        (m, "PointVector")
    ;

    /// Export Vector3i
    py::class_<Vector3i, Unique<Vector3i>>(m,"Vector3i")
        .def(py::init<std::array<IndexType,3>>())
        .def(py::init<IndexType, IndexType, IndexType>())
        .def("X", [](const Vector3i& self){return self[0];} )
        .def("Y", [](const Vector3i& self){return self[1];} )
        .def("Z", [](const Vector3i& self){return self[2];} )
        .def("__getitem__",  [](const Vector3i &v, IndexType i){return v[i];} )
        ;

    /// Export Integration Points
    py::class_<IntegrationPointType, Unique<IntegrationPointType>>(m, "IntegrationPoint")
        .def(py::init<double, double, double, double>())
        .def("X", [](const IntegrationPointType& self){return self[0];} )
        .def("Y", [](const IntegrationPointType& self){return self[1];} )
        .def("Z", [](const IntegrationPointType& self){return self[2];} )
        .def("__getitem__",  [](const IntegrationPointType &v, IndexType i){return v[i];} )
        .def("Weight", &IntegrationPointType::Weight)
        .def("SetWeight", &IntegrationPointType::SetWeight)
    ;

    /// Export IntegrationPoint Vector
    py::bind_vector<IntegrationPointVectorType, Unique<IntegrationPointVectorType>>
        (m, "IntegrationPointVector")
    ;

    /// Export BoundaryIntegrationPoint
    py::class_<BoundaryIntegrationPointType, Unique<BoundaryIntegrationPointType>, IntegrationPointType>(m, "BoundaryIntegrationPoint")
        .def(py::init<double, double, double, double, const std::array<double,3>& >())
        .def("Normal", &BoundaryIntegrationPointType::Normal )
    ;

    /// Export BoundaryIntegrationPoint Vector
    py::bind_vector<BoundaryIpVectorType, Unique<BoundaryIpVectorType>>
        (m, "BoundaryIPVector")
    ;

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
    py::class_<ElementType, Unique<ElementType>>(m,"Element")
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
    py::class_<ElementVectorPtrType, Unique<ElementVectorPtrType>>(m, "ElementVector")
        .def("__getitem__", [](const ElementVectorPtrType &v, const IndexType i) { return &(*v[i]); })
        .def("__len__", [](const ElementVectorPtrType &v) { return v.size(); })
        .def("__iter__", [](ElementVectorPtrType &v) {
            return py::make_iterator( dereference_iterator(v.begin()), dereference_iterator(v.end()) );
        }, py::keep_alive<0, 1>() )
    ;

    /// Export BackgroundGrid
    py::class_<BackgroundGrid<ElementType>, Unique<BackgroundGrid<ElementType>>>(m, "BackgroundGrid")
        .def(py::init<const Settings&>())
        .def("GetElements", &BackgroundGrid<ElementType>::GetElements, py::return_value_policy::reference_internal)
        .def("NumberOfActiveElements", &BackgroundGrid<ElementType>::NumberOfActiveElements)
        .def("GetConditions", &BackgroundGrid<ElementType>::GetConditions, py::return_value_policy::reference_internal)
        .def("NumberOfConditions", &BackgroundGrid<ElementType>::NumberOfConditions)
    ;

    /// Export Condition Segment
    py::class_<ConditionSegmentType, Unique<ConditionSegmentType>>(m,"ConditionSegment")
        .def("GetTriangleMesh", &ConditionSegmentType::GetTriangleMesh , py::return_value_policy::reference_internal )
    ;

    // Export ConditionSegment Vector
    py::class_<ConditionSegmentVectorPtrType, Unique<ConditionSegmentVectorPtrType>>(m, "ConditionSegmentVector")
        .def("__getitem__", [](const ConditionSegmentVectorPtrType &v, const IndexType i) { return &(*v[i]); })
        .def("__len__", [](const ConditionSegmentVectorPtrType &v) { return v.size(); })
        .def("__iter__", [](ConditionSegmentVectorPtrType &v) {
            return py::make_iterator( dereference_iterator(v.begin()), dereference_iterator(v.end()) );
        }, py::keep_alive<0, 1>() )
    ;

    /// Export Condition
    py::class_<ConditionType, Unique<ConditionType>>(m,"Condition")
        .def("is_weak_condition", [](const ConditionType& rCondition)->bool { return true; } )
        .def("GetSettings", &ConditionType::GetSettings, py::return_value_policy::reference_internal)
        .def("GetSegments", &ConditionType::GetSegments, py::return_value_policy::reference_internal)
        .def("NumberOfSegments", &ConditionType::NumberOfSegments )
        .def("__len__", [](const ConditionType &v) { return v.NumberOfSegments(); })
        .def("__iter__", [](ConditionType &v) {
            return py::make_iterator( v.SegmentsBegin(), v.SegmentsEnd() );
        }, py::keep_alive<0, 1>())
    ;

    // Export Condition Vector
    py::class_<ConditionVectorPtrType, Unique<ConditionVectorPtrType>>(m, "ConditionVector")
        .def("__getitem__", [](const ConditionVectorPtrType &v, const IndexType i)
            { return &(*v[i]); }, py::return_value_policy::reference_internal)
        .def("__len__", [](const ConditionVectorPtrType &v) { return v.size(); })
        .def("__iter__", [](ConditionVectorPtrType &v) {
            return py::make_iterator( dereference_iterator(v.begin()), dereference_iterator(v.end()) );
        })
    ;

    /// Export Integration Points 1D vector. Just a: (std::vector<std::array<double,2>>)
    py::bind_vector<IntegrationPoint1DVectorType, Unique<IntegrationPoint1DVectorType>>
        (m, "IntegrationPoint1DVector")
    ;

    /// Export IntegrationPointFactory1D (mainly for Testing in py)
    py::class_<IntegrationPointFactory1D>(m,"IntegrationPointFactory1D")
        .def_static("GetGGQ", &IntegrationPointFactory1D::GetGGQ, py::return_value_policy::move)
    ;

    /// Export QuESo
    py::class_<EmbeddedModel>(m,"EmbeddedModel")
        .def(py::init<const Settings&>())
        .def("CreateAllFromSettings", &EmbeddedModel::CreateAllFromSettings)
        .def("GetElements", &EmbeddedModel::GetElements, py::return_value_policy::reference_internal)
        .def("GetConditions", &EmbeddedModel::GetConditions, py::return_value_policy::reference_internal )
        .def("GetSettings", &EmbeddedModel::GetSettings, py::return_value_policy::reference_internal)
        .def("GetModelInfo", &EmbeddedModel::GetModelInfo, py::return_value_policy::reference_internal)
    ;

} // End AddContainersToPython

} // End namespace Python
} // End namespace queso

