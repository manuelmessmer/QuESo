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

// Note: PYBIND11_MAKE_OPAQUE can not be captured within namespace
typedef std::vector<queso::PointType> PointVectorType;
typedef queso::IntegrationPoint IntegrationPointType;
typedef queso::BoundaryIntegrationPoint BoundaryIntegrationPointType;
typedef queso::Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;
typedef ElementType::IntegrationPoint1DVectorType IntegrationPoint1DVectorType;

typedef ElementType::IntegrationPointVectorType IntegrationPointVectorType;
typedef std::vector<BoundaryIntegrationPointType> BoundaryIpVectorType;

typedef std::vector<queso::Unique<ElementType>> ElementVectorPtrType;
typedef std::vector<queso::Condition> ConditionVectorType;

PYBIND11_MAKE_OPAQUE(PointVectorType);
PYBIND11_MAKE_OPAQUE(BoundaryIpVectorType);
PYBIND11_MAKE_OPAQUE(IntegrationPointVectorType);
PYBIND11_MAKE_OPAQUE(ElementVectorPtrType);
PYBIND11_MAKE_OPAQUE(ConditionVectorType);
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

    /// Export PointType
    py::class_<Vector3i, Unique<Vector3i>>(m,"Vector3i")
        .def(py::init<std::array<IndexType,3>>())
        .def(py::init<IndexType, IndexType, IndexType>())
        .def("X", [](const Vector3i& self){return self[0];} )
        .def("Y", [](const Vector3i& self){return self[1];} )
        .def("Z", [](const Vector3i& self){return self[2];} )
        .def("__getitem__",  [](const Vector3i &v, IndexType i){return v[i];} )
        ;

    /// Export PointVector
    py::bind_vector<PointVectorType, Unique<PointVectorType>>
        (m, "PointVector")
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

    /// Export TriangleMesh
    py::class_<TriangleMesh, Unique<TriangleMesh>>(m,"TriangleMesh")
        .def(py::init<>())
        .def("Center", &TriangleMesh::Center)
        .def("Normal", [](TriangleMesh& self, IndexType Id){
            return self.Normal(Id); })
        .def("Area", [](TriangleMesh& self, IndexType Id){
            return self.Area(Id); })
        .def("GetIntegrationPointsGlobal", [](TriangleMesh& self, IndexType Id, IndexType Method){
            return self.pGetIPsGlobal<BoundaryIntegrationPointType>(Id, Method); })
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

    /// Export Element Container
    py::class_<BackgroundGrid<ElementType>, Unique<BackgroundGrid<ElementType>>>(m, "BackgroundGrid")
        .def(py::init<const SettingsBaseType&>())
        .def("__len__", [](const BackgroundGrid<ElementType> &rGrid) { return rGrid.size(); })
        .def("__iter__", [](BackgroundGrid<ElementType> &rGrid) {
            return py::make_iterator( rGrid.begin(), rGrid.end() );
        }, py::keep_alive<0, 1>())
    ;

    /// Export Condition
    py::class_<Condition, Unique<Condition>>(m,"Condition")
        .def("IsWeakCondition", [](const Condition& rCondition)->bool { return true; } )
        .def("GetTriangleMesh", &Condition::GetConformingMesh , py::return_value_policy::reference_internal )
        .def("GetSettings", &Condition::GetSettings)
    ;

    /// Export Condition Vector
    py::class_<ConditionVectorType, Unique<ConditionVectorType>>(m, "ConditionVector")
        .def(py::init<>())
        .def("__len__", [](const ConditionVectorType &v) { return v.size(); })
        .def("__iter__", [](ConditionVectorType &v) {
            return py::make_iterator( v.begin(), v.end() );
        }, py::keep_alive<0, 1>())
    ;

    /// Export Integration Points 1D vector. Just a: (std::vector<std::array<double,2>>)
    py::bind_vector<IntegrationPoint1DVectorType, Unique<IntegrationPoint1DVectorType>>
        (m, "IntegrationPoint1DVector")
    ;

    /// Export IntegrationPointFactory1D (mainly for Testing in py)
    py::class_<IntegrationPointFactory1D>(m,"IntegrationPointFactory1D")
        .def_static("GetGGQ", &IntegrationPointFactory1D::GetGGQ, py::return_value_policy::move)
    ;

} // End AddContainersToPython

} // End namespace Python
} // End namespace queso

