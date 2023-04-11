// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <iostream>
#include <vector>
//// Project includes
#include "TIBRA_main.h"
#include "containers/element.hpp"
#include "containers/element_container.hpp"
#include "containers/triangle_mesh.hpp"
#include "containers/integration_point.hpp"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "utilities/mesh_utilities.h"
#include "io/io_utilities.h"

// Note: PYBIND11_MAKE_OPAQUE can not be captured within namespace
typedef std::vector<tibra::PointType> PointVectorType;
typedef std::vector<std::array<double,2>> IntegrationPoint1DVectorType;
typedef std::vector<tibra::IntegrationPoint> IntegrationPointVectorType;
typedef std::vector<std::shared_ptr<tibra::Element>> ElementVectorPtrType;
typedef std::vector<tibra::BoundaryIntegrationPoint> BoundaryIpVectorType;

PYBIND11_MAKE_OPAQUE(PointVectorType);
PYBIND11_MAKE_OPAQUE(BoundaryIpVectorType);
PYBIND11_MAKE_OPAQUE(IntegrationPoint1DVectorType);
PYBIND11_MAKE_OPAQUE(IntegrationPointVectorType);
PYBIND11_MAKE_OPAQUE(ElementVectorPtrType);

namespace tibra {
namespace Python {
namespace py = pybind11;

template <class T> class ptr_wrapper
{
    public:
        ptr_wrapper() : ptr(nullptr) {}
        ptr_wrapper(const T* ptr, std::size_t size) : ptr(ptr), size(size) {}
        ptr_wrapper(const ptr_wrapper& other) : ptr(other.ptr), size(other.size) {}
        const T& operator* () const { return *ptr; }
        const T* operator->() const { return  ptr; }
        const T* get() const { return ptr; }
        std::size_t get_size() const {return size; };
        void destroy() { delete ptr; }
        const T& operator[](std::size_t idx) const { return ptr[idx]; }
    private:
        const T* ptr;
        std::size_t size;
};

IntegrationMethodType GetIntegrationMethodFromString(const std::string& rValue){
    if( rValue == "Gauss" )
        return IntegrationMethod::Gauss;
    else if( rValue == "Gauss_Reduced1" )
        return IntegrationMethod::Gauss_Reduced1;
    else if( rValue == "Gauss_Reduced2" )
        return IntegrationMethod::Gauss_Reduced2;
    else if( rValue == "GGQ_Optimal")
        return IntegrationMethod::GGQ_Optimal;
    else if( rValue == "GGQ_Reduced1")
        return IntegrationMethod::GGQ_Reduced1;
    else if( rValue == "GGQ_Reduced2")
        return IntegrationMethod::GGQ_Reduced2;
    else
        TIBRA_ERROR("Parameters::GetIntegrationMethodFromString") << "Integration Method: " + rValue + " not available! \n";
}

PYBIND11_MODULE(TIBRA_Application,m) {

    m.doc() = "This is a Python binding for TIBRA";

    /// Required for GetPostMeshPointsRaw()
    py::class_<ptr_wrapper<double>>(m,"pdouble")
        .def(py::init<>())
        .def("__len__", [](const ptr_wrapper<double> &v) { return v.get_size(); })
        .def("__getitem__",  [](const ptr_wrapper<double> &v, unsigned int i){return v[i];} )
        .def("__iter__", [](ptr_wrapper<double> &v) {
            return py::make_iterator(v.get(), v.get() + v.get_size()) ;
        }, py::keep_alive<0, 1>())
        ;

    /// Export PointType
    py::class_<PointType, std::shared_ptr<PointType>>(m,"Point")
        .def(py::init<std::array<double,3>>())
        .def(py::init<double, double, double>())
        .def("GetX", static_cast< double (PointType::*)() const>(&PointType::X)) // Return const version of X()
        .def("GetY", static_cast< double (PointType::*)() const>(&PointType::Y)) // Return const version of Y()
        .def("GetZ", static_cast< double (PointType::*)() const>(&PointType::Z)) // Return const version of Z()
        .def("__getitem__",  [](const PointType &v, IndexType i){return v[i];} )
        ;

    /// Export PointType
    py::class_<Vector3i, std::shared_ptr<Vector3i>>(m,"Vector3i")
        .def(py::init<std::array<IndexType,3>>())
        .def(py::init<IndexType, IndexType, IndexType>())
        .def("GetX", static_cast< IndexType (Vector3i::*)() const>(&Vector3i::X)) // Return const version of X()
        .def("GetY", static_cast< IndexType (Vector3i::*)() const>(&Vector3i::Y)) // Return const version of Y()
        .def("GetZ", static_cast< IndexType (Vector3i::*)() const>(&Vector3i::Z)) // Return const version of Z()
        .def("__getitem__",  [](const Vector3i &v, IndexType i){return v[i];} )
        ;

    /// Export PointVector
    py::bind_vector<PointVectorType,Unique<PointVectorType>>
        (m, "PointVector")
    ;

    /// Export Integration Points 1D vector. Just a: (std::vector<std::array<double,2>>)
    py::bind_vector<IntegrationPoint1DVectorType,Unique<IntegrationPoint1DVectorType>>
        (m, "IntegrationPoint1DVector")
    ;

    /// Export Integration Points
    py::class_<IntegrationPoint, std::shared_ptr<IntegrationPoint>, PointType>(m, "IntegrationPoint")
        .def(py::init<double, double, double, double>())
        .def("GetWeight", &IntegrationPoint::GetWeight)
        .def("SetWeight", &IntegrationPoint::SetWeight)
    ;

    /// Export IntegrationPoint Vector
    py::bind_vector<IntegrationPointVectorType,std::shared_ptr<IntegrationPointVectorType>>
        (m, "IntegrationPointVector")
    ;

    /// Export BoundaryIntegrationPoint
    py::class_<BoundaryIntegrationPoint, std::shared_ptr<BoundaryIntegrationPoint>, IntegrationPoint>(m, "BoundaryIntegrationPoint")
        .def(py::init<double, double, double, double, const std::array<double,3>& >())
        .def("Normal", &BoundaryIntegrationPoint::Normal )
    ;

    /// Export BoundaryIntegrationPoint Vector
    py::bind_vector<BoundaryIpVectorType,Unique<BoundaryIpVectorType>>
        (m, "BoundaryIPVector")
    ;

    /// Export TriangleMesh
    py::class_<TriangleMesh, Unique<TriangleMesh>>(m,"TriangleMesh")
        .def(py::init<>())
        .def("Center", &TriangleMesh::Center)
        .def("Normal", &TriangleMesh::Normal)
        .def("Area", [](TriangleMesh& self, IndexType Id){
            return self.Area(Id);
        })
        .def("GetIntegrationPointsGlobal", [](TriangleMesh& self, IndexType Id, IndexType Method){
            return self.pGetIPsGlobal(Id, Method);
        })
        .def("NumOfTriangles", &TriangleMesh::NumOfTriangles)
        .def("P1", &TriangleMesh::P1)
        .def("P2", &TriangleMesh::P2)
        .def("P3", &TriangleMesh::P3)
        .def("Reserve", &TriangleMesh::Reserve)
    ;

    /// Export MeshUtilities
    py::class_<MeshUtilities, Unique<MeshUtilities>>(m,"MeshUtilities")
        .def_static("Append", [](TriangleMesh& rMesh, const TriangleMesh& rNewMesh){
            return MeshUtilities::Append(rMesh, rNewMesh);
        })
    ;

    /// Export Element
    py::class_<Element, std::shared_ptr<Element>>(m,"Element")
        .def("GetIntegrationPoints",  static_cast< const IntegrationPointVectorType& (Element::*)() const>(&Element::GetIntegrationPoints)
            ,py::return_value_policy::reference_internal ) // Export const version
        .def("GetTriangleMesh", [](const Element& rElement){
            return rElement.pGetTrimmedDomain()->GetTriangleMesh();
        }, py::return_value_policy::reference_internal)
        .def("GetBCTriangleMesh", [](const Element& rElement, std::function<bool(double, double,double)> &IsInDomain){
            return rElement.pGetTrimmedDomain()->pGetTriangleMesh(IsInDomain);
        })
        .def("GetLowerBoundParam", &Element::GetLowerBoundParam)
        .def("GetUpperBoundParam", &Element::GetUpperBoundParam)
        .def("GetNumberBoundaryTriangles", [](const Element& rElement ){
            return rElement.pGetTrimmedDomain()->GetTriangleMesh().NumOfTriangles();
        })
        .def("ID", &Element::GetId)
        .def("IsTrimmed", &Element::IsTrimmed)
    ;

    /// Export Element Vector
    py::class_<ElementVectorPtrType>(m, "ElementVector")
        .def(py::init<>())
        .def("__len__", [](const ElementVectorPtrType &v) { return v.size(); })
        .def("__iter__", [](ElementVectorPtrType &v) {
            return py::make_iterator( v.begin(), v.end() );
        }, py::keep_alive<0, 1>())
    ;

    /// Export enum IntegrationMethod
    py::enum_<IntegrationMethod>(m, "IntegrationMethod")
        .value("Gauss", IntegrationMethod::Gauss)
        .value("Gauss_Reduced1", IntegrationMethod::Gauss_Reduced1)
        .value("Gauss_Reduced2", IntegrationMethod::Gauss_Reduced2)
        .value("GGQ_Optimal", IntegrationMethod::GGQ_Optimal)
        .value("GGQ_Reduced1", IntegrationMethod::GGQ_Reduced1)
        .value("GGQ_Reduced2", IntegrationMethod::GGQ_Reduced2)
        .export_values()
    ;

    /// Export IntegrationPointFactory1D (mainly for Testing in py)
    py::class_<IntegrationPointFactory1D, std::shared_ptr<IntegrationPointFactory1D>>(m,"IntegrationPointFactory1D")
        .def_static("GetGGQ", &IntegrationPointFactory1D::GetGGQ, py::return_value_policy::move)
    ;

    /// Export Parameters
    py::class_<Parameters,std::shared_ptr<Parameters>>(m,"Parameters")
        .def(py::init<>())
        .def("Set",[](Parameters& rParams, const std::string& rName, const PointType& rValue){
            rParams.Set(rName, rValue); })
        .def("Set",[](Parameters& rParams, const std::string& rName, const std::array<double,3>& rValue){
            rParams.Set(rName, PointType(rValue)); })
        .def("Set",[](Parameters& rParams, const std::string& rName, const Vector3i& rValue){
            rParams.Set(rName, rValue); })
        .def("Set",[](Parameters& rParams, const std::string& rName, const std::array<IndexType,3>& rValue){
            rParams.Set(rName, Vector3i(rValue)); })
        .def("Set",[](Parameters& rParams, const std::string& rName, const bool rValue){
            rParams.Set(rName, rValue); })
        .def("Set",[](Parameters& rParams, const std::string& rName, const unsigned long rValue){
            rParams.Set(rName, rValue); })
        .def("Set",[](Parameters& rParams, const std::string& rName, const double rValue){
            rParams.Set(rName, rValue); })
        .def("Set",[](Parameters& rParams, const std::string& rName, const IntegrationMethodType& rValue){
            rParams.Set(rName, rValue); })
        .def("Set",[](Parameters& rParams, const std::string& rName, const std::string& rValue){
            // Allow integration method to be string. JSON is parsed as string.
            if( rName == "integration_method" ){
                rParams.Set(rName, GetIntegrationMethodFromString(rValue)); }
            else {
                rParams.Set(rName, rValue); } })
        .def("EchoLevel", &Parameters::EchoLevel)
        // Return std::array<type,3> types. Easier to handle in python.
        .def("LowerBound", []( const Parameters& rParams ) { return rParams.LowerBound().Coordinates(); })
        .def("UpperBound", []( const Parameters& rParams ) { return rParams.UpperBound().Coordinates(); })
        .def("Order", []( const Parameters& rParams ) { return rParams.Order().Coordinates(); })
        .def("NumberOfElements", []( const Parameters& rParams ) { return rParams.NumberOfElements().Coordinates(); })
        ;

    /// Export TIBRA
    py::class_<TIBRA,std::shared_ptr<TIBRA>>(m,"TIBRA")
        .def(py::init<const Parameters&>())
        .def("GetElements",  &TIBRA::GetElements, py::return_value_policy::reference_internal )
        .def("ReadWritePostMesh", &TIBRA::ReadWritePostMesh )
        .def("GetPostMeshPoints", [](const TIBRA& v){
            auto& mesh = v.GetPostMesh();
            return  mesh.GetVertices();
        });
    ;

    /// Export free floating function
    /// @todo wrap this in class
    m.def("WriteDisplacementToVTK", &IO::WriteDisplacementToVTK);
}

}// End namespace Python
}// End namespace tibra
