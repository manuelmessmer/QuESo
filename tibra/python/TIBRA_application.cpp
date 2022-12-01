// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <iostream>
#include <vector>

// Project includes
#include "TIBRA_main.hpp"
#include "containers/element.h"
#include "containers/element_container.h"
#include "containers/triangle_mesh.h"
#include "containers/integration_point.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"
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
        .def(py::init<double, double, double>())
        .def("GetX", static_cast< double (PointType::*)() const>(&PointType::X)) // Return const version of X()
        .def("GetY", static_cast< double (PointType::*)() const>(&PointType::Y)) // Return const version of Y()
        .def("GetZ", static_cast< double (PointType::*)() const>(&PointType::Z)) // Return const version of Z()
        .def("__getitem__",  [](const PointType &v, IndexType i){return v[i];} )
        ;

    /// Export PointVector
    py::bind_vector<PointVectorType,std::unique_ptr<PointVectorType>>
        (m, "PointVector")
    ;

    /// Export Integration Points 1D vector. Just a: (std::vector<std::array<double,2>>)
    py::bind_vector<IntegrationPoint1DVectorType,std::unique_ptr<IntegrationPoint1DVectorType>>
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
    py::bind_vector<BoundaryIpVectorType,std::unique_ptr<BoundaryIpVectorType>>
        (m, "BoundaryIPVector")
    ;

    /// Export TriangleMesh
    py::class_<TriangleMesh, std::unique_ptr<TriangleMesh>>(m,"TriangleMesh")
        .def(py::init<>())
        .def("Center", &TriangleMesh::Center)
        .def("Normal", &TriangleMesh::Normal)
        .def("Area", [](TriangleMesh& self, IndexType Id){
            return self.Area(Id);
        })
        .def("GetIntegrationPointsGlobal", [](TriangleMesh& self, IndexType Id, IndexType Method){
            return self.GetIPsGlobal(Id, Method);
        })
        .def("Append", [](TriangleMesh& self, TriangleMesh& rOthers){
            return self.Append(rOthers);
        })
        .def("NumOfTriangles", &TriangleMesh::NumOfTriangles)
        .def("P1", &TriangleMesh::P1)
        .def("P2", &TriangleMesh::P2)
        .def("P3", &TriangleMesh::P3)
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
    py::enum_<IntegrationPointFactory1D::IntegrationMethod>(m, "IntegrationMethod")
        .value("Gauss", IntegrationPointFactory1D::IntegrationMethod::Gauss)
        .value("ReducedGauss1", IntegrationPointFactory1D::IntegrationMethod::ReducedGauss1)
        .value("ReducedGauss2", IntegrationPointFactory1D::IntegrationMethod::ReducedGauss2)
        .value("ReducedExact", IntegrationPointFactory1D::IntegrationMethod::ReducedExact)
        .value("ReducedOrder1", IntegrationPointFactory1D::IntegrationMethod::ReducedOrder1)
        .value("ReducedOrder2", IntegrationPointFactory1D::IntegrationMethod::ReducedOrder2)
        .export_values()
    ;

    /// Export IntegrationPointFactory1D (mainly for Testing in py)
    py::class_<IntegrationPointFactory1D, std::shared_ptr<IntegrationPointFactory1D>>(m,"IntegrationPointFactory1D")
        .def_static("GetGGQ", &IntegrationPointFactory1D::GetGGQ, py::return_value_policy::move)
    ;

    /// Export TIBRA
    py::class_<TIBRA,std::shared_ptr<TIBRA>>(m,"TIBRA")
        .def(py::init<const std::string, std::array<double,3>, std::array<double,3>, std::array<IndexType,3>, std::array<IndexType,3>, double, int, double, double, std::string, int>())
        .def(py::init<const std::string, std::array<double,3>, std::array<double,3>, std::array<IndexType,3>, std::array<IndexType,3>, double, int, double, double, std::string, int, bool>())
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
