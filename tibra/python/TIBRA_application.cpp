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
#include "geometries/element.h"
#include "geometries/element_container.h"
#include "geometries/triangle_3d_3n.h"
#include "geometries/integration_point.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "io/io_utilities.h"

typedef std::vector<std::array<double,2>> IntegrationPoint1DVectorType;
typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
typedef std::vector<std::shared_ptr<Element>> ElementVectorPtrType;
typedef std::vector<Triangle3D3N> TriangleVectorType;

PYBIND11_MAKE_OPAQUE(IntegrationPoint1DVectorType);
PYBIND11_MAKE_OPAQUE(IntegrationPointVectorType);
PYBIND11_MAKE_OPAQUE(ElementVectorPtrType);
PYBIND11_MAKE_OPAQUE(TriangleVectorType);

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


//ptr_wrapper<double> get_ptr(const TIBRA& v) { return const_cast<double*>(v.GetMeshPoints()); }
// double array[3] = { 3.14, 2.18, -1 };
// static ptr_wrapper<double> get_ptr() { return ptr_wrapper<double>(array, 3); }

PYBIND11_MODULE(TIBRA_Application,m) {

    m.doc() = "This is a Python binding for TIBRA";

    py::class_<ptr_wrapper<double>>(m,"pdouble")
        .def(py::init<>())
        .def("__len__", [](const ptr_wrapper<double> &v) { return v.get_size(); })
        .def("__getitem__",  [](const ptr_wrapper<double> &v, unsigned int i){return v[i];} )
        .def("__iter__", [](ptr_wrapper<double> &v) {
            return py::make_iterator(v.get(), v.get() + v.get_size()) ;
        }, py::keep_alive<0, 1>())
        ;

    py::class_<IntegrationPoint, std::shared_ptr<IntegrationPoint>>(m, "IntegrationPoint")
        .def(py::init<double, double, double, double>())
        .def("GetX", &IntegrationPoint::X)
        .def("GetY", &IntegrationPoint::Y)
        .def("GetZ", &IntegrationPoint::Z)
        .def("Coordinates", &IntegrationPoint::Coordinates)
        .def("GetWeight", &IntegrationPoint::GetWeight)
        .def("SetWeight", &IntegrationPoint::SetWeight)
    ;

    py::bind_vector<IntegrationPointVectorType,std::shared_ptr<IntegrationPointVectorType>>
        (m, "VectorOfIntegrationPoints")
    ;

    py::class_<Triangle3D3N, std::shared_ptr<Triangle3D3N>>(m,"Triangle3D3N")
        .def(py::init<std::array<double, 3>, std::array<double, 3>, std::array<double, 3>,std::array<double, 3>>())
        .def("Center", &Triangle3D3N::Center)
        .def("Jacobian", &Triangle3D3N::Jacobian)
        .def("Normal", &Triangle3D3N::Normal)
        .def("Area", &Triangle3D3N::Area)
        .def("GetIntegrationPointsGlobal", [](Triangle3D3N& self, IndexType Method){
            IntegrationPointVectorPtrType ptr;
            ptr = self.GetIntegrationPointsGlobal(Method);

            return *ptr;
        })
        .def("P1", &Triangle3D3N::P1)
        .def("P2", &Triangle3D3N::P2)
        .def("P3", &Triangle3D3N::P3)
    ;

    py::bind_vector<TriangleVectorType,std::shared_ptr<TriangleVectorType>>
        (m, "VectorOfTriangles")
    ;

    py::class_<ElementVectorPtrType>(m, "ElementVector")
        .def(py::init<>())
        .def("__len__", [](const ElementVectorPtrType &v) { return v.size(); })
        .def("__iter__", [](ElementVectorPtrType &v) {
            return py::make_iterator( v.begin(), v.end() );
        }, py::keep_alive<0, 1>())
        ;


    py::class_<Element, std::shared_ptr<Element>>(m,"Element")
        .def("GetIntegrationPointsTrimmed",  &Element::GetIntegrationPointsTrimmed, py::return_value_policy::reference_internal )
        .def("GetIntegrationPointsInside",  &Element::GetIntegrationPointsInside, py::return_value_policy::reference_internal )
        .def("GetIntegrationPointsFictitious",  &Element::GetIntegrationPointsFictitious, py::return_value_policy::reference_internal )
        .def("GetTriangles", &Element::GetTriangles, py::return_value_policy::reference_internal)
        .def("GetNeumannTriangles", &Element::GetNeumannTriangles, py::return_value_policy::reference_internal)
        .def("GetDirichletTriangles", &Element::GetDirichletTriangles, py::return_value_policy::reference_internal)
        .def("GetLocalLowerPoint", &Element::GetLocalLowerPoint)
        .def("GetLocalUpperPoint", &Element::GetLocalUpperPoint)
        .def("GetNumberBoundaryTriangles", &Element::GetNumberBoundaryTriangles)
        .def("ID", &Element::GetId)
        .def("IsTrimmed", &Element::IsTrimmed)
    ;

    py::bind_vector<IntegrationPoint1DVectorType,std::unique_ptr<IntegrationPoint1DVectorType>>
        (m, "VectorOfIntegrationPoints1D")
    ;

    py::enum_<IntegrationPointFactory1D::IntegrationMethod>(m, "IntegrationMethod")
        .value("Gauss", IntegrationPointFactory1D::IntegrationMethod::Gauss)
        .value("ReducedGauss1", IntegrationPointFactory1D::IntegrationMethod::ReducedGauss1)
        .value("ReducedGauss2", IntegrationPointFactory1D::IntegrationMethod::ReducedGauss2)
        .value("ReducedExact", IntegrationPointFactory1D::IntegrationMethod::ReducedExact)
        .value("ReducedOrder1", IntegrationPointFactory1D::IntegrationMethod::ReducedOrder1)
        .value("ReducedOrder2", IntegrationPointFactory1D::IntegrationMethod::ReducedOrder2)
        .export_values()
    ;

    py::class_<IntegrationPointFactory1D, std::shared_ptr<IntegrationPointFactory1D>>(m,"IntegrationPointFactory1D")
        .def_static("GetGGQ", &IntegrationPointFactory1D::GetGGQ, py::return_value_policy::move)
    ;

    py::class_<TIBRA,std::shared_ptr<TIBRA>>(m,"TIBRA")
        .def(py::init<const std::string, std::array<double, 3>, std::array<double, 3>, std::array<int, 3>, std::array<int, 3>, double, int, double, double, std::string, int>())
        .def(py::init<const std::string, std::array<double, 3>, std::array<double, 3>, std::array<int, 3>, std::array<int, 3>, double, int, double, double, std::string, int, bool>())
        .def("GetElements",  &TIBRA::GetElements, py::return_value_policy::reference_internal )
        .def("ReadWritePostMesh", &TIBRA::ReadWritePostMesh )
        .def("GetPostMeshPointsRaw", [](const TIBRA& v){

            auto& mesh = v.GetPostMesh();
            const std::size_t num_p = mesh.num_vertices();
            const auto raw_ptr = mesh.points().data()->cartesian_begin();
            return  ptr_wrapper<double>(raw_ptr, num_p*3);
        });
    ;


    m.def("WriteDisplacementToVTK", &IO::WriteDisplacementToVTK);
}

}// namespace Python
