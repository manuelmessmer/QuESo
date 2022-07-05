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
#include "TrIGA.hpp"
#include "geometries/element.h"
#include "geometries/element_container.h"
#include "geometries/triangle_3d_3n.h"
#include "geometries/integration_point.h"

typedef std::vector<IntegrationPoint> IntegrationPointType;
typedef std::vector<std::shared_ptr<Element>> ElementVectorPtrType;
typedef std::vector<Triangle3D3N> TriangleVectorType;

PYBIND11_MAKE_OPAQUE(IntegrationPointType);
PYBIND11_MAKE_OPAQUE(ElementVectorPtrType);
PYBIND11_MAKE_OPAQUE(TriangleVectorType);

namespace Python {

namespace py = pybind11;

static const std::vector<std::array<double, 2>>& GetIntegrationPointsexport( int PolynomialDegree, int NumberKnotSpans  ){

    return IntegrationPointFactory::GetIntegrationPoints(PolynomialDegree, NumberKnotSpans, IntegrationPointFactory::IntegrationMethod::ReducedOrder2 );
}

PYBIND11_MODULE(TrIGA_Application,m) {

    m.doc() = "This is a Python binding for the STLEmbedder";

    py::class_<IntegrationPoint, std::shared_ptr<IntegrationPoint>>(m, "IntegrationPoint")
        .def(py::init<double, double, double, double>())
        .def("GetX", &IntegrationPoint::X)
        .def("GetY", &IntegrationPoint::Y)
        .def("GetZ", &IntegrationPoint::Z)
        .def("Coordinates", &IntegrationPoint::Coordinates)
        .def("GetWeight", &IntegrationPoint::GetWeight)
        .def("SetWeight", &IntegrationPoint::SetWeight)
    ;

    py::bind_vector<IntegrationPointType,std::shared_ptr<IntegrationPointType>>
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
        .def("MfIterationsResidual", &Element::MfIterationsResidual)
        .def("MfIterationsPoints", &Element::MfIterationsPoints)
        .def("ID", &Element::GetId)
        .def("IsTrimmed", &Element::IsTrimmed)
    ;


    py::class_<TrIGA,std::shared_ptr<TrIGA>>(m,"TrIGA")
        .def(py::init<const std::string, std::array<double, 3>, std::array<double, 3>, std::array<int, 3>, std::array<int, 3>, double, int, double, double, std::string, int>())
        .def(py::init<const std::string, std::array<double, 3>, std::array<double, 3>, std::array<int, 3>, std::array<int, 3>, double, int, double, double, std::string, int, bool>())
        .def("GetElements",  &TrIGA::ExportElements, py::return_value_policy::reference_internal )
    ;


    m.def("GetIntegrationPoints", &GetIntegrationPointsexport )
    ;
}

}// namespace Python
