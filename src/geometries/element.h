// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_H
#define ELEMENT_H

// VTK includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataNormals.h>

// Project includes
#include "geometries/integration_point.h"
#include "geometries/triangle_3d_3n.h"
#include "utilities/parameters.h"
#include "utilities/inside_test.h"
#include "utilities/mesh_utilities.h"

// External includes
#include "omp.h"
#include <memory>

class InsideTest;
class Parameters;

class Element
{

public:

    // Typedefs
    typedef std::shared_ptr<IntegrationPoint> IntegrationPointPtrType;
    typedef std::vector<IntegrationPointPtrType> IntegrationPointPtrVectorType;
    typedef std::vector<std::array<double, 2>> IntegrationPoint1DVectorType;
    typedef std::vector<Triangle3D3N> TriangleVectorType;
    typedef std::array<double,3> PointType;
    typedef std::size_t IndexType;

    typedef vtkSmartPointer<vtkPolyData> SurfaceMeshPtrType;
    typedef InsideTest InsideTestType;
    typedef std::vector<std::array<double,3>> BoundingBox;

    // Constructor
    Element(std::size_t ID, PointType PointLocalA, PointType PointLocalB, const Parameters& rParam) :
        mElementId(ID), mLocalLowerPoint(PointLocalA), mLocalUpperPoint(PointLocalB), mParameters(rParam)
    {
        mIsTrimmed = false;
        mSurfaceMeshSetFlag = false;
    }

    // Delete copy constructor
    Element(Element const& rOther) = delete;

    void SetIsTrimmed(bool value){
        mIsTrimmed = value;
    }

    void SetId(std::size_t value){
        mElementId = value;
    }

    const std::size_t GetId() const{
        return mElementId;
    }

    bool IsTrimmed(){
        return mIsTrimmed;
    }

    int GetNumberBoundaryTriangles(){
        return GetTriangles().size();
    }

    IntegrationPointPtrVectorType& GetIntegrationPointsTrimmed(){
        return mIntegrationPointsTrimmed;
    }
    IntegrationPointPtrVectorType& GetIntegrationPointsInside(){
        return mIntegrationPointsInside;
    }
    IntegrationPointPtrVectorType& GetIntegrationPointsFictitious(){
        return mIntegrationPointsFictitious;
    }

    TriangleVectorType& GetTriangles();

    TriangleVectorType& GetNeumannTriangles(std::function<bool(double, double,double)> &is_neumann){

        auto& triangles = this->GetTriangles();

        mNeumannTriangles.clear();
        for( auto triangle : triangles){
            if( is_neumann(triangle.P1()[0], triangle.P1()[1], triangle.P1()[2]) ){
                if( is_neumann(triangle.P2()[0], triangle.P2()[1], triangle.P2()[2]) ){
                    if( is_neumann(triangle.P3()[0], triangle.P3()[1], triangle.P3()[2]) ){
                        mNeumannTriangles.push_back(triangle);
                    }
                }
            }
        }

        return mNeumannTriangles;
    }


    TriangleVectorType& GetDirichletTriangles(std::function<bool(double, double,double)> &is_dirichlet){

        auto& triangles = this->GetTriangles();

        mDirichletTriangles.clear();
        for( auto triangle : triangles){
            if( is_dirichlet(triangle.P1()[0], triangle.P1()[1], triangle.P1()[2]) ){
                if( is_dirichlet(triangle.P2()[0], triangle.P2()[1], triangle.P2()[2]) ){
                    if( is_dirichlet(triangle.P3()[0], triangle.P3()[1], triangle.P3()[2]) ){
                        mDirichletTriangles.push_back(triangle);
                    }
                }
            }
        }

        return mDirichletTriangles;
    }

    PointType GetLocalUpperPoint(){
        return mLocalUpperPoint;
    }

    PointType GetLocalLowerPoint(){
        return mLocalLowerPoint;
    }

    IntegrationPoint1DVectorType& IntegrationPoints1D(int i){
        if(i ==0)
            return mIntegrationPointsX;
        if(i==1)
            return mIntegrationPointsY;

        return mIntegrationPointsZ;
    }

    IntegrationPoint1DVectorType& IntegrationPointsY(){
        return mIntegrationPointsY;
    }
    IntegrationPoint1DVectorType& IntegrationPointsZ(){
        return mIntegrationPointsZ;
    }

    void pSetSurfaceMesh(SurfaceMeshPtrType& pSurfaceMesh );

    SurfaceMeshPtrType pGetSurfaceMesh();

    bool IsPointInTrimmedDomain(PointType& rTestPoint);

    BoundingBox ComputeTrimmedBoundingBox();

    void ClearSurfaceMesh();

    void SetNeighbourCoefficient(double value, int direction){
        mNumberOfNeighbours[direction] = value;
    }

    double NeighbourCoefficient(){
        return mNumberOfNeighbours[0]*mNumberOfNeighbours[1]*mNumberOfNeighbours[2];
    }

    void SetVisited(bool value){
        mIsVisited = value;
    }

    bool IsVisited(){
        return mIsVisited;
    }

    // Information to keep track of iterations of point elimination algorithm. Only for publication.
    void AddResidual(double value){
        mResidual.push_back(value);
    }
    void AddNumberIps(int value){
        mNumberIntegrationPoints.push_back(value);
    }
    std::vector<double> MfIterationsResidual(){
        return mResidual;
    }
    std::vector<int> MfIterationsPoints(){
        return mNumberIntegrationPoints;
    }

private:
    IntegrationPointPtrVectorType mIntegrationPointsTrimmed;
    IntegrationPointPtrVectorType mIntegrationPointsInside;
    IntegrationPointPtrVectorType mIntegrationPointsFictitious;

    TriangleVectorType mTriangles{};
    TriangleVectorType mNeumannTriangles;
    TriangleVectorType mDirichletTriangles;

    const Parameters& mParameters;
    PointType mLocalUpperPoint;
    PointType mLocalLowerPoint;

    IntegrationPoint1DVectorType mIntegrationPointsX;
    IntegrationPoint1DVectorType mIntegrationPointsY;
    IntegrationPoint1DVectorType mIntegrationPointsZ;

    SurfaceMeshPtrType mpSurfaceMesh;
    bool mSurfaceMeshSetFlag;

    std::unique_ptr<InsideTestType> mpInsideTest;

    std::size_t mElementId;
    bool mIsTrimmed;

    std::array<double, 3> mNumberOfNeighbours{};
    bool mIsVisited{};

    // Information to keep track of iterations of point elimination algorithm. Only for publication.
    std::vector<double> mResidual{};
    std::vector<int> mNumberIntegrationPoints{};
};

#endif