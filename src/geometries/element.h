// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_H
#define ELEMENT_H

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

// Project includes
#include "geometries/integration_point.h"
#include "geometries/triangle_3d_3n.h"
#include "utilities/parameters.h"
#include "utilities/mapping_utilities.h"

// External includes
#include "omp.h"
#include <memory>


class Element
{

public:

    // Typedefs
    typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
    typedef std::vector<std::array<double, 2>> IntegrationPoint1DVectorType;
    typedef std::vector<Triangle3D3N> TriangleVectorType;
    typedef std::array<double,3> PointType;
    typedef std::size_t IndexType;

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef K::Vector_3 Vector;
    typedef CGAL::Surface_mesh<Point_3> SurfaceMeshType;
    typedef std::unique_ptr<SurfaceMeshType> SurfaceMeshPtrType;
    typedef CGAL::Side_of_triangle_mesh<SurfaceMeshType, K> InsideTestType;
    typedef std::vector<std::array<double,3>> BoundingBox;
    typedef typename boost::graph_traits<CGAL::Surface_mesh<Point_3> >::vertex_descriptor vertex_descriptor;
    typedef CGAL::Surface_mesh<Point_3>::Property_map<vertex_descriptor, CGAL::Surface_mesh<Point_3>::Point> PositionType;

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
        return mpSurfaceMesh->number_of_faces();
    }

    const Parameters& GetParameters() const {
        return mParameters;
    }

    IntegrationPointVectorType& GetIntegrationPointsTrimmed(){
        return mIntegrationPointsTrimmed;
    }
    IntegrationPointVectorType& GetIntegrationPointsInside(){
        return mIntegrationPointsInside;
    }
    IntegrationPointVectorType& GetIntegrationPointsFictitious(){
        return mIntegrationPointsFictitious;
    }
    TriangleVectorType& GetTriangles(){
        if( !mSurfaceMeshSetFlag ){
            throw  std::runtime_error("Element :: Surface Mesh Pointer has not been set" );
        }
        mTriangles.reserve(mpSurfaceMesh->number_of_faces());

        for(auto fd : faces(*mpSurfaceMesh))
        {
            Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(fd, *mpSurfaceMesh);
            std::array<double,3> tmp_normal = {n[0], n[1], n[2]};

            std::vector<std::array<double,3>> tmp_coordinates(3);

            PositionType positions = mpSurfaceMesh->points();
            int index = 0;
            for( auto vh : vertices_around_face(halfedge(fd,*mpSurfaceMesh),*mpSurfaceMesh) ){
            std::array<double,3> tmp_array; // = {p.x(), p.y(), p.z()};
            tmp_array[0] = positions[vh].x();
            tmp_array[1] = positions[vh].y();
            tmp_array[2] = positions[vh].z();
            tmp_coordinates[index] = tmp_array;
            index++;
            }
            if( tmp_coordinates.size() == 3){
                mTriangles.push_back(Triangle3D3N(tmp_coordinates[0], tmp_coordinates[1], tmp_coordinates[2], tmp_normal));
            }
            else {
                throw std::runtime_error("Facet is not a triangle!");
            }
        }

        this->ClearSurfaceMesh();

        return mTriangles;
    }

    TriangleVectorType& GetNeumannTriangles(std::function<bool(double, double,double)> &is_neumann){
        if( !mSurfaceMeshSetFlag && mTriangles.size() == 0){
            throw  std::runtime_error("Element :: No boundary mesh stored." );
        }
        mNeumannTriangles.clear();
        if( !mSurfaceMeshSetFlag ){
            for( auto triangle : mTriangles){
                if( is_neumann(triangle.P1()[0], triangle.P1()[1], triangle.P1()[2]) ){
                    if( is_neumann(triangle.P2()[0], triangle.P2()[1], triangle.P2()[2]) ){
                        if( is_neumann(triangle.P3()[0], triangle.P3()[1], triangle.P3()[2]) ){
                            mNeumannTriangles.push_back(triangle);
                        }
                    }
                }
            }

        }
        else {
            for(auto fd : faces(*mpSurfaceMesh)) {
                Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(fd, *mpSurfaceMesh);
                std::array<double,3> tmp_normal = {n[0], n[1], n[2]};

                std::vector<std::array<double,3>> tmp_coordinates(3);
                PositionType positions = mpSurfaceMesh->points();
                int index = 0;
                for( auto vh : vertices_around_face(halfedge(fd,*mpSurfaceMesh),*mpSurfaceMesh) ){
                    std::array<double,3> tmp_array;
                    tmp_array[0] = positions[vh].x();
                    tmp_array[1] = positions[vh].y();
                    tmp_array[2] = positions[vh].z();
                    tmp_coordinates[index] = tmp_array;
                    index++;
                }
                if( is_neumann(tmp_coordinates[0][0], tmp_coordinates[0][1], tmp_coordinates[0][2]) ){
                    if( is_neumann(tmp_coordinates[1][0], tmp_coordinates[1][1], tmp_coordinates[1][2]) ){
                        if( is_neumann(tmp_coordinates[2][0], tmp_coordinates[2][1], tmp_coordinates[2][2]) ){
                            mNeumannTriangles.push_back(Triangle3D3N(tmp_coordinates[0], tmp_coordinates[1], tmp_coordinates[2], tmp_normal));
                        }
                    }
                }
            }
        }

        return mNeumannTriangles;
    }

    TriangleVectorType& GetDirichletTriangles(std::function<bool(double, double,double)> &is_dirichlet){
        if( !mSurfaceMeshSetFlag && mTriangles.size() == 0){
            throw  std::runtime_error("Element :: No boundary mesh stored." );
        }
        mDirichletTriangles.clear();
        if( !mSurfaceMeshSetFlag ){
            for( auto triangle : mTriangles){
                if( is_dirichlet(triangle.P1()[0], triangle.P1()[1], triangle.P1()[2]) ){
                    if( is_dirichlet(triangle.P2()[0], triangle.P2()[1], triangle.P2()[2]) ){
                        if( is_dirichlet(triangle.P3()[0], triangle.P3()[1], triangle.P3()[2]) ){
                            mDirichletTriangles.push_back(triangle);
                        }
                    }
                }
            }
        }
        else {
            for(auto fd : faces(*mpSurfaceMesh)) {
                Vector n = CGAL::Polygon_mesh_processing::compute_face_normal(fd, *mpSurfaceMesh);
                std::array<double,3> tmp_normal = {n[0], n[1], n[2]};

                std::vector<std::array<double,3>> tmp_coordinates(3);
                PositionType positions = mpSurfaceMesh->points();
                int index = 0;
                for( auto vh : vertices_around_face(halfedge(fd,*mpSurfaceMesh),*mpSurfaceMesh) ){
                    std::array<double,3> tmp_array;
                    tmp_array[0] = positions[vh].x();
                    tmp_array[1] = positions[vh].y();
                    tmp_array[2] = positions[vh].z();
                    tmp_coordinates[index] = tmp_array;
                    index++;
                }
                if( is_dirichlet(tmp_coordinates[0][0], tmp_coordinates[0][1], tmp_coordinates[0][2]) ){
                    if( is_dirichlet(tmp_coordinates[1][0], tmp_coordinates[1][1], tmp_coordinates[1][2]) ){
                        if( is_dirichlet(tmp_coordinates[2][0], tmp_coordinates[2][1], tmp_coordinates[2][2]) ){
                            mDirichletTriangles.push_back(Triangle3D3N(tmp_coordinates[0], tmp_coordinates[1], tmp_coordinates[2], tmp_normal));
                        }
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

    PointType GetGlobalUpperPoint(){
        return MappingUtilities::FromLocalToGlobalSpace(mLocalUpperPoint, mParameters.PointA(), mParameters.PointB());
    }

    PointType GetGlobalLowerPoint(){
        return MappingUtilities::FromLocalToGlobalSpace(mLocalLowerPoint, mParameters.PointA(), mParameters.PointB());
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

    void pSetSurfaceMesh(SurfaceMeshPtrType& pSurfaceMesh ){
        mpSurfaceMesh = std::move(pSurfaceMesh);
        mpInsideTest = std::make_unique<InsideTestType>(*mpSurfaceMesh);
        mSurfaceMeshSetFlag = true;
    }

    SurfaceMeshType& GetSurfaceMesh(){
        if( !mSurfaceMeshSetFlag ){
            throw  std::runtime_error("Element :: Surface Mesh Pointer has not been set" );
        }
        return *mpSurfaceMesh;
    }

    void ClearSurfaceMesh(){
        mpInsideTest.reset();
        mpSurfaceMesh.reset();
        mSurfaceMeshSetFlag = false;
    }

    // TODO: Use inside_test.h and remove corresponding include.
    bool IsPointInTrimmedDomain(PointType& rTestPoint){
        if( !mSurfaceMeshSetFlag ){
            throw  std::runtime_error("Element :: Surface Mesh Pointer has not been set" );
        }
        Point_3 point_3(rTestPoint[0], rTestPoint[1], rTestPoint[2]);
        CGAL::Bounded_side res = (*mpInsideTest)(point_3);
        if (res == CGAL::ON_BOUNDED_SIDE) { return true; }
        if (res == CGAL::ON_BOUNDARY) { return false; }

        return false;
    }

    BoundingBox ComputeTrimmedBoundingBox(){
        if( !mSurfaceMeshSetFlag ){
            throw  std::runtime_error("Element :: Surface Mesh Pointer has not been set" );
        }
        BoundingBox bounding_box = { {1e10, 1e10, 1e10},
                                 {-1e10, -1e10, -1e10} };

        typedef SurfaceMeshType::Vertex_index vertex_descriptor;
        typedef typename boost::property_map<SurfaceMeshType, CGAL::vertex_point_t>::const_type VPMap;

        VPMap vpmap = get(CGAL::vertex_point, *mpSurfaceMesh);
        for(vertex_descriptor v : vertices(*mpSurfaceMesh))
        {
            const Point_3& p = get(vpmap, v);
            std::array<double, 3> tmp_coordinate = {CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())};
            // Loop over all 3 dimensions
            for( IndexType i = 0; i < 3; ++i){
                if( tmp_coordinate[i] < bounding_box[0][i] ){ // Find min values
                    bounding_box[0][i] = tmp_coordinate[i];
                }
                if( tmp_coordinate[i] > bounding_box[1][i] ){ // Find max values
                    bounding_box[1][i] = tmp_coordinate[i];
                }
            }
        }
        return bounding_box;
    }

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
    IntegrationPointVectorType mIntegrationPointsTrimmed;
    IntegrationPointVectorType mIntegrationPointsInside;
    IntegrationPointVectorType mIntegrationPointsFictitious;

    TriangleVectorType mTriangles;
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