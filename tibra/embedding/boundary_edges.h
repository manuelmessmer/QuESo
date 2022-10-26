// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef EDGE_2D_INCLUDE_H
#define EDGE_2D_INCLUDE_H

/// External libraries
#include <cstddef>
#include <array>
#include <algorithm>

#include "geometries/boundary_integration_point.h"
#include "io/io_utilities.h"

///Project includes

///@name TIBRA Classes
///@{

/**
 * @class  BoundaryEdges
 * @author Manuel Messmer
 * @brief
*/
class BoundaryEdges {

public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;
    typedef std::array<double, 2> Point2DType;
    typedef std::array<double, 3> Point3DType;
    typedef std::vector<BoundaryIntegrationPoint> BoundaryIPVectorType;
    typedef std::unique_ptr<BoundaryIPVectorType> BoundaryIPVectorPtrType;

    /// Default constructor.
    BoundaryEdges(IndexType PlaneIndex, const Point3DType& LowerBound, const Point3DType& UpperBound) : mLowerBound(LowerBound), mUpperBound(UpperBound) {
        if( PlaneIndex == 2 ){
            DIRINDEX1 = 0;
            DIRINDEX2 = 1;
            DIRINDEX3 = 2;
        }
        else if( PlaneIndex == 1 ){
            DIRINDEX1 = 0;
            DIRINDEX2 = 2;
            DIRINDEX3 = 1;
        }
        else if( PlaneIndex == 0 ){
            DIRINDEX1 = 1;
            DIRINDEX2 = 2;
            DIRINDEX3 = 0;
        }
        else {
            throw std::runtime_error("BoundaryEdges :: Contructor :: Wrong PlaneIndex.");
        }
    }
    /**
     * @class  Edge
     * @author Manuel Messmer
     * @brief
    */
    class Egde2D {

    public:
        ///@name Type Definitions
        ///@{

        Egde2D( IndexType V1, IndexType V2 ) :
            mV1(V1), mV2(V2), mIsVisited(false)
        {
        }

        IndexType V1() const{
            return mV1;
        }
        IndexType V2() const{
            return mV2;
        }

        std::vector<Point2DType>& GetSplitPoints()  {
            return mSplitPoints;
        }

        const std::vector<Point2DType>& GetSplitPoints() const {
            return mSplitPoints;
        }

        void AddSplitPoint(const Point2DType& rPoint){
            mSplitPoints.push_back(rPoint);
        }

        bool IsVisited(){
            return mIsVisited;
        }

        void SetVisited(bool Value){
            mIsVisited = Value;
        }

    private:
        IndexType mV1;
        IndexType mV2;
        std::vector<Point2DType> mSplitPoints{};
        bool mIsVisited;
    };

    ///@}
    void Reserve(IndexType Size){
        mEdgesPositiveOriented.reserve(Size);
        mEdgesNegativeOriented.reserve(Size);
    }

    void InsertEdges( const std::vector<std::array<Point3DType,2>> rEdges, const Point3DType& rNormal ) {
        for( auto& edge : rEdges ){
            // Get vertices of edge
            const auto& V1 = edge[0];
            const auto& V2 = edge[1];

            Point2DType v1{};
            Point2DType v2{};
            // Make sure edge is oriented along DIRINDEX1 such that x2 > x1
            if( V2[DIRINDEX1] > V1[DIRINDEX1] ){
                v1 = {V1[DIRINDEX1], V1[DIRINDEX2]};
                v2 = {V2[DIRINDEX1], V2[DIRINDEX2]};
            }
            else {
                v1 = {V2[DIRINDEX1], V2[DIRINDEX2]};
                v2 = {V1[DIRINDEX1], V1[DIRINDEX2]};
            }
            // Positive oriented
            if( rNormal[DIRINDEX2] > 1e-10 ){
                // Get unique vertex index.
                auto index_1 = InsertUniqueVertex(v1, true);
                auto index_2 = InsertUniqueVertex(v2, true);
                mEdgesPositiveOriented.push_back( Egde2D(index_1, index_2) );
            } // Negative oriented
            else if ( rNormal[DIRINDEX2] < -1e-10 ) {
                // Get unique vertex index.
                auto index_1 = InsertUniqueVertex(v1, false);
                auto index_2 = InsertUniqueVertex(v2, false);
                mEdgesNegativeOriented.push_back( Egde2D(index_1, index_2) );
            }
        }
        // Do not add vertical edges.
    }

    IndexType InsertUniqueVertex(const Point2DType& rPoint, bool Positive){
        // TODO: Improve this.
        if( Positive ){ // Positive insertion
            auto res = std::find_if( mVerticesPositive.begin(), mVerticesPositive.end(),
                [&rPoint](const auto& x) { return  (rPoint[0] > x[0] - 1e-8 && rPoint[0] < x[0] + 1e-8
                                                 && rPoint[1] > x[1] - 1e-8 && rPoint[1] < x[1] + 1e-8) ;});

            if( res != mVerticesPositive.end() ){
                IndexType index = std::distance(mVerticesPositive.begin(), res);
                return index;
            } else {
                IndexType index = mVerticesPositive.size();
                mVerticesPositive.push_back(rPoint);
                return index;
            }
        } else { // Negative insertion
            auto res = std::find_if( mVerticesNegative.begin(), mVerticesNegative.end(),
                [&rPoint](const auto& x) { return  (rPoint[0] > x[0] - 1e-8 && rPoint[0] < x[0] + 1e-8
                                                 && rPoint[1] > x[1] - 1e-8 && rPoint[1] < x[1] + 1e-8) ;});

            if( res != mVerticesNegative.end() ){
                IndexType index = std::distance(mVerticesNegative.begin(), res);
                return index;
            } else {
                IndexType index = mVerticesNegative.size();
                mVerticesNegative.push_back(rPoint);
                return index;
            }
        }
    }

    IndexType InsertVertex(const Point2DType& rPoint, bool Positive){
        //TODO: Improve this.
        if( Positive ){ // Positive insertion
            IndexType index = mVerticesPositive.size();
            mVerticesPositive.push_back(rPoint);
            return index;
        } else { // Negative insertion
            IndexType index = mVerticesNegative.size();
            mVerticesNegative.push_back(rPoint);
            return index;
        }
    }

    const Point2DType& V1(IndexType i,  bool Positive){
        if( Positive ){
            return mVerticesPositive[mEdgesPositiveOriented[i].V1()];
        }
        else {
            return mVerticesNegative[mEdgesNegativeOriented[i].V1()];
        }
    }

    const Point2DType& V2(IndexType i,  bool Positive){
        if( Positive ){
            return mVerticesPositive[mEdgesPositiveOriented[i].V2()];
        }
        else {
            return mVerticesNegative[mEdgesNegativeOriented[i].V2()];
        }
    }

    int FindIntersectingEdge(const Point2DType& rPoint, bool Positive){
        for( int i = 0; i < GetNumberEdges(Positive); ++i ){
            if( rPoint[DIRINDEX1] > V1(i,Positive)[DIRINDEX1] && rPoint[DIRINDEX1] < V2(i, Positive)[DIRINDEX1] ){
                return i;
            }
        }

        return -1;
    }

    void ConstructIPSbyRayShooting( BoundaryIPVectorPtrType& pIPs, double PlanePosition,
                const Point3DType& rNormal, bool ShootDownwards ){

        auto& r_edges_origin = ShootDownwards ? mEdgesPositiveOriented : mEdgesNegativeOriented;
        auto& r_edges_dest = ShootDownwards ? mEdgesNegativeOriented : mEdgesPositiveOriented;

        bool orientation_origin = ShootDownwards;
        bool orientation_dest = !ShootDownwards;

        // First we shoot rays from upper edges downwards.
        for( int i = 0; i < r_edges_origin.size(); ++i){

            if( !r_edges_origin[i].IsVisited() ) {

                const auto& v1_up = V1(i, orientation_origin);
                const auto& v2_up = V2(i, orientation_origin);

                // Get center
                Point2DType c_positive = { 0.5*(v1_up[0]+v2_up[0]), 0.5*(v1_up[1]+v2_up[1]) };
                int edge_id_dest = FindIntersectingEdge(c_positive, orientation_dest);

                double ray_width = std::abs(v2_up[0] - v1_up[0]);
                double ray_x = c_positive[0];
                double ray_y_start = std::min(v1_up[1], v2_up[1]);
                double ray_y_end = mLowerBound[DIRINDEX2];

                if( std::abs(v1_up[1] - v2_up[1]) > 1e-8 ){ // If edge is inclined
                    double weight = 0.5*std::abs( (v2_up[0]-v1_up[0]) * (v2_up[1]-v1_up[1]) );
                    Point3DType tmp_point = {0.0, 0.0, 0.0};
                    tmp_point[DIRINDEX3] = PlanePosition;

                    if( v1_up[1] > v2_up[1] ){

                        tmp_point[DIRINDEX1] = 1.0/3.0 * (2*v1_up[0]+v2_up[0] );
                        tmp_point[DIRINDEX2] = 1.0/3.0 * (v1_up[1]+2*v2_up[1] );
                        pIPs->push_back( BoundaryIntegrationPoint(tmp_point[0], tmp_point[1], tmp_point[2], weight, rNormal) );
                    }
                    else {
                        tmp_point[DIRINDEX1] = 1.0/3.0 * (v1_up[0]+2*v2_up[0] );
                        tmp_point[DIRINDEX2] = 1.0/3.0 * (2*v1_up[1]+v2_up[1] );
                        pIPs->push_back( BoundaryIntegrationPoint(tmp_point[0], tmp_point[1], tmp_point[2], weight, rNormal) );
                    }
                }

                if( edge_id_dest > -1){ // Lower edge is found.
                    r_edges_dest[edge_id_dest].SetVisited(true);
                    const auto& v1_low = V1(edge_id_dest, false);
                    const auto& v2_low = V2(edge_id_dest, false);
                    ray_y_end = std::max(v1_low[1], v2_low[1] );
                    if( std::abs(v1_low[1] - v2_low[1]) > 1e-8 ){ // If edge is inclined
                        double weight = 0.5*std::abs( (v2_low[0]-v1_low[0]) * (v2_low[1]-v1_low[1]) );
                        Point3DType tmp_point = {0.0, 0.0, 0.0};
                        tmp_point[DIRINDEX3] = PlanePosition;

                        if( v1_low[1] > v2_low[1] ){

                            tmp_point[DIRINDEX1] = 1.0/3.0 * (v1_low[0]+2*v2_low[0] );
                            tmp_point[DIRINDEX2] = 1.0/3.0 * (2*v1_low[1]+v2_low[1] );
                            pIPs->push_back( BoundaryIntegrationPoint(tmp_point[0], tmp_point[1], tmp_point[2], weight, rNormal) );
                        }
                        else {
                            tmp_point[DIRINDEX1] = 1.0/3.0 * (2*v1_low[0]+v2_low[0] );
                            tmp_point[DIRINDEX2] = 1.0/3.0 * (v1_low[1]+2*v2_low[1] );
                            pIPs->push_back( BoundaryIntegrationPoint(tmp_point[0], tmp_point[1], tmp_point[2], weight, rNormal) );
                        }
                    }
                }

                double delta_y = std::abs((ray_y_end-ray_y_start))/ (10.0);
                double weight = delta_y*ray_width;
                Point3DType tmp_point = {0.0, 0.0, 0.0};
                tmp_point[DIRINDEX3] = PlanePosition;
                double ray_y = ray_y_start - 0.5*delta_y;
                for( int i = 0; i < 10.0; ++i){
                    tmp_point[DIRINDEX1] = ray_x;
                    tmp_point[DIRINDEX2] = ray_y;
                    pIPs->push_back( BoundaryIntegrationPoint(tmp_point[0], tmp_point[1], tmp_point[2], weight, rNormal) );
                    ray_y -= delta_y;
                }
            } // End if( !r_edges_origin[i].IsVisited() )
        } // for( int i = 0; i < r_edges_origin.size(); ++i){
    }

    BoundaryIPVectorPtrType ConstructIPS(bool upper_bound){
        auto p_boundary_ips = std::make_unique<BoundaryIPVectorType>();

        Point3DType normal = {0.0, 0.0, 0.0};
        double plane_position = 0.0;
        if( upper_bound ){
            normal[DIRINDEX3] = 1.0;
            plane_position = mUpperBound[DIRINDEX3];
        } else {
            normal[DIRINDEX3] = -1.0;
            plane_position = mLowerBound[DIRINDEX3];
        }

        bool shoot_downwards = true;
        ConstructIPSbyRayShooting(p_boundary_ips, plane_position, normal, shoot_downwards);

        // When unvisited negative ones have to shoot upwards.

        IO::WritePointsToVTK(p_boundary_ips, "test_points.vtk", true);
        return std::move(p_boundary_ips);
    }

    BoundaryIPVectorPtrType pGetBoundaryIPs(){
        RefineEdges();

        // Actually refine edges.
        return ConstructIPS(false);
    }

    void RefineEdges(){
        IndexType num_positive_edges = mEdgesPositiveOriented.size();
        IndexType num_negative_edges = mEdgesNegativeOriented.size();

        // std::cout << "positive: " << num_positive_edges << std::endl;
        // std::cout << "negative: " << num_negative_edges << std::endl;

        if( num_positive_edges == 0 ||  num_negative_edges == 0 ){
            return;
        }
        //std::cout << "00000000000000000: " << std::endl;
        // Must ne verticess!!

        for( int i = 0; i < mVerticesPositive.size(); ++i){
            const auto& v = mVerticesPositive[i];
            SetSplitPoint(v, mEdgesNegativeOriented, true);
        }

        for( int i = 0; i < mVerticesNegative.size(); ++i){
            const auto& v = mVerticesNegative[i];
            SetSplitPoint(v, mEdgesPositiveOriented, false);
        }

        //std::cout << "before: " << std::endl;
        SplitEdgesAtSplitPoint(mEdgesPositiveOriented, true);
        SplitEdgesAtSplitPoint(mEdgesNegativeOriented, false);
        //std::cout << "after: " << std::endl;

        // std::cout << "positive: " << mEdgesPositiveOriented.size() << std::endl;
        // std::cout << "negative: " << mEdgesNegativeOriented.size() << std::endl;
    }

    void SplitEdgesAtSplitPoint(std::vector<Egde2D>& rEdges, bool Positive){
        std::vector<Egde2D> new_edges{};
        for( int i = 0; i < rEdges.size(); ++i){
            const auto& edge = rEdges[i];
            auto split_points = edge.GetSplitPoints();
            const IndexType num_split_points = split_points.size();
            if( num_split_points > 0 ){

                std::sort(split_points.begin(), split_points.end(), [](const auto& point_a, const auto& point_b) -> bool {
                    return point_a[0] > point_b[0];
                });

                IndexType index1 = InsertVertex( this->V1(i, Positive), Positive );
                IndexType index2 = InsertVertex( split_points[0], Positive );
                new_edges.push_back( Egde2D( index1, index2) );

                for( int j = 0; j < num_split_points-1; ++j){
                    index1 = InsertVertex( split_points[j], Positive );
                    index2 = InsertVertex( split_points[j+1], Positive );
                    new_edges.push_back( Egde2D( index1, index2) );
                }

                index1 = InsertVertex( split_points[num_split_points-1], Positive );
                index2 = InsertVertex(  this->V2(i, Positive), Positive );
                new_edges.push_back( Egde2D( index1, index2) );
            }
            else {
                new_edges.push_back(edge);
            }
        }

        //TODO: Make ptr and swap!!!
        rEdges.clear();
        rEdges.insert( rEdges.begin(), new_edges.begin(), new_edges.end() );
    }

    void SetSplitPoint(const Point2DType& rPoint, std::vector<Egde2D>& rEdges, bool Positive){
        // Take bool from Edge.
        double current_distance = 1e15;
        int edge_id = -1;
        Point2DType intersection_point{};
        for( IndexType i = 0; i < rEdges.size(); ++i){
            const auto& v1 = V1(i, !Positive);
            const auto& v2 = V2(i, !Positive);

            double pos1 = rPoint[0];
            double pos2 = 0.0;
            if( pos1 > v1[0]+1e-8 &&  pos1 < v2[0]-1e-8 ){
                // Note: x2 > x1. See AddEdge()
                // y1 + (y2-y1)/(x2-x1) * (x-x1)
                pos2 = v1[1] + (v2[1] - v1[1]) / (v2[0] - v1[0]) * (rPoint[0] - v1[0]);
                double distance = Positive ? (rPoint[1] - pos2) : (pos2-rPoint[1]);

                if( (distance < current_distance ) && distance > 0 ){
                    current_distance = pos2 - rPoint[0];
                    intersection_point[0] = pos1;
                    intersection_point[1] = pos2;
                    edge_id = i;
                }
            }
        }
        if( edge_id > -1){
            rEdges[edge_id].AddSplitPoint(intersection_point);
        }
    }


    const std::vector<Egde2D>& GetEdges(bool Positive){
        if( Positive )
            return mEdgesPositiveOriented;
        else
            return mEdgesNegativeOriented;
    }

    IndexType GetNumberEdges(bool Positive){
        if( Positive )
            return mEdgesPositiveOriented.size();
        else
            return mEdgesNegativeOriented.size();
    }

private:
    std::vector<Egde2D> mEdgesPositiveOriented{};
    std::vector<Egde2D> mEdgesNegativeOriented{};
    std::vector<Point2DType> mVerticesPositive{};
    std::vector<Point2DType> mVerticesNegative{};

    IndexType DIRINDEX1; // In plane
    IndexType DIRINDEX2; // In plane
    IndexType DIRINDEX3; // Out of plane

    Point3DType mLowerBound;
    Point3DType mUpperBound;
};

#endif