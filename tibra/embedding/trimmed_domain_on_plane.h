// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef EDGE_2D_INCLUDE_H
#define EDGE_2D_INCLUDE_H

//// STL includes
#include <cstddef>
#include <array>
#include <algorithm>
//// Project includes
#include "containers/boundary_integration_point.h"
#include "io/io_utilities.h"

///Project includes

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  TrimmedDomainOnPlane
 * @author Manuel Messmer
 * @brief Protype: Projects the trimmed domain onto a plane of the AABB and provides functions to construct
 *        integration points on trimmed domain. Still requires validation. There are also still some bugs.
*/
class TrimmedDomainOnPlane {

public:
    ///@name Type Definitions
    ///@{
    typedef std::array<double, 2> Point2DType;
    typedef Vector3d Point3DType;
    typedef std::vector<BoundaryIntegrationPoint> BoundaryIPVectorType;
    typedef std::unique_ptr<BoundaryIPVectorType> BoundaryIPVectorPtrType;

    /**
     * @class  Edge2D
     * @author Manuel Messmer
     * @brief Simple edge with 2 vertices in 2D.
    */
    class Egde2D {

    public:
        ///@name Type Definitions
        ///@{

        Egde2D( IndexType V1, IndexType V2 ) :
            mV1(V1), mV2(V2), mIsVisited(false)
        {
        }

        ///@}
        ///@name Type Definitions
        ///@{
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
        ///@}

    private:
        ///@name Private Members
        ///@{
        IndexType mV1;
        IndexType mV2;
        std::vector<Point2DType> mSplitPoints{};
        bool mIsVisited;
        ///@}
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    TrimmedDomainOnPlane(IndexType PlaneIndex, bool UpperBoundary, const Point3DType& LowerBound, const Point3DType& UpperBound)
            : mUpperBoundary(UpperBoundary), mLowerBound(LowerBound), mUpperBound(UpperBound)
    {
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
            throw std::runtime_error("TrimmedDomainOnPlane :: Contructor :: Wrong PlaneIndex.");
        }
    }

    ///@}
    ///@name Public Operations
    ///@{

    void Reserve(IndexType Size){
        mEdgesPositiveOriented.reserve(Size);
        mEdgesNegativeOriented.reserve(Size);
        mVerticesPositive.reserve(Size);
        mVerticesNegative.reserve(Size);
    }

    ///@brief Insert new edges to trimmed domain on plane. Check for dublicate nodes. Only inserts new node if unique.
    ///       Only inserts edges if std::abs(normal[DIRINDEX2]) > 1e-10. All other edges have no contribution for constant terms
    ///       of moment fitting equation.
    ///@param rEdged
    ///@param rNormal
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

    ///@brief Returns boundary points of trimmed domain on plane.
    ///@param [out] pIPs
    ///@param state_a_c Must be true if point a is inside.
    ///@param state_b_d Must be true if point a is inside.
    ///@note Param mapping:
    //
    //     a_______b                 y
    //     /      /|                ´|`
    //  c /_____d/ |<-- plane lower  |-->x
    //    |     |  /                /
    //    |     |</-- plane upper  Z
    //    |_____|/
    //
    void pGetBoundaryIPs(BoundaryIPVectorPtrType& pIPs, bool state_a_c, bool state_b_d){
        RefineEdges(state_a_c, state_b_d);
        ConstructIPS(pIPs);
    }
    ///@}

private:
    ///@name Public Operations
    ///@{

    ///@brief Inserts point to container if point is not contained in container.
    ///@param rPoint NewPoint
    ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
    ///@param Return vertex index in container.
    IndexType InsertUniqueVertex(const Point2DType& rPoint, bool Positive){

        auto& r_vertices = GetVertices(Positive);
        // Find new vertex.
        auto res = std::find_if( r_vertices.begin(), r_vertices.end(),
            [&rPoint](const auto& x) { return  (rPoint[0] > x[0] - 1e-8 && rPoint[0] < x[0] + 1e-8
                                                && rPoint[1] > x[1] - 1e-8 && rPoint[1] < x[1] + 1e-8) ;});

        if( res != r_vertices.end() ){ // Vertex already exists
            IndexType index = std::distance(r_vertices.begin(), res);
            return index;
        } else { // Add new vertex
            IndexType index = r_vertices.size();
            r_vertices.push_back(rPoint);
            return index;
        }
    }

    ///@brief Fast insertion of point to container. Check if point is already contained is omittede.
    ///@param rPoint NewPoint
    ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
    ///@param Return vertex index in container.
    IndexType InsertVertex(const Point2DType& rPoint, bool Positive){

        auto& r_vertices = GetVertices(Positive);
        IndexType index = r_vertices.size();
        r_vertices.push_back(rPoint);
        return index;
    }

    ///@brief Get Vertex V1 by EdgeId
    ///@param EdgeId
    ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
    const Point2DType& V1byEdgeId(IndexType EdgeId,  bool Positive){
        if( Positive ){
            return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V1()];
        }
        else {
            return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V1()];
        }
    }

    ///@brief Refine edges on trimmed domain on plane. See "details" for more information.
    ///@param state_a Positive if point a is inside domain.
    ///@param state_b Positive if point b is inside domain.
    ///@note See pGetBoundaryIPs for mapping/indexing of points a and b.
    ///@details 1.) Maps all negative oriented points onto positive oriented edges if overlap and vice versa.
    ///             Such that the vertices of postive and negative oriented edges are at same position in: DIRINDEX1.
    ///                 x--^-----x---^--x   positive oriented edges
    ///                 |  |     |   |                                     DIRINDEX2
    ///                 |  |     |   |                                         ^
    ///            x----v--x-----v---x      negative oriented edges            |---> DIRINDEX1
    ///
    ///         2.) Checks if positive oriented edges span the entiere AABB (we consider only the part that is inside the domain) in DIRINDEX1.
    ///             If not new edges are introduced at top of AABB (UpperBound[DIRINDEX2]) to fill gaps.
    ///         a             b      a        _____b
    ///         |  out  /     |      |  out  /     |  new edges introduced since Point b is inside domain.
    ///         |------/      | ---> |------/      |
    ///         |  inside     |      |  inside     |
    ///        LB             UB     LB            UP
    ///
    ///         |  in   /     |      |  in   /     |  no new edges introduced since Point b is outside.
    ///         |------/      | ---> |------/      |       LB - lower bound of AABB in DIRINDEX1
    ///         |  outsid     |      |  outside    |       UP - upper bound of AABB in DIRINDEX1
    ///        LB             UB     LB            UP
    ///
    ///@todo    Not enough, we need to find all intersection on top plane and then shoot only downwards.
    void RefineEdges(bool state_a, bool state_b){
        IndexType num_positive_edges = mEdgesPositiveOriented.size();
        IndexType num_negative_edges = mEdgesNegativeOriented.size();


        if( num_positive_edges == 0 ||  num_negative_edges == 0 ){
            return;
        }

        /// Split negative edges at positions of positive points.
        for( int i = 0; i < mVerticesPositive.size(); ++i){
            const auto& v = mVerticesPositive[i];
            SetSplitPoint(v, mEdgesNegativeOriented, true);
        }
        SplitEdgesAtSplitPoint(mEdgesNegativeOriented, false);


        /// Split positive edges at positions of negative points.
        for( int i = 0; i < mVerticesNegative.size(); ++i){
            const auto& v = mVerticesNegative[i];
            SetSplitPoint(v, mEdgesPositiveOriented, false);
        }
        SplitEdgesAtSplitPoint(mEdgesPositiveOriented, true);

        ///TODO: Not complete yet. Find all intersection with top edge from a-b.
        // Search from left boundary
        if( state_a ){
            const double corner_point_x = mLowerBound[DIRINDEX1];
            double min_distance = 1e10;
            for(int i = 0; i < mVerticesPositive.size(); ++i ){
                double distance = std::abs(mVerticesPositive[i][0] - corner_point_x);
                if( distance < min_distance ){
                    min_distance = distance;
                }
            }
            if( min_distance > 1e-8 ){
                double point_y = mUpperBound[DIRINDEX2];
                const IndexType num_edges = 5;
                double delta_x = min_distance/num_edges;
                for(IndexType i = 0; i < 5; ++i ){
                    double point_x_1 = corner_point_x + i*delta_x;
                    IndexType index1 = InsertVertex( {point_x_1, point_y}, true );
                    double point_x_2 = corner_point_x + (i+1)*delta_x;
                    IndexType index2 = InsertVertex( {point_x_2, point_y}, true );
                    mEdgesPositiveOriented.push_back( Egde2D( index1, index2) );
                }
            }
        }

        // Search from right boundary
        if( state_b ){
            const double corner_point_x = mUpperBound[DIRINDEX1];
            double min_distance = 1e10;
            for(int i = 0; i < mVerticesPositive.size(); ++i ){
                double distance = std::abs(mVerticesPositive[i][0] - corner_point_x);
                if( distance < min_distance ){
                    min_distance = distance;
                }
            }
            if( min_distance > 1e-8 ){
                double point_y = mUpperBound[DIRINDEX2];
                const IndexType num_edges = 5;
                double delta_x = min_distance/num_edges;
                for(IndexType i = 0; i < 5; ++i ){
                    double point_x_1 = corner_point_x - i*delta_x;
                    IndexType index1 = InsertVertex( {point_x_1, point_y}, true );
                    double point_x_2 = corner_point_x - (i+1)*delta_x;
                    IndexType index2 = InsertVertex( {point_x_2, point_y}, true );
                    mEdgesPositiveOriented.push_back( Egde2D( index1, index2) );
                }
            }
        }
    }

    ///@brief Get Vertex V2 by EdgeId
    ///@param EdgeId
    ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
    const Point2DType& V2byEdgeId(IndexType EdgeId,  bool Positive){
        if( Positive ){
            return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V2()];
        }
        else {
            return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V2()];
        }
    }

    ///@brief Return EdgeId that overlaps rPoint in DIRINDEX1 direction. Return -1 if no edge is found.
    ///       Found if: rPoint[DIRINDEX1] > EgdeV1[DIRINDEX1] and rPoint[DIRINDEX1] < EgdeV2[DIRINDEX1]
    ///@param rPoint
    ///@param Positive Orientation of edges to be searched.
    ///@param int
    int FindIntersectingEdge(const Point2DType& rPoint, bool Positive){
        for( int i = 0; i < GetNumberEdges(Positive); ++i ){
            if( rPoint[DIRINDEX1] > V1byEdgeId(i,Positive)[DIRINDEX1] && rPoint[DIRINDEX1] < V2byEdgeId(i, Positive)[DIRINDEX1] ){
                return i;
            }
        }
        return -1;
    }

    ///@brief Create integration points by ray shooting.
    ///@details Shoot ray from each positive oriented edge downards in negative DIRINDEX2 direction. End point is defined by intersecting
    ///         negative oriented edge or mLowerBound[DIRINDEX2].
    ///@param [out] pIPs (new points are pushed_back(). Old points are preserved.)
    ///@param PlanePosition Distance to origin in DIRINDEX3 direction.
    ///@param rNormal Normal vector of plane
    void ConstructIPSbyRayShooting( BoundaryIPVectorPtrType& pIPs, double PlanePosition,
                const Point3DType& rNormal ){

        auto& r_edges_origin = mEdgesPositiveOriented;
        auto& r_edges_dest = mEdgesNegativeOriented;

        bool orientation_origin = true;
        bool orientation_dest = false;

        // First we shoot rays from upper edges downwards.
        for( int i = 0; i < r_edges_origin.size(); ++i){

            //TODO: I think visited is not required.
            //TODO: Improve this.
            if( !r_edges_origin[i].IsVisited() ) {

                r_edges_origin[i].SetVisited(true);

                const auto& v1_up = V1byEdgeId(i, orientation_origin);
                const auto& v2_up = V2byEdgeId(i, orientation_origin);

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
                    const auto& v1_low = V1byEdgeId(edge_id_dest, false);
                    const auto& v2_low = V2byEdgeId(edge_id_dest, false);
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

    ///@brief Constructs normal vector and plane position and calls ConstructIPSbyRayShooting().
    ///@param [out] pIPs
    void ConstructIPS(BoundaryIPVectorPtrType& pIPs){
        Point3DType normal = {0.0, 0.0, 0.0};
        double plane_position = 0.0;
        if( mUpperBoundary ){
            normal[DIRINDEX3] = 1.0;
            plane_position = mUpperBound[DIRINDEX3];
        } else {
            normal[DIRINDEX3] = -1.0;
            plane_position = mLowerBound[DIRINDEX3];
        }

        ConstructIPSbyRayShooting(pIPs, plane_position, normal);
    }

    ///@brief Splits all edges at their defined split point.
    ///@param rEdges to be split.
    ///@param Positive Orientation of edges,
    ///@todo See TODO in definition.
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

                IndexType index1 = InsertVertex( V1byEdgeId(i, Positive), Positive );
                IndexType index2 = InsertVertex( split_points[0], Positive );
                new_edges.push_back( Egde2D( index1, index2) );

                for( int j = 0; j < num_split_points-1; ++j){
                    index1 = InsertVertex( split_points[j], Positive );
                    index2 = InsertVertex( split_points[j+1], Positive );
                    new_edges.push_back( Egde2D( index1, index2) );
                }

                index1 = InsertVertex( split_points[num_split_points-1], Positive );
                index2 = InsertVertex(  V2byEdgeId(i, Positive), Positive );
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

    ///@brief Find intersecting point on Edge with x=Point[DIRINDEX1] and mark as split point.
    ///@param rPoint Potential split point.
    ///@param rEdges Edges to be searched and marked.
    ///@param Positive Orientation of edges,
    void SetSplitPoint(const Point2DType& rPoint, std::vector<Egde2D>& rEdges, bool Positive){
        // Take bool from Edge.
        double current_distance = 1e15;
        int edge_id = -1;
        Point2DType intersection_point{};
        for( IndexType i = 0; i < rEdges.size(); ++i){
            const auto& v1 = V1byEdgeId(i, !Positive);
            const auto& v2 = V2byEdgeId(i, !Positive);

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


    ///@brief Returns Edges container.  (const version)
    ///@param Positive orientation of edges
    ///@return const std::vector<Egde2D>&
    const std::vector<Egde2D>& GetEdges(bool Positive) const{
        if( Positive )
            return mEdgesPositiveOriented;
        else
            return mEdgesNegativeOriented;
    }

    ///@brief Returns Edges container.  (non const version)
    ///@param Positive orientation of edges
    ///@return std::vector<Egde2D>&
    std::vector<Egde2D>& GetEdges(bool Positive){
        if( Positive )
            return mEdgesPositiveOriented;
        else
            return mEdgesNegativeOriented;
    }

    ///@brief Returns vertices container.  (const version)
    ///@param Positive orientation of vertices
    ///@return const std::vector<Point2DType>&
    const std::vector<Point2DType>& GetVertices(bool Positive) const{
        if( Positive )
            return mVerticesPositive;
        else
            return mVerticesNegative;
    }

    ///@brief Returns vertices container.  (non-const version)
    ///@param Positive orientation of vertices
    ///@return std::vector<Point2DType>&
    std::vector<Point2DType>& GetVertices(bool Positive){
        if( Positive )
            return mVerticesPositive;
        else
            return mVerticesNegative;
    }

    ///@brief Returns numver of edges.
    ///@param Positive Orientation of edges.
    ///@return IndexType
    IndexType GetNumberEdges(bool Positive){
        if( Positive )
            return mEdgesPositiveOriented.size();
        else
            return mEdgesNegativeOriented.size();
    }

    ///@}
    ///@name Private Members
    ///@{

    std::vector<Egde2D> mEdgesPositiveOriented{};
    std::vector<Egde2D> mEdgesNegativeOriented{};
    std::vector<Point2DType> mVerticesPositive{};
    std::vector<Point2DType> mVerticesNegative{};

    IndexType DIRINDEX1; // In plane 1
    IndexType DIRINDEX2; // In plane 2
    IndexType DIRINDEX3; // Out of plane

    Point3DType mLowerBound; // Lower bound AABB
    Point3DType mUpperBound; // Upper bound AABB

    bool mUpperBoundary; //Is current plane upper bound?
    ///@}
}; // End TrimmedDomainOnPlane

///@} // End Tibra Classes

} // End namespace tibra
#endif