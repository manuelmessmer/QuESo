// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef EDGE_2D_INCLUDE_H
#define EDGE_2D_INCLUDE_H

//// STL includes
#include <cstddef>
#include <array>
#include <algorithm>
#include <iomanip>
//// Project includes
#include "containers/boundary_integration_point.h"
#include "utilities/tolerances.h"
#include "embedding/polygon.h"
#include "io/io_utilities.h"

/// Project includes

namespace tibra
{

    ///@name TIBRA Classes
    ///@{

    /**
     * @class  TrimmedDomainOnPlane
     * @author Manuel Messmer
     * @brief Protype: Projects the trimmed domain onto a plane of the AABB and provides functions to construct
     *        integration points on trimmed domain. Still requires validation. There are also still some bugs.
     */
    class TrimmedDomainOnPlane
    {

    public:
        ///@name Type Definitions
        ///@{
        typedef std::array<double, 2> Point2DType;
        typedef Vector3d Point3DType;
        typedef std::vector<BoundaryIntegrationPoint> BoundaryIPVectorType;
        typedef std::unique_ptr<BoundaryIPVectorType> BoundaryIPVectorPtrType;
        typedef std::unique_ptr<TriangleMesh> TriangleMeshPtrType;

        enum Orientation
        {
            Positive,
            Negative,
            Vertical
        };
        typedef enum Orientation OrientationType;
        /**
         * @class  Edge2D
         * @author Manuel Messmer
         * @brief Simple edge with 2 vertices in 2D.
         */
        class Egde2D
        {

        public:
            ///@name Type Definitions
            ///@{

            Egde2D(IndexType V1, IndexType V2) : mV1(V1), mV2(V2), mIsVisited(false)
            {
            }

            ///@}
            ///@name Type Definitions
            ///@{
            IndexType V1() const
            {
                return mV1;
            }
            IndexType V2() const
            {
                return mV2;
            }

            std::vector<Point2DType> &GetSplitPoints()
            {
                return mSplitPoints;
            }

            const std::vector<Point2DType> &GetSplitPoints() const
            {
                return mSplitPoints;
            }

            void ClearSplitPoints()
            {
                mSplitPoints.clear();
            }

            void AddSplitPoint(const Point2DType &rPoint)
            {
                mSplitPoints.push_back(rPoint);
            }

            bool IsVisited()
            {
                return mIsVisited;
            }

            void SetVisited(bool Value)
            {
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
        TrimmedDomainOnPlane(IndexType PlaneIndex, bool UpperBoundary, const Point3DType &LowerBound, const Point3DType &UpperBound, const TrimmedDomainBase *pTrimmedDomain)
            : mUpperBoundary(UpperBoundary), mLowerBound(LowerBound), mUpperBound(UpperBound), mpTrimmedDomain(pTrimmedDomain)
        {
            if (PlaneIndex == 2)
            {
                DIRINDEX1 = 0;
                DIRINDEX2 = 1;
                DIRINDEX3 = 2;
            }
            else if (PlaneIndex == 1)
            {
                DIRINDEX1 = 2;
                DIRINDEX2 = 0;
                DIRINDEX3 = 1;
            }
            else if (PlaneIndex == 0)
            {
                DIRINDEX1 = 1;
                DIRINDEX2 = 2;
                DIRINDEX3 = 0;
            }
            else
            {
                throw std::runtime_error("TrimmedDomainOnPlane :: Contructor :: Wrong PlaneIndex.");
            }
        }

        ///@}
        ///@name Public Operations
        ///@{

        void Reserve(IndexType Size)
        {
            mEdgesPositiveOriented.reserve(Size);
            mEdgesNegativeOriented.reserve(Size);
            mEdgesVertical.reserve(Size);
            mVerticesPositive.reserve(Size);
            mVerticesNegative.reserve(Size);
            mVerticesVertical.reserve(Size);
        }

        void InsertEdge(const std::array<Point3DType, 2> rEdge, const Point3DType &rNormal)
        {
            // Get vertices of edge
            const auto &V1 = rEdge[0];
            const auto &V2 = rEdge[1];

            Point2DType v1{};
            Point2DType v2{};
            // Make sure edge is oriented along DIRINDEX1 such that x2 > x1
            bool vertical_flag = false;
            if (std::abs(V2[DIRINDEX1] - V1[DIRINDEX1]) < EPS3)
                vertical_flag = true;

            if (V2[DIRINDEX1] > V1[DIRINDEX1])
            {
                v1 = {V1[DIRINDEX1], V1[DIRINDEX2]};
                v2 = {V2[DIRINDEX1], V2[DIRINDEX2]};
            }
            else
            {
                v1 = {V2[DIRINDEX1], V2[DIRINDEX2]};
                v2 = {V1[DIRINDEX1], V1[DIRINDEX2]};
            }
            // Positive oriented
            if ((rNormal[DIRINDEX2] > EPS3) && (!vertical_flag))
            {
                // Get unique vertex index.GetVertexID
                auto indices = GetVertexID(v1, v2, Orientation::Positive);
                if (indices.first != indices.second)
                {
                    InsertVertex(v1, indices.first, Orientation::Positive);
                    InsertVertex(v2, indices.second, Orientation::Positive);
                    mEdgesPositiveOriented.push_back(Egde2D(indices.first, indices.second));
                }
            } // Negative oriented
            else if ((rNormal[DIRINDEX2] < -EPS3) && (!vertical_flag))
            {

                // Get unique vertex index.
                auto indices = GetVertexID(v1, v2, Orientation::Negative);
                if (indices.first != indices.second)
                {
                    InsertVertex(v1, indices.first, Orientation::Negative);
                    InsertVertex(v2, indices.second, Orientation::Negative);
                    mEdgesNegativeOriented.push_back(Egde2D(indices.first, indices.second));
                }
            }
            else
            {

                auto index_1 = InsertVertex(v1, Orientation::Vertical);
                auto index_2 = InsertVertex(v2, Orientation::Vertical);
                mEdgesVertical.push_back(Egde2D(index_1, index_2));
            }
        }

        ///@brief Insert new edges to trimmed domain on plane. Check for dublicate nodes. Only inserts new node if unique.
        ///       Only inserts edges if std::abs(normal[DIRINDEX2]) > EPS2. All other edges have no contribution for constant terms
        ///       of moment fitting equation.
        ///@param rEdged
        ///@param rNormal
        void InsertEdges(const std::vector<std::array<Point3DType, 2>> rEdges, const Point3DType &rNormal)
        {
            for (auto &edge : rEdges)
            {
                InsertEdge(edge, rNormal);
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
        TriangleMeshPtrType pGetTriangulation(const TriangleMesh &rTriangleMesh)
        {
            CollectEdgesOnPlane(rTriangleMesh);
            // std::cout << "waf: " << std::endl;
            // std::cout << mEdgesNegativeOriented.size() << std::endl;
            // std::cout << mEdgesPositiveOriented.size() << std::endl;
            // std::cout << mEdgesVertical.size() << std::endl;
            RefineEdges();
            return Triangulate();
        }
        ///@}

    private:
        ///@name Public Operations
        ///@{
        void CollectEdgesOnPlane(const TriangleMesh &rTriangleMesh)
        {
            for (IndexType triangle_id = 0; triangle_id < rTriangleMesh.NumOfTriangles(); ++triangle_id)
            {
                const auto &P1 = rTriangleMesh.P1(triangle_id);
                const auto &P2 = rTriangleMesh.P2(triangle_id);
                const auto &P3 = rTriangleMesh.P3(triangle_id);
                const auto &normal = rTriangleMesh.Normal(triangle_id);

                const double area = rTriangleMesh.Area(triangle_id);
                if (IsPointOnPlane(P1) && IsPointOnPlane(P2))
                {
                    std::array<PointType, 2> new_edge = {P1, P2};
                    InsertEdge(new_edge, normal);
                }
                else if (IsPointOnPlane(P2) && IsPointOnPlane(P3))
                {
                    std::array<PointType, 2> new_edge = {P2, P3};
                    InsertEdge(new_edge, normal);
                }
                else if (IsPointOnPlane(P3) && IsPointOnPlane(P1))
                {
                    std::array<PointType, 2> new_edge = {P3, P1};
                    InsertEdge(new_edge, normal);
                }
            }
        }

        inline bool IsPointOnPlane(const PointType &rPoint, double Tolerance = EPS1) const
        {
            if (mUpperBoundary)
            {
                if ((rPoint[DIRINDEX3] < mUpperBound[DIRINDEX3] + Tolerance) && (rPoint[DIRINDEX3] > mUpperBound[DIRINDEX3] - Tolerance))
                {
                    return true;
                }
                return false;
            }
            else
            {
                if ((rPoint[DIRINDEX3] < mLowerBound[DIRINDEX3] + Tolerance) && (rPoint[DIRINDEX3] > mLowerBound[DIRINDEX3] - Tolerance))
                {
                    return true;
                }
                return false;
            }
        }

        ///@brief Inserts point to container if point is not contained in container.
        ///@param rPoint NewPoint
        ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
        ///@param Return vertex index in container.
        // IndexType GetVertexID(const Point2DType& rPoint, Orientation orientation){

        //     auto& r_vertices = GetVertices(orientation);
        //     // Find new vertex.
        //     auto res = std::find_if( r_vertices.begin(), r_vertices.end(),
        //         [&rPoint](const auto& x) { return  (rPoint[0] > x[0] - EPS2 && rPoint[0] < x[0] + EPS2
        //                                             && rPoint[1] > x[1] - EPS2 && rPoint[1] < x[1] + EPS2) ;});

        //     if( res != r_vertices.end() ){ // Vertex already exists
        //         IndexType index = std::distance(r_vertices.begin(), res);
        //         return index;
        //     } else { // Add new vertex
        //         IndexType index = r_vertices.size();
        //         //r_vertices.push_back(rPoint);
        //         return index;
        //     }
        // }

        ///@brief Inserts point to container if point is not contained in container.
        ///@param rPoint NewPoint
        ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
        ///@param Return vertex index in container.
        std::pair<IndexType, IndexType> GetVertexID(const Point2DType &rV1, const Point2DType &rV2, Orientation orientation)
        {

            auto &r_vertices = GetVertices(orientation);
            auto &r_edges = GetEdges(orientation);

            // If vertices are the same.
            if (rV1[0] > rV2[0] - EPS1 && rV1[0] < rV2[0] + EPS1 && rV1[1] > rV2[1] - EPS1 && rV1[1] < rV2[1] + EPS1)
            {
                return std::make_pair<IndexType, IndexType>(0UL, 0UL);
            }

            // Find new vertex.
            auto res1 = std::find_if(r_vertices.begin(), r_vertices.end(),
                                     [&rV1](const auto &x)
                                     { return (rV1[0] > x[0] - EPS1 && rV1[0] < x[0] + EPS1 && rV1[1] > x[1] - EPS1 && rV1[1] < x[1] + EPS1); });

            IndexType index_1 = 0UL;
            IndexType index_2 = 0UL;
            IndexType v1_is_new = 0UL;
            if (res1 != r_vertices.end())
            { // Vertex already exists
                index_1 = std::distance(r_vertices.begin(), res1);
            }
            else
            { // Add new vertex
                index_1 = r_vertices.size();
                v1_is_new = 1UL;
            }
            auto res_edges = std::find_if(r_edges.begin(), r_edges.end(),
                                          [index_1](const auto &edge)
                                          { return (edge.V1() == index_1); });

            if (res_edges == r_edges.end())
            {
                auto res2 = std::find_if(r_vertices.begin(), r_vertices.end(),
                                         [&rV2](const auto &x)
                                         { return (rV2[0] > x[0] - EPS1 && rV2[0] < x[0] + EPS1 && rV2[1] > x[1] - EPS1 && rV2[1] < x[1] + EPS1); });
                if (res2 != r_vertices.end())
                { // Vertex already exists
                    index_2 = std::distance(r_vertices.begin(), res2);
                }
                else
                { // Add new vertex
                    index_2 = r_vertices.size() + v1_is_new;
                }
                auto res2_edges = std::find_if(r_edges.begin(), r_edges.end(),
                                               [index_2](const auto &edge)
                                               { return (edge.V2() == index_2); });

                if (res2_edges == r_edges.end())
                {
                    return std::make_pair(index_1, index_2);
                }
            }
            return std::make_pair<IndexType, IndexType>(0UL, 0UL);
        }

        ///@brief Fast insertion of point to container. Check if point is already contained is omittede.
        ///@param rPoint NewPoint
        ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
        ///@param Return vertex index in container.
        IndexType InsertVertex(const Point2DType &rPoint, Orientation orientation)
        {

            auto &r_vertices = GetVertices(orientation);
            IndexType index = r_vertices.size();
            r_vertices.push_back(rPoint);
            return index;
        }

        ///@brief Fast insertion of point to container. Check if point is already contained is omittede.
        ///@param rPoint NewPoint
        ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
        ///@param Return vertex index in container.
        void InsertVertex(const Point2DType &rPoint, IndexType NewIndex, Orientation orientation)
        {

            auto &r_vertices = GetVertices(orientation);
            IndexType index = r_vertices.size();
            if (NewIndex < index)
            {
                r_vertices[NewIndex] = rPoint;
            }
            else if (NewIndex == index)
            {
                r_vertices.push_back(rPoint);
            }
            else
            {
                throw std::runtime_error(" Problem!!!!!");
            }
        }

        ///@brief Get Vertex V1 by EdgeId
        ///@param EdgeId
        ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
        const Point2DType &V1byEdgeIdBool(IndexType EdgeId, bool Positive)
        {
            if (Positive)
            {
                return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V1()];
            }
            else
            {
                return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V1()];
            }
        }

        const Point2DType &V1byEdgeId(IndexType EdgeId, OrientationType orientation)
        {
            switch (orientation)
            {
            case Orientation::Positive:
                return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V1()];
            case Orientation::Negative:
                return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V1()];
            case Orientation::Vertical:
                return mVerticesVertical[mEdgesVertical[EdgeId].V1()];
            default:
                throw std::runtime_error("V1byEdgeId");
            }
        }

        ///@brief Get Vertex V2 by EdgeId
        ///@param EdgeId
        ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
        const Point2DType &V2byEdgeIdBool(IndexType EdgeId, bool Positive)
        {
            if (Positive)
            {
                return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V2()];
            }
            else
            {
                return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V2()];
            }
        }

        ///@brief Get Vertex V2 by EdgeId
        ///@param EdgeId
        ///@param Positive Orientation of corresponding edge. Upper edge: Normal[DIRINDEX2]>0 -> Positive.
        const Point2DType &V2byEdgeId(IndexType EdgeId, Orientation orientation)
        {
            switch (orientation)
            {
            case Orientation::Positive:
                return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V2()];
            case Orientation::Negative:
                return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V2()];
            case Orientation::Vertical:
                return mVerticesVertical[mEdgesVertical[EdgeId].V2()];
            default:
                throw std::runtime_error("V2byEdgeId");
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
        void RefineEdges()
        {
            IndexType num_positive_edges = mEdgesPositiveOriented.size();
            IndexType num_negative_edges = mEdgesNegativeOriented.size();

            // Remove duplicates
            // TODO: this can probably be remove again.
            // mEdgesNegativeOriented.erase(std::unique(mEdgesNegativeOriented.begin(), mEdgesNegativeOriented.end(),
            //     [](const auto& edge1, const auto& edge2) { return (edge1.V1() == edge2.V1()) && (edge1.V2() == edge2.V2());  }) , mEdgesNegativeOriented.end());

            // mEdgesPositiveOriented.erase(std::unique(mEdgesPositiveOriented.begin(), mEdgesPositiveOriented.end(),
            //     [](const auto& edge1, const auto& edge2) { return (edge1.V1() == edge2.V1()) && (edge1.V2() == edge2.V2());  }) , mEdgesPositiveOriented.end());
            // mEdgesVertical.erase(std::unique(mEdgesVertical.begin(), mEdgesVertical.end(),
            //     [](const auto& edge1, const auto& edge2) { return (edge1.V1() == edge2.V1()) && (edge1.V2() == edge2.V2());  }) , mEdgesVertical.end());

            // for( auto& edge : mEdgesPositiveOriented ){
            //     std::cout << edge.V1() << ", " << edge.V2() << std::endl;
            // }

            std::vector<double> intersected_vertex_ids{};

            intersected_vertex_ids.push_back(mLowerBound[DIRINDEX1]);
            FindAllIntersectionWithUpperEdge(intersected_vertex_ids, Orientation::Positive);
            FindAllIntersectionWithUpperEdge(intersected_vertex_ids, Orientation::Negative);
            FindAllIntersectionWithUpperEdge(intersected_vertex_ids, Orientation::Vertical);

            // intersected_vertex_ids.push_back( -1.5 );
            intersected_vertex_ids.push_back(mUpperBound[DIRINDEX1]);

            // Sort
            std::sort(intersected_vertex_ids.begin(), intersected_vertex_ids.end());

            // Remove duplicates
            intersected_vertex_ids.erase(std::unique(intersected_vertex_ids.begin(), intersected_vertex_ids.end(),
                                                     [](const auto &v1, const auto &v2)
                                                     { return std::abs(v1 - v2) < EPS2; }),
                                         intersected_vertex_ids.end());

            double plane_position;
            if (mUpperBoundary)
            {
                plane_position = mUpperBound[DIRINDEX3];
            }
            else
            {
                plane_position = mLowerBound[DIRINDEX3];
            }

            for (IndexType i = 0; i < intersected_vertex_ids.size() - 1; ++i)
            {
                double v_left = intersected_vertex_ids[i];
                double v_right = intersected_vertex_ids[i + 1];
                if ((v_right - v_left) > EPS3)
                {
                    double center = v_left + 0.5 * (v_right - v_left);

                    Point3DType new_point{};
                    new_point[DIRINDEX1] = center;
                    new_point[DIRINDEX2] = mUpperBound[DIRINDEX2];
                    new_point[DIRINDEX3] = plane_position;
                    if (mpTrimmedDomain->IsInsideTrimmedDomain(new_point))
                    {
                        auto indices = GetVertexID({v_left, mUpperBound[DIRINDEX2]}, {v_right, mUpperBound[DIRINDEX2]}, Orientation::Positive);
                        if (indices.first != indices.second)
                        {
                            InsertVertex({v_left, mUpperBound[DIRINDEX2]}, indices.first, Orientation::Positive);
                            InsertVertex({v_right, mUpperBound[DIRINDEX2]}, indices.second, Orientation::Positive);
                            mEdgesPositiveOriented.push_back(Egde2D(indices.first, indices.second));
                        }
                    }
                }
            }

            /// Split negative edges at positions of positive points.
            for (int i = 0; i < mVerticesPositive.size(); ++i)
            {
                const auto &v = mVerticesPositive[i];
                SetSplitPoint(v, mEdgesNegativeOriented, Orientation::Positive);
            }

            /// Split positive edges at positions of negative points.
            for (int i = 0; i < mVerticesNegative.size(); ++i)
            {
                const auto &v = mVerticesNegative[i];
                SetSplitPoint(v, mEdgesPositiveOriented, Orientation::Negative);
            }

            SplitEdgesAtSplitPoint(mEdgesNegativeOriented, Orientation::Negative);
            SplitEdgesAtSplitPoint(mEdgesPositiveOriented, Orientation::Positive);
        }

        void FindAllIntersectionWithUpperEdge(std::vector<double> &rVertices, Orientation orientation, double Tolerance = EPS0)
        {
            for (IndexType edge_id = 0; edge_id < GetNumberEdges(orientation); ++edge_id)
            {
                const auto &v1 = V1byEdgeId(edge_id, orientation);
                const auto &v2 = V2byEdgeId(edge_id, orientation);
                bool v1_on_edge = false;
                if (v1[1] < mUpperBound[DIRINDEX2] + Tolerance && v1[1] > mUpperBound[DIRINDEX2] - Tolerance)
                {
                    rVertices.push_back(v1[0]);
                    v1_on_edge = true;
                }
                bool v2_on_edge = false;
                if (v2[1] < mUpperBound[DIRINDEX2] + Tolerance && v2[1] > mUpperBound[DIRINDEX2] - Tolerance)
                {
                    v2_on_edge = true;
                    rVertices.push_back(v2[0]);
                }
            }
        }

        ///@brief Return EdgeId that overlaps rPoint in DIRINDEX1 direction. Return -1 if no edge is found.
        ///       Found if: rPoint[DIRINDEX1] > EgdeV1[DIRINDEX1] and rPoint[DIRINDEX1] < EgdeV2[DIRINDEX1]
        ///@param rPoint
        ///@param Positive Orientation of edges to be searched.
        ///@param int
        int FindIntersectingEdge(const Point2DType &rPoint, Orientation orientation)
        {
            double min_distance = std::numeric_limits<double>::max();
            IndexType found_id = -1;
            for (int i = 0; i < GetNumberEdges(orientation); ++i)
            {
                const auto &v1 = V1byEdgeId(i, orientation);
                const auto &v2 = V2byEdgeId(i, orientation);
                if (rPoint[0] > v1[0] - EPS3 && rPoint[0] < v2[0] + EPS3)
                {
                    double distance = rPoint[1] - 0.5 * (v2[1] + v1[1]);

                    if (distance > 0.0 && distance < min_distance)
                    {
                        found_id = i;
                        min_distance = distance;
                    }
                }
            }
            return found_id;
        }

        ///@brief Create integration points by ray shooting.
        ///@details Shoot ray from each positive oriented edge downards in negative DIRINDEX2 direction. End point is defined by intersecting
        ///         negative oriented edge or mLowerBound[DIRINDEX2].
        ///@param [out] pIPs (new points are pushed_back(). Old points are preserved.)
        ///@param PlanePosition Distance to origin in DIRINDEX3 direction.
        ///@param rNormal Normal vector of plane
        TriangleMeshPtrType Triangulate(double PlanePosition, const Point3DType &rNormal)
        {
            auto p_new_mesh = std::make_unique<TriangleMesh>();
            auto &r_edges_origin = mEdgesPositiveOriented;
            auto &r_edges_dest = mEdgesNegativeOriented;

            auto orientation_origin = Orientation::Positive;
            auto orientation_dest = Orientation::Negative;

            p_new_mesh->Reserve(5 * r_edges_origin.size());
            if (r_edges_origin.size() > 0)
            {
                const auto &test = V2byEdgeId(r_edges_origin.size() - 1, orientation_origin);
            }
            // First we shoot rays from upper edges downwards.
            for (int i = 0; i < r_edges_origin.size(); ++i)
            {

                // TODO: I think visited is not required.
                // TODO: Improve this.
                if (!r_edges_origin[i].IsVisited())
                {

                    r_edges_origin[i].SetVisited(true);

                    const auto &v1_up = V1byEdgeId(i, orientation_origin);
                    const auto &v2_up = V2byEdgeId(i, orientation_origin);

                    // Get center
                    Point2DType c_positive = {0.5 * (v1_up[0] + v2_up[0]), 0.5 * (v1_up[1] + v2_up[1])};
                    int edge_id_dest = FindIntersectingEdge(c_positive, orientation_dest);
                    bool skip = false;
                    if (edge_id_dest > -1)
                    { // Lower edge is found.
                        r_edges_dest[edge_id_dest].SetVisited(true);
                        const auto &v1_low = V1byEdgeId(edge_id_dest, Orientation::Negative);
                        const auto &v2_low = V2byEdgeId(edge_id_dest, Orientation::Negative);

                        if (std::abs(v1_low[1] - v1_up[1]) < EPS2)
                        {
                            Point3DType tmp_point = {0.0, 0.0, 0.0};
                            tmp_point[DIRINDEX3] = PlanePosition;
                            tmp_point[DIRINDEX1] = v1_low[0];
                            tmp_point[DIRINDEX2] = v1_low[1];
                            IndexType v1 = p_new_mesh->AddVertex(tmp_point);
                            tmp_point[DIRINDEX1] = v2_low[0];
                            tmp_point[DIRINDEX2] = v2_low[1];
                            IndexType v2 = p_new_mesh->AddVertex(tmp_point);
                            tmp_point[DIRINDEX1] = v2_up[0];
                            tmp_point[DIRINDEX2] = v2_up[1];
                            IndexType v3 = p_new_mesh->AddVertex(tmp_point);
                            if (mUpperBoundary)
                                p_new_mesh->AddTriangle(Vector3i(v1, v2, v3));
                            else
                                p_new_mesh->AddTriangle(Vector3i(v2, v1, v3));
                            p_new_mesh->AddNormal(rNormal);

                            skip = true;
                        }
                        else if (std::abs(v2_low[1] - v2_up[1]) < EPS2)
                        {
                            Point3DType tmp_point = {0.0, 0.0, 0.0};
                            tmp_point[DIRINDEX3] = PlanePosition;
                            tmp_point[DIRINDEX1] = v1_low[0];
                            tmp_point[DIRINDEX2] = v1_low[1];
                            IndexType v1 = p_new_mesh->AddVertex(tmp_point);
                            tmp_point[DIRINDEX1] = v2_low[0];
                            tmp_point[DIRINDEX2] = v2_low[1];
                            IndexType v2 = p_new_mesh->AddVertex(tmp_point);
                            tmp_point[DIRINDEX1] = v1_up[0];
                            tmp_point[DIRINDEX2] = v1_up[1];
                            IndexType v3 = p_new_mesh->AddVertex(tmp_point);
                            if (mUpperBoundary)
                                p_new_mesh->AddTriangle(Vector3i(v1, v2, v3));
                            else
                                p_new_mesh->AddTriangle(Vector3i(v2, v1, v3));
                            p_new_mesh->AddNormal(rNormal);

                            skip = true;
                        }
                    }

                    if (!skip)
                    {

                        std::array<Point3DType, 4> corner_points{};
                        corner_points[0][DIRINDEX1] = v1_up[0];
                        corner_points[0][DIRINDEX2] = v1_up[1];
                        corner_points[0][DIRINDEX3] = PlanePosition;

                        corner_points[3][DIRINDEX1] = v2_up[0];
                        corner_points[3][DIRINDEX2] = v2_up[1];
                        corner_points[3][DIRINDEX3] = PlanePosition;

                        if (edge_id_dest > -1)
                        {
                            const auto &v1_low = V1byEdgeId(edge_id_dest, Orientation::Negative);
                            const auto &v2_low = V2byEdgeId(edge_id_dest, Orientation::Negative);

                            corner_points[1][DIRINDEX1] = v1_low[0];
                            corner_points[1][DIRINDEX2] = v1_low[1];
                            corner_points[1][DIRINDEX3] = PlanePosition;

                            corner_points[2][DIRINDEX1] = v2_low[0];
                            corner_points[2][DIRINDEX2] = v2_low[1];
                            corner_points[2][DIRINDEX3] = PlanePosition;
                        }
                        else
                        {
                            corner_points[1][DIRINDEX1] = v1_up[0];
                            corner_points[1][DIRINDEX2] = mLowerBound[DIRINDEX2];
                            corner_points[1][DIRINDEX3] = PlanePosition;

                            corner_points[2][DIRINDEX1] = v2_up[0];
                            corner_points[2][DIRINDEX2] = mLowerBound[DIRINDEX2];
                            corner_points[2][DIRINDEX3] = PlanePosition;
                        }

                        Polygon<4> polygon(rNormal);
                        if (mUpperBoundary)
                        {
                            polygon.AddVertex(corner_points[0]);
                            polygon.AddVertex(corner_points[1]);
                            polygon.AddVertex(corner_points[2]);
                            polygon.AddVertex(corner_points[3]);
                        }
                        else
                        {
                            polygon.AddVertex(corner_points[0]);
                            polygon.AddVertex(corner_points[3]);
                            polygon.AddVertex(corner_points[2]);
                            polygon.AddVertex(corner_points[1]);
                        }
                        const auto p_new_triangles = polygon.pGetTriangleMesh();
                        p_new_mesh->Append(*p_new_triangles);
                    }
                } // End Skip
            }     // for( int i = 0; i < r_edges_origin.size(); ++i){

            return std::move(p_new_mesh);
        }

        ///@brief Constructs normal vector and plane position and calls ConstructIPSbyRayShooting().
        ///@param [out] pIPs
        TriangleMeshPtrType Triangulate()
        {
            Point3DType normal = {0.0, 0.0, 0.0};
            double plane_position = 0.0;
            if (mUpperBoundary)
            {
                normal[DIRINDEX3] = 1.0;
                plane_position = mUpperBound[DIRINDEX3];
            }
            else
            {
                normal[DIRINDEX3] = -1.0;
                plane_position = mLowerBound[DIRINDEX3];
            }

            return Triangulate(plane_position, normal);
        }

        ///@brief Splits all edges at their defined split point.
        ///@param rEdges to be split.
        ///@param Positive Orientation of edges,
        ///@todo See TODO in definition.
        void SplitEdgesAtSplitPoint(std::vector<Egde2D> &rEdges, Orientation orientation)
        {
            std::vector<Egde2D> new_edges{};
            for (int i = 0; i < rEdges.size(); ++i)
            {
                const auto &edge = rEdges[i];
                auto split_points = edge.GetSplitPoints();
                const IndexType num_split_points = split_points.size();
                if (num_split_points > 0)
                {

                    std::sort(split_points.begin(), split_points.end(), [](const auto &point_a, const auto &point_b) -> bool
                              { return point_a[0] < point_b[0]; });

                    IndexType index1 = InsertVertex(V1byEdgeId(i, orientation), orientation);
                    IndexType index2 = InsertVertex(split_points[0], orientation);
                    new_edges.push_back(Egde2D(index1, index2));

                    for (int j = 0; j < num_split_points - 1; ++j)
                    {
                        index1 = InsertVertex(split_points[j], orientation);
                        index2 = InsertVertex(split_points[j + 1], orientation);
                        new_edges.push_back(Egde2D(index1, index2));
                    }

                    index1 = InsertVertex(split_points[num_split_points - 1], orientation);
                    index2 = InsertVertex(V2byEdgeId(i, orientation), orientation);
                    new_edges.push_back(Egde2D(index1, index2));
                }
                else
                {
                    new_edges.push_back(edge);
                }
            }

            // TODO: Make ptr and swap!!!
            rEdges.clear();
            rEdges.insert(rEdges.begin(), new_edges.begin(), new_edges.end());
        }

        ///@brief Find intersecting point on Edge with x=Point[DIRINDEX1] and mark as split point.
        ///@param rPoint Potential split point.
        ///@param rEdges Edges to be searched and marked.
        ///@param Positive Orientation of edges,
        void SetSplitPoint(const Point2DType &rPoint, std::vector<Egde2D> &rEdges, Orientation Orientation)
        {
            // Take bool from Edge.
            double current_distance = 1e15;
            int edge_id = -1;
            Point2DType intersection_point{};
            bool is_positive = IsPositive(Orientation);
            auto &wrong_egdes = GetEdges(Orientation);
            for (IndexType i = 0; i < wrong_egdes.size(); ++i)
            {
                const auto &v1 = V1byEdgeIdBool(i, is_positive);
                const auto &v2 = V2byEdgeIdBool(i, is_positive);

                double pos1 = rPoint[0];
                double pos2 = 0.0;
                // Keep looser tolerance here.
                if (pos1 > v1[0] - EPS2 && pos1 < v2[0] + EPS2)
                {
                    pos2 = v1[1] + (v2[1] - v1[1]) / (v2[0] - v1[0]) * (rPoint[0] - v1[0]);
                    double distance = is_positive ? (rPoint[1] - pos2) : (pos2 - rPoint[1]);

                    if ((distance < current_distance) && distance > EPS2)
                    {
                        current_distance = distance;
                    }
                }
            }

            for (IndexType i = 0; i < rEdges.size(); ++i)
            {
                const auto &v1 = V1byEdgeIdBool(i, !is_positive);
                const auto &v2 = V2byEdgeIdBool(i, !is_positive);

                double pos1 = rPoint[0];
                double pos2 = 0.0;

                if (pos1 >= v1[0] + EPS2 && pos1 <= v2[0] - EPS2)
                {
                    pos2 = v1[1] + (v2[1] - v1[1]) / (v2[0] - v1[0]) * (rPoint[0] - v1[0]);
                    double distance = is_positive ? (rPoint[1] - pos2) : (pos2 - rPoint[1]);

                    if ((distance < current_distance) && distance > EPS2)
                    {
                        current_distance = distance;
                        intersection_point[0] = pos1;
                        intersection_point[1] = pos2;
                        edge_id = i;
                    }
                }
            }
            if (edge_id > -1)
            {
                rEdges[edge_id].AddSplitPoint(intersection_point);
            }
        }

        bool IsPositive(Orientation Orientation)
        {
            switch (Orientation)
            {
            case Orientation::Positive:
                return true;
            case Orientation::Negative:
                return false;
            case Orientation::Vertical:
                throw std::runtime_error("error");
            default:
                throw std::runtime_error("error");
            }
        }

        ///@brief Returns Edges container.  (const version)
        ///@param Positive orientation of edges
        ///@return const std::vector<Egde2D>&
        const std::vector<Egde2D> &GetEdges(Orientation Orientation) const
        {
            switch (Orientation)
            {
            case Orientation::Positive:
                return mEdgesPositiveOriented;
            case Orientation::Negative:
                return mEdgesNegativeOriented;
            case Orientation::Vertical:
                return mEdgesVertical;
            default:
                throw std::runtime_error("error");
            }
        }

        ///@brief Returns Edges container.  (non const version)
        ///@param Positive orientation of edges
        ///@return std::vector<Egde2D>&
        std::vector<Egde2D> &GetEdges(Orientation Orientation)
        {
            switch (Orientation)
            {
            case Orientation::Positive:
                return mEdgesPositiveOriented;
            case Orientation::Negative:
                return mEdgesNegativeOriented;
            case Orientation::Vertical:
                return mEdgesVertical;
            default:
                throw std::runtime_error("error");
            }
        }

        ///@brief Returns vertices container.  (const version)
        ///@param Positive orientation of vertices
        ///@return const std::vector<Point2DType>&
        const std::vector<Point2DType> &GetVertices(Orientation Orientation) const
        {
            switch (Orientation)
            {
            case Orientation::Positive:
                return mVerticesPositive;
            case Orientation::Negative:
                return mVerticesNegative;
            case Orientation::Vertical:
                return mVerticesVertical;
            default:
                throw std::runtime_error("Error");
                break;
            }
        }

        ///@brief Returns vertices container.  (non-const version)
        ///@param Positive orientation of vertices
        ///@return std::vector<Point2DType>&
        std::vector<Point2DType> &GetVertices(Orientation orientation)
        {
            switch (orientation)
            {
            case Orientation::Positive:
                return mVerticesPositive;
            case Orientation::Negative:
                return mVerticesNegative;
            case Orientation::Vertical:
                return mVerticesVertical;
            default:
                throw std::runtime_error("Error");
                break;
            }
        }

        ///@brief Returns numver of edges.
        ///@param Positive Orientation of edges.
        ///@return IndexType
        IndexType GetNumberEdges(Orientation orientation)
        {
            switch (orientation)
            {
            case Orientation::Positive:
                return mEdgesPositiveOriented.size();
            case Orientation::Negative:
                return mEdgesNegativeOriented.size();
            case Orientation::Vertical:
                return mEdgesVertical.size();
            default:
                throw std::runtime_error("Error");
                break;
            }
        }

        ///@}
        ///@name Private Members
        ///@{

        std::vector<Egde2D> mEdgesPositiveOriented{};
        std::vector<Egde2D> mEdgesNegativeOriented{};
        std::vector<Egde2D> mEdgesVertical{};
        std::vector<Point2DType> mVerticesPositive{};
        std::vector<Point2DType> mVerticesNegative{};
        std::vector<Point2DType> mVerticesVertical{};

        IndexType DIRINDEX1; // In plane 1
        IndexType DIRINDEX2; // In plane 2
        IndexType DIRINDEX3; // Out of plane

        Point3DType mLowerBound; // Lower bound AABB
        Point3DType mUpperBound; // Upper bound AABB

        bool mUpperBoundary; // Is current plane upper bound?
        const TrimmedDomainBase *mpTrimmedDomain;
        ///@}
    }; // End TrimmedDomainOnPlane

    ///@} // End Tibra Classes

} // End namespace tibra
#endif