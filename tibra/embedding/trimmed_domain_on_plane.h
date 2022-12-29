// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIMMED_DOMAIN_ON_PLANE_INCLUDE_H
#define TRIMMED_DOMAIN_ON_PLANE_INCLUDE_H

//// STL includes
#include <cstddef>
#include <array>
#include <algorithm>
#include <iomanip>
#include <set>
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
     * @brief Projects the trimmed domain onto a plane of the AABB and provides functions to construct
     *        integration points on trimmed domain.
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

        struct PointComparison {
            bool operator() (const std::vector<Point2DType>::iterator& rLhs, const std::vector<Point2DType>::iterator& rRhs) const {
                if( std::abs( (*rLhs)[0] - (*rRhs)[0] ) < EPS1 ){ // If equal
                    return (*rLhs)[1] < (*rRhs)[1] - EPS1;
                }
                else {
                    return (*rLhs)[0] < (*rRhs)[0] - EPS1;
                }
            }
        };

        typedef std::set<std::vector<Point2DType>::iterator, PointComparison> Point2DSetType;
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

            Egde2D(IndexType V1, IndexType V2) : mV1(V1), mV2(V2)
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

            void Set(IndexType V1, IndexType V2){
                mV1 = V1;
                mV2 = V2;
            }

            const std::vector<Point2DType> &GetSplitPoints() const
            {
                return mSplitPoints;
            }

            std::vector<Point2DType>& GetSplitPoints()
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

            ///@}

        private:
            ///@name Private Members
            ///@{
            IndexType mV1;
            IndexType mV2;
            std::vector<Point2DType> mSplitPoints{};
            ///@}
        };

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor.
        TrimmedDomainOnPlane(IndexType PlaneIndex, bool UpperBoundary, const Point3DType &LowerBound, const Point3DType &UpperBound, const TrimmedDomainBase *pTrimmedDomain)
            : mUpperBoundary(UpperBoundary), mLowerBound(LowerBound), mUpperBound(UpperBound), mpTrimmedDomain(pTrimmedDomain)
        {
            // Orientation
            //
            //      _______                  DIRINDEX2
            //     /      /|                ´|`
            //    /_____ / |<-- lower plane  |-->DIRINDEX1
            //    |     |  /                /
            //    |     |</-- upper plane  DIRINDEX3
            //    |_____|/
            //
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

        ///@brief Returns a triangulated mesh of trimmed domain.
        ///@param rTriangleMesh Input mesh. Must hold the edges on planes.
        ///@return TriangleMeshPtrType.
        ///@details see CollectEdgesOnPlane(), CloseContourEdges(), Triangulate().
        //
        //     a_______b                 y
        //     /      /|                ´|`
        //  c /_____d/ |<-- plane lower  |-->x
        //    |     |  /                /
        //    |     |</-- plane upper  Z
        //    |_____|/
        //
        TriangleMeshPtrType pGetTriangulation(const TriangleMesh &rTriangleMesh) {
            CollectEdgesOnPlane(rTriangleMesh);
            CloseContourEdges();
            return TriangulateDomain();
        }

        /// @brief Reserve all containers.
        /// @param NewSize
        void Reserve(IndexType NewSize) {
            mEdgesPositiveOriented.reserve(NewSize);
            mEdgesNegativeOriented.reserve(NewSize);
            mEdgesVertical.reserve(NewSize);
            mVerticesPositive.reserve(NewSize);
            mVerticesNegative.reserve(NewSize);
            mVerticesVertical.reserve(NewSize);
        }

    private:

        ///@}
        ///@name Private Operations
        ///@{

        /// @brief Coolects all edges and stores them in local containers. Sorts edges corresponding to their orientation.
        ///        Positive oriented, Negative oriented, and vertically oriented.
        /// @details
        ///             ^
        ///        -----|----- Positive oriented.
        ///
        ///        -----|----- Negative oriented.
        ///             V
        /// @param rTriangleMesh Clipped triangle mesh inside trimmed domain.
        void CollectEdgesOnPlane(const TriangleMesh &rTriangleMesh);

        ///@brief Refine edges on trimmed domain on plane. See @details for more information.
        ///@details 1.) Checks if positive oriented edges span the entire AABB (we consider only the part that is inside the domain) in DIRINDEX1.
        ///             If not new edges are introduced at top of AABB (UpperBound[DIRINDEX2]) to fill gaps.
        ///        UB2                  UP2       _____
        ///         |  out  /     |      |  out  /     |  new edges introduced since Point b is inside domain.
        ///         |------/      | ---> |------/      |
        ///         |  inside     |      |  inside     |
        ///        LB12          UB1    LB12           UP1
        ///
        ///         UB2
        ///         |  in   /     |      |  in   /     |  no new edges introduced since Point b is outside.
        ///         |------/      | ---> |------/      |       LB1 - lower bound of AABB in DIRINDEX1
        ///         |  outsid     |      |  outside    |       UP1 - upper bound of AABB in DIRINDEX1
        ///        LB12          UB1    LB12           UP1
        ///
        ///         2.) Maps all negative oriented points onto positive oriented edges if overlap and vice versa.
        ///             Such that the vertices of postive and negative oriented edges are at same position in: DIRINDEX1.
        ///                 x--^-----x---^--x   positive oriented edges
        ///                 |  |     |   |                                     DIRINDEX2
        ///                 |  |     |   |                                         ^
        ///            x----v--x-----v---x      negative oriented edges            |---> DIRINDEX1
        void CloseContourEdges();

        ///@brief Triangulates the trimmed domain on plane (see @details for more information).
        ///@details Assumes that the vertices of positive and negative oriented edges coincide along DIRINDEX1. See RefineEdges().
        ///
        ///          x___x   Positive Oriented      DIRINDEX2
        ///         /|   |                              ^
        ///   x___x/ |   |                              |----> DIRINDEX1
        ///   |   |  |   |
        ///   |   |  |   |
        ///   x___x__x___x    Negative Oriented      x-Vertices
        ///
        /// Each vertical stripe is individually triangulated. First we generate a polygon with 4 vertices and the triangulate it.
        /// If vertices of positve and negative edges coincide, only triangle are created.
        ///
        ///          x___x   Positive Oriented
        ///         /|   |
        ///       x/_|___|   Negative Oriented       x-Vertices
        ///       Tri Poly
        ///
        ///@return TriangleMeshPtrType.
        TriangleMeshPtrType TriangulateDomain();

        /// @brief Returns intersections (only values in DIRINDEX1 direction) with upper bound (mUpperBound[DIRINDEX2]).
        /// @param [out] rVertices
        /// @param Orientation
        /// @param Tolerance: Default EPS0.
        void FindAllIntersectionWithUpperBound(std::vector<double> &rVertices, OrientationType Orientation, double Tolerance = EPS0);

        ///@brief Return EdgeId that overlaps rPoint in DIRINDEX1 direction. If multiple overlap, closest intersection is returned.
        ///       Return -1 if no edge is found.
        ///       Found if: rPoint[DIRINDEX1] > EgdeV1[DIRINDEX1] and rPoint[DIRINDEX1] < EgdeV2[DIRINDEX1]
        ///@param rPoint
        ///@param Positive Orientation of edges to be searched.
        ///@param int
        int FindIntersectingEdge(const Point2DType &rPoint, OrientationType Orientation);

        ///@brief Find intersecting point on edge with Orientation==OrientationDest with x=Point[DIRINDEX1] and mark as split point.
        ///@param rPoint Potential split point.
        ///@param OrientationDest Oreintation of edges to split.
        void SetSplitPoint(const Point2DType &rPoint, Orientation OrientationDest);

        ///@brief Splits all edges at their defined split point.
        ///@param Orientation Orientation of edges.
        void SplitEdgesAtSplitPoint(OrientationType Orientation);

        ////////////////////////
        /// Setter Functions ///
        ////////////////////////

        /// @brief Inserts edge into local containers.
        /// @param rV1 Point1 (3D-Point).
        /// @param rV2 Point2 (3D-Point).
        /// @param rNormal Corresponsing normal direction.
        void InsertEdge(const Point3DType& rV1, const Point3DType& rV2, const Point3DType &rNormal);

        /// @brief Inserts edge into local containers.
        /// @param rV1 Point1 (2D-Point).
        /// @param rV2 Point1 (2D-Point).
        /// @param rOrientation Current orientation.
        /// @return True if edge is inserted.
        bool InsertEdge(const Point2DType& rV1, const Point2DType& rV2, OrientationType rOrientation );

        ///@brief Returns vertex IDs of point container, for two points.
        ///       Returns std::pair(0,0), if rV1==rV2 (according to PointComparison()).
        ///       If rV1 or rV2 are new points, a new index is created.
        ///@details std::vector<Point2DType> are the actual vertex containers. This function returns, the index of those vectors.
        ///         std::set<std::vector<Point2DType>::iterator> is used to allow a fast search of duplicate vertices.
        ///@param rV1 Vertex1.
        ///@param rV2 Vertex2.
        ///@param Orientation Current orientation.
        std::pair<IndexType, IndexType> GetUniqueVertexIDs(const Point2DType &rV1, const Point2DType &rV2, OrientationType Orientation);

        ///@brief Simple and fast insertion of vertex to container. Check if point is already contained is omitted.
        ///@brief If uniqueness check of point is required, use: GetUniqueVertexIDs() + InsertVertex(Point2DType, IndexType, Orientation).
        ///@param rPoint NewPoint
        ///@param OrientationType Current orientation.
        ///@param Return Vertex index in container.
        IndexType InsertVertex(const Point2DType &rPoint, OrientationType Orientation);

        ///@brief Inserts point to vertex container by index. If index is known through GetUniqueVertexIDs, this function is used
        ///       to insert point at corresponding position to std::vector<Point2DType> containers.
        ///@param rPoint NewPoint
        ///@param NewIndex Index of new point.
        ///@param OrientationType Current orientation.
        void InsertVertex(const Point2DType &rPoint, IndexType NewIndex, OrientationType Orientation);

        ////////////////////////
        /// Setter Functions ///
        ////////////////////////

        /// @brief Returns vertex V1 of edge.
        /// @param EdgeId
        /// @param Orientation
        /// @return const Point2DType&
        const Point2DType& V1byEdgeId(IndexType EdgeId, OrientationType Orientation);

        /// @brief Returns vertex V2 of edge.
        /// @param EdgeId
        /// @param Orientation
        /// @return const Point2DType&
        const Point2DType& V2byEdgeId(IndexType EdgeId, OrientationType Orientation);

        ///@brief Returns Edges container. (const version)
        ///@param Orientation orientation of edges
        ///@return const std::vector<Egde2D>&
        const std::vector<Egde2D>& GetEdges(OrientationType Orientation) const;

        ///@brief Returns Edges container.  (non const version)
        ///@param Orientation orientation of edges
        ///@return std::vector<Egde2D>&
        std::vector<Egde2D>& GetEdges(OrientationType Orientation);

        ///@brief Returns vertices container. (const version)
        ///@param Orientation Orientation of vertices
        ///@return const std::vector<Point2DType>&
        const std::vector<Point2DType>& GetVertices(OrientationType Orientation) const;

        ///@brief Returns vertices container.  (non-const version)
        ///@param Orientation Orientation of vertices
        ///@return std::vector<Point2DType>&
        std::vector<Point2DType>& GetVertices(OrientationType Orientation);

        ///@brief Returns vertices set.  (non-const version)
        ///@param Orientation Orientation of vertices
        ///@return Point2DSetType&
        Point2DSetType& GetVerticesSet(OrientationType Orientation);

        ///@brief Returns number of edges.
        ///@param Positive Orientation of edges.
        ///@return IndexType
        IndexType GetNumberEdges(OrientationType Orientation);

        /// @brief Returns offset of plane to origin.
        /// @return double
        double GetPlanePosition() const;

        ///@}
        ///@name Private Members
        ///@{

        /// Edges container
        std::vector<Egde2D> mEdgesPositiveOriented{};
        std::vector<Egde2D> mEdgesNegativeOriented{};
        std::vector<Egde2D> mEdgesVertical{};

        /// Vertices container
        std::vector<Point2DType> mVerticesPositive{};
        std::vector<Point2DType> mVerticesNegative{};
        std::vector<Point2DType> mVerticesVertical{};

        /// Vertices sets. Used for fast search, if new vertex already exists.
        Point2DSetType mVerticesSetPositive{};
        Point2DSetType mVerticesSetNegative{};
        Point2DSetType mVerticesSetVertical{};

        /// Plane direction indices
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
#endif // TRIMMED_DOMAIN_ON_PLANE_INCLUDE_H