// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIMMED_DOMAIN_ON_PLANE_INCLUDE_H
#define TRIMMED_DOMAIN_ON_PLANE_INCLUDE_H

//// STL includes
#include <set>
//// Project includes
#include "define.hpp"
#include "embedding/brep_operator_base.h"
#include "embedding/trimmed_domain_base.h"

namespace tibra
{
///@name TIBRA Classes
///@{

/**
 * @class  TrimmedDomainOnPlane
 * @author Manuel Messmer
 * @brief Projects the trimmed domain onto a plane of the AABB and provides functions to triangluate
 *        the trimmed domain.
 * @details 1. Stores vertices and edges in different containers according to the
 *             orientation of the corresponding edge (see: CollectEdgesOnPlane(rTriangleMesh)).
 *             Edges are retrieved from EdgesOnPlane storef in rTriangleMesh).
 *          2. New edges are introduced, if neccessary, to close polyline (see. CloseContourEdges()).
 *          3. Domain on plane is triangulated (see: TriangulateDomain()).
 *          All three functions are called in pGetTriangulation().
 */
class TrimmedDomainOnPlane
{
public:
    ///@name Enum's
    ///@{

    enum Orientation {Positive, Negative, Vertical};

    ///@}
    ///@name Type Definitions
    ///@{

    typedef enum Orientation OrientationType;
    typedef std::array<double, 2> Point2DType;
    typedef Vector3d Point3DType;
    typedef std::vector<BoundaryIntegrationPoint> BoundaryIPVectorType;
    typedef Unique<BoundaryIPVectorType> BoundaryIPVectorPtrType;
    typedef Unique<TriangleMesh> TriangleMeshPtrType;


    struct PointComparison {
        PointComparison(double Tolerance) : mTolerance(Tolerance){}
        bool operator() (const std::vector<Point2DType>::iterator& rLhs, const std::vector<Point2DType>::iterator& rRhs) const {
            if( std::abs( (*rLhs)[0] - (*rRhs)[0] ) < mTolerance ){ // If equal
                return (*rLhs)[1] < (*rRhs)[1] - mTolerance;
            }
            else {
                return (*rLhs)[0] <= (*rRhs)[0] - mTolerance;
            }
        }
        double mTolerance;
    };
    typedef std::set<std::vector<Point2DType>::iterator, PointComparison> Point2DSetType;

    /**
     * @class  Edge2D
     * @author Manuel Messmer
     * @brief Simple edge with 2 vertices in 2D.
     */
    class Egde2D {
    public:
        ///@name Life Cycle
        ///@{
        /// Contructor
        Egde2D(IndexType V1, IndexType V2, const Point2DType& rNormal) : mV1(V1), mV2(V2), mNormal(rNormal) {}
        ///@}
        ///@name Public Operations
        ///@{
        IndexType V1() const {return mV1;}
        IndexType V2() const {return mV2;}
        void Set(IndexType V1, IndexType V2){mV1 = V1; mV2 = V2;}
        void SetVerticesOnUpperBoundary(bool V1, bool V2){ mVertexOnUpperBoundary = std::make_pair(V1, V2); }
        std::pair<bool, bool> GetVerticesOnUpperBoundary() const { return mVertexOnUpperBoundary; }
        const std::vector<Point2DType> &GetSplitPoints() const {return mSplitPoints;}
        std::vector<Point2DType>& GetSplitPoints(){return mSplitPoints;}
        void ClearSplitPoints(){mSplitPoints.clear();}
        void AddSplitPoint(const Point2DType &rPoint){mSplitPoints.push_back(rPoint);}
        void Reserve(IndexType NewSize){mSplitPoints.reserve(NewSize);}
        Point2DType Normal() const {return mNormal; }
        ///@}
    private:
        ///@name Private Members
        ///@{
        IndexType mV1, mV2;
        Point2DType mNormal;
        std::vector<Point2DType> mSplitPoints{};
        std::pair<bool, bool> mVertexOnUpperBoundary{};
        ///@}
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    TrimmedDomainOnPlane(IndexType PlaneIndex, bool UpperBoundary, const Point3DType &LowerBound, const Point3DType &UpperBound, const TrimmedDomainBase *pTrimmedDomain, bool Switch)
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
        mPlaneIndex = PlaneIndex;
        if (PlaneIndex == 2) {
            if( Switch ){
                DIRINDEX1 = 0;
                DIRINDEX2 = 1;
            } else {
                DIRINDEX1 = 1;
                DIRINDEX2 = 0;
            }
            DIRINDEX3 = 2;
        }
        else if (PlaneIndex == 1) {
            if( Switch ){
                DIRINDEX1 = 2;
                DIRINDEX2 = 0;
            } else {
                DIRINDEX1 = 2;
                DIRINDEX2 = 0;
            }
            DIRINDEX3 = 1;
        }
        else if (PlaneIndex == 0) {
            if( Switch ){
                DIRINDEX1 = 1;
                DIRINDEX2 = 2;
            } else {
                DIRINDEX1 = 2;
                DIRINDEX2 = 1;
            }
            DIRINDEX3 = 0;
        }
        else {
            TIBRA_ERROR("TrimmedDomainOnPlane::Constructor") << "Wrong PlaneIndex.\n";
        }
        const auto delta = (mUpperBound - mLowerBound);
        const double max_extension = std::max(delta[0], std::max(delta[1], delta[2]));

        mSnapTolerance = std::max( std::max(delta[0], std::max(delta[1], delta[2]))*SNAPTOL, SNAPTOL);

        mVerticesSetPositive = MakeUnique<Point2DSetType>(PointComparison(mSnapTolerance));
        mVerticesSetNegative = MakeUnique<Point2DSetType>(PointComparison(mSnapTolerance));
        mVerticesSetVertical = MakeUnique<Point2DSetType>(PointComparison(mSnapTolerance));

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
    TriangleMeshPtrType pGetTriangulation(const TriangleMesh &rTriangleMesh, const BRepOperatorBase* pOperator) {
        CollectEdgesOnPlane(rTriangleMesh);
        CloseContourEdges(pOperator);
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
    ///        Positive oriented, Negative oriented, and Vertically oriented.
    /// @details
    ///             ^                                 DIRINDEX2
    ///        -----|----- Positive oriented.             ^
    ///                                                   |---> DIRINDEX2
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
    void CloseContourEdges(const BRepOperatorBase* pOperator);

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
    TriangleMeshPtrType TriangulateDomain() const;

    /// @brief Returns intersections (only values in DIRINDEX1 direction) with upper bound (mUpperBound[DIRINDEX2]).
    /// @param [out] rVertices
    /// @param Orientation
    /// @param Tolerance: Default EPS0.
    void FindAllIntersectionWithUpperBound(std::vector<Egde2D> &rEdges, std::vector<Point2DType> &rPoints, OrientationType Orientation);

    ///@brief Return EdgeId that overlaps rPoint in DIRINDEX1 direction. If multiple overlap, closest intersection is returned.
    ///       Return -1 if no edge is found.
    ///       Found if: rPoint[DIRINDEX1] > EgdeV1[DIRINDEX1] and rPoint[DIRINDEX1] < EgdeV2[DIRINDEX1]
    ///@param rPoint
    ///@param Positive Orientation of edges to be searched.
    ///@param int
    int FindIntersectingEdge(const Point2DType &rV1, const Point2DType &rV2, const Point2DType &rNormal, OrientationType Orientation) const;

    ///@brief Return EdgeId that overlaps rPoint in DIRINDEX1 direction. If multiple overlap, closest intersection is returned.
    ///       Return -1 if no edge is found.
    ///       Found if: rPoint[DIRINDEX1] > EgdeV1[DIRINDEX1] and rPoint[DIRINDEX1] < EgdeV2[DIRINDEX1]
    ///@param rPoint
    ///@param Positive Orientation of edges to be searched.
    ///@param int
    bool DoesIntersectEdge(const Point2DType &rV1, const Point2DType &rV2, OrientationType Orientation) const;

    ///@brief Find intersecting point on edge with Orientation==OrientationDest with x=Point[DIRINDEX1] and mark as split point.
    ///@param rPoint Potential split point.
    ///@param OrientationDest Oreintation of edges to split.
    void SetSplitPoint(const Point2DType &rPoint, OrientationType OrientationDest);

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
    bool InsertEdge(const Point2DType& rV1, const Point2DType& rV2, const Point2DType& rNormal, OrientationType rOrientation );

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
    /// Getter Functions ///
    ////////////////////////

    ///@brief Returns vertex IDs of point container, for two points.
    ///       Returns std::pair(0,0), if rV1==rV2 (according to PointComparison()).
    ///       If rV1 or rV2 are new points, a new index is created.
    ///@details std::vector<Point2DType> are the actual vertex containers. This function returns, the index of those vectors.
    ///         std::set<std::vector<Point2DType>::iterator> is used to allow a fast search of duplicate vertices.
    ///@param rV1 Vertex1.
    ///@param rV2 Vertex2.
    ///@param Orientation Current orientation.
    std::pair<IndexType, IndexType> GetUniqueVertexIDs(const Point2DType &rV1, const Point2DType &rV2, OrientationType Orientation) const;

    /// @brief Returns vertex V1 of edge.
    /// @param EdgeId
    /// @param Orientation
    /// @return const Point2DType&
    const Point2DType& V1byEdgeId(IndexType EdgeId, OrientationType Orientation) const;

    /// @brief Returns vertex V2 of edge.
    /// @param EdgeId
    /// @param Orientation
    /// @return const Point2DType&
    const Point2DType& V2byEdgeId(IndexType EdgeId, OrientationType Orientation) const;

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

    ///@brief Returns vertices set.  (const version)
    ///@param Orientation Orientation of vertices
    ///@return Point2DSetType&
    const Point2DSetType& GetVerticesSet(OrientationType Orientation) const;

    ///@brief Returns number of edges.
    ///@param Positive Orientation of edges.
    ///@return IndexType
    IndexType GetNumberEdges(OrientationType Orientation) const;

    /// @brief Returns offset of plane to origin.
    /// @return double
    double GetPlanePosition() const;

    void AddPositiveTouch(Egde2D* pEdge, std::vector<std::pair<double, bool>>& rVertices);

    void AddPositive(Egde2D* pEdge, std::vector<std::pair<double, bool>>& rVertices);
    void AddNegative(Egde2D* pEdge, std::vector<std::pair<double, bool>>& rVertices);
    void AddVertical(Egde2D* pEdge, std::vector<std::pair<double, bool>>& rVertices);
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
    Unique<Point2DSetType> mVerticesSetPositive;
    Unique<Point2DSetType> mVerticesSetNegative;
    Unique<Point2DSetType> mVerticesSetVertical;

    /// Plane direction indices
    IndexType mPlaneIndex;
    IndexType DIRINDEX1; // In plane 1
    IndexType DIRINDEX2; // In plane 2
    IndexType DIRINDEX3; // Out of plane

    Point3DType mLowerBound; // Lower bound AABB
    Point3DType mUpperBound; // Upper bound AABB

    bool mUpperBoundary; // Is current plane upper bound?
    const TrimmedDomainBase *mpTrimmedDomain;

    double mSnapTolerance;
    double mTolerance2;
    ///@}
}; // End TrimmedDomainOnPlane

///@} // End Tibra Classes

} // End namespace tibra
#endif // TRIMMED_DOMAIN_ON_PLANE_INCLUDE_H