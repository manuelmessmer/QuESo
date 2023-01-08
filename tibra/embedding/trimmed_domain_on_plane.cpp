// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes

//// Project includes
#include "embedding/trimmed_domain_on_plane.h"
#include "embedding/polygon.h"

namespace tibra
{

typedef TrimmedDomainOnPlane::Point2DType Point2DType;
typedef TrimmedDomainOnPlane::TriangleMeshPtrType TriangleMeshPtrType;
typedef TrimmedDomainOnPlane::Egde2D Egde2D;
typedef TrimmedDomainOnPlane::Point2DType Point2DType;
typedef TrimmedDomainOnPlane::Point2DSetType Point2DSetType;


void TrimmedDomainOnPlane::CollectEdgesOnPlane(const TriangleMesh &rTriangleMesh)
{
    const auto& edges_on_planes = rTriangleMesh.GetEdgesOnPlanes();
    const IndexType plane_index = DIRINDEX3*2UL + static_cast<IndexType>(mUpperBoundary);

    const auto& r_vertices = rTriangleMesh.GetVertices();
    const auto& edges_on_plane = edges_on_planes[plane_index];
    // Should only require 2UL*edges_on_plane.size(), however small buffer is used.
    // Note, initial capacity must be large enough. New allocation is not allowed and will crash.
    Reserve( std::max(3UL*edges_on_plane.size(), 10UL) );
    for( const auto& edge : edges_on_plane ){
        const IndexType vertex_index_1 = std::get<0>(edge);
        const IndexType vertex_index_2 = std::get<1>(edge);
        const IndexType normal_index = std::get<2>(edge);

        const auto& P1 = r_vertices[vertex_index_1];
        const auto& P2 = r_vertices[vertex_index_2];
        const auto& normal = rTriangleMesh.Normal(normal_index);
        InsertEdge(P1, P2, normal);
    }
}

void TrimmedDomainOnPlane::CloseContourEdges() {
    /// 1.Checks if positive oriented edges span the entire AABB

    /// Find all intersections of Edges with upper bound of AABB.
    ///     UP2 x1-------x2---x3
    ///         | out   /     |      Point  x2 will be found.
    ///         |------/      |      Points x1 and x3 are always added.
    ///     LB2 |  inside     |
    ///        LB1           UP1
    /// LB1 - lower bound of AABB in DIRINDEX1
    /// UP1 - upper bound of AABB in DIRINDEX1
    /// We store only the value in DIRINDEX1 in interected_vertices, since the position in
    /// DIRINDEX2 is always mUpperBound[DIRINDEX2].
    std::vector<double> interected_vertices{};
    interected_vertices.reserve(20);
    interected_vertices.push_back(mLowerBound[DIRINDEX1]);
    FindAllIntersectionWithUpperBound(interected_vertices, Orientation::Positive);
    FindAllIntersectionWithUpperBound(interected_vertices, Orientation::Negative);
    FindAllIntersectionWithUpperBound(interected_vertices, Orientation::Vertical);
    interected_vertices.push_back(mUpperBound[DIRINDEX1]);
    // Sort vertices
    std::sort(interected_vertices.begin(), interected_vertices.end());
    // Remove duplicates
    interected_vertices.erase(std::unique(interected_vertices.begin(), interected_vertices.end(),
        [](const auto &v1, const auto &v2){ return std::abs(v1 - v2) < EPS1; }), interected_vertices.end());

    /// Loop over intersected vertices.
    /// We find the center between x_i and x_i+1 and check if it is inside or outside.
    /// If center is inside, then a new edge bewteen x_i and x_i+1 is drawn.
    ///
    ///     UP2 x1--c1--x2-c2-x3
    ///         | out   /     |      c1 is outside -> no new line.
    ///         |------/      |      c2 is inside -> new edge from x2-x3.
    ///     LB2 |  inside     |      c1 is center between x1 and x2, etc.
    ///        LB1           UP1
    const double plane_position = GetPlanePosition();
    for (IndexType i = 0; i < interected_vertices.size() - 1; ++i) {
        const double v_left = interected_vertices[i];
        const double v_right = interected_vertices[i + 1];
        if ((v_right - v_left) > EPS3)
        {
            double center = v_left + 0.5 * (v_right - v_left);
            Point3DType new_point{}; // We need 3D point, to use IsInsideTrimmedDomain().
            new_point[DIRINDEX1] = center;
            new_point[DIRINDEX2] = mUpperBound[DIRINDEX2];
            new_point[DIRINDEX3] = plane_position;
            if (mpTrimmedDomain->IsInsideTrimmedDomain(new_point)) {
                InsertEdge({v_left, mUpperBound[DIRINDEX2]}, {v_right, mUpperBound[DIRINDEX2]}, Orientation::Positive);
            }
        }
    }

    /// 2. We split all positive edges at all negative oriented vertices and vice versa.
    ///                 x--^-----x---^--x   positive oriented edges
    ///                 |  |     |   |                                     DIRINDEX2
    ///                 |  |     |   |                                         ^
    ///            x----v--x-----v---x      negative oriented edges            |---> DIRINDEX1

    /// Set split points at negative oriented vertices.
    for (int i = 0; i < mVerticesPositive.size(); ++i) {
        SetSplitPoint(mVerticesPositive[i], Orientation::Negative);
    }

    /// Set split points at positive oriented vertices.
    for (int i = 0; i < mVerticesNegative.size(); ++i) {
        SetSplitPoint(mVerticesNegative[i], Orientation::Positive);
    }

    /// Split edges
    SplitEdgesAtSplitPoint(Orientation::Negative);
    SplitEdgesAtSplitPoint(Orientation::Positive);
}

TriangleMeshPtrType TrimmedDomainOnPlane::TriangulateDomain() const
{
    Point3DType normal = {0.0, 0.0, 0.0};
    double plane_position = GetPlanePosition();
    if (mUpperBoundary) {
        normal[DIRINDEX3] = 1.0;
    }
    else {
        normal[DIRINDEX3] = -1.0;
    }

    auto orientation_origin = Orientation::Positive;
    auto orientation_dest = Orientation::Negative;

    auto &r_edges_origin = GetEdges(orientation_origin);
    auto &r_edges_dest = GetEdges(orientation_dest);

    ///Instantiate new mesh ptr.
    auto p_new_mesh = MakeUnique<TriangleMesh>();
    p_new_mesh->Reserve(5 * r_edges_origin.size());

    // Loop over positive oriented edges.
    for (int i = 0; i < r_edges_origin.size(); ++i) {
        const auto &v1_up = V1byEdgeId(i, orientation_origin);
        const auto &v2_up = V2byEdgeId(i, orientation_origin);

        // Get center
        Point2DType c_positive = {0.5 * (v1_up[0] + v2_up[0]), 0.5 * (v1_up[1] + v2_up[1])};
        int edge_id_dest = FindIntersectingEdge(c_positive, orientation_dest);
        bool skip = false;
        if (edge_id_dest > -1) { // Lower edge is found.
            const auto &v1_low = V1byEdgeId(edge_id_dest, Orientation::Negative);
            const auto &v2_low = V2byEdgeId(edge_id_dest, Orientation::Negative);

            /// If v1 of lower and upper edge coincide.
            /// Add:    /|
            ///        / |
            ///         \|
            if (std::abs(v1_low[1] - v1_up[1]) < EPS2) {
                Point3DType tmp_point = {0.0, 0.0, 0.0};
                tmp_point[DIRINDEX3] = plane_position;
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
                p_new_mesh->AddNormal(normal);

                skip = true; // Skip polygon construction.
            }

            else if (std::abs(v2_low[1] - v2_up[1]) < EPS2) {
                /// If v2 of lower and upper edge coincide.
                /// Add:  |\
                ///       | \
                ///       |/
                Point3DType tmp_point = {0.0, 0.0, 0.0};
                tmp_point[DIRINDEX3] = plane_position;
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
                p_new_mesh->AddNormal(normal);

                skip = true; // Skip polygon construction.
            }
        }

        if (!skip) {
            /// Counter clock-wise orientation. Get first and last point of polygon.
            std::array<Point3DType, 4> corner_points{};
            corner_points[0][DIRINDEX1] = v1_up[0];
            corner_points[0][DIRINDEX2] = v1_up[1];
            corner_points[0][DIRINDEX3] = plane_position;

            corner_points[3][DIRINDEX1] = v2_up[0];
            corner_points[3][DIRINDEX2] = v2_up[1];
            corner_points[3][DIRINDEX3] = plane_position;

            /// Second and third point depends if lower edge is found or not.
            if (edge_id_dest > -1) {
                const auto &v1_low = V1byEdgeId(edge_id_dest, Orientation::Negative);
                const auto &v2_low = V2byEdgeId(edge_id_dest, Orientation::Negative);

                corner_points[1][DIRINDEX1] = v1_low[0];
                corner_points[1][DIRINDEX2] = v1_low[1];
                corner_points[1][DIRINDEX3] = plane_position;

                corner_points[2][DIRINDEX1] = v2_low[0];
                corner_points[2][DIRINDEX2] = v2_low[1];
                corner_points[2][DIRINDEX3] = plane_position;
            }
            else {
                corner_points[1][DIRINDEX1] = v1_up[0];
                corner_points[1][DIRINDEX2] = mLowerBound[DIRINDEX2];
                corner_points[1][DIRINDEX3] = plane_position;

                corner_points[2][DIRINDEX1] = v2_up[0];
                corner_points[2][DIRINDEX2] = mLowerBound[DIRINDEX2];
                corner_points[2][DIRINDEX3] = plane_position;
            }

            /// Instantiate new polygon.
            Polygon<4> polygon(normal);
            /// Orientation of triangle depends whether we are on upper or lower bound of AABB.
            if (mUpperBoundary) {
                polygon.AddVertex(corner_points[0]);
                polygon.AddVertex(corner_points[1]);
                polygon.AddVertex(corner_points[2]);
                polygon.AddVertex(corner_points[3]);
            }
            else {
                polygon.AddVertex(corner_points[0]);
                polygon.AddVertex(corner_points[3]);
                polygon.AddVertex(corner_points[2]);
                polygon.AddVertex(corner_points[1]);
            }
            polygon.AddToTriangleMesh(*p_new_mesh.get());
        }
    }
    return std::move(p_new_mesh);
}

void TrimmedDomainOnPlane::FindAllIntersectionWithUpperBound(std::vector<double> &rVertices, OrientationType Orientation, double Tolerance ) {
    for (IndexType edge_id = 0; edge_id < GetNumberEdges(Orientation); ++edge_id) {
        const auto &v1 = V1byEdgeId(edge_id, Orientation);
        const auto &v2 = V2byEdgeId(edge_id, Orientation);
        bool v1_on_edge = false;
        if ( std::abs(v1[1] - mUpperBound[DIRINDEX2]) < Tolerance )
        {
            rVertices.push_back(v1[0]);
            v1_on_edge = true;
        }
        bool v2_on_edge = false;
        if (std::abs(v2[1] - mUpperBound[DIRINDEX2]) < Tolerance )
        {
            v2_on_edge = true;
            rVertices.push_back(v2[0]);
        }
    }
}

int TrimmedDomainOnPlane::FindIntersectingEdge(const Point2DType &rPoint, OrientationType Orientation) const {
    double min_distance = MAXD;
    IndexType found_id = -1;
    for (int i = 0; i < GetNumberEdges(Orientation); ++i) {
        const auto &v1 = V1byEdgeId(i, Orientation);
        const auto &v2 = V2byEdgeId(i, Orientation);
        if (rPoint[0] > v1[0] && rPoint[0] < v2[0] ) {
            const double distance = rPoint[1] - 0.5 * (v2[1] + v1[1]);
            // Distance must be larger than 0.0 and large than alreade found min_distance.
            if (distance > 0.0 && distance < min_distance) {
                found_id = i;
                min_distance = distance;
            }
        }
    }
    return found_id;
}

void TrimmedDomainOnPlane::SetSplitPoint(const Point2DType &rPoint, OrientationType OrientationDest)
{
    // Take bool from Edge.
    double current_distance = 1e15;
    int edge_id = -1;
    Point2DType intersection_point{};
    bool is_positive = OrientationDest == Orientation::Positive;
    OrientationType orientation_origin = is_positive ? Orientation::Negative : Orientation::Positive;

    auto& r_edges_dest = GetEdges(OrientationDest);
    auto& r_edges_origin = GetEdges(orientation_origin);

    for (IndexType i = 0; i < r_edges_origin.size(); ++i) {
        const auto &v1 = V1byEdgeId(i, orientation_origin);
        const auto &v2 = V2byEdgeId(i, orientation_origin);

        double pos1 = rPoint[0];
        double pos2 = 0.0;
        // Keep looser tolerance here.
        if (pos1 > v1[0] - EPS2 && pos1 < v2[0] + EPS2)
        {
            pos2 = v1[1] + (v2[1] - v1[1]) / (v2[0] - v1[0]) * (rPoint[0] - v1[0]);
            double distance = (!is_positive) ? (rPoint[1] - pos2) : (pos2 - rPoint[1]);

            if ((distance < current_distance) && distance > EPS2)
            {
                current_distance = distance;
            }
        }
    }

    for (IndexType i = 0; i < r_edges_dest.size(); ++i)
    {
        const auto &v1 = V1byEdgeId(i, OrientationDest);
        const auto &v2 = V2byEdgeId(i, OrientationDest);

        double pos1 = rPoint[0];
        double pos2 = 0.0;

        if (pos1 >= v1[0] + EPS2 && pos1 <= v2[0] - EPS2)
        {
            pos2 = v1[1] + (v2[1] - v1[1]) / (v2[0] - v1[0]) * (rPoint[0] - v1[0]);
            double distance = (!is_positive) ? (rPoint[1] - pos2) : (pos2 - rPoint[1]);

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
        r_edges_dest[edge_id].AddSplitPoint(intersection_point);
    }
}

void TrimmedDomainOnPlane::SplitEdgesAtSplitPoint(OrientationType Orientation)
{
    auto& r_edges = GetEdges(Orientation);

    IndexType pos = 0;
    IndexType size = r_edges.size();
    while( pos < size ) {
        auto& edge = r_edges[pos];
        auto& split_points = edge.GetSplitPoints();

        const IndexType num_split_points = split_points.size();
        if (num_split_points > 0) {
            std::sort(split_points.begin(), split_points.end(), [](auto &point_a, auto &point_b) -> bool
                { return point_a[0] < point_b[0]; });

            /// Insert egde (vertex V1 + first split point)
            IndexType index1 = edge.V1();
            IndexType index2 = InsertVertex(split_points[0], Orientation);
            r_edges.push_back(Egde2D(index1, index2));

            /// Insert egdes (only split points.)
            for (int j = 0; j < num_split_points - 1; ++j) {
                index1 = InsertVertex(split_points[j], Orientation);
                index2 = InsertVertex(split_points[j + 1], Orientation);
                r_edges.push_back(Egde2D(index1, index2));
            }

            /// Insert egde (last split point + vertex V2)
            index1 = InsertVertex(split_points[num_split_points - 1], Orientation);
            index2 = edge.V2();
            r_edges.push_back(Egde2D(index1, index2));

            split_points.clear();

            // Remove original edge
            r_edges.erase( r_edges.begin() + pos );
            --size;
        }
        else {
            // Keep original edge as is and increment to next one.
            ++pos;
        }
    }
}

////////////////////////
/// Setter Functions ///
////////////////////////

bool TrimmedDomainOnPlane::InsertEdge(const Point2DType& rV1, const Point2DType& rV2, OrientationType Orientation ){
    // Get unique vertex indices.
    auto indices = GetUniqueVertexIDs(rV1, rV2, Orientation);
    // Only insert if indices are not the same.
    if (indices.first != indices.second) {
        InsertVertex(rV1, indices.first, Orientation);
        InsertVertex(rV2, indices.second, Orientation);
        auto& r_edges = GetEdges(Orientation);
        r_edges.push_back(Egde2D(indices.first, indices.second));
        return true;
    }
    return false;
}

void TrimmedDomainOnPlane::InsertEdge(const Point3DType& rV1, const Point3DType& rV2, const Point3DType &rNormal)
{
    // Insantiate 2D points.
    Point2DType v1{};
    Point2DType v2{};

    // Make sure edge is oriented along DIRINDEX1 such that x2 > x1.
    if (rV2[DIRINDEX1] > rV1[DIRINDEX1]) {
        v1 = {rV1[DIRINDEX1], rV1[DIRINDEX2]};
        v2 = {rV2[DIRINDEX1], rV2[DIRINDEX2]};
    }
    else {
        v1 = {rV2[DIRINDEX1], rV2[DIRINDEX2]};
        v2 = {rV1[DIRINDEX1], rV1[DIRINDEX2]};
    }

    if ((rNormal[DIRINDEX2] > EPS3) ) { // Positive oriented
        InsertEdge(v1, v2, Orientation::Positive);
    }
    else if ((rNormal[DIRINDEX2] < -EPS3) ) { // Negative oriented
        InsertEdge(v1, v2, Orientation::Negative);
    }
    else { // Vertical oriented
        InsertEdge(v1, v2, Orientation::Vertical);
    }
}

IndexType TrimmedDomainOnPlane::InsertVertex(const Point2DType &rPoint, OrientationType Orientation) {
    auto &r_vertices = GetVertices(Orientation);
    IndexType index = r_vertices.size();
    r_vertices.push_back(rPoint);
    return index;
}

void TrimmedDomainOnPlane::InsertVertex(const Point2DType &rPoint, IndexType NewIndex, OrientationType Orientation) {
    auto &r_vertices = GetVertices(Orientation);
    auto& r_vertices_set = GetVerticesSet(Orientation);
    IndexType index = r_vertices.size();
    if (NewIndex == index) {
        r_vertices.push_back(rPoint);
        auto res = r_vertices_set.insert(--r_vertices.end());
        if( !res.second )
            throw std::runtime_error("TrimmedDomainOnPlane :: InsertVertex :: Vetrex already exists.");
    }
    else if(NewIndex > index) {
        throw std::runtime_error("TrimmedDomainOnPlane :: InsertVertex :: Given index out of range.");
    }
}

////////////////////////
/// Getter Functions ///
////////////////////////

std::pair<IndexType, IndexType> TrimmedDomainOnPlane::GetUniqueVertexIDs(const Point2DType &rV1, const Point2DType &rV2, OrientationType Orientation) const
{
    const auto& r_vertices = GetVertices(Orientation);
    const auto& r_vertices_set = GetVerticesSet(Orientation);
    const auto& r_edges = GetEdges(Orientation);

    // Instaniate a tmp_vector, as the following functio call require iterators.
    std::vector<Point2DType> tmp_vector = {rV1, rV2};
    // If points are the same.
    if( !PointComparison()(tmp_vector.begin(), ++tmp_vector.begin())  &&
            !PointComparison()(++tmp_vector.begin(), tmp_vector.begin()) ){
        return std::make_pair<IndexType, IndexType>(0UL, 0UL);
    }

    // r_vertices_set is only used to have a fast search here.
    const auto v1_res = r_vertices_set.find(tmp_vector.begin());
    const auto v2_res = r_vertices_set.find(++tmp_vector.begin());

    // If same points are found and both are not r_vertices_set.end().
    if( v1_res == v2_res && v1_res != r_vertices_set.end() ){
        return std::make_pair<IndexType, IndexType>(0UL, 0UL);
    }

    IndexType index_1 = 0UL;
    IndexType index_2 = 0UL;
    IndexType v1_is_new = 0UL;
    if (v1_res != r_vertices_set.end()) { // Vertex 1 already exists
        index_1 = std::distance<std::vector<Point2DType>::const_iterator>(r_vertices.cbegin(), (*v1_res));
    }
    else { // Add new vertex 1
        index_1 = r_vertices.size();
        v1_is_new = 1UL;
    }

    if (v2_res != r_vertices_set.end()) { // Vertex 2 already exists
        index_2 = std::distance<std::vector<Point2DType>::const_iterator>(r_vertices.cbegin(), (*v2_res));
    }
    else { // Add new vertex 2
        index_2 = r_vertices.size() + v1_is_new;
    }

    return std::make_pair(index_1, index_2);
}

const Point2DType& TrimmedDomainOnPlane::V1byEdgeId(IndexType EdgeId, OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V1()];
    case Orientation::Negative:
        return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V1()];
    case Orientation::Vertical:
        return mVerticesVertical[mEdgesVertical[EdgeId].V1()];
    default:
        throw std::runtime_error("TrimmedDomainOnPlane :: V1byEdgeId :: Given Orientation not available.");
    }
}

const Point2DType& TrimmedDomainOnPlane::V2byEdgeId(IndexType EdgeId, OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mVerticesPositive[mEdgesPositiveOriented[EdgeId].V2()];
    case Orientation::Negative:
        return mVerticesNegative[mEdgesNegativeOriented[EdgeId].V2()];
    case Orientation::Vertical:
        return mVerticesVertical[mEdgesVertical[EdgeId].V2()];
    default:
        throw std::runtime_error("TrimmedDomainOnPlane :: V2byEdgeId :: Given Orientation not available.");
    }
}

const std::vector<Egde2D>& TrimmedDomainOnPlane::GetEdges(OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mEdgesPositiveOriented;
    case Orientation::Negative:
        return mEdgesNegativeOriented;
    case Orientation::Vertical:
        return mEdgesVertical;
    default:
        throw std::runtime_error("TrimmedDomainOnPlane :: GetEdges :: Orientation not valid.");
    }
}

std::vector<Egde2D>& TrimmedDomainOnPlane::GetEdges(OrientationType Orientation) {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mEdgesPositiveOriented;
    case Orientation::Negative:
        return mEdgesNegativeOriented;
    case Orientation::Vertical:
        return mEdgesVertical;
    default:
        throw std::runtime_error("TrimmedDomainOnPlane :: GetEdges :: Orientation not valid.");
    }
}

const std::vector<Point2DType>& TrimmedDomainOnPlane::GetVertices(OrientationType Orientation) const
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
        throw std::runtime_error("TrimmedDomainOnPlane :: GetVertices :: Orientation not valid.");
        break;
    }
}

std::vector<Point2DType>& TrimmedDomainOnPlane::GetVertices(OrientationType Orientation)
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
        throw std::runtime_error("TrimmedDomainOnPlane :: GetVertices :: Orientation not valid.");
        break;
    }
}

Point2DSetType& TrimmedDomainOnPlane::GetVerticesSet(OrientationType Orientation) {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mVerticesSetPositive;
    case Orientation::Negative:
        return mVerticesSetNegative;
    case Orientation::Vertical:
        return mVerticesSetVertical;
    default:
        throw std::runtime_error("TrimmedDomainOnPlane :: GetVerticesSet :: Orientation not valid.");
        break;
    }
}

const Point2DSetType& TrimmedDomainOnPlane::GetVerticesSet(OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mVerticesSetPositive;
    case Orientation::Negative:
        return mVerticesSetNegative;
    case Orientation::Vertical:
        return mVerticesSetVertical;
    default:
        throw std::runtime_error("TrimmedDomainOnPlane :: GetVerticesSet :: Orientation not valid.");
        break;
    }
}

IndexType TrimmedDomainOnPlane::GetNumberEdges(OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return mEdgesPositiveOriented.size();
    case Orientation::Negative:
        return mEdgesNegativeOriented.size();
    case Orientation::Vertical:
        return mEdgesVertical.size();
    default:
        throw std::runtime_error("TrimmedDomainOnPlane :: GetNumberEdges :: Orientation not valid.");
        break;
    }
}

double TrimmedDomainOnPlane::GetPlanePosition() const {
    return mUpperBoundary ? mUpperBound[DIRINDEX3] : mLowerBound[DIRINDEX3];
}



} // End namespace tibra