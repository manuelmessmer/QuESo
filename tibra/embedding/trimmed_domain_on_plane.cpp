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


void TrimmedDomainOnPlane::AddPositiveTouch(Egde2D* pEdge, std::vector<std::pair<double, bool>>& rVertices){
    auto status = pEdge->GetVerticesOnUpperBoundary();
    if( status.first ){
        auto v1 = mVerticesPositive[pEdge->V1()];
        rVertices.push_back( std::make_pair(v1[0], true)  );
        rVertices.push_back( std::make_pair(v1[0], false)  );

    }
    else if( status.second ){
        auto normal = pEdge->Normal();
        auto v2 = mVerticesPositive[pEdge->V2()];
        rVertices.push_back( std::make_pair(v2[0], true)  );
        rVertices.push_back( std::make_pair(v2[0], false)  );
    }
}

void TrimmedDomainOnPlane::AddPositive(Egde2D* pEdge, std::vector<std::pair<double, bool>>& rVertices){
    auto status = pEdge->GetVerticesOnUpperBoundary();
    if( status.first ){
        auto normal = pEdge->Normal();
        auto v1 = mVerticesPositive[pEdge->V1()];
        if (normal[0] > -ZEROTOL) {
            rVertices.push_back( std::make_pair(v1[0], false)  );
        }
    }
    else if( status.second ){
        auto normal = pEdge->Normal();
        auto v2 = mVerticesPositive[pEdge->V2()];
        if ( normal[0] < ZEROTOL ){
            rVertices.push_back( std::make_pair(v2[0], true)  );
        }
    }
}

void TrimmedDomainOnPlane::AddNegative(Egde2D* pEdge, std::vector<std::pair<double, bool>>& rVertices){
    auto status = pEdge->GetVerticesOnUpperBoundary();
    if( status.first ){
        auto normal = pEdge->Normal();
        auto v1 = mVerticesNegative[pEdge->V1()];
        if( normal[0] < ZEROTOL ) {
            rVertices.push_back( std::make_pair(v1[0], true)  );
        }
    }
    else if( status.second ){
        auto normal = pEdge->Normal();
        auto v2 = mVerticesNegative[pEdge->V2()];
        if(normal[0] > -ZEROTOL) {
            rVertices.push_back( std::make_pair(v2[0], false)  );
        }
    }
}

void TrimmedDomainOnPlane::AddVertical(Egde2D* pEdge, std::vector<std::pair<double, bool>>& rVertices){
    auto status = pEdge->GetVerticesOnUpperBoundary();
    if( status.first || status.second ){
        auto normal = pEdge->Normal();
        auto v1 = mVerticesVertical[pEdge->V1()];
        if( normal[0] < ZEROTOL ) {
            rVertices.push_back( std::make_pair(v1[0], true)  );
        } else {
            rVertices.push_back( std::make_pair(v1[0], false)  );
        }
    }
}

void TrimmedDomainOnPlane::CloseContourEdges(const BRepOperatorBase* pOperator) {
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

    std::vector<double> new_positive{};

    std::vector<Egde2D> intersect_edges_positive{};
    std::vector<Point2DType> intersecting_points{};
    std::vector<Egde2D> intersect_edges_negative{};
    std::vector<Egde2D> intersect_edges_vertical{};
    FindAllIntersectionWithUpperBound(intersect_edges_positive, intersecting_points, Orientation::Positive);
    FindAllIntersectionWithUpperBound(intersect_edges_negative, intersecting_points, Orientation::Negative);
    FindAllIntersectionWithUpperBound(intersect_edges_vertical, intersecting_points, Orientation::Vertical);

    // for( auto edge : intersect_edges_vertical){
    //     std::cout << mVerticesVertical[edge.V1()][0] << ", " << mVerticesVertical[edge.V2()][0] << std::endl;
    // }
    // std::cout << "waf: " << std::endl;
    // std::cout << "intersect_edges_positive: " << intersect_edges_positive.size() << std::endl;
    // std::cout << "intersect_edges_negative: " << intersect_edges_negative.size() << std::endl;
    // std::cout << "intersect_edges_vertical: " << intersect_edges_vertical.size() << std::endl;

    // std::cout << "nenenennenen" << std::endl;
    // for( auto edge : intersect_edges_positive){
    //     std::cout << edge.Normal()[0] << std::endl;
    //     auto status = edge.GetVerticesOnUpperBoundary();
    //     std::cout << status.first << ", " << status.second << std::endl;
    //     std::cout << std::setprecision(16) << mVerticesPositive[edge.V1()][0] << ", " << mVerticesPositive[edge.V2()][0] << std::endl;
    // }



    intersect_edges_vertical.erase(
            std::unique(intersect_edges_vertical.begin(), intersect_edges_vertical.end(), [this](const auto &rValue1, const auto &rValue2){
                const bool normal_1_is_neg = rValue1.Normal()[0] < 0.0;
                const bool normal_2_is_neg = rValue2.Normal()[0] < 0.0;

                return normal_1_is_neg == normal_2_is_neg && std::abs(this->mVerticesVertical[rValue1.V1()][0] - (this->mVerticesVertical[rValue2.V1()][0])) < mSnapTolerance; }), intersect_edges_vertical.end() );

    std::vector<Egde2D> intersect_edges_vertical_copy( intersect_edges_vertical.begin(), intersect_edges_vertical.end() );
    intersect_edges_vertical.clear();

    for( auto r_edge : intersect_edges_vertical_copy){
        IndexType count = std::count_if(intersect_edges_vertical_copy.begin(), intersect_edges_vertical_copy.end(), [&r_edge, this](const auto& rValue)
            { return std::abs(this->mVerticesVertical[r_edge.V1()][0] - (this->mVerticesVertical[rValue.V1()][0])) < mSnapTolerance; });

        if(count == 1 )
            intersect_edges_vertical.push_back(r_edge);
    }

    intersect_edges_positive.erase(
            std::unique(intersect_edges_positive.begin(), intersect_edges_positive.end(), [this](const auto &rValue1, const auto &rValue2){
                return std::abs(this->mVerticesPositive[rValue1.V1()][0] - (this->mVerticesPositive[rValue2.V1()][0])) < mSnapTolerance
                    && std::abs(this->mVerticesPositive[rValue1.V2()][0] - (this->mVerticesPositive[rValue2.V2()][0])) < mSnapTolerance; }), intersect_edges_positive.end() );

    for( auto edge : intersect_edges_positive ){
        auto v1 = mVerticesPositive[edge.V1()];
        auto found1 = std::find_if(intersect_edges_vertical.begin(), intersect_edges_vertical.end(), [v1,this](const auto& rValue){
            return std::abs(v1[0]- this->mVerticesVertical[rValue.V1()][0]) < mSnapTolerance;
        }  );
        auto status = edge.GetVerticesOnUpperBoundary();
        if( found1 != intersect_edges_vertical.end() ){
            status.first = false;
        }
        auto v2 = mVerticesPositive[edge.V2()];
        auto found2 = std::find_if(intersect_edges_vertical.begin(), intersect_edges_vertical.end(), [v2,this](const auto& rValue){
            return std::abs(v2[0]- this->mVerticesVertical[rValue.V1()][0]) < mSnapTolerance;
        }  );
        if( found2 != intersect_edges_vertical.end() ){
            status.second = false;
        }
        edge.SetVerticesOnUpperBoundary(status.first, status.second);
    }

    intersect_edges_positive.erase(std::remove_if(intersect_edges_positive.begin(), intersect_edges_positive.end(), [](const auto& rValue) {
                auto status = rValue.GetVerticesOnUpperBoundary();
                return !status.first && !status.second; }), intersect_edges_positive.end());

    intersect_edges_negative.erase(
            std::unique(intersect_edges_negative.begin(), intersect_edges_negative.end(), [this](const auto &rValue1, const auto &rValue2){
                return std::abs(this->mVerticesNegative[rValue1.V1()][0] - (this->mVerticesNegative[rValue2.V1()][0])) < mSnapTolerance
                    && std::abs(this->mVerticesNegative[rValue1.V2()][0] - (this->mVerticesNegative[rValue2.V2()][0])) < mSnapTolerance;
                }), intersect_edges_negative.end() );

    intersect_edges_negative.erase(std::remove_if(intersect_edges_negative.begin(), intersect_edges_negative.end(), [](const auto& rValue) {
                auto status = rValue.GetVerticesOnUpperBoundary();
                return status.first && status.second; }), intersect_edges_negative.end());

    // std::cout << "waf: " << std::endl;
    // std::cout << "intersect_edges_positive: " << intersect_edges_positive.size() << std::endl;
    // std::cout << "intersect_edges_negative: " << intersect_edges_negative.size() << std::endl;
    // std::cout << "intersect_edges_vertical: " << intersect_edges_vertical.size() << std::endl;

    // /// BEGIN: PLOT VALUES
    // std::cout << "waf: " << std::endl;
    // std::cout << "intersect_edges_positive: " << intersect_edges_positive.size() << std::endl;
    // std::cout << "intersect_edges_negative: " << intersect_edges_negative.size() << std::endl;
    // std::cout << "intersect_edges_vertical: " << intersect_edges_vertical.size() << std::endl;
    // for( auto edge : intersect_edges_vertical){
    //     std::cout << edge.Normal()[0] << std::endl;
    //     std::cout << mVerticesVertical[edge.V1()][0] << ", " << mVerticesVertical[edge.V2()][0] << std::endl;
    // }
    // // for( auto edge : intersect_edges_negative){
    // //     std::cout << mVerticesNegative[edge.V1()][0] << ", " << mVerticesNegative[edge.V2()][0] << std::endl;
    // // }
    // for( auto edge : intersect_edges_positive){
    //     // std::cout << edge.Normal()[0] << std::endl;
    //     // auto status = edge.GetVerticesOnUpperBoundary();
    //     // std::cout << status.first << ", " << status.second << std::endl;
    //     std::cout << std::setprecision(16) << mVerticesPositive[edge.V1()][0] << ", " << mVerticesPositive[edge.V2()][0] << std::endl;
    // }
    //TIBRA_ERROR( "waf") << "fuckii\n";
    /// END: PLOT VALUES

    for( auto& edge : intersect_edges_positive){
       auto status = edge.GetVerticesOnUpperBoundary();
       if( status.first && status.second ){
            if( std::abs(mVerticesPositive[edge.V1()][0]-mVerticesPositive[edge.V2()][0]) < mSnapTolerance ){
                edge.SetVerticesOnUpperBoundary(true, false);
            }
       }
    }
    for( auto& edge : intersect_edges_negative){
       auto status = edge.GetVerticesOnUpperBoundary();
       if( status.first && status.second ){
            if( std::abs(mVerticesNegative[edge.V1()][0]-mVerticesNegative[edge.V2()][0]) < mSnapTolerance ){
                edge.SetVerticesOnUpperBoundary(true, false);
            }
       }
    }

    // Remove intermediates
    Point2DType corner_left = { mLowerBound[DIRINDEX1], mUpperBound[DIRINDEX2] };
    Point2DType corner_right = { mUpperBound[DIRINDEX1], mUpperBound[DIRINDEX2] };
    std::vector<Point2DType> tmp_vector = {corner_left, corner_right};
    const auto v1_res = mVerticesSetPositive->find(tmp_vector.begin());
    const auto v2_res = mVerticesSetPositive->find(++tmp_vector.begin());


    bool add_corner_left = v1_res == mVerticesSetPositive->end();
    bool add_corner_right = v2_res == mVerticesSetPositive->end();;
    // for( auto& r_vertex : mVerticesPositive ){
    //     if( std::abs(r_vertex[0]-mLowerBound[DIRINDEX1]) < ZEROTOL )
    //         add_corner_left = false;
    //     if( std::abs(r_vertex[0]-mUpperBound[DIRINDEX1]) < ZEROTOL ){
    //         add_corner_right = false;
    //     }
    // }

    // std::cout << "add_corner_left: " << add_corner_left << std::endl;
    // std::cout << "add_corner_right: " << add_corner_right << std::endl;

    std::vector<std::pair<double, bool>> interected_vertices{};
    if( add_corner_left ){
        interected_vertices.push_back( std::make_pair(mLowerBound[DIRINDEX1], true) );
    }

    const IndexType size_positive = intersect_edges_positive.size();
    const IndexType size_negative = intersect_edges_negative.size();
    const IndexType size_vertical = intersect_edges_vertical.size();
    IndexType pos_positive = 0;
    IndexType pos_negative = 0;
    IndexType pos_vertical = 0;
    while( pos_vertical < size_vertical || pos_negative < size_negative || pos_positive < size_positive){

        Egde2D* edge_positive = nullptr;
        Egde2D* edge_negative = nullptr;
        Egde2D* edge_vertical = nullptr;
        if( pos_positive < size_positive )
            edge_positive = &intersect_edges_positive[pos_positive];
        if( pos_negative < size_negative )
            edge_negative = &intersect_edges_negative[pos_negative];
        if( pos_vertical < size_vertical)
            edge_vertical = &intersect_edges_vertical[pos_vertical];

        IndexType size = interected_vertices.size();
        double left_bound = size > 0 ? interected_vertices[size-1].first : mLowerBound[DIRINDEX1];
        double distance_pos = MAXD;
        bool double_vertex_edge = false;
        if( edge_positive ){
            auto status = edge_positive->GetVerticesOnUpperBoundary();
            if( status.first && status.second ){
                double_vertex_edge = true;
            }
            if( status.first ){
                distance_pos = mVerticesPositive[edge_positive->V1()][0] - left_bound;
            }
            else {
                double_vertex_edge = false;
                distance_pos = mVerticesPositive[edge_positive->V2()][0] - left_bound;
            }

        }

        double distance_neg = MAXD;
        if( edge_negative ){
            auto status = edge_negative->GetVerticesOnUpperBoundary();
            if( status.first )
                distance_neg = mVerticesNegative[edge_negative->V1()][0] - left_bound;
            else
                distance_neg = mVerticesNegative[edge_negative->V2()][0] - left_bound;
        }

        double distance_ver = MAXD;
        if( edge_vertical ){
            auto status = edge_vertical->GetVerticesOnUpperBoundary();
            if( status.first )
                distance_ver = mVerticesVertical[edge_vertical->V1()][0] - left_bound;
            else
                distance_ver = mVerticesVertical[edge_vertical->V2()][0] - left_bound;
        }

        if( std::abs(distance_pos) < mSnapTolerance ){
            const IndexType size = interected_vertices.size();
            if( (!(add_corner_left && size == 1)) &&  (size > 0) ){
                interected_vertices[size-1].second = true;
                double val = interected_vertices[size-1].first;
                interected_vertices.push_back( std::make_pair(val, false)  );
            }
            if( double_vertex_edge ){
                auto status = edge_positive->GetVerticesOnUpperBoundary();
                edge_positive->SetVerticesOnUpperBoundary(false, status.second);
            }
            else {
                ++pos_positive;
            }
        }
        else if( std::abs(distance_ver) < mSnapTolerance ){
            IndexType size = interected_vertices.size();
            if( (!(add_corner_left && size == 1)) &&  (size > 0) ){
                interected_vertices[size-1].second = true;
                double val = interected_vertices[size-1].first;
                interected_vertices.push_back( std::make_pair(val, false)  );
            }
            ++pos_vertical;
        }
        else if( std::abs(distance_neg) < mSnapTolerance ){
            IndexType size = interected_vertices.size();
            // if( (!(add_corner_left && size == 1)) &&  (size > 0) ){
            //     interected_vertices[size-1].second = true;
            //     double val = interected_vertices[size-1].first;
            //     interected_vertices.push_back( std::make_pair(val, false)  );
            // }
            ++pos_negative;
        }
        else if( (distance_pos+distance_ver) < 0.1*MAXD && std::abs(distance_pos-distance_ver) <= mSnapTolerance ){ // Vert and pos are the same
            if( double_vertex_edge ){
                auto status = edge_positive->GetVerticesOnUpperBoundary();
                edge_positive->SetVerticesOnUpperBoundary(false, status.second);
            }
            else {
                ++pos_positive;
            }
            ++pos_vertical;
        } else if ( (distance_pos+distance_neg) < 0.1*MAXD &&  std::abs(distance_pos-distance_neg) <= mSnapTolerance ) {
            ++pos_positive;
            ++pos_negative;
        } else if ( (distance_ver+distance_neg) < 0.1*MAXD && std::abs(distance_ver-distance_neg) <= mSnapTolerance ) {
            //AddVertical(edge_vertical, interected_vertices);
            //++pos_vertical;
            ++pos_negative;
        } else {
            /// ADD POSITIVE
            if( distance_pos < distance_neg && distance_pos < distance_ver ){
                AddPositive(edge_positive, interected_vertices);
                if( double_vertex_edge ){
                    auto status = edge_positive->GetVerticesOnUpperBoundary();
                    edge_positive->SetVerticesOnUpperBoundary(false, status.second);
                }
                else {
                    ++pos_positive;
                }
            } // ADD Negative
            else if( distance_neg < distance_pos && distance_neg < distance_ver ){
                AddNegative(edge_negative, interected_vertices);
                ++pos_negative;
            } // ADD VERTICAL
            else if( distance_ver < distance_neg && distance_ver < distance_pos ){
                AddVertical(edge_vertical, interected_vertices);
                ++pos_vertical;
            }
        }

    }

    if( add_corner_right ){
        IndexType size = interected_vertices.size();
        if( size == 0 || std::abs(interected_vertices[size-1].first - mUpperBound[DIRINDEX1]) > mSnapTolerance ){
            interected_vertices.push_back( std::make_pair(mUpperBound[DIRINDEX1], false) );
        } else {
            add_corner_right = false;
        }
    }

    // std::cout << "star: -------------------" << interected_vertices.size() << std::endl;
    // std::cout << "positive size: " << mEdgesPositiveOriented.size() << std::endl;
    // for( IndexType i = 0; i < interected_vertices.size(); ++i){
    //     std::cout << std::setprecision(16) << interected_vertices[i].first << ", " << interected_vertices[i].second << std::endl;
    // }

    // std::cout << "mEdgesPositiveOriented: " << mEdgesPositiveOriented.size() << std::endl;

    // interected_vertices.clear();
    /// Loop over intersected vertices.
    /// We find the center between x_i and x_i+1 and check if it is inside or outside.
    /// If center is inside, then a new edge bewteen x_i and x_i+1 is drawn.
    ///
    ///     UP2 x1--c1--x2-c2-x3
    ///         | out   /     |      c1 is outside -> no new line.
    ///         |------/      |      c2 is inside -> new edge from x2-x3.
    ///     LB2 |  inside     |      c1 is center between x1 and x2, etc.
    ///        LB1           UP1

    double y_max = mLowerBound[DIRINDEX2];
    for( auto& vertex : mVerticesNegative ){
        if( vertex[1] > y_max ){
            y_max = vertex[1];
        }
    }
    for( auto& vertex : mVerticesPositive ){
        if( vertex[1] > y_max ){
            y_max = vertex[1];
        }
    }
    const double plane_position = GetPlanePosition();
    //std::cout << "COrner: << " << add_corner_left << add_corner_right << std::endl;
    if( add_corner_left && add_corner_right && interected_vertices.size() == 2 ){
        //if( mEdgesPositiveOriented.size() == 0 ){
            const double v_left = interected_vertices[0].first;
            const double v_right = interected_vertices[1].first;

            double center = v_left + 0.5 * (v_right - v_left);
            Point3DType new_point{}; // We need 3D point, to use IsInsideTrimmedDomain().
            new_point[DIRINDEX1] = center;
            new_point[DIRINDEX2] = 0.5*(mUpperBound[DIRINDEX2]+y_max);
            new_point[DIRINDEX3] = plane_position;
            Point2DType normal = {0, 1};

            // std::cout << "new_point: " << new_point << std::endl;
            // std::cout << mpTrimmedDomain->IsInsideTrimmedDomain(new_point) << std::endl;
            // std::cout << pOperator->IsInside(new_point) << std::endl;

            bool Success = true;
            if ( mpTrimmedDomain->IsInsideTrimmedDomain(new_point, Success) ) { //mpTrimmedDomain->IsInsideTrimmedDomain(new_point)) {
                InsertEdge({v_left, mUpperBound[DIRINDEX2]}, {v_right, mUpperBound[DIRINDEX2]}, normal, Orientation::Positive);
            }
            if( !Success ){
                if( pOperator->IsInside(new_point) ){
                     InsertEdge({v_left, mUpperBound[DIRINDEX2]}, {v_right, mUpperBound[DIRINDEX2]}, normal, Orientation::Positive);
                }
            }
        //}

    }
    else if( interected_vertices.size() > 1 ){
        for (IndexType i = 0; i < interected_vertices.size() - 1; ++i) {
            const double v_left = interected_vertices[i].first;
            const double v_right = interected_vertices[i + 1].first;
            Point2DType normal = {0, 1};
            if( interected_vertices[i].second && !interected_vertices[i + 1].second  ){
                InsertEdge({v_left, mUpperBound[DIRINDEX2]}, {v_right, mUpperBound[DIRINDEX2]}, normal, Orientation::Positive);
            }
        }
    }

    // std::cout << mEdgesPositiveOriented.size() << std::endl;
    // for( auto edge : mEdgesPositiveOriented){
    //     std::cout << std::setprecision(16) << mVerticesPositive[edge.V1()][0] << ", " << mVerticesPositive[edge.V2()][0] << std::endl;
    // }
    // std::cout << mEdgesNegativeOriented.size() << std::endl;
    // for( auto edge : mEdgesNegativeOriented){
    //     std::cout << std::setprecision(16) << mVerticesNegative[edge.V1()][0] << ", " << mVerticesNegative[edge.V2()][0] << std::endl;
    // }

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
    // std::cout << "Wieso geth das nciht??? \n";
    // std::cout << mEdgesPositiveOriented.size() << std::endl;
    // for( auto edge : mEdgesPositiveOriented){
    //     std::cout << std::setprecision(16) << mVerticesPositive[edge.V1()][0] << ", " << mVerticesPositive[edge.V2()][0] << std::endl;
    // }
    // std::cout << mEdgesNegativeOriented.size() << std::endl;
    // for( auto edge : mEdgesNegativeOriented){
    //     std::cout << std::setprecision(16) << mVerticesNegative[edge.V1()][0] << ", " << mVerticesNegative[edge.V2()][0] << std::endl;
    // }

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
        const auto& r_edge = r_edges_origin[i];

        int edge_id_dest = FindIntersectingEdge(v1_up, v2_up, r_edge.Normal(), orientation_dest);

        bool skip = false;
        if (edge_id_dest > -1) { // Lower edge is found.
            const auto &v1_low = V1byEdgeId(edge_id_dest, Orientation::Negative);
            const auto &v2_low = V2byEdgeId(edge_id_dest, Orientation::Negative);

            /// If v1 of lower and upper edge coincide.
            /// Add:    /|
            ///        / |
            ///         \|
            if (std::abs(v1_low[1] - v1_up[1]) < mSnapTolerance) {
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

            else if (std::abs(v2_low[1] - v2_up[1]) < mSnapTolerance) {
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
    //MeshUtilities::Refine(*p_new_mesh, 10);
    return std::move(p_new_mesh);
}

void TrimmedDomainOnPlane::FindAllIntersectionWithUpperBound(std::vector<Egde2D> &rEdges, std::vector<Point2DType> &rPoints, OrientationType Orientation ) {

    auto& r_edges = GetEdges(Orientation);
    for (IndexType edge_id = 0; edge_id < GetNumberEdges(Orientation); ++edge_id) {
        const auto &v1 = V1byEdgeId(edge_id, Orientation);
        const auto &v2 = V2byEdgeId(edge_id, Orientation);
        bool v1_on_edge = false;
        if ( std::abs(v1[1] - mUpperBound[DIRINDEX2]) < 10.0*mSnapTolerance )
        {
            rPoints.push_back(v1);
            v1_on_edge = true;
        }
        bool v2_on_edge = false;
        if (std::abs(v2[1] - mUpperBound[DIRINDEX2]) < 10.0*mSnapTolerance )
        {

            rPoints.push_back(v2);
            v2_on_edge = true;

        }

        if( v1_on_edge || v2_on_edge ){
            auto new_edge = r_edges[edge_id];
            new_edge.SetVerticesOnUpperBoundary(v1_on_edge, v2_on_edge);
            rEdges.push_back( new_edge );
        }
    }

    const auto& r_vertices = GetVertices(Orientation);
    // for( auto edge : rEdges){
    //     std::cout << edge.GetVerticesOnUpperBoundary().first << ", " << edge.GetVerticesOnUpperBoundary().second << std::endl;
    //     std::cout << r_vertices[edge.V1()][DIRINDEX1] << ", " << r_vertices[edge.V2()][DIRINDEX1] << std::endl;
    // }
    // Sort Edges from left to right
    const IndexType dir_index_1 = DIRINDEX1;
    std::sort( rEdges.begin(), rEdges.end(), [&r_vertices, dir_index_1](const auto& rLHs, const auto& rRHs){
        auto status_1 = rLHs.GetVerticesOnUpperBoundary();
        auto status_2 = rRHs.GetVerticesOnUpperBoundary();
        double value_left = status_1.second ? r_vertices[rLHs.V2()][0] : r_vertices[rLHs.V1()][0];
        double value_right = status_2.second ? r_vertices[rRHs.V2()][0] : r_vertices[rRHs.V1()][0];

        return value_left < value_right;
    }  );


    // for( auto edge : rEdges){
    //     std::cout << edge.GetVerticesOnUpperBoundary().first << ", " << edge.GetVerticesOnUpperBoundary().second << std::endl;
    //     std::cout << r_vertices[edge.V1()][DIRINDEX1] << ", " << r_vertices[edge.V2()][DIRINDEX1] << std::endl;
    // }

}

bool TrimmedDomainOnPlane::DoesIntersectEdge(const Point2DType &rV1, const Point2DType &rV2, OrientationType Orientation) const {
    // Get center
    Point2DType c_positive = {0.5 * (rV1[0] + rV2[0]), 0.5 * (rV1[1] + rV2[1])};
    double min_distance = MAXD;

    for (int i = 0; i < GetNumberEdges(Orientation); ++i) {
        const auto &v1 = V1byEdgeId(i, Orientation);
        const auto &v2 = V2byEdgeId(i, Orientation);

        if (c_positive[0] > v1[0]-mSnapTolerance && c_positive[0] < v2[0]+mSnapTolerance ) {
            if( v1[1] > mUpperBound[DIRINDEX2]-mSnapTolerance &&  v2[1] > mUpperBound[DIRINDEX2]-mSnapTolerance )
                return true;

        }
    }

    return false;
}

int TrimmedDomainOnPlane::FindIntersectingEdge(const Point2DType &rV1, const Point2DType &rV2, const Point2DType &rNormal, OrientationType Orientation) const {
    // Get center
    Point2DType c_positive = {0.5 * (rV1[0] + rV2[0]), 0.5 * (rV1[1] + rV2[1])};
    double min_distance = MAXD;
    IndexType found_id = -1;

    for (int i = 0; i < GetNumberEdges(Orientation); ++i) {
        const auto &v1 = V1byEdgeId(i, Orientation);
        const auto &v2 = V2byEdgeId(i, Orientation);

        if (c_positive[0] > v1[0] && c_positive[0] < v2[0] ) {
            // TODO: add comment and use normal in which one is on top of other one!
            if( std::abs(rV1[0] - v1[0]) < 1000*mSnapTolerance && std::abs(rV2[0] - v2[0]) < 1000*mSnapTolerance) {
                Point2DType c_negative = {0.5 * (v1[0] + v2[0]), 0.5 * (v1[1] + v2[1])};
                // Vector from center to center
                Point2DType c_pos_c_neg = {c_negative[0] - c_positive[0], c_negative[1] - c_positive[1]  };
                // Dot produced
                double value = c_pos_c_neg[0]*rNormal[0] + c_pos_c_neg[1]*rNormal[1];

                bool is_inside = true;
                if( std::abs(value) <= ZEROTOL ){
                    Point3DType test_point{};
                    test_point[DIRINDEX1] = c_positive[0];
                    test_point[DIRINDEX2] = c_positive[1];
                    double plane_position = GetPlanePosition();

                    test_point[DIRINDEX3] = mUpperBoundary ? plane_position + 100.0*mSnapTolerance : plane_position - 100.0*mSnapTolerance;

                    bool is_inside = mpTrimmedDomain->IsInsideTrimmedDomain(test_point);
                    //std::cout << test_point << ": " << is_inside << std::endl;
                    if ( !is_inside ){
                        const double distance = c_positive[1] - 0.5*(v1[1] + v2[1]);
                        // Distance must be larger than 0.0 and large than alreade found min_distance.
                        if (distance > -ZEROTOL && distance < min_distance) {
                            found_id = i;
                            min_distance = distance;
                        }
                    }
                }
                else if( value < ZEROTOL ){
                    const double distance = c_positive[1] - 0.5*(v1[1] + v2[1]);
                    // Distance must be larger than 0.0 and large than alreade found min_distance.

                    if (distance >= -ZEROTOL && distance < min_distance) {
                        found_id = i;
                        min_distance = distance;
                    }
                }
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
        if (pos1 > v1[0] - mSnapTolerance && pos1 < v2[0] + mSnapTolerance)
        {
            pos2 = v1[1] + (v2[1] - v1[1]) / (v2[0] - v1[0]) * (rPoint[0] - v1[0]);
            double distance = (!is_positive) ? (rPoint[1] - pos2) : (pos2 - rPoint[1]);
            if ((distance < current_distance) && distance > mSnapTolerance)
            {
                current_distance = distance;
            }
        }
    }

    double current_lenght = MAXD;

    for (IndexType i = 0; i < r_edges_dest.size(); ++i)
    {
        const auto &v1 = V1byEdgeId(i, OrientationDest);
        const auto &v2 = V2byEdgeId(i, OrientationDest);

        double pos1 = rPoint[0];
        double pos2 = 0.0;

        if (pos1 >= v1[0] + mSnapTolerance && pos1 <= v2[0] - mSnapTolerance)
        {
            pos2 = v1[1] + (v2[1] - v1[1]) / (v2[0] - v1[0]) * (rPoint[0] - v1[0]);
            double distance = (!is_positive) ? (rPoint[1] - pos2) : (pos2 - rPoint[1]);

            if ((distance < current_distance-ZEROTOL) && distance > mSnapTolerance)
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

        auto normal = edge.Normal();
        const IndexType num_split_points = split_points.size();
        if (num_split_points > 0) {
            std::sort(split_points.begin(), split_points.end(), [](auto &point_a, auto &point_b) -> bool
                { return point_a[0] < point_b[0]; });

            /// Insert egde (vertex V1 + first split point)
            IndexType index1 = edge.V1();
            IndexType index2 = InsertVertex(split_points[0], Orientation);
            r_edges.push_back(Egde2D(index1, index2, normal) );

            /// Insert egdes (only split points.)
            for (int j = 0; j < num_split_points - 1; ++j) {
                index1 = InsertVertex(split_points[j], Orientation);
                index2 = InsertVertex(split_points[j + 1], Orientation);
                r_edges.push_back(Egde2D(index1, index2, normal));
            }

            /// Insert egde (last split point + vertex V2)
            index1 = InsertVertex(split_points[num_split_points - 1], Orientation);
            index2 = edge.V2();
            r_edges.push_back(Egde2D(index1, index2, normal));

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

bool TrimmedDomainOnPlane::InsertEdge(const Point2DType& rV1, const Point2DType& rV2, const Point2DType& rNormal, OrientationType Orientation ){
    // Get unique vertex indices.
    auto indices = GetUniqueVertexIDs(rV1, rV2, Orientation);
    // Only insert if indices are not the same.
    if (indices.first != indices.second) {
        InsertVertex(rV1, indices.first, Orientation);
        InsertVertex(rV2, indices.second, Orientation);
        auto& r_edges = GetEdges(Orientation);
        r_edges.push_back(Egde2D(indices.first, indices.second, rNormal));
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

    Point2DType normal;
    normal[0] = rNormal[DIRINDEX1];
    normal[1] = rNormal[DIRINDEX2];
    if( std::abs( v1[0] - v2[0] ) > mSnapTolerance ){ //|| std::abs( v1[1] - v2[1] ) > mSnapTolerance ){
        if ((rNormal[DIRINDEX2] > ZEROTOL) ) { // Positive oriented
            InsertEdge(v1, v2, normal, Orientation::Positive);
        }
        else if ((rNormal[DIRINDEX2] < -ZEROTOL) ) { // Negative oriented
            InsertEdge(v1, v2, normal, Orientation::Negative);
        }
        else {
            // Vertical oriented
            if( std::abs( v1[1] - v2[1] ) > mSnapTolerance ){
                InsertEdge(v1, v2, normal, Orientation::Vertical);
            }
        }
    } else {
        // Vertical oriented
        if( std::abs( v1[1] - v2[1] ) > mSnapTolerance ){
            InsertEdge(v1, v2, normal, Orientation::Vertical);
        }
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
        //TIBRA_ERROR_IF("TrimmedDomainOnPlane::InsertVertex", !res.second) << "Vetrex already exists.\n";
    }
    else if(NewIndex > index) {
        TIBRA_ERROR("TrimmedDomainOnPlane::InsertVertex") << "Given index out of range.\n";
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
    if( !PointComparison(mSnapTolerance)(tmp_vector.begin(), ++tmp_vector.begin())  &&
            !PointComparison(mSnapTolerance)(++tmp_vector.begin(), tmp_vector.begin()) ){
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
        TIBRA_ERROR("TrimmedDomainOnPlane::V1byEdgeId") << "Given Orientation not available.\n";
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
        TIBRA_ERROR("TrimmedDomainOnPlane::V2byEdgeId") << "Given Orientation not available.\n";
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
        TIBRA_ERROR("TrimmedDomainOnPlane::GetEdges") << "Given Orientation not available.\n";
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
        TIBRA_ERROR("TrimmedDomainOnPlane::GetEdges") << "Given Orientation not available.\n";
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
        TIBRA_ERROR("TrimmedDomainOnPlane::GetVertices") << "Given Orientation not available.\n";
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
        TIBRA_ERROR("TrimmedDomainOnPlane::GetVertices") << "Given Orientation not available.\n";
    }
}

Point2DSetType& TrimmedDomainOnPlane::GetVerticesSet(OrientationType Orientation) {
    switch (Orientation)
    {
    case Orientation::Positive:
        return *mVerticesSetPositive.get();
    case Orientation::Negative:
        return *mVerticesSetNegative.get();
    case Orientation::Vertical:
        return *mVerticesSetVertical.get();
    default:
        TIBRA_ERROR("TrimmedDomainOnPlane::GetVerticesSet") << "Given Orientation not available.\n";
        break;
    }
}

const Point2DSetType& TrimmedDomainOnPlane::GetVerticesSet(OrientationType Orientation) const {
    switch (Orientation)
    {
    case Orientation::Positive:
        return *mVerticesSetPositive.get();
    case Orientation::Negative:
        return *mVerticesSetNegative.get();
    case Orientation::Vertical:
        return *mVerticesSetVertical.get();
    default:
        TIBRA_ERROR("TrimmedDomainOnPlane::GetVerticesSet") << "Given Orientation not available.\n";
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
        TIBRA_ERROR("TrimmedDomainOnPlane::GetNumberEdges") << "Given Orientation not available.\n";
        break;
    }
}

double TrimmedDomainOnPlane::GetPlanePosition() const {
    return mUpperBoundary ? mUpperBound[DIRINDEX3] : mLowerBound[DIRINDEX3];
}



} // End namespace tibra