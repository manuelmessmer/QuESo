// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#include "embedding/clipper.h"
#include "utilities/utilities.h"

typedef std::size_t IndexType;
typedef std::array<double, 3> PointType;

// Function Definitions Polygon
template<IndexType SIZE>
void Polygon<SIZE>::AddVertex(const PointType& rPoint) {
    if( mNumVertices >= SIZE ){
        throw std::runtime_error("Polygon :: AddVertex :: Size of Polygon is exceeded.");
    }
    mVertices[mNumVertices] = rPoint;
    ++mNumVertices;
}

template<IndexType SIZE>
IndexType Polygon<SIZE>::NumVertices() const {
    return mNumVertices;
}

template<IndexType SIZE>
const PointType& Polygon<SIZE>::GetVertex(IndexType i) const {
    if( i >= mNumVertices ){
        throw std::runtime_error("Polygon :: GetVertex :: Size of Polygon is exceeded.");
    }
    return mVertices[i];
}

template<IndexType SIZE>
const PointType& Polygon<SIZE>::operator[] (IndexType i) const {
    return mVertices[i];
}

template<IndexType SIZE>
std::unique_ptr<std::vector<std::array<PointType, 3>>> Polygon<SIZE>::pGetTriangles() const {

    auto p_new_triangles = std::make_unique<std::vector<std::array<PointType, 3>>>();
    p_new_triangles->reserve(mNumVertices);

    if(mNumVertices < 3){
        throw std::runtime_error("Obacht");
    }

    if( mNumVertices == 3 ){

        p_new_triangles->push_back({ mVertices[0], mVertices[1], mVertices[2]} );

        return std::move(p_new_triangles);
    }

    // Compute mean of vertices
    PointType centroid = {0.0, 0.0, 0.0};
    const double inv_num_vertices = 1.0/mNumVertices;

    for( IndexType i = 0 ; i < mNumVertices; ++i){
        centroid[0] += inv_num_vertices*mVertices[i][0];
        centroid[1] += inv_num_vertices*mVertices[i][1];
        centroid[2] += inv_num_vertices*mVertices[i][2];
    }

    for( IndexType i = 0 ; i < mNumVertices-1; ++i){
        p_new_triangles->push_back( {mVertices[i], mVertices[i+1], centroid} );
    }
    p_new_triangles->push_back( {mVertices[mNumVertices-1], mVertices[0], centroid} );

    return std::move(p_new_triangles);
}

template<IndexType SIZE>
const PointType& Polygon<SIZE>::GetLastVertex() const{
    return mVertices[mNumVertices-1];
}

template<IndexType SIZE>
void Polygon<SIZE>::Clear(){
    std::fill(mVertices.begin(), mVertices.end(), PointType{});
    mNumVertices = 0;
}

// Explicit instantiation Polygon
template class Polygon<9>;

// Function Definitions Clipper
typedef Clipper::PolygonType PolygonType;

std::unique_ptr<PolygonType> Clipper::ClipTriangle(const PointType& rV1, const PointType& rV2, const PointType& rV3,
                const PointType& rLowerBound, const PointType& rUpperBound){

    std::unique_ptr<PolygonType> p_current_poly = std::make_unique<PolygonType>();
    std::unique_ptr<PolygonType> p_prev_poly = std::make_unique<PolygonType>();

    const PointType x_values = {rV1[0], rV2[0], rV3[0]};
    const PointType y_values = {rV1[1], rV2[1], rV3[1]};
    const PointType z_values = {rV1[2], rV2[2], rV3[2]};

    auto x_min_max = std::minmax_element(x_values.begin(), x_values.end());
    auto y_min_max = std::minmax_element(y_values.begin(), y_values.end());
    auto z_min_max = std::minmax_element(z_values.begin(), z_values.end());

    PointType min_tri = {*x_min_max.first, *y_min_max.first, *z_min_max.first};
    PointType max_tri = {*x_min_max.second, *y_min_max.second, *z_min_max.second};

    p_current_poly->AddVertex(rV1);
    p_current_poly->AddVertex(rV2);
    p_current_poly->AddVertex(rV3);

    /// Loop ober all three dimensions.
    for(IndexType dimension = 0; dimension < 3; ++dimension)
    {
        if(max_tri[dimension] > rLowerBound[dimension])
        {
            utilities::swap_unique(p_prev_poly, p_current_poly);
            // Plane is oriented in negative direction.
            Plane plane(2 * dimension + 0, rLowerBound[dimension]);
            ClipPolygonByPlane(p_prev_poly.get(), p_current_poly.get(), plane);
        }

        if(min_tri[dimension] < rUpperBound[dimension])
        {
            utilities::swap_unique(p_prev_poly, p_current_poly);
            // Plane is positive in negative direction.
            Plane plane(2 * dimension + 1, rUpperBound[dimension]);
            ClipPolygonByPlane(p_prev_poly.get(), p_current_poly.get(), plane);
        }
    }
    return p_current_poly;
}

void Clipper::ClipPolygonByPlane(const PolygonType* pPrevPoly,
                         PolygonType* pCurrentPoly,
                         const Plane& rPlane) {

    pCurrentPoly->Clear();
    int numVerts = pPrevPoly->NumVertices();

    if(numVerts == 0) {
        throw std::runtime_error("Clipper :: ClipAxisByPlane :: pPrevPoly contains 0 vertices.");
    }

    // Initialize point a with the last vertex of the polygon
    const PointType* a = &pPrevPoly->GetLastVertex();
    int aSide = ClassifyPointToPlane(*a, rPlane);

    for(int i = 0; i < numVerts; ++i)
    {
        const PointType& b = pPrevPoly->GetVertex(i);
        int bSide = ClassifyPointToPlane(b, rPlane);

        switch(bSide)
        {
        case IN_FRONT_OF_PLANE:
            if(aSide == BEHIND_PLANE) {
                pCurrentPoly->AddVertex(FindIntersectionPointOnPlane(*a, b, rPlane));
            }
            break;
        case ON_PLANE:
            if(aSide == BEHIND_PLANE) {
                pCurrentPoly->AddVertex(b);
            }
            break;
        case BEHIND_PLANE:
            switch(aSide)
            {
            case IN_FRONT_OF_PLANE:
                pCurrentPoly->AddVertex(FindIntersectionPointOnPlane(*a, b, rPlane));
                pCurrentPoly->AddVertex(b);
                break;
            case ON_PLANE:
                pCurrentPoly->AddVertex(*a);
                pCurrentPoly->AddVertex(b);
                break;
            case BEHIND_PLANE:
                pCurrentPoly->AddVertex(b);
                break;
            }
            break;
            }

        // swap a and b
        a = &b;
        aSide = bSide;
    }
}

IndexType Clipper::ClassifyPointToPlane(const PointType& rPoint,
                            const Plane& rPlane,
                            const double Eps) {

    // Note: Plane is axis aligned. Distance is directly given by coordinate.
    const bool is_negative_oriented = rPlane.IsNegativeOriented();
    const IndexType index = rPlane.GetIndexOfNormalDirection();

    const double dist = is_negative_oriented ? rPlane.mPosition - rPoint[index] : rPoint[index] - rPlane.mPosition;

    if(dist > Eps) {
        return IN_FRONT_OF_PLANE;
    }
    if(dist < -Eps) {
        return BEHIND_PLANE;
    }

    return ON_PLANE;
}

PointType Clipper::FindIntersectionPointOnPlane(const PointType& rA,
                                                const PointType& rB,
                                                const Plane& rPlane) {

    // Need to find a parameter t for the point pt, such that,
    // * 0 <= t <= 1
    // * pt = rA + t * (rB-rA)
    // * pt[index]  == val

    IndexType index = rPlane.GetIndexOfNormalDirection();
    double t = (rPlane.mPosition - rA[index]) / (rB[index] - rA[index]);

    if( !(0.0 <= t && t <= 1.0) ){
        throw std::runtime_error("Error");
    }

    PointType ret = {rA[0] + t*(rB[0] - rA[0]),
                     rA[1] + t*(rB[1] - rA[1]),
                     rA[2] + t*(rB[2] - rA[2]) };

    auto status = ClassifyPointToPlane(ret, rPlane);

    if( status != ON_PLANE){
        throw std::runtime_error("Intersection Point is not on Boundary");
    }

    return ret;
}
