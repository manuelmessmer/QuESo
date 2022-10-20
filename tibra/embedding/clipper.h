// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CLIPPER_INCLUDE_H
#define CLIPPER_INCLUDE_H

/// External libraries
#include <cstddef>
#include <array>

///@name TIBRA Classes
///@{

template<std::size_t SIZE>
class Polygon {

public:
    typedef std::size_t IndexType;
    typedef std::array<double,3> PointType;

    void AddVertex(const PointType& rPoint){
        if( mNumVertices >= SIZE ){
            throw std::runtime_error("Polygon :: Size of Polygon is exceeded.");
        }
        mVertices[mNumVertices] = rPoint;
        ++mNumVertices;
    }

    IndexType NumVertices() const{
        return mNumVertices;
    }

    const PointType& GetVertex(IndexType i) const{
        if( i >= mNumVertices ){
            throw std::runtime_error("Polygon :: Size of Polygon is exceeded.");
        }
        return mVertices[i];
    }

    const PointType& operator[] (IndexType i) const
    {
        return mVertices[i];
    }

    void Clear(){
        std::fill(mVertices.begin(), mVertices.end(), 0.0);
        mNumVertices = 0;
    }

private:
    std::array<PointType, SIZE> mVertices{}; // Keep static array to be fast.
    IndexType mNumVertices = 0;
};

/**
 * @class  Clipper
 * @author Manuel Messmer
 * @brief  Provides functions to clip triangle by cube.
*/
class Clipper  {

public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;
    typedef std::array<double,3> PointType;
    typedef Polygon<6> PolygonType; // Intersections contains maximum 6 vertices.

    enum Side{ON_POSITIVE_SIDE, ON_NEGATIVE_SIDE, ON_BOUNDARY};

    ///@}
    ///@name Operations
    ///@{

    ///@brief Clips triangle by AABB. Expects that input triangle does intersect one of the faces of the bounding box.
    ///       Meaning, triangle is not fully outside, neither fully contained inside the AABB.
    ///@param rV1 Vertex 1 of Triangle
    ///@param rV2 Vertex 2 of Triangle
    ///@param rV3 Vertex 3 of Triangle
    ///@param rLowerBound Lower bound of AABB.
    ///@param rUpperBound Upper bound of AABB.
    ///@return std::unique_ptr<Polygon> (Will contain maximal 6 vertices).
    std::unique_ptr<Polygon<6>> Clip(const PointType& rV1, const PointType& rV2, const PointType& rV3,
                 const PointType& rLowerBound, const PointType& rUpperBound){


        PolygonType poly[2] = {PolygonType{}, PolygonType{}};
        PolygonType* currentPoly = &poly[0];
        PolygonType* prevPoly = &poly[1];

        const PointType x_values = {rV1[0], rV2[0], rV3[0]};
        const PointType y_values = {rV1[1], rV2[1], rV3[1]};
        const PointType z_values = {rV1[2], rV2[2], rV3[2]};

        auto x_min_max = std::minmax_element(x_values.begin(), x_values.end());
        auto y_min_max = std::minmax_element(y_values.begin(), y_values.end());
        auto z_min_max = std::minmax_element(z_values.begin(), z_values.end());

        PointType min_tri = {*x_min_max.first, *y_min_max.first, *z_min_max.first};
        PointType max_tri = {*x_min_max.second, *y_min_max.second, *z_min_max.second};

        //Loop through the planes of the bbox and clip the vertices
        for(IndexType dim = 0; dim < 3; ++dim)
        {
            // Optimization note: we should be able to save some work based on
            // the clipping plane and the triangle's bounding box
            if(max_tri[dim] > rLowerBound[dim])
            {
                swap(prevPoly, currentPoly);
                ClipAxisByPlane(*prevPoly, *currentPoly, 2 * dim + 0, rLowerBound[dim]);
            }

            if(min_tri[dim] < rUpperBound[dim])
            {
                swap(prevPoly, currentPoly);
                ClipAxisByPlane(prevPoly, currentPoly, 2 * dim + 1, rUpperBound[dim]);
            }
        }

        return std::make_unique<Polygon<6>>(*currentPoly);
    }

    ///@}

private:

    template <typename T>
    inline void swap(T& a, T& b)
    {
        T tmp = a;
        a = b;
        b = tmp;
    }

    void ClipAxisByPlane(const PolygonType* rPrevPoly,
                         PolygonType* rCurrentPoly,
                         IndexType Index,
                         double Val)
    {

        rCurrentPoly->Clear();
        int numVerts = rPrevPoly->NumVertices();

        if(numVerts == 0)
        {
            return;
        }

        // Initialize point a with the last vertex of the polygon
        const PointType* a = &(*rPrevPoly)[numVerts - 1];
        int aSide = ClassifyPointAxisPlane(*a, Index, Val);

        for(int i = 0; i < numVerts; ++i)
        {
            const PointType* b = &(*rPrevPoly)[i];
            int bSide = ClassifyPointAxisPlane(*b, Index, Val);

            switch(bSide)
            {
            case ON_POSITIVE_SIDE:
            if(aSide == ON_NEGATIVE_SIDE)
            {
                rCurrentPoly->AddVertex(FindIntersectionPoint(*a, *b, Index, Val));
            }
            break;
            case ON_BOUNDARY:
            if(aSide == ON_NEGATIVE_SIDE)
            {
                rCurrentPoly->AddVertex(*b);
            }
            break;
            case ON_NEGATIVE_SIDE:
            switch(aSide)
            {
            case ON_POSITIVE_SIDE:
                rCurrentPoly->AddVertex(FindIntersectionPoint(*a, *b, Index, Val));
                rCurrentPoly->AddVertex(*b);
                break;
            case ON_BOUNDARY:
                rCurrentPoly->AddVertex(*a);
                rCurrentPoly->AddVertex(*b);
                break;
            case ON_NEGATIVE_SIDE:
                rCurrentPoly->AddVertex(*b);
                break;
            }
            break;
            }

            // swap a and b
            a = b;
            aSide = bSide;
        }
    }

    int ClassifyPointAxisPlane(const PointType& rPoint,
                               IndexType Index,
                               double Val,
                               const double Eps = 1e-8)
    {
    // Note: we are exploiting the fact that the planes are axis aligned
    // So the dot product is +/- the given coordinate.
    // In general, we would need to call distance(pt, plane) here
    bool is_even = Index % 2;
    const double dist = is_even ? Val - rPoint[Index / 2] : rPoint[Index / 2] - Val;

    if(dist > Eps)
    {
        return ON_POSITIVE_SIDE;
    }
    if(dist < -Eps)
    {
        return ON_NEGATIVE_SIDE;
    }

    return ON_BOUNDARY;
    }

    // /*!
    // * \brief Finds the clipping intersection point between points a and b.
    // *
    // * \param [in] a The point behind the plane
    // * \param [in] b The point in front of the plane
    // * \param [in] index The index of the axis aligned plane.
    // * \param [in] val The plane's coordinate with respect to the given axis
    // * \return The point between a and b whose corresponding coordinate is val
    // *
    // * \see classifyPointAxisPlane for description of how index maps to coordinates.
    // */
    PointType FindIntersectionPoint(const PointType& rA,
                                    const PointType& rB,
                                    IndexType Index,
                                    double Val)
    {

    // Need to find a parameter t for the point pt, such that,
    // * 0 <= t <= 1
    // * pt = a + t * (b-a)
    // * pt[ index/2]  == val

    double t = (Val - rA[Index / 2]) / (rB[Index / 2] - rA[Index / 2]);

    if( !(0.0 <= t && t <= 1.0) ){
        throw std::runtime_error("Error");
    }

    PointType ret = {rA[0] + t*(rB[0] - rA[0]),
                     rA[1] + t*(rB[1] - rA[1]),
                     rA[2] + t*(rB[2] - rA[2]) };

    auto status = ClassifyPointAxisPlane(ret, Index, Val);
    if( status == ON_BOUNDARY){
        throw std::runtime_error("Point is in boudnary.");
    }

    return ret;
    }
};
///@}

#endif