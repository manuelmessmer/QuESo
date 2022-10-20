// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef CLIPPER_INCLUDE_H
#define CLIPPER_INCLUDE_H

/// External libraries
#include <cstddef>
#include <array>

#include "geometries/triangle_mesh.h"

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

    std::vector<std::array<PointType, 3>> GetTriangles() const{
        std::vector<std::array<PointType, 3>> new_triangles;
        new_triangles.reserve(mNumVertices);
        std::cout << "mNumVertices: " << mNumVertices << std::endl;
        if(mNumVertices < 3){
            throw std::runtime_error("Obacht");
        }

        if( mNumVertices == 3 ){

            new_triangles.push_back({ mVertices[0], mVertices[1], mVertices[2]} );

            return new_triangles;
        }
        // TODO: Improve this
        PointType centroid = {0.0, 0.0, 0.0};
        const double inv_num_vertices = 1.0/mNumVertices;
        for( IndexType i = 0 ; i < mNumVertices; ++i){
            centroid[0] += inv_num_vertices*mVertices[i][0];
            centroid[1] += inv_num_vertices*mVertices[i][1];
            centroid[2] += inv_num_vertices*mVertices[i][2];
        }

        for( IndexType i = 0 ; i < mNumVertices-1; ++i){
            new_triangles.push_back( {mVertices[i], mVertices[i+1], centroid} );
        }
        new_triangles.push_back( {mVertices[mNumVertices-1], mVertices[0], centroid} );

        return new_triangles;
    }

    void Clear(){
        std::fill(mVertices.begin(), mVertices.end(), PointType{});
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

    ///@brief Clips triangle mesh by AABB.
    ///@param rV1 Vertex 1 of Triangle
    ///@param rV2 Vertex 2 of Triangle
    ///@param rV3 Vertex 3 of Triangle
    ///@param rLowerBound Lower bound of AABB.
    ///@param rUpperBound Upper bound of AABB.
    ///@return std::unique_ptr<Polygon> (Will contain maximal 6 vertices).
    static std::unique_ptr<TriangleMesh> ClipTriangleMesh(const TriangleMesh& rTriangleMesh,
                 const PointType& rLowerBound, const PointType& rUpperBound){

        TriangleMesh new_mesh{};
        std::map<IndexType, IndexType> index_map{};
        IndexType vertex_count = 0;
        for( IndexType triangle_id = 0; triangle_id < rTriangleMesh.NumOfTriangles(); ++triangle_id ){
            const auto& P1 = rTriangleMesh.P1(triangle_id);
            const auto& P2 = rTriangleMesh.P2(triangle_id);
            const auto& P3 = rTriangleMesh.P3(triangle_id);

            if(    IsContained(P1, rLowerBound, rUpperBound )
                && IsContained(P2, rLowerBound, rUpperBound )
                && IsContained(P3, rLowerBound, rUpperBound ) ){ // Triangle is fully contained, does not need to be clipped.

                new_mesh.AddVertex(P1);
                new_mesh.AddVertex(P2);
                new_mesh.AddVertex(P3);

                // Copy triangles and normals.
                new_mesh.AddTriangle({vertex_count, vertex_count+1, vertex_count+2 });
                vertex_count += 3;
                new_mesh.AddNormal( rTriangleMesh.Normal(triangle_id) );
            }
            else { // Triangle needs to be clipped.
                //throw std::runtime_error("Wht the fuck!!");
                auto polygon = ClipTriangle(P1, P2, P3, rLowerBound, rUpperBound);
                for( auto triangle : polygon->GetTriangles() ){
                    new_mesh.AddVertex(triangle[0]);
                    new_mesh.AddVertex(triangle[1]);
                    new_mesh.AddVertex(triangle[2]);

                    // Copy triangles and normals.
                    new_mesh.AddTriangle({vertex_count, vertex_count+1, vertex_count+2 });
                    vertex_count += 3;
                    new_mesh.AddNormal( rTriangleMesh.Normal(triangle_id) );
                }
            }

        }
        std::cout << "vertices: " << new_mesh.NumOfVertices() << std::endl;
        return std::make_unique<TriangleMesh>(new_mesh);
    }

    static bool IsContained(const PointType& rPoint, const PointType& rLowerBound, const PointType& rUpperBound){
        if(    rPoint[0] < rLowerBound[0]
            || rPoint[0] > rUpperBound[0]
            || rPoint[1] < rLowerBound[1]
            || rPoint[1] > rUpperBound[1]
            || rPoint[2] < rLowerBound[2]
            || rPoint[2] > rUpperBound[2] )
        {
            return false;
        }

        return true;
    }

    ///@brief Clips triangle by AABB. Expects that input triangle does intersect one of the faces of the bounding box.
    ///       Meaning, triangle is not fully outside, neither fully contained inside the AABB.
    ///@param rV1 Vertex 1 of Triangle
    ///@param rV2 Vertex 2 of Triangle
    ///@param rV3 Vertex 3 of Triangle
    ///@param rLowerBound Lower bound of AABB.
    ///@param rUpperBound Upper bound of AABB.
    ///@return std::unique_ptr<Polygon> (Will contain maximal 6 vertices).
    static std::unique_ptr<Polygon<6>> ClipTriangle(const PointType& rV1, const PointType& rV2, const PointType& rV3,
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

        currentPoly->AddVertex(rV1);
        currentPoly->AddVertex(rV2);
        currentPoly->AddVertex(rV3);

        std::cout << "current POly: " << currentPoly->NumVertices() << std::endl;
        std::cout << "orev POly: " << prevPoly->NumVertices() << std::endl;
        //swap(prevPoly, currentPoly);
        std::cout << "current POly: " << currentPoly->NumVertices() << std::endl;
        std::cout << "orev POly: " << prevPoly->NumVertices() << std::endl;
        //Loop through the planes of the bbox and clip the vertices
        for(IndexType dim = 0; dim < 3; ++dim)
        {
            // Optimization note: we should be able to save some work based on
            // the clipping plane and the triangle's bounding box
            if(max_tri[dim] > rLowerBound[dim])
            {
                swap(prevPoly, currentPoly);
                ClipAxisByPlane(prevPoly, currentPoly, 2 * dim + 0, rLowerBound[dim]);
            }

            if(min_tri[dim] < rUpperBound[dim])
            {
                swap(prevPoly, currentPoly);
                ClipAxisByPlane(prevPoly, currentPoly, 2 * dim + 1, rUpperBound[dim]);
            }
        }
        std::cout << "in the end: " << currentPoly->NumVertices() << std::endl;
        return std::make_unique<Polygon<6>>(*currentPoly);
    }


    ///@}

private:

    template <typename T>
    static inline void swap(T& a, T& b)
    {
        T tmp = a;
        a = b;
        b = tmp;
    }

    static void ClipAxisByPlane(const PolygonType* rPrevPoly,
                         PolygonType* rCurrentPoly,
                         IndexType Index,
                         double Val)
    {

        rCurrentPoly->Clear();
        int numVerts = rPrevPoly->NumVertices();
        std::cout << "in here: " << numVerts << std::endl;
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

            std::cout << "aSide: " << aSide << std::endl;
            std::cout << "bSide: " << bSide << std::endl;
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
        std::cout << "hhhhhhh: " << rCurrentPoly->NumVertices() << std::endl;
    }

    static int ClassifyPointAxisPlane(const PointType& rPoint,
                               IndexType Index,
                               double Val,
                               const double Eps = 1e-10)
    {
    // Note: we are exploiting the fact that the planes are axis aligned
    // So the dot product is +/- the given coordinate.
    // In general, we would need to call distance(pt, plane) here
    bool is_even =  (Index & 1) == 0;

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
    static PointType FindIntersectionPoint(const PointType& rA,
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
    // if( status == ON_BOUNDARY){
    //     throw std::runtime_error("Point is on boudnary.");
    // }

    return ret;
    }
};
///@}

#endif