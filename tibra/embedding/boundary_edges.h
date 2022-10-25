// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef EDGE_2D_INCLUDE_H
#define EDGE_2D_INCLUDE_H

/// External libraries
#include <cstddef>
#include <array>

///Project includes

///@name TIBRA Classes
///@{

/**
 * @class  BoundaryEdges
 * @author Manuel Messmer
 * @brief Hallo
 * DIRINDEX1 plane index 1.
 * DIRINDEX2 plane index 2.
*/
template<std::size_t DIRINDEX1, std::size_t DIRINDEX2>
class BoundaryEdges {

public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;
    typedef std::array<double, 2> Point2DType;
    typedef std::array<double, 3> Point3DType;

    /// Default constructor.
    BoundaryEdges() {
        if( DIRINDEX1 == DIRINDEX2){
            throw std::runtime_error("BoundaryEdges :: Contructor :: Template arguments must not be the same.");
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

        Egde2D( const Point2DType& V1, const Point2DType& V2 ) :
            mV1(V1), mV2(V2)
        {}

        const Point2DType& V1(){
            return V1;
        }
        const Point2DType& V2(){
            return V2;
        }

    private:
        Point2DType mV1{};
        Point2DType mV2{};
    };

    ///@}
    void Reserve(IndexType Size){
        mEdgesPositiveOriented.reserve(Size);
        mEdgesNegativeOriented.reserve(Size);
    }

    void AddEdge( const Point3DType& V1, const Point3DType& V2, const Point3DType& rNormal ) {
        Point2DType v1 = {V1[DIRINDEX1], V1[DIRINDEX2]};
        Point2DType v2 = {V2[DIRINDEX1], V2[DIRINDEX2]};

        bool is_positive_oriented = rNormal[DIRINDEX1] > 0.0;

        if( is_positive_oriented ){
            mEdgesPositiveOriented.push_back( Egde2D(v1, v2) );
        } else {
            mEdgesNegativeOriented.push_back( Egde2D(v1, v2) );
        }
    }

    void RefineEdges(){
        IndexType num_positive_edges = mEdgesPositiveOriented.size();
        IndexType num_negative_edges = mEdgesNegativeOriented.size();
        if( num_positive_edges == 0 ||  num_negative_edges == 0 ){
            return;
        }

        std::map<IndexType, double> index_map_1{};
        for( int i = 0; i < num_positive_edges; ++i){

            const auto& v1 = mEdgesPositiveOriented[i].V1();
            const auto& v1 = mEdgesPositiveOriented[i].V2();
        }


    }

    Point2DType FindIntersection(Point2DType, const std::vector<Egde2D>& rEdges) const{

    }

private:
    std::vector<Egde2D> mEdgesPositiveOriented{};
    std::vector<Egde2D> mEdgesNegativeOriented{};
};

#endif