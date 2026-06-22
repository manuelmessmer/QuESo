//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// Project includes
#include "queso/containers/boundary_integration_point.hpp"
#include "queso/containers/trimmed_element.hpp"
#include "queso/containers/integration_point.hpp"
#include "queso/embedding/octree.h"
#include "queso/embedding/trimmed_domain.h"

namespace queso {

///
/// Function Definition of Octree<TOperator>::Node
///

template<typename TOperator>
void Octree<TOperator>::Node::Refine(IndexType MinLevel, IndexType MaxLevel, const TOperator* pOperator){
    const PointType& r_lower_bound_xyz = mBoundsXYZ.lower;
    const PointType& r_upper_bound_xyz = mBoundsXYZ.upper;
    const PointType& r_lower_bound_uvw = mBoundsUVW.lower;
    const PointType& r_upper_bound_uvw = mBoundsUVW.upper;
    if( this->IsLeaf() ){
        const auto delta_xyz = (0.5 * (r_upper_bound_xyz - r_lower_bound_xyz));
        const PointType delta_x_xyz{delta_xyz[0], 0.0, 0.0};
        const PointType delta_y_xyz{0.0, delta_xyz[1], 0.0};
        const PointType delta_z_xyz{0.0, 0.0, delta_xyz[2]};
        const auto delta_uvw = (0.5 * (r_upper_bound_uvw - r_lower_bound_uvw));
        const PointType delta_x_uvw{delta_uvw[0], 0.0, 0.0};
        const PointType delta_y_uvw{0.0, delta_uvw[1], 0.0};
        const PointType delta_z_uvw{0.0, 0.0, delta_uvw[2]};

        //       d_________c
        //      /        /|                 y
        //     /        / |                ´|`
        //   h/_______g/  |b                |-->x
        //    |  a     |  /                /
        //    |        | /                Z
        //    |________|/
        //    e        f
        //
        if( (mLevel < MinLevel) || (mLevel < MaxLevel && mStatus == IntersectionState::trimmed) ){
            // Corner a
            CreateNewNode(MinLevel, MaxLevel, 0, MakeBox(r_lower_bound_xyz, (r_lower_bound_xyz + delta_xyz) ),
                                                 MakeBox(r_lower_bound_uvw, (r_lower_bound_uvw + delta_uvw) ), pOperator);
            // Corner b (a+delta_x)
            CreateNewNode(MinLevel, MaxLevel, 1, MakeBox((r_lower_bound_xyz + delta_x_xyz), ((r_lower_bound_xyz + delta_x_xyz) + delta_xyz) ),
                                                 MakeBox((r_lower_bound_uvw + delta_x_uvw), ((r_lower_bound_uvw + delta_x_uvw) + delta_uvw) ), pOperator);
            // Corner c (g-delta_z)
            CreateNewNode(MinLevel, MaxLevel, 2, MakeBox(((r_upper_bound_xyz - delta_z_xyz) - delta_xyz), (r_upper_bound_xyz - delta_z_xyz) ),
                                                 MakeBox(((r_upper_bound_uvw - delta_z_uvw) - delta_uvw), (r_upper_bound_uvw - delta_z_uvw) ), pOperator);
            // Corner d (a+delta_y)
            CreateNewNode(MinLevel, MaxLevel, 3, MakeBox((r_lower_bound_xyz + delta_y_xyz), ((r_lower_bound_xyz + delta_y_xyz) + delta_xyz) ),
                                                 MakeBox((r_lower_bound_uvw + delta_y_uvw), ((r_lower_bound_uvw + delta_y_uvw) + delta_uvw) ), pOperator);
            // Corner e (a+delta_z)
            CreateNewNode(MinLevel, MaxLevel, 4, MakeBox((r_lower_bound_xyz + delta_z_xyz), ((r_lower_bound_xyz + delta_z_xyz) + delta_xyz) ),
                                                 MakeBox((r_lower_bound_uvw + delta_z_uvw), ((r_lower_bound_uvw + delta_z_uvw) + delta_uvw) ), pOperator);
            // Corner f (g-delta_y)
            CreateNewNode(MinLevel, MaxLevel, 5, MakeBox(((r_upper_bound_xyz - delta_y_xyz) - delta_xyz), (r_upper_bound_xyz - delta_y_xyz) ),
                                                 MakeBox(((r_upper_bound_uvw - delta_y_uvw) - delta_uvw), (r_upper_bound_uvw - delta_y_uvw) ), pOperator);
            // Corner g
            CreateNewNode(MinLevel, MaxLevel, 6, MakeBox((r_upper_bound_xyz - delta_xyz), r_upper_bound_xyz),
                                                 MakeBox((r_upper_bound_uvw - delta_uvw), r_upper_bound_uvw), pOperator);
            // Corner h (g-delta_x)
            CreateNewNode(MinLevel, MaxLevel, 7, MakeBox(((r_upper_bound_xyz - delta_x_xyz) - delta_xyz), (r_upper_bound_xyz - delta_x_xyz) ),
                                                 MakeBox(((r_upper_bound_uvw - delta_x_uvw) - delta_uvw), (r_upper_bound_uvw - delta_x_uvw) ), pOperator);

        }
    }
    else {
        for( IndexType i = 0; i < 8UL; ++i){
            if( mChildren[i] ){
                mChildren[i]->Refine(MinLevel, MaxLevel, pOperator);
            }
        }
    }

}

template<typename TOperator>
void Octree<TOperator>::Node::NumberOfLeafs(IndexType& rValue) const {
    if( this->IsLeaf() ){
        ++rValue;
    } else {
        for( IndexType i = 0; i < 8UL; ++i){
            if( mChildren[i] ){ // If not nullptr
                mChildren[i]->NumberOfLeafs(rValue);
            }
        }
    }
}

template<typename TOperator>
void Octree<TOperator>::Node::NumberOfNodes(IndexType& rValue) const {
    ++rValue;
    if( !this->IsLeaf() ){
        for( IndexType i = 0; i < 8UL; ++i){
            if( mChildren[i] ){ // If not nullptr
                mChildren[i]->NumberOfNodes(rValue);
            }
        }
    }
}

template<typename TOperator>
void Octree<TOperator>::Node::CreateNewNode(IndexType MinLevel, IndexType MaxLevel, IndexType ChildIndex, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW, const TOperator* pOperator){
    const auto status = pOperator->template GetIntersectionState<CoordinateSpace::global>(rBoundsXYZ.lower, rBoundsXYZ.upper);
    if( status != IntersectionState::outside ){
        mChildren[ChildIndex] = MakeUnique<Node>(rBoundsXYZ, rBoundsUVW, status, mLevel+1);
        ++mNumChildren;
        mChildren[ChildIndex]->Refine(MinLevel, MaxLevel, pOperator);
    }
}

template<typename TOperator>
bool Octree<TOperator>::Node::IsLeaf() const {
    return (mNumChildren == 0UL);
}

///
/// Function Definition of Octree<TOperator>
///

template<typename TOperator>
void Octree<TOperator>::Refine(IndexType MinLevel, IndexType MaxLevel){
    QuESo_ERROR_IF( MinLevel > MaxLevel ) << "MinLevel must be smaller/equal than MaxLevel. "
        << "Given MinLevel: " << MinLevel << ", MaxLevel: " << MaxLevel << ".\n";
    mMinLevel = MinLevel;
    mMaxLevel = MaxLevel;
    mpRoot->Refine(MinLevel, MaxLevel, mpOperator);
}


template<typename TOperator>
SizeType Octree<TOperator>::NumberOfLeafs() const{
    SizeType number_of_leafs = 0UL;
    mpRoot->NumberOfLeafs(number_of_leafs);
    return number_of_leafs;
}

template<typename TOperator>
SizeType Octree<TOperator>::NumberOfNodes() const{
    SizeType number_of_nodes = 0UL;
    mpRoot->NumberOfNodes(number_of_nodes);
    return number_of_nodes;
}

// Explicit class instantiation
template class Octree<TrimmedElement<IntegrationPoint, BoundaryIntegrationPoint>>;

} // End namespace queso
