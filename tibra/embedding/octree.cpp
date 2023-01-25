// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// Project includes
#include "embedding/octree.h"
#include "embedding/trimmed_domain_base.h"

namespace tibra {

///
/// Function Definition of Octree<TOperator>::Node
///

template<typename TOperator>
void Octree<TOperator>::Node::Refine(IndexType MinLevel, IndexType MaxLevel, const TOperator* pOperator){
    if( this->IsLeaf() ){
        const auto delta = (mUpperBound - mLowerBound) * 0.5;
        const PointType delta_x(delta[0], 0.0, 0.0);
        const PointType delta_y(0.0, delta[1], 0.0);
        const PointType delta_z(0.0, 0.0, delta[2]);
        //       d_________c
        //      /        /|                 y
        //     /        / |                ´|`
        //   h/_______g/  |b                |-->x
        //    |  a     |  /                /
        //    |        | /                Z
        //    |________|/
        //    e        f
        //
        if( (mLevel < MinLevel) || (mLevel < MaxLevel && mStatus == IntersectionStatus::Trimmed) ){
            CreateNewNode(MinLevel, MaxLevel, 0, mLowerBound, mLowerBound+delta, pOperator);                 // Corner a
            CreateNewNode(MinLevel, MaxLevel, 1, mLowerBound+delta_x, mLowerBound+delta+delta_x, pOperator); // Corner b (a+delta_x)
            CreateNewNode(MinLevel, MaxLevel, 2, mUpperBound-delta-delta_z, mUpperBound-delta_z, pOperator); // Corner c (g-delta_z)
            CreateNewNode(MinLevel, MaxLevel, 3, mLowerBound+delta_y, mLowerBound+delta+delta_y, pOperator); // Corner d (a+delta_y)
            CreateNewNode(MinLevel, MaxLevel, 4, mLowerBound+delta_z, mLowerBound+delta+delta_z, pOperator); // Corner e (a+delta_z)
            CreateNewNode(MinLevel, MaxLevel, 5, mUpperBound-delta-delta_y, mUpperBound-delta_y, pOperator); // Corner f (g-delta_y)
            CreateNewNode(MinLevel, MaxLevel, 6, mUpperBound-delta, mUpperBound, pOperator);                 // Corner g
            CreateNewNode(MinLevel, MaxLevel, 7, mUpperBound-delta-delta_x, mUpperBound-delta_x, pOperator); // Corner h (g-delta_x)
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
void Octree<TOperator>::Node::GetIntegrationPoints(IntegrationPointVectorType* pPoints, const PointType& GlobalLowerBound, const PointType& GlobalUpperBound, const Vector3i& rOrder, const TOperator* pOperator) {
    if( this->IsLeaf() ){
        const auto lower_bound_param = Mapping::GlobalToParam(mLowerBound, GlobalLowerBound, GlobalUpperBound);
        const auto upper_bound_param = Mapping::GlobalToParam(mUpperBound, GlobalLowerBound, GlobalUpperBound);

        IntegrationPointVectorType integration_points_tmp{};
        SingleElement::AssembleIPs(integration_points_tmp, lower_bound_param, upper_bound_param, rOrder);
        if( mStatus == IntersectionStatus::Inside )
            pPoints->insert(pPoints->end(), integration_points_tmp.begin(), integration_points_tmp.end());
        else {
            for( auto& point : integration_points_tmp){
                const auto tmp_point = Mapping::ParamToGlobal(point, GlobalLowerBound, GlobalUpperBound);
                if( pOperator->IsInsideTrimmedDomain( tmp_point ) ){
                    pPoints->push_back(point);
                }
            }
        }
    } else {
        for( IndexType i = 0; i < 8UL; ++i){
            if( mChildren[i] ){ // If not nullptr
                mChildren[i]->GetIntegrationPoints(pPoints, GlobalLowerBound, GlobalUpperBound, rOrder, pOperator);
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
void Octree<TOperator>::Node::CreateNewNode(IndexType MinLevel, IndexType MaxLevel, IndexType i, const PointType& rLowerBound, const PointType& rUpperBound, const TOperator* pOperator){
    const auto status = pOperator->GetIntersectionState(rLowerBound, rUpperBound);
    if( status != IntersectionStatus::Outside ){
        mChildren[i] = MakeUnique<Node>(rLowerBound, rUpperBound, status, mLevel+1);
        ++mNumChildren;
        mChildren[i]->Refine(MinLevel, MaxLevel, pOperator);
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
    TIBRA_ERROR_IF("Octree :: Constructor", MinLevel > MaxLevel ) << "MinLevel must be smaller/equal than MaxLevel. "
        << "Given MinLevel: " << MinLevel << ", MaxLevel: " << MaxLevel << ".\n";
    mpRoot->Refine(MinLevel, MaxLevel, mpOperator);
}


template<typename TOperator>
SizeType Octree<TOperator>::NumberOfLeafs() const{
    SizeType number_of_leafs = 0UL;
    mpRoot->NumberOfLeafs(number_of_leafs);
    return number_of_leafs;
}

template<typename TOperator>
typename Octree<TOperator>::IntegrationPointVectorPtrType Octree<TOperator>::pGetIntegrationPoints(const Vector3i& rOrder){
    auto p_points = MakeUnique<IntegrationPointVectorType>();
    p_points->reserve(NumberOfLeafs());
    mpRoot->GetIntegrationPoints(p_points.get(), mrParameters.LowerBound(), mrParameters.UpperBound(), rOrder, mpOperator );
    return std::move(p_points);
}

template<typename TOperator>
void Octree<TOperator>::AddIntegrationPoints(IntegrationPointVectorType& rPoints, const Vector3i& rOrder){
    rPoints.reserve(NumberOfLeafs());
    mpRoot->GetIntegrationPoints(&rPoints, mrParameters.LowerBound(), mrParameters.UpperBound(), rOrder, mpOperator );
}

// Explicit instantiation Octree
template class Octree<TrimmedDomainBase>;

} // End namespace tibra