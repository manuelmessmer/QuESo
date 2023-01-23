// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef OCTREE_INCLUDE_H
#define OCTREE_INCLUDE_H

//// STL includes
#include <array>

//// Project includes
#include "embedding/trimmed_domain_base.h"
#include "quadrature/single_element.h"
#include "define.hpp"

namespace tibra {

///@name TIBRA Classes
///@{



/**
 * @class  Octree
 * @author Manuel Messmer
 * @tparam TOperator: Required member operations:
 *                   -GetIntersectionState(const PointType& rLowerBound, const PointType& rUpperBound)
 *                   -IsInside(const PointType& rPoint)
*/
template<typename TOperator>
class Octree {
private:
    typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
    typedef Unique<IntegrationPointVectorType> IntegrationPointVectorPtrType;

    class Node {
        public:
            Node( const PointType& rLowerBound, const PointType& rUpperBound, IntersectionStatusType Status, IndexType Level = 0UL) :
                mLowerBound(rLowerBound), mUpperBound(rUpperBound), mStatus(Status), mLevel(Level)
            {
                mChildren = {nullptr};
                mNumChildren = 0;
            }

            void Refine(IndexType MinLevel, IndexType MaxLevel, TOperator* pOperator){
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
                    for( IndexType i = 0; i < 8; ++i){
                        if( mChildren[i] ){
                            mChildren[i]->Refine(MinLevel, MaxLevel, pOperator);
                        }
                    }
                }

            }

            void GetIntegrationPoints(IntegrationPointVectorType* pPoints, const PointType& GlobalLowerBound, const PointType& GlobalUpperBound, TOperator* pOperator) {
                if( this->IsLeaf() ){
                    auto lower_bound_param = Mapping::GlobalToParam(mLowerBound, GlobalLowerBound, GlobalUpperBound);
                    auto upper_bound_param = Mapping::GlobalToParam(mUpperBound, GlobalLowerBound, GlobalUpperBound);

                    IntegrationPointVectorType r_integration_points{};
                    SingleElement::AssembleIPs(r_integration_points, lower_bound_param, upper_bound_param, {2Ul, 2Ul, 2Ul});
                    if( mStatus == IntersectionStatus::Inside )
                        pPoints->insert(pPoints->end(), r_integration_points.begin(), r_integration_points.end());
                    else {
                        for( auto& point : r_integration_points){
                            auto tmp_point = Mapping::ParamToGlobal(point, GlobalLowerBound, GlobalUpperBound);
                            if( pOperator->IsInsideTrimmedDomain( tmp_point ) ){
                                pPoints->push_back(point);
                            }
                        }
                    }
                } else {
                    for( IndexType i = 0; i < 8; ++i){
                        if( mChildren[i] ){
                            mChildren[i]->GetIntegrationPoints(pPoints, GlobalLowerBound, GlobalUpperBound, pOperator);
                        }
                    }
                }
            }

            void Count(IndexType& rValue) {
                if( this->IsLeaf() ){
                    ++rValue;
                } else {
                    for( IndexType i = 0; i < 8; ++i){
                        if( mChildren[i] ){
                            mChildren[i]->Count(rValue);
                        }
                    }
                }
            }

            void Volume(double& volume){
                if( IsLeaf() ){
                    const auto delta = (mUpperBound - mLowerBound);
                    volume += delta[0]*delta[1]*delta[2];
                    TIBRA_INFO << delta[0]*delta[1]*delta[2] << std::endl;
                } else {
                    for( IndexType i = 0; i < 8; ++i){
                        if( mChildren[i] ){
                            mChildren[i]->Volume(volume);
                        }
                    }
                }
            }

        private:
            void CreateNewNode(IndexType MinLevel, IndexType MaxLevel, IndexType i, const PointType& rLowerBound, const PointType& rUpperBound, TOperator* pOperator){
                const auto status = pOperator->GetIntersectionState(rLowerBound, rUpperBound, 1e-14);
                if( status != IntersectionStatus::Outside ){
                    mChildren[i] = MakeUnique<Node>(rLowerBound, rUpperBound, status, mLevel+1);
                    ++mNumChildren;
                    mChildren[i]->Refine(MinLevel, MaxLevel, pOperator);
                }
            }

            bool IsLeaf(){
                return (mNumChildren == 0);
            }

            PointType mLowerBound;
            PointType mUpperBound;
            IntersectionStatus mStatus;
            std::array<Unique<Node>, 8> mChildren{};
            SizeType mLevel;
            SizeType mNumChildren;
    };

public:
    ///@name Type Definitions
    ///@{

    Octree(TOperator* pOperator, const PointType& rLowerBound, const PointType& rUpperBound, const Parameters& rParameters)
        : mpOperator(pOperator), mrParameters(rParameters)
    {
        mpRoot = MakeUnique<Node>(rLowerBound, rUpperBound, IntersectionStatus::Trimmed, 0UL);
    }

    void Refine(IndexType MinLevel, IndexType MaxLevel){
        TIBRA_ERROR_IF("Octree :: Constructor", MinLevel > MaxLevel ) << "MinLevel must be smaller/equal than MaxLevel. "
            << "Given MinLevel: " << MinLevel << ", MaxLevel: " << MaxLevel << ".\n";
        mpRoot->Refine(MinLevel, MaxLevel, mpOperator);
    }

    SizeType NumberOfLeafs(){
        SizeType number_of_leafs = 0UL;
        mpRoot->Count(number_of_leafs);
        return number_of_leafs;

    }

    IntegrationPointVectorPtrType GetIntegrationPoints(){
        auto p_points = MakeUnique<IntegrationPointVectorType>();
        p_points->reserve(NumberOfLeafs());
        mpRoot->GetIntegrationPoints(p_points.get(), mrParameters.LowerBound(), mrParameters.UpperBound(), mpOperator );
        return std::move(p_points);
    }

    double Volume(){
        double volume = 0.0;
        mpRoot->Volume(volume);
        return volume;
    }

private:

    Unique<Node> mpRoot;
    const Parameters& mrParameters;
    TOperator* mpOperator;

protected:

};

} // End tibr namespace

#endif // OCTREE_INCLUDE_H