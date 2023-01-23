// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef OCTREE_INCLUDE_H
#define OCTREE_INCLUDE_H

//// STL includes
#include <array>

//// Project includes
#include "embedding/trimmed_domain_base.h"
#include "define.hpp"

namespace tibra {

///@name TIBRA Classes
///@{

class Node {
    public:
        Node( const PointType& rLowerBound, const PointType& rUpperBound, IndexType Level = 0UL) :
            mLowerBound(rLowerBound), mUpperBound(rUpperBound), mLevel(Level)
        {
            mNumChildren = 0;
        }

        void Refine(IndexType MinLevel, IndexType MaxLevel){
            if( this->IsLeaf() ){
                if( mLevel < MinLevel ){
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
                    mChildren[0] = MakeUnique<Node>(mLowerBound, mLowerBound+delta);                 // Corner a
                    mChildren[1] = MakeUnique<Node>(mLowerBound+delta_x, mLowerBound+delta+delta_x); // Corner b (a+delta_x)
                    mChildren[2] = MakeUnique<Node>(mUpperBound-delta-delta_z, mUpperBound-delta_z); // Corner c (g-delta_z)
                    mChildren[3] = MakeUnique<Node>(mLowerBound+delta_y, mLowerBound+delta+delta_y); // Corner d (a+delta_y)
                    mChildren[4] = MakeUnique<Node>(mLowerBound+delta_z, mLowerBound+delta+delta_z); // Corner e (a+delta_z)
                    mChildren[5] = MakeUnique<Node>(mUpperBound-delta-delta_z, mUpperBound-delta_z); // Corner f (g-delta_z)
                    mChildren[7] = MakeUnique<Node>(mUpperBound-delta, mUpperBound);                 // Corner g
                    mChildren[6] = MakeUnique<Node>(mUpperBound-delta-delta_x, mUpperBound-delta_x); // Corner h (g-delta_x)
                }
            }
            else {
                for( IndexType i = 0; i < 8; ++i){
                    if( mChildren[i] ){
                        mChildren[i]->Refine(MinLevel, MaxLevel);
                    }
                }
            }
        }

        void Volume(double& volume){
            if( IsLeaf() ){
                const auto delta = (mUpperBound - mLowerBound);
                volume += delta[0]*delta[1]*delta[2];
            } else {
                for( IndexType i = 0; i < 8; ++i){
                    if( mChildren[i] ){
                        mChildren[i]->Volume(volume);
                    }
                }
            }
        }

        bool IsLeaf(){
            return (mNumChildren == 0);
        }

        IndexType Level(){
            return mLevel;
        }

    private:
        PointType mLowerBound;
        PointType mUpperBound;
        std::array<Unique<Node>, 8> mChildren{};
        SizeType mLevel;
        SizeType mNumChildren;
};

/**
 * @class  Octree
 * @author Manuel Messmer
*/
class Octree {

public:
    ///@name Type Definitions
    ///@{

    Octree(TrimmedDomainBase* pTrimmedDomain) : mpTrimmedDomain(pTrimmedDomain)
    {
        const auto bounding_box = pTrimmedDomain->GetBoundingBoxOfTrimmedDomain();
        mpRoot = MakeUnique<Node>(bounding_box.first, bounding_box.second);
    }

    void Refine(IndexType MinLevel, IndexType MaxLevel){
        mpRoot->Refine(MinLevel, MaxLevel);
    }

    double Volume(){
        double volume = 0.0;
        mpRoot->Volume(volume);
        return volume;
    }

private:


    TrimmedDomainBase* mpTrimmedDomain;
    Unique<Node> mpRoot;

protected:

};

} // End tibr namespace

#endif // OCTREE_INCLUDE_H