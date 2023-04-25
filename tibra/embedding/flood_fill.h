// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef FLOOD_FILL_INCLUDE_H
#define FLOOD_FILL_INCLUDE_H

//// STL includes
#include <array>
#include <algorithm>
#include <stack>
#include <unordered_set>

//// Project includes
#include "define.hpp"
#include "utilities/timer.hpp"
#include "utilities/vector_matrix_id_utilities.h"
#include "utilities/mesh_utilities.h"
#include "containers/element_container.hpp"
#include "embedding/brep_operator.h"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  FloodFill
 * @author Manuel Messmer
*/
class FloodFill {
public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<std::pair<IndexType, IntersectionStatusType>> StatusVectorType;
    typedef std::vector<IntersectionStatusType> StatusVector2Type;
    typedef std::stack<IndexType> IndexStackType;
    typedef std::vector< std::tuple<int, int, int> > GroupVectorType;
                       // Partition index, Indices, Count
    typedef std::tuple<IndexType, std::set<IndexType>, int > GroupSetType;
    typedef std::vector<GroupSetType> GroupVectorSetType;
    typedef std::pair<PointType, PointType> BoundingBoxType;
    typedef std::pair<Vector3i, Vector3i> PartitionBoxType;

    ///@}
    ///@name Life cycle
    ///@{
    FloodFill(BRepOperatorBase* pBrepOperator, const Parameters& Parameters) :
        mpBrepOperator(pBrepOperator), mIdMapper(Parameters.NumberOfElements()), mLowerBound(Parameters.LowerBound()), mUpperBound(Parameters.UpperBound()),
        mNumberOfElements( Parameters.NumberOfElements() )
    {
        // Obtain discretization of background mesh.
        mDelta[0] = std::abs(mUpperBound[0] - mLowerBound[0]) / (mNumberOfElements[0]);
        mDelta[1] = std::abs(mUpperBound[1] - mLowerBound[1]) / (mNumberOfElements[1]);
        mDelta[2] = std::abs(mUpperBound[2] - mLowerBound[2]) / (mNumberOfElements[2]);
    }

    ///@}
    ///@name Operations
    ///@{


// bool AABB_primitive::intersect(const AABB_primitive &aabb) const  {
//     for (unsigned int i = 0; i < 3; ++i) {
//         if (aabb.upperBound[i] < lowerBound[i] || aabb.lowerBound[i] > upperBound[i] ) {
//             return false;
//         }
//     }
//     return true;
// }
    bool Touch(const PartitionBoxType& rBox1, const PartitionBoxType& rBox2) const  {
    for (unsigned int i = 0; i < 3; ++i) {
        if (rBox1.second[i] < rBox2.first[i]-1 || rBox1.first[i] > rBox2.second[i]+1 ) {
            return false;
        }
    }
    return true;
}

    int GetNextIndex(IndexType Direction, IndexType Index) const {
        auto indices = mIdMapper.GetMatrixIndicesFromVectorIndex(Index);
        Vector3i next_indices = indices;
        switch(Direction){
            case 0:
                if( indices[0] < mNumberOfElements[0]-1 ){ next_indices[0] += 1; }
                else { return -1; }
                break;
            case 1:
                if( indices[0] > 0 ){ next_indices[0] -= 1; }
                else { return -1; }
                break;
            case 2:
                if( indices[1] < mNumberOfElements[1]-1 ){ next_indices[1] += 1; }
                else { return -1; }
                break;
            case 3:
                if( indices[1] > 0 ){ next_indices[1] -= 1; }
                else { return -1; }
                break;
            case 4:
                if( indices[2] < mNumberOfElements[2]-1 ){ next_indices[2] += 1; }
                else { return -1; }
                break;
            case 5:
                if( indices[2] > 0 ){ next_indices[2] -= 1; }
                else { return -1; }
                break;
            default:
                TIBRA_ERROR("error") << "error\n";
        }

        return mIdMapper.GetVectorIndexFromMatrixIndices(next_indices[0], next_indices[1], next_indices[2]);
    }

    void PartitionedFill(IndexType PartitionIndex, GroupVectorSetType& rGroupVectorSet, PartitionBoxType rPartition, std::vector<bool>& rVisited, StatusVector2Type& rStates) const;

    void SinglePartitionFill(IndexType Index, GroupSetType& rGroupSet, PartitionBoxType rPartition, std::vector<bool>& rVisited, StatusVector2Type& rStates) const;

    int FillDirectionNew(IndexType Direction, IndexStackType& rStack, GroupSetType& rGroupSet, std::vector<bool>& rVisited, PartitionBoxType& rPartition, StatusVector2Type& rStates ) const;

    void MergeGroups(GroupVectorSetType& rMergedGroups, GroupVectorSetType& rGroupVectorSet, StatusVector2Type& rStates) const;

    void GroupFill(std::vector<std::vector<std::set<IndexType>>>& rBoundaryIndics, GroupVectorSetType& rMergedGroups, IndexType GroupIndex, std::vector<PartitionBoxType>& rPartitionBox, GroupVectorSetType& rGroupVectorSet, std::vector<bool>& rVisited, StatusVector2Type& rStates) const;

    Unique<StatusVector2Type> ClassifyElements() const;

    int Fill(IndexType index, GroupVectorType& rGroups, IndexType& rGroupId) const;

private:

    int FillDirection(IndexType Direction, IndexStackType& rStack, GroupVectorType& rGroups, IndexType& rGroupId) const;

    std::pair<PointType, PointType> GetBoundingBoxFromIndex(IndexType Index) const;

    std::pair<PointType, PointType> GetBoundingBoxFromIndex(Vector3i Indices) const{
        const PointType indices_d( Indices[0], Indices[1], Indices[2] );
        return std::make_pair( mLowerBound + indices_d * mDelta,
                               mLowerBound + (indices_d+1.0) * mDelta );
    }

    std::pair<PointType, PointType> GetBoundingBoxFromIndex(IndexType i, IndexType j, IndexType k) const{
        const PointType indices_d( i, j, k);
        return std::make_pair( mLowerBound + indices_d * mDelta,
                               mLowerBound + (indices_d+1.0) * mDelta );
    }

private:
    const BRepOperatorBase* mpBrepOperator;
    VectorMatrixIdUtilities mIdMapper;
    const PointType mLowerBound;
    const PointType mUpperBound;
    const Vector3i mNumberOfElements;
    PointType mDelta;

}; // End class FloodFill

///@} TIBRA Classes
} // End tibra namespace

#endif // FLOOD_FILL_INCLUDE_H