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

    typedef std::vector<IntersectionStatusType> StatusVectorType;
    typedef std::stack<IndexType> IndexStackType;
    typedef std::vector< std::tuple<int, int, int> > GroupVectorType;
    typedef std::vector<bool> BoolVectorType;
                       // Partition index, Indices, Count
    typedef std::tuple<IndexType, std::set<IndexType>, int > GroupSetType;
    typedef std::vector<GroupSetType> GroupVectorSetType;
    typedef std::pair<PointType, PointType> BoundingBoxType;
    typedef std::pair<Vector3i, Vector3i> PartitionBoxType;
    typedef std::pair<IndexType, IndexType> Partition1DBoxType;
    typedef std::vector<std::vector<std::set<IndexType>>> BoundaryIndicesVectorType;

    ///@}
    ///@name Life cycle
    ///@{
    FloodFill(BRepOperatorBase* pBrepOperator, const Parameters& Parameters) :
        mpBrepOperator(pBrepOperator), mMapper(Parameters), mLowerBound(Parameters.LowerBound()), mUpperBound(Parameters.UpperBound()),
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

    Unique<StatusVectorType> ClassifyElements() const;

protected:
    std::pair<Unique<StatusVectorType>, Unique<GroupVectorSetType>> ClassifyElementsForTest() const;
private:

    Unique<StatusVectorType> ClassifyElements(GroupVectorSetType& rGroupsOutput) const;

    bool Touch(const Partition1DBoxType& rBox1, const Partition1DBoxType& rBox2) const  {

        if ( static_cast<int>(rBox1.second) < static_cast<int>(rBox2.first)-1
                || rBox1.first > rBox2.second+1 ) {
            return false;
        }

        return true;
    }

    int GetIsInsideCount( IndexType Index, IndexType NextIndex, const PointType& rLowerOffset, const PointType& rUpperOffset ) const{
        auto box_current = mMapper.GetBoundingBoxFromIndex(Index);
        auto box_next = mMapper.GetBoundingBoxFromIndex(NextIndex);
        const PointType center_box = (box_current.first + box_current.second)*0.5;
        if( mpBrepOperator->OnBoundedSideOfClippedSection(center_box, box_next.first + rLowerOffset , box_next.second + rUpperOffset) ) {
            return 1;
        } else {
            return -1;
        }
    }

    int GetNextIndex( IndexType Direction, IndexType Index ) const {
        PointType dummy_1, dummy_2;
        PartitionBoxType partition_box = std::make_pair(Vector3i(0, 0, 0), mNumberOfElements);
        return GetNextIndex(Direction, Index, partition_box, dummy_1, dummy_2);
    }

    int GetNextIndex( IndexType Direction, IndexType Index, PointType& rLowerBoundOffset, PointType& rUpperBoundOffset ) const {
        PartitionBoxType partition_box = std::make_pair(Vector3i(0, 0, 0), mNumberOfElements);
        return GetNextIndex(Direction, Index, partition_box, rLowerBoundOffset, rUpperBoundOffset);
    }

    int GetNextIndex(IndexType Direction, IndexType Index, const PartitionBoxType& rPartition, PointType& rLowerBoundOffset, PointType& rUpperBoundOffset) const {
        auto indices = mMapper.GetMatrixIndicesFromVectorIndex(Index);
        Vector3i next_indices = indices;
        rLowerBoundOffset = {0.0, 0.0, 0.0};
        rUpperBoundOffset = {0.0, 0.0, 0.0};
        double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);
        switch(Direction){
            case 0:
                if( indices[0] < rPartition.second[0]-1 ){ next_indices[0] += 1; }
                else { return -1; }
                rLowerBoundOffset = {-tolerance, 0.0, 0.0};
                break;
            case 1:
                if( indices[0] > rPartition.first[0] ){ next_indices[0] -= 1; }
                else { return -1; }
                rUpperBoundOffset = {tolerance, 0.0, 0.0};
                break;
            case 2:
                if( indices[1] < rPartition.second[1]-1 ){ next_indices[1] += 1; }
                else { return -1; }
                rLowerBoundOffset = {0.0, -tolerance, 0.0};
                break;
            case 3:
                if( indices[1] > rPartition.first[1] ){ next_indices[1] -= 1; }
                else { return -1; }
                rUpperBoundOffset = {0.0, tolerance, 0.0};
                break;
            case 4:
                if( indices[2] < rPartition.second[2]-1 ){ next_indices[2] += 1; }
                else { return -1; }
                rLowerBoundOffset = {0.0, 0.0, -tolerance};
                break;
            case 5:
                if( indices[2] > rPartition.first[2] ){ next_indices[2] -= 1; }
                else { return -1; }
                rUpperBoundOffset = {0.0, 0.0, tolerance};
                break;
            default:
                TIBRA_ERROR("FloodFill::FillDirection") << " Direction is out-of-range.\n";
        }

        return mMapper.GetVectorIndexFromMatrixIndices(next_indices[0], next_indices[1], next_indices[2]);
    }

    void PartitionedFill(IndexType PartitionIndex, GroupVectorSetType& rGroupVectorSet, PartitionBoxType rPartition, BoolVectorType& rVisited, StatusVectorType& rStates) const;

    void SinglePartitionFill(IndexType Index, GroupSetType& rGroupSet, PartitionBoxType rPartition, BoolVectorType& rVisited, StatusVectorType& rStates) const;

    int FillDirection(IndexType Direction, IndexStackType& rStack, GroupSetType& rGroupSet, BoolVectorType& rVisited, PartitionBoxType& rPartition, StatusVectorType& rStates ) const;

    void MergeGroups(IndexType PartitionDir, std::vector<PartitionBoxType>& rPartitions, GroupVectorSetType& rMergedGroups, GroupVectorSetType& rGroupVectorSet, StatusVectorType& rStates) const;

    void GroupFill(IndexType PartitionDir, BoundaryIndicesVectorType& rBoundaryIndics, GroupVectorSetType& rMergedGroups, IndexType GroupIndex, GroupVectorSetType& rGroupVectorSet, BoolVectorType& rVisited, StatusVectorType& rStates) const;


    const BRepOperatorBase* mpBrepOperator;
    Mapper mMapper;
    const PointType mLowerBound;
    const PointType mUpperBound;
    const Vector3i mNumberOfElements;
    PointType mDelta;

}; // End class FloodFill

///@} TIBRA Classes
} // End tibra namespace

#endif // FLOOD_FILL_INCLUDE_H