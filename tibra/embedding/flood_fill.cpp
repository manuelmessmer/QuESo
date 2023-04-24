// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// External includes
#include <omp.h>
#include <map>
//// Project includes
#include "embedding/flood_fill.h"


namespace tibra {

typedef FloodFill::StatusVectorType StatusVectorType;
typedef FloodFill::StatusVector2Type StatusVector2Type;

void FloodFill::PartitionedFill(GroupVectorSetType& rGroupVectorSet, PartitionBoxType rPartition, std::vector<bool>& rVisited, StatusVector2Type& rStates) const {

    for( IndexType i = rPartition.first[0]; i < rPartition.second[0]; ++i ){
        for( IndexType j = rPartition.first[1]; j < rPartition.second[1]; ++j ) {
            for( IndexType k = rPartition.first[2]; k < rPartition.second[2]; ++k ) {
                const IndexType index = mIdMapper.GetVectorIndexFromMatrixIndices(i, j, k);
                if( !rVisited[index] ) { // Unvisited
                    GroupSetType new_group;
                    SinglePartitionFill(index, new_group, rPartition, rVisited, rStates);
                    if( new_group.first.size() > 0 ){
                        #pragma omp critical
                        rGroupVectorSet.push_back(new_group);
                    }
                }
            }
        }
    }
}

void FloodFill::SinglePartitionFill(IndexType Index, GroupSetType& rGroupSet, PartitionBoxType rPartition, std::vector<bool>& rVisited, StatusVector2Type& rStates) const {

    // Set index as visited
    rVisited[Index] = true;
    // Get indices
    const auto indices = mIdMapper.GetMatrixIndicesFromVectorIndex(Index);
    const auto box = GetBoundingBoxFromIndex(indices);
    if( mpBrepOperator->IsTrimmed(box.first, box.second) ){
        rStates[Index] = IntersectionStatus::Trimmed;
    } else {

        rGroupSet.first.insert(Index);
        // If box is not trimmed, run flood fill.
        IndexStackType index_stack;

        index_stack.push( Index );
        //std::cout << "start flodd 2D "  << std::endl;
        std::array<int, 6> new_indices;
        while( !index_stack.empty() ){
            for( IndexType i = 0; i < 6; ++i){
                new_indices[i] = FillDirectionNew(i, index_stack, rGroupSet, rVisited, rPartition, rStates );
            }

            index_stack.pop();
            for( IndexType i = 0; i < 6; ++i){
                if( new_indices[i] > -1 ){
                    index_stack.push(new_indices[i]);
                }
            }
        }
    }

}

int FloodFill::FillDirectionNew(IndexType Direction, IndexStackType& rStack, GroupSetType& rGroupSet, std::vector<bool>& rVisited, PartitionBoxType& rPartition, StatusVector2Type& rStates ) const{
    bool local_end = false;
    const IndexType index = rStack.top();
    const Vector3i indices = mIdMapper.GetMatrixIndicesFromVectorIndex(index);
    const IndexType max_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];
    Vector3i next_indices = indices;

    PointType lower_perturb = {0.0, 0.0, 0.0};
    PointType upper_perturb = {0.0, 0.0, 0.0};
    double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);

    switch(Direction){
        case 0:
            if( indices[0] < rPartition.second[0]-1 ){ next_indices[0] += 1; }
            else { return -1; }
            lower_perturb = {-tolerance, 0.0, 0.0};
            break;
        case 1:
            if( indices[0] > rPartition.first[0] ){ next_indices[0] -= 1; }
            else { return -1; }
            upper_perturb = {tolerance, 0.0, 0.0};
            break;
        case 2:
            if( indices[1] < rPartition.second[1]-1 ){ next_indices[1] += 1; }
            else { return -1; }
            lower_perturb = {0.0, -tolerance, 0.0};
            break;
        case 3:
            if( indices[1] > rPartition.first[1] ){ next_indices[1] -= 1; }
            else { return -1; }
            upper_perturb = {0.0, tolerance, 0.0};
            break;
        case 4:
            if( indices[2] < rPartition.second[2]-1 ){ next_indices[2] += 1; }
            else { return -1; }
            lower_perturb = {0.0, 0.0, -tolerance};
            break;
        case 5:
            if( indices[2] > rPartition.first[2] ){ next_indices[2] -= 1; }
            else { return -1; }
            upper_perturb = {0.0, 0.0, tolerance};
            break;
        default:
            TIBRA_ERROR("error") << "error\n";
    }


    const IndexType next_index = mIdMapper.GetVectorIndexFromMatrixIndices(next_indices[0], next_indices[1], next_indices[2]);

    //Check if out-of-range
    if( next_index >= max_num_elements || next_index < 0 ) {
        return -1;
    }

    auto box_next = GetBoundingBoxFromIndex(next_indices);

    if( mpBrepOperator->IsTrimmed(box_next.first + lower_perturb , box_next.second + upper_perturb) ){
        auto box_current = GetBoundingBoxFromIndex(indices);
        PointType center_box = (box_current.first + box_current.second)*0.5;
        if( mpBrepOperator->OnBoundedSideOfClippedSection(center_box, box_next.first + lower_perturb , box_next.second + upper_perturb) ) {
            rGroupSet.second += 1;
        } else {
            rGroupSet.second -= 1;
        }
        return -1;
    }
    // Already classified
    if( rVisited[next_index] ) {
        return -1;
    }

    // Is trimmed
    if( mpBrepOperator->IsTrimmed(box_next.first, box_next.second) ){
        rVisited[next_index] = true;
        rStates[next_index] = IntersectionStatus::Trimmed;

        return -1;
    }

    rVisited[next_index] = true;
    rGroupSet.first.insert(next_index);

    return next_index;
}

Unique<StatusVector2Type> FloodFill::ClassifyElements() const {

        Timer timer{};
        Parameters params( {Component("lower_bound", mLowerBound),
                            Component("upper_bound", mUpperBound),
                            Component("number_of_elements", mNumberOfElements) });
        ElementContainer elcont(params);

        IndexType max_num_elements = Math::Max(mNumberOfElements);
        StatusVector2Type states(mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2]);
        std::fill(states.begin(), states.end(), IntersectionStatus::NotVisited );


        IndexType total_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];


        GroupVectorSetType groups;
        PartitionBoxType partition_box = std::make_pair(Vector3i(0, 0, 0), mNumberOfElements );
        std::vector<bool> visited(total_num_elements, false);
        PartitionedFill(groups, partition_box, visited, states);

        int id = 0;
        int cc = 0;
        for( auto& r_group : groups ){
            int inside_count = r_group.second;
            auto state = (inside_count > 0) ? IntersectionStatus::Inside : IntersectionStatus::Outside;
            //std::cout << "Size: " << r_group.first.size() << ", Count: " << inside_count << std::endl;
            for( auto group_it = r_group.first.begin(); group_it !=  r_group.first.end(); ++group_it){
                // auto box = GetBoundingBoxFromIndex(*group_it);
                // auto lower_bound_param = Mapping::GlobalToParam(box.first, mLowerBound,mUpperBound);
                // auto upper_bound_param = Mapping::GlobalToParam(box.second, mLowerBound, mUpperBound);
                // Shared<Element> new_el = MakeShared<Element>(id,  lower_bound_param,upper_bound_param, params );
                // if( cc > 0  ){
                //     elcont.AddElement(new_el);
                //     id++;
                // }
                states[*group_it] = state;
            }
            // cc++;
            // IO::WriteElementsToVTK(elcont, "Elements.vtk", true);
            //TIBRA_ERROR("test") << "waf\n";

        }

        return MakeUnique<StatusVector2Type>(states);
    }

    int FloodFill::Fill(IndexType index, GroupVectorType& rGroups, IndexType& rGroupId) const {

        const auto indices = mIdMapper.GetMatrixIndicesFromVectorIndex(index);

        std::get<0>(rGroups[index]) = indices[0];
        std::get<1>(rGroups[index]) = rGroupId;
        //std::cout << "index " << index << std::endl;
        //rStates[index].second = IntersectionStatus::PartOfNewGroup;
        const auto box = GetBoundingBoxFromIndex(indices);
        if( mpBrepOperator->IsTrimmed(box.first, box.second) ){
            std::get<0>(rGroups[index]) = indices[0];
            std::get<1>(rGroups[index]) = -1;
        } else {
            // If box is not trimmed, run flood fill.
            IndexStackType index_stack;

            index_stack.push( index );
            //std::cout << "start flodd 2D "  << std::endl;
            std::array<int, 6> new_indices;
            while( !index_stack.empty() ){
                for( IndexType i = 2; i < 6; ++i){
                    new_indices[i] = FillDirection(i, index_stack, rGroups, rGroupId );
                }

                index_stack.pop();
                for( IndexType i = 2; i < 6; ++i){
                    if( new_indices[i] > -1 ){
                        index_stack.push(new_indices[i]);
                    }
                }
            }
            ++rGroupId;

        }

        return 1;
    }


    int FloodFill::FillDirection(IndexType Direction, IndexStackType& rStack, GroupVectorType& rGroups, IndexType& rGroupId ) const{
        bool local_end = false;
        const IndexType index = rStack.top();
        const Vector3i indices = mIdMapper.GetMatrixIndicesFromVectorIndex(index);
        const IndexType max_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];
        Vector3i next_indices = indices;

        PointType lower_perturb = {0.0, 0.0, 0.0};
        PointType upper_perturb = {0.0, 0.0, 0.0};
        double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);

        switch(Direction){
            case 0:
                if( indices[0] < mNumberOfElements[0]-1 ){ next_indices[0] += 1; }
                else { return -1; }
                lower_perturb = {-tolerance, 0.0, 0.0};
                break;
            case 1:
                if( indices[0] > 0 ){ next_indices[0] -= 1; }
                else { return -1; }
                upper_perturb = {tolerance, 0.0, 0.0};
                break;
            case 2:
                if( indices[1] < mNumberOfElements[1]-1 ){ next_indices[1] += 1; }
                else { return -1; }
                lower_perturb = {0.0, -tolerance, 0.0};
                break;
            case 3:
                if( indices[1] > 0 ){ next_indices[1] -= 1; }
                else { return -1; }
                upper_perturb = {0.0, tolerance, 0.0};
                break;
            case 4:
                if( indices[2] < mNumberOfElements[2]-1 ){ next_indices[2] += 1; }
                else { return -1; }
                lower_perturb = {0.0, 0.0, -tolerance};
                break;
            case 5:
                if( indices[2] > 0 ){ next_indices[2] -= 1; }
                else { return -1; }
                upper_perturb = {0.0, 0.0, tolerance};
                break;
            default:
                TIBRA_ERROR("error") << "error\n";
        }


        const IndexType next_index = mIdMapper.GetVectorIndexFromMatrixIndices(next_indices[0], next_indices[1], next_indices[2]);

        //Check if out-of-range
        if( next_index >= max_num_elements || next_index < 0 ) {
            return -1;
        }

        auto box_next = GetBoundingBoxFromIndex(next_indices);

        if( mpBrepOperator->IsTrimmed(box_next.first + lower_perturb , box_next.second + upper_perturb) ){
            auto box_current = GetBoundingBoxFromIndex(indices);
            PointType center_box = (box_current.first + box_current.second)*0.5;
            if( mpBrepOperator->OnBoundedSideOfClippedSection(center_box, box_next.first + lower_perturb , box_next.second + upper_perturb) ) {
                std::get<2>(rGroups[index]) += 1;
            } else {
                std::get<2>(rGroups[index]) -= 1;
            }
            return -1;
        }
        // Already classified
        if( std::get<0>(rGroups[next_index]) > -1 ) {
            return -1;
        }

        // Is trimmed
        if( mpBrepOperator->IsTrimmed(box_next.first, box_next.second) ){
            std::get<0>(rGroups[next_index]) = indices[0];
            std::get<1>(rGroups[next_index]) = -1;

            return -1;
        }

        std::get<0>(rGroups[next_index]) = indices[0];
        std::get<1>(rGroups[next_index]) = rGroupId;
        return next_index;
    }

    // int FloodFill::FillDirection(IndexType Direction, IndexStackType& stack, StatusVectorType& rStates, int& rInsideCount) const{
    //     bool local_end = false;
    //     IndexType index = stack.top();

    //     int next_index = -1;
    //     PointType lower_perturb = {0.0, 0.0, 0.0};
    //     PointType upper_perturb = {0.0, 0.0, 0.0};
    //     double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);
    //     switch(Direction){
    //         case 0:
    //             next_index = mIdMapper.GetNextIndexX(index, local_end);
    //             lower_perturb = {-tolerance, 0.0, 0.0};
    //             break;
    //         case 1:
    //             next_index = mIdMapper.GetPreviousIndexX(index, local_end);
    //             upper_perturb = {tolerance, 0.0, 0.0};
    //             break;
    //         case 2:
    //             next_index = mIdMapper.GetNextIndexY(index, local_end);
    //             lower_perturb = {0.0, -tolerance, 0.0};
    //             break;
    //         case 3:
    //             next_index = mIdMapper.GetPreviousIndexY(index, local_end);
    //             upper_perturb = {0.0, tolerance, 0.0};
    //             break;
    //         case 4:
    //             next_index = mIdMapper.GetNextIndexZ(index, local_end);
    //             lower_perturb = {0.0, 0.0, -tolerance};
    //             break;
    //         case 5:
    //             next_index = mIdMapper.GetPreviousIndexZ(index, local_end);
    //             upper_perturb = {0.0, 0.0, tolerance};
    //             break;
    //         default:
    //             TIBRA_ERROR("error") << "error\n";
    //     }

    //     // Check if out-of-range
    //     if( next_index >= rStates.size() || next_index < 0 ) {
    //         return -1;
    //     }

    //     // Already classified
    //     if( rStates[next_index] != IntersectionStatus::NotVisited ) {
    //         return -1;
    //     }


    //     auto indices =  mIdMapper.GetMatrixIndicesFromVectorIndex(next_index);
    //     PointType indices_d( indices[0], indices[1], indices[2] );
    //     PointType lower_bound = mLowerBound + indices_d * mDelta;
    //     PointType upper_bound = mLowerBound + (indices_d+1.0) * mDelta;

    //     if( mpBrepOperator->IsTrimmed(lower_bound + lower_perturb , upper_bound + upper_perturb) ){
    //         auto indices_tmp =  mIdMapper.GetMatrixIndicesFromVectorIndex(index);
    //         PointType indices_d_tmp( indices_tmp[0], indices_tmp[1], indices_tmp[2] );
    //         PointType lower_bound_tmp = mLowerBound + indices_d_tmp * mDelta;
    //         PointType upper_bound_tmp = mLowerBound + (indices_d_tmp+1.0) * mDelta;
    //         PointType center_box = (upper_bound_tmp + lower_bound_tmp)*0.5;

    //         if( mpBrepOperator->OnBoundedSideOfClippedSection(center_box, lower_bound + lower_perturb , upper_bound + upper_perturb) ) {
    //             rInsideCount += 1;
    //         } else {
    //             rInsideCount -= 1;
    //         }
    //         return -1;
    //     }

    //     if( mpBrepOperator->IsTrimmed(lower_bound, upper_bound) ){
    //         rStates[next_index] = IntersectionStatus::Trimmed;
    //         return -1;
    //     }

    //     if( mpBrepOperator->IsTrimmed(lower_bound + lower_perturb , upper_bound + upper_perturb) ){
    //         rStates[next_index] = IntersectionStatus::OppositeGroup;
    //         return -1;
    //     }

    //     rStates[next_index] = IntersectionStatus::PartOfNewGroup;
    //     return next_index;
    // }

    std::pair<PointType, PointType> FloodFill::GetBoundingBoxFromIndex(IndexType Index) const{
        const auto indices =  mIdMapper.GetMatrixIndicesFromVectorIndex(Index);
        const PointType indices_d( indices[0], indices[1], indices[2] );
        return std::make_pair( mLowerBound + indices_d * mDelta,
                               mLowerBound + (indices_d+1.0) * mDelta );
    }


} // End namespace tibra