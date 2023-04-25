// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// External includes
#include <omp.h>
#include <map>
#include <cstdlib>
//// Project includes
#include "embedding/flood_fill.h"


namespace tibra {

typedef FloodFill::StatusVectorType StatusVectorType;
typedef FloodFill::StatusVector2Type StatusVector2Type;

Unique<StatusVector2Type> FloodFill::ClassifyElements() const {

    Timer timer{};

    StatusVector2Type states(mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2]);
    std::fill(states.begin(), states.end(), IntersectionStatus::Outside );


    IndexType total_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];


    GroupVectorSetType groups;
    PartitionBoxType partition_box = std::make_pair(Vector3i(0, 0, 0), mNumberOfElements );

    std::vector<PartitionBoxType> partition;

    IndexType num_threads = 1;
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }

    num_threads = 11;
    IndexType size = std::max<IndexType>( std::ceil( static_cast<double>(mNumberOfElements[0]) / static_cast<double>(num_threads) ), 1);
    for( IndexType i = 0; i < num_threads; ++i ){
        //std::cout << i*5 << ", " << 5*(i+1)+2 << std::endl;
        IndexType index_1 = size*i;
        if( index_1 < mNumberOfElements[0] ){
            IndexType index_2 = std::min<IndexType>(size*(i+1), mNumberOfElements[0]);
            partition.push_back( std::make_pair(Vector3i(index_1, 0, 0), Vector3i(index_2, mNumberOfElements[1], mNumberOfElements[2]) ) );
        }
    }

    std::vector<bool> visited(total_num_elements, false);
    #pragma omp parallel for firstprivate(visited)
    for( IndexType i = 0; i < partition.size(); ++i){
        PartitionedFill(i, groups, partition[i], visited, states);
    }

    GroupVectorSetType merged_groups;
    MergeGroups( merged_groups, groups, states );


    for( auto& r_group : merged_groups ){
        int inside_count = std::get<2>(r_group);
        auto state = (inside_count > 0) ? IntersectionStatus::Inside : IntersectionStatus::Outside;
        for( auto group_it = std::get<1>(r_group).begin(); group_it !=  std::get<1>(r_group).end(); ++group_it){
            states[*group_it] = state;
        }
    }


    timer.Reset();
    #pragma omp parallel for
    for( IndexType i = 0; i < total_num_elements; ++i){
        auto indices = mIdMapper.GetMatrixIndicesFromVectorIndex(i);
        auto box = GetBoundingBoxFromIndex(i);
        mpBrepOperator->GetIntersectionState(box.first, box.second);
    }
    //std::cout << "Total Classical: " << timer.Measure() << std::endl;

    return MakeUnique<StatusVector2Type>(states);
}


void FloodFill::PartitionedFill(IndexType PartitionIndex, GroupVectorSetType& rGroupVectorSet, PartitionBoxType rPartition, std::vector<bool>& rVisited, StatusVector2Type& rStates) const {

    for( IndexType i = rPartition.first[0]; i < rPartition.second[0]; ++i ){
        for( IndexType j = rPartition.first[1]; j < rPartition.second[1]; ++j ) {
            for( IndexType k = rPartition.first[2]; k < rPartition.second[2]; ++k ) {
                const IndexType index = mIdMapper.GetVectorIndexFromMatrixIndices(i, j, k);
                if( !rVisited[index] ) { // Unvisited
                    GroupSetType new_group;
                    std::get<0>(new_group) = PartitionIndex;
                    SinglePartitionFill(index, new_group, rPartition, rVisited, rStates);
                    if( std::get<1>(new_group).size() > 0 ){
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

        std::get<1>(rGroupSet).insert(Index);
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
            std::get<2>(rGroupSet) += 1;
        } else {
            std::get<2>(rGroupSet) -= 1;
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
    std::get<1>(rGroupSet).insert(next_index);

    return next_index;
}




    void FloodFill::MergeGroups(GroupVectorSetType& rMergedGroups, GroupVectorSetType& rGroupVectorSet, StatusVector2Type& rStates) const {

        const IndexType num_groups = rGroupVectorSet.size();
        std::vector<bool> visited(num_groups, false);

        double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);

        // Compute bounding box for each group
        std::vector<PartitionBoxType> group_bounding_boxes(num_groups);
        std::fill(group_bounding_boxes.begin(), group_bounding_boxes.end(),
            std::make_pair( mNumberOfElements,  Vector3i(0, 0, 0)) );
        std::vector<std::vector<std::set<IndexType>>> group_boundary_indices(num_groups);
        for( IndexType i = 0; i < num_groups; ++i){
            group_boundary_indices[i].resize(6);
        }

        #pragma omp parallel for
        for( IndexType group_index = 0; group_index < num_groups;  ++group_index){
            const auto& index_set = std::get<1>(rGroupVectorSet[group_index]);
            auto& group_bounding_box = group_bounding_boxes[group_index];
            auto& boundary_indices = group_boundary_indices[group_index];

            for(auto& index : index_set ){
                const auto indices = mIdMapper.GetMatrixIndicesFromVectorIndex(index);

                if( indices[0] >= group_bounding_box.second[0] ){
                    auto box = GetBoundingBoxFromIndex(indices[0]+1, indices[1], indices[2]);
                    if( !mpBrepOperator->IsTrimmed(box.first - PointType(tolerance, 0.0, 0.0), box.second) ) {
                        if( indices[0] > group_bounding_box.second[0] ){
                            group_bounding_box.second[0] = indices[0];
                            boundary_indices[0].clear();
                        }
                        boundary_indices[0].insert(index);
                    }
                }
                if( indices[0] <= group_bounding_box.first[0] ){
                    auto box = GetBoundingBoxFromIndex(indices[0]-1, indices[1], indices[2]);
                    if( !mpBrepOperator->IsTrimmed(box.first, box.second + PointType(tolerance, 0.0, 0.0) ) ){
                        if( indices[0] < group_bounding_box.first[0] ){
                            group_bounding_box.first[0] = indices[0];
                            boundary_indices[1].clear();
                        }
                        boundary_indices[1].insert(index);
                    }
                }
                if( indices[1] >= group_bounding_box.second[1] ){
                    auto box = GetBoundingBoxFromIndex(indices[0], indices[1]+1, indices[2]);
                    if( !mpBrepOperator->IsTrimmed(box.first - PointType(0.0, tolerance, 0.0), box.second) ) {
                        if( indices[1] > group_bounding_box.second[1] ){
                            group_bounding_box.second[1] = indices[1];
                            boundary_indices[2].clear();
                        }
                        boundary_indices[2].insert(index);
                    }
                }

                if( indices[1] <= group_bounding_box.first[1] ){
                    auto box = GetBoundingBoxFromIndex(indices[0], indices[1]-1, indices[2]);
                    if( !mpBrepOperator->IsTrimmed(box.first, box.second + PointType(0.0, tolerance, 0.0)) ){
                        if( indices[1] < group_bounding_box.first[1] ){
                            group_bounding_box.first[1] = indices[1];
                            boundary_indices[3].clear();
                        }
                        boundary_indices[3].insert(index);
                    }
                }

                if( indices[2] >= group_bounding_box.second[2] ){
                    auto box = GetBoundingBoxFromIndex(indices[0], indices[1], indices[2]+1);
                    if( !mpBrepOperator->IsTrimmed(box.first - PointType(0.0, 0.0, tolerance), box.second) ) {
                        if( indices[2] > group_bounding_box.second[2] ){
                            group_bounding_box.second[2] = indices[2];
                            boundary_indices[4].clear();
                        }
                        boundary_indices[4].insert(index);
                    }
                }

                if( indices[2] <= group_bounding_box.first[2] ){
                    auto box = GetBoundingBoxFromIndex(indices[0], indices[1], indices[2]-1);
                    if( !mpBrepOperator->IsTrimmed(box.first, box.second  + PointType(0.0, 0.0, tolerance)) ) {
                        if( indices[2] < group_bounding_box.first[2] ){
                            group_bounding_box.first[2] = indices[2];
                            boundary_indices[5].clear();
                        }
                        boundary_indices[5].insert(index);
                    }
                }



            }
        }

        //std::cout << "GroupFill: " << std::endl;

        // Test if nothing is found
        for( IndexType group_index = 0; group_index < num_groups; ++group_index){
            if( !visited[group_index] ){
                GroupFill(group_boundary_indices, rMergedGroups, group_index, group_bounding_boxes, rGroupVectorSet, visited, rStates);
            }
        }

    }

    void FloodFill::GroupFill(std::vector<std::vector<std::set<IndexType>>>& rBoundaryIndices, GroupVectorSetType& rMergedGroups, IndexType GroupIndex, std::vector<PartitionBoxType>& rPartitionBox, GroupVectorSetType& rGroupVectorSet, std::vector<bool>& rVisited, StatusVector2Type& rStates) const {

        // Set index as visited
        rVisited[GroupIndex] = true;

        // If box is not trimmed, run flood fill.
        IndexStackType index_stack;
        index_stack.push( GroupIndex );

        rMergedGroups.push_back( rGroupVectorSet[GroupIndex] );
        while( !index_stack.empty() ){
            IndexType current_group_index = index_stack.top();
            index_stack.pop();
            IndexType partition_index = std::get<0>(rGroupVectorSet[current_group_index]);

            auto& current_boundary_indices = rBoundaryIndices[current_group_index];

           // std::vector<PartitionBoxType>
            auto& current_bounding_box = rPartitionBox[current_group_index];
            //std::cout << current_bounding_box.first << "; " << current_bounding_box.second << std::endl;
            for( IndexType i = 0; i < rGroupVectorSet.size(); ++i ){
                auto other_partition_index = std::get<0>(rGroupVectorSet[i]);
                if( !rVisited[i] && std::abs( static_cast<int>(partition_index-other_partition_index)) <=1 ){

                    auto& other_bounding_box = rPartitionBox[i];
                    //if( Touch(current_bounding_box, other_bounding_box) ){
                        auto& other_boundary_indices = rBoundaryIndices[i];

                        bool are_boundary = false;
                        std::array<IndexType, 6> map_direction = {1, 0, 3, 2, 5, 4};
                        for( IndexType direction = 0; direction < 2; ++direction){
                            if( are_boundary ){
                                break;
                            }
                            for( auto& iii : current_boundary_indices[direction] ) { // -1{
                                auto next_index = GetNextIndex(direction, iii);
                                if( other_boundary_indices[map_direction[direction]].find(next_index) != other_boundary_indices[map_direction[direction]].end() ){
                                    are_boundary = true;
                                    break;
                                }
                            }
                        }

                        if( are_boundary ) {
                            auto& merged_group = rMergedGroups[rMergedGroups.size()-1];
                            auto &new_set = std::get<1>(rGroupVectorSet[i]);

                            std::get<1>(merged_group).insert( new_set.begin(), new_set.end() );
                            std::get<2>(merged_group) +=  std::get<2>(rGroupVectorSet[i]);

                            for( IndexType direction = 0; direction < 2; ++direction){
                                current_boundary_indices[0].insert(other_boundary_indices[direction].begin(),
                                    other_boundary_indices[direction].end() );
                            }
                            rVisited[i] = true;
                            index_stack.push(i);
                        }
                    //}
                }

            }



        }
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