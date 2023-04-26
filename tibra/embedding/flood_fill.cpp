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
typedef FloodFill::StatusVectorType StatusVectorType;

Unique<StatusVectorType> FloodFill::ClassifyElements() const {

    Timer timer{};
    const IndexType total_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];

    // Set all elements as outside.
    StatusVectorType states(total_num_elements);
    std::fill(states.begin(), states.end(), IntersectionStatus::Outside );

    // Get max num elements. Partition will happen along side with max_elements.
    IndexType max_num_elements_per_dir = std::max<IndexType>( std::max<IndexType>(
            mNumberOfElements[0], mNumberOfElements[1] ), mNumberOfElements[2] );

    // Direction along we apply partition.
    IndexType dir_index;
    std::array<IndexType, 2> directions;
    if( mNumberOfElements[0] == max_num_elements_per_dir ){
        dir_index = 0;
        directions[0] = 0; // +x
        directions[1] = 1; // -x
    } else if ( mNumberOfElements[1] == max_num_elements_per_dir ) {
        dir_index = 1;
        directions[0] = 2; // +y
        directions[1] = 3; // -y
    } else if ( mNumberOfElements[2] == max_num_elements_per_dir ) {
        dir_index = 2;
        directions[0] = 4; // +z
        directions[1] = 5; // -z
    }

    std::vector<PartitionBoxType> partition;
    dir_index = 2;

    IndexType num_threads = 1;
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }
    num_threads = 10;
    //std::cout << "num Threads: " << num_threads << std::endl;
    const IndexType partition_size = std::max<IndexType>( std::ceil( static_cast<double>(mNumberOfElements[dir_index]) /
        static_cast<double>(num_threads) ), 1);

    Vector3i partition_lower_bound(0, 0, 0);
    Vector3i partition_upper_bound = mNumberOfElements;
    for( IndexType i = 0; i < num_threads; ++i ){
        if( partition_size*i < mNumberOfElements[dir_index] ){
            partition_lower_bound[dir_index] = partition_size*i;
            partition_upper_bound[dir_index] = std::min<IndexType>(partition_size*(i+1), mNumberOfElements[dir_index]);
            partition.push_back( std::make_pair(partition_lower_bound,  partition_upper_bound) );
        }
    }

    GroupVectorSetType groups;
    BoolVectorType visited(total_num_elements, false);
    #pragma omp parallel for firstprivate(visited) schedule(static, 1)
    for( IndexType i = 0; i < partition.size(); ++i){
        PartitionedFill(i, groups, partition[i], visited, states);
    }


    GroupVectorSetType groups2;
    BoolVectorType visited2(total_num_elements, false);

    auto partition_2 = std::make_pair(Vector3i(0, 0, 0), mNumberOfElements);
    PartitionedFill(1000, groups2, partition_2, visited2, states);


    std::cout << "Num Groups: " << groups.size() << std::endl;
    GroupVectorSetType merged_groups;
    MergeGroups( dir_index, merged_groups, groups, states );

    std::cout << "Merged groups: " << merged_groups.size() << std::endl;
    if( merged_groups.size()  != groups2.size() ){
        TIBRA_ERROR("errorrrrrrrrrrr") << "Waf\n";
    }
    for( auto& r_group : merged_groups ){
        int inside_count = std::get<2>(r_group);
        //std::cout << std::get<1>(r_group).size() << ", count: " << inside_count << std::endl;
        auto state = (inside_count > 0) ? IntersectionStatus::Inside : IntersectionStatus::Outside;
        for( auto group_it = std::get<1>(r_group).begin(); group_it !=  std::get<1>(r_group).end(); ++group_it){
            states[*group_it] = state;
        }
    }

    //std::cout << "Total: " << timer.Measure() << std::endl;
    timer.Reset();
    #pragma omp parallel for
    for( IndexType i = 0; i < total_num_elements; ++i){
        auto indices = mIdMapper.GetMatrixIndicesFromVectorIndex(i);
        auto box = GetBoundingBoxFromIndex(i);
        mpBrepOperator->GetIntersectionState(box.first, box.second);
    }
    //std::cout << "Total Classical: " << timer.Measure() << std::endl;

    return MakeUnique<StatusVectorType>(states);
}


void FloodFill::PartitionedFill(IndexType PartitionIndex, GroupVectorSetType& rGroupVectorSet, PartitionBoxType rPartition, BoolVectorType& rVisited, StatusVectorType& rStates) const {
    // Loop through current partition
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

void FloodFill::SinglePartitionFill(IndexType Index, GroupSetType& rGroupSet, PartitionBoxType rPartition, BoolVectorType& rVisited, StatusVectorType& rStates) const {

    // Set Index as visited
    rVisited[Index] = true;
    const auto box = GetBoundingBoxFromIndex(Index);
    // Only start filling of current element is not trimmed.
    if( mpBrepOperator->IsTrimmed(box.first, box.second) ){
        rStates[Index] = IntersectionStatus::Trimmed;
    } else {
        std::get<1>(rGroupSet).insert(Index);

        // If box is not trimmed, run flood fill.
        IndexStackType index_stack;
        index_stack.push( Index );
        std::array<int, 6> new_indices;
        while( !index_stack.empty() ){
            /// 0:+x, 1:-x, 2:+y, 3:-y, 4:+z, 5:-z
            for( IndexType i = 0; i < 6; ++i){
                new_indices[i] = FillDirection(i, index_stack, rGroupSet, rVisited, rPartition, rStates );
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

int FloodFill::FillDirection(IndexType Direction, IndexStackType& rStack, GroupSetType& rGroupSet,
        BoolVectorType& rVisited, PartitionBoxType& rPartition, StatusVectorType& rStates ) const {

    const IndexType index = rStack.top();
    const Vector3i indices = mIdMapper.GetMatrixIndicesFromVectorIndex(index);
    const IndexType max_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];
    Vector3i next_indices = indices;

    PointType lower_perturb = {0.0, 0.0, 0.0};
    PointType upper_perturb = {0.0, 0.0, 0.0};
    double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);

    int next_index = GetNextIndex(Direction, index, rPartition, lower_perturb, upper_perturb);

    // Check if out-of-range
    if( next_index >= max_num_elements || next_index < 0 ) {
        return -1;
    }

    // If next box is trimmed, add inside count.
    // Note that next box is slightly shifted towards original box to capture trims directly at the boundary.
    auto box_next = GetBoundingBoxFromIndex(next_index);
    if( mpBrepOperator->IsTrimmed(box_next.first + lower_perturb , box_next.second + upper_perturb) ){
        std::get<2>(rGroupSet) += GetIsInsideCount(index, next_index, lower_perturb, upper_perturb);
        return -1;
    }

    // Already visited.
    if( rVisited[next_index] ) {
        return -1;
    }

    // Is trimmed.
    if( mpBrepOperator->IsTrimmed(box_next.first, box_next.second) ){
        rVisited[next_index] = true;
        rStates[next_index] = IntersectionStatus::Trimmed;
        return -1;
    }

    // Add next_index to current group.
    rVisited[next_index] = true;
    std::get<1>(rGroupSet).insert(next_index);

    return next_index;
}




void FloodFill::MergeGroups(IndexType PartitionDir, GroupVectorSetType& rMergedGroups, GroupVectorSetType& rGroupVectorSet, StatusVectorType& rStates) const {

    const IndexType num_groups = rGroupVectorSet.size();
    BoolVectorType visited(num_groups, false);

    double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);

    const std::array<IndexType, 2> walk_directions = {2*PartitionDir, (2*PartitionDir)+1};
    // Compute bounding box for each group (only along PartitionDir).
    std::vector<std::pair<IndexType, IndexType>> group_bounding_boxes(num_groups);
    std::fill(group_bounding_boxes.begin(), group_bounding_boxes.end(),
        std::make_pair( mNumberOfElements[PartitionDir]+1,  0 ) );

    BoundaryIndicesVectorType group_boundary_indices(num_groups);
    for( IndexType i = 0; i < num_groups; ++i){
        group_boundary_indices[i].resize(2);
    }

    #pragma omp parallel for
    for( IndexType group_index = 0; group_index < num_groups;  ++group_index){
        auto& group_set = rGroupVectorSet[group_index];
        const auto& index_set = std::get<1>(group_set);
        auto& group_bounding_box = group_bounding_boxes[group_index];
        auto& boundary_indices = group_boundary_indices[group_index];

        for(auto& index : index_set ){
            const auto indices = mIdMapper.GetMatrixIndicesFromVectorIndex(index);
            PointType lower_offset;
            PointType upper_offset;

            if( indices[PartitionDir] >= std::get<1>( group_bounding_box ) ){
                auto next_index = GetNextIndex(walk_directions[0], index, lower_offset, upper_offset );
                auto box_next = GetBoundingBoxFromIndex(next_index);
                if( !mpBrepOperator->IsTrimmed(box_next.first + lower_offset, box_next.second + upper_offset ) ) {
                    if( indices[PartitionDir] > std::get<1>( group_bounding_box ) ){
                        std::get<1>( group_bounding_box ) = indices[PartitionDir];
                        boundary_indices[1].clear();
                    }
                    boundary_indices[1].insert(index);
                } else {
                    // Slight error here. This should acutally only be applied to all boundary indices.
                    // This would require another loop.
                    std::get<2>(group_set) += GetIsInsideCount(index, next_index, lower_offset, upper_offset);
                }
            }
            if( indices[PartitionDir] <= std::get<0>( group_bounding_box ) ){
                auto next_index = GetNextIndex(walk_directions[1], index, lower_offset, upper_offset );
                auto box_next = GetBoundingBoxFromIndex(next_index);
                if( !mpBrepOperator->IsTrimmed(box_next.first + lower_offset, box_next.second + upper_offset ) ){
                    if( indices[PartitionDir] < std::get<0>( group_bounding_box ) ){
                        std::get<0>( group_bounding_box ) = indices[PartitionDir];
                        boundary_indices[0].clear();
                    }
                    boundary_indices[0].insert(index);
                } else {
                    std::get<2>(group_set) += GetIsInsideCount(index, next_index, lower_offset, upper_offset);
                }
            }
        }
    }

    for( IndexType group_index = 0; group_index < num_groups; ++group_index){
        if( !visited[group_index] ){
            GroupFill(PartitionDir, group_boundary_indices, rMergedGroups, group_index, group_bounding_boxes, rGroupVectorSet, visited, rStates);
        }
    }

}

    void FloodFill::GroupFill(IndexType PartitionDir, BoundaryIndicesVectorType& rBoundaryIndices, GroupVectorSetType& rMergedGroups, IndexType GroupIndex, std::vector<Partition1DBoxType>& rPartition1DBoxex, GroupVectorSetType& rGroupVectorSet, BoolVectorType& rVisited, StatusVectorType& rStates) const {

        const std::array<IndexType, 2> walk_directions = {2*PartitionDir, (2*PartitionDir)+1};

        // Set index as visited
        rVisited[GroupIndex] = true;

        // If box is not trimmed, run flood fill.
        IndexStackType index_stack;
        index_stack.push( GroupIndex );

        rMergedGroups.push_back( rGroupVectorSet[GroupIndex] );
        const std::array<IndexType, 6> map_direction = {1, 0, 3, 2, 5, 4};
        while( !index_stack.empty() ){
            IndexType current_group_index = index_stack.top();
            index_stack.pop();
            IndexType partition_index = std::get<0>(rGroupVectorSet[current_group_index]);

            auto& current_boundary_indices = rBoundaryIndices[current_group_index];

            // std::vector<PartitionBoxType>
            auto& current_partition_box = rPartition1DBoxex[current_group_index];
            for( IndexType i = 0; i < rGroupVectorSet.size(); ++i ){
                auto other_partition_index = std::get<0>(rGroupVectorSet[i]);
                if( !rVisited[i] && std::abs( static_cast<int>(partition_index-other_partition_index)) <=1 ){

                    auto& other_partition_box = rPartition1DBoxex[i];

                    //if( Touch(other_partition_box, current_partition_box) ){
                        auto& other_boundary_indices = rBoundaryIndices[i];
                        bool are_boundary = false;

                        for( auto direction : walk_directions){
                            if( are_boundary ){
                                break;
                            }

                            for( auto& iii : current_boundary_indices[map_direction[direction] % 2] ) { // -1{
                                auto next_index = GetNextIndex(direction, iii);
                                if( other_boundary_indices[direction % 2].find(next_index) != other_boundary_indices[direction % 2].end() ){
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

                            rVisited[i] = true;
                            index_stack.push(i);
                        }
                    //}
                }

            }

        }
    }



    std::pair<PointType, PointType> FloodFill::GetBoundingBoxFromIndex(IndexType Index) const{
        const auto indices =  mIdMapper.GetMatrixIndicesFromVectorIndex(Index);
        const PointType indices_d( indices[0], indices[1], indices[2] );
        return std::make_pair( mLowerBound + indices_d * mDelta,
                               mLowerBound + (indices_d+1.0) * mDelta );
    }


} // End namespace tibra