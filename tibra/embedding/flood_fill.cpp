// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// External includes
#include <omp.h>
#include <map>
//// Project includes
#include "embedding/flood_fill.h"


namespace tibra {

typedef FloodFill::StatusVectorType StatusVectorType;

Unique<StatusVectorType> FloodFill::ClassifyElements() const {

        Timer timer{};
        Parameters params( {Component("lower_bound", mLowerBound),
                            Component("upper_bound", mUpperBound),
                            Component("number_of_elements", mNumberOfElements) });
        ElementContainer elcont(params);

        IndexType max_num_elements = Math::Max(mNumberOfElements);
        StatusVectorType states(mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2]);
        std::fill(states.begin(), states.end(), std::make_pair( 0, IntersectionStatus::NotVisited) );


        IndexType total_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];
        GroupVectorType groups(total_num_elements);

        std::fill(groups.begin(), groups.end(), std::make_tuple<int, int, int>(-1, -1, 0));

        std::vector<IndexType> max_group_ids(mNumberOfElements[0]);
        #pragma omp parallel for
        for( IndexType i = 0; i < mNumberOfElements[0]; ++i ){
            IndexType group_id = i * mNumberOfElements[1]*mNumberOfElements[2]; // Ensure unique ID.
            for( IndexType j = 0; j < mNumberOfElements[1]; ++j ) {
                for( IndexType k = 0; k < mNumberOfElements[2]; ++k ) {
                    const IndexType index = mIdMapper.GetVectorIndexFromMatrixIndices(i, j, k);
                    if( std::get<0>(groups[index]) <= -1 ) { // Unvisited
                        Fill(index, groups, group_id);
                        if ( std::get<1>(groups[index]) == -1 ){
                            states[index].second = IntersectionStatus::Trimmed;
                        }
                    }
                    max_group_ids[i] = group_id;
                }
            }
        }

        double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);
        PointType offset(-tolerance, 0.0, 0.0);
        // Group_id and count.
        std::map<IndexType, int> group_inside_count;

        //#pragma omp parallel for
        for( IndexType j = 0; j < mNumberOfElements[1]; ++j ) {
            for( IndexType k = 0; k < mNumberOfElements[2]; ++k ) {
                for( IndexType i = 0; i < mNumberOfElements[0]; ++i ){
                    const IndexType index = mIdMapper.GetVectorIndexFromMatrixIndices(i, j, k);
                    const int group_id = std::get<1>(groups[index]);
                    if( group_id > -1 ){
                        if (group_inside_count.find(group_id) == group_inside_count.end()) {
                            group_inside_count.insert(std::make_pair(group_id, std::get<2>(groups[index])) );
                        } else {
                            group_inside_count[group_id] += std::get<2>(groups[index]);
                        }
                    }
                    if( i < mNumberOfElements[0]-1) {
                        const auto box = GetBoundingBoxFromIndex(index);
                        const PointType center = (box.first + box.second) * 0.5;
                        const IndexType index_next = mIdMapper.GetVectorIndexFromMatrixIndices(i+1, j, k);
                        const auto box_next = GetBoundingBoxFromIndex(index_next);
                        const int group_id_next = std::get<1>(groups[index_next]);

                        if( group_id > -1 && group_id_next > -1 ){
                            if ( !mpBrepOperator->IsTrimmed(box_next.first+offset, box_next.second) ){

                                // std::get<1>(groups[index_next]) = group_id;
                                // muss fuer alle auf der ebenen j, k passieren.
                                // make group id reference to something..

                                /// Das Funktioneirt noch nciht!!

                                // if( mpBrepOperator->OnBoundedSideOfClippedSection(center, box_next.first+offset, box_next.second+offset) ){
                                //     /// reduced next
                                //     if (group_inside_count.find(group_id) == group_inside_count.end()) {
                                //         group_inside_count.insert(std::make_pair(group_id, 1));
                                //     } else {
                                //         group_inside_count[group_id] += 1;
                                //     }
                                // } else {
                                //     /// Update next
                                //     if (group_inside_count.find(group_id) == group_inside_count.end()) {
                                //         group_inside_count.insert(std::make_pair(group_id, -1));
                                //     } else {
                                //         group_inside_count[group_id] -= 1;
                                //     }
                                // }
                            } else {
                                 /// Das Funktioneirt noch nciht!!
                            }
                        }
                    }
                }
            }
        }

        // for( auto it = group_inside_count.begin(); it != group_inside_count.end(); ++it){
        //     std::cout << it->first << ", " << it->second << std::endl;
        // }

        #pragma omp parallel for
        for( IndexType index = 0; index < total_num_elements; ++index) {
            const int group_id = std::get<1>(groups[index]);
            if( group_id > -1 ){
                const int inside_count = group_inside_count[group_id];
                IntersectionStatus status = (inside_count > 0) ? IntersectionStatus::Inside : IntersectionStatus::Outside;
                states[index].second = status;
            } else {
                states[index].second = IntersectionStatus::Trimmed;
            }
        }

        //std::cout << "Timer: " << timer.Measure() << std::endl;
        return MakeUnique<StatusVectorType>(states);
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