// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

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
        std::fill(states.begin(), states.end(), IntersectionStatus::NotVisited );

        #pragma omp parallel for
        for( IndexType i = 0; i < mNumberOfElements[0]; ++i ){
            for( IndexType j = 0; j < mNumberOfElements[1]; ++j ) {
                for( IndexType k = 0; k < mNumberOfElements[2]; ++k ) {
                    const IndexType index = mIdMapper.GetVectorIndexFromMatrixIndices(i, j, k);
                    if( states[index] == IntersectionStatus::NotVisited ) {
                        states[index] = IntersectionStatus::PartOfNewGroup;
                        const auto box = GetBoundingBoxFromIndex(i,j,k);
                        if( mpBrepOperator->IsTrimmed(box.first, box.second) ){
                            states[index] = IntersectionStatus::Trimmed;
                        } else {
                            // If box is not trimmed, run flood fill.
                            IndexStackType index_stack;

                            index_stack.push( index );
                            int inside_count = Fill2D(index_stack, states);

                            IntersectionStatus current_state = (inside_count > 0) ? IntersectionStatus::Inside : IntersectionStatus::Outside;
                            IntersectionStatus opposite_state = (current_state == IntersectionStatus::Inside) ? IntersectionStatus::Outside : IntersectionStatus::Inside;
                            IndexType tmp_count = 0;
                            //for( IndexType i_ = 0; i_ < mNumberOfElements[0]; ++i_ ) {
                                for( IndexType j_ = 0; j_ < mNumberOfElements[1]; ++j_ ) {
                                    for( IndexType k_ = 0; k_ < mNumberOfElements[2]; ++k_ ) {
                                        const auto index_ = mIdMapper.GetVectorIndexFromMatrixIndices(i, j_, k_);
                                        if( states[index_] == IntersectionStatus::PartOfNewGroup ) {
                                            states[index_] = current_state;
                                            // auto bounding_box = GetBoundingBoxFromIndex(j);
                                            // auto low_el = Mapping::GlobalToParam(bounding_box.first, mLowerBound, mUpperBound);
                                            // auto up_el = Mapping::GlobalToParam(bounding_box.second, mLowerBound, mUpperBound);

                                            // Shared<Element> el_ptr = MakeShared<Element>(j, low_el, up_el, params);
                                            // elcont.AddElement(el_ptr);
                                            tmp_count++;
                                        }
                                        if( states[index_] == IntersectionStatus::OppositeGroup ) {
                                            states[index_] = opposite_state;
                                        }
                                    }
                                }
                            //}
                            // std::cout << "Group Size: " << tmp_count << std::endl;
                            // std::cout << "inside_count: " << inside_count << std::endl;
                            // std::cout << "----------------: " << inside_count << std::endl;
                        }
                    }
                }
            }

            // Merge
        }
        //std::cout << "Timer: " << timer.Measure() << std::endl;

        // Start inner test!
        // TIBRA_INFO << "Start Inner Test\n";
        // for( IndexType i = 0; i < states.size(); ++i){
        //     auto box = GetBoundingBoxFromIndex(i);
        //     auto status_ref = mpBrepOperator->GetIntersectionState(box.first, box.second);

        //     if( status_ref != states[i] ){
        //         auto cube = MeshUtilities::pGetCuboid(box.first, box.second);
        //         IO::WriteMeshToSTL(*cube, "cube.stl", true);
        //         std::cout << "i: " << i << std::endl;
        //         std::cout << "status_ref: " << status_ref << ", " << states[i] << std::endl;
        //         TIBRA_ERROR("waf") << "test\n";
        //     }

        // }
        return MakeUnique<StatusVectorType>(states);
    }

    int FloodFill::Fill2D( IndexStackType& rIndexStack, StatusVectorType& rStates ) const {
        int inside_count = 0;
        std::array<int, 6> new_indices;
        while( !rIndexStack.empty() ){

            for( IndexType i = 2; i < 6; ++i){
                new_indices[i] = FillDirection(i, rIndexStack, rStates, inside_count );
            }

            rIndexStack.pop();
            for( IndexType i = 2; i < 6; ++i){
                if( new_indices[i] > -1 ){
                    rIndexStack.push(new_indices[i]);
                }
            }
        }

        return inside_count;
    }


    int FloodFill::Fill( IndexStackType& rIndexStack, StatusVectorType& rStates ) const {
        int inside_count = 0;
        std::array<int, 6> new_indices;
        while( !rIndexStack.empty() ){

            for( IndexType i = 0; i < 6; ++i){
                new_indices[i] = FillDirection(i, rIndexStack, rStates, inside_count );
            }

            rIndexStack.pop();
            for( IndexType i = 0; i < 6; ++i){
                if( new_indices[i] > -1 ){
                    rIndexStack.push(new_indices[i]);
                }
            }
        }

        return inside_count;
    }

    // int FloodFill::Fill( IndexStackType& rIndexStack, StatusVectorType& rStates ) const {
    //     int inside_count = 0;
    //     std::array<int, 6> new_indices;
    //     while( !rIndexStack.empty() ){

    //         for( IndexType i = 0; i < 6; ++i){
    //             new_indices[i] = FillDirection(i, rIndexStack, rStates, inside_count );
    //         }

    //         rIndexStack.pop();
    //         for( IndexType i = 0; i < 6; ++i){
    //             if( new_indices[i] > -1 ){
    //                 rIndexStack.push(new_indices[i]);
    //             }
    //         }
    //     }
    //     return inside_count;
    // }



    int FloodFill::FillDirection(IndexType Direction, IndexStackType& stack, StatusVectorType& rStates, int& rInsideCount) const{
        bool local_end = false;
        const IndexType index = stack.top();
        const Vector3i indices = mIdMapper.GetMatrixIndicesFromVectorIndex(index);
        Vector3i next_indices = indices;

        PointType lower_perturb = {0.0, 0.0, 0.0};
        PointType upper_perturb = {0.0, 0.0, 0.0};
        double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);
        //std::cout << mNumberOfElements << std::endl;
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
        if( next_index >= rStates.size() || next_index < 0 ) {
            return -1;
        }
        // Already classified
        if( rStates[next_index] != IntersectionStatus::NotVisited ) {
            return -1;
        }

        auto box_next = GetBoundingBoxFromIndex(next_indices);

        if( mpBrepOperator->IsTrimmed(box_next.first + lower_perturb , box_next.second + upper_perturb) ){
            auto box_current = GetBoundingBoxFromIndex(indices);
            PointType center_box = (box_current.first + box_current.second)*0.5;
            if( mpBrepOperator->OnBoundedSideOfClippedSection(center_box, box_next.first + lower_perturb , box_next.second + upper_perturb) ) {
                rInsideCount += 1;
            } else {
                rInsideCount -= 1;
            }
            return -1;
        }

        if( mpBrepOperator->IsTrimmed(box_next.first, box_next.second) ){
            rStates[next_index] = IntersectionStatus::Trimmed;
            return -1;
        }

        if( mpBrepOperator->IsTrimmed(box_next.first + lower_perturb , box_next.second + upper_perturb) ){
            rStates[next_index] = IntersectionStatus::OppositeGroup;
            return -1;
        }

        rStates[next_index] = IntersectionStatus::PartOfNewGroup;
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