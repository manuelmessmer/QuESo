// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef FLOOD_FILL_INCLUDE_H
#define FLOOD_FILL_INCLUDE_H

//// STL includes
#include <array>
#include <algorithm>
#include <stack>

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

    typedef std::vector<IntersectionStatusType> StatusVectorType;
    typedef std::stack<IndexType> IndexStackType;

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
    void ClassifyElements( ) {

        StatusVectorType states(mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2]);
        std::fill(states.begin(), states.end(), IntersectionStatus::NotVisited );

        Timer timer{};


        for( IndexType i = 0; i < states.size(); ++i) {
            if( states[i] == IntersectionStatus::NotVisited ) {
                const auto box = GetBoundingBoxFromIndex(i);
                if( mpBrepOperator->IsTrimmed(box.first, box.second) ){
                    states[i] = IntersectionStatus::Trimmed;
                } else {
                    // If box is not trimmed, run flood fill.
                    std::stack<IndexType> index_stack;
                    index_stack.push(i);
                    states[i] = IntersectionStatus::PartOfNewGroup;
                    double criterion = Fill(index_stack, states);

                    IntersectionStatus current_state = (criterion > 0) ? IntersectionStatus::Inside : IntersectionStatus::Outside;
                    IntersectionStatus opposite_state = (current_state == IntersectionStatus::Inside) ? IntersectionStatus::Outside : IntersectionStatus::Inside;
                    int tmp_count = 0;
                    for( IndexType j = 0; j < states.size(); ++j) {
                        if( states[j] == IntersectionStatus::PartOfNewGroup ) {
                            states[j] = current_state;
                        }
                        if( states[j] == IntersectionStatus::OppositeGroup ) {
                            states[j] = opposite_state;
                        }
                    }
                }
            }
        }

        std::cout << "Timer: " << timer.Measure() << std::endl;

        #pragma omp parallel for
        for( IndexType i = 0; i < states.size(); ++i){
            auto box = GetBoundingBoxFromIndex(i);
            auto status_ref = mpBrepOperator->GetIntersectionState(box.first, box.second);

            if( status_ref != states[i] ){
                auto cube = MeshUtilities::pGetCuboid(box.first, box.second);
                IO::WriteMeshToSTL(*cube, "cube.stl", true);
                std::cout << "i: " << i << std::endl;
                std::cout << "status_ref: " << status_ref << ", " << states[i] << std::endl;
                TIBRA_ERROR("waf") << "test\n";
            }

        }
    }

    double Fill( IndexStackType& rIndexStack, StatusVectorType& rStates ){
        double criterion = 0.0;
        std::array<int, 6> new_indices;
        while( !rIndexStack.empty() ){

            for( IndexType i = 0; i < 6; ++i){
                new_indices[i] = FillDirection(i, rIndexStack, rStates, criterion );
            }

            rIndexStack.pop();
            for( IndexType i = 0; i < 6; ++i){
                if( new_indices[i] > -1 ){
                    rIndexStack.push(new_indices[i]);
                }
            }
        }
        return criterion;
    }

private:

    int FillDirection(IndexType Direction, IndexStackType& stack, StatusVectorType& rStates, double& rCriterion){
        bool local_end = false;
        IndexType index = stack.top();

        int next_index = -1;
        PointType perturb;
        double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);
        switch(Direction){
            case 0:
                next_index = mIdMapper.GetNextIndexX(index, local_end);
                perturb = {-tolerance, 0.0, 0.0};
                break;
            case 1:
                next_index = mIdMapper.GetPreviousIndexX(index, local_end);
                perturb = {tolerance, 0.0, 0.0};
                break;
            case 2:
                next_index = mIdMapper.GetNextIndexY(index, local_end);
                perturb = {0.0, -tolerance, 0.0};
                break;
            case 3:
                next_index = mIdMapper.GetPreviousIndexY(index, local_end);
                perturb = {0.0, tolerance, 0.0};
                break;
            case 4:
                next_index = mIdMapper.GetNextIndexZ(index, local_end);
                perturb = {0.0, 0.0, -tolerance};
                break;
            case 5:
                next_index = mIdMapper.GetPreviousIndexZ(index, local_end);
                perturb = {0.0, 0.0, tolerance};
                break;
            default:
                TIBRA_ERROR("error") << "error\n";
        }

        // Check if out-of-range
        if( next_index >= rStates.size() || next_index < 0 ) {
            return -1;
        }

        // Already classified
        if(    rStates[next_index] == IntersectionStatus::Trimmed
            || rStates[next_index] == IntersectionStatus::PartOfNewGroup
            || rStates[next_index] == IntersectionStatus::Outside
            || rStates[next_index] == IntersectionStatus::Inside ) {
            return -1;
        }


        auto indices =  mIdMapper.GetMatrixIndicesFromVectorIndex(next_index);
        PointType indices_d( indices[0], indices[1], indices[2] );
        PointType lower_bound = mLowerBound + indices_d * mDelta;
        PointType upper_bound = mLowerBound + (indices_d+1.0) * mDelta;

        if( mpBrepOperator->IsTrimmed(lower_bound, upper_bound) ){
            auto indices_tmp =  mIdMapper.GetMatrixIndicesFromVectorIndex(index);
            PointType indices_d_tmp( indices_tmp[0], indices_tmp[1], indices_tmp[2] );
            PointType lower_bound_tmp = mLowerBound + indices_d_tmp * mDelta;
            PointType upper_bound_tmp = mLowerBound + (indices_d_tmp+1.0) * mDelta;
            PointType center_box = (upper_bound_tmp + lower_bound_tmp)*0.5;

            auto p_mesh = mpBrepOperator->pClipTriangleMesh(lower_bound, upper_bound);

            for( IndexType k = 0; k < p_mesh->NumOfTriangles(); ++k ){
                auto center_tri = p_mesh->Center(k);
                auto normal = p_mesh->Normal(k);
                double area = p_mesh->Area(k);
                auto c_c = center_tri - center_box;

                //#pragma omp atomic update
                rCriterion += area / c_c.Norm() * Math::Dot(normal, c_c);
            }
            rStates[next_index] = IntersectionStatus::Trimmed;

            return -1;
        }

        if( mpBrepOperator->IsTrimmed(lower_bound + perturb , upper_bound + perturb) ){
            rStates[next_index] = IntersectionStatus::OppositeGroup;
            return -1;
        }

        rStates[next_index] = IntersectionStatus::PartOfNewGroup;
        return next_index;
    }

    std::pair<PointType, PointType> GetBoundingBoxFromIndex(IndexType Index){
        const auto indices =  mIdMapper.GetMatrixIndicesFromVectorIndex(Index);
        const PointType indices_d( indices[0], indices[1], indices[2] );
        return std::make_pair( mLowerBound + indices_d * mDelta,
                               mLowerBound + (indices_d+1.0) * mDelta );
    }
    // int FillPreviousX(IndexStackType& stack, StatusVectorType& rStates){
    //     bool local_end = false;
    //     IndexType index = stack.top();


    //     int next_index = mIdMapper.GetPreviousIndexX(index, local_end);
    //     if( next_index >= rStates.size() || next_index < 0 ){ return -1; }
    //     if( rStates[next_index] == IntersectionStatus::Trimmed || rStates[next_index] == IntersectionStatus::Outside) {return -1; }
    //     rStates[next_index] = IntersectionStatus::Outside;
    //     return next_index;
    // }

    // int FillNextY(IndexStackType& stack, StatusVectorType& rStates){
    //     bool local_end = false;
    //     IndexType index = stack.top();


    //     int next_index = mIdMapper.GetNextIndexY(index, local_end);
    //     if( next_index >= rStates.size() || next_index < 0 ){ return -1; }
    //     if( rStates[next_index] == IntersectionStatus::Trimmed || rStates[next_index] == IntersectionStatus::Outside) {return -1; }
    //     rStates[next_index] = IntersectionStatus::Outside;
    //     return next_index;
    // }

    // int FillPreviousY(IndexStackType& stack, StatusVectorType& rStates){
    //     bool local_end = false;
    //     IndexType index = stack.top();


    //     int next_index = mIdMapper.GetPreviousIndexY(index, local_end);
    //     if( next_index >= rStates.size() || next_index < 0 ){ return -1; }
    //     if( rStates[next_index] == IntersectionStatus::Trimmed || rStates[next_index] == IntersectionStatus::Outside) {return -1; }
    //     return next_index;
    // }

    // int FillNextZ(IndexStackType& stack, StatusVectorType& rStates){
    //     bool local_end = false;
    //     IndexType index = stack.top();


    //     int next_index = mIdMapper.GetNextIndexZ(index, local_end);
    //     if( next_index >= rStates.size() || next_index < 0 ){ return -1; }
    //     if( rStates[next_index] == IntersectionStatus::Trimmed || rStates[next_index] == IntersectionStatus::Outside) {return -1; }
    //     rStates[next_index] = IntersectionStatus::Outside;
    //     return next_index;
    // }

    // int FillPreviousZ(IndexStackType& stack, StatusVectorType& rStates){
    //     bool local_end = false;
    //     IndexType index = stack.top();


    //     int next_index = mIdMapper.GetPreviousIndexZ(index, local_end);
    //     if( next_index >= rStates.size() || next_index < 0 ){ return -1; }
    //     if( rStates[next_index] == IntersectionStatus::Trimmed || rStates[next_index] == IntersectionStatus::Outside) {return -1; }
    //     rStates[next_index] = IntersectionStatus::Outside;
    //     return next_index;
    // }

private:
    BRepOperatorBase* mpBrepOperator;
    VectorMatrixIdUtilities mIdMapper;
    const PointType mLowerBound;
    const PointType mUpperBound;
    const Vector3i mNumberOfElements;
    PointType mDelta;

}; // End class FloodFill

///@} TIBRA Classes
} // End tibra namespace

#endif // FLOOD_FILL_INCLUDE_H