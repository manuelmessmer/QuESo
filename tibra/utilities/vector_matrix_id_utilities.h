// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef VECTOR_MATRIX_ID_UTILITIES_INCLUDE_H
#define VECTOR_MATRIX_ID_UTILITIES_INCLUDE_H


//// Project incudes
#include "define.hpp"

namespace tibra {

///@name  TIBRA Classes
///@{

/**
 * @class  VectorMatrixIdUtilities
 * @author Manuel Messmer
 * @brief  Provides function to map beetwen matrix and vector indices.
**/
class VectorMatrixIdUtilities {

public:

    ///@name  Life cycle
    ///@{
    VectorMatrixIdUtilities(const Vector3i& rNumberOfElements) : mNumberOfElements(rNumberOfElements)
    {
    }

    ///@}
    ///@name  Operation
    ///@{
    IndexType GetNextIndexX(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[0] < mNumberOfElements[0]-1 ){
            local_end = false;
        } else {
            local_end = true;
        }
        return i + 1;
    }

    IndexType GetNextIndexY(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[1] < mNumberOfElements[1]-1) {
            indices[1] += 1;
            local_end = false;
        }
        else if( indices[0] < mNumberOfElements[0]-1){
            indices[0] += 1;
            indices[1] = 0;
            local_end = true;
        }
        else {
            indices[2] += 1;
            indices[1] = 0;
            indices[0] = 0;
            local_end = true;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    IndexType GetNextIndexZ(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[2] < mNumberOfElements[2]-1) {
            indices[2] += 1;
            local_end = false;
        }
        else if( indices[0] < mNumberOfElements[0]-1){
            indices[0] += 1;
            indices[2] = 0;
            local_end = true;
        }
        else {
            indices[1] += 1;
            indices[2] = 0;
            indices[0] = 0;
            local_end = true;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    IndexType GetPreviousIndexX(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[0] > 0) {
            local_end = false;
        } else {
            local_end = true;
        }
        return i - 1;
    }

    IndexType GetPreviousIndexY(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[1] > 0) {
            indices[1] -= 1;
            local_end = false;
        }
        else if( indices[0] > 0){
            indices[0] -= 1;
            indices[1] = mNumberOfElements[1]-1;
            local_end = true;
        }
        else {
            indices[2] -= 1;
            indices[1] = mNumberOfElements[1]-1;
            indices[0] = mNumberOfElements[0]-1;
            local_end = true;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    IndexType GetPreviousIndexZ(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[2] > 0 ){
            indices[2] -= 1;
            local_end = false;
        }
        else if( indices[0] > 0){
            indices[0] -= 1;
            indices[2] = mNumberOfElements[2]-1;
            local_end = true;
        }
        else {
            indices[1] -= 1;
            indices[2] = mNumberOfElements[2]-1;
            indices[0] = mNumberOfElements[0]-1;
            local_end = true;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    std::array<IndexType,3> GetMatrixIndicesFromVectorIndex(const IndexType Index) noexcept
    {
        std::array<IndexType,3> result;
        const IndexType index_in_row_column_plane = Index % (mNumberOfElements[0]*mNumberOfElements[1]);
        result[0] = index_in_row_column_plane % mNumberOfElements[0]; // row
        result[1] = index_in_row_column_plane / mNumberOfElements[0]; // column
        result[2] = Index / (mNumberOfElements[0]*mNumberOfElements[1]);   // depth

        return result;
    }

    IndexType GetVectorIndexFromMatrixIndices(const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex) noexcept
    {
        return DepthIndex * (mNumberOfElements[1]*mNumberOfElements[0]) + ColumnIndex * mNumberOfElements[0] + RowIndex;
    }

public:
    const Vector3i mNumberOfElements;

}; // End class Parameters
///@} End TIBRA classes

} // End namespace tibra

#endif // VECTOR_MATRIX_ID_UTILITIES_INCLUDE_H