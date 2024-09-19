//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef GRID_INDEXER_INCLUDE_HPP
#define GRID_INDEXER_INCLUDE_HPP

//// STL includes
#include <cassert>

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/includes/settings.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  GridIndexer
 * @author Manuel Messmer
 * @brief  Provides interface to quickly access elements in the background grid.
*/
class GridIndexer {
public:
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rSettings
    GridIndexer( const SettingsBaseType& rSettings ) :
        mBoundXYZ( std::make_pair(rSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz),
                                  rSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz)) ),
        mBoundUVW( std::make_pair(rSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw),
                                  rSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw)) ),
        mNumberOfElements(rSettings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::number_of_elements) ),
        mGlobalPartition( std::make_pair(Vector3i({0, 0, 0}), Vector3i({mNumberOfElements[0]-1, mNumberOfElements[1]-1, mNumberOfElements[2]-1}) )),
        mBSplineMesh( rSettings[MainSettings::background_grid_settings].GetValue<BackgroundGridType>(BackgroundGridSettings::grid_type) ==  BackgroundGridType::b_spline_grid )
    {
    }

    ///@}
    ///@name Public Operations
    ///@{

    enum class IndexState {okay, local_end, out_of_bounds};

    /// @brief Maps univariate index to trivariate indices i -> (i,j,k).
    /// @see GetVectorIndexFromMatrixIndices.
    /// @param Index
    /// @return Vector3i.
    inline Vector3i GetMatrixIndicesFromVectorIndex(const IndexType Index) const {
        Vector3i result;
        const IndexType index_in_row_column_plane = Index % (mNumberOfElements[0]*mNumberOfElements[1]);
        result[0] = index_in_row_column_plane % mNumberOfElements[0]; // row
        result[1] = index_in_row_column_plane / mNumberOfElements[0]; // column
        result[2] = Index / (mNumberOfElements[0]*mNumberOfElements[1]);   // depth

        return result;
    }

    /// @brief Maps trivariate indices to univariate index (i,j,k) -> i.
    /// @see GetMatrixIndicesFromVectorIndex.
    /// @param RowIndex
    /// @param ColumnIndex
    /// @param DepthIndex
    /// @return IndexType.
    inline IndexType GetVectorIndexFromMatrixIndices(const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex ) const {
        return DepthIndex * (mNumberOfElements[1]*mNumberOfElements[0]) + ColumnIndex * mNumberOfElements[0] + RowIndex;
    }

    /// @brief Maps trivariate indices to univariate index (i,j,k) -> i.
    /// @see GetMatrixIndicesFromVectorIndex.
    /// @param rIndices
    /// @return IndexType.
    inline IndexType GetVectorIndexFromMatrixIndices(const Vector3i& rIndices ) const {
        return rIndices[2] * (mNumberOfElements[1]*mNumberOfElements[0]) + rIndices[1] * mNumberOfElements[0] + rIndices[0];
    }

    /// @brief Creates bounding box in physical space from given index.
    /// @param Index
    /// @return BoundingBoxType.
    inline BoundingBoxType GetBoundingBoxXYZFromIndex(IndexType Index) const {
        const auto indices = GetMatrixIndicesFromVectorIndex(Index);
        return GetBoundingBoxFromIndex(indices[0], indices[1], indices[2], mBoundXYZ.first, mBoundXYZ.second);
    }

    /// @brief Creates bounding box in physical space from given indices.
    /// @param Indices
    /// @return BoundingBoxType.
    inline BoundingBoxType GetBoundingBoxXYZFromIndex(const Vector3i& rIndices) const {
        return GetBoundingBoxFromIndex(rIndices[0], rIndices[1], rIndices[2], mBoundXYZ.first, mBoundXYZ.second);
    }

    /// @brief Creates bounding box in physical space from given indices.
    /// @param i
    /// @param j
    /// @param k
    /// @return BoundingBoxType
    inline BoundingBoxType GetBoundingBoxXYZFromIndex(IndexType i, IndexType j, IndexType k) const {
        return GetBoundingBoxFromIndex(i, j, k, mBoundXYZ.first, mBoundXYZ.second);
    }

    /// @brief Creates bounding box in parametric space from given index.
    /// @param Index
    /// @return BoundingBoxType.
    inline BoundingBoxType GetBoundingBoxUVWFromIndex(IndexType Index) const {
        if( mBSplineMesh ) {
            const auto indices = GetMatrixIndicesFromVectorIndex(Index);
            return GetBoundingBoxFromIndex(indices[0], indices[1], indices[2], mBoundUVW.first, mBoundUVW.second);
        }
        return mBoundUVW;
    }

    /// @brief Creates bounding box in parametric space from given indices.
    /// @param Indices
    /// @return BoundingBoxType.
    inline BoundingBoxType GetBoundingBoxUVWFromIndex(const Vector3i& rIndices) const {
        if( mBSplineMesh ) {
            return GetBoundingBoxFromIndex(rIndices[0], rIndices[1], rIndices[2], mBoundUVW.first, mBoundUVW.second);
        }
        return mBoundUVW;
    }

    /// @brief Creates bounding box in parametric space from given indices.
    /// @param i
    /// @param j
    /// @param k
    /// @return BoundingBoxType
    inline BoundingBoxType GetBoundingBoxUVWFromIndex(IndexType i, IndexType j, IndexType k) const {
        if( mBSplineMesh ) {
            return GetBoundingBoxFromIndex(i, j, k, mBoundUVW.first, mBoundUVW.second);
        }
        return mBoundUVW;
    }

    /// @brief Returns global number of elements (including inactive elements).
    /// @return IndexType.
    inline IndexType NumberOfElements() const {
        return mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];
    }

    /// @brief Returns next index in given direction.
    /// @param Index current index.
    /// @param Direction Move Direction: 0:+x, 1:-x, 2:+y, 3:-y, 4:+z, 5:-z
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetNextIndex(IndexType Index, IndexType Direction, bool& LocalEnd) const {
        switch(Direction){
            case 0:
                return GetNextIndexX(Index, LocalEnd);
            case 1:
                return GetPreviousIndexX(Index, LocalEnd);
            case 2:
                return GetNextIndexY(Index, LocalEnd);
            case 3:
                return GetPreviousIndexY(Index, LocalEnd);
            case 4:
                return GetNextIndexZ(Index, LocalEnd);
            case 5:
                return GetPreviousIndexZ(Index, LocalEnd);
            default:
                assert(false);
        }
    }

    /// @brief Returns next index in given direction and partition.
    /// @param Index current index.
    /// @param rPartition active partition.
    /// @param Direction Move Direction: 0:+x, 1:-x, 2:+y, 3:-y, 4:+z, 5:-z
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetNextIndex(IndexType Index, IndexType Direction, const PartitionBoxType& rPartition, bool& LocalEnd) const {
        switch(Direction){
            case 0:
                return GetNextIndexX(Index, rPartition, LocalEnd);
            case 1:
                return GetPreviousIndexX(Index, rPartition, LocalEnd);
            case 2:
                return GetNextIndexY(Index, rPartition, LocalEnd);
            case 3:
                return GetPreviousIndexY(Index, rPartition, LocalEnd);
            case 4:
                return GetNextIndexZ(Index, rPartition, LocalEnd);
            case 5:
                return GetPreviousIndexZ(Index, rPartition, LocalEnd);
            default:
                assert(false);
        }
    }

    /// @brief Returns next index in x-direction.
    /// @param i current index.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetNextIndexX(IndexType i,  bool& rLocalEnd) const {
        return GetNextIndexX(i, mGlobalPartition, rLocalEnd);
    }

    /// @brief Returns next index in x-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetNextIndexX(IndexType i, const PartitionBoxType& rPartition, bool& rLocalEnd) const {
        rLocalEnd = false;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[0] < rPartition.second[0]) {
            indices[0] += 1;
        }
        else if (indices[1] < rPartition.second[1]) {
            indices[0] = rPartition.first[0];
            indices[1] += 1;
            rLocalEnd = true;
        }
        else if (indices[2] < rPartition.second[2]) {
            indices[2] += 1;
            indices[1] = rPartition.first[1];
            indices[0] = rPartition.first[0];
            rLocalEnd = true;
        } else {
            assert(false); // Should not be reached.
        }

        return GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);
    }

    /// @brief Returns next index in y-direction.
    /// @param i current index.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetNextIndexY(IndexType i,  bool& rLocalEnd) const {
        return GetNextIndexY(i, mGlobalPartition, rLocalEnd);
    }

    /// @brief Returns next index in y-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetNextIndexY(IndexType i, const PartitionBoxType& rPartition, bool& rLocalEnd) const {
        rLocalEnd = false;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[1] < rPartition.second[1]) {
            indices[1] += 1;
        }
        else if( indices[0] < rPartition.second[0]){
            indices[0] += 1;
            indices[1] = rPartition.first[1];
            rLocalEnd = true;
        }
        else if( indices[2] < rPartition.second[2] ) {
            indices[2] += 1;
            indices[1] = rPartition.first[1];
            indices[0] = rPartition.first[0];
            rLocalEnd = true;
        } else {
            assert(false); // Should not be reached.
        }

        return GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);
    }

    /// @brief Returns next index in z-direction.
    /// @param i current index.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetNextIndexZ(IndexType i,  bool& rLocalEnd) const {
        return GetNextIndexZ(i, mGlobalPartition, rLocalEnd);
    }

    /// @brief Returns next index in z-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetNextIndexZ(IndexType i, const PartitionBoxType& rPartition, bool& rLocalEnd) const {
        rLocalEnd = false;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[2] < rPartition.second[2]) {
            indices[2] += 1;
        }
        else if( indices[0] < rPartition.second[0]){
            indices[0] += 1;
            indices[2] = rPartition.first[2];
            rLocalEnd = true;
        }
        else if( indices[1] < rPartition.second[1]) {
            indices[1] += 1;
            indices[2] = rPartition.first[2];
            indices[0] = rPartition.first[0];
            rLocalEnd = true;
        } else {
            assert(false); // Should not be reached.
        }

        return GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);
    }


    /// @brief Returns previous index in x-direction.
    /// @param i current index..
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetPreviousIndexX(IndexType i, bool& rLocalEnd) const {
        return GetPreviousIndexX(i, mGlobalPartition, rLocalEnd);
    }

    /// @brief Returns previous index in x-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetPreviousIndexX(IndexType i, const PartitionBoxType& rPartition, bool& rLocalEnd) const {
        rLocalEnd = false;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[0] > rPartition.first[0] ) {
            indices[0] -= 1;
        }
        else if( indices[1] > rPartition.first[1] ){
            indices[0] = rPartition.second[0];
            indices[1] -= 1;
            rLocalEnd = true;
        }
        else if( indices[2] > rPartition.first[2] ){
            indices[2] -= 1;
            indices[1] = rPartition.second[1];
            indices[0] = rPartition.second[0];
            rLocalEnd = true;
        } else {
            assert(false); // Should not be reached.
        }

        return GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);
    }

    /// @brief Returns previous index in y-direction.
    /// @param i current index.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetPreviousIndexY(IndexType i, bool& rLocalEnd) const {
        return GetPreviousIndexY(i, mGlobalPartition, rLocalEnd);
    }

    /// @brief Returns previous index in y-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetPreviousIndexY(IndexType i, const PartitionBoxType& rPartition, bool& rLocalEnd) const {
        rLocalEnd = false;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[1] > rPartition.first[1] ) {
            indices[1] -= 1;
        }
        else if( indices[0] > rPartition.first[0] ){
            indices[0] -= 1;
            indices[1] = rPartition.second[1];
            rLocalEnd = true;
        }
        else if( indices[2] > rPartition.first[2] ){
            indices[2] -= 1;
            indices[1] = rPartition.second[1];
            indices[0] = rPartition.second[0];
            rLocalEnd = true;
        } else {
            assert(false); // Should not be reached.
        }

        return GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);
    }

    /// @brief Returns previous index in z-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetPreviousIndexZ(IndexType i, bool& rLocalEnd) const {
        return GetPreviousIndexZ(i, mGlobalPartition, rLocalEnd);
    }

    /// @brief Returns previous index in z-direction.
    /// @param i current index.
    /// @param[out] LocalEnd True if current index is a local end.
    /// @return IndexType
    inline IndexType GetPreviousIndexZ(IndexType i, const PartitionBoxType& rPartition, bool& rLocalEnd) const {
        rLocalEnd = false;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[2] > rPartition.first[2] ){
            indices[2] -= 1;
        }
        else if( indices[0] > rPartition.first[0]){
            indices[0] -= 1;
            indices[2] = rPartition.second[2];
            rLocalEnd = true;
        }
        else if( indices[1] > rPartition.first[1]){
            indices[1] -= 1;
            indices[2] = rPartition.second[2];
            indices[0] = rPartition.second[0];
            rLocalEnd = true;
        } else {
            assert(false);
        }

        return GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);
    }

    ///@}
private:

    inline BoundingBoxType GetBoundingBoxFromIndex(IndexType i, IndexType j, IndexType k, const PointType& rLowerBound, const PointType& rUpperBound) const {
        const PointType indices_d{ static_cast<double>(i), static_cast<double>(j), static_cast<double>(k) };
        PointType delta;
        delta[0] = std::abs(rUpperBound[0] - rLowerBound[0]) / (mNumberOfElements[0]);
        delta[1] = std::abs(rUpperBound[1] - rLowerBound[1]) / (mNumberOfElements[1]);
        delta[2] = std::abs(rUpperBound[2] - rLowerBound[2]) / (mNumberOfElements[2]);
        return std::make_pair( Math::Add( rLowerBound, Math::MultElementWise(delta, indices_d)),
                            Math::Add( rLowerBound, Math::MultElementWise(delta, Math::Add({1.0, 1.0, 1.0}, indices_d))) );

    }
    ///@name Private Members
    ///@{

    const BoundingBoxType mBoundXYZ;
    const BoundingBoxType mBoundUVW;
    const Vector3i mNumberOfElements;
    const PartitionBoxType mGlobalPartition;
    const bool mBSplineMesh;

    ///@}
}; // End class GridIndexer.
///@} End queso classes.

} // End namespace queso

#endif // GRID_INDEXER_INCLUDE_HPP