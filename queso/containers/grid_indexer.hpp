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

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/containers/dictionary.hpp"

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

    enum class IndexInfo {middle, local_end, global_end};
    enum class Direction {x_forward, x_backward, y_forward, y_backward, z_forward, z_backward};

    using IndexReturnType = std::pair<IndexType, IndexInfo>;
    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;

    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rSettings
    GridIndexer( const MainDictionaryType& rSettings ) :
        mBoundXYZ( std::make_pair(rSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_xyz),
                                  rSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_xyz)) ),
        mBoundUVW( std::make_pair(rSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::lower_bound_uvw),
                                  rSettings[MainSettings::background_grid_settings].GetValue<PointType>(BackgroundGridSettings::upper_bound_uvw)) ),
        mNumberOfElements(rSettings[MainSettings::background_grid_settings].GetValue<Vector3i>(BackgroundGridSettings::number_of_elements) ),
        mGlobalPartition( std::make_pair(Vector3i({0, 0, 0}), Vector3i({mNumberOfElements[0]-1, mNumberOfElements[1]-1, mNumberOfElements[2]-1}) )),
        mBSplineMesh( rSettings[MainSettings::background_grid_settings].GetValue<GridType>(BackgroundGridSettings::grid_type) ==  GridType::b_spline_grid )
    {
    }

    ///@}
    ///@name Public Operations
    ///@{

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
    /// @param Direction Move Direction.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetNextIndex(IndexType Index, Direction Dir) const {
        return GetNextIndex(Index, Dir, mGlobalPartition);
    }

    /// @brief Returns next index in given direction and partition.
    /// @param Index current index.
    /// @param Direction Move Direction.
    /// @param rPartition active partition.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetNextIndex(IndexType Index, Direction Dir, const PartitionBoxType& rPartition) const {
        switch(Dir){
            case Direction::x_forward:
                return GetNextIndexX(Index, rPartition);
            case Direction::x_backward:
                return GetPreviousIndexX(Index, rPartition);
            case Direction::y_forward:
                return GetNextIndexY(Index, rPartition);
            case Direction::y_backward:
                return GetPreviousIndexY(Index, rPartition);
            case Direction::z_forward:
                return GetNextIndexZ(Index, rPartition);
            case Direction::z_backward:
                return GetPreviousIndexZ(Index, rPartition);
            default:
                QuESo_ASSERT(false, "Should no be reached");
                return std::make_pair(0, IndexInfo::middle);
        }
    }

    /// @brief Returns next index in x-direction.
    /// @param i current index.
    /// @return IndexType
    inline IndexReturnType GetNextIndexX(IndexType i) const {
        return GetNextIndexX(i, mGlobalPartition);
    }

    /// @brief Returns next index in x-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetNextIndexX(IndexType i, const PartitionBoxType& rPartition) const {
        IndexInfo index_info = IndexInfo::middle;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[0] < rPartition.second[0]) {
            indices[0] += 1;
        }
        else if (indices[1] < rPartition.second[1]) {
            indices[0] = rPartition.first[0];
            indices[1] += 1;
            index_info = IndexInfo::local_end;
        }
        else if (indices[2] < rPartition.second[2]) {
            indices[2] += 1;
            indices[1] = rPartition.first[1];
            indices[0] = rPartition.first[0];
            index_info = IndexInfo::local_end;
        } else {
            index_info = IndexInfo::global_end;
            return std::make_pair(i, index_info);
        }
        return std::make_pair(GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]), index_info);
    }

    /// @brief Returns next index in y-direction.
    /// @param i current index.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetNextIndexY(IndexType i) const {
        return GetNextIndexY(i, mGlobalPartition);
    }

    /// @brief Returns next index in y-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetNextIndexY(IndexType i, const PartitionBoxType& rPartition) const {
        IndexInfo index_info = IndexInfo::middle;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[1] < rPartition.second[1]) {
            indices[1] += 1;
        }
        else if( indices[0] < rPartition.second[0]){
            indices[0] += 1;
            indices[1] = rPartition.first[1];
            index_info = IndexInfo::local_end;
        }
        else if( indices[2] < rPartition.second[2] ) {
            indices[2] += 1;
            indices[1] = rPartition.first[1];
            indices[0] = rPartition.first[0];
            index_info = IndexInfo::local_end;
        } else {
            index_info = IndexInfo::global_end;
            return std::make_pair(i, index_info);
        }
        return std::make_pair(GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]), index_info);
    }

    /// @brief Returns next index in z-direction.
    /// @param i current index.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetNextIndexZ(IndexType i) const {
        return GetNextIndexZ(i, mGlobalPartition);
    }

    /// @brief Returns next index in z-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetNextIndexZ(IndexType i, const PartitionBoxType& rPartition) const {
        IndexInfo index_info = IndexInfo::middle;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[2] < rPartition.second[2]) {
            indices[2] += 1;
        }
        else if( indices[0] < rPartition.second[0]){
            indices[0] += 1;
            indices[2] = rPartition.first[2];
            index_info = IndexInfo::local_end;
        }
        else if( indices[1] < rPartition.second[1]) {
            indices[1] += 1;
            indices[2] = rPartition.first[2];
            indices[0] = rPartition.first[0];
            index_info = IndexInfo::local_end;
        } else {
            index_info = IndexInfo::global_end;
            return std::make_pair(i, index_info);
        }
        return std::make_pair(GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]), index_info);
    }


    /// @brief Returns previous index in x-direction.
    /// @param i current index..
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetPreviousIndexX(IndexType i) const {
        return GetPreviousIndexX(i, mGlobalPartition);
    }

    /// @brief Returns previous index in x-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetPreviousIndexX(IndexType i, const PartitionBoxType& rPartition) const {
        IndexInfo index_info = IndexInfo::middle;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[0] > rPartition.first[0] ) {
            indices[0] -= 1;
        }
        else if( indices[1] > rPartition.first[1] ){
            indices[0] = rPartition.second[0];
            indices[1] -= 1;
            index_info = IndexInfo::local_end;
        }
        else if( indices[2] > rPartition.first[2] ){
            indices[2] -= 1;
            indices[1] = rPartition.second[1];
            indices[0] = rPartition.second[0];
            index_info = IndexInfo::local_end;
        } else {
            index_info = IndexInfo::global_end;
            return std::make_pair(i, index_info);
        }
        return std::make_pair(GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]), index_info);
    }

    /// @brief Returns previous index in y-direction.
    /// @param i current index.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetPreviousIndexY(IndexType i) const {
        return GetPreviousIndexY(i, mGlobalPartition);
    }

    /// @brief Returns previous index in y-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetPreviousIndexY(IndexType i, const PartitionBoxType& rPartition) const {
        IndexInfo index_info = IndexInfo::middle;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[1] > rPartition.first[1] ) {
            indices[1] -= 1;
        }
        else if( indices[0] > rPartition.first[0] ){
            indices[0] -= 1;
            indices[1] = rPartition.second[1];
            index_info = IndexInfo::local_end;
        }
        else if( indices[2] > rPartition.first[2] ){
            indices[2] -= 1;
            indices[1] = rPartition.second[1];
            indices[0] = rPartition.second[0];
            index_info = IndexInfo::local_end;
        } else {
            index_info = IndexInfo::global_end;
            return std::make_pair(i, index_info);
        }
        return std::make_pair(GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]), index_info);
    }

    /// @brief Returns previous index in z-direction.
    /// @param i current index.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetPreviousIndexZ(IndexType i) const {
        return GetPreviousIndexZ(i, mGlobalPartition);
    }

    /// @brief Returns previous index in z-direction.
    /// @param i current index.
    /// @param rPartition active partition.
    /// @return std::pair<IndexType, IndexInfo>: IndexInfo relates to the CURRENT Index and can be: middle, local_end, global_end.
    inline IndexReturnType GetPreviousIndexZ(IndexType i, const PartitionBoxType& rPartition) const {
        IndexInfo index_info = IndexInfo::middle;
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        if( indices[2] > rPartition.first[2] ){
            indices[2] -= 1;
        }
        else if( indices[0] > rPartition.first[0]){
            indices[0] -= 1;
            indices[2] = rPartition.second[2];
            index_info = IndexInfo::local_end;
        }
        else if( indices[1] > rPartition.first[1]){
            indices[1] -= 1;
            indices[2] = rPartition.second[2];
            indices[0] = rPartition.second[0];
            index_info = IndexInfo::local_end;
        } else {
            index_info = IndexInfo::global_end;
            return std::make_pair(i, index_info);
        }
        return std::make_pair(GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]), index_info);
    }

    /// @brief Returns true if current index is a local end w.r.t. the given direction.
    /// @param i current index
    /// @param Direction Move Direction.
    /// @return bool
    bool IsEnd(IndexType i, Direction Dir) const {
        return IsEnd(i, Dir, mGlobalPartition);
    }

    /// @brief Returns true if current index is a local end w.r.t. the given direction.
    /// @param i current index
    /// @param Direction Move Direction.
    /// @param rPartition active Partition
    /// @return bool
    bool IsEnd(IndexType i, Direction Dir, const PartitionBoxType& rPartition ) const {
        auto indices = GetMatrixIndicesFromVectorIndex(i);
        switch( Dir )
        {
        case Direction::x_forward:
            return (indices[0] == (rPartition.second[0]));
        case Direction::x_backward:
            return (indices[0] == (rPartition.first[0]));
        case Direction::y_forward:
            return (indices[1] == (rPartition.second[1]));
        case Direction::y_backward:
            return (indices[1] == (rPartition.first[1]));
        case Direction::z_forward:
            return (indices[2] == (rPartition.second[2]));
        case Direction::z_backward:
            return (indices[2] == (rPartition.first[2]));
        default:
            QuESo_ASSERT(false, "Should no be reached");
            return false;
        }
    }

    /// @brief Helper function that returns an array of all valid directions.
    /// @return std::array<Direction, 6>
    static constexpr std::array<Direction, 6> GetDirections() {
        return {Direction::x_forward, Direction::x_backward,
                Direction::y_forward, Direction::y_backward,
                Direction::z_forward, Direction::z_backward};
    }

    /// @brief Helper function that returns an array of all valid directions in reverse order.
    /// @return std::array<Direction, 6>
    static constexpr Direction ReverseDirection(Direction Dir) {
        constexpr std::array<Direction, 6> map_direction = {
            Direction::x_backward, Direction::x_forward,
            Direction::y_backward, Direction::y_forward,
            Direction::z_backward, Direction::z_forward
        };
        return map_direction[static_cast<IndexType>(Dir)];
    }
    ///@}
private:

    /// Helper function.
    inline BoundingBoxType GetBoundingBoxFromIndex(IndexType i, IndexType j, IndexType k, const PointType& rLowerBound, const PointType& rUpperBound) const {
        const PointType indices_d{ static_cast<double>(i), static_cast<double>(j), static_cast<double>(k) };
        PointType delta;
        delta[0] = std::abs(rUpperBound[0] - rLowerBound[0]) / static_cast<double>(mNumberOfElements[0]);
        delta[1] = std::abs(rUpperBound[1] - rLowerBound[1]) / static_cast<double>(mNumberOfElements[1]);
        delta[2] = std::abs(rUpperBound[2] - rLowerBound[2]) / static_cast<double>(mNumberOfElements[2]);
        return std::make_pair( Math::Add( rLowerBound, Math::MultElementWise(delta, indices_d)),
                               Math::Add( rLowerBound, Math::MultElementWise(delta, Math::Add({1.0, 1.0, 1.0}, indices_d))) );

    }

    ///@name Private Members
    ///@{

    BoundingBoxType mBoundXYZ;
    BoundingBoxType mBoundUVW;
    Vector3i mNumberOfElements;
    PartitionBoxType mGlobalPartition;
    bool mBSplineMesh;

    ///@}
}; // End class GridIndexer.
///@} End queso classes.


/// Output stream function
inline std::ostream& operator<<(std::ostream& os, GridIndexer::IndexInfo p) {
    switch(p) {
        case GridIndexer::IndexInfo::middle:
            return (os << "middle");
        case GridIndexer::IndexInfo::local_end:
            return (os << "local_end");
        case GridIndexer::IndexInfo::global_end:
            return (os << "global_end");
        default:
            return os;
    }
}

} // End namespace queso

#endif // GRID_INDEXER_INCLUDE_HPP
