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

#pragma once

//// STL includes
#include <array>
#include <cmath>
#include <ostream>
#include <utility>

//// Project includes
#include "queso/containers/dictionary.hpp"
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo classes
///@{

/// @class  GridIndexer
/// @author Manuel Messmer
/// @brief  Provides fast zero-based indexing, traversal, and cell-bound queries for the background grid.
class GridIndexer
{
public:
    ///@name Type Definitions
    ///@{
    enum class IndexInfo { middle, local_end, global_end, _end };
    enum class Direction { x_forward, x_backward, y_forward, y_backward, z_forward, z_backward, _end };

    using IndexReturnType = std::pair<IndexType, IndexInfo>;
    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructs the indexer from the background-grid settings.
    /// @param rSettings
    explicit GridIndexer(const MainDictionaryType& rSettings)
        : mBoundXYZ(MakeBox(
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                  BackgroundGridSettings::lower_bound_xyz
              ),
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                  BackgroundGridSettings::upper_bound_xyz
              )
          )),
          mBoundUVW(MakeBox(
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                  BackgroundGridSettings::lower_bound_uvw
              ),
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<PointType>(
                  BackgroundGridSettings::upper_bound_uvw
              )
          )),
          mNumberOfElementsX(
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<Vector3i>(
                  BackgroundGridSettings::number_of_elements
              )[0]
          ),
          mNumberOfElementsY(
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<Vector3i>(
                  BackgroundGridSettings::number_of_elements
              )[1]
          ),
          mNumberOfElementsZ(
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<Vector3i>(
                  BackgroundGridSettings::number_of_elements
              )[2]
          ),
          mNumberOfElementsXY(mNumberOfElementsX * mNumberOfElementsY),
          mDeltaXYZ(Delta(mBoundXYZ.lower, mBoundXYZ.upper)), mDeltaUVW(Delta(mBoundUVW.lower, mBoundUVW.upper)),
          mGlobalPartition(
              std::make_pair(
                  Vector3i({ 0, 0, 0 }),
                  Vector3i({ mNumberOfElementsX - 1, mNumberOfElementsY - 1, mNumberOfElementsZ - 1 })
              )
          ),
          mIsBSplineGrid(
              rSettings[MainSettings::background_grid_settings].GetRequiredValue<GridType>(
                  BackgroundGridSettings::grid_type
              )
              == GridType::b_spline_grid
          )
    {}

    ///@}
    ///@name Public Operations
    ///@{

    /// @brief Maps a zero-based linear index to zero-based matrix indices.
    /// @see GetVectorIndexFromMatrixIndices.
    /// @param Index Linear index.
    /// @return Matrix indices as {row, column, depth}, corresponding to the {x, y, z} cell indices.
    [[nodiscard]] Vector3i GetMatrixIndicesFromVectorIndex(IndexType Index) const noexcept
    {
        const IndexType depth_index = Index / mNumberOfElementsXY;
        const IndexType index_in_row_column_plane = Index - depth_index * mNumberOfElementsXY;
        Vector3i result{ index_in_row_column_plane % mNumberOfElementsX,  // row
                         index_in_row_column_plane / mNumberOfElementsX,  // column
                         depth_index };  // depth

        return result;
    }

    /// @brief Maps zero-based matrix indices to a zero-based linear index.
    /// @see GetMatrixIndicesFromVectorIndex.
    /// @param RowIndex Zero-based row index, corresponding to the x-direction cell index.
    /// @param ColumnIndex Zero-based column index, corresponding to the y-direction cell index.
    /// @param DepthIndex Zero-based depth index, corresponding to the z-direction cell index.
    /// @return IndexType.
    [[nodiscard]] IndexType
        GetVectorIndexFromMatrixIndices(IndexType RowIndex, IndexType ColumnIndex, IndexType DepthIndex) const noexcept
    { return DepthIndex * mNumberOfElementsXY + ColumnIndex * mNumberOfElementsX + RowIndex; }

    /// @brief Maps zero-based matrix indices to a zero-based linear index.
    /// @see GetMatrixIndicesFromVectorIndex.
    /// @param rIndices Matrix indices as {row, column, depth}, corresponding to the {x, y, z} cell indices.
    /// @return IndexType.
    [[nodiscard]] IndexType GetVectorIndexFromMatrixIndices(const Vector3i& rIndices) const noexcept
    { return rIndices[2] * mNumberOfElementsXY + rIndices[1] * mNumberOfElementsX + rIndices[0]; }

    /// @brief Creates the physical-space cell bounding box for a zero-based linear index.
    /// @param Index Linear index.
    /// @return BoundingBoxType.
    [[nodiscard]] BoundingBoxType GetBoundingBoxXYZFromIndex(IndexType Index) const noexcept
    {
        const auto indices = GetMatrixIndicesFromVectorIndex(Index);
        return GetBoundingBoxFromIndex(indices[0], indices[1], indices[2], mBoundXYZ.lower, mDeltaXYZ);
    }

    /// @brief Creates the physical-space cell bounding box for zero-based matrix indices.
    /// @param rIndices Matrix indices as {row, column, depth}, corresponding to the {x, y, z} cell indices.
    /// @return BoundingBoxType.
    [[nodiscard]] BoundingBoxType GetBoundingBoxXYZFromIndex(const Vector3i& rIndices) const noexcept
    { return GetBoundingBoxFromIndex(rIndices[0], rIndices[1], rIndices[2], mBoundXYZ.lower, mDeltaXYZ); }

    /// @brief Creates the physical-space cell bounding box for zero-based matrix indices.
    /// @param RowIndex Zero-based row index, corresponding to the x-direction cell index.
    /// @param ColumnIndex Zero-based column index, corresponding to the y-direction cell index.
    /// @param DepthIndex Zero-based depth index, corresponding to the z-direction cell index.
    /// @return BoundingBoxType
    [[nodiscard]] BoundingBoxType
        GetBoundingBoxXYZFromIndex(IndexType RowIndex, IndexType ColumnIndex, IndexType DepthIndex) const noexcept
    { return GetBoundingBoxFromIndex(RowIndex, ColumnIndex, DepthIndex, mBoundXYZ.lower, mDeltaXYZ); }

    /// @brief Creates the parametric-space bounding box for a zero-based linear index.
    /// @details B-spline grids return a per-cell UVW box. FE grids return the full stored UVW bound.
    /// @param Index Linear index.
    /// @return BoundingBoxType.
    [[nodiscard]] BoundingBoxType GetBoundingBoxUVWFromIndex(IndexType Index) const noexcept
    {
        if (mIsBSplineGrid) {
            const auto indices = GetMatrixIndicesFromVectorIndex(Index);
            return GetBoundingBoxFromIndex(indices[0], indices[1], indices[2], mBoundUVW.lower, mDeltaUVW);
        }
        return mBoundUVW;
    }

    /// @brief Creates the parametric-space bounding box for zero-based matrix indices.
    /// @details B-spline grids return a per-cell UVW box. FE grids return the full stored UVW bound.
    /// @param rIndices Matrix indices as {row, column, depth}, corresponding to the {x, y, z} cell indices.
    /// @return BoundingBoxType.
    [[nodiscard]] BoundingBoxType GetBoundingBoxUVWFromIndex(const Vector3i& rIndices) const noexcept
    {
        if (mIsBSplineGrid) {
            return GetBoundingBoxFromIndex(rIndices[0], rIndices[1], rIndices[2], mBoundUVW.lower, mDeltaUVW);
        }
        return mBoundUVW;
    }

    /// @brief Creates the parametric-space bounding box for zero-based matrix indices.
    /// @details B-spline grids return a per-cell UVW box. FE grids return the full stored UVW bound.
    /// @param RowIndex Zero-based row index, corresponding to the x-direction cell index.
    /// @param ColumnIndex Zero-based column index, corresponding to the y-direction cell index.
    /// @param DepthIndex Zero-based depth index, corresponding to the z-direction cell index.
    /// @return BoundingBoxType
    [[nodiscard]] BoundingBoxType
        GetBoundingBoxUVWFromIndex(IndexType RowIndex, IndexType ColumnIndex, IndexType DepthIndex) const noexcept
    {
        if (mIsBSplineGrid) {
            return GetBoundingBoxFromIndex(RowIndex, ColumnIndex, DepthIndex, mBoundUVW.lower, mDeltaUVW);
        }
        return mBoundUVW;
    }

    /// @brief Returns the global number of grid elements, including inactive elements.
    /// @return IndexType.
    [[nodiscard]] IndexType NumberOfElements() const noexcept
    { return mNumberOfElementsX * mNumberOfElementsY * mNumberOfElementsZ; }

    /// @brief Returns the next linear index in the given traversal direction.
    /// @param Index Current zero-based linear index.
    /// @param Dir Traversal direction.
    /// @return Pair of next index and status. The status describes whether the current index is middle, local end,
    ///         or global end.
    [[nodiscard]] IndexReturnType GetNextIndex(IndexType Index, Direction Dir) const noexcept(NOTDEBUG)
    { return GetNextIndex(Index, Dir, mGlobalPartition); }

    /// @brief Returns the next linear index in the given traversal direction within a partition.
    /// @param Index Current zero-based linear index.
    /// @param Dir Traversal direction.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return Pair of next index and status. The status describes whether the current index is middle, local end,
    ///         or global end.
    [[nodiscard]] IndexReturnType
        GetNextIndex(IndexType Index, Direction Dir, const PartitionBoxType& rPartition) const noexcept(NOTDEBUG)
    {
        switch (Dir) {
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
        case Direction::_end:
            Unreachable("Direction::_end is a sentinel, not a runtime value.");
        }
        Unreachable("Invalid Direction value.");
    }

    /// @brief Returns the next linear index in the compile-time traversal direction.
    /// @tparam TDir Traversal direction.
    /// @param Index Current zero-based linear index.
    /// @return Pair of next index and status for the current index.
    template<Direction TDir>
    [[nodiscard]] IndexReturnType GetNextIndex(IndexType Index) const noexcept
    { return GetNextIndex<TDir>(Index, mGlobalPartition); }

    /// @brief Returns the next linear index in the compile-time traversal direction within a partition.
    /// @tparam TDir Traversal direction.
    /// @param Index Current zero-based linear index.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return Pair of next index and status for the current index.
    template<Direction TDir>
    [[nodiscard]] IndexReturnType GetNextIndex(IndexType Index, const PartitionBoxType& rPartition) const noexcept
    {
        if constexpr (TDir == Direction::x_forward) {
            return GetNextIndexX(Index, rPartition);
        } else if constexpr (TDir == Direction::x_backward) {
            return GetPreviousIndexX(Index, rPartition);
        } else if constexpr (TDir == Direction::y_forward) {
            return GetNextIndexY(Index, rPartition);
        } else if constexpr (TDir == Direction::y_backward) {
            return GetPreviousIndexY(Index, rPartition);
        } else if constexpr (TDir == Direction::z_forward) {
            return GetNextIndexZ(Index, rPartition);
        } else if constexpr (TDir == Direction::z_backward) {
            return GetPreviousIndexZ(Index, rPartition);
        } else {
            static_assert(always_false_v<TDir>, "Unsupported traversal direction.");
        }
    }

    /// @brief Returns whether the current index lies on the global end plane for the given direction.
    /// @param i Current zero-based linear index.
    /// @param Dir Traversal direction.
    /// @return bool
    [[nodiscard]] bool IsEnd(IndexType i, Direction Dir) const noexcept(NOTDEBUG)
    { return IsEnd(i, Dir, mGlobalPartition); }

    /// @brief Returns whether the current index lies on the partition end plane for the given direction.
    /// @param i Current zero-based linear index.
    /// @param Dir Traversal direction.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return bool
    [[nodiscard]] bool IsEnd(IndexType i, Direction Dir, const PartitionBoxType& rPartition) const noexcept(NOTDEBUG)
    {
        switch (Dir) {
        case Direction::x_forward:
            return IsEndXForward(i, rPartition);
        case Direction::x_backward:
            return IsEndXBackward(i, rPartition);
        case Direction::y_forward:
            return IsEndYForward(i, rPartition);
        case Direction::y_backward:
            return IsEndYBackward(i, rPartition);
        case Direction::z_forward:
            return IsEndZForward(i, rPartition);
        case Direction::z_backward:
            return IsEndZBackward(i, rPartition);
        case Direction::_end:
            Unreachable("Direction::_end is a sentinel, not a runtime value.");
        }
        Unreachable("Invalid Direction value.");
    }

    /// @brief Returns whether the current index lies on the global end plane for the compile-time direction.
    /// @tparam TDir Traversal direction.
    /// @param i Current zero-based linear index.
    /// @return bool
    template<Direction TDir>
    [[nodiscard]] bool IsEnd(IndexType i) const noexcept
    { return IsEnd<TDir>(i, mGlobalPartition); }

    /// @brief Returns whether the current index lies on the partition end plane for the compile-time direction.
    /// @tparam TDir Traversal direction.
    /// @param i Current zero-based linear index.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return bool
    template<Direction TDir>
    [[nodiscard]] bool IsEnd(IndexType i, const PartitionBoxType& rPartition) const noexcept
    {
        if constexpr (TDir == Direction::x_forward) {
            return IsEndXForward(i, rPartition);
        } else if constexpr (TDir == Direction::x_backward) {
            return IsEndXBackward(i, rPartition);
        } else if constexpr (TDir == Direction::y_forward) {
            return IsEndYForward(i, rPartition);
        } else if constexpr (TDir == Direction::y_backward) {
            return IsEndYBackward(i, rPartition);
        } else if constexpr (TDir == Direction::z_forward) {
            return IsEndZForward(i, rPartition);
        } else if constexpr (TDir == Direction::z_backward) {
            return IsEndZBackward(i, rPartition);
        } else {
            static_assert(always_false_v<TDir>, "Unsupported traversal direction.");
        }
    }

    /// @brief Returns the opposite traversal direction.
    /// @param Dir Traversal direction. Must not be Direction::_end.
    /// @return Direction.
    [[nodiscard]] static constexpr Direction ReverseDirection(Direction Dir) noexcept
    {
        QuESo_ASSERT(Dir != Direction::_end, "Direction::_end is a sentinel, not a runtime value.");
        constexpr std::array<Direction, 6> map_direction = { Direction::x_backward, Direction::x_forward,
                                                             Direction::y_backward, Direction::y_forward,
                                                             Direction::z_backward, Direction::z_forward };
        return map_direction[static_cast<IndexType>(Dir)];
    }
    ///@}

private:
    /// @brief Returns the next index in x-forward traversal order within a partition.
    /// @param i Current zero-based linear index.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return Pair of next index and status for the current index.
    [[nodiscard]] IndexReturnType GetNextIndexX(IndexType i, const PartitionBoxType& rPartition) const noexcept
    {
        const IndexType row_index = i % mNumberOfElementsX;
        if (row_index < rPartition.second[0]) { return { i + 1, IndexInfo::middle }; }

        const IndexType column_index = (i / mNumberOfElementsX) % mNumberOfElementsY;
        if (column_index < rPartition.second[1]) {
            const IndexType depth_index = i / mNumberOfElementsXY;
            return { GetVectorIndexFromMatrixIndices(rPartition.first[0], column_index + 1, depth_index),
                     IndexInfo::local_end };
        }

        const IndexType depth_index = i / mNumberOfElementsXY;
        if (depth_index < rPartition.second[2]) {
            return { GetVectorIndexFromMatrixIndices(rPartition.first[0], rPartition.first[1], depth_index + 1),
                     IndexInfo::local_end };
        }

        return { i, IndexInfo::global_end };
    }

    /// @brief Returns the next index in y-forward traversal order within a partition.
    /// @param i Current zero-based linear index.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return Pair of next index and status for the current index.
    [[nodiscard]] IndexReturnType GetNextIndexY(IndexType i, const PartitionBoxType& rPartition) const noexcept
    {
        const IndexType column_index = (i / mNumberOfElementsX) % mNumberOfElementsY;
        if (column_index < rPartition.second[1]) { return { i + mNumberOfElementsX, IndexInfo::middle }; }

        const IndexType row_index = i % mNumberOfElementsX;
        if (row_index < rPartition.second[0]) {
            const IndexType depth_index = i / mNumberOfElementsXY;
            return { GetVectorIndexFromMatrixIndices(row_index + 1, rPartition.first[1], depth_index),
                     IndexInfo::local_end };
        }

        const IndexType depth_index = i / mNumberOfElementsXY;
        if (depth_index < rPartition.second[2]) {
            return { GetVectorIndexFromMatrixIndices(rPartition.first[0], rPartition.first[1], depth_index + 1),
                     IndexInfo::local_end };
        }

        return { i, IndexInfo::global_end };
    }

    /// @brief Returns the next index in z-forward traversal order within a partition.
    /// @param i Current zero-based linear index.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return Pair of next index and status for the current index.
    [[nodiscard]] IndexReturnType GetNextIndexZ(IndexType i, const PartitionBoxType& rPartition) const noexcept
    {
        const IndexType depth_index = i / mNumberOfElementsXY;
        if (depth_index < rPartition.second[2]) { return { i + mNumberOfElementsXY, IndexInfo::middle }; }

        const IndexType row_index = i % mNumberOfElementsX;
        if (row_index < rPartition.second[0]) {
            const IndexType column_index = (i / mNumberOfElementsX) % mNumberOfElementsY;
            return { GetVectorIndexFromMatrixIndices(row_index + 1, column_index, rPartition.first[2]),
                     IndexInfo::local_end };
        }

        const IndexType column_index = (i / mNumberOfElementsX) % mNumberOfElementsY;
        if (column_index < rPartition.second[1]) {
            return { GetVectorIndexFromMatrixIndices(rPartition.first[0], column_index + 1, rPartition.first[2]),
                     IndexInfo::local_end };
        }

        return { i, IndexInfo::global_end };
    }


    /// @brief Returns the next index in x-backward traversal order within a partition.
    /// @param i Current zero-based linear index.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return Pair of next index and status for the current index.
    [[nodiscard]] IndexReturnType GetPreviousIndexX(IndexType i, const PartitionBoxType& rPartition) const noexcept
    {
        const IndexType row_index = i % mNumberOfElementsX;
        if (row_index > rPartition.first[0]) { return { i - 1, IndexInfo::middle }; }

        const IndexType column_index = (i / mNumberOfElementsX) % mNumberOfElementsY;
        if (column_index > rPartition.first[1]) {
            const IndexType depth_index = i / mNumberOfElementsXY;
            return { GetVectorIndexFromMatrixIndices(rPartition.second[0], column_index - 1, depth_index),
                     IndexInfo::local_end };
        }

        const IndexType depth_index = i / mNumberOfElementsXY;
        if (depth_index > rPartition.first[2]) {
            return { GetVectorIndexFromMatrixIndices(rPartition.second[0], rPartition.second[1], depth_index - 1),
                     IndexInfo::local_end };
        }

        return { i, IndexInfo::global_end };
    }

    /// @brief Returns the next index in y-backward traversal order within a partition.
    /// @param i Current zero-based linear index.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return Pair of next index and status for the current index.
    [[nodiscard]] IndexReturnType GetPreviousIndexY(IndexType i, const PartitionBoxType& rPartition) const noexcept
    {
        const IndexType column_index = (i / mNumberOfElementsX) % mNumberOfElementsY;
        if (column_index > rPartition.first[1]) { return { i - mNumberOfElementsX, IndexInfo::middle }; }

        const IndexType row_index = i % mNumberOfElementsX;
        if (row_index > rPartition.first[0]) {
            const IndexType depth_index = i / mNumberOfElementsXY;
            return { GetVectorIndexFromMatrixIndices(row_index - 1, rPartition.second[1], depth_index),
                     IndexInfo::local_end };
        }

        const IndexType depth_index = i / mNumberOfElementsXY;
        if (depth_index > rPartition.first[2]) {
            return { GetVectorIndexFromMatrixIndices(rPartition.second[0], rPartition.second[1], depth_index - 1),
                     IndexInfo::local_end };
        }

        return { i, IndexInfo::global_end };
    }

    /// @brief Returns the next index in z-backward traversal order within a partition.
    /// @param i Current zero-based linear index.
    /// @param rPartition Inclusive zero-based partition bounds.
    /// @return Pair of next index and status for the current index.
    [[nodiscard]] IndexReturnType GetPreviousIndexZ(IndexType i, const PartitionBoxType& rPartition) const noexcept
    {
        const IndexType depth_index = i / mNumberOfElementsXY;
        if (depth_index > rPartition.first[2]) { return { i - mNumberOfElementsXY, IndexInfo::middle }; }

        const IndexType row_index = i % mNumberOfElementsX;
        if (row_index > rPartition.first[0]) {
            const IndexType column_index = (i / mNumberOfElementsX) % mNumberOfElementsY;
            return { GetVectorIndexFromMatrixIndices(row_index - 1, column_index, rPartition.second[2]),
                     IndexInfo::local_end };
        }

        const IndexType column_index = (i / mNumberOfElementsX) % mNumberOfElementsY;
        if (column_index > rPartition.first[1]) {
            return { GetVectorIndexFromMatrixIndices(rPartition.second[0], column_index - 1, rPartition.second[2]),
                     IndexInfo::local_end };
        }

        return { i, IndexInfo::global_end };
    }

    /// @brief Returns whether the current index lies on the x-forward end plane.
    [[nodiscard]] bool IsEndXForward(IndexType i, const PartitionBoxType& rPartition) const noexcept
    { return (i % mNumberOfElementsX) == rPartition.second[0]; }

    /// @brief Returns whether the current index lies on the x-backward end plane.
    [[nodiscard]] bool IsEndXBackward(IndexType i, const PartitionBoxType& rPartition) const noexcept
    { return (i % mNumberOfElementsX) == rPartition.first[0]; }

    /// @brief Returns whether the current index lies on the y-forward end plane.
    [[nodiscard]] bool IsEndYForward(IndexType i, const PartitionBoxType& rPartition) const noexcept
    { return ((i / mNumberOfElementsX) % mNumberOfElementsY) == rPartition.second[1]; }

    /// @brief Returns whether the current index lies on the y-backward end plane.
    [[nodiscard]] bool IsEndYBackward(IndexType i, const PartitionBoxType& rPartition) const noexcept
    { return ((i / mNumberOfElementsX) % mNumberOfElementsY) == rPartition.first[1]; }

    /// @brief Returns whether the current index lies on the z-forward end plane.
    [[nodiscard]] bool IsEndZForward(IndexType i, const PartitionBoxType& rPartition) const noexcept
    { return (i / mNumberOfElementsXY) == rPartition.second[2]; }

    /// @brief Returns whether the current index lies on the z-backward end plane.
    [[nodiscard]] bool IsEndZBackward(IndexType i, const PartitionBoxType& rPartition) const noexcept
    { return (i / mNumberOfElementsXY) == rPartition.first[2]; }

    /// @brief Creates a cell bounding box from zero-based matrix indices and a precomputed grid spacing.
    [[nodiscard]] BoundingBoxType GetBoundingBoxFromIndex(
        IndexType RowIndex,
        IndexType ColumnIndex,
        IndexType DepthIndex,
        const PointType& rLowerBound,
        const PointType& rDelta
    ) const noexcept
    {
        const PointType lower{ rLowerBound[0] + rDelta[0] * static_cast<double>(RowIndex),
                               rLowerBound[1] + rDelta[1] * static_cast<double>(ColumnIndex),
                               rLowerBound[2] + rDelta[2] * static_cast<double>(DepthIndex) };
        const PointType upper{ rLowerBound[0] + rDelta[0] * static_cast<double>(RowIndex + 1),
                               rLowerBound[1] + rDelta[1] * static_cast<double>(ColumnIndex + 1),
                               rLowerBound[2] + rDelta[2] * static_cast<double>(DepthIndex + 1) };
        return MakeBox(lower, upper);
    }

    /// @brief Calculates the grid spacing for the configured number of elements.
    [[nodiscard]] PointType Delta(const PointType& rLowerBound, const PointType& rUpperBound) const noexcept
    {
        return PointType{ std::abs(rUpperBound[0] - rLowerBound[0]) / static_cast<double>(mNumberOfElementsX),
                          std::abs(rUpperBound[1] - rLowerBound[1]) / static_cast<double>(mNumberOfElementsY),
                          std::abs(rUpperBound[2] - rLowerBound[2]) / static_cast<double>(mNumberOfElementsZ) };
    }

    ///@name Private Members
    ///@{

    BoundingBoxType mBoundXYZ;
    BoundingBoxType mBoundUVW;

    IndexType mNumberOfElementsX;
    IndexType mNumberOfElementsY;
    IndexType mNumberOfElementsZ;
    IndexType mNumberOfElementsXY;

    PointType mDeltaXYZ;
    PointType mDeltaUVW;

    PartitionBoxType mGlobalPartition;

    bool mIsBSplineGrid;
    ///@}
};  // class GridIndexer
///@}


/// @brief Writes a GridIndexer::IndexInfo value to an output stream.
inline std::ostream& operator<<(std::ostream& os, GridIndexer::IndexInfo p)
{
    switch (p) {
    case GridIndexer::IndexInfo::middle:
        return (os << "middle");
    case GridIndexer::IndexInfo::local_end:
        return (os << "local_end");
    case GridIndexer::IndexInfo::global_end:
        return (os << "global_end");
    case GridIndexer::IndexInfo::_end:
        Unreachable("IndexInfo::_end is a sentinel, not a runtime value.");
    default:
        return os;
    }
}

}  // namespace queso
