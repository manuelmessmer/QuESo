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
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <ostream>
#include <utility>

//// Project includes
#include "queso/containers/dictionary.hpp"
#include "queso/containers/geometry_tolerance.hpp"
#include "queso/includes/define.hpp"

namespace queso {

/// @brief Dense identity of one physical background-grid face.
struct GridFaceId
{
    IndexType value;

    friend bool operator==(const GridFaceId&, const GridFaceId&) = default;
};

/// @brief Selects an adjacent cell by its position relative to a physical grid face.
/// @details For an internal face, the positive side is the canonical side under the grid's lower-inclusive,
///          upper-exclusive convention. An exterior face canonically belongs to its only adjacent cell, regardless
///          of side.
enum class GridFaceSide : std::uint8_t { negative, positive };

/// @brief Returns the opposite adjacent side of a physical grid face.
/// @param Side Side to reverse.
/// @return Opposite side.
[[nodiscard]] constexpr GridFaceSide Opposite(GridFaceSide Side) noexcept
{ return Side == GridFaceSide::negative ? GridFaceSide::positive : GridFaceSide::negative; }

/// @brief Writes a GridFaceSide value to an output stream.
inline std::ostream& operator<<(std::ostream& rStream, GridFaceSide Side)
{ return rStream << (Side == GridFaceSide::negative ? "negative" : "positive"); }

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
    /// @param rSettings Settings providing grid bounds, element counts, and grid type.
    explicit GridIndexer(const MainDictionaryType& rSettings) : GridIndexer(ValidateAndExtract(rSettings))
    {}

    ///@}
    ///@name Cell Indexing and Geometry
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

    /// @brief Returns the number of cells along the x, y, and z axes.
    /// @return Cell counts as {x, y, z}.
    [[nodiscard]] Vector3i ElementCounts() const noexcept
    { return { mNumberOfElementsX, mNumberOfElementsY, mNumberOfElementsZ }; }

    ///@}
    ///@name Face Topology
    ///@{

    /// @brief Returns the number of physical faces normal to the x, y, and z axes.
    /// @return Face counts as {x-normal, y-normal, z-normal}.
    [[nodiscard]] Vector3i FaceCounts() const noexcept
    {
        return { (mNumberOfElementsX + 1) * mNumberOfElementsY * mNumberOfElementsZ,
                 mNumberOfElementsX * (mNumberOfElementsY + 1) * mNumberOfElementsZ,
                 mNumberOfElementsXY * (mNumberOfElementsZ + 1) };
    }

    /// @brief Returns the total number of physical grid faces, including exterior faces.
    /// @return Number of dense face IDs.
    [[nodiscard]] IndexType NumberOfFaces() const noexcept
    {
        const Vector3i counts = FaceCounts();
        return counts[0] + counts[1] + counts[2];
    }

    /// @brief Resolves a cell-local direction to one stable physical grid-face identity.
    /// @details Opposite directions from neighboring cells resolve to the same face. The returned identity contains
    ///          no state, geometry, activity, or parent assignment.
    /// @param CellIndex Zero-based adjacent cell index.
    /// @param FaceDirection Cell-local face direction.
    /// @return Dense physical face identity.
    [[nodiscard]] GridFaceId GetFace(IndexType CellIndex, Direction FaceDirection) const noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(CellIndex < NumberOfElements(), "Cell index is out-of-bounds.");
        QuESo_ASSERT(FaceDirection != Direction::_end, "Direction::_end is not a physical grid face.");

        const Vector3i indices = GetMatrixIndicesFromVectorIndex(CellIndex);
        switch (FaceDirection) {
        case Direction::x_forward:
            return MakeFaceId(0, indices[0] + 1, indices[1], indices[2]);
        case Direction::x_backward:
            return MakeFaceId(0, indices[0], indices[1], indices[2]);
        case Direction::y_forward:
            return MakeFaceId(1, indices[0], indices[1] + 1, indices[2]);
        case Direction::y_backward:
            return MakeFaceId(1, indices[0], indices[1], indices[2]);
        case Direction::z_forward:
            return MakeFaceId(2, indices[0], indices[1], indices[2] + 1);
        case Direction::z_backward:
            return MakeFaceId(2, indices[0], indices[1], indices[2]);
        case Direction::_end:
            Unreachable("Direction::_end is not a physical grid face.");
        }
        Unreachable("Invalid grid-face direction.");
    }

    /// @brief Returns the cell adjacent to a physical face on the selected side.
    /// @details Internal faces have cells on both sides. Exterior faces return `std::nullopt` for the side outside
    ///          the grid. This operation performs topology navigation only and does not consider cell activity.
    /// @param Face Physical grid face.
    /// @param Side Negative or positive adjacent side.
    /// @return Adjacent cell index, or `std::nullopt` when that side lies outside the grid.
    [[nodiscard]] std::optional<IndexType> GetAdjacentCell(GridFaceId Face, GridFaceSide Side) const noexcept(NOTDEBUG)
    {
        const auto [axis, indices] = GetFaceIndices(Face);
        Vector3i cell_indices = indices;
        if (Side == GridFaceSide::negative) {
            if (indices[axis] == 0) { return std::nullopt; }
            --cell_indices[axis];
        } else {
            const Vector3i counts = ElementCounts();
            if (indices[axis] == counts[axis]) { return std::nullopt; }
        }
        return GetVectorIndexFromMatrixIndices(cell_indices);
    }

    /// @brief Returns which side of a physical face contains the supplied adjacent cell.
    /// @param Face Physical grid face.
    /// @param CellIndex Zero-based cell index adjacent to the face.
    /// @return Negative or positive adjacent side.
    [[nodiscard]] GridFaceSide GetAdjacentSide(GridFaceId Face, IndexType CellIndex) const noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(CellIndex < NumberOfElements(), "Cell index is out-of-bounds.");
        if (GetAdjacentCell(Face, GridFaceSide::negative) == CellIndex) { return GridFaceSide::negative; }
        QuESo_ASSERT(
            GetAdjacentCell(Face, GridFaceSide::positive) == CellIndex, "Cell is not adjacent to the physical face."
        );
        return GridFaceSide::positive;
    }

    /// @brief Returns the canonical adjacent cell of a physical grid face.
    /// @details An internal face canonically selects its positive-side cell. An exterior face canonically selects its
    ///          only adjacent cell, including a positive exterior face whose cell lies on the negative side. Cell
    ///          activity does not affect this fixed topology rule.
    /// @param Face Physical grid face.
    /// @return Canonical adjacent cell index.
    [[nodiscard]] IndexType GetCanonicalAdjacentCell(GridFaceId Face) const noexcept(NOTDEBUG)
    {
        if (const auto cell = GetAdjacentCell(Face, GridFaceSide::positive)) { return *cell; }
        const auto cell = GetAdjacentCell(Face, GridFaceSide::negative);
        QuESo_ASSERT(cell.has_value(), "A physical grid face must have an adjacent cell.");
        return *cell;
    }

    /// @brief Returns the coordinate axis normal to a physical grid face.
    /// @param Face Physical grid face.
    /// @return Axis in `[0, 3)`.
    [[nodiscard]] IndexType GetFaceAxis(GridFaceId Face) const noexcept(NOTDEBUG)
    { return GetFaceIndices(Face).first; }

    /// @brief Returns the exact physical plane coordinate of a grid face.
    /// @param Face Physical grid face.
    /// @return Plane coordinate along `GetFaceAxis(Face)`.
    [[nodiscard]] double GetFaceCoordinate(GridFaceId Face) const noexcept(NOTDEBUG)
    {
        const auto [axis, indices] = GetFaceIndices(Face);
        return mBoundXYZ.lower[axis] + mDeltaXYZ[axis] * static_cast<double>(indices[axis]);
    }

    ///@}
    ///@name Spatial Lookup and Bounds
    ///@{

    /// @brief Returns the unique half-open cell containing a point within the closed grid AABB.
    /// @details Points on internal planes select the positive-side cell. Points on a positive exterior plane select
    ///          the last cell because no positive-side cell exists.
    /// @param rPoint Point within the closed physical grid bounds.
    /// @return Zero-based containing cell index.
    [[nodiscard]] IndexType GetContainingCell(PointView rPoint) const noexcept(NOTDEBUG)
    {
        const Vector3i counts = ElementCounts();
        Vector3i indices{};
        for (IndexType axis = 0; axis < 3; ++axis) {
            QuESo_ASSERT(
                rPoint[axis] >= mBoundXYZ.lower[axis] && rPoint[axis] <= mBoundXYZ.upper[axis],
                "Point is outside the closed background-grid bounds."
            );
            const double position = SnapGridPosition((rPoint[axis] - mBoundXYZ.lower[axis]) / mDeltaXYZ[axis]);
            indices[axis] = std::min(static_cast<IndexType>(std::floor(position)), counts[axis] - 1);
        }
        return GetVectorIndexFromMatrixIndices(indices);
    }

    /// @brief Returns the physical bounds of the complete background grid.
    /// @return Global-space grid bounds.
    [[nodiscard]] const BoundingBoxType& BoundsXYZ() const noexcept
    { return mBoundXYZ; }

    /// @brief Returns the authoritative cell-scale geometry tolerance.
    /// @return Immutable geometry tolerance resolved from cell size and grid coordinate magnitude.
    [[nodiscard]] const GeometryTolerance& GetGeometryTolerance() const noexcept
    { return mGeometryTolerance; }

    ///@}
    ///@name Traversal
    ///@{

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
    ///@name Configuration
    ///@{

    /// @brief Validated values extracted from the background-grid settings.
    struct Configuration
    {
        BoundingBoxType bounds_xyz;  ///< Physical grid bounds.
        BoundingBoxType bounds_uvw;  ///< Parametric grid bounds.
        Vector3i number_of_elements;  ///< Cell counts along x, y, and z.
        bool is_b_spline_grid;  ///< Whether UVW bounds are subdivided per cell.
    };

    /// @brief Initializes the indexer from values extracted and validated by ValidateAndExtract().
    /// @param rConfiguration Validated grid configuration.
    explicit GridIndexer(const Configuration& rConfiguration)
        : mBoundXYZ(rConfiguration.bounds_xyz), mBoundUVW(rConfiguration.bounds_uvw),
          mNumberOfElementsX(rConfiguration.number_of_elements[0]),
          mNumberOfElementsY(rConfiguration.number_of_elements[1]),
          mNumberOfElementsZ(rConfiguration.number_of_elements[2]),
          mNumberOfElementsXY(mNumberOfElementsX * mNumberOfElementsY),
          mDeltaXYZ(Delta(mBoundXYZ.lower, mBoundXYZ.upper)), mDeltaUVW(Delta(mBoundUVW.lower, mBoundUVW.upper)),
          mGeometryTolerance(
              GeometryTolerance::FromScale(
                  { .length_scale = std::max({ mDeltaXYZ[0], mDeltaXYZ[1], mDeltaXYZ[2] }),
                    .coordinate_scale = std::max(
                        { 1.0,
                          std::abs(mBoundXYZ.lower[0]),
                          std::abs(mBoundXYZ.lower[1]),
                          std::abs(mBoundXYZ.lower[2]),
                          std::abs(mBoundXYZ.upper[0]),
                          std::abs(mBoundXYZ.upper[1]),
                          std::abs(mBoundXYZ.upper[2]) }
                    ) }
              )
          ),
          mGlobalPartition(
              std::make_pair(
                  Vector3i({ 0, 0, 0 }),
                  Vector3i({ mNumberOfElementsX - 1, mNumberOfElementsY - 1, mNumberOfElementsZ - 1 })
              )
          ),
          mIsBSplineGrid(rConfiguration.is_b_spline_grid)
    {
        const double tolerance = mGeometryTolerance.SnapDistance();
        QuESo_ERROR_IF(
            mDeltaXYZ[0] <= 2.0 * tolerance || mDeltaXYZ[1] <= 2.0 * tolerance || mDeltaXYZ[2] <= 2.0 * tolerance
        ) << "Background-grid cells must be larger than twice their geometry snap distance.\n";
    }

    /// @brief Extracts and validates the settings required to initialize the indexer.
    /// @param rSettings Background-grid settings.
    /// @return Validated configuration read exactly once from the settings.
    [[nodiscard]] static Configuration ValidateAndExtract(const MainDictionaryType& rSettings)
    {
        const auto& r_grid_settings = rSettings[MainSettings::background_grid_settings];
        const Vector3i counts = r_grid_settings.GetRequiredValue<Vector3i>(BackgroundGridSettings::number_of_elements);
        const PointType lower_xyz =
            r_grid_settings.GetRequiredValue<PointType>(BackgroundGridSettings::lower_bound_xyz);
        const PointType upper_xyz =
            r_grid_settings.GetRequiredValue<PointType>(BackgroundGridSettings::upper_bound_xyz);
        const PointType lower_uvw =
            r_grid_settings.GetRequiredValue<PointType>(BackgroundGridSettings::lower_bound_uvw);
        const PointType upper_uvw =
            r_grid_settings.GetRequiredValue<PointType>(BackgroundGridSettings::upper_bound_uvw);
        const GridType grid_type = r_grid_settings.GetRequiredValue<GridType>(BackgroundGridSettings::grid_type);
        for (IndexType axis = 0; axis < 3; ++axis) {
            QuESo_ERROR_IF(counts[axis] == 0)
                << "Background-grid element counts must be positive. Axis: " << axis << ".\n";
            QuESo_ERROR_IF(
                !std::isfinite(lower_xyz[axis]) || !std::isfinite(upper_xyz[axis]) || lower_xyz[axis] >= upper_xyz[axis]
            ) << "Background-grid XYZ bounds must be finite and strictly increasing. Axis: "
              << axis << ".\n";
        }
        return { .bounds_xyz = MakeBox(lower_xyz, upper_xyz),
                 .bounds_uvw = MakeBox(lower_uvw, upper_uvw),
                 .number_of_elements = counts,
                 .is_b_spline_grid = grid_type == GridType::b_spline_grid };
    }

    ///@}
    ///@name Numerical Helpers
    ///@{

    /// @brief Normalizes a dimensionless grid position affected only by floating-point roundoff.
    /// @details Values within a machine-epsilon-scale neighborhood of an integer are replaced by that integer. This is
    ///          not a geometric snapping tolerance.
    /// @param Position Position measured in cell widths from the grid lower bound.
    /// @return Position with near-integer roundoff removed.
    [[nodiscard]] static double SnapGridPosition(double Position) noexcept
    {
        const double nearest = std::round(Position);
        const double tolerance = 16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(Position));
        return std::abs(Position - nearest) <= tolerance ? nearest : Position;
    }

    ///@}
    ///@name Face Indexing Helpers
    ///@{

    /// @brief Creates a dense physical-face ID from its axis-specific integer coordinates.
    /// @param Axis Face-normal axis in `[0, 3)`.
    /// @param X Face-grid x coordinate.
    /// @param Y Face-grid y coordinate.
    /// @param Z Face-grid z coordinate.
    /// @return Dense physical-face identity.
    [[nodiscard]] GridFaceId MakeFaceId(IndexType Axis, IndexType X, IndexType Y, IndexType Z) const noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(Axis < 3, "Grid-face axis is out-of-bounds.");
        const Vector3i face_counts = FaceCounts();
        if (Axis == 0) {
            return { Z * (mNumberOfElementsX + 1) * mNumberOfElementsY + Y * (mNumberOfElementsX + 1) + X };
        }
        if (Axis == 1) {
            return { face_counts[0] + Z * mNumberOfElementsX * (mNumberOfElementsY + 1) + Y * mNumberOfElementsX + X };
        }
        return { face_counts[0] + face_counts[1] + Z * mNumberOfElementsXY + Y * mNumberOfElementsX + X };
    }

    /// @brief Decodes a dense physical-face ID into its normal axis and axis-specific integer coordinates.
    /// @param Face Physical grid-face identity.
    /// @return Pair containing the normal axis and face-grid coordinates.
    [[nodiscard]] std::pair<IndexType, Vector3i> GetFaceIndices(GridFaceId Face) const noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(Face.value < NumberOfFaces(), "Grid-face ID is out-of-bounds.");
        const Vector3i face_counts = FaceCounts();
        if (Face.value < face_counts[0]) {
            const IndexType plane_size = (mNumberOfElementsX + 1) * mNumberOfElementsY;
            const IndexType z = Face.value / plane_size;
            const IndexType in_plane = Face.value - z * plane_size;
            return { 0, { in_plane % (mNumberOfElementsX + 1), in_plane / (mNumberOfElementsX + 1), z } };
        }
        const IndexType y_offset = Face.value - face_counts[0];
        if (y_offset < face_counts[1]) {
            const IndexType plane_size = mNumberOfElementsX * (mNumberOfElementsY + 1);
            const IndexType z = y_offset / plane_size;
            const IndexType in_plane = y_offset - z * plane_size;
            return { 1, { in_plane % mNumberOfElementsX, in_plane / mNumberOfElementsX, z } };
        }
        const IndexType z_offset = y_offset - face_counts[1];
        const IndexType z = z_offset / mNumberOfElementsXY;
        const IndexType in_plane = z_offset - z * mNumberOfElementsXY;
        return { 2, { in_plane % mNumberOfElementsX, in_plane / mNumberOfElementsX, z } };
    }

    ///@}
    ///@name Traversal Helpers
    ///@{

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

    ///@}
    ///@name Cell Geometry Helpers
    ///@{

    /// @brief Creates a cell bounding box from zero-based matrix indices and a precomputed grid spacing.
    /// @param RowIndex Zero-based x cell index.
    /// @param ColumnIndex Zero-based y cell index.
    /// @param DepthIndex Zero-based z cell index.
    /// @param rLowerBound Lower bound of the complete grid.
    /// @param rDelta Cell size along each axis.
    /// @return Cell bounding box.
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
    /// @param rLowerBound Lower grid bound.
    /// @param rUpperBound Upper grid bound.
    /// @return Cell size along x, y, and z.
    [[nodiscard]] PointType Delta(const PointType& rLowerBound, const PointType& rUpperBound) const noexcept
    {
        return PointType{ std::abs(rUpperBound[0] - rLowerBound[0]) / static_cast<double>(mNumberOfElementsX),
                          std::abs(rUpperBound[1] - rLowerBound[1]) / static_cast<double>(mNumberOfElementsY),
                          std::abs(rUpperBound[2] - rLowerBound[2]) / static_cast<double>(mNumberOfElementsZ) };
    }

    ///@}
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

    GeometryTolerance mGeometryTolerance;

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
