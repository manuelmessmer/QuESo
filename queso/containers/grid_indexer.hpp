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

#ifndef MAPPING_UTILITIES_INCLUDE_H
#define MAPPING_UTILITIES_INCLUDE_H

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
        mBSplineMesh( rSettings[MainSettings::background_grid_settings].GetValue<BackgroundGridType>(BackgroundGridSettings::grid_type) ==  BackgroundGridType::b_spline_grid )
    {
    }

    ///@}
    ///@name Public Operations
    ///@{


    /// @brief Maps univariate index to trivariate indices i -> (i,j,k).
    /// @see GetVectorIndexFromMatrixIndices.
    /// @param Index
    /// @return Vector3i.
    Vector3i GetMatrixIndicesFromVectorIndex(const IndexType Index) const {
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
    IndexType GetVectorIndexFromMatrixIndices(const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex ) const {
        return DepthIndex * (mNumberOfElements[1]*mNumberOfElements[0]) + ColumnIndex * mNumberOfElements[0] + RowIndex;
    }

    /// @brief Maps trivariate indices to univariate index (i,j,k) -> i.
    /// @see GetMatrixIndicesFromVectorIndex.
    /// @param rIndices
    /// @return IndexType.
    IndexType GetVectorIndexFromMatrixIndices(const Vector3i& rIndices ) const {
        return rIndices[2] * (mNumberOfElements[1]*mNumberOfElements[0]) + rIndices[1] * mNumberOfElements[0] + rIndices[0];
    }

    /// @brief Creates bounding box in physical space from given index.
    /// @param Index
    /// @return BoundingBoxType.
    BoundingBoxType GetBoundingBoxXYZFromIndex(IndexType Index) const {
        const auto indices = GetMatrixIndicesFromVectorIndex(Index);
        return GetBoundingBoxFromIndex(indices[0], indices[1], indices[2], mBoundXYZ.first, mBoundXYZ.second);
    }

    /// @brief Creates bounding box in physical space from given indices.
    /// @param Indices
    /// @return BoundingBoxType.
    BoundingBoxType GetBoundingBoxXYZFromIndex(const Vector3i& rIndices) const {
        return GetBoundingBoxFromIndex(rIndices[0], rIndices[1], rIndices[2], mBoundXYZ.first, mBoundXYZ.second);
    }

    /// @brief Creates bounding box in physical space from given indices.
    /// @param i
    /// @param j
    /// @param k
    /// @return BoundingBoxType
    BoundingBoxType GetBoundingBoxXYZFromIndex(IndexType i, IndexType j, IndexType k) const {
        return GetBoundingBoxFromIndex(i, j, k, mBoundXYZ.first, mBoundXYZ.second);
    }

    /// @brief Creates bounding box in parametric space from given index.
    /// @param Index
    /// @return BoundingBoxType.
    BoundingBoxType GetBoundingBoxUVWFromIndex(IndexType Index) const {
        if( mBSplineMesh ) {
            const auto indices = GetMatrixIndicesFromVectorIndex(Index);
            return GetBoundingBoxFromIndex(indices[0], indices[1], indices[2], mBoundUVW.first, mBoundUVW.second);
        }
        return mBoundUVW;
    }

    /// @brief Creates bounding box in parametric space from given indices.
    /// @param Indices
    /// @return BoundingBoxType.
    BoundingBoxType GetBoundingBoxUVWFromIndex(const Vector3i& rIndices) const {
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
    BoundingBoxType GetBoundingBoxUVWFromIndex(IndexType i, IndexType j, IndexType k) const {
        if( mBSplineMesh ) {
            return GetBoundingBoxFromIndex(i, j, k, mBoundUVW.first, mBoundUVW.second);
        }
        return mBoundUVW;
    }

    /// @brief Returns global number of elements (including inactive elements).
    /// @return IndexType.
    IndexType NumberOfElements() const {
        return mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];
    }

    ///@}
private:

    BoundingBoxType GetBoundingBoxFromIndex(IndexType i, IndexType j, IndexType k, const PointType& rLowerBound, const PointType& rUpperBound) const {
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
    const bool mBSplineMesh;

    ///@}
}; // End class GridIndexer.
///@} End queso classes.

} // End namespace queso

#endif