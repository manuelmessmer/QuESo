// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// STL includes
#include <cstdlib>

// Project includes
#include "utilities/mapping_utilities.h"

namespace tibra {

// Static member operations for Mapping
PointType Mapping::PointFromGlobalToParam( const PointType& rGlobalCoord, const PointType& rLowerPoint, const PointType& rUpperPoint){
    return PointType( (rGlobalCoord[0] - rLowerPoint[0]) / std::abs(rUpperPoint[0] - rLowerPoint[0]),
                      (rGlobalCoord[1] - rLowerPoint[1]) / std::abs(rUpperPoint[1] - rLowerPoint[1]),
                      (rGlobalCoord[2] - rLowerPoint[2]) / std::abs(rUpperPoint[2] - rLowerPoint[2]) );
}

PointType Mapping::PointFromParamToGlobal( const PointType& rLocalCoord, const PointType& rLowerPoint, const PointType& rUpperPoint ) {
    return PointType( (rLocalCoord[0] * std::abs(rUpperPoint[0] - rLowerPoint[0]) ) + rLowerPoint[0],
                      (rLocalCoord[1] * std::abs(rUpperPoint[1] - rLowerPoint[1]) ) + rLowerPoint[1],
                      (rLocalCoord[2] * std::abs(rUpperPoint[2] - rLowerPoint[2]) ) + rLowerPoint[2] );
}

Vector3i Mapping::GetMatrixIndicesFromVectorIndex(const IndexType Index, const Vector3i& rNumberOfElements) {
    Vector3i result;
    const IndexType index_in_row_column_plane = Index % (rNumberOfElements[0]*rNumberOfElements[1]);
    result[0] = index_in_row_column_plane % rNumberOfElements[0]; // row
    result[1] = index_in_row_column_plane / rNumberOfElements[0]; // column
    result[2] = Index / (rNumberOfElements[0]*rNumberOfElements[1]);   // depth

    return result;
}

IndexType Mapping::GetVectorIndexFromMatrixIndices(const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex, const Vector3i& rNumberOfElements) {
    return DepthIndex * (rNumberOfElements[1]*rNumberOfElements[0]) + ColumnIndex * rNumberOfElements[0] + RowIndex;
}

IndexType Mapping::GetVectorIndexFromMatrixIndices(const Vector3i& rIndices, const Vector3i& rNumberOfElements) {
    return rIndices[2] * (rNumberOfElements[1]*rNumberOfElements[0]) + rIndices[1] * rNumberOfElements[0] + rIndices[0];
}
std::pair<PointType, PointType> Mapping::GetBoundingBoxFromIndex(IndexType Index, const PointType& rLowerBound, const PointType& rUpperBound, const Vector3i& rNumberOfElements)  {
    const auto indices = GetMatrixIndicesFromVectorIndex(Index, rNumberOfElements);
    return GetBoundingBoxFromIndex( indices[0], indices[1], indices[2], rLowerBound, rUpperBound, rNumberOfElements);
}

std::pair<PointType, PointType> Mapping::GetBoundingBoxFromIndex(const Vector3i& rIndices, const PointType& rLowerBound, const PointType& rUpperBound, const Vector3i& rNumberOfElements )  {
    return GetBoundingBoxFromIndex( rIndices[0], rIndices[1], rIndices[2], rLowerBound, rUpperBound, rNumberOfElements);
}

std::pair<PointType, PointType> Mapping::GetBoundingBoxFromIndex(IndexType i, IndexType j, IndexType k, const PointType& rLowerBound, const PointType& rUpperBound, const Vector3i& rNumberOfElements)  {
    const PointType indices_d( static_cast<double>(i), static_cast<double>(j), static_cast<double(k) );
    PointType delta;
    delta[0] = std::abs(rUpperBound[0] - rLowerBound[0]) / (rNumberOfElements[0]);
    delta[1] = std::abs(rUpperBound[1] - rLowerBound[1]) / (rNumberOfElements[1]);
    delta[2] = std::abs(rUpperBound[2] - rLowerBound[2]) / (rNumberOfElements[2]);
    return std::make_pair( rLowerBound + indices_d * delta,
                            rLowerBound + (indices_d+1.0) * delta );
}

// Member operations for Mapper
PointType Mapper::PointFromGlobalToParam( const PointType& rGlobalCoord ) const {
    return Mapping::PointFromGlobalToParam(rGlobalCoord, mLowerBound, mUpperBound);
}

PointType Mapper::PointFromParamToGlobal( const PointType& rLocalCoord ) const {
    return Mapping::PointFromParamToGlobal(rLocalCoord, mLowerBound, mUpperBound);
}

Vector3i Mapper::GetMatrixIndicesFromVectorIndex(const IndexType Index ) const {
    return Mapping::GetMatrixIndicesFromVectorIndex(Index, mNumberOfElements);
}

IndexType Mapper::GetVectorIndexFromMatrixIndices(const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex ) const {
    return Mapping::GetVectorIndexFromMatrixIndices(RowIndex, ColumnIndex, DepthIndex, mNumberOfElements );
}

IndexType Mapper::GetVectorIndexFromMatrixIndices(const Vector3i& rIndices ) const {
    return Mapping::GetVectorIndexFromMatrixIndices(rIndices, mNumberOfElements );
}

std::pair<PointType, PointType> Mapper::GetBoundingBoxFromIndex(IndexType Index) const {
    return Mapping::GetBoundingBoxFromIndex(Index, mLowerBound, mUpperBound, mNumberOfElements);
}

std::pair<PointType, PointType> Mapper::GetBoundingBoxFromIndex(const Vector3i& Indices) const {
    return Mapping::GetBoundingBoxFromIndex(Indices, mLowerBound, mUpperBound, mNumberOfElements);
}

std::pair<PointType, PointType> Mapper::GetBoundingBoxFromIndex(IndexType i, IndexType j, IndexType k) const {
    return Mapping::GetBoundingBoxFromIndex(i, j, k, mLowerBound, mUpperBound, mNumberOfElements);
}

} // End namespace tibra