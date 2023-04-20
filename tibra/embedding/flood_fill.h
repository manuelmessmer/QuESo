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
    Unique<StatusVectorType> ClassifyElements() const;

    int Fill( IndexStackType& rIndexStack, StatusVectorType& rStates ) const;

private:

    int FillDirection(IndexType Direction, IndexStackType& stack, StatusVectorType& rStates, int& rInsideCount) const;

    std::pair<PointType, PointType> GetBoundingBoxFromIndex(IndexType Index) const;

private:
    const BRepOperatorBase* mpBrepOperator;
    VectorMatrixIdUtilities mIdMapper;
    const PointType mLowerBound;
    const PointType mUpperBound;
    const Vector3i mNumberOfElements;
    PointType mDelta;

}; // End class FloodFill

///@} TIBRA Classes
} // End tibra namespace

#endif // FLOOD_FILL_INCLUDE_H