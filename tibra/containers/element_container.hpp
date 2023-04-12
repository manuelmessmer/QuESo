// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_CONTAINER_INCLUDE_H
#define ELEMENT_CONTAINER_INCLUDE_H

//// STL includes
#include <cstring>
#include <sstream>
#include <stdexcept>
//// Project includes
#include "containers/element.hpp"
#include "utilities/parameters.h"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  ElementContainer
 * @author Manuel Messmer
 * @brief  Stores elements in vector and provides fast access via Id map.
 * @note Only active elements/knot spans are stored.
*/
class ElementContainer {

public:
    ///@name Type Defintitions
    ///@{
    typedef std::shared_ptr<Element> ElementPtrType;
    typedef std::vector<ElementPtrType> ElementVectorPtrType;
    typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
    typedef Unique<IntegrationPointVectorType> IntegrationPointVectorPtrType;
    typedef std::unordered_map<IndexType, IndexType> ElementIdMapType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ElementContainer(const Parameters& rParameters){
        mNumberOfElements = rParameters.NumberOfElements();
        mLastElementId = 0;
    }

    ///@}
    ///@name Operations
    ///@{

    const Element& GetElement(std::size_t id) const{
        return *pGetElement(id);
    }

    const ElementPtrType pGetElement(std::size_t id) const {
        auto found_key = mElementIdMap.find(id);
        if( found_key == mElementIdMap.end() )
            TIBRA_ERROR("ElementContainer::pGetElement") << "ID does not exist.\n";
        return mElements[found_key->second];
    }

    const ElementPtrType pGetElement(std::size_t id, bool& found) const{
        auto found_key = mElementIdMap.find(id);
        found = false;
        if( found_key != mElementIdMap.end() ){
            found = true;
            return mElements[found_key->second];
        }
        return nullptr;
    }

    ElementVectorPtrType::iterator begin(){
        return mElements.begin();
    }

    ElementVectorPtrType::const_iterator begin() const {
        return mElements.begin();
    }

    ElementVectorPtrType::iterator end(){
        return mElements.end();
    }

    ElementVectorPtrType::const_iterator end() const {
        return mElements.end();
    }

    const ElementVectorPtrType& GetElements() const{
        return mElements;
    }

    void AddElement(const ElementPtrType& rElement){
        const int current_id = rElement->GetId();
        auto found_key = mElementIdMap.find(current_id);
        if( found_key == mElementIdMap.end() ){
            // critical section
            if( rElement->GetId() > static_cast<IndexType>(mLastElementId) ){
                mLastElementId = rElement->GetId();
            }
            mElementIdMap.insert(std::pair<IndexType, IndexType>(rElement->GetId(), mElements.size()));
            mElements.push_back(std::move(rElement));
        }
        else {
            TIBRA_ERROR("ElementContainer::AddElement") << "ID already exists.\n";
        }
    }

    std::size_t size() const {
        return mElements.size();
    }

    void reserve(std::size_t new_capacity){
        mElements.reserve(new_capacity);
        mElementIdMap.reserve(new_capacity);
    }

    ElementPtrType pGetNextElementInX(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        local_end = false;
        next_id = id + 1;
        int next_index = id + 1;
        auto indices = GetMatrixIndicesFromVectorIndex(id);
        if( indices[0] == mNumberOfElements[0]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);
        if( found == false){  // Element is not found
            local_end = true;
        }
        return found_element;
    }

    ElementPtrType pGetNextElementInY(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        // Make sure current element exists
        // TODO:: if id >= mLastElement error
        local_end = false;
        int next_index = GetNextIndexY(id, local_end);
        next_id = next_index;

        auto indices = GetMatrixIndicesFromVectorIndex(next_id-1); // Matrix starts with 0. Here Id's start with 1.
        if( indices[1] == mNumberOfElements[1]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    ElementPtrType pGetNextElementInZ(std::size_t id, std::size_t& next_id, bool& found, bool& local_end){
        local_end = false;

        int next_index = GetNextIndexZ(id, local_end);
        next_id = next_index;
        auto indices = GetMatrixIndicesFromVectorIndex(next_id-1); // Matrix starts with 0. Here Id's start with 1.
        if( indices[2] == mNumberOfElements[2]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    ElementPtrType pGetPreviousElementInX(std::size_t id, std::size_t& next_id, bool& found, bool& local_end){
        local_end = false;
        next_id = id - 1;
        int next_index = id - 1;
        auto indices = GetMatrixIndicesFromVectorIndex(id);
        if( indices[0] == mNumberOfElements[0]-1) {
            local_end = true;
        }

        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }
        return found_element;
    }

    ElementPtrType pGetPreviousElementInY(std::size_t id, std::size_t& next_id, bool& found, bool& local_end){
        // Make sure current element exists
        // TODO:: if id >= mLastElement error
        local_end = false;
        int next_index = GetPreviousIndexY(id, local_end);
        next_id = next_index;

        auto indices = GetMatrixIndicesFromVectorIndex(next_id-1); // Matrix starts with 0. Here Id's start with 1.
        if( indices[1] == mNumberOfElements[1]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    ElementPtrType pGetPreviousElementInZ(std::size_t id, std::size_t& next_id, bool& found, bool& local_end){
        local_end = false;

        int next_index = GetPreviousIndexZ(id, local_end);
        next_id = next_index;

        auto indices = GetMatrixIndicesFromVectorIndex(next_id-1); // Matrix starts with 0. Here Id's start with 1.
        if( indices[2] == mNumberOfElements[2]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    bool IsLast(std::size_t id, std::size_t direction ){
        auto indices = GetMatrixIndicesFromVectorIndex(id-1);

        switch( direction )
        {
        case 0: // Forward x
            return (indices[0] == (mNumberOfElements[0]-1));
        case 1: // Backward X
            return (indices[0] == 0);
        case 2: // Forward Y
            return (indices[1] == (mNumberOfElements[1]-1));
        case 3: // Backward Y
            return (indices[1] == 0);
        case 4: // Forward Z
            return (indices[2] == (mNumberOfElements[2]-1));
        case 5: // Backward Z
            return (indices[2] == 0);
        default:
            TIBRA_ERROR("ElementContainer::IsLast") << "There are only 6 different directions! \n";
        }
    }

    const IntegrationPointVectorPtrType pGetPoints(const char* type) const {
        IntegrationPointVectorPtrType points = MakeUnique<IntegrationPointVectorType>();
        const auto begin_el_itr_ptr = this->begin();
        for( IndexType i = 0; i < this->size(); ++i){
            auto el_itr = *(begin_el_itr_ptr + i);
            IntegrationPointVectorType points_tmp;
            if( std::strcmp(type,"Trimmed") == 0 ){
                if( el_itr->IsTrimmed() )
                    points_tmp = el_itr->GetIntegrationPoints();
            }
            else if( std::strcmp(type,"Inside") == 0 ){
                if( !el_itr->IsTrimmed() )
                    points_tmp = el_itr->GetIntegrationPoints();
            }
            else if( std::strcmp(type,"All") == 0 ){
                points_tmp = el_itr->GetIntegrationPoints();
            }
            else {
                TIBRA_ERROR("ElementContainer::pGetPoints") << "Given type '" << type << "' not available.\n";
            }
            points->insert(points->end(), points_tmp.begin(), points_tmp.end());
        }
        return points;
    }

    double GetVolumeOfAllIPs(){
        const auto p_points = pGetPoints("All");
        double weight = 0.0;
        const IndexType num_points = p_points->size();
        const auto it_begin = p_points->begin();
        #pragma omp parallel for reduction(+ : weight)
        for( int i = 0; i < static_cast<int>(num_points); ++i ){
            auto it = it_begin + i;
            weight += it->GetWeight();
        }

        const auto parameters = (*this->begin())->GetParameters();
        const auto jacobian = parameters.UpperBound() - parameters.LowerBound();
        return weight*(jacobian[0]*jacobian[1]*jacobian[2]);
    }

private:

    IndexType GetNextIndexX(IndexType i){
        return i + 1;
    }

    IndexType GetNextIndexY(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[1] < mNumberOfElements[1]-1) {
            indices[1] += 1;
        }
        else if( indices[0] < mNumberOfElements[0]-1){
            indices[0] += 1;
            indices[1] = 0;
        }
        else {
            indices[2] += 1;
            indices[1] = 0;
            indices[0] = 0;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    IndexType GetNextIndexZ(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[2] < mNumberOfElements[2]-1) {
            indices[2] += 1;
        }
        else if( indices[0] < mNumberOfElements[0]-1){
            indices[0] += 1;
            indices[2] = 0;
        }
        else {
            indices[1] += 1;
            indices[2] = 0;
            indices[0] = 0;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    IndexType GetPreviousIndexX(IndexType i){
        return i - 1;
    }

    IndexType GetPreviousIndexY(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[1] > 0) {
            indices[1] -= 1;
        }
        else if( indices[0] > 0){
            indices[0] -= 1;
            indices[1] = mNumberOfElements[1]-1;
        }
        else {
            indices[2] -= 1;
            indices[1] = mNumberOfElements[1]-1;
            indices[0] = mNumberOfElements[0]-1;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    IndexType GetPreviousIndexZ(IndexType i, bool& local_end){
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[2] > 0 ){
            indices[2] -= 1;
        }
        else if( indices[0] > 0){
            indices[0] -= 1;
            indices[2] = mNumberOfElements[2]-1;;
        }
        else {
            indices[1] -= 1;
            indices[2] = mNumberOfElements[2]-1;;
            indices[0] = mNumberOfElements[0]-1;;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    inline std::array<IndexType,3> GetMatrixIndicesFromVectorIndex(
        const IndexType Index) noexcept
    {
        std::array<IndexType,3> result;
        const IndexType index_in_row_column_plane = Index % (mNumberOfElements[0]*mNumberOfElements[1]);
        result[0] = index_in_row_column_plane % mNumberOfElements[0]; // row
        result[1] = index_in_row_column_plane / mNumberOfElements[0]; // column
        result[2] = Index / (mNumberOfElements[0]*mNumberOfElements[1]);   // depth

        return result;
    }

    inline IndexType GetVectorIndexFromMatrixIndices(
        const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex) noexcept
    {
        return DepthIndex * (mNumberOfElements[1]*mNumberOfElements[0]) + ColumnIndex * mNumberOfElements[0] + RowIndex;
    }

    ///@}
    ///@name Private member variables
    ///@{

    int mLastElementId;
    ElementVectorPtrType mElements{};
    ElementIdMapType mElementIdMap{};
    Vector3i mNumberOfElements{};
    ///@}
}; // End class Element container
///@} // End TIBRA classes

} // End namespace tibra
#endif // ELEMENT_CONTAINER_INCLUDE_H