// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_CONTAINER_INCLUDE_H
#define ELEMENT_CONTAINER_INCLUDE_H

// External includes
#include <cstring>
#include <sstream>
#include <stdexcept>

// Project includes
#include "geometries/element.h"
#include "utilities/parameters.h"


class ElementContainer {

public:

    // Typedefs
    typedef std::shared_ptr<Element> ElementPtrType;
    typedef std::vector<ElementPtrType> ElementVectorPtrType;
    typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
    typedef std::unique_ptr<IntegrationPointVectorType> IntegrationPointVectorPtrType;
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef std::unordered_map<IndexType, IndexType> ElementHashMap;

    // Constructor
    ElementContainer(const Parameters& rParameters){
        mNumberOfElements = rParameters.NumberOfElements();
        mLastElementId = 0;
    }

    Element& GetElement(std::size_t id){
        return *pGetElement(id);
    }

    ElementPtrType pGetElement(std::size_t id) {
        auto found_key = mElementHashMap.find(id);
        if( found_key == mElementHashMap.end() )
            throw std::runtime_error("ID does not exist");

        return mElements[found_key->second];
    }

    ElementPtrType pGetElement(std::size_t id, bool& found){
        auto found_key = mElementHashMap.find(id);
        found = false;
        if( found_key != mElementHashMap.end() ){
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

    ElementVectorPtrType& GetElements(){
        return mElements;
    }

    void AddElement(ElementPtrType& rElement){
        const int current_id = rElement->GetId();
        auto found_key = mElementHashMap.find(current_id);
        if( found_key == mElementHashMap.end() ){
            // critical section
            if( rElement->GetId() > mLastElementId){
                mLastElementId = rElement->GetId();
            }
            mElementHashMap.insert(std::pair<IndexType, IndexType>(rElement->GetId(), mElements.size()));
            mElements.push_back(std::move(rElement));
        }
        else {
            throw std::runtime_error("ID already exists");
        }
    }

    std::size_t size() const {
        return mElements.size();
    }

    void reserve(std::size_t new_capacity){
        mElements.reserve(new_capacity);
        mElementHashMap.reserve(new_capacity);
    }

    ElementPtrType pGetNextElementInX(std::size_t id, std::size_t& next_id, bool& found, bool& local_end){
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

    ElementPtrType pGetNextElementInY(std::size_t id, std::size_t& next_id, bool& found, bool& local_end){
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
            throw std::invalid_argument("Element Container: There are only 6 different directions!" );
        }
    }

    const IntegrationPointVectorPtrType pGetPoints(const char* type) const {
        IntegrationPointVectorPtrType points = std::make_unique<IntegrationPointVectorType>();
        const auto begin_el_itr_ptr = this->begin();
        for( int i = 0; i < this->size(); ++i){
            auto el_itr = *(begin_el_itr_ptr + i);
            IntegrationPointVectorType points_tmp;
            if( std::strcmp(type,"Trimmed") == 0 || std::strcmp(type,"All") == 0){
                points_tmp = el_itr->GetIntegrationPointsTrimmed();
            }
            else if( std::strcmp(type,"Inside") == 0 || std::strcmp(type,"All") == 0){
                points_tmp = el_itr->GetIntegrationPointsInside();
            }
            else {
                std::stringstream error_message;
                error_message << "Element Container: Given type '" << type << "' not available.";
                throw std::invalid_argument(error_message.str());
            }
            points->insert(points->end(), points_tmp.begin(), points_tmp.end());
        }
        return std::move(points);
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

    int mLastElementId;
    ElementVectorPtrType mElements{};
    ElementHashMap mElementHashMap{};
    std::array<int,3> mNumberOfElements{};
};

#endif