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

#ifndef BACKGROUND_GRID_INCLUDE_HPP
#define BACKGROUND_GRID_INCLUDE_HPP

//// STL includes
#include <cstring>
#include <sstream>
#include <stdexcept>
//// Project includes
#include "queso/containers/grid_indexer.hpp"
#include "queso/containers/element.hpp"
#include "queso/containers/condition.hpp"
#include "queso/includes/settings.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  BackgroundGrid
 * @author Manuel Messmer
 * @brief  Stores elements in vector and provides fast access via Id map.
 * @note Only active elements/knot spans are stored.
*/
template<typename TElementType>
class BackgroundGrid {

public:

    ///@name Type Defintitions
    ///@{

    typedef TElementType ElementType;
    typedef typename ElementType::IntegrationPointType IntegrationPointType;
    typedef typename ElementType::BoundaryIntegrationPointType BoundaryIntegrationPointType;

    typedef Unique<ElementType> ElementPtrType;
    typedef std::vector<ElementPtrType> ElementContainerType;
    typedef std::vector<IntegrationPointType> IntegrationPointVectorType;
    typedef Unique<IntegrationPointVectorType> IntegrationPointVectorPtrType;
    typedef std::unordered_map<IndexType, IndexType> ElementIdMapType;

    typedef Condition ConditionType;
    typedef Unique<Condition> ConditionPtrType;
    typedef std::vector<std::vector<ConditionPtrType>> ConditionContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    BackgroundGrid(const SettingsBaseType& rSettings) :
        mGridIndexer(rSettings), mLastElementId(0)
    {
    }
    /// Destructor
    ~BackgroundGrid() = default;
    /// Copy Constructor
    BackgroundGrid(BackgroundGrid const& rOther) = delete;
    /// Assignement Operator
    BackgroundGrid& operator=(BackgroundGrid const& rOther) = delete;
    /// Move constructor
    BackgroundGrid(BackgroundGrid&& rOther) = delete;
    /// Move assignement operator
    BackgroundGrid& operator=(BackgroundGrid&& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns reference to element with given Id.
    /// @param ElementId
    ///@return const ElementType&
    const ElementType& GetElement(IndexType ElementId) const{
        return *pGetElement(ElementId);
    }

    /// @brief Returns raw pointer to element with given Id (const version).
    /// @param ElementId
    ///@return const ElementType*
    const ElementType* pGetElement(IndexType ElementId) const {
        auto found_key = mElementIdMap.find(ElementId);
        if( found_key != mElementIdMap.end() ){
            return mElements[found_key->second].get();
        }
        return nullptr;
    }

    /// @brief Returns raw pointer to element with given Id (non-const version).
    /// @param ElementId
    ///@return const ElementType*
    ElementType* pGetElement(IndexType ElementId){
        auto found_key = mElementIdMap.find(ElementId);
        if( found_key != mElementIdMap.end() ){
            return mElements[found_key->second].get();
        }
        return nullptr;
    }

    /// @brief Returns all active elements.
    /// @return const std::vector<Unique<Element>>&
    const ElementContainerType& GetElements() const{
        return mElements;
    }

    /// @brief Moves element into container.
    void AddElement(ElementPtrType& pElement){
        const IndexType current_id = pElement->GetId();
        const auto found_key = mElementIdMap.find(current_id);
        if( found_key == mElementIdMap.end() ){
            if( pElement->GetId() > static_cast<IndexType>(mLastElementId) ){
                mLastElementId = pElement->GetId();
            }
            mElementIdMap.insert(std::pair<IndexType, IndexType>(pElement->GetId(), mElements.size()));
            mElements.push_back(std::move(pElement));
        }
        else {
            QuESo_ERROR << "ID already exists.\n";
        }
    }

    /// @brief Returns number of stored elements.
    /// @return IndexType.
    IndexType size() const {
        return mElements.size();
    }

    /// @brief Adjust capacity of element containers.
    void reserve(IndexType NewCapacity){
        mElements.reserve(NewCapacity);
        mElementIdMap.reserve(NewCapacity);
    }

    /// @brief Returns raw ptr to next element in x-direction.
    /// @param CurrentId element id to start from.
    /// @param[out] rNextId
    /// @param[out] rIsEnd true if CurrentId is a local or global end. That means we can not move towards the given direction.
    /// @return ElementType*
    ElementType* pGetNextElementInX(IndexType CurrentId, IndexType& rNextId, bool& rIsEnd) {
        const auto [next_index, index_info] = mGridIndexer.GetNextIndexX(CurrentId-1);
        rIsEnd = ( index_info != GridIndexer::IndexInfo::middle );
        rNextId = next_index + 1;
        ElementType* p_element = pGetElement(rNextId);
        if( !p_element ){ // Element not found. This also indicates an end.
            rIsEnd = true;
        }
        return p_element;
    }

    /// @brief Returns raw ptr to next element in y-direction.
    /// @param CurrentId element id to start from.
    /// @param[out] rNextId
    /// @param[out] rIsEnd true if CurrentId is a local or global end. That means we can not move towards the given direction.
    /// @return ElementType*
    ElementType* pGetNextElementInY(IndexType CurrentId, IndexType& rNextId, bool& rIsEnd) {
        const auto [next_index, index_info] = mGridIndexer.GetNextIndexY(CurrentId-1);
        rIsEnd = ( index_info != GridIndexer::IndexInfo::middle );
        rNextId = next_index + 1;
        ElementType* p_element = pGetElement(rNextId);
        if( !p_element ){ // Element not found. This also indicates an end.
            rIsEnd = true;
        }
        return p_element;
    }

    /// @brief Returns raw ptr to next element in z-direction.
    /// @param CurrentId element id to start from.
    /// @param[out] rNextId
    /// @param[out] rIsEnd true if CurrentId is a local or global end. That means we can not move towards the given direction.
    /// @return ElementType*
    ElementType* pGetNextElementInZ(IndexType CurrentId, IndexType& rNextId, bool& rIsEnd) {
        const auto [next_index, index_info] = mGridIndexer.GetNextIndexZ(CurrentId-1);
        rIsEnd = ( index_info != GridIndexer::IndexInfo::middle );
        rNextId = next_index + 1;
        ElementType* p_element = pGetElement(rNextId);
        if( !p_element ){ // Element not found. This also indicates an end.
            rIsEnd = true;
        }
        return p_element;
    }

    /// @brief Returns raw ptr to previous element in x-direction.
    /// @param CurrentId element id to start from.
    /// @param[out] rNextId
    /// @param[out] rIsEnd true if CurrentId is a local or global end. That means we can not move towards the given direction.
    /// @return ElementType*
    ElementType* pGetPreviousElementInX(IndexType CurrentId, IndexType& rNextId, bool& rIsEnd) {
        const auto [next_index, index_info] = mGridIndexer.GetPreviousIndexX(CurrentId-1);
        rIsEnd = ( index_info != GridIndexer::IndexInfo::middle );
        rNextId = next_index + 1;
        ElementType* p_element = pGetElement(rNextId);
        if( !p_element ){ // Element not found. This also indicates an end.
            rIsEnd = true;
        }
        return p_element;
    }

    /// @brief Returns raw ptr to previous element in y-direction.
    /// @param CurrentId element id to start from.
    /// @param[out] rNextId
    /// @param[out] rIsEnd true if CurrentId is a local or global end. That means we can not move towards the given direction.
    /// @return ElementType*
    ElementType* pGetPreviousElementInY(IndexType CurrentId, IndexType& rNextId, bool& rIsEnd) {
        const auto [next_index, index_info] = mGridIndexer.GetPreviousIndexY(CurrentId-1);
        rIsEnd = ( index_info != GridIndexer::IndexInfo::middle );
        rNextId = next_index + 1;
        ElementType* p_element = pGetElement(rNextId);
        if( !p_element ){ // Element not found. This also indicates an end.
            rIsEnd = true;
        }
        return p_element;
    }

    /// @brief Returns raw ptr to previous element in z-direction.
    /// @param CurrentId element id to start from.
    /// @param[out] rNextId
    /// @param[out] rIsEnd true if CurrentId is a local or global end. That means we can not move towards the given direction.
    /// @return ElementType*
    ElementType* pGetPreviousElementInZ(IndexType CurrentId, IndexType& rNextId, bool& rIsEnd) {
        const auto [next_index, index_info] = mGridIndexer.GetPreviousIndexZ(CurrentId-1);
        rIsEnd = ( index_info != GridIndexer::IndexInfo::middle );
        rNextId = next_index + 1;
        ElementType* p_element = pGetElement(rNextId);
        if( !p_element ){ // Element not found. This also indicates an end.
            rIsEnd = true;
        }
        return p_element;
    }

    /// @brief Returns true if CurrentId is a local or global end. That means we can not move towards the given direction.
    /// @param CurrentId
    /// @param Direction – Move Direction: 0:+x, 1:-x, 2:+y, 3:-y, 4:+z, 5:-z
    /// @return bool
    bool IsEnd(IndexType CurrentId, IndexType Direction ){
        return mGridIndexer.IsEnd(CurrentId-1, Direction);
    }

    /// @brief Returns the volume stored on all integration points.
    /// @return double.
    double GetVolumeOfAllIPs() const {
        double volume = 0.0;
        const auto el_it_ptr_begin = this->begin();
        #pragma omp parallel for reduction(+ : volume)
        for( int i = 0; i < static_cast<int>(this->size()); ++i ){
            const auto& el_ptr = (*(el_it_ptr_begin + i));
            const double det_j = el_ptr->DetJ();
            const auto& r_points = el_ptr->GetIntegrationPoints();
            for( const auto& r_point : r_points ){
                volume += r_point.Weight()*det_j;
            }
        }
        return volume;
    }

    /// @todo Remove this function
    const IntegrationPointVectorPtrType pGetPoints(const char* type) const {
        IntegrationPointVectorPtrType points = MakeUnique<IntegrationPointVectorType>();
        const auto begin_el_itr_ptr = this->begin();
        for( IndexType i = 0; i < this->size(); ++i){
            const auto& el_ptr = *(begin_el_itr_ptr + i);
            IntegrationPointVectorType points_tmp;
            if( std::strcmp(type,"Trimmed") == 0 ){
                if( el_ptr->IsTrimmed() )
                    points_tmp = el_ptr->GetIntegrationPoints();
            }
            else if( std::strcmp(type,"Inside") == 0 ){
                if( !el_ptr->IsTrimmed() )
                    points_tmp = el_ptr->GetIntegrationPoints();
            }
            else if( std::strcmp(type,"All") == 0 ){
                points_tmp = el_ptr->GetIntegrationPoints();
            }
            else {
                QuESo_ERROR << "Given type '" << type << "' not available.\n";
            }
            points->insert(points->end(), points_tmp.begin(), points_tmp.end());
        }
        return points;
    }

    ///@}
    ///@name Get Iterators
    ///@{

    DereferenceIterator<typename std::vector<std::unique_ptr<ElementType>>::iterator> begin() {
        return dereference_iterator(mElements.begin());
    }

    DereferenceIterator<typename std::vector<std::unique_ptr<ElementType>>::const_iterator> begin() const {
        return dereference_iterator(mElements.begin());
    }

    DereferenceIterator<typename std::vector<std::unique_ptr<ElementType>>::iterator> end() {
        return dereference_iterator(mElements.end());
    }

    DereferenceIterator<typename std::vector<std::unique_ptr<ElementType>>::const_iterator> end() const {
        return dereference_iterator(mElements.end());
    }

    RawPointerIterator<typename std::vector<std::unique_ptr<ElementType>>::iterator> begin_to_ptr() {
        return raw_pointer_iterator(mElements.begin());
    }

    RawPointerIterator<typename std::vector<std::unique_ptr<ElementType>>::const_iterator> begin_to_ptr() const {
        return raw_pointer_iterator(mElements.begin());
    }

    RawPointerIterator<typename std::vector<std::unique_ptr<ElementType>>::iterator> end_to_ptr() {
        return raw_pointer_iterator(mElements.end());
    }

    RawPointerIterator<typename std::vector<std::unique_ptr<ElementType>>::const_iterator> end_to_ptr() const {
        return raw_pointer_iterator(mElements.end());
    }

private:

    ///@}
    ///@name Private member variables
    ///@{

    GridIndexer mGridIndexer;
    int mLastElementId;
    ElementContainerType mElements;
    ElementIdMapType mElementIdMap;

    ConditionContainerType mConditions;

    ///@}
}; // End class BackgroundGrid
///@} // End QuESo classes

} // End namespace queso
#endif // BACKGROUND_GRID_INCLUDE_HPP