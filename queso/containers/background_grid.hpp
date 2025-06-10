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

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/containers/grid_indexer.hpp"
#include "queso/containers/element.hpp"
#include "queso/containers/condition.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  BackgroundGrid
 * @author Manuel Messmer
 * @brief  Stores active elements and conditions.
 *         Indexing is based on GridIndexer.
 * @tparam TElementType must provide the following typedefs: IntegrationPointType, BoundaryIntegrationPointType.
 * @see    grid_indexer.hpp
**/
template<typename TElementType>
class BackgroundGrid {

public:

    ///@name Type Defintitions
    ///@{

    using ElementType = TElementType;
    using IntegrationPointType = typename ElementType::IntegrationPointType;
    using BoundaryIntegrationPointType = typename ElementType::BoundaryIntegrationPointType;

    using ElementPtrType = Unique<ElementType>;
    using ElementContainerType = std::vector<ElementPtrType>;
    using ElementIdMapType = std::unordered_map<IndexType, IndexType>;

    using ConditionType = Condition<ElementType>;
    using ConditionPtrType = Unique<ConditionType>;
    using ConditionContainerType = std::vector<ConditionPtrType>;

    using MainDictionaryType = Dictionary<key::MainValuesTypeTag>;
    using Direction = GridIndexer::Direction;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rSettings
    BackgroundGrid(const MainDictionaryType& rSettings) :
        mGridIndexer(rSettings)
    {
    }

    /// Destructor
    ~BackgroundGrid() = default;
    /// Copy Constructor
    BackgroundGrid(BackgroundGrid const& rOther) = delete;
    /// Assignement Operator
    BackgroundGrid& operator=(BackgroundGrid const& rOther) = delete;
    /// Move constructor
    BackgroundGrid(BackgroundGrid&& rOther) noexcept = delete;
    /// Move assignement operator
    BackgroundGrid& operator=(BackgroundGrid&& rOther) noexcept = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns reference to element with given Id.
    /// @param ElementId
    /// @return const ElementType&
    const ElementType& GetElement(IndexType ElementId) const{
        return *pGetElement(ElementId);
    }

    /// @brief Returns raw pointer to element with given Id (const version).
    /// @param ElementId
    /// @return const ElementType*
    const ElementType* pGetElement(IndexType ElementId) const {
        auto found_key = mElementIdMap.find(ElementId);
        if( found_key != mElementIdMap.end() ){
            return mElements[found_key->second].get();
        }
        return nullptr;
    }

    /// @brief Returns raw pointer to element with given Id (non-const version).
    /// @param ElementId
    /// @return const ElementType*
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

    /// @brief Returns all conditions.
    /// @return const std::vector<Unique<Element>>&
    const ConditionContainerType& GetConditions() const{
        return mConditions;
    }

    /// @brief Adds a new condition to the background grid.
    /// @param ConditionPtrType
    void AddCondition(ConditionPtrType& pCondition){
        mConditions.push_back(std::move(pCondition));
    }

    /// @todo Add the following functions for invalid elements.
    // /// @brief Returns raw pointer to a potentially stored invalid element with given Id (non-const version).
    // /// @param ElementId
    // ///@return const ElementType*
    // const ElementType* pGetInvalidElement(IndexType ElementId) const{
    //     auto found_key = mInvalidElementIdMap.find(ElementId);
    //     if( found_key != mInvalidElementIdMap.end() ){
    //         return mInvalidElements[found_key->second].get();
    //     }
    //     return nullptr;
    // }

    // /// @brief Moves element into container.
    // void AddInvalidElement(ElementPtrType& pElement){
    //     const IndexType current_id = pElement->GetId();
    //     const auto found_key = mInvalidElementIdMap.find(current_id);
    //     if( found_key == mInvalidElementIdMap.end() ){
    //         mInvalidElementIdMap.insert(std::pair<IndexType, IndexType>(pElement->GetId(), mInvalidElements.size()));
    //         mInvalidElements.push_back(std::move(pElement));
    //     }
    //     else {
    //         QuESo_ERROR << "ID already exists.\n";
    //     }
    // }

    /// @brief Adds element to container. Element is moved into container.
    /// @param pElement
    void AddElement(ElementPtrType& pElement){
        const IndexType current_id = pElement->GetId();
        const auto found_key = mElementIdMap.find(current_id);
        if( found_key == mElementIdMap.end() ){
            mElementIdMap.insert(std::pair<IndexType, IndexType>(pElement->GetId(), mElements.size()));
            mElements.push_back(std::move(pElement));
        }
        else {
            QuESo_ERROR << "ID already exists.\n"; // Make this assert instead?
        }
    }

    /// @brief Returns number of stored conditions.
    /// @return IndexType.
    IndexType NumberOfConditions() const {
        return mConditions.size();
    }

    /// @brief Returns number of stored elements.
    /// @return IndexType.
    IndexType NumberOfActiveElements() const {
        return mElements.size();
    }

    /// @brief Adjusts capacity of element container.
    /// @param NewCapacity.
    void ReserveElements(IndexType NewCapacity){
        mElements.reserve(NewCapacity);
        mElementIdMap.reserve(NewCapacity);
    }

    /// @brief Adjusts capacity of condition container.
    /// @param NewCapacity.
    void ReserveConditions(IndexType NewCapacity){
        mConditions.reserve(NewCapacity);
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

    /// @brief Returns raw ptr to next element.
    /// @param CurrentId element id to start from.
    /// @param Dir Move direction.
    /// @param[out] rNextId
    /// @param[out] rIsEnd true if CurrentId is a local or global end. That means we can not move towards the given direction.
    /// @return ElementType*
    ElementType* pGetNextElement(IndexType CurrentId, Direction Dir, IndexType& rNextId, bool& rIsEnd ) {
        switch( Dir )
        {
        case Direction::x_forward:
            return pGetNextElementInX(CurrentId, rNextId, rIsEnd);
        case Direction::x_backward:
            return pGetPreviousElementInX(CurrentId, rNextId, rIsEnd);
        case Direction::y_forward:
            return pGetNextElementInY(CurrentId, rNextId, rIsEnd);
        case Direction::y_backward:
            return pGetPreviousElementInY(CurrentId, rNextId, rIsEnd);
        case Direction::z_forward:
            return pGetNextElementInZ(CurrentId, rNextId, rIsEnd);
        case Direction::z_backward:
            return pGetPreviousElementInZ(CurrentId, rNextId, rIsEnd);
        default:
            QuESo_ASSERT(false, "There are only 6 different directions.\n");
            return nullptr;
        }
    }

    /// @brief Returns true if CurrentId is a local or global end. That means we can not move towards the given direction.
    /// @param CurrentId
    /// @param Direction – Move Direction.
    /// @return bool
    bool IsEnd(IndexType CurrentId, Direction Dir ){
        return mGridIndexer.IsEnd(CurrentId-1, Dir);
    }

    /// @brief Returns the volume stored on all integration points.
    /// @return double.
    /// @todo Remove this function
    double GetVolumeOfAllIPs() const {
        double volume = 0.0;
        const auto el_it_ptr_begin = ElementsBegin();
        #pragma omp parallel for reduction(+ : volume)
        for( int i = 0; i < static_cast<int>(NumberOfActiveElements()); ++i ){
            const auto& el_ptr = *(el_it_ptr_begin + i);
            const double det_j = el_ptr->DetJ();
            const auto& r_points = el_ptr->GetIntegrationPoints();
            for( const auto& r_point : r_points ){
                volume += r_point.Weight()*det_j;
            }
        }
        return volume;
    }
    /// @brief Returns total number of integration points.
    /// @return IndexType
    /// @todo Remove this function
    IndexType NumberOfIntegrationPoints() const {
        IndexType number = 0;
        const auto el_it_ptr_begin = ElementsBegin();
        #pragma omp parallel for reduction(+ : number)
        for( int i = 0; i < static_cast<int>(NumberOfActiveElements()); ++i ){
            const auto& el_ptr = *(el_it_ptr_begin + i);
            const auto& r_points = el_ptr->GetIntegrationPoints();
            number += r_points.size();
        }
        return number;
    }

    ///@}
    ///@name Get Iterators
    ///@{

    //////////////////
    //// Elements ////
    //////////////////

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    auto ElementsBegin() {
        return dereference_iterator<typename ElementContainerType::iterator>(mElements.begin());
    }

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    auto ElementsBegin() const {
        return dereference_iterator<typename ElementContainerType::const_iterator>(mElements.begin());
    }

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    auto ElementsEnd() {
        return dereference_iterator<typename ElementContainerType::iterator>(mElements.end());
    }

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    auto ElementsEnd() const {
        return dereference_iterator<typename ElementContainerType::const_iterator>(mElements.end());
    }

    /// @brief Returns range of dereferenced iterators for the elements container (non-const version).
    /// @return Range.
    auto Elements() {
        using IteratorType = DereferenceIterator<typename ElementContainerType::iterator>;
        return Range<IteratorType>{ ElementsBegin(), ElementsEnd() };
    }

    /// @brief Returns range of dereferenced iterators for the elements container (const version).
    /// @return Range.
    auto Elements() const {
        using IteratorType = DereferenceIterator<typename ElementContainerType::const_iterator>;
        return Range<IteratorType>{ ElementsBegin(), ElementsEnd() };
    }

    ////////////////////
    //// Conditions ////
    ////////////////////

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    DereferenceIterator<typename ConditionContainerType::iterator> ConditionsBegin() {
        return dereference_iterator(mConditions.begin());
    }

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    DereferenceIterator<typename ConditionContainerType::const_iterator> ConditionsBegin() const {
        return dereference_iterator(mConditions.begin());
    }

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    DereferenceIterator<typename ConditionContainerType::iterator> ConditionsEnd() {
        return dereference_iterator(mConditions.end());
    }

    /// @brief Returns dereferenced iterator. This means iterator does not point to UniquePtr,
    ///        but directly to the actual object.
    /// @return DereferenceIterator
    DereferenceIterator<typename ConditionContainerType::const_iterator> ConditionsEnd() const {
        return dereference_iterator(mConditions.end());
    }

    /// @brief Returns range of dereferenced iterators for the condition container (non-const version).
    /// @return Range.
    auto Conditions() {
        using IteratorType = DereferenceIterator<typename ConditionContainerType::iterator>;
        return Range<IteratorType>{ ConditionsBegin(), ConditionsEnd() };
    }

    /// @brief Returns range of dereferenced iterators for the condition container (const version).
    /// @return Range.
    auto Conditions() const {
        using IteratorType = DereferenceIterator<typename ConditionContainerType::const_iterator>;
        return Range<IteratorType>{ ConditionsBegin(), ConditionsEnd() };
    }

private:

    ///@}
    ///@name Private member variables
    ///@{

    GridIndexer mGridIndexer;

    ElementContainerType mElements;
    ElementIdMapType mElementIdMap;

    // ElementContainerType mInvalidElements;
    // ElementIdMapType mInvalidElementIdMap;

    ConditionContainerType mConditions;

    ///@}
}; // End class BackgroundGrid
///@} // End QuESo classes

} // End namespace queso
#endif // BACKGROUND_GRID_INCLUDE_HPP