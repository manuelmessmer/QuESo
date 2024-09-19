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
 * @todo Refactor. Store elements as unique_ptr
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

    ///@}
    ///@name Operations
    ///@{

    const ElementType& GetElement(std::size_t id) const{
        return *pGetElement(id);
    }

    const ElementType* pGetElement(std::size_t id) const {
        auto found_key = mElementIdMap.find(id);
        if( found_key == mElementIdMap.end() )
            QuESo_ERROR << "ID does not exist.\n";
        return mElements[found_key->second].get();
    }

    ElementType* pGetElement(std::size_t id, bool& found){
        auto found_key = mElementIdMap.find(id);
        found = false;
        if( found_key != mElementIdMap.end() ){
            found = true;
            return mElements[found_key->second].get();
        }
        return nullptr;
    }

    const ElementContainerType& GetElements() const{
        return mElements;
    }

    void AddElement(ElementPtrType& rElement){
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
            QuESo_ERROR << "ID already exists.\n";
        }
    }

    std::size_t size() const {
        return mElements.size();
    }

    void reserve(std::size_t new_capacity){
        mElements.reserve(new_capacity);
        mElementIdMap.reserve(new_capacity);
    }

    ElementType* pGetNextElementInX(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        local_end = false;
        next_id = mGridIndexer.GetNextIndexX(id-1, local_end)+1;
        auto found_element = pGetElement(next_id, found);
        if( found == false){  // Element is not found
            local_end = true;
        }
        return found_element;
    }

    ElementType* pGetNextElementInY(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        local_end = false;
        next_id = mGridIndexer.GetNextIndexY(id-1, local_end)+1;
        auto found_element = pGetElement(next_id, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    ElementType* pGetNextElementInZ(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        local_end = false;
        next_id = mGridIndexer.GetNextIndexZ(id-1, local_end)+1;
        auto found_element = pGetElement(next_id, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    ElementType* pGetPreviousElementInX(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        local_end = false;
        next_id = mGridIndexer.GetPreviousIndexX(id-1, local_end)+1;
        auto found_element = pGetElement(next_id, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }
        return found_element;
    }

    ElementType* pGetPreviousElementInY(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        // Make sure current element exists
        // TODO:: if id >= mLastElement error
        local_end = false;
        next_id = mGridIndexer.GetPreviousIndexY(id-1, local_end)+1;
        auto found_element = pGetElement(next_id, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    ElementType* pGetPreviousElementInZ(std::size_t id, std::size_t& next_id, bool& found, bool& local_end){
        local_end = false;

        next_id = mGridIndexer.GetPreviousIndexZ(id-1, local_end)+1;
        auto found_element = pGetElement(next_id, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    bool IsLast(std::size_t id, std::size_t direction ){
        return mGridIndexer.IsLocalEnd(id-1, direction);
    }

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

    double GetVolumeOfAllIPs(){
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