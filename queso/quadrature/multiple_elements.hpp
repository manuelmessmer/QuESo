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

#ifndef MULITPLE_ELEMENTS_INCLUDE_HPP
#define MULITPLE_ELEMENTS_INCLUDE_HPP

//// STL includes
#include <queue>

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/containers/element.hpp"
#include "queso/containers/background_grid.hpp"
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace queso {

///@name QuESo Classes
///@{

/// @class  QuadratureMultipleElements.
/// @author Manuel Messmer
/// @brief  Provides assembly operations for tensor-product quadrature rules that can be used for multiple non-trimmed elements.
///         Available quadrature rules: {GGQ_Optimal, GGQ_Reduced1, GGQ_Reduced2}.
/// @details Implements algorithm from Section 3.2.1 in 10.1016/j.cma.2022.115584.
template<typename TElementType>
class QuadratureMultipleElements {

public:
    ///@name Type Defintitions
    ///@{

    using ElementType = TElementType;
    using BackgroundGridType = BackgroundGrid<ElementType>;
    using Direction = GridIndexer::Direction;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Assembles integration points for all non-trimmed elements using tensor-product quadrature rules.
    ///        These rules require less integration points by exploiting the coninuity property of B-Splines and NURBS.
    /// @details Implements algorithm from Section 3.2.1 in 10.1016/j.cma.2022.115584.
    /// @param rGrid
    /// @param rIntegrationOrder
    /// @param Method Integration method - Options: {GGQ_Optimal, GGQ_Reduced1, GGQ_Reduced2}.
    static void AssembleIPs(BackgroundGridType& rGrid, const Vector3i& rIntegrationOrder, IntegrationMethodType Method) {
        using ElementVectorType = std::vector<ElementType*>;

        // Initialize
        ComputeNeighborCoefficients(rGrid);

        // Create priority queue for all non-trimmed elements. The element with the highest neighbor coefficient is
        // always at front.
        std::priority_queue<ElementType*, ElementVectorType, CompareByCoefficient> unvisited_elements;
        for( auto& r_el : rGrid.Elements() ) {
            r_el.SetValue(ElementValues::is_visited, false);
            if( !r_el.IsTrimmed() ){
                unvisited_elements.push(&r_el);
            }
        }

        // Repeat until all elements are visited.
        while( !unvisited_elements.empty()  ){
            ElementType* p_max_el = unvisited_elements.top();
            unvisited_elements.pop();
            if( p_max_el->template GetValue<bool>(ElementValues::is_visited) ) {
                continue; // Skip already visited elements.
            }
            p_max_el->SetValue(ElementValues::is_visited, true);

            // Repetively try to expand current_box.
            ElementVectorType current_box{};
            current_box.reserve(10000);
            current_box.push_back(p_max_el);
            Vector3i current_box_sizes = {1, 1, 1};
            BoundingBoxType current_box_bounds = p_max_el->GetBoundsUVW();
            while( TryToExpandCurrentBox(rGrid, current_box, current_box_sizes, current_box_bounds) ){
            }

            // Assemble integration points.
            AssembleGGQRulesOnBox(current_box, current_box_sizes, current_box_bounds, rIntegrationOrder, Method);
        }
    }

private:
    ///@}
    ///@name Type Defintions
    ///@{

    using ElementVectorType = std::vector<ElementType*>;

    /// Helper struct to define the order in the priority queue.
    struct CompareByCoefficient {
        bool operator()(const ElementType* pLHS, const ElementType* pRHS) const {
            double l_value = pLHS->template GetValue<double>(ElementValues::neighbor_coefficient);
            double r_value = pRHS->template GetValue<double>(ElementValues::neighbor_coefficient);
            if(std::abs(l_value - r_value) < ZEROTOL) {
                return (pLHS->GetId() > pRHS->GetId());
            }
            return l_value < r_value;
        }
    };

    ///@}
    ///@name Private Operations
    ///@{

    /// @brief Computes the neighbor coefficients for each full element in the given grid.
    /// @details Implements algorithm from Fig. 8 in 10.1016/j.cma.2022.115584.
    /// @param rElements
    static void ComputeNeighborCoefficients(BackgroundGridType& rElements) {
        // Loop over all forward directions.
        constexpr std::array<Direction, 3> directions = {Direction::x_forward, Direction::y_forward, Direction::z_forward};
        for( auto dir : directions ){
            bool local_end = false;
            IndexType current_id = 1;
            IndexType next_id = 0;
            ElementVectorType neighbors;
            neighbors.reserve(20);
            // Check if element with index=1 is part of rElements.
            const auto first_element = rElements.pGetElement(1);
            if( first_element && !first_element->IsTrimmed() )
                neighbors.push_back(first_element);

            IndexType el_counter = 1;
            // Loop until all elements in rElements have beend visited/found
            while( el_counter < rElements.NumberOfActiveElements() ){
                ElementType* neighbour = rElements.pGetNextElement(current_id, dir, next_id, local_end);
                if( neighbour ){
                    el_counter++;
                    if(neighbour->IsTrimmed()){
                        local_end = true;
                    }
                    else {
                        neighbors.push_back(neighbour);
                    }
                }
                if( local_end ){
                    AssignNeighborCoefficients(neighbors);
                }
                current_id = next_id;
            }
            AssignNeighborCoefficients(neighbors);
        }
    }

    /// @brief Assigns the neighbor coefficients.
    /// @details Implements algorithm from Fig. 8 in 10.1016/j.cma.2022.115584.
    /// @param rNeighbors (Neighboring elements in a row/column).
    static void AssignNeighborCoefficients(ElementVectorType& rNeighbors) {
        const IndexType number_neighbours = rNeighbors.size();
        if( number_neighbours > 1) {
            const auto el_it_begin = rNeighbors.begin();
            for(IndexType i = 0; i < number_neighbours; ++i){
                auto el_ptr = *(el_it_begin + static_cast<std::ptrdiff_t>(i));
                const double old_value = el_ptr->template GetValue<double>(ElementValues::neighbor_coefficient);
                el_ptr->SetValue(ElementValues::neighbor_coefficient, old_value*LinearFunction(i, number_neighbours));
            }
        }
        rNeighbors.clear();
    }

    /// @brief Expands the current box while preserving the box-like shape.
    /// @details Implements algorithm from Fig. 9 in 10.1016/j.cma.2022.115584.
    /// @param rGrid
    /// @param [out] rCurrentBox
    /// @param [out] rBoxSize of current box.
    /// @param [out] rBoxBounds of current box.
    /// @return true if expansion was successful.
    static bool TryToExpandCurrentBox(BackgroundGridType& rGrid,
                                      ElementVectorType& rCurrentBox,
                                      Vector3i& rBoxSize,
                                      BoundingBoxType& rBoxBounds) {
        // Find neighbors of rCurrentBox.
        std::array<ElementVectorType,6> neighbours{}; // <- List of neighbors for each direction.
        std::array<double, 6> neighbour_coeffs = {0.0}; // <- Coefficients in each direction.
        for( const auto* p_element : rCurrentBox ){
            IndexType current_id = p_element->GetId();
            for( auto direction : GridIndexer::GetDirections() ){
                const IndexType dir_index = static_cast<IndexType>(direction);
                ElementType* p_neighbour = pNextElement(rGrid, current_id, direction );
                if( p_neighbour ){
                    const bool is_visited = p_neighbour->template GetValue<bool>(ElementValues::is_visited);
                    if( !is_visited && !p_neighbour->IsTrimmed() ){
                        double coeff = p_neighbour->template GetValue<double>(ElementValues::neighbor_coefficient);
                        neighbour_coeffs[dir_index] += coeff; // <- Sum up coefficients.
                        neighbours[dir_index].push_back(p_neighbour);
                    }
                }
            }
        }

        // Lets move towards the direction with the highest coefficient.
        // There are six possible directions, e.i., six attempts.
        for( IndexType attempts = 0; attempts < 6; ++attempts ){
            IndexType move_dir_index = static_cast<IndexType>(std::distance(neighbour_coeffs.begin(),
                std::max_element(neighbour_coeffs.begin(), neighbour_coeffs.end()))); // <- Index of move direction.
            const auto& candidates_to_add =  neighbours[move_dir_index]; // <- Neighbors with highest coefficients.

            // Get the number of neighbors on the plane orthogonal to the move direction.
            IndexType dimension_index = move_dir_index / 2;
            constexpr std::array<std::pair<IndexType, IndexType>, 3> dimension_index_map = {{
                {1, 2}, // for dimension_index 0
                {0, 2}, // for dimension_index 1
                {0, 1}  // for dimension_index 2
            }};
            const auto [i, j] = dimension_index_map[dimension_index];
            IndexType required_number_neighbors = rBoxSize[i] * rBoxSize[j];

            // Add candidates if their inclusion preserves the box shape.
            if( candidates_to_add.size() == required_number_neighbors ){
                for( auto* p_element : candidates_to_add ) {
                    p_element->SetValue(ElementValues::is_visited, true);
                    rCurrentBox.push_back(p_element);
                    // Update bounds of current box
                    const auto& r_el_bounds = p_element->GetBoundsUVW();
                    for (IndexType k = 0; k < 3; ++k) {
                        rBoxBounds.first[k]  = std::min(rBoxBounds.first[k], r_el_bounds.first[k]);
                        rBoxBounds.second[k] = std::max(rBoxBounds.second[k], r_el_bounds.second[k]);
                    }
                }
                rBoxSize[dimension_index]++;
                return true;
            }
            else { // No valid move direction.
                neighbour_coeffs[move_dir_index] = 0.0; // <- Exclude this one from the potential move directions
                // in next iteration.
            }
        }
        return false;
    }

    /// @brief Assembles the GGQ rules on the given box.
    /// @param rElements Elements within the box.
    /// @param rNumberKnotspans Dimensions of the box.
    /// @param rBounds Bounds of the box in parametric space.
    /// @param rIntegrationOrder
    /// @param Method Integration method - Options: {GGQ_Optimal, GGQ_Reduced1, GGQ_Reduced2}.
    static void AssembleGGQRulesOnBox(ElementVectorType& rElements,
                                      const Vector3i& rNumberKnotspans,
                                      const BoundingBoxType& rBounds,
                                      const Vector3i& rIntegrationOrder,
                                      IntegrationMethodType Method) {

        // Loop over all elements
        for( auto* p_el : rElements ){

            // Local lower and upper points
            const auto lower_point_param = p_el->GetBoundsUVW().first;
            const auto upper_point_param = p_el->GetBoundsUVW().second;

            std::array<std::vector<std::array<double,2>>, 3> tmp_integration_points{};

            for( IndexType direction = 0; direction < 3; ++direction){
                const double distance_global = rBounds.second[direction] - rBounds.first[direction];
                const double length_global = std::abs(rBounds.second[direction] - rBounds.first[direction]);

                const IntegrationPointFactory1D::Ip1DVectorPtrType p_ggq_points =
                    IntegrationPointFactory1D::GetGGQ(rIntegrationOrder[direction], rNumberKnotspans[direction], Method);
                const IntegrationPointFactory1D::Ip1DVectorType& r_ggq_points = *p_ggq_points;

                for( IndexType j = 0; j < r_ggq_points.size(); ++j){
                    const double position = rBounds.first[direction] + distance_global* (r_ggq_points)[j][0];
                    const double weight = length_global *  (r_ggq_points)[j][1];
                    std::array<double, 2> tmp_point = {position, weight};
                    if( lower_point_param[direction]-EPS3 <= position && position < upper_point_param[direction]-EPS3){
                        if( lower_point_param[direction]+EPS2 > position) {
                            // Make sure point is clearly inside one element.
                            tmp_point[0] = lower_point_param[direction]+2*EPS3;
                        }
                        tmp_integration_points[direction].push_back(tmp_point);
                    }
                }
            }
            const SizeType PointsInU = tmp_integration_points[0].size();
            const SizeType PointsInV = tmp_integration_points[1].size();
            const SizeType PointsInW = tmp_integration_points[2].size();

            auto& r_integration_points = p_el->GetIntegrationPoints();
            r_integration_points.reserve(PointsInU * PointsInV * PointsInW);

            for (SizeType u = 0; u < PointsInU; ++u) {
                for (SizeType v = 0; v < PointsInV; ++v) {
                    for( SizeType w = 0; w < PointsInW; ++w) {
                        const double weight = tmp_integration_points[0][u][1]*tmp_integration_points[1][v][1]*tmp_integration_points[2][w][1];
                        r_integration_points.emplace_back(tmp_integration_points[0][u][0],
                                                          tmp_integration_points[1][v][0],
                                                          tmp_integration_points[2][w][0],
                                                          weight);
                    }
                }
            }
        }
    }

    /// @brief Helper function to move within the background grid.
    /// @param rElements
    /// @param CurrentId
    /// @param Dir
    /// @return ElementType*
    static ElementType* pNextElement(BackgroundGridType& rElements, IndexType CurrentId, Direction Dir ) {
        IndexType dummy_next_id;
        bool dummy_local_end;

        return rElements.pGetNextElement(CurrentId, Dir, dummy_next_id, dummy_local_end);
    }

    /// @brief Helper function to compute the neighbor coefficient using a linear function.
    /// @param X Value
    /// @param NumNeighbors
    /// @return double.
    static double LinearFunction(IndexType X, IndexType NumNeighbors) {
        const double center = static_cast<double>(NumNeighbors-1) / 2.0;
        const double delta = std::abs(center - static_cast<double>(X));

        const double value = (1.0 - (0.9/center)*delta)*static_cast<double>(NumNeighbors);

        QuESo_ERROR_IF(value < ZEROTOL) << "Value too low\n";

        return value;
    }

    ///@}

}; // End Class QuadratureMultipleElements.

///@} End QuESo Classes.

} // End namespace queso.

#endif // MULITPLE_ELEMENTS_INCLUDE_HPP
