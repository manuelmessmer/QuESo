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
#include "queso/containers/background_grid.hpp"
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  QuadratureMultipleElements. Provides assembly opeartions for tensor-product quadrature rules of multiple non-trimmed elements.
 * @author Manuel Messmer
 * @brief  Provides assembly for 3D quadrature rules.
 * @details Available Quadrature rules:
 *          {GGQ_Optimal, GGQ_Reduced1, GGQ_Reduced2}
 * @todo   This class requires major refactoring. Implementation should be similar to FloodFill.
*/
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

    struct CompareByCoefficient {
        bool operator()(const ElementType* a, const ElementType* b) const {
            auto av = a->template GetValue<double>(ElementValues::neighbor_coefficient);
            auto bv = b->template GetValue<double>(ElementValues::neighbor_coefficient);
            if (std::abs(av - bv) < ZEROTOL) return a->GetId() > b->GetId();
            return av < bv; // max-heap behavior (largest comes first)
        }
    };

    /// @brief Assemble tensor product quadrature rules.
    /// @param rGrid
    /// @param rNumberOfElements
    /// @param rIntegrationOrder
    /// @param Method Integration method
    static void AssembleIPs(BackgroundGridType& rGrid, const Vector3i& rNumberOfElements, const Vector3i& rIntegrationOrder, IntegrationMethodType Method) {
        // Loop over all 3 space dimensions

        ComputeNeighborCoefficients(rGrid);

        std::priority_queue<ElementType*, std::vector<ElementType*>, CompareByCoefficient> unvisited_elements;
        for( auto& r_el : rGrid.Elements() ) {
            r_el.SetValue(ElementValues::is_visited, false);
            if( !r_el.IsTrimmed() ){
                unvisited_elements.push(&r_el);
            }
        }

        std::vector<ElementType*> current_box{};
        current_box.reserve(10000);

        // Repeat until all elements are visited.
        while( !unvisited_elements.empty()  ){

            auto* p_max_el = unvisited_elements.top();
            unvisited_elements.pop();

            if( p_max_el->template GetValue<bool>(ElementValues::is_visited) ) {
                continue;
            }
            p_max_el->SetValue(ElementValues::is_visited, true);

            // Repetively increase current_box if possible.
            current_box.push_back(p_max_el);
            Vector3i current_box_sizes = {1, 1, 1};
            BoundingBoxType current_box_bounds = p_max_el->GetBoundsUVW();
            while( true ){
                // Find neighbors of current_box.
                std::array<std::vector<ElementType*>,6> neighbours{}; // <- List of neighbors for each direction.
                std::array<double, 6> neighbour_coeffs = {0.0}; // <- Coefficients in each direction.
                for( auto* p_element : current_box ){
                    IndexType current_id = p_element->GetId();
                    for( auto direction : GridIndexer::GetDirections() ){
                        const IndexType dir_index = static_cast<IndexType>(direction);
                        const auto p_neighbour = NextElement(rGrid, current_id, direction );
                        if( p_neighbour ){
                            const bool is_visited = p_neighbour->template GetValue<bool>(ElementValues::is_visited);
                            if( !is_visited && !p_neighbour->IsTrimmed() ){
                                const double coeff = p_neighbour->template GetValue<double>(ElementValues::neighbor_coefficient);
                                neighbour_coeffs[dir_index] += coeff;
                                neighbours[dir_index].push_back(p_neighbour);
                            }
                        }
                    }
                }

                // Lets move towards the direction with the highest coefficient.
                bool successful_move = false;
                IndexType tried_directions = 0;
                while( !successful_move && tried_directions++ < 6 ){

                    IndexType move_dir_index = std::distance(neighbour_coeffs.begin(),
                        std::max_element(neighbour_coeffs.begin(), neighbour_coeffs.end())); // <- Index of move direction.
                    const auto& candidates_to_add =  neighbours[move_dir_index]; // <- Neighbors with highest coefficients.

                    // Get the number of neighbors on the plane orthogonal to the move direction.
                    IndexType dimension_index = move_dir_index / 2;
                    constexpr std::array<std::pair<IndexType, IndexType>, 3> dimension_index_map = {{
                        {1, 2}, // for dimension_index 0
                        {0, 2}, // for dimension_index 1
                        {0, 1}  // for dimension_index 2
                    }};
                    const auto [i, j] = dimension_index_map[dimension_index];
                    IndexType required_number_neighbors = current_box_sizes[i] * current_box_sizes[j];

                    // Add candidates if their inclusion preserves the box shape.
                    if( candidates_to_add.size() == required_number_neighbors){
                        for (auto* p_element : candidates_to_add) {
                            p_element->SetValue(ElementValues::is_visited, true);
                            current_box.push_back(p_element);
                            // Update bounds of current box
                            const auto& r_el_bounds = p_element->GetBoundsUVW();
                            for (IndexType i = 0; i < 3; ++i) {
                                current_box_bounds.first[i]  = std::min(current_box_bounds.first[i], r_el_bounds.first[i]);
                                current_box_bounds.second[i] = std::max(current_box_bounds.second[i], r_el_bounds.second[i]);
                            }
                        }
                        successful_move = true;
                        current_box_sizes[dimension_index]++;
                    }
                    else { // No valid move direction.
                        neighbour_coeffs[move_dir_index] = 0.0;
                    }
                }

                if( tried_directions >= 6){
                    break;
                }

            }

            StoreIntegrationPoints( current_box, current_box_sizes, current_box_bounds, rIntegrationOrder, Method);
            current_box.clear();
        }
    }

private:

    static void ComputeNeighborCoefficients(BackgroundGridType& rElements) {

        constexpr std::array<Direction, 3> directions = {Direction::x_forward, Direction::y_forward, Direction::z_forward};
        for( auto dir : directions ){
            bool local_end = false;
            IndexType current_id = 1;
            IndexType next_id = 0;
            std::vector<ElementType*> neighbors;
            neighbors.reserve(20);
            // Check if element with index=1 is part of rElements.
            const auto first_element = rElements.pGetElement(1);
            if( first_element && !first_element->IsTrimmed() )
                neighbors.push_back(first_element);

            IndexType el_counter = 1;
            // Loop until all elements in rElements have beend visited/found
            while( el_counter < rElements.NumberOfActiveElements() ){
                ElementType* neighbour = rElements.GetNextElement(current_id, dir, next_id, local_end);

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

    ///@todo Refactor this
    static void AssignNeighborCoefficients(std::vector<ElementType*>& rElements) {
        const auto element_it_begin = rElements.begin();
        const IndexType number_neighbours = rElements.size();

        for(IndexType i = 0; i < number_neighbours; ++i){
            auto el_ptr = *(element_it_begin + i);
            const double old_value = el_ptr->template GetValue<double>(ElementValues::neighbor_coefficient);
            if( number_neighbours > 1)
                el_ptr->SetValue(ElementValues::neighbor_coefficient, old_value*LinearFunction(i, number_neighbours));
        }
        rElements.clear();
    }

    static ElementType* NextElement(BackgroundGridType& rElements, IndexType CurrentId, Direction Dir ) {
        IndexType dummy_next_id;
        bool dummy_local_end;

        return rElements.GetNextElement(CurrentId, Dir, dummy_next_id, dummy_local_end);
    }

    static void StoreIntegrationPoints(std::vector<ElementType*>& rElements, const Vector3i& rNumberKnotspans,
            const BoundingBoxType& rBounds, const Vector3i& rIntegrationOrder, IntegrationMethodType Method) {

        // Loop over all elements
        for( auto* p_el : rElements ){

            // Local lower and upper points
            const auto lower_point_param = p_el->GetBoundsUVW().first;
            const auto upper_point_param = p_el->GetBoundsUVW().second;

            std::array<std::vector<std::array<double,2>>, 3> tmp_integration_points{};

            for( int direction = 0; direction < 3; ++direction){
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

    static double LinearFunction(IndexType X, IndexType NumNeighbors) {
        const double center = static_cast<double>(NumNeighbors-1) / 2.0;
        const double delta = std::abs(center - X);

        const double value = (1.0 - (0.9/center)*delta)*static_cast<double>(NumNeighbors);

        QuESo_ERROR_IF(value < ZEROTOL) << "Value too low\n";

        return value;
    }
}; // End Class QuadratureMultipleElements

///@}

} // End namespace queso

#endif // MULITPLE_ELEMENTS_INCLUDE_HPP