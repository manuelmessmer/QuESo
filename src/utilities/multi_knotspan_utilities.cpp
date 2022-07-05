#include "utilities/mapping_utilities.h"
#include "utilities/multi_knotspan_utilities.h"
#include "utilities/integration_points/integration_points_factory.h"

typedef Element::IntegrationPoint1DVectorType IntegrationPoint1DVectorType;

void MultiKnotspanUtilities::ComputeIntegrationPoints(ElementContainer& rElements, const Parameters& rParameters){

    // Loop over all 3 space dimensions
    // i = 0: x
    // i = 1: y
    // i = 2: z
    for( IndexType i = 0; i <= 2; ++i){
        bool local_end = false;
        bool found = false;
        IndexType current_id = 1;
        IndexType next_id = 0;
        ElementContainer::ElementVectorPtrType door_to_door_neighbours;
        door_to_door_neighbours.reserve(20);
        // Check if element with index=1 is part of rElements.
        auto first_element = rElements.pGetElement(1, found);
        if( found && !first_element->IsTrimmed() )
            door_to_door_neighbours.push_back(first_element);

        if( rParameters.NumberOfElements()[i] == 1){
            if(door_to_door_neighbours.size() > 0 ){
                    Store1DIntegrationPoints(door_to_door_neighbours, i, rParameters);
            }
            door_to_door_neighbours.clear();
        }
        std::size_t active_element_counter = 1;
        // Loop until all elements in rElements have beend visited/found
        while( active_element_counter < rElements.size() ){
            ElementContainer::ElementPtrType neighbour;

            if(i == 0)
                neighbour = rElements.pGetNextElementInX(current_id, next_id, found, local_end);
            else if(i==1)
                neighbour = rElements.pGetNextElementInY(current_id, next_id, found, local_end);
            else
                neighbour = rElements.pGetNextElementInZ(current_id, next_id, found, local_end);

            if( found ){
                active_element_counter++;
                if(neighbour->IsTrimmed()){ // Is trimmed can only be evaluated if element is found (otherwise nullptr)
                    local_end = true;
                }
                else {
                    door_to_door_neighbours.push_back(neighbour);
                }
            }

            if( door_to_door_neighbours.size() >= 40 ){ // Maximum number of points is reached or element is trimmed
                 local_end = true;
            }

            if( local_end ){
                if(door_to_door_neighbours.size() > 0 ){
                    Store1DIntegrationPoints(door_to_door_neighbours, i, rParameters);
                }
                door_to_door_neighbours.clear();
            }
            current_id = next_id;
        }
        if(door_to_door_neighbours.size() > 0 ){
            Store1DIntegrationPoints(door_to_door_neighbours, i, rParameters);
        }
        door_to_door_neighbours.clear();
    }

    // Assemble integration points on element level.
    Assemble1DIntegrationPoints(rElements, rParameters);

}

void MultiKnotspanUtilities::Store1DIntegrationPoints(ElementContainer::ElementVectorPtrType& rElements, IndexType direction, const Parameters& rParameters){

    const int number_of_neighbours = rElements.size();
    const auto global_lower_point_param = rElements[0]->GetLocalLowerPoint();
    const auto global_upper_point_param = rElements[rElements.size()-1]->GetLocalUpperPoint();

    const double distance_global = global_upper_point_param[direction] - global_lower_point_param[direction];
    const double length_global = std::abs(global_upper_point_param[direction] - global_lower_point_param[direction]);
    std::cout << "length global: " << length_global << std::endl;
    const int polynomial_degree = rParameters.Order()[direction];
    // Loop over all elements
    const auto element_it_begin = rElements.begin();

    // const auto intial_lower_point_param = (*element_it_begin)->GetLocalLowerPoint();
    // const auto intial_upper_point_param = (*element_it_begin)->GetLocalUpperPoint();
    // const auto length_knot_span = intial_upper_point_param[direction] - intial_lower_point_param[direction];

    // int extenstion = 0;
    // if (number_of_neighbours%2 == 0 ){
    //     std::cout << "True " << std::endl;
    //     extenstion = 21;
    //     distance_global += extenstion*length_knot_span;
    //     global_lower_point_param[direction] -= 11*length_knot_span;
    //     global_upper_point_param[direction] += 10*length_knot_span;
    // }
    // else{
    //     extenstion = 20;
    //     distance_global += extenstion*length_knot_span;
    //     global_lower_point_param[direction] -= 10*length_knot_span;
    //     global_upper_point_param[direction] += 10*length_knot_span;
    // }

    // length_global *= (extenstion+number_of_neighbours)/ ( (double) number_of_neighbours );


    for( int i = 0; i < rElements.size(); ++i){
        auto element_it = *(element_it_begin + i);

        // Local lower and upper points
        const auto lower_point_param = element_it->GetLocalLowerPoint();
        const auto upper_point_param = element_it->GetLocalUpperPoint();

        const double distance = upper_point_param[direction] - lower_point_param[direction];
        const double length = std::abs(upper_point_param[direction] - lower_point_param[direction]);

        const auto method = IntegrationPointFactory::ReducedOrder2;
        const std::vector<std::array<double, 2>>& integration_point_list =
            IntegrationPointFactory::GetIntegrationPoints(polynomial_degree, number_of_neighbours, method);

        IntegrationPoint1DVectorType& integration_points_1d = element_it->IntegrationPoints1D(direction);
        for( int j = 0; j < integration_point_list.size(); ++j){
            double position = global_lower_point_param[direction] + distance_global* integration_point_list[j][0];
            double weight = length_global *  integration_point_list[j][1];
            std::array<double, 2> tmp_point = {position, weight};
            // Keep this in mind.
            if( position > lower_point_param[direction] && position <= upper_point_param[direction]){
                integration_points_1d.push_back(tmp_point);
            }
        }
    }
}

void MultiKnotspanUtilities::Assemble1DIntegrationPoints(ElementContainer& rElements, const Parameters& rParameter){

    double total_weight_inside = 0.0;
    int number_elements_inside = 0;

    auto elemet_it_begin = rElements.begin();
    for( int i = 0; i < rElements.size(); ++i){
        auto element_it = elemet_it_begin + i;
        if( !(*element_it)->IsTrimmed() ){
            number_elements_inside++;
            auto integration_points_x = (*element_it)->IntegrationPoints1D(0);
            auto integration_points_y = (*element_it)->IntegrationPoints1D(1);
            auto integration_points_z = (*element_it)->IntegrationPoints1D(2);

            const int PointsInU = integration_points_x.size();
            const int PointsInV = integration_points_y.size();
            const int PointsInW = integration_points_z.size();

            for (SizeType u = 0; u < PointsInU; ++u) {
                for (SizeType v = 0; v < PointsInV; ++v) {
                    for( SizeType w = 0; w < PointsInW; ++w) {
                        const double weight = integration_points_x[u][1]*integration_points_y[v][1]*integration_points_z[w][1];
                        total_weight_inside += weight;
                        (*element_it)->GetIntegrationPointsInside().push_back(
                                                        IntegrationPoint( integration_points_x[u][0],
                                                                          integration_points_y[v][0],
                                                                          integration_points_z[w][0],
                                                                          weight ) );
                    }
                }
            }
        }
    }

    // Apply correction
    const auto number_elements = rParameter.NumberOfElements();
    const double targeted_weight = ((double)number_elements_inside) / ((double)(number_elements[0]*number_elements[1]*number_elements[2]));
    double correction_factor = (double)targeted_weight/total_weight_inside;
    //correction_factor = correction_factor+0.5*(1-correction_factor);
    std::cout << "Number el inside: " << number_elements_inside << std::endl;
    std::cout << "Correction factor: " << correction_factor << std::endl;
    for( int i = 0; i < rElements.size(); ++i){
        auto element_it = elemet_it_begin + i;
        if( !(*element_it)->IsTrimmed() ){
            auto& integration_points = (*element_it)->GetIntegrationPointsInside();
            for(auto& point : integration_points){
                point.SetWeight(point.GetWeight()*correction_factor);
            }
        }
    }
}