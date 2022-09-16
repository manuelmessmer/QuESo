// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#include "utilities/ggq_utilities.h"

GGQRule::IntegrationPoint1DVectorPtrType GGQRule::GetGGQ_Rule(SizeType p, SizeType r, SizeType e, double a, double b ){
//void GGQRule::GetGGQ_Rule(SizeType p, SizeType r, SizeType e, double a, double b ){

    const SizeType n = (p+1)*2 + (e-1)*(p-r) - p - 1;
    const SizeType m = std::ceil(n/2.0);

    std::cout << "n " << n << std::endl;
    std::cout << "m " << m << std::endl;

    std::vector<std::array<double,2>> points(m);

    // std::vector<double> weights[m];

    // std::vector<double> bnodes{};
    // std::vector<double> bweigths{};
    // std::vector<double> inodes{};
    // std::vector<double> iweigths{};
    // std::vector<double> mnodes{};
    // std::vector<double> mweights{};

    std::vector<std::vector<std::array<double, 2>>> base_points{};

    // Get correct base rule points
    if( p == 4 && r == 0 ){
        if( e % 2 == 0 ){
            // 'quadrule_p4_r0/quadrule_even_r0_p4.mat'
            base_points = S_4_0_even;
        }
        else {
            base_points = S_4_0_odd;
            // 'quadrule_p4_r0/quadrule_odd_r0_p4.mat'
        }
    }

    SizeType m1 = base_points[0].size(); // boundary nodes
    SizeType m2 = base_points[1].size(); // internal nodes
    SizeType m3 = base_points[2].size(); // center nodes

    if( 2*m == n ){
        // TODO..
        std::cout << "2*m = n" << std::endl;
    }
    else { // For odd number of nodes
        if( m > 2*(m1+m3)-1 ){
            std::cout << "in here?? " << std::endl;
            const SizeType z = std::ceil(0.5*m);
            std::copy_n(base_points[0].begin(), m1, points.begin());
            SizeType left = std::ceil(base_points[0].back()[0]);
            const SizeType right = std::ceil(base_points[1].back()[0]);
            SizeType ii = m1;
            while( ii < z ){
                std::copy_n(base_points[1].begin(), m2, points.begin()+ii);
                std::for_each(points.begin()+ii, points.begin()+ii+m2, [left](auto& rValue) { rValue[0] +=left;});
                left += right;
                ii += m2;
            }

            std::copy_n(base_points[2].begin(), m3, points.begin()+z-m3);
            std::for_each(points.begin()+z-m3, points.begin()+z, [e](auto& rValue) { rValue[0] +=0.5*e;});

            std::reverse_copy(points.begin(), points.begin()+z, points.end()- z );
            std::for_each(points.end()-z, points.end(), [e](auto& rValue) { rValue[0] = e - rValue[0];});
            std::cout << "Sucess?? " << std::endl;
        }
        else {
            // Here compute precomputed rule!!
        }
    }

    const double h = (b-a) / e;
    std::for_each(points.begin(), points.end(), [a, h](auto& rValue) { rValue[0] = a + h*rValue[0];
                                                                        rValue[1] *= h; });
    int count = 1;
    for( auto point : points){
        count++;
        std::cout << count << ": " << point[0] << ", " << point[1] << std::endl;
    }

    return std::make_unique<IntegrationPoint1DVectorType>(points);
}