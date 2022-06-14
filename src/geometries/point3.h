// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef Point3_INCLUDE_H
#define Point3_INCLUDE_H

#include <array>

class Point3 : public std::array<double,3>
{
public:
    typedef std::array<double, 3> BaseType;

    // Default constructor
    Point3()
    {}

    // Constructor
    Point3(double x, double y, double z)
    {
        this->data()[0] = x;
        this->data()[1] = y;
        this->data()[2] = z;
    }

    // Assignement operator
    Point3& operator=(const Point3& rOther)
    {
        std::array<double,3>::operator=(rOther);
        return *this;
    }

    BaseType& Coordinates(){
        return *this;
    }

    double& X(){
        return this->operator[](0);
    }

    double& Y(){
        return this->operator[](1);
    }

    double& Z(){
        return this->operator[](2);
    }

    double& operator [] (std::size_t i){
        return this->data()[i];
    }
};

#endif // Point3_INCLUDE_H