// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef POINT_INCLUDE_H
#define POINT_INCLUDE_H

#include <array>

///@TODO: Add copy constructor, move etc. and const getter functions.
class Point : public std::array<double,3>
{
public:
    typedef std::array<double, 3> BaseType;

    // Default constructor
    Point()
    {}

    // Constructor
    Point(double x, double y, double z)
    {
        this->data()[0] = x;
        this->data()[1] = y;
        this->data()[2] = z;
    }

    // Assignement operator
    Point& operator=(const Point& rOther)
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

    double X() const{
        return this->operator[](0);
    }

    double Y() const{
        return this->operator[](1);
    }

    double Z() const{
        return this->operator[](2);
    }

    double& operator [] (std::size_t i){
        return this->data()[i];
    }

    double operator [] (std::size_t i) const{
        return this->data()[i];
    }
};

#endif // POINT_INCLUDE_H