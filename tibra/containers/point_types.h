// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef VECTOR_3_INCLUDE_H
#define VECTOR_3_INCLUDE_H

#include <array>

template<typename type>
class Vector3 : public std::array<type,3>
{
public:
    typedef std::array<type,3> BaseType;

    /// Default constructor
    Vector3()
    {}

    /// Constructor
    Vector3(type x, type y, type z)
    {
        this->data()[0] = x;
        this->data()[1] = y;
        this->data()[2] = z;
    }

    /// Destructor
    ~Vector3() = default;

    /// Copy Constructor
    Vector3(const BaseType& rOther) : BaseType(rOther)
    {
    }

    /// Copy Constructor
    Vector3(const Vector3& rOther) : BaseType(rOther)
    {
    }

    /// Copy Assignement
    Vector3& operator=(const Vector3& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }


    BaseType& Coordinates(){
        return *this;
    }

    type& X(){
        return this->operator[](0);
    }

    type& Y(){
        return this->operator[](1);
    }

    type& Z(){
        return this->operator[](2);
    }

    type X() const{
        return this->operator[](0);
    }

    type Y() const{
        return this->operator[](1);
    }

    type Z() const{
        return this->operator[](2);
    }

    type& operator [] (std::size_t i){
        return this->data()[i];
    }

    type operator [] (std::size_t i) const{
        return this->data()[i];
    }
};

typedef std::size_t  SizeType;
typedef std::size_t  IndexType;

typedef Vector3<double> PointType;
typedef Vector3<double> Vector3d;
typedef Vector3<IndexType> Vector3i;


#endif // VECTOR_3_INCLUDE_H