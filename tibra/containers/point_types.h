// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef POINT_TYPES_INCLUDE_H
#define POINT_TYPES_INCLUDE_H

//// STL includes
#include <array>
#include <cmath>
#include <iostream>

namespace tibra {

/// Forward declaration
template<typename type>
class Vector3;

/// Global Typedefs to be used in TIBRA
typedef std::size_t  SizeType;
typedef std::size_t  IndexType;

typedef Vector3<double> PointType;
typedef Vector3<double> Vector3d;
typedef Vector3<IndexType> Vector3i;

///@name TIBRA Classes
///@{

/**
 * @class  Vector3
 * @author Manuel Messmer
*/
template<typename type>
class Vector3 : public std::array<type,3>
{
public:
    ///@name Type defintion
    ///@{
    typedef std::array<type,3> BaseType;

    ///@}
    ///@name Life cycle
    ///@{

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

    /// Copy Constructor from BaseType
    Vector3(const BaseType& rOther) : BaseType(rOther)
    {
    }

    /// Copy Constructor from Vector3
    Vector3(const Vector3& rOther) : BaseType(rOther)
    {
    }

    /// Copy Assignement from Vector3
    Vector3& operator=(const Vector3& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

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

    Vector3 operator+ (const Vector3& rOther) const {
        return Vector3(this->data()[0] + rOther[0], this->data()[1] + rOther[1], this->data()[2] + rOther[2]);
    }

    Vector3 operator- (const Vector3& rOther) const {
        return Vector3(this->data()[0] - rOther[0], this->data()[1] - rOther[1], this->data()[2] - rOther[2]);
    }

    Vector3 operator*= (const type& rValue)  {
        this->data()[0] *= rValue;
        this->data()[1] *= rValue;
        this->data()[2] *= rValue;
        return *this;
    }

    Vector3 operator/= (const type& rValue)  {
        this->data()[0] /= rValue;
        this->data()[1] /= rValue;
        this->data()[2] /= rValue;
        return *this;
    }

    double Norm(){
        return std::sqrt( this->data()[0]*this->data()[0]
                         +this->data()[1]*this->data()[1]
                         +this->data()[2]*this->data()[2] );
    }

    void PrintInfo(std::ostream& rOStream) const {
        rOStream << '(' << this->data()[0] << ", " << this->data()[1] << ", " << this->data()[2] << ')';
    }
    ///@}
}; // End Vector3 class
///@} // End TIBRA classes

template<typename type>
std::ostream& operator<<(std::ostream& rOStream, const Vector3<type>& rThis)  {
    rThis.PrintInfo(rOStream);
    return rOStream;
}

} // End namespace tibra

#endif // POINT_TYPES_INCLUDE_H