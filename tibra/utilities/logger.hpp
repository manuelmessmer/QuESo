// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef LOGGER_INCLUDE_HPP
#define LOGGER_INCLUDE_HPP

//// STL includes
#include <stdexcept>
#include <iostream>
#include <sstream>

namespace tibra {

///@name TIBRA Classes
///@{

///
/**
 * @class  Logger
 * @author Manuel Messmer
*/
class Logger {
public:
    ///@}
    ///@name Life cycle
    ///@{

    /// Default Constructor
    Logger() {}

    /// Constructor
    Logger(std::string rWhere ) {
        AppendMessage( rWhere );
        AppendMessage( " :: " );

    }

    /// Destructor
    ~Logger() {
        std::cout << mMessage << '\n';
    }

    ///@}
    ///@name Operations
    ///@{
    void AppendMessage(std::string const& rMessage) {
		mMessage.append(rMessage);
	}

    Logger& operator << (const char * rString) {
        AppendMessage(rString);
        return *this;
    }
    ///@}

private:
    ///@}
    ///@name Private Members
    ///@{
    std::string mMessage;

    ///@}
};

#define TIBRA_INFO Logger()
#define TIBRA_INFO_IF(Conditional) if(Conditional) Logger()

///
/**
 * @class  Exception
 * @author Manuel Messmer
*/
class Exception : public std::exception {
public:
    ///@}
    ///@name Life cycle
    ///@{

    /// Constructor
    Exception(std::string rWhere ) {
        AppendMessage( "Where: " );
        AppendMessage( rWhere );
        AppendMessage( " :: Message: " );
    }

    ///@}
    ///@name Operations
    ///@{

    void AppendMessage(std::string const& rMessage) {
		mMessage.append(rMessage);
	}

    Exception& operator << (const char * rString) {
        AppendMessage(rString);
        return *this;
    }

    const char* what() const noexcept override {
        return mMessage.c_str();
    }

private:
    ///@}
    ///@name Private Members
    ///@{
    std::string mMessage;
    ///@}
};

#define TIBRA_ERROR(Where) throw Exception(Where)
#define TIBRA_ERROR_IF(Where, Conditional) if(Conditional) throw Exception(Where)}

///@} // End TIBRA Classes
} // End namespace tibra

#endif // LOGGER_INCLUDE_HPP