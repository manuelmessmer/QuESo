// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef LOGGER_INCLUDE_HPP
#define LOGGER_INCLUDE_HPP

//// STL includes
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace queso {

///@name QuESo Classes
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
        AppendMessage( " -- " );

    }

    /// Destructor
    ~Logger() {
        std::cout << mMessage;
    }

    ///@}
    ///@name Operations
    ///@{

    Logger& operator << (std::ostream& (*pf)(std::ostream&)) {
       	std::stringstream buffer;
		pf(buffer);

        AppendMessage(buffer.str());

        return *this;
    }

    template<class StreamValueType>
    Logger& operator << (StreamValueType const& rValue)
    {
        std::stringstream buffer;
        buffer << std::setprecision(10) << rValue;

        AppendMessage(buffer.str());

        return *this;
    }

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

#define QuESo_INFO Logger()
#define QuESo_INFO_IF(Conditional) if(Conditional) Logger()

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
        AppendMessage( "\n Where: " );
        AppendMessage( rWhere );
        AppendMessage( " -- What: " );
    }

    ///@}
    ///@name Operations
    ///@{

    template<class StreamValueType>
    Exception& operator << (StreamValueType const& rValue)
    {
        std::stringstream buffer;
        buffer << rValue;

        AppendMessage(buffer.str());

        return *this;
    }

    Exception& operator << (std::ostream& (*pf)(std::ostream&)) {
       	std::stringstream buffer;
		pf(buffer);

        AppendMessage(buffer.str());

        return *this;
    }

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

#define QuESo_ERROR(Where) throw Exception(Where)
#define QuESo_ERROR_IF(Where, Conditional) if(Conditional) throw Exception(Where)

///@} // End QuESo Classes
} // End namespace queso

#endif // LOGGER_INCLUDE_HPP