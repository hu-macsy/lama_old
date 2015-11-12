#pragma once

// base class
#include <scai/common/exception/Exception.hpp>

// std
#include <string>

namespace scai
{

namespace common
{

class COMMON_DLL_IMPORTEXPORT NotSupportedValueTypeException : public Exception
{
public:
    /**
     * @brief The default constructor creates an NotSupportedValueTypeException with no message.
     */
	NotSupportedValueTypeException( );

    /**
     * @brief This constructor creates an NotSupportedValueTypeException with the passed message.
     *
     * @param[in] message  the message to assign to this.
     */
	NotSupportedValueTypeException( const std::string& message );

    /**
     * @brief The destructor destroys this NotSupportedValueTypeException.
     */
    virtual ~NotSupportedValueTypeException() throw ();

    virtual const char* what() const throw();

protected:

    std::string mMessage;
};

} /* end namespace common */

} /* end namespace scai */

