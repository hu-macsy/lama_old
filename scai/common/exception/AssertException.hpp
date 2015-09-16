#pragma once

// base class
#include <scai/common/exception/Exception.hpp>

// std
#include <string>

namespace scai
{

namespace common
{

class COMMON_DLL_IMPORTEXPORT AssertException : public Exception
{
public:
    /**
     * @brief The default constructor creates an AssertException with no message.
     */
	AssertException( );

    /**
     * @brief This constructor creates an AssertException with the passed message.
     *
     * @param[in] message  the message to assign to this.
     */
	AssertException( const std::string& message );

    /**
     * @brief The destructor destroys this AssertException.
     */
    virtual ~AssertException() throw ();
};

} /* end namespace common */

} /* end namespace scai */

