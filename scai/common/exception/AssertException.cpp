/*
 * AssertException.cpp
 *
 *  Created on: Aug 31, 2015
 *      Author: eschricker
 */

// hpp
#include <scai/common/exception/AssertException.hpp>

namespace scai
{

namespace common
{

AssertException::AssertException() :
		mMessage( "AssertException" )
{
}

AssertException::AssertException( const std::string& message ) :
		mMessage( message )
{
	mMessage += "@AssertException";
}

AssertException::~AssertException() throw ()
{
}

const char* AssertException::what() const throw ()
{
    return mMessage.c_str();
}

} /* end namespace common */

} /* end namespace scai */


