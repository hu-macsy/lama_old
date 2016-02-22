/*
 * NotSupportedValueTypeException.cpp
 *
 *  Created on: Aug 31, 2015
 *      Author: eschricker
 */

// hpp
#include <scai/common/exception/NotSupportedValueTypeException.hpp>

namespace scai
{

namespace common
{

NotSupportedValueTypeException::NotSupportedValueTypeException() :
		mMessage( "NotSupportedValueTypeException" )
{
}

NotSupportedValueTypeException::NotSupportedValueTypeException( const std::string& message ) :
		mMessage( message )
{
	mMessage += "@NotSupportedValueTypeException";
}

NotSupportedValueTypeException::~NotSupportedValueTypeException() throw ()
{
}

const char* NotSupportedValueTypeException::what() const throw ()
{
    return mMessage.c_str();
}

} /* end namespace common */

} /* end namespace scai */


