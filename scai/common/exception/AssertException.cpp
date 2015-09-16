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

AssertException::AssertException()
{
}

AssertException::AssertException( const std::string& message ) : Exception( message )
{
}

AssertException::~AssertException() throw ()
{
}

} /* end namespace common */

} /* end namespace scai */


