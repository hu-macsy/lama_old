/**
 * @file Exception.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Implementation of methods for class Exception.
 * @author Jiri Kraus, Thomas Brandes
 * @date 01.03.2011
 */

#include <scai/common/exception/Exception.hpp>

// std
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <sstream>

#ifndef _WIN32
#include <execinfo.h>
#endif //_WIND32

#ifdef __GNUC__
#include <cxxabi.h>
#include <cstring>
#endif // __GNUC__

namespace scai
{

namespace common
{

Exception::Exception()
{
}

Exception::Exception( const std::string& message ) : mMessage( message )
{
}

Exception::~Exception() throw ()
{
}

const char* Exception::what() const throw ()
{
    return mMessage.c_str();
}

#ifdef __GNUC__

void Exception::addCallStack( std::ostream& output )
{
    const size_t maxDepth = 20;
    void* stackAddrs[maxDepth];
    size_t stackDepth = backtrace( stackAddrs, maxDepth );
    char** stackStrings = backtrace_symbols( stackAddrs, stackDepth );

    for ( size_t i = 1; i < stackDepth; i++ )
    {
        output << "    stack[" << i << "] : " << demangle( stackStrings[i] ) << std::endl;
    }

    free( stackStrings ); // malloc()ed by backtrace_symbols
}

std::string Exception::demangle( const char* functionName )
{
    // We need a copy of functionName for demangling it.
    std::vector<char> fName( strlen( functionName ) + 1 );
    strcpy( fName.data(), functionName );
    std::string demangledString;
    char* begin = 0;
    char* end = 0;

    // find the parentheses and address offset surrounding the mangled name
    for ( char* j = fName.data(); *j; ++j )
    {
        if ( *j == '(' )
        {
            begin = j;
        }
        else if ( *j == '+' )
        {
            end = j;
        }
    }

    if ( begin && end )
    {
        begin++;
        *end = '\0';
        // found our mangled name, now in [begin, end)
        int status;
        char* ret = abi::__cxa_demangle( begin, 0, 0, &status );

        if ( status == 0 )
        {
            // return value may be a realloc() of the input
            demangledString = ret;
        }
        else
        {
            // demangling failed, just pretend it's a C function with no args
            //std::strncpy(function, begin, sz);
            //std::strncat(function, "()", sz);
            //function[sz-1] = ' ';
            demangledString = functionName;
        }

        if ( ret )
        {
            free( ret );
        }
    }
    else
    {
        // didn't find the mangled name, just print the whole line
        demangledString = functionName;
    }

    return demangledString;
}

#else

void Exception::addCallStack( std::ostream& )
{
}

std::string Exception::demangle( const char* functionName )
{
    return functionName;
}

#endif

}  /* end namespace common */

} /* end namespace scai */
