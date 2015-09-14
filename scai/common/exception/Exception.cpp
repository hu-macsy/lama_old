/**
 * @file Exception.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Implementation of methods for class Exception.
 * @author Jiri Kraus
 * @date 01.03.2011
 */

// hpp
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

void Exception::addCallStack( std::ostringstream& output )
{
    const size_t maxDepth = 20;

    void *stackAddrs[maxDepth];

    size_t stackDepth = backtrace( stackAddrs, maxDepth );
    char** stackStrings = backtrace_symbols( stackAddrs, stackDepth );

    for( size_t i = 1; i < stackDepth; i++ )
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
    for( char *j = fName.data(); *j; ++j )
    {
        if( *j == '(' )
        {
            begin = j;
        }
        else if( *j == '+' )
        {
            end = j;
        }
    }

    if( begin && end )
    {
        begin++;
        *end = '\0';
        // found our mangled name, now in [begin, end)

        int status;
        char *ret = abi::__cxa_demangle( begin, 0, 0, &status );

        if( status == 0 )
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

        if( ret )
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

void Exception::addCallStack( std::ostringstream& )
{
}

std::string Exception::demangle( const char* functionName )
{
    return functionName;
}

#endif

}  /* end namespace common */

} /* end namespace scai */
