/**
 * @file SourceLocation.cpp
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
 * @brief SourceLocation.cpp
 * @author Thomas Brandes
 * @date 01.03.2011
 */

// hpp
#include <logging/SourceLocation.hpp>

#include <string.h>

namespace logging
{

SourceLocation::SourceLocation( const char* const filename, const char* const funcname, const int line )
    : mFileName( filename ), mFuncName( funcname ), mLine( line )
{
    const char* shortname = strrchr( filename, '/' );

    if ( shortname )
    {
        mFileName = shortname + 1;
    }
}

std::ostream& operator<<( std::ostream& os, const SourceLocation& loc )
{
    os << loc.mFileName << "::" << loc.mLine << ", funct=" << loc.mFuncName;
    return os;
}

}
