/**
 * @file SourceLocation.hpp
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
 * @brief Definition of a struct that describes a source code location.
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

#include <ostream>

#include <scai/common/config.hpp>

namespace scai
{

namespace logging
{

/** SourceLocation is a structure containing file, line and function info;
 * it specifies a location in a source file, e.g. for logging statement
 *
 */

struct COMMON_DLL_IMPORTEXPORT SourceLocation
{

    const char* mFileName; //!< Name of the source file
    const char* mFuncName; //!< Name of the function
    int mLine; //!< Line number of location in source file

    /** Constructor of a location */

    SourceLocation( const char* const filename, const char* const funcname, const int line );

};

} /* end namespace logging */

} /* end namespace scai */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& os, const scai::logging::SourceLocation& loc );

#if !defined(LOG4SCAI_LOCATION)
#if defined(_MSC_VER)
#if _MSC_VER >= 1300
#define __LOG4SCAI_FUNC__ __FUNCSIG__
#endif
#else
#if defined(__GNUC__)
// this macro gives the whole signature
// #define __LOG4SCAI_FUNC__ __PRETTY_FUNCTION__
// this macro is standardized
#define __LOG4SCAI_FUNC__ __func__
#endif
#endif

#if !defined(__LOG4SCAI_FUNC__)
#define __LOG4SCAI_FUNC__ ""
#endif

#define LOG4SCAI_LOCATION scai::logging::SourceLocation(__FILE__, __LOG4SCAI_FUNC__, __LINE__)
#endif
