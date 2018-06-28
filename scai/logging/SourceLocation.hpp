/**
 * @file SourceLocation.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Definition of a struct that describes a source code location.
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

// for dll import
#include <scai/common/config.hpp>

#include <ostream>

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

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& os, const SourceLocation& loc );

} /* end namespace logging */

} /* end namespace scai */

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
