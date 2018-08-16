/**
 * @file PardisoError.hpp
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
 * @brief Error handling for Pardiso
 * @author Lauretta Schubert
 * @date 18.07.2016
 */

#pragma once

namespace scai
{

namespace sparsekernel
{

/** Function that translates enum cublasStatus to strings. */

COMMON_DLL_IMPORTEXPORT const char* pardisoErrorString( int error );

} /* end namespace sparsekernel */

} /* end namespace scai */

#define SCAI_PARDISO_ERROR_CHECK( error, msg )                                          \
    {                                                                                   \
        if ( error != 0 )                                                               \
        {                                                                               \
            std::ostringstream errorStr;                                                \
            errorStr << "Pardiso error in line " << __LINE__;                           \
            errorStr << " of file " << __FILE__ << std::endl;                           \
            errorStr << "  Msg  : " << msg << std::endl;                                \
            errorStr << "  Error: ";                                                    \
            errorStr << scai::sparsekernel::pardisoErrorString( error );                \
            errorStr << ", pardisoError = " << error << "\n";                           \
            scai::common::Exception::addCallStack( errorStr );                          \
            throw scai::common::Exception( errorStr.str() );                            \
        }                                                                               \
    }
