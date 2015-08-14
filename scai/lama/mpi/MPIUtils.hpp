/**
 * @file MPIUtils.hpp
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
 * @brief Utility macros for MPI calls
 * @author Thomas Brandes
 * @date 23.03.2011
 * @since 1.0.0
 */
#ifndef LAMA_MPIUTILS_HPP_
#define LAMA_MPIUTILS_HPP_

#include <scai/lama/mpi/MPIException.hpp>

#include <sstream>
#include <mpi.h>

#ifdef LAMA_CHECK_ASSERTS

#define LAMA_MPICALL( logger, exp, msg)                                             \
    {                                                                               \
        SCAI_LOG_TRACE( logger, "MPI call " << msg );                               \
        int status = exp;                                                           \
        SCAI_LOG_TRACE( logger, "MPI call " << msg  << ", status = " << status );   \
        if ( status != MPI_SUCCESS )                                                \
        {                                                                           \
            std::ostringstream errorStr;                                            \
            errorStr << "MPI error in line " << __LINE__ ;                          \
            errorStr << " of file " << __FILE__ << ": ";                            \
            errorStr << msg<< "\n";                                                 \
            common::Exception::addCallStack( errorStr );                           \
            fprintf( stderr, "%s\n", errorStr.str().c_str() );                      \
            throw MPIException( errorStr.str(), status );                           \
        }                                                                           \
    }

#else

#define LAMA_MPICALL( logger, exp, msg) exp;

#endif // LAMA_CHECK_ASSERTS

#endif // LAMA_MPIUTILS_HPP_
