/**
 * @file MPIUtils.hpp
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
 * @brief Utility macros for MPI calls
 * @author Thomas Brandes
 * @date 23.03.2011
 */

#pragma once

#include <mpi.h> //Intel MPI need mpi.h to be included before stdio.h so this header comes first

// local library
#include <scai/dmemo/mpi/MPIException.hpp>

// std
#include <sstream>
#include <iostream>

#ifdef SCAI_CHECK_ASSERTS

#define SCAI_MPICALL( logger, exp, msg)                                             \
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
            common::Exception::addCallStack( errorStr );                            \
            throw MPIException( errorStr.str(), status );                           \
        }                                                                           \
    }

#else

#define SCAI_MPICALL( logger, exp, msg) exp;

#endif // SCAI_CHECK_ASSERTS
