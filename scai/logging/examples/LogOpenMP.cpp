/**
 * @file LogOpenMP.cpp
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
 * @brief Example program for logging in OpenMP Threads
 *        /
 * @author Thomas Brandes
 * @date 21.08.2015
 */

#include <scai/logging.hpp>
#include <omp.h>

SCAI_LOG_DEF_LOGGER( myLogger, "LogOpenMP" )

int main( int, char** )
{
    SCAI_LOG_THREAD( "main" )
    // Here we name the OpenMP Threads
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        if ( thread_id > 0 )
        {
            SCAI_LOG_THREAD( "OMP_Thread_" << thread_id )
        }
    }
    const int N = 20;
    #pragma omp parallel for schedule( static, 2 )

    for ( int i = 0; i < N; ++i )
    {
        // logging per thread shows exactly which thread executes which iteration
        SCAI_LOG_INFO( myLogger, "executes iteration " << i << " of " << N )
    }
}
