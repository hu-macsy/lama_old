/**
 * @file scai/tracing/test/simpleTracing.cpp
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
 * @brief simple executable that is used by the tracing tests
 * @author Jan Ecker
 * @date 18.02.2016
 */

#include <scai/tracing.hpp>
#include <scai/logging.hpp>
#include <cstdio>
#include <scai/common/OpenMP.hpp>

void subA( int& X )
{
    SCAI_REGION( "A" )
    ++X;
}

void subB( int& X )
{
    SCAI_REGION( "B" )
    X++;
}

int main()
{
    int X = 0;
    SCAI_REGION( "main" )
    #pragma omp parallel for

    for ( int i = 0; i < 10000; ++i )
    {
#ifndef UNNAMED_THREADS
        SCAI_LOG_THREAD( omp_get_thread_num() );
#endif
        SCAI_REGION_START( "main.loopA" )

        for ( int j = 0; j < 30; ++ j )
        {
            subA( X );
        }

        SCAI_REGION_END( "main.loopA" )
        SCAI_REGION_START( "main.loopB" )

        for ( int j = 0; j < 20; ++ j )
        {
            subB( X );
        }

        SCAI_REGION_END( "main.loopB" )
    }
}
