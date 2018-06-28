/**
 * @file tracing/examples/TraceCalls.cpp
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
 * @brief ToDo: Missing description in ./tracing/examples/TraceCalls.cpp
 * @author Thomas Brandes
 * @date 18.06.2015
 */

#include <scai/tracing.hpp>
#include <scai/common/Walltime.hpp>

#include <cstdio>

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
    double time = scai::common::Walltime::get();
    int X = 0;
    SCAI_REGION( "main" )

    for ( int i = 0; i < 10000; ++i )
    {
        for ( int j = 0; j < 30; ++ j )
        {
            subA( X );
        }

        for ( int j = 0; j < 20; ++ j )
        {
            subB( X );
        }
    }

    time = scai::common::Walltime::get() - time;
    printf( "X = %d, number of calls, time = %f s\n", X, time );
}
