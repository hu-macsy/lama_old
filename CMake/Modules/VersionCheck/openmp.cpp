/**
 * @file VersionCheck/openmp.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Checking the OpenMP version
 * @author Eric Schricker
 * @date 07.04.2016
 */

#include <stdio.h>
#include <omp.h>

#ifndef _OPENMP
    #define _OPENMP -1
#endif

int main(int argc, char *argv[])
{
    double version = 0.0;

    if ( _OPENMP >= 200505 )
    {
        version = 2.5;
    }

    if ( _OPENMP >= 200805 )
    {
        version = 3.0;
    }

    if ( _OPENMP >= 201107 )
    {
        version = 3.1;
    }

    if ( _OPENMP >= 201307 )
    {
        version = 4.0;
    }

    if ( _OPENMP >= 201511 )
    {
        version = 4.5;
    }

    if ( version == 0.0 )
    {
        fprintf( stderr, "Unknown OpenMP-Version, _OPENMP = %d\n", _OPENMP );
        return -1;
    }
    
    printf("%3.1f", version );
}
