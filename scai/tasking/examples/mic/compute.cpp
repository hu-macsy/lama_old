/**
 * @file tasking/examples/mic/compute.cpp
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
 * @brief ToDo: Missing description in ./tasking/examples/mic/compute.cpp
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <cstdlib>

__declspec( target(mic, host) )
void compute( size_t NSIZE )
{
    double* a = new double[ NSIZE * NSIZE ];
    double* b = new double[ NSIZE * NSIZE ];
    double* c = new double[ NSIZE * NSIZE ];

    for ( int i = 0; i < NSIZE; ++i )
    {
        for ( int j = 0; j < NSIZE; ++j )
        {
            a[ i * NSIZE + j ] = 1;
            b[ i * NSIZE + j ] = 1;
            c[ i * NSIZE + j ] = 0;
        }
    }

    #pragma omp parallel for

    for ( int i = 0; i < NSIZE; ++i )
    {
        for ( int j = 0; j < NSIZE; ++j )
        {
            double tmp = 0;

            for ( int k = 0; k < NSIZE; ++k )
            {
                tmp += a[ i * NSIZE + k ] * b[ k * NSIZE + j ];
            }

            c[ i * NSIZE + j ] = tmp;
        }
    }

    delete[] a;
    delete[] b;
    delete[] c;
}

