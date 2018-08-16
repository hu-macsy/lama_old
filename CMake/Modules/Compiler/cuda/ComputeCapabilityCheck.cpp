/**
 * @file ComputeCapabilityCheck.cpp
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
 * @brief ComputeCapabilityCheck.cpp
 * @author Jan Ecker
 * @date 10.01.2014
 */

#include <stdio.h>
#include <cuda_runtime.h>
#include <stdlib.h>

void printVersion ( int major, int minor )
{
    printf( "%d%d", major, minor );
}

int main()
{
    int deviceCount;
    struct cudaDeviceProp properties;

    if ( cudaGetDeviceCount( &deviceCount ) != cudaSuccess )
    {
        return 1;
    }

    if ( getenv( "SCAI_DEVICE" ) )
    {
        char *pEnd;
        long int cuda_device = strtol( getenv( "SCAI_DEVICE" ), &pEnd, 10 );

        if ( deviceCount > cuda_device )
        {
            cudaGetDeviceProperties( &properties, cuda_device );
            printVersion( properties.major, properties.minor );
            return 0;  // success
        }
    }
    else
    {
        int max_major = -1;
        int max_minor = -1;
        for ( int device = 0; device < deviceCount; ++device )
        {
            int major = -1, minor = -1;
            cudaGetDeviceProperties( &properties, device );
            if ( properties.major != 9999 ) // 9999 means emulation only
            { 
                major = properties.major;
                minor = properties.minor;
                if ( major > max_major )
                {
                    max_major = major;
                    max_minor = minor;
                }
                else if ( minor > max_minor )
                {
                    max_minor = minor;
                }
            }
        }
        if ( max_major != -1 && max_minor != -1 )
        {
            printVersion( max_major, max_minor );
            return 0; // success
        }
    }

    return 1; // failure
}
