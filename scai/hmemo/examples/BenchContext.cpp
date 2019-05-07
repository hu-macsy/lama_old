/**
 * @file BenchContext.cpp
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
 * @brief Benchmarking of memory transfers HOST <-> context and different contextes
 * @author Thomas Brandes
 * @date 14.09.2015
 */

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>
#include <scai/logging.hpp>

using namespace std;
using namespace scai::hmemo;
using scai::common::ContextType;

/** number of iterations for the different problem cases */
static int ITER_VEC[] = { 1000, 1000,  1000,    300,     100,      50,       20,        10 };

/** number of array entries for the different problem cases 
 *
 *  Decrease values if they do not fit in memory.
 */
static int N_VEC[]    = {    1,  100, 10000, 100000, 1000000, 8000000, 16000000, 128000000 };

/**
 *  Benchmark memory transfer between host and device
 */
void benchHostGPU( ContextPtr host, ContextPtr device )
{
    int NCASES = sizeof( ITER_VEC ) / sizeof( int );

    for ( int k = 0; k < NCASES; ++k )
    {
        int ITER = ITER_VEC[k];
        int N    = N_VEC[k];

        cout << "Case " << k << ": N = " << N << ", ITER = " << ITER << endl;

        HArray<double> array( N, double( 0 ), device );

        double time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            // write only access invalidates all data
            {
                WriteOnlyAccess<double> hostWrite( array, host, N );
            }
            {
                ReadAccess<double> devRead( array, device );
            }
        }

        double t1 = scai::common::Walltime::get() - time;

        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            {
                WriteOnlyAccess<double> devWrite( array, device, N );
            }
            {
                ReadAccess<double> hostRead( array, host );
            }
        }

        double t2 = scai::common::Walltime::get() - time;
        double Bytes = static_cast<double>( N ) * sizeof( double ) * ITER;
        double GByte = 1024.0 * 1024.0 * 1024.0;
        double GBytes1 = ( Bytes / GByte ) / t1;
        double GBytes2 = ( Bytes / GByte ) / t2;
        cout << "Case " << k << ": N = " << N << ", ITER = " << ITER << ", Bytes = " << Bytes << endl;
        cout << "Host -> " << *device << " : " << t1 << " s, is " << GBytes1 << " GByte/s" << endl;
        cout << *device << " -> Host : " << t2 << " s, is " << GBytes2 << " GByte/s" << endl;
        cout << endl;
    }
}

/**
 *  Benchmark memory transfer from device gpu1 to device gpu2
 */
void benchInterGPU( ContextPtr gpu1, ContextPtr gpu2 )
{
    int NCASES = sizeof( ITER_VEC ) / sizeof( int );

    for ( int k = 0; k < NCASES; ++k )
    {
        int ITER = ITER_VEC[k];
        int N    = N_VEC[k];
        cout << "Case " << k << ": N = " << N << ", ITER = " << ITER << endl;

        HArray<double> array( N, gpu1 );  // allocate on gpu1

        {
            WriteAccess<double> wGPU2( array, gpu2 );
        }
        {
            ReadAccess<double> wGPU1( array, gpu1 );
        }

        double time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            // write access invalidates all data
            {
                WriteAccess<double> write( array, gpu1 );
            }
            // read access: transfer data from gpu1 to gpu2 
            {
                ReadAccess<double> read( array, gpu2 );
            }
        }

        time = scai::common::Walltime::get() - time;

        double Bytes = static_cast<double>( N ) * sizeof( double ) * ITER;
        double GBytes = Bytes / ( 1024.0 * 1024.0 * 1024.0 * time );
        cout << "Case " << k << ": N = " << N << ", ITER = " << ITER << ", Bytes = " << Bytes << endl;
        cout << "Transfer " << *gpu1 << " -> " << *gpu2 << ": " << GBytes << " GByte/s" << endl;
        cout << endl;
    }
}

#define MAX_DEVICES 8

int main()
{
    ContextPtr host = Context::getHostPtr();

    ContextPtr devices[MAX_DEVICES];

    int countDevices = 0;

    for ( int i = 0; i < MAX_DEVICES; ++i )
    {
        try
        {
            ContextPtr cuda = Context::getContextPtr( ContextType::CUDA, i );
            devices[ countDevices++ ] = cuda;
            std::cout << "Context " << *cuda << " created." << std::endl;
        }
        catch ( std::exception& )
        {
            std::cout << "Failed to create CUDA device " << i << std::endl;
        }
    }

    std::cout << "Benchmark host - device memory transfer for " << countDevices << " devices." << std::endl;

    for ( int i = 0; i < countDevices; ++ i )
    {
        std::cout << "CUDA bench test, context = " << *devices[i] << std::endl;
        benchHostGPU( host, devices[i] );
    }

    std::cout << "Benchmark device - device memory transfer for " << countDevices << " devices." << std::endl;

    for ( int i = 0; i < countDevices; ++i )
    {
        for ( int j = 0; j < countDevices; ++j )
        {
            if ( i == j )
            {
                continue;
            }

            benchInterGPU( devices[i], devices[j] );
        }
    }
}
