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
 * @brief Benchmarking of memory transfers HOST <-> CUDA devices and
 *        between the different CUDA devices
 * @author Thomas Brandes
 * @date 14.09.2015
 */

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>
#include <scai/logging.hpp>

using namespace scai;
using namespace hmemo;

using namespace std;

using scai::common::ContextType;

/**
 *  Benchmark memory transfer from one context to another context 
 *
 *  @param[in] ctx1 source context
 *  @param[in] ctx2 target context
 *  @param[in] N    number of (double) array elements to transfer
 *  @param[in] nIterWarmup  number of transfers before timing is done
 *  @param[in] nIterBench   number of transfers that are measured
 */
void benchContextTransfer( ContextPtr ctx1, ContextPtr ctx2, const IndexType N,
                           const IndexType nIterWarmup, const IndexType nIterBench)
{
    // first touch on CUDA context to allow fast memory transfer

    ContextPtr firstContext = ctx1->getType() == ContextType::Host ? ctx2 : ctx1;

    HArray<double> array( N, firstContext );  // allocate on CUDA / CUDAHost context

    for ( int i = 0; i < nIterWarmup; ++i )
    {
        {
            // write access on ctx1 invalidates data on ctx2
            WriteAccess<double> wGPU2( array, ctx1 );
        }
        {
            // read access implies transfer ctx1 -> ctx2
            ReadAccess<double> wGPU1( array, ctx2 );
        }
    }

    double time = scai::common::Walltime::get();

    for ( int i = 0; i < nIterBench; ++i )
    {
        {
            WriteAccess<double> write( array, ctx1 );
        }
        // read access: transfer data from ctx1 to ctx2 
        {
            ReadAccess<double> read( array, ctx2 );
        }
    }

    time = scai::common::Walltime::get() - time;

    double msgLength  = static_cast<double>( N ) * sizeof( double );
    double totalBytes = msgLength * nIterBench;
    double GByte = 1024.0 * 1024.0 * 1024.0;
    double GBytes = ( totalBytes / GByte ) / time;
    cout << "Transfer " << *ctx1 << " -> " << *ctx2 << ": " << GBytes << " GByte/s" << endl;
}

#define MAX_DEVICES 8

/* ----------------------------------------------------------------------------- */

void printHelp( const char* argv[] )
{
    std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
    std::cout << "   -m  <message_size>" << std::endl;
    std::cout << "   -x  <warmup-iterations>" << std::endl;
    std::cout << "   -i  <bench-iterations>" << std::endl;
    std::cout << "   -d  <number-of-devices>, max = " << MAX_DEVICES << std::endl;
    exit( -1 );
}

/* ----------------------------------------------------------------------------- */

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    // parse remaining arguments

    int iarg = 1;

    // default: 256 MByte 

    IndexType msgLength = 256 * 1024 * 1024;   //!<  message-size 

    int nIterBench = 100;
    int nIterWarmup = 2;

    int nTryDevices = MAX_DEVICES;

    while ( iarg < argc )
    {
        if ( strcmp( argv[iarg], "-m" ) == 0 )
        {
             iarg++;

             if ( iarg  < argc )
             {
                 // take lenght in MByte

                 msgLength = atoi( argv[iarg] ) * 1024 * 1024;
                 iarg++;
             }
        }
        else if ( strcmp( argv[iarg], "-x" ) == 0 )
        {
             iarg++;

             if ( iarg  < argc )
             {
                 nIterWarmup = atoi( argv[iarg] );
                 iarg++;
             }
        }
        else if ( strcmp( argv[iarg], "-i" ) == 0 )
        {
             iarg++;

             if ( iarg  < argc )
             {
                 nIterBench = atoi( argv[iarg] );
                 iarg++;
             }
        }
        else if ( strcmp( argv[iarg], "-d" ) == 0 )
        {
             iarg++;

             if ( iarg  < argc )
             {
                 nTryDevices = atoi( argv[iarg] );

                 if ( nTryDevices > MAX_DEVICES )
                 {
                     nTryDevices = MAX_DEVICES;
                 }

                 iarg++;
             }
        }
        else
        {
            printHelp( argv );
        }
    }

    std::cout << "Benchmark context transfer" << std::endl;
    std::cout << "==========================" << std::endl;

    std::cout << "msg size = " << ( double( msgLength ) / double( 1024 * 1024 ) ) << " MByte"
              << ", #iter warmup = " << nIterWarmup
              << ", #iter bench = " << nIterBench << std::endl;

    std::cout << std::endl;

    IndexType N = msgLength / sizeof( double );

    ContextPtr host = Context::getHostPtr();

    ContextPtr devices[MAX_DEVICES];

    int countDevices = 0;

    for ( int i = 0; i < nTryDevices; ++i )
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

    std::cout << "Benchmark host -> device memory transfer for " << countDevices << " devices." << std::endl;

    for ( int i = 0; i < countDevices; ++ i )
    {
        benchContextTransfer( host, devices[i], N, nIterWarmup, nIterBench );
    }

    std::cout << "Benchmark device -> host memory transfer for " << countDevices << " devices." << std::endl;

    for ( int i = 0; i < countDevices; ++ i )
    {
        benchContextTransfer( devices[i], host, N, nIterWarmup, nIterBench );
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

            benchContextTransfer( devices[i], devices[j], N, nIterWarmup, nIterBench );
        }
    }
}
