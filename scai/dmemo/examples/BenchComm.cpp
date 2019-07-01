/**
 * @file BenchComm.cpp
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
 * @brief Benchmark program to measure communication bandwidth between two processes
 * @author Thomas Brandes
 * @date 24.04.2019
 */

#include <iostream>
#include <iomanip>

#include <scai/dmemo.hpp>
#include <scai/dmemo/SingleDistribution.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <mpi.h>

using namespace scai;
using namespace dmemo;
using namespace hmemo;

/* ----------------------------------------------------------------------------- */

void benchBW( const Communicator& comm, const PartitionId source, const PartitionId target, 
              IndexType m, int nIter, bool isBench )
{
    IndexType N = m / sizeof( double );

    // std::cout << "Bench comm " << source << " -> " << target << std::endl;

    PartitionId rank = comm.getRank();

    // build communication plans for a pair communication between source and target

    CommunicationPlan sendPlan;
    CommunicationPlan recvPlan;

    if ( rank == source )
    {
        sendPlan.defineBySingleEntry( N, target );
    }
    if ( rank == target )
    {
        recvPlan.defineBySingleEntry( N, source );
    }

    ContextPtr ctx = Context::getContextPtr(); // take default context, can be set by SCAI_CONTEXT

    HArray<double> sourceArray( ctx );

    if ( rank == source )
    {
        auto wArray = hostWriteOnlyAccess( sourceArray, N );

        for ( IndexType i = 0; i < N; ++i )
        {
            wArray[i] = double( i ) / double( i + 1 );
        }
    }

    HArray<double> targetArray( ctx );
    
    comm.synchronize();

    double time = common::Walltime::get();

    for ( int iter = 0; iter < nIter; ++iter )
    {
        comm.exchangeByPlan( targetArray, recvPlan, sourceArray, sendPlan );
    
        if ( rank == target )
        {
            WriteAccess<double> wA( targetArray, ctx );
        }

        comm.exchangeByPlan( sourceArray, sendPlan, targetArray, recvPlan );

        // make sure that the data is really where it is needed

        if ( rank == source )
        {
            WriteAccess<double> wA( sourceArray, ctx );
        }
    }

    comm.synchronize();

    time = common::Walltime::get() - time;

    if ( rank == source && isBench )
    {
        double bandwidthGB = 2.0 * double( N ) * sizeof( double ) * double( nIter ) / ( 1024.0 * 1024.0 * 1024.0 * time );

        std::cout << "Processor " << source << " -> " << target 
                  << ": bandwidth = " << bandwidthGB << " GB/s" << std::endl;
    }
}

/* ----------------------------------------------------------------------------- */

void benchBiBW( const Communicator& comm, const PartitionId p1, const PartitionId p2, 
                IndexType m, int nIter, bool isBench )
{
    IndexType N = m / sizeof( double );

    PartitionId rank = comm.getRank();

    // build communication plans for a pair communication between p1 and p2

    CommunicationPlan sendPlan;
    CommunicationPlan recvPlan;

    if ( rank == p1 )
    {
        sendPlan.defineBySingleEntry( N, p2 );
        recvPlan.defineBySingleEntry( N, p2 );
    }
    if ( rank == p2 )
    {
        sendPlan.defineBySingleEntry( N, p1 );
        recvPlan.defineBySingleEntry( N, p1 );
    }

    ContextPtr ctx = Context::getContextPtr(); // take default context, can be set by SCAI_CONTEXT

    HArray<double> sourceArray( ctx );

    if ( rank == p1 || rank == p2 )
    {
        auto wArray = hostWriteOnlyAccess( sourceArray, N );

        for ( IndexType i = 0; i < N; ++i )
        {
            wArray[i] = double( i ) / double( i + 1 );
        }
    }

    HArray<double> targetArray( ctx );
    
    comm.synchronize();

    double time = common::Walltime::get();

    for ( int iter = 0; iter < nIter; ++iter )
    {
        comm.exchangeByPlan( targetArray, recvPlan, sourceArray, sendPlan );
    
        if ( rank == p1 || rank == p2 )
        {
            WriteAccess<double> wA( targetArray, ctx );
        }

        targetArray.swap( sourceArray );
    }

    comm.synchronize();

    time = common::Walltime::get() - time;

    if ( rank == p1 && isBench )
    {
        double bandwidthGB = 2.0 * double( N ) * sizeof( double ) * double( nIter ) / ( 1024.0 * 1024.0 * 1024.0 * time );

        std::cout << "Processor " << p1 << " <-> " << p2 
                  << ": bidirectional bandwidth = " << bandwidthGB << " GB/s" << std::endl;
    }
}

/* ----------------------------------------------------------------------------- */

void printHelp( const char* argv[] )
{
    std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
    std::cout << "   -m  <message_size>" << std::endl;          
    std::cout << "   -x  <warmup-iterations>" << std::endl;
    std::cout << "   -i  <bench-iterations>" << std::endl;
    std::cout << "   -2  # will measure bi-directional bandwidth" << std::endl;
    exit( -1 );
}

/* ----------------------------------------------------------------------------- */

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );
 
    // parse remaining arguments

    int iarg = 1;

    // default: 256 MByte 

    IndexType msgLength = sizeof( double ) * 32 * 1024 * 1024;   //!<  message-size 

    int nIterBench = 100;
    int nIterWarmup = 2;

    bool bidirectional = false;

    while ( iarg < argc )
    {
        if ( strcmp( argv[iarg], "-2" ) == 0 )
        {
            bidirectional = true;
            iarg++;
        }
        else if ( strcmp( argv[iarg], "-m" ) == 0 )
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
        else
        {
            printHelp( argv );
        }
    }

    ContextPtr ctx = Context::getContextPtr();

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    // std::cout << *comm << ": run comm bench with this context: " << *ctx << std::endl;

    if ( comm->getRank() == 0 )
    {
        std::cout << "Bench, msg length = " << ( double( msgLength) / double( 1024 * 1024 ) ) << " MByte"
                  << ", #iter ( warmup ) = " << nIterWarmup
                  << ", #iter ( bench ) = " << nIterBench << std::endl;
    }

    PartitionId size = comm->getSize();

    for ( PartitionId p1 = 0; p1 < size; p1 ++ )
    {
        for ( PartitionId p2 = 0; p2 < size; p2 ++ )
        {
            if ( p1 == p2 ) continue;

            benchBW( *comm, p1, p2, msgLength, nIterWarmup, false ); 
        } 
    } 

    if ( comm->getRank() == 0 )
    {
        std::cout << "Warmup finished, now start timing." << std::endl;
    }

    // bench phase 

    for ( PartitionId p1 = 0; p1 < size; p1 ++ )
    {
        for ( PartitionId p2 = 0; p2 < size; p2 ++ )
        {
            if ( p1 == p2 ) continue;

            if ( bidirectional )
            {
                benchBiBW( *comm, p1, p2, msgLength, nIterBench, true );  
            }
            else
            {
                benchBW( *comm, p1, p2, msgLength, nIterBench, true );  
            }
        } 
    } 

    // MPI Finalize must be called with actual context

    comm->finalize();
}
