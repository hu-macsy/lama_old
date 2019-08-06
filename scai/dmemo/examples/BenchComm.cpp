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

using namespace scai;
using namespace dmemo;
using namespace hmemo;

/* ----------------------------------------------------------------------------- */

void bench( const Communicator& comm, const PartitionId source, const PartitionId target )
{
    const IndexType N = 1024 * 1024;
    const IndexType NITER = 100;

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

    for ( IndexType iter = 0; iter < NITER; ++iter )
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

    double bandwidthMB = 2.0 * double( N ) * sizeof( double ) * double( NITER ) / ( 1024.0 * 1024.0 * time );

    if ( rank == source )
    {
        std::cout << "Processor " << source << " -> " << target 
                  << ": bandwidth = " << bandwidthMB << " MB/s" << std::endl;
    }
}

int main( int argc, const char* argv[] )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    common::Settings::setRank( comm->getNodeRank() );

    common::Settings::parseArgs( argc, argv );

    ContextPtr ctx = Context::getContextPtr();

    std::cout << *comm << ": run comm bench with this context: " << *ctx << std::endl;

    PartitionId size = comm->getSize();

    for ( PartitionId p1 = 0; p1 < size; p1 ++ )
    {
        for ( PartitionId p2 = 0; p2 < size; p2 ++ )
        {
            if ( p1 == p2 ) continue;

            bench( *comm, p1, p2 );
        } 
    } 
}
