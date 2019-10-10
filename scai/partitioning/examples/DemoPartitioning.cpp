/**
 * @file DemoPartitioning.cpp
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
 * @brief Demo program for using partitioning
 * @author Thomas Brandes
 * @date 23.08.2017
 */

#include <scai/partitioning.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/freeFunction.hpp>
#include <scai/lama/io/PartitionIO.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/RedistributePlan.hpp>

#define HOST_PRINT( comm, msg )            \
    {                                      \
        int myRank = comm->getRank();      \
        if ( myRank == 0 )                 \
        {                                  \
            std::cout << msg << std::endl; \
        }                                  \
    }

using namespace scai;
using namespace dmemo;
using namespace lama;
using namespace partitioning;

typedef DefaultReal ValueType;

template<typename ValueType>
static void statistic( const SparseMatrix<ValueType>& matrix, const char* name )
{
    auto plan = matrix.getHaloExchangePlan();

    auto comm = matrix.getRowDistribution().getCommunicatorPtr();

    IndexType size = plan.getHaloSize();

    IndexType min   = comm->min( size );
    IndexType max   = comm->max( size );
    IndexType total = comm->sum( size );

    HOST_PRINT( comm, "Communication[" << name << "]: min = " << min << ", max = " << max << " avg = " << ( total / comm->getSize()  ) )

    size = matrix.getRowDistribution().getLocalSize();

    min   = comm->min( size );
    max   = comm->max( size );
    total = comm->sum( size );

    HOST_PRINT( comm, "#rows[" << name << "]: min = " << min << ", max = " << max << " avg = " << ( total / comm->getSize()  ) )
 
    size = matrix.getLocalStorage().getNumValues() + matrix.getHaloStorage().getNumValues();

    min   = comm->min( size );
    max   = comm->max( size );
    total = comm->sum( size );

    HOST_PRINT( comm, "#nnz[" << name << "]: min = " << min << ", max = " << max << " avg = " << ( total / comm->getSize()  ) )
}

int main( int narg, const char* argv[] )
{
    if ( narg < 2 )
    {
        std::cout << "call " << argv[0] << " <matrixFileName>" << std::endl;
    }

    std::string fileName = argv[1];

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    SCAI_REGION( "main.DemoPartitioning" )

    auto A = read<CSRSparseMatrix<ValueType>>( fileName );

    HOST_PRINT( comm, ": Read matrix A = " << A )

    auto blockDist = blockDistribution( A.getNumRows(), comm );
    auto Ab = distribute<CSRSparseMatrix<ValueType>>( A, blockDist, blockDist );

    statistic( Ab, "block" );

    const char PARMETIS_ID[] = "PARMETIS";

    if ( Partitioning::canCreate( PARMETIS_ID ) )
    {
        const char DIST_FILE_NAME[] = "distParMetis.txt";

        Ab.redistribute( blockDist, noDistribution( Ab.getNumColumns() ) );

        PartitioningPtr thePartitioning = Partitioning::create( PARMETIS_ID );

        std::cout << *comm << ":Parmetis partitioning, A = " << Ab << std::endl;

        scai::hmemo::HArray<IndexType> newLocalOwners;

        float procWeight = 1.0f;

        thePartitioning->squarePartitioning(newLocalOwners, Ab, procWeight);

        auto plan = dmemo::redistributePlanByNewOwners(newLocalOwners, blockDist);

        auto dist = plan.getTargetDistributionPtr();

        Ab.redistribute( dist, dist );
        statistic( Ab, PARMETIS_ID );

        HOST_PRINT( comm, "Write distribution to file " << DIST_FILE_NAME )
        PartitionIO::write( *dist, DIST_FILE_NAME );
    }

    const char METIS_ID[] = "METIS";

    if ( Partitioning::canCreate( METIS_ID ) )
    {
        const char DIST_FILE_NAME[] = "distMetis.txt";

        PartitioningPtr thePartitioning = Partitioning::create( "METIS" );

        DistributionPtr dist( thePartitioning->partitionIt( comm, A, 1.0f ) );

        A.redistribute( dist, dist );

        statistic( A, METIS_ID );

        HOST_PRINT( comm, "Write distribution to file " << DIST_FILE_NAME )

        PartitionIO::write( *dist, DIST_FILE_NAME );
    }
}
