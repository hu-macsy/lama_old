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
 * @brief Demo program for using partitioning with weights
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

template<typename ValueType>
static void statistic( const DenseVector<ValueType>& vector, const char* name )
{
    auto comm = vector.getDistribution().getCommunicatorPtr();

    ValueType localSum = utilskernel::HArrayUtils::sum( vector.getLocalValues() );

    ValueType min   = comm->min( localSum );
    ValueType max   = comm->max( localSum );
    ValueType total = comm->sum( localSum );

    HOST_PRINT( comm, "Weights[" << name << "]: min = " << min << ", max = " << max << " avg = " << ( total / comm->getSize()  ) )
}

void printHelp( CommunicatorPtr comm, const char* cmd )
{
    HOST_PRINT( comm, cmd << " <matrixFileName> [-w <weightsFileName>] [-p <partitioningTool>]" )
}

int main( int narg, const char* argv[] )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( narg < 2 )
    {
        HOST_PRINT( comm, "not enough arguments" )
        printHelp( comm, argv[0] );
        return 1;
    }

    bool hasVertexWeights = false;

    std::string matrixFileName = argv[1];
    std::string weightsFileName;                 // default is no weights
    std::string PARTITIONING_ID = "PARMETIS";    // default partitioning tool

    int k = 2;

    while ( k < narg )
    {
        if ( strcmp( argv[k], "-w" ) == 0 )
        {
            if ( k + 1 < narg )
            {
                weightsFileName = argv[k + 1];
                hasVertexWeights = true;
            }
            else
            {
                HOST_PRINT( comm, "-w must be followed by filename" )
                printHelp( comm, argv[0] );
                return 1;
            }

            k += 2;
        }
        else if ( strcmp( argv[k], "-p" ) == 0 )
        {
            if ( k + 1 < narg )
            {
                PARTITIONING_ID = argv[k + 1];
            }
            else
            {
                HOST_PRINT( comm, "-p must be followed by partitioning tool" )
                printHelp( comm, argv[0] );
                return 1;
            }

            k += 2;
        }
        else
        {
            HOST_PRINT( comm, "illegal arg " << argv[k] )
            printHelp( comm, argv[0] );
            return 1;
        }
    }
    
    SCAI_REGION( "main.DemoPartitioning" )

    auto A = read<CSRSparseMatrix<ValueType>>( matrixFileName );

    SCAI_ASSERT_EQ_ERROR( A.getNumRows(), A.getNumColumns(), "not square matrix" )

    HOST_PRINT( comm, "Read matrix A = " << A )

    auto oldDist = blockDistribution( A.getNumRows(), comm );
    A.redistribute( oldDist, oldDist );

    statistic( A, "BLOCK" );

    DenseVector<float> W;
  
    if ( hasVertexWeights )
    {
        W.readFromFile( weightsFileName );
        HOST_PRINT( comm, "Read weights = " << W );
        SCAI_ASSERT_EQ_ERROR( W.size(), A.getNumRows(), "#weights does not match matrix size" )
        W.redistribute( oldDist );
        statistic( W, "block" );
    }

    if ( Partitioning::canCreate( PARTITIONING_ID.c_str() ) )
    {
        std::string DIST_FILE_NAME = "dist_";
        DIST_FILE_NAME += PARTITIONING_ID;
        DIST_FILE_NAME += ".txt";

        PartitioningPtr thePartitioning = Partitioning::create( PARTITIONING_ID.c_str() );

        scai::hmemo::HArray<IndexType> newLocalOwners;

        float procWeight = 1.0f;

        if ( W.size() >  0 )
        {
            HOST_PRINT( comm, PARTITIONING_ID << ": square partitioning with vertex weights" )
            thePartitioning->squarePartitioningW( newLocalOwners, A, W.getLocalValues(), procWeight);
        }
        else
        {
            HOST_PRINT( comm, PARTITIONING_ID << ": square partitioning, NO vertex weights" )
            thePartitioning->squarePartitioning( newLocalOwners, A, procWeight);
        }

        HOST_PRINT( comm, PARTITIONING_ID << ": new owners computed, now redistribute" )

        auto plan = dmemo::redistributePlanByNewOwners( newLocalOwners, oldDist );

        auto dist = plan.getTargetDistributionPtr();

        A.redistribute( dist, dist );
        statistic( A, PARTITIONING_ID.c_str() );

        if ( hasVertexWeights )
        {
            W.redistribute( dist );
            statistic( W, PARTITIONING_ID.c_str() );
        }

        DenseVector<IndexType> newOwners( oldDist, newLocalOwners );

        newOwners += 3;

        HOST_PRINT( comm, "Write distribution to file " << DIST_FILE_NAME )

        newOwners.writeToFile( DIST_FILE_NAME );
    }
    else
    { 
        HOST_PRINT( comm, "ERROR: partitioning " << PARTITIONING_ID.c_str() << " not supported"  )
        return 1;
    }
}
