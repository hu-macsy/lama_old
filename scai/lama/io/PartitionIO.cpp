/**
 * @file lama/io/PartitionIO.cpp
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
 * @brief IO support for partitioned read/write of vectors, matrices, distributions
 * @author Thomas Brandes
 * @date 19.06.2016
 */

#include <scai/lama/io/PartitionIO.hpp>
#include <scai/lama/io/FileIO.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>

#include <iostream>

using namespace std;

namespace scai
{

using namespace dmemo;
using utilskernel::HArrayUtils;

namespace lama
{

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( PartitionIO::logger, "PartitionIO" )

/* --------------------------------------------------------------------------------- */

void PartitionIO::getPartitionFileName( string& fileName, bool& isPartitioned, const PartitionId rank, const PartitionId size )
{
    size_t pos = fileName.find( "%r" );

    isPartitioned = false;

    if ( pos == string::npos )
    {
        return;
    }

    ostringstream rankStr;

    if ( size > 1 )
    {
        rankStr << rank << "." << size;
        isPartitioned = true;
    }

    fileName.replace( pos, 2, rankStr.str() );
}

/* --------------------------------------------------------------------------------- */

void PartitionIO::getPartitionFileName( string& fileName, bool& isPartitioned, const Communicator& comm )
{
    size_t pos = fileName.find( "%r" );

    isPartitioned = false;

    if ( pos == string::npos )
    {
        return;
    }
    else
    {
        ostringstream rankStr;

        if ( comm.getSize() > 1 )
        {
            rankStr << comm.getRank() << "." << comm.getSize();
            isPartitioned = true;
        }

        fileName.replace( pos, 2, rankStr.str() );
    }
}

/* --------------------------------------------------------------------------------- */

void PartitionIO::getSingleFileName( string& fileName )
{
    size_t pos = fileName.find( "%r" );

    if ( pos != string::npos )
    {
        fileName.replace( pos, 2, "" );
    }
}

/* --------------------------------------------------------------------------------- */

bool PartitionIO::isPartitionFileName( const string& fileName )
{
    size_t pos = fileName.find( "%r" );

    if ( pos == string::npos )
    {
        return false;
    }
    else
    {
        return true;
    }
}

/* --------------------------------------------------------------------------------- */

int PartitionIO::removeFile( const std::string& fileName, const dmemo::Communicator& comm )
{
    int rc = 0;

    std::string pFileName = fileName;

    bool isPartitioned;

    getPartitionFileName( pFileName, isPartitioned, comm );

    if ( isPartitioned )
    {
        // partitoned file name, each processor deletes its part, sync results

        rc = FileIO::removeFile( pFileName );

        bool okay = rc == 0;

        okay = comm.all( okay );

        if ( !okay )
        {
            rc = -1;
        }
    }
    else
    {
        // serial file name, only root processor deletes and bcasts the result

        const PartitionId root = 0;

        if ( comm.getRank() == root )
        {
            rc = FileIO::removeFile( pFileName );
        }

        comm.bcast( &rc, 1, root );
    }

    return rc;
}

/* --------------------------------------------------------------------------------- */

bool PartitionIO::fileExists( const std::string& fileName, const dmemo::Communicator& comm )
{
    bool exists = true;

    std::string pFileName = fileName;

    bool isPartitioned;

    getPartitionFileName( pFileName, isPartitioned, comm );

    if ( isPartitioned )
    {
        // partitoned file name, each processor deletes its part, sync results

        exists = FileIO::fileExists( pFileName );

        exists = comm.all( exists );
    }
    else
    {
        int iExists = 0;   // help variable as not bcast<bool> available

        const PartitionId root = 0;

        if ( comm.getRank() == root )
        {
            exists = FileIO::fileExists( pFileName );
            iExists = exists ? 1 : 0 ;
        }

        comm.bcast( &iExists, 1, root );

        exists = iExists == 1;
    }

    return exists;
}

/* --------------------------------------------------------------------------------- */

DistributionPtr PartitionIO::readSDistribution( const string& inFileName, CommunicatorPtr comm )
{
    using hmemo::HArray;

    SCAI_LOG_INFO( logger, "read distribution from one single file " << inFileName )

    HArray<IndexType> owners;
    HArray<IndexType> localSizes( { 0 } );  // at least one entry

    typedef enum
    {
        FAIL,    //!< read was not successul, all processes will throw an exception
        BLOCKED, //!< owners are ascending, a general block distribution is constructed
        GENERAL  //!< owners are arbitrary, a general distribution is construcuted
    } Status;

    Status status = FAIL;

    if ( comm->getRank() == MASTER )
    {
        try
        {
            FileIO::read( owners, inFileName );

            // Bucketsort of owners gives info about illegal values

            HArrayUtils::bucketCount( localSizes, owners, comm->getSize() );

            IndexType nLegalValues = HArrayUtils::sum( localSizes );

            SCAI_LOG_INFO( logger, *comm << ": read owners = " << owners << ", #legal values = " << nLegalValues )

            SCAI_ASSERT_EQ_ERROR( nLegalValues, owners.size(),
                                  *comm << ": mapping file " << inFileName << " contains illegal owners" )
        }
        catch ( common::Exception& e )
        {
            SCAI_LOG_ERROR( logger, "Reading distribution from file " << inFileName << " failed: " << e.what() )
        }

        bool isAscending = HArrayUtils::isSorted( owners, common::CompareOp::LE );

        if ( isAscending )
        {
            status = BLOCKED;
        }
        else
        {
            status = GENERAL;
        }
    }

    // now broadcast status

    comm->bcast( reinterpret_cast<int*>( &status ), 1, MASTER );

    if ( status == FAIL )
    {
        COMMON_THROWEXCEPTION( "Reading distribution failed" )
    }

    IndexType globalSize = owners.size();
    comm->bcast( &globalSize, 1, MASTER );

    DistributionPtr dist;

    if ( status == BLOCKED )
    {
        IndexType localSize;
        hmemo::ReadAccess<IndexType> rSizes( localSizes );
        comm->scatter( &localSize, 1, MASTER, rSizes.get() );
        dist.reset( new GenBlockDistribution ( globalSize, localSize, comm ) );
    }
    else
    {
        // general distribution can be

        dist.reset( new GeneralDistribution( owners, comm ) );
    }

    return dist;
}

/* --------------------------------------------------------------------------------- */

DistributionPtr PartitionIO::readPDistribution( const string& inFileName, CommunicatorPtr comm )
{
    hmemo::HArray<IndexType> owners;

    SCAI_LOG_INFO( logger, *comm << ", read partitioned distribution from " << inFileName )

    hmemo::HArray<IndexType> myIndexes;

    // Some logic needed here to guarantee that all processors will throw an exception
    // if any read fails

    bool errorFlag = false;

    try
    {
        FileIO::read( myIndexes, inFileName );
    }
    catch ( common::Exception& e )
    {
        errorFlag = true;
    }

    errorFlag = comm->any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "Could not read partitioned distribution" )
    }

    IndexType globalSize = comm->sum( myIndexes.size() );

    DistributionPtr dist( new GeneralDistribution( globalSize, myIndexes, comm ) );

    return dist;
}

/* --------------------------------------------------------------------------------- */

DistributionPtr PartitionIO::readDistribution( const string& inFileName, CommunicatorPtr comm )
{
    bool isPartitioned = false;

    string mapFileName = inFileName;

    getPartitionFileName( mapFileName, isPartitioned, *comm );

    if ( isPartitioned )
    {
        return readPDistribution( mapFileName, comm );
    }
    else
    {
        return readSDistribution( mapFileName, comm );
    }
}

/* --------------------------------------------------------------------------------- */

void PartitionIO::write( const Distribution& distribution, const string& fileName )
{
    string distFileName = fileName;

    bool writePartitions;

    getPartitionFileName( distFileName, writePartitions, distribution.getReduceCommunicator() );

    SCAI_LOG_INFO( logger, distribution.getReduceCommunicator() << ": write ( partitioned = " << writePartitions << " ) to " << distFileName )

    if ( writePartitions )
    {
        writePDistribution( distribution, distFileName );
    }
    else
    {
        writeSDistribution( distribution, distFileName );
    }
}

/* --------------------------------------------------------------------------------- */

void PartitionIO::writeSDistribution( const Distribution& distribution, const string& fileName )
{
    using namespace hmemo;

    CommunicatorPtr comm = distribution.getTargetCommunicatorPtr();

    PartitionId rank = comm->getRank();

    SCAI_LOG_INFO( logger, "write distribution to single file " << fileName )

    HArray<IndexType> indexes;

    if ( rank == MASTER )
    {
        // we need the owners only on the host processor
        // indexes = 0, 1, 2, ..., globalSize - 1

        utilskernel::HArrayUtils::setOrder( indexes, distribution.getGlobalSize() );
    }

    HArray<IndexType> owners;

    // Note: gather all owners on master process

    distribution.allOwners( owners, MASTER );

    SCAI_LOG_INFO( logger, *comm << ", owner computation for " << distribution << " finished, owners = " << owners )

    int errorFlag = 0;

    if ( rank == MASTER )
    {
        SCAI_LOG_INFO( logger, *comm << ", MASTER, write distribution to " << fileName )

        try
        {
            FileIO::write( owners, fileName );
        }
        catch ( common::Exception& e )
        {
            SCAI_LOG_ERROR( logger, "master process could not write owner file " << fileName << "\n" << e.what() )
            errorFlag = 1;
        }
    }

    // bcast error status, implies also synchronization

    comm->bcast( &errorFlag, 1, MASTER );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "Error @ writing distrubution to " << fileName )
    }
}

/* --------------------------------------------------------------------------------- */

void PartitionIO::writePDistribution( const Distribution& distribution, const string& fileName )
{
    SCAI_LOG_INFO( logger, distribution.getTargetCommunicator() << ": write distribution to partition file " << fileName )

    // each processor writes a file with its global indexes

    using namespace hmemo;

    HArray<IndexType> myGlobalIndexes;

    distribution.getOwnedIndexes( myGlobalIndexes );

    bool errorFlag = false;

    try
    {
        FileIO::write( myGlobalIndexes, fileName );
    }
    catch ( common::Exception& e )
    {
        errorFlag = true;
    }

    // in case of error: all processors should throw an exception

    const Communicator& comm = distribution.getTargetCommunicator();

    errorFlag = comm.any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "Could not write partitioned distribution" )
    }
}

/* --------------------------------------------------------------------------------- */

}  // namespace

}  // namespace
