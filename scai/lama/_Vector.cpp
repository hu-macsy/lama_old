/**
 * @file _Vector.cpp
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
 * @brief Implementations of methods for class _Vector.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

// hpp
#include <scai/lama/_Vector.hpp>

// local library

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/Distribution.hpp>

#include <scai/lama/matrix/_Matrix.hpp>
#include <scai/lama/io/PartitionIO.hpp>

// tracing
#include <scai/tracing.hpp>

// std
#include <ostream>

namespace scai
{

using namespace common;
using namespace hmemo;
using namespace dmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( _Vector::logger, "Vector" )

/* ---------------------------------------------------------------------------------------*/
/*    Factory to create a vector                                                          */
/* ---------------------------------------------------------------------------------------*/

_Vector* _Vector::getVector( const VectorKind kind, const common::ScalarType type )
{
    return _Vector::create( VectorCreateKeyType( kind, type ) );
}

/* ---------------------------------------------------------------------------------------*/
/*    Constructor / Destructor                                                            */
/* ---------------------------------------------------------------------------------------*/

_Vector::_Vector( const IndexType size, hmemo::ContextPtr context ) :

    Distributed( DistributionPtr( new NoDistribution( size ) ) ),
    mContext( context )
{
    if ( !mContext )
    {
        mContext = Context::getHostPtr();
    }

    SCAI_LOG_INFO( logger, "Vector(" << size << "), replicated, on " << *mContext )
}

_Vector::_Vector( DistributionPtr distribution, hmemo::ContextPtr context )
    : Distributed( distribution ), mContext( context )
{
    if ( !mContext )
    {
        mContext = Context::getHostPtr();
    }

    SCAI_LOG_INFO( logger,
                   "_Vector(" << distribution->getGlobalSize() << ") with " << getDistribution() << " constructed" )
}

_Vector::_Vector( const _Vector& other )
    : Distributed( other ), mContext( other.getContextPtr() )
{
    SCAI_ASSERT_ERROR( mContext, "NULL context not allowed" )
    SCAI_LOG_INFO( logger, "_Vector(" << other.getDistribution().getGlobalSize() << "), distributed, copied" )
}

_Vector::~_Vector()
{
    SCAI_LOG_DEBUG( logger, "~_Vector(" << getDistribution().getGlobalSize() << ")" )
}

/* ---------------------------------------------------------------------------------------*/
/*    Reading vector from a file, only host reads                                         */
/* ---------------------------------------------------------------------------------------*/

void _Vector::readFromSingleFile( const std::string& fileName )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const PartitionId MASTER = 0;
    const PartitionId myRank = comm->getRank();

    // this is a bit tricky stuff, but it avoids an additional copy from array -> vector

    IndexType globalSize;

    if ( myRank == MASTER )
    {
        SCAI_LOG_INFO( logger, *comm << ": read local array from file " << fileName )

        globalSize = readLocalFromFile( fileName );
    }

    comm->bcast( &globalSize, 1, MASTER );

    if ( myRank != MASTER )
    {
        clearValues();
    }

    DistributionPtr distribution( new SingleDistribution( globalSize, comm, MASTER ) );

    setDistributionPtr( distribution );

    SCAI_LOG_INFO( logger, "readFromSingleFile, vector = " << *this )
}

/* ---------------------------------------------------------------------------------------*/
/*    Reading vector from a file, every processor reads its partition                     */
/* ---------------------------------------------------------------------------------------*/

void _Vector::readFromSingleFile( const std::string& fileName, const DistributionPtr distribution )
{
    if ( distribution.get() == NULL )
    {
        SCAI_LOG_INFO( logger, "readFromSingleFile( " << fileName << ", master only" )
        readFromSingleFile( fileName );
        return;
    }

    const IndexType n = distribution->getBlockDistributionSize();

    if ( n == invalidIndex )
    {
        SCAI_LOG_INFO( logger, "readFromSingleFile( " << fileName << " ), master only + redistribute" )
        readFromSingleFile( fileName );
        redistribute( distribution );
        return;
    }

    // we have a block distribution, so every processor reads its own part

    IndexType first = 0;

    if ( n > 0 )
    {
        first = distribution->local2global( 0 );   // first global index
    }

    bool error = false;

    SCAI_LOG_INFO( logger, "readFromSingleFile( " << fileName << " ), block dist = " << *distribution 
                           << ", read my block, first = " << first << ", n = " << n )

    try
    {
        IndexType localSize = readLocalFromFile( fileName, first, n );
        error = localSize != distribution->getLocalSize();
    }
    catch ( Exception& ex )
    {
        SCAI_LOG_ERROR( logger, ex.what() )
        error = true;
    }

    error = distribution->getCommunicator().any( error );

    if ( error )
    {
        COMMON_THROWEXCEPTION( "readFromSingleFile " << fileName << " failed, dist = " << *distribution )
    }

    setDistributionPtr( distribution );
}

/* ---------------------------------------------------------------------------------*/

void _Vector::readFromPartitionedFile( const std::string& myPartitionFileName, DistributionPtr dist )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    bool errorFlag = false;

    IndexType localSize = 0;

    try
    {
        localSize = readLocalFromFile( myPartitionFileName );

        if ( dist.get() )
        {
            // size of storage must match the local size of distribution

            SCAI_ASSERT_EQUAL( dist->getLocalSize(), localSize, "serious mismatch: local matrix in " << myPartitionFileName << " has illegal local size" )
        }
    }
    catch ( common::Exception& e )
    {
        SCAI_LOG_ERROR( logger, *comm << ": failed to read " << myPartitionFileName << ": " << e.what() )
        errorFlag = true;
    }

    errorFlag = comm->any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "error reading partitioned matrix" )
    }

    DistributionPtr vectorDist = dist;

    if ( !vectorDist.get() )
    {
        // we have no distribution so assume a general block distribution

        IndexType globalSize = comm->sum( localSize );

        vectorDist.reset( new GenBlockDistribution( globalSize, localSize, comm ) );
    }

    setDistributionPtr( vectorDist );   // distribution matches size of local part
}

/* ---------------------------------------------------------------------------------*/

void _Vector::readFromFile( const std::string& vectorFileName, const std::string& distributionFileName )
{
    // read the distribution

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr distribution;

    if ( distributionFileName == "BLOCK" )
    {
        // for a single file we set a BlockDistribution

        if ( vectorFileName.find( "%r" ) == std::string::npos )
        {
            PartitionId root = 0;

            IndexType numRows = invalidIndex;

            if ( comm->getRank() == root )
            {
                numRows = FileIO::getStorageSize( vectorFileName );
            }

            comm->bcast( &numRows, 1, root );

            distribution.reset( new BlockDistribution( numRows, comm ) );
        }

        // for a partitioned file general block distribution is default
    }
    else
    {
        distribution = PartitionIO::readDistribution( distributionFileName, comm );
    }

    readFromFile( vectorFileName, distribution );
}

/* ---------------------------------------------------------------------------------*/

void _Vector::readFromFile( const std::string& fileName, DistributionPtr distribution )
{
    SCAI_LOG_INFO( logger, *this << ": readFromFile( " << fileName << " )" )

    std::string newFileName = fileName;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();  // take default

    if ( distribution.get() )
    {
        comm = distribution->getCommunicatorPtr();
    }

    bool isPartitioned;

    PartitionIO::getPartitionFileName( newFileName, isPartitioned, *comm );

    if ( !isPartitioned )
    {
        readFromSingleFile( newFileName, distribution );
    }
    else
    {
        readFromPartitionedFile( newFileName, distribution );
    }
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::assign( const _HArray& localValues, DistributionPtr dist )
{
    SCAI_ASSERT_EQ_ERROR( localValues.size(), dist->getLocalSize(), "Mismatch local size of vecotr" )

    setDistributionPtr( dist );
    setDenseValues( localValues );
}

void _Vector::assign( const _HArray& globalValues )
{
    SCAI_LOG_INFO( logger, "assign vector with globalValues = " << globalValues )

    setDistributionPtr( DistributionPtr( new NoDistribution( globalValues.size() ) ) );
    setDenseValues( globalValues );
}

void _Vector::assignDistribute( const _HArray& globalValues, DistributionPtr dist )
{
    assign( globalValues );
    redistribute( dist );
}

void _Vector::assignDistribute( const _Vector& other, DistributionPtr distribution )
{
    assign( other );
    redistribute( distribution );
}

/* ---------------------------------------------------------------------------------------*/
/*   writeToFile                                                                          */
/* ---------------------------------------------------------------------------------------*/

void _Vector::writeToSingleFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::ScalarType dataType,
    const FileIO::FileMode fileMode
) const
{
    if ( getDistribution().isReplicated() )
    {
        // make sure that only one processor writes to file

        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        if ( comm->getRank() == 0 )
        {
            writeLocalToFile( fileName, fileType, dataType, fileMode );
        }

        // synchronization to avoid that other processors start with
        // something that might depend on the finally written file

        comm->synchronize();
    }
    else
    {
        // writing a distributed vector into a single file requires redistributon

        std::unique_ptr<_Vector> repV( copy() );
        repV->replicate();
        repV->writeToSingleFile( fileName, fileType, dataType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::replicate()
{
    if ( getDistribution().isReplicated() )
    {
        return;
    }

    redistribute( DistributionPtr( new NoDistribution( size() ) ) );
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::writeToPartitionedFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::ScalarType dataType,
    const FileIO::FileMode fileMode ) const
{
    bool errorFlag = false;

    try
    {
        writeLocalToFile( fileName, fileType, dataType, fileMode );
    }
    catch ( common::Exception& e )
    {
        errorFlag = true;
    }

    const Communicator& comm = getDistribution().getCommunicator();

    errorFlag = comm.any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "Partitioned IO of vector failed" )
    }
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::writeToFile(
    const std::string& fileName,
    const std::string& fileType,               /* = "", take IO type by suffix   */
    const common::ScalarType dataType, /* = UNKNOWN, take defaults of IO type */
    const FileIO::FileMode fileMode            /* = DEFAULT_MODE */
) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": writeToFile( " << fileName << ", fileType = " << fileType << ", dataType = " << dataType << " )" )

    std::string newFileName = fileName;

    bool writePartitions;

    const Communicator& comm = getDistribution().getCommunicator();

    PartitionIO::getPartitionFileName( newFileName, writePartitions, comm );

    if ( !writePartitions )
    {
        writeToSingleFile( newFileName, fileType, dataType, fileMode );
    }
    else
    {
        // matrix_%r.mtx -> matrix_0.4.mtx,  ..., matrix_3.4.mtxt
        writeToPartitionedFile( newFileName, fileType, dataType, fileMode );
    }
}

/* ---------------------------------------------------------------------------------------*/
/*   Miscellaneous                                                                        */
/* ---------------------------------------------------------------------------------------*/

void _Vector::swapVector( _Vector& other )
{
    // swaps only on this base class, not whole vectors
    mContext.swap( other.mContext );
    Distributed::swap( other );
}

void _Vector::writeAt( std::ostream& stream ) const
{
    stream << "Vector(" << getDistributionPtr()->getGlobalSize() << ")";
}

void _Vector::setContextPtr( ContextPtr context )
{
    SCAI_ASSERT_DEBUG( context, "NULL context invalid" )

    if ( mContext->getType() != context->getType() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": new context = " << *context << ", old context = " << *mContext )
    }

    mContext = context;
}

void _Vector::prefetch() const
{
    prefetch( mContext );
}

/* ---------------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
