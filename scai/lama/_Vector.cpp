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
#include <scai/dmemo/GeneralDistribution.hpp>
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

    Distributed( std::make_shared<NoDistribution>( size ) ),
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

void _Vector::readFromFile( const std::string& fileName )
{
    std::string suffix = FileIO::getSuffix( fileName );

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        std::unique_ptr<FileIO> file( FileIO::create( suffix ) );

        file->open( fileName.c_str(), "r" );
        readFromFile( *file );
        file->close();
    }
    else
    {
        COMMON_THROWEXCEPTION( "File : " << fileName << ", unknown suffix" )
    }
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::readFromFile( const std::string& fileName, DistributionPtr distribution )
{
    // ToDo: distribution might be used for collective mode if it is a block distribution
    readFromFile( fileName );
    redistribute( distribution );
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::readFromFile( const std::string& vectorFileName, const std::string& distributionFileName )
{
    readFromFile( vectorFileName );

    DistributionPtr distribution;

    if ( distributionFileName == "BLOCK" )
    {
        std::string kind = getDistribution().getKind();  // kind of the current distribution

        if ( kind == "BLOCK" ||  kind == "GEN_BLOCK" )
        {
            return;   // already done.
        }

        distribution = dmemo::blockDistribution( this->size() );
    }
    else
    {
        DenseVector<IndexType> owners;

        owners.readFromFile( distributionFileName );

        SCAI_ASSERT_EQ_ERROR( owners.size(), this->size(), 
                              "size of distribution does not match size of vector" )

        distribution = generalDistributionByNewOwners( owners.getDistribution(), owners.getLocalValues() );
    }

    redistribute( distribution );
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

void _Vector::replicate()
{
    if ( getDistribution().isReplicated() )
    {
        return;
    }

    redistribute( getDistribution().toReplicatedDistribution() );
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::writeToFileMaster( FileIO& file ) const
{
    // SINGLE mode: only master process writes to file the whole storage

    CommunicatorPtr comm = file.getCommunicatorPtr();

    const Distribution& dist = getDistribution();

    if ( dist.isMasterDistributed( comm ) )
    {
        SCAI_LOG_DEBUG( logger, *this << ": MASTER writes it" )

        if ( comm->getRank() == 0 )
        {
            writeLocalToFile( file );
        }

        // add a barrier ??, might be better to have it on the close that is called by all processors
    }
    else
    {
        SCAI_LOG_DEBUG( logger, *this << ": master mode, replicate it" )

        DistributionPtr masterDist = dist.toMasterDistribution( comm );

        std::unique_ptr<_Vector> vecMaster( copy() );
        vecMaster->redistribute( masterDist );
        vecMaster->writeToFile( file );
    }
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::writeToFileBlocked( FileIO& file ) const
{
    // this vector must be block distributed onto the processors assigned to the file

    CommunicatorPtr comm = file.getCommunicatorPtr();

    const Distribution& dist = getDistribution();

    SCAI_LOG_DEBUG( logger, "matrix distribution = " << dist << ", file comm = " << *comm )

    if ( dist.isBlockDistributed( comm ) )
    {
        SCAI_LOG_INFO( logger, *this << ": independent/collective write blocked, write local data" )
        writeLocalToFile( file );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, *this << ": independent/collective write, not blocked yet" )

        DistributionPtr blockDist = dist.toBlockDistribution( comm );

        SCAI_ASSERT_DEBUG( blockDist->isBlockDistributed( comm ), *blockDist << " no block distribution" )

        std::unique_ptr<_Vector> vecBlockDistributed( copy() );
        vecBlockDistributed->redistribute( blockDist );
        vecBlockDistributed->writeToFile( file );
    }
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::writeToFile( FileIO& file ) const
{
    // Note: this method must be called by all processors of the file communicator

    SCAI_LOG_INFO( logger, *this << ": writeToFile( file = " << file << " )" )

    if ( file.getDistributedIOMode() == DistributedIOMode::MASTER )
    {
        // master processor will write the whole data

        writeToFileMaster( file );
    }
    else
    {
        // all processors will write the 'block-distributed' data

        writeToFileBlocked( file );
    }
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::writeToFile(
    const std::string& fileName,
    const FileMode fileMode,
    const common::ScalarType dataType,
    const common::ScalarType indexType  ) const
{
    std::string suffix = FileIO::getSuffix( fileName );

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        std::unique_ptr<FileIO> file( FileIO::create( suffix ) );

        file->setMode( fileMode );
        file->setDataType( dataType );
        file->setIndexType( indexType );

        file->open( fileName.c_str(), "w" );
        writeToFile( *file );
        file->close();
    }
    else
    {
        COMMON_THROWEXCEPTION( "File : " << fileName << ", unknown suffix" )
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
