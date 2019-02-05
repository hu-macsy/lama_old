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

void _Vector::readFromFile( const std::string& fileName, const std::string& fileType )
{
    std::string suffix = fileType;

    if ( suffix == "" )
    {
        suffix = FileIO::getSuffix( fileName );
    }

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

    redistribute( dmemo::noDistribution( size() ) );
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::writeToFile( FileIO& file ) const
{
    SCAI_LOG_INFO( logger, *this << ": writeToFile( file = " << file << " )" )

    if ( file.getDistributedIOMode() == DistributedIOMode::SINGLE )
    {
        // SINGLE mode: only master process writes to file

        if ( getDistribution().isReplicated() )
        {
            const Communicator& comm = *Communicator::getCommunicatorPtr();

            if ( comm.getRank() == 0 )
            {
                writeLocalToFile( file );
            }
        }
        else
        {
            std::unique_ptr<_Vector> repV( copy() );
            repV->replicate();
            repV->writeToFile( file );
        }
    }
    else
    {
        // INDEPENDENT / COLLECTIVE mode: write local parts of block distributed vector

        if ( getDistribution().getBlockDistributionSize() == invalidIndex )
        {
            std::unique_ptr<_Vector> vBlockDistributed( copy() );
            vBlockDistributed->redistribute( blockDistribution( vBlockDistributed->size() ) );
            vBlockDistributed->writeToFile( file );
        }
        else
        {
            writeLocalToFile( file );
        }
    }
}

/* ---------------------------------------------------------------------------------------*/

void _Vector::writeToFile(
    const std::string& fileName,
    const std::string& fileType,               /* = "", take IO type by suffix   */
    const common::ScalarType dataType, /* = UNKNOWN, take defaults of IO type */
    const FileMode fileMode            /* = DEFAULT_MODE */
) const
{
    std::string suffix = fileType;

    if ( suffix == "" )
    {
        suffix = FileIO::getSuffix( fileName );
    }

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        std::unique_ptr<FileIO> file( FileIO::create( suffix ) );

        if ( dataType != common::ScalarType::UNKNOWN )
        {
            // overwrite the default settings

            file->setDataType( dataType );
        }

        if ( fileMode != FileMode::DEFAULT )
        {
            // overwrite the default settings

            file->setMode( fileMode );
        }

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
