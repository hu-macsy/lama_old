/**
 * @file _MatrixStorage.cpp
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
 * @brief Implementation of methods for common base class of all matrix storage formats.
 * @author Thomas Brandes
 * @date 27.04.2011
 */

// hpp
#include <scai/lama/storage/_MatrixStorage.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/hmemo/Context.hpp>

#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>

#include <scai/tracing.hpp>

namespace scai
{

using namespace hmemo;
using namespace dmemo;

using common::BinaryOp;

using utilskernel::UtilKernelTrait;
using utilskernel::OpenMPUtils;

namespace lama
{

SCAI_LOG_DEF_LOGGER( _MatrixStorage::logger, "MatrixStorage" )

/* ---------------------------------------------------------------------------------- */
/*   Constructors                                                                     */
/* ---------------------------------------------------------------------------------- */

_MatrixStorage::_MatrixStorage( IndexType numRows, IndexType numColumns, ContextPtr ctx ) :

    mNumRows( numRows ),
    mNumColumns( numColumns ),
    mRowIndexes(), 
    mCompressThreshold( 0.0f ), 
    mContext( ctx )

{
    SCAI_ASSERT_ERROR( ctx, "Null context not allowed for storage" )
    SCAI_LOG_DEBUG( logger, "constructed MatrixStorage( " << mNumRows << " x " << mNumColumns << " ) @ " << *mContext )
}

_MatrixStorage::_MatrixStorage( const _MatrixStorage& other ) :

    mNumRows( other.mNumRows ),
    mNumColumns( other.mNumColumns ),
    mRowIndexes(),
    mCompressThreshold( other.mCompressThreshold ),
    mContext( other.mContext )

{
    SCAI_LOG_DEBUG( logger, "_MatrixStorage( " << mNumRows << " x " << mNumColumns << " copied" )
}

_MatrixStorage::_MatrixStorage( _MatrixStorage&& other ) noexcept :

    mNumRows( other.mNumRows ),
    mNumColumns( other.mNumColumns ),
    mRowIndexes( std::move( other.mRowIndexes ) ),
    mCompressThreshold( other.mCompressThreshold ),
    mContext( other.mContext )

{
    // reset the shape of the moved input storage

    other.mNumRows = 0;
    other.mNumColumns = 0;

    SCAI_LOG_DEBUG( logger, "_MatrixStorage( " << mNumRows << " x " << mNumColumns << " moved" )
}

_MatrixStorage::~_MatrixStorage()
{
    SCAI_LOG_DEBUG( logger, "~_MatrixStorage" )
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::setDimension( const IndexType numRows, const IndexType numColumns )
{
    // in any case set dimensions
    mNumRows = numRows;
    mNumColumns = numColumns;
    mRowIndexes.clear();
    // but do not reset threshold
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::splitUp( IndexType& numRows, IndexType& numColumns )
{
    numRows = mNumRows;
    numColumns = mNumColumns;

    mNumRows = 0;
    mNumColumns = 0;
    mRowIndexes.clear();
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::resetNumColumns( const IndexType numColumns )
{
    SCAI_ASSERT_GE_ERROR( numColumns, mNumColumns, "Number of columns cannot be decreased" )

    if ( numColumns != mNumColumns )
    {
        SCAI_ASSERT_NE_ERROR( getFormat(), Format::DENSE, "number of columns cannot be changed for dense storage." )
    }

    mNumColumns = numColumns;
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::writeAt( std::ostream& stream ) const
{
    stream << " MatrixStorage: (" << mNumRows << " x " << mNumColumns << ")";
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::setCompressThreshold( float ratio )
{
    if ( ratio < 0.0f || ratio > 1.0f )
    {
        COMMON_THROWEXCEPTION( "Illegal threshold " << ratio << ", must be from 0.0 to 1.0" )
    }

    mCompressThreshold = ratio;
    SCAI_LOG_INFO( logger, "set compress threshold, ratio = " << ratio << " : " << *this )
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::swap( _MatrixStorage& other )
{
    std::swap( mNumRows, other.mNumRows );
    std::swap( mNumColumns, other.mNumColumns );
    mRowIndexes.swap( other.mRowIndexes );
    std::swap( mCompressThreshold, other.mCompressThreshold );
    std::swap( mContext, other.mContext );
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::_assignTranspose( const _MatrixStorage& other )
{
    // make it safe also for other == &this
    IndexType tmpNumRows = other.mNumRows;
    mNumRows = other.mNumColumns;
    mNumColumns = tmpNumRows;
    mRowIndexes.clear();
    // remains unchanged: mCompressThreshold = other.mCompressThreshold;
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::_assign( const _MatrixStorage& other )
{
    mNumRows = other.mNumRows;
    mNumColumns = other.mNumColumns;
    mRowIndexes.clear();
    // remains unchanged: mCompressThreshold = other.mCompressThreshold;
}

/* --------------------------------------------------------------------------- */

_MatrixStorage& _MatrixStorage::operator=( const _MatrixStorage& other )
{
    assign( other ); // assign can deal with all kind of storage formats/types
    return *this;
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::moveImpl( _MatrixStorage&& other )
{
    SCAI_LOG_DEBUG( logger, "move assignment: other = " << other << " is moved to this = " << *this )

    mNumRows    = other.mNumRows;
    mNumColumns = other.mNumColumns;

    // reset the shape of the moved input storage

    other.mNumRows    = 0;
    other.mNumColumns = 0;

    mContext = other.mContext;

    mRowIndexes = std::move( other.mRowIndexes );
    mCompressThreshold = other.mCompressThreshold;
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::setContextPtr( ContextPtr context )
{
    if ( context.get() != mContext.get() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": new location = " << *context << ", old location = " << *mContext )
    }

    mContext = context;
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::localize( const _MatrixStorage& global, const Distribution& rowDist )
{
    SCAI_ASSERT_EQUAL_ERROR( getNumColumns(), global.getNumColumns() )
    SCAI_ASSERT_EQUAL_ERROR( global.getNumRows(), rowDist.getGlobalSize() )
    COMMON_THROWEXCEPTION( "No default implementation for localize available, matrix = " << *this )
}

/* --------------------------------------------------------------------------- */

IndexType _MatrixStorage::getNumValues() const
{
    // Default implementation builds sum of row sizes, derived classes have more efficient routines

    HArray<IndexType> sizes;

    buildCSRSizes( sizes );

    IndexType numValues = utilskernel::HArrayUtils::reduce( sizes, BinaryOp::ADD );

    return numValues;
}

/* ========================================================================= */
/*       File I/O                                                            */
/* ========================================================================= */

void _MatrixStorage::writeToFile(
    const std::string& fileName,
    const FileMode mode,
    const common::ScalarType dataType,
    const common::ScalarType indexType ) const
{
    SCAI_REGION( "Storage.writeToFile" )

    std::string suffix = FileIO::getSuffix( fileName );

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        std::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );

        fileIO->setDataType( dataType );
        fileIO->setIndexType( indexType );
        fileIO->setMode( mode );

        SCAI_LOG_INFO( logger, "write matrix storage to file, FileIO = " << *fileIO << ", storage = " << *this )

        fileIO->open( fileName.c_str(), "w", DistributedIOMode::INDEPENDENT );
        fileIO->writeStorage( *this );
        fileIO->close();
    }
    else
    {
        COMMON_THROWEXCEPTION( "writeToFile " << fileName << ", unknown file type " << suffix )
    }
}

/* ---------------------------------------------------------------------------------- */

size_t _MatrixStorage::_getMemoryUsage() const
{
    size_t memoryUsage = 0;

    memoryUsage += 2 * sizeof( IndexType );                  // mNumRows, mNumColumns
    memoryUsage += sizeof( float );                          // mCompressThreshold
    memoryUsage += sizeof( IndexType ) * mRowIndexes.size(); // mRowIndexes

    return memoryUsage;
}

} /* end namespace lama */

} /* end namespace scai */
