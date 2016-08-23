/**
 * @file SegmentData.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Implementation of methods for class SegmentData
 * @author Thomas Brandes
 * @date 06.05.2014
 */

// hpp
#include <scai/dmemo/gpi/SegmentData.hpp>

// GASPI logging
#include <scai/dmemo/gpi/GPIUtils.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/assert.hpp>

using namespace std;

namespace scai
{

namespace dmemo
{

/* ---------------------------------------------------------------------------------- */
/*           Logger                                                                   */
/* ---------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( SegmentData::logger, "SegmentData" )

/* ---------------------------------------------------------------------------------- */
/*           Constructor                                                              */
/* ---------------------------------------------------------------------------------- */

SegmentData::SegmentData( common::scalar::ScalarType stype, const GPICommunicator* comm, const IndexType size )
{
    mScalarType = stype;
    mTypeSize = common::typeSize( stype );
    mComm = comm;  // GPI communicator not really needed
    reserve( size );
}

void SegmentData::reserve( const IndexType size )
{
    SCAI_REGION( "GASPI.SegmentData.reserve" )
    GPIMemManager::getSegmentData( mId, mPtr, mOffsetBytes, size * mTypeSize );
    SCAI_ASSERT_EQ_ERROR( 0, mOffsetBytes % mTypeSize, "unaligned data" )
    mData  = mPtr;  // mData will allow for pointer arithmetic helpful for arrays
    mSize = size;
    SCAI_LOG_DEBUG( logger, "SegmentData ( size = " << size << " ), uses id = " << static_cast<int>( mId )
                    << ", offset = " << mOffsetBytes << ", mPtr = " << mPtr )
    mReleaseFlag = true;
}

SegmentData::SegmentData( common::scalar::ScalarType stype, const GPICommunicator* comm, const IndexType size, void* data )
{
    mScalarType = stype;
    mTypeSize = common::typeSize( stype );

    mComm = comm;  // GPI communicator not really needed
    bool found = GPIMemManager::findSegment( mId, mOffsetBytes, data );
    SCAI_LOG_DEBUG( logger, "SegmentData( size = " << size << ", data = " << data
                    << ", found = " << found )

    if ( found )
    {
        SCAI_LOG_DEBUG( logger, "no reservation of SegmentData required, size = " << size )
        // @ToDo: verify that segment is sufficient
        mPtr  = data;
        mData = data;
        mSize = size;
        mReleaseFlag = false; // do not release at end
    }
    else
    {
        reserve( size );
    }
}

void SegmentData::assign( const void* values, const IndexType n )
{
    void* data = get();

    if ( values == data )
    {
        SCAI_LOG_INFO( logger, "skipped assign of size " << n )
        // might happen if SegmentData has not been reserved here
        return;
    }

    memcpy( data, values, n * mTypeSize );
}

void SegmentData::copyTo( void* values, const IndexType n ) const
{
    const void* data = get();

    if ( values == data )
    {
        SCAI_LOG_INFO( logger, "skipped copyTo of size " << n )
        // might happen if SegmentData has not been reserved here
        return;
    }

    memcpy( values, data, n * mTypeSize );
}

void SegmentData::release()
{
    if ( mComm )
    {
        if ( mReleaseFlag )
        {
            GPIMemManager::releaseSegmentData( mId, mOffsetBytes );
        }

        mComm = NULL;
        mPtr  = NULL;
    }
}

void SegmentData::writeAt( std::ostream& stream ) const
{
    stream << "Segment<" << mScalarType<< ">( id = " << static_cast<int>( mId )
           << ", offs = " << mOffsetBytes << ", ptr = " << mPtr << " )";
}

/* --------------------------------------------------------------------------- */

} // namespace dmemo

} // namespace scai
