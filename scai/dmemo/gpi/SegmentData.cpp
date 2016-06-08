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

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, SegmentData<ValueType>::logger, "SegmentData" )

/* ---------------------------------------------------------------------------------- */
/*           Constructor                                                              */
/* ---------------------------------------------------------------------------------- */

template<typename T>
SegmentData<T>::SegmentData( const GPICommunicator* comm, const IndexType size )
{
    mComm = comm;  // GPI communicator not really needed

    reserve( size );
}

template<typename T>
void SegmentData<T>::reserve( const IndexType size )
{
    SCAI_REGION( "GASPI.SegmentData.reserve" )

    GPIMemManager::getSegmentData( mId, mPtr, mOffsetBytes, size * sizeof( T ) );

    SCAI_ASSERT_EQUAL_ERROR( 0, mOffsetBytes % sizeof( T ) )

    mData  = static_cast<T*>( mPtr );  // mData will allow for pointer arithmetic helpful for arrays

    mSize = size;

    SCAI_LOG_DEBUG( logger, "SegmentData ( size = " << size << " ), uses id = " << static_cast<int>( mId )
                    << ", offset = " << mOffsetBytes << ", mPtr = " << mPtr )

    mReleaseFlag = true;
}

template<typename T>
SegmentData<T>::SegmentData( const GPICommunicator* comm, const IndexType size, T* data )
{
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

template<typename T>
void SegmentData<T>::assign( const T values[], const IndexType n )
{
    T* data = get();

    if ( values == data )
    {
        SCAI_LOG_INFO( logger, "skipped assign of size " << n )

        // might happen if SegmentData has not been reserved here

        return;
    }

    memcpy( data, values, n * sizeof( T ) );
}

template<typename T>
void SegmentData<T>::copyTo( T values[], const IndexType n ) const
{
    const T* data = get();

    if ( values == data )
    {
        SCAI_LOG_INFO( logger, "skipped copyTo of size " << n )

        // might happen if SegmentData has not been reserved here

        return;
    }

    memcpy( values, data, n * sizeof( T ) );
}

template<typename T>
void SegmentData<T>::release()
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

template<typename T>
void SegmentData<T>::writeAt( std::ostream& stream ) const
{
    stream << "Segment( id = " << static_cast<int>( mId )
           << ", offs = " << mOffsetBytes << ", ptr = " << mPtr << " )";
}

/* --------------------------------------------------------------------------- */

// Note: char is not one of the SCAI supported array types

template class COMMON_DLL_IMPORTEXPORT SegmentData<char> ;

// generic template instantiation for the supported data types

SCAI_COMMON_INST_CLASS( SegmentData, SCAI_ARITHMETIC_HOST_CNT, SCAI_ARITHMETIC_HOST )

} // namespace dmemo

} // namespace scai
