/**
 * @file SegmentData.cpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Implementation of methods for class SegmentData
 * @author Thomas Brandes
 * @date 06.05.2014
 * @since 1.1.0
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

SCAI_COMMON_INST_CLASS( SegmentData, ARITHMETIC_HOST_CNT, ARITHMETIC_HOST )

} // namespace dmemo

} // namespace scai
