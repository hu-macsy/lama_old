/**
 * @file MICMemory.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Contains the implementation of the class MICMemory.
 * @author Thomas Brandes
 * @date 15.07.2011
 * @since 1.1.0
 */

// hpp

#include <scai/lama/mic/MICMemory.hpp>
#include <scai/common/exception/Exception.hpp>
#include <scai/tracing.hpp>

#include <scai/lama/mic/MICContext.hpp>
#include <scai/tasking/NoSyncToken.hpp>

#include <cstring>

namespace scai
{

namespace hmemo
{

/**  static variables *****************************************************/

SCAI_LOG_DEF_LOGGER( MICMemory::logger, "MIC.Memory" )

/**  constructor  *********************************************************/

MICMemory::MICMemory( common::shared_ptr<const MICContext> micContext )

    : Memory( memtype::MICMemory ),
      mMICContext( micContext )

{
    SCAI_ASSERT( micContext, "NULL context for MIC memory" )

    SCAI_LOG_DEBUG( logger, "construct MICMemory for context " << micContext )
}

/**  destructor   *********************************************************/

MICMemory::~MICMemory()
{
    SCAI_LOG_INFO( logger, "~MICMemory: " << *this )
}

/**  const methods ********************************************************/

int MICMemory::getDeviceNr() const
{
    return mMICContext->getDeviceNr();
}

/* ----------------------------------------------------------------------------- */

void MICMemory::writeAt( std::ostream& stream ) const
{
    stream << "MICMemory( @ " << *mMICContext << " )";
}

/* ----------------------------------------------------------------------------- */

void* MICMemory::allocate( const size_t size ) const
{
    void* pointer = NULL;

    SCAI_LOG_INFO( logger, "allocate, init pointer = " << pointer )

    int deviceNr = mMICContext->getDeviceNr();

    // allocate data on the MIC device and return the pointer

    #pragma offload target( mic : deviceNr ), in( size ), out( pointer )
    {
        pointer = ::malloc( size );
    }

    SCAI_LOG_INFO( logger, "allocated " << size << " bytes on device, ptr = " << pointer )

    return pointer;
}

/* ----------------------------------------------------------------------------- */

void MICMemory::free( void* pointer, const size_t size ) const
{
    SCAI_LOG_INFO( logger, *this << ": free " << size << " bytes, ptr = " << pointer )

    int deviceNr = mMICContext->getDeviceNr();

    #pragma offload target( mic : deviceNr ), in( pointer )
    {
        ::free( pointer );
    }
}

/* ----------------------------------------------------------------------------- */

void MICMemory::memcpy( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "memcpy " << size << " bytes on MIC device, " << dst << " <- " << src )

    const size_t dst_ptr = (size_t) dst;
    const size_t src_ptr = (size_t) src;

    int deviceNr = mMICContext->getDeviceNr();

    if( ( size & 7 ) == 0 && ( dst_ptr & 7 ) == 0 && ( src_ptr & 7 ) == 0 )
    {
        // SCAI_REGION( "MIC.memcpy8" )

        #pragma offload target( mic : deviceNr ) in( src_ptr, dst_ptr, size )
        {
            double* dst = ( double* ) dst_ptr;
            const double* src = ( const double* ) src_ptr;

            #pragma omp parallel for
            for ( int i = 0; i < ( size >> 3 ); ++ i )
            {
               dst[i] = src[i];
            }
        }
    }
    else if( ( size & 3 ) == 0 && ( dst_ptr & 3 ) == 0 && ( src_ptr & 3 ) == 0 )
    {
        // SCAI_REGION( "MIC.memcpy4" )

        #pragma offload target( mic : deviceNr ) in( src_ptr, dst_ptr, size )
        {
            float* dst = ( float* ) dst_ptr;
            const float* src = ( const float* ) src_ptr;

            #pragma omp parallel for
            for ( int i = 0; i < ( size >> 2 ); ++ i )
            {
               dst[i] = src[i];
            }
        }
    }
    else
    {
        // SCAI_REGION( "MIC.memcpy1" )

        #pragma offload target( mic : deviceNr ) in( src_ptr, dst_ptr, size )
        {
            void* dst = ( void* ) dst_ptr;
            const void* src = ( const void* ) src_ptr;

            ::memcpy( dst, src, size );
        }
    }
}

/* ----------------------------------------------------------------------------- */

bool MICMemory::canCopyFrom( const Memory& other ) const
{
    bool supported = false;

    MemoryType otherType = other.getType();

    if ( otherType == memtype::HostMemory )
    {
        // MICDevice -> Host is supported

        supported = true;
    }

    SCAI_LOG_INFO( logger, "canCopyFrom " << other << " to this " << *this << ", supported = " << supported )

    return supported;
}

/* ----------------------------------------------------------------------------- */

bool MICMemory::canCopyTo( const Memory& other ) const
{
    bool supported = false;

    MemoryType otherType = other.getType();

    if ( otherType == memtype::HostMemory )
    {
        // MICMemory -> HostMemory is supported

        supported = true;
    }

    SCAI_LOG_INFO( logger, "canCopyTo " << other << " from this " << *this << ", supported = " << supported )

    return supported;
}

/* ----------------------------------------------------------------------------- */

void MICMemory::memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const
{
    if ( srcMemory.getType() == memtype::HostMemory )
    {
        memcpyFromHost( dst, src, size );
    }
    else
    {
        SCAI_LOG_ERROR( logger, "copy from " << srcMemory << " to " << *this << " not supported" )
        COMMON_THROWEXCEPTION( "copy from " << srcMemory << " to " << *this << " not supported" )
    }
}

/* ----------------------------------------------------------------------------- */

void MICMemory::memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const
{
    if ( dstMemory.getType() == memtype::HostMemory )
    {
        memcpyToHost( dst, src, size );
    }
    else
    {
        SCAI_LOG_ERROR( logger, "copy to " << dstMemory << " from " << *this << " not supported" )
        COMMON_THROWEXCEPTION( "copy to " << dstMemory << " from " << *this << " not supported" )
    }
}

/* ----------------------------------------------------------------------------- */

void MICMemory::memcpyToHost( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "memcpy " << size << " bytes from MIC to Host" )

    // SCAI_REGION( "MIC.memcpyToHost" )

    uint8_t* dst8 = ( uint8_t* ) dst;

    int deviceNr = mMICContext->getDeviceNr();

    #pragma offload target( mic : deviceNr ) out( dst8 : length( size ) ), in( size ), in( src )
    {
        const uint8_t* src8 = ( const uint8_t* ) src;
        ::memcpy( dst8, src8, size );
    }
}

/* ----------------------------------------------------------------------------- */

void MICMemory::memcpyFromHost( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "memcpy " << size << " bytes on MIC device, dst = " << dst << ", src = " << src )

    // SCAI_REGION( "MIC.memcpyFromHost" )

    uint8_t* src8 = ( uint8_t* ) src;

    int deviceNr = mMICContext->getDeviceNr();

    #pragma offload target( mic : deviceNr ) in( src8 : length( size ) ), in( size ), in( dst )
    {
        uint8_t* dst8 = ( uint8_t* ) dst;
        ::memcpy( dst8, src8, size );
    }
}

/* ----------------------------------------------------------------------------- */

tasking::SyncToken* MICMemory::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    memcpy( dst, src, size );
    return new tasking::NoSyncToken();
}

/* ----------------------------------------------------------------------------- */

ContextPtr MICMemory::getContextPtr() const
{
    return mMICContext;
}

} /* end namespace hmemo */

} /* end namespace scai */
