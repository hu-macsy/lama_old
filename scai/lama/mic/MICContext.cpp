/**
 * @file MICContext.cpp
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
 * @brief Contains the implementation of the class MICContext.
 * @author Thomas Brandes
 * @date 15.07.2011
 * @since 1.1.0
 */

// hpp
#include <scai/lama/mic/MICContext.hpp>

#include <scai/lama/exception/LAMAAssert.hpp>
#include <scai/lama/NoSyncToken.hpp>
#include <scai/lama/mic/MICSyncToken.hpp>

#include <scai/tracing.hpp>
#include <omp.h>
#include <mkl.h>

namespace lama
{

/**  static variables *****************************************************/

SCAI_LOG_DEF_LOGGER( MICContext::logger, "MIC.Context" )

int MICContext::currentDeviceNr = -1;

int MICContext::numUsedDevices = 0;

/**  constructor  *********************************************************/

__attribute__( ( target( mic ) ) ) bool chk_target00()
{
    bool retval;

#ifdef __MIC__
    retval = true;
#else
    retval = false;
#endif

    return retval;
}

MICContext::MICContext( int deviceNr )
    : Context( MIC ), mDeviceNr( deviceNr )
{
    SCAI_LOG_INFO( logger, "construct MICContext, device nr = = " << deviceNr )

    int numDevices = 0;

    SCAI_LOG_INFO( logger, "Checking for Intel(R) MIC Architecture (Target CPU) devices" )

#ifdef __INTEL_OFFLOAD
    numDevices = _Offload_number_of_devices();
#endif

    if( numDevices < 1 )
    {
        COMMON_THROWEXCEPTION( "No mic devices available" )
    }

    // ToDo: allow for any device

    mDeviceNr = 0;

    SCAI_LOG_INFO( logger, "Using device " << mDeviceNr << " of " << numDevices << " installed devices" )

    bool targetOK = false;

#pragma offload target( mic: mDeviceNr ) in( mDeviceNr )
    targetOK = chk_target00();

    if( !targetOK )
    {
        COMMON_THROWEXCEPTION( "Could not offload to device " << mDeviceNr )
    }

    int numCores;

#pragma offload target( mic: mDeviceNr ) out( numCores )
    {
        #pragma omp parallel
        {
            #pragma omp single
            {
                numCores = omp_get_num_threads();
            }
        }
    }

    SCAI_LOG_INFO( logger, "Uses " << numCores << " threads for parallel execution" );

    mNumThreads = numCores;
}

/**  destructor   *********************************************************/

MICContext::~MICContext()
{
    SCAI_LOG_INFO( logger, "~MICContext: " << *this )
}

/* ----------------------------------------------------------------------------- */

void MICContext::writeAt( std::ostream& stream ) const
{
    stream << "MICContext(" << mDeviceNr << ": " << mNumThreads << " Threads)";
}

/* ----------------------------------------------------------------------------- */

void MICContext::disable( const char* file, int line ) const
{
    currentDeviceNr = -1;
}

/* ----------------------------------------------------------------------------- */

void MICContext::enable( const char* file, int line ) const
{
    SCAI_LOG_DEBUG( logger, *this << ": enable" )
    currentDeviceNr = mDeviceNr;
}

/* ----------------------------------------------------------------------------- */

bool MICContext::canUseData( const Context& other ) const
{
    // same object by pointer can always use same data.

    if( this == &other )
    {
        return true;
    }

    // MIC device can use only data on same MIC device

    if( other.getType() == MIC )
    {
        const MICContext& otherMIC = static_cast<const MICContext&>( other );
        return otherMIC.mDeviceNr == mDeviceNr;
    }

    return false;
}

/* ----------------------------------------------------------------------------- */

void* MICContext::allocate( const size_t size ) const
{
    void* pointer = NULL;

    SCAI_LOG_INFO( logger, "allocate, init pointer = " << pointer )

#pragma offload target( mic : mDeviceNr ), in( size ), out( pointer )
    {
        pointer = ::malloc( size );
    }

    SCAI_LOG_INFO( logger, "allocated " << size << " bytes on device, ptr = " << pointer )

    return pointer;
}

/* ----------------------------------------------------------------------------- */

void MICContext::allocate( ContextData& contextData, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "context allocate " << size << " bytes on MIC device" )
    contextData.pointer = allocate( size );
}

/* ----------------------------------------------------------------------------- */

void MICContext::free( void* pointer, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "free " << size << " bytes on MIC device, ptr = " << pointer )

#pragma offload target( mic : mDeviceNr ), in( pointer )
    {
        ::free( pointer );
    }
}

/* ----------------------------------------------------------------------------- */

void MICContext::free( ContextData& contextData ) const
{
    LAMA_ASSERT_EQUAL_ERROR( contextData.context->getType(), getType() )
    free( contextData.pointer, contextData.size );
}

/* ----------------------------------------------------------------------------- */

void MICContext::memcpy( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "memcpy " << size << " bytes on MIC device, " << dst << " <- " << src )

    const size_t dst_ptr = (size_t) dst;
    const size_t src_ptr = (size_t) src;

    if( ( size & 7 ) == 0 && ( dst_ptr & 7 ) == 0 && ( src_ptr & 7 ) == 0 )
    {
        LAMA_REGION( "MIC.memcpy8" )

#pragma offload target( mic : mDeviceNr ) in( src_ptr, dst_ptr, size )
        {
            double* dst = (double*) dst_ptr;
            const double* src = (const double*) src_ptr;

            #pragma omp parallel for

            for( int i = 0; i < ( size >> 3 ); ++i )
            {
                dst[i] = src[i];
            }
        }
    }
    else if( ( size & 3 ) == 0 && ( dst_ptr & 3 ) == 0 && ( src_ptr & 3 ) == 0 )
    {
        LAMA_REGION( "MIC.memcpy4" )

#pragma offload target( mic : mDeviceNr ) in( src_ptr, dst_ptr, size )
        {
            float* dst = (float*) dst_ptr;
            const float* src = (const float*) src_ptr;

            #pragma omp parallel for

            for( int i = 0; i < ( size >> 2 ); ++i )
            {
                dst[i] = src[i];
            }
        }
    }
    else
    {
        LAMA_REGION( "MIC.memcpy1" )

#pragma offload target( mic : mDeviceNr ) in( src_ptr, dst_ptr, size )
        {
            void* dst = (void*) dst_ptr;
            const void* src = (const void*) src_ptr;

            ::memcpy( dst, src, size );
        }
    }
}

/* ----------------------------------------------------------------------------- */

void MICContext::memcpyToHost( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "memcpy " << size << " bytes from MIC to Host" )

    LAMA_REGION( "MIC.memcpyToHost" )

    uint8_t* dst8 = (uint8_t*) dst;

#pragma offload target( mic : mDeviceNr ) out( dst8 : length( size ) ), in( size ), in( src )
    {
        const uint8_t* src8 = (const uint8_t*) src;
        ::memcpy( dst8, src8, size );
    }
}

/* ----------------------------------------------------------------------------- */

void MICContext::memcpyFromHost( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "memcpy " << size << " bytes on MIC device, dst = " << dst << ", src = " << src )

    LAMA_REGION( "MIC.memcpyFromHost1" )

    uint8_t* src8 = (uint8_t*) src;

#pragma offload target( mic : mDeviceNr ) in( src8 : length( size ) ), in( size ), in( dst )
    {
        uint8_t* dst8 = (uint8_t*) dst;
        ::memcpy( dst8, src8, size );
    }
}

/* ----------------------------------------------------------------------------- */

SyncToken* MICContext::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    memcpy( dst, src, size );
    return new NoSyncToken();
}

/* ----------------------------------------------------------------------------- */

bool MICContext::cancpy( const ContextData& dst, const ContextData& src ) const
{
    const Context::ContextType dstType = dst.context->getType();
    const Context::ContextType srcType = src.context->getType();

    return ( srcType == Host && dstType == MIC ) || ( srcType == MIC && dstType == Host )
           || ( srcType == MIC && dstType == MIC );
}

/* ----------------------------------------------------------------------------- */

void MICContext::memcpy( ContextData& dst, const ContextData& src, const size_t size ) const
{
    const Context::ContextType dstType = dst.context->getType();
    const Context::ContextType srcType = src.context->getType();

    if( srcType == Host && dstType == MIC )
    {
        SCAI_LOG_INFO( logger, "copy from Host to MIC" )
        memcpyFromHost( dst.pointer, src.pointer, size );
        SCAI_LOG_INFO( logger, "ready copy from Host to MIC" )
    }
    else if( srcType == MIC && dstType == Host )
    {
        SCAI_LOG_INFO( logger, "copy from MIC to Host" )
        memcpyToHost( dst.pointer, src.pointer, size );
        SCAI_LOG_INFO( logger, "ready copy from MIC to Host" )
    }
    else if( srcType == MIC && dstType == MIC )
    {
        SCAI_LOG_INFO( logger, "copy from MIC to MIC" )
        memcpy( dst.pointer, src.pointer, size );
        SCAI_LOG_INFO( logger, "ready copy from MIC to MIC" )
    }
    else
    {
        COMMON_THROWEXCEPTION( "Can not copy from "<< *src.context << " to " << *dst.context )
    }
}

/* ----------------------------------------------------------------------------- */

SyncToken* MICContext::memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const
{
    memcpy( dst, src, size );

    return new NoSyncToken();
}

/* ----------------------------------------------------------------------------- */

SyncToken* MICContext::getSyncToken() const
{
    return new MICSyncToken( shared_from_this() );
}

/* ----------------------------------------------------------------------------- */

int MICContext::getCurrentDevice()
{
    // ToDo: get current device, make sure that access has been enabled

    return 0;
}

/* ----------------------------------------------------------------------------- */

} //namespace lama

