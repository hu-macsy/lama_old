/**
 * @file MICMemory.cpp
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
 * @brief Contains the implementation of the class MICMemory.
 * @author Thomas Brandes
 * @date 15.07.2011
 */

// hpp
#include <scai/hmemo/mic/MICMemory.hpp>

// local library
#include <scai/hmemo/mic/MICContext.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/common/macros/assert.hpp>

// std
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

void MICMemory::memset( void* dst, const int val, const size_t size ) const
{
    const size_t dst_ptr = (size_t) dst;

    int deviceNr = mMICContext->getDeviceNr();

#pragma offload target( mic : deviceNr ) in( dst_ptr, size, val )
    {
        void* dst = ( void* ) dst_ptr;

        ::memset( dst, val, size );
    }
}

/* ----------------------------------------------------------------------------- */

bool MICMemory::canCopyFrom( const Memory& other ) const
{
    bool supported = false;

    memtype::MemoryType otherType = other.getType();

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

    memtype::MemoryType otherType = other.getType();

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
