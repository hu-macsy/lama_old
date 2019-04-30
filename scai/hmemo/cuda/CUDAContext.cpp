/**
 * @file CUDAContext.cpp
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
 * @brief Contains the implementation of the class CUDAContext.
 * @author Thomas Brandes, Jiri Kraus
 * @date 15.07.2011
 */

// local library
#include <scai/hmemo/cuda/CUDAContext.hpp>
#include <scai/hmemo/cuda/CUDAHostMemory.hpp>
#include <scai/hmemo/cuda/CUDAMemory.hpp>

#include <scai/hmemo/ContextAccess.hpp>

// internal scai libraries

#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>
#include <scai/tasking/cuda/CUDAStreamPool.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/macros/assert.hpp>

// std
#include <memory>

namespace scai
{

using tasking::SyncToken;
using tasking::StreamType;
using tasking::CUDAStreamSyncToken;
using tasking::CUDAStreamPool;

namespace hmemo
{

/* static variables *****************************************************/

SCAI_LOG_DEF_LOGGER( CUDAContext::logger, "Context.CUDAContext" )

/* constructor  *********************************************************/

CUDAContext::CUDAContext( int deviceNr ) :

    Context( common::ContextType::CUDA ),
    CUDACtx( deviceNr )

{
    SCAI_LOG_INFO( logger, "construct CUDAContext, device nr = = " << deviceNr )
    // Note: logging is safe as CUDA context is always created after static initializations
    {
        char deviceName[256];
        SCAI_CUDA_DRV_CALL( cuDeviceGetName( deviceName, 256, getCUdevice() ), "cuDeviceGetName" );
        mDeviceName = deviceName; // save it as string member variable for output
        SCAI_LOG_DEBUG( logger, "got device " << mDeviceName )
    }
    // count used devices to shutdown CUDA if no more device is used
    SCAI_LOG_DEBUG( logger, *this << " constructed, is disabled" )
}

/*  destructor  ---------------------------------------------------------------- */

CUDAContext::~CUDAContext()
{
    SCAI_LOG_INFO( logger, "~CUDAContext: " << *this )
}

/* ----------------------------------------------------------------------------- */

bool CUDAContext::isEqual( const Context& other ) const
{
    if ( &other == this )
    {
        return true;
    }

    if ( other.getType() != common::ContextType::CUDA )
    {
        return false;
    }

    const CUDAContext& otherCUDA = static_cast<const CUDAContext&>( other );
   
    return getDeviceNr() == otherCUDA.getDeviceNr();
}

/* ----------------------------------------------------------------------------- */

MemoryPtr CUDAContext::getLocalMemoryPtr() const
{
    MemoryPtr memory;

    if ( mMemory.expired() )
    {
        memory.reset( new CUDAMemory( shared_from_this() ) );
        mMemory = memory; // save it here as a weak pointer to avoid cycles
    }
    else
    {
        // the last memory instance is still valid, so we return just shared pointer to it
        memory = mMemory.lock();
    }

    return memory;
}

/* ----------------------------------------------------------------------------- */

MemoryPtr CUDAContext::getHostMemoryPtr() const
{
    MemoryPtr memory;

    if ( mHostMemory.expired() )
    {
        memory.reset( new CUDAHostMemory( shared_from_this() ) );
        mHostMemory = memory; // save it here as a weak pointer to avoid cycles
    }
    else
    {
        // the last host context instance is still valid, so we return just shared pointer to it
        memory = mHostMemory.lock();
    }

    return memory;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::writeAt( std::ostream& stream ) const
{
    stream << "CUDAContext(" << getDeviceNr() << ": " << mDeviceName << ")";
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::disable( const char* file, int line ) const
{
    Context::disable( file, line ); // call routine of base class
    common::CUDAAccess::disable( this );
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::enable( const char* file, int line ) const
{
    Context::enable( file, line ); // call routine of base class
    common::CUDAAccess::enable( *this );
}

/* ----------------------------------------------------------------------------- */

bool CUDAContext::canUseMemory( const Memory& other ) const
{
    bool canUse = false;

    // CUDA device can use only data on same CUDA device

    if ( other.getType() == MemoryType::CUDAMemory )
    {
        const CUDAMemory* otherCUDAMem = dynamic_cast<const CUDAMemory*>( &other );
        SCAI_ASSERT( otherCUDAMem, "serious type mismatch" )
        canUse = otherCUDAMem->getDeviceNr() == getDeviceNr();
    }

    // Zero-Copy: we can use CUDA Host memory

    if ( other.getType() == MemoryType::CUDAHostMemory )
    {
        const CUDAHostMemory* otherCUDAHostMem = dynamic_cast<const CUDAHostMemory*>( &other );
        SCAI_ASSERT( otherCUDAHostMem, "serious type mismatch" )
        canUse = otherCUDAHostMem->getCUDAContext().getDeviceNr() == getDeviceNr();
    }

    SCAI_LOG_DEBUG( logger, *this << ": " << ( canUse ? "can use " : "can't use " )
                    << other )
    return canUse;
}

/* ----------------------------------------------------------------------------- */

CUDAStreamSyncToken* CUDAContext::getComputeSyncToken() const
{
    // ToDo: A possible problem might be that this CUDAContext is deleted before
    // synchronization has taken place. Solution: add a dummy routine where
    // one argument is bind to this context.
    return new CUDAStreamSyncToken( *this, StreamType::ComputeStream );
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAContext::getSyncToken() const
{
    return new CUDAStreamSyncToken( *this, StreamType::ComputeStream );
}

/* ----------------------------------------------------------------------------- */

CUDAStreamSyncToken* CUDAContext::getTransferSyncToken() const
{
    return new CUDAStreamSyncToken( *this, StreamType::TransferStream );
}

/* ----------------------------------------------------------------------------- */

#define SCAI_DEFAULT_DEVICE_NUMBER -1
#define SCAI_MAX_CUDA_DEVICES 8

static int getDefaultDeviceNr()
{
    int device = 0;
    common::Settings::getEnvironment( device, "SCAI_DEVICE" );
    return device;
}

/* ----------------------------------------------------------------------------- */
/*      Factory::Register - create( int )                                        */
/* ----------------------------------------------------------------------------- */

static std::weak_ptr<CUDAContext> mCUDAContext[SCAI_MAX_CUDA_DEVICES];

ContextPtr CUDAContext::create( int deviceNr )
{
    SCAI_LOG_INFO( logger, "create CUDAContext( device = " << deviceNr << " )" )

    int cudaDeviceNr = deviceNr;

    if ( cudaDeviceNr == SCAI_DEFAULT_DEVICE_NUMBER )
    {
        cudaDeviceNr = getDefaultDeviceNr();
    }

    SCAI_ASSERT_VALID_INDEX_ERROR( cudaDeviceNr, SCAI_MAX_CUDA_DEVICES, 
                                   "illegal device specified, MAX_DEVICES = " << SCAI_MAX_CUDA_DEVICES )

    std::shared_ptr<CUDAContext> context = std::shared_ptr<CUDAContext>();

    if ( mCUDAContext[cudaDeviceNr].expired() )
    {
        // create a new context for the device and return the shared pointer
        context = std::shared_ptr<CUDAContext>( new CUDAContext( cudaDeviceNr ) );
        // we keep a weak pointer so that we can return
        mCUDAContext[cudaDeviceNr] = context;
    }
    else
    {
        // the weak pointer to the device is still okay, so return a shared pointer for it
        context = mCUDAContext[cudaDeviceNr].lock();
    }

    return context;
}

/* ----------------------------------------------------------------------------- */

} /* end namespace hmemo */

} /* end namespace scai */
