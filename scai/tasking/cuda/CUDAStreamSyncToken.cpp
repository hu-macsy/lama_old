/**
 * @file CUDAStreamSyncToken.cpp
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
 * @brief Implementation of methods for class CUDAStreamSyncToken.
 * @author Jiri Kraus
 * @date 20.05.2011
 */

// hpp
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <scai/tasking/cuda/CUDAStreamPool.hpp>

// internal scai libraries

#include <scai/common/macros/assert.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

namespace scai
{

using common::CUDACtx;
using common::CUDAAccess;

namespace tasking
{

CUDAStreamSyncToken::CUDAStreamSyncToken( const CUDACtx& cuda, const StreamType type ) : mCUDA( cuda )
{
    // take a CUDA stream from a pool, that might be allocated at first use
    mStream = CUDAStreamPool::getPool( mCUDA ).reserveStream( type );
    mEvent = 0;
    SCAI_LOG_DEBUG( logger, "StreamSyncToken for " << mCUDA.getDeviceNr() << " generated, stream = " << mStream )
}

CUDAStreamSyncToken::~CUDAStreamSyncToken()
{
    wait();
    CUDAStreamPool::getPool( mCUDA ).releaseStream( mStream );
}

void CUDAStreamSyncToken::wait()
{
    if ( isSynchronized() )
    {
        return;
    }

    // Note: do not throw exceptions here

    SCAI_LOG_DEBUG( logger, "wait on CUDA stream synchronization" )
    {
        common::CUDAAccess tmpAccess( mCUDA );

        if ( mEvent != 0 )
        {
            SCAI_CUDA_DRV_CALL_NOTHROW( cuEventSynchronize( mEvent ), "cuEventSynchronize( " << mEvent << " ) failed." )
            SCAI_CUDA_DRV_CALL_NOTHROW( cuEventDestroy( mEvent ), "cuEventDestroy( " << mEvent << " ) failed." );
            mEvent = 0;
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "synchronize with stream " << mStream );
            SCAI_CUDA_DRV_CALL_NOTHROW( cuStreamSynchronize( mStream ), "cuStreamSynchronize( " << mStream << " ) failed." );
            SCAI_LOG_DEBUG( logger, "synchronized with stream " << mStream );
        }

        // finally called functions might also need the context, e.g. unbindTexture
        setSynchronized();
    }
}

bool CUDAStreamSyncToken::probe() const
{
    if ( isSynchronized() )
    {
        return true;
    }

    return probeEvent( mEvent );
}

cudaStream_t CUDAStreamSyncToken::getCUDAStream() const
{
    return ( cudaStream_t ) mStream;
}

bool CUDAStreamSyncToken::probeEvent( const CUevent& stopEvent ) const
{
    SCAI_ASSERT( stopEvent != 0, "probe on invalid event" )
    common::CUDAAccess tmpAccess ( mCUDA );
    CUresult result = cuEventQuery( stopEvent );

    if ( result != CUDA_SUCCESS && result != CUDA_ERROR_NOT_READY )
    {
        SCAI_CUDA_DRV_CALL( result, "cuEventQuery failed for CUevent " << stopEvent << '.' );
    }

    return result == CUDA_SUCCESS;
}

bool CUDAStreamSyncToken::queryEvent( const CUevent event ) const
{
    common::CUDAAccess cudaAccess ( mCUDA );
    CUresult result = cuEventQuery( event );

    if ( result != CUDA_SUCCESS || result != CUDA_ERROR_NOT_READY )
    {
        SCAI_CUDA_DRV_CALL( result, "cuEventQuery failed for CUevent " << event << '.' );
    }

    return result == CUDA_SUCCESS;
}

void CUDAStreamSyncToken::synchronizeEvent( const CUevent event ) const
{
    common::CUDAAccess cudaAccess ( mCUDA );
    SCAI_CUDA_DRV_CALL( cuEventSynchronize( event ), "cuEventSynchronize failed for CUevent " << event << '.' )
}

CUDAStreamSyncToken* CUDAStreamSyncToken::getCurrentSyncToken()
{
    SyncToken* syncToken = SyncToken::getCurrentSyncToken();

    if ( syncToken == NULL )
    {
        return NULL;
    }

    // make a dynamic CAST
    CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );

    // If the current sync token is not a CUDA stream token it is very likely an error

    if ( cudaStreamSyncToken == NULL )
    {
        SCAI_LOG_ERROR( logger, "Current sync token = " << *syncToken << " not CUDAStreamSyncToken as expected" )
    }

    // But might not be too serious so probably NULL results in synchronous execution
    return cudaStreamSyncToken;
}

void CUDAStreamSyncToken::writeAt( std::ostream& stream ) const
{
    stream << "CUDAStreamSyncToken( ";

    stream << ", mStream = " << mStream;

    stream << ", synchronized = " << isSynchronized() << ")";
}

} /* end namespace tasking */

} /* end namespace scai */
