/**
 * @file CUDAStreamPool.cpp
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
 * @brief Managing multiple streams for a CUDA context
 * @author Thomas Brandes
 * @date 09.07.2016
 */

#include <scai/tasking/cuda/CUDAStreamPool.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/bind.hpp>

#include <map>

namespace scai
{

namespace tasking
{

/* -----------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( CUDAStreamPool::logger, "CUDAStreamPool" )

/* -----------------------------------------------------------------------------*/

CUstream CUDAStreamPool::reserveStream( const StreamType type )
{
    if ( type == ComputeStream )
    {
        mComputeReservations++;

        SCAI_LOG_INFO( logger, "reserved computed stream " << mComputeStream << ", #reservations = " << mComputeReservations )

        return mComputeStream;
    }
    else
    {
        mTransferReservations++;

        SCAI_LOG_INFO( logger, "reserved transfer stream " << mTransferStream << ", #reservations = " << mTransferReservations )

        return mTransferStream;
    }
}

/* -----------------------------------------------------------------------------*/

void CUDAStreamPool::releaseStream( CUstream stream )
{
    if ( stream == mComputeStream )
    {
        SCAI_ASSERT_GT( mComputeReservations, 0, "Compute stream not reserved" )

        mComputeReservations--;

        SCAI_LOG_INFO( logger, "released computed stream " << mComputeStream << ", #reservations = " << mComputeReservations )
    }
    else if ( stream == mTransferStream )
    {
        SCAI_ASSERT_GT( mTransferReservations, 0, "Memory transfer stream not reserved" )

        mTransferReservations--;

        SCAI_LOG_INFO( logger, "released transfer stream " << mTransferStream << ", #reservations = " << mTransferReservations )
    }
    else
    {
        COMMON_THROWEXCEPTION( "CUstream " << stream << " not of this stream pool." )
    }
}

/* -----------------------------------------------------------------------------*/

CUDAStreamPool::CUDAStreamPool( const common::CUDACtx& cuda ) : mCUDA( cuda )
{
    SCAI_LOG_INFO( logger, "CUDAStreamPool( device = " << mCUDA.getDeviceNr() << " )" )

    common::CUDAAccess tmpAccess( mCUDA ); 

    int flags = 0; // must be 0 by specification of CUDA driver API

    SCAI_CUDA_DRV_CALL( cuStreamCreate( &mTransferStream, flags ), "cuStreamCreate for transfer failed" )
    SCAI_CUDA_DRV_CALL( cuStreamCreate( &mComputeStream, flags ), "cuStreamCreate for compute failed" );

    mComputeReservations = 0;
    mTransferReservations = 0;
}

/* -----------------------------------------------------------------------------*/

CUDAStreamPool::~CUDAStreamPool()
{
    SCAI_LOG_INFO( logger, "~CUDAStreamPool( device = " << mCUDA.getDeviceNr() << " )" )

    common::CUDAAccess tmpAccess( mCUDA ); 

    // No exceptions in destructor !!

    SCAI_ASSERT_EQUAL( mComputeReservations, 0, "Not all compute streams released" )
    SCAI_ASSERT_EQUAL( mTransferReservations, 0, "Not all transfer streams released" )

    SCAI_CUDA_DRV_CALL( cuStreamSynchronize( mComputeStream ), "cuStreamSynchronize for compute failed" )
    SCAI_CUDA_DRV_CALL( cuStreamDestroy( mComputeStream ), "cuStreamDestroy for compute failed" )
    SCAI_CUDA_DRV_CALL( cuStreamSynchronize( mTransferStream ), "cuStreamSynchronize for transfer failed" )
    SCAI_CUDA_DRV_CALL( cuStreamDestroy( mTransferStream ), "cuStreamDestroy for transfer failed" )
}

typedef std::map<CUcontext, CUDAStreamPool*>  PoolMap;

/* -----------------------------------------------------------------------------*/

static PoolMap& getPoolMap()
{
    static PoolMap* thePoolMap = NULL;

    if ( thePoolMap == NULL )
    {
        thePoolMap = new PoolMap();
    }

    return *thePoolMap;
}

/* -----------------------------------------------------------------------------*/

CUDAStreamPool& CUDAStreamPool::getPool( const common::CUDACtx& cuda )
{
    PoolMap& poolMap = getPoolMap();

    PoolMap::iterator it = poolMap.find( cuda.getCUcontext() );
   
    if ( it == poolMap.end() )
    {
        CUDAStreamPool* pool = new CUDAStreamPool( cuda );

        // map takes ownership of pool 

        poolMap.insert( std::pair<CUcontext, CUDAStreamPool*>( cuda.getCUcontext(), pool) );
  
        // ATTENTION / pitfall
        // Pool must be freed before cuda is destroyed 
        // solution: Add shutdown routine to the CUDA device

        common::CUDACtx& cuda1 = const_cast< common::CUDACtx& >( cuda );

        cuda1.addShutdown( common::bind( &CUDAStreamPool::freePool, common::cref( cuda ) ) );

        return *pool;
    }
    else
    {
        return *it->second;
    }
}

/* -----------------------------------------------------------------------------*/

void CUDAStreamPool::freePool( const common::CUDACtx& cuda )
{
    SCAI_LOG_INFO( logger, "freePool: pool for CUDA device " << cuda.getDeviceNr() )

    PoolMap& poolMap = getPoolMap();

    PoolMap::iterator it = poolMap.find( cuda.getCUcontext() );

    if ( it != poolMap.end() )
    {
        delete it->second;

        poolMap.erase( it );
    }
    else
    {
        SCAI_LOG_INFO( logger, "freePool: pool for CUDA device " << cuda.getDeviceNr() << " released twice or never created" )
    }
}

/* -----------------------------------------------------------------------------*/

} /* end namespace tasking */

} /* end namespace scai */