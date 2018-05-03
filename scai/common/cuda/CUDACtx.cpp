/**
 * @file CUDACtx.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Implementation of methods for CUDA device.
 * @author Thomas Brandes
 * @date 08.03.2016
 */

#include <scai/common/cuda/CUDACtx.hpp>

#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAException.hpp>

#include <iostream>

#include <functional>

namespace scai
{

namespace common
{

static bool cudaInitialized = 0;

/* --------------------------------------------------------------------- */

CUDACtx::CUDACtx( int deviceNr )
{
    if ( !cudaInitialized )
    {
        unsigned int flags = 0;    // must be set to zero
        SCAI_CUDA_DRV_CALL( cuInit( flags ), "cuInit failed, probably no GPU devices available" )
        cudaInitialized = true;
    }

    mDeviceNr = deviceNr;   // no alternative is taken here
    SCAI_CUDA_DRV_CALL( cuDeviceGet( &mCUdevice, mDeviceNr ),
                        "cuDeviceGet device = " << mDeviceNr << " failed, probably not available" )

    SCAI_CUDA_DRV_CALL( cuCtxCreate( &mCUcontext, CU_CTX_SCHED_SPIN | CU_CTX_MAP_HOST, mCUdevice ),
                        "cuCtxCreate for " << mDeviceNr )

    SCAI_CUBLAS_CALL( cublasCreate( &mcuBLASHandle ), "Initialization of cuBLAS library" );

    SCAI_CUSPARSE_CALL( cusparseCreate( &mcuSparseHandle ), "Initialization of cuSparse library" );

#if ( CUDART_VERSION >= 7050 )
    SCAI_CUSOLVER_CALL( cusolverDnCreate( &mcuSolverDnHandle ), "Initialization of cuSolverDn library" )
    SCAI_CUSOLVER_CALL( cusolverSpCreate( &mcuSolverSpHandle ), "Initialization of cuSolverSp library" )
#endif
    CUcontext tmp; // temporary for last context, not necessary to save it
    SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( &tmp ), "could not pop context" )
}

/* --------------------------------------------------------------------- */

cusparseHandle_t CUDACtx::getcuSparseHandle() const
{
    return mcuSparseHandle;
}

/* --------------------------------------------------------------------- */

cublasHandle_t CUDACtx::getcuBLASHandle() const
{
    return mcuBLASHandle;
}

/* --------------------------------------------------------------------- */
#if ( CUDART_VERSION >= 7050 )

cusolverDnHandle_t CUDACtx::getcuSolverDnHandle() const
{
    return mcuSolverDnHandle;
}

/* --------------------------------------------------------------------- */

cusolverSpHandle_t CUDACtx::getcuSolverSpHandle() const
{
    return mcuSolverSpHandle;
}

#endif
/* --------------------------------------------------------------------- */

CUDACtx::~CUDACtx()
{
    // call added shutdown routines
    for ( size_t i = 0; i < mShutdownFunctions.size(); ++i )
    {
        mShutdownFunctions[i]();
    }

    // do not throw exceptions in destructor
    CUresult res = cuCtxPushCurrent( mCUcontext );

    if ( res == CUDA_ERROR_DEINITIALIZED )
    {
        // this might happen with other software, e.g. VampirTrace
        std::cerr << "Warn: CUDA already deinitialized" << std::endl;
        return;
    }
    else if ( res != CUDA_SUCCESS )
    {
        std::cerr << "Error: push context for device " << mDeviceNr << " failed." << std::endl;
        return;
    }

    // Be careful: cublasDestroy should be called within the current CUDA context

    if ( mcuBLASHandle )
    {
        cublasStatus_t error = cublasDestroy( mcuBLASHandle );

        if ( error != CUBLAS_STATUS_SUCCESS )
        {
            std::cerr << "Warn: could not destroy cublas handle, status = " << error << std::endl;
        }

        mcuBLASHandle = 0;
    }

    // Be careful: cusparseDestroy should be called within the current CUDA context

    if ( mcuSparseHandle )
    {
        cusparseStatus_t error = cusparseDestroy( mcuSparseHandle );

        if ( error != CUSPARSE_STATUS_SUCCESS )
        {
            std::cerr << "Warn: could not destroy cusparse handle, status = " << error << std::endl;
        }

        mcuSparseHandle = 0;
    }

#if ( CUDART_VERSION >= 7050 )
    // Be careful: cusolverDnDestroy should be called within the current CUDA context

    if ( mcuSolverDnHandle )
    {
        cusolverStatus_t error = cusolverDnDestroy( mcuSolverDnHandle );

        if ( error != CUSOLVER_STATUS_SUCCESS )
        {
            std::cerr << "Warn: could not destroy cusolverDn handle, status = " << error << std::endl;
        }

        mcuSolverDnHandle = 0;
    }

    // Be careful: cusolverSpDestroy should be called within the current CUDA context

    if ( mcuSolverSpHandle )
    {
        cusolverStatus_t error = cusolverSpDestroy( mcuSolverSpHandle );

        if ( error != CUSOLVER_STATUS_SUCCESS )
        {
            std::cerr << "Warn: could not destroy cusolverSp handle, status = " << error << std::endl;
        }

        mcuSolverSpHandle = 0;
    }

#endif

    SCAI_CUDA_DRV_CALL_NOTHROW( cuCtxDestroy( mCUcontext ),
                                "Error: destroy context for device " << mDeviceNr << " failed." )
}

/* ----------------------------------------------------------------------- */

void CUDACtx::addShutdown( std::function<void()> routine )
{
    mShutdownFunctions.push_back( routine );
}

}  // namespace common

}  // namespace scai
