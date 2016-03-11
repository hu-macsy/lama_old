/**
 * @file CUDADevice.cpp
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
 * @brief Implementation of methods for CUDA device.
 * @author: Thomas Brandes
 * @date 08.03.2016
 **/

#include <scai/common/cuda/CUDADevice.hpp>

#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <iostream>

namespace scai
{

namespace common
{

static bool cudaInitialized = 0;

/* --------------------------------------------------------------------- */

CUDADevice::CUDADevice( int deviceNr )
{
    if ( !cudaInitialized )
    {
        unsigned int flags = 0;    // must be set to zero

        SCAI_CUDA_DRV_CALL( cuInit( flags ), "cuInit failed, probably no GPU devices available" )
        cudaInitialized = true;
    }

    mDeviceNr = deviceNr;   // no alternative is taken here

    SCAI_CUDA_DRV_CALL( cuDeviceGet( &mCUdevice, mDeviceNr ), 
                        "cuDeviceGet device = " << mDeviceNr << " failed, probably not available" );
        
    SCAI_CUDA_DRV_CALL( cuCtxCreate( &mCUcontext, CU_CTX_SCHED_SPIN | CU_CTX_MAP_HOST, mCUdevice ),
                        "cuCtxCreate for " << mDeviceNr )
        
    SCAI_CUBLAS_CALL( cublasCreate( &mcuBLASHandle ), "Initialization of cuBLAS library" );

    CUcontext tmp; // temporary for last context, not necessary to save it

    SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( &tmp ), "could not pop context" )
}

/* --------------------------------------------------------------------- */

cublasHandle_t CUDADevice::getcuBLASHandle() const
{
    return mcuBLASHandle;
}

/* --------------------------------------------------------------------- */

CUDADevice::~CUDADevice()
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

    res = cuCtxDestroy( mCUcontext );

    if ( res != CUDA_SUCCESS )
    {
        std::cerr << "Error: destroy context for device " << mDeviceNr << " failed." << std::endl;
        return;
    }
}

/* ----------------------------------------------------------------------- */

void CUDADevice::addShutdown( common::function<void()> routine )
{
    mShutdownFunctions.push_back( routine );
}

}  // namespace common

}  // namespace scai
