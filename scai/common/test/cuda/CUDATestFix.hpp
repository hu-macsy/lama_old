/**
 * @file CUDATestFix.hpp
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
 * @brief Fixture for CUDA Tests
 * @author: Thomas Brandes
 * @date 08.03.2016
 **/

#pragma once

#include <scai/common/Settings.hpp>

#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

using namespace scai;

/* --------------------------------------------------------------------- */

struct CUDAFix
{   
    CUDAFix()
    {
        unsigned int flags = 0;    // must be set to zero

        SCAI_CUDA_DRV_CALL( cuInit( flags ), "cuInit failed, probably no GPU devices available" )

        int deviceNr = 0;  // take this as default

        scai::common::Settings::getEnvironment( deviceNr, "SCAI_DEVICE" );  // can be set
        
        SCAI_CUDA_DRV_CALL( cuDeviceGet( &mCUdevice, deviceNr ), "cuDeviceGet device " << deviceNr );
        
        SCAI_CUDA_DRV_CALL( cuCtxCreate( &mCUcontext, CU_CTX_SCHED_SPIN | CU_CTX_MAP_HOST, mCUdevice ),
                            "cuCtxCreate for " << deviceNr )
        
        CUcontext tmp; // temporary for last context, not necessary to save it

        SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( &tmp ), "could not pop context" )
    }

    ~CUDAFix()
    {
        SCAI_CUDA_DRV_CALL( cuCtxPushCurrent( mCUcontext ), "push context failed" );
        SCAI_CUDA_DRV_CALL( cuCtxDestroy( mCUcontext ), "cuCtxDestroy failed" )
    }

    CUdevice mCUdevice;
    CUcontext mCUcontext;
};

