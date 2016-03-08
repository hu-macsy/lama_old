/**
 * @file CUDAAccessTest.cpp
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
 * @brief Test of class CUDAAccess
 * @author: Thomas Brandes
 * @date 08.03.2016
 **/

#include <boost/test/unit_test.hpp>

#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/test/cuda/CUDAKernel.hpp>

#include <iostream>

using namespace scai;

/* --------------------------------------------------------------------- */

struct CUDAFix
{   
    CUDAFix()
    {
        std::cout << "Init cuda" << std::endl;

        unsigned int flags = 0;    // must be set to zero

        SCAI_CUDA_DRV_CALL( cuInit( flags ), "cuInit failed, probably no GPU devices available" )

        int deviceNr = 0;
        
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

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CommonCUDATest );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( accessTest )
{
    CUDAFix myCuda;

    const int N = 100;

    CUdeviceptr pointer = 0;

    size_t size = sizeof( float ) * N;

    // driver call without access must throw exception

    BOOST_CHECK_THROW( 
        {
            SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, size ), "cuMemAlloc( size = " << size << " ) failed." )
        }, 
        scai::common::Exception )

    {
        // create an access for CUDA calls 

        scai::common::CUDAAccess tmpAccess( myCuda.mCUcontext );

        SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, size ), "cuMemAlloc( size = " << size << " ) failed." )
    
        float* fpointer = reinterpret_cast<float*>( pointer );
    
        init( fpointer, N, 3.0 );

        float s = sum( fpointer, N );

        BOOST_CHECK_CLOSE( 3.0f * N, s, 0.01 );

        SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << pointer << " ) failed" )
    }

    BOOST_CHECK_THROW( 
        {
            SCAI_CHECK_CUDA_ACCESS
        }, 
        scai::common::Exception )
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

