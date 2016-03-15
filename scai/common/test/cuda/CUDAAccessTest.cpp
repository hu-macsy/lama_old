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

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>

#include <scai/common/test/cuda/CUDAKernel.hpp>

#include <iostream>


/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDAAccessTest );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( accessTest )
{
    int nr = 0;

    scai::common::Settings::getEnvironment( nr, "SCAI_DEVICE" );

    scai::common::CUDACtx myCuda( nr );

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

        scai::common::CUDAAccess tmpAccess( myCuda );

        SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, size ), "cuMemAlloc( size = " << size << " ) failed." )
    
        float* fpointer = reinterpret_cast<float*>( pointer );
    
        init( fpointer, N, 3.0 );

        float s = sum( fpointer, N );

        BOOST_CHECK_CLOSE( 3.0f * N, s, 0.01 );

        SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << pointer << " ) failed" )
    }

    // checkt that access is released at end of previous scope

    BOOST_CHECK_THROW( 
        {
            SCAI_CHECK_CUDA_ACCESS
        }, 
        scai::common::Exception )
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getCurrentDeviceTest )
{
    using namespace scai::common;

    int nr = 0;

    Settings::getEnvironment( nr, "SCAI_DEVICE" );

    CUDACtx myCuda( nr );

    BOOST_CHECK_THROW( 
        {
            CUDAAccess::getCurrentCUDACtx();
        }, 
        Exception );

    {
        CUDAAccess access( myCuda );
        const CUDACtx& current = CUDAAccess::getCurrentCUDACtx();
        BOOST_CHECK_EQUAL( &myCuda, &current );
    }

    BOOST_CHECK_THROW( 
        {
            CUDAAccess::getCurrentCUDACtx();
        }, 
        Exception );

}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( stackAccessTest )
{
    using namespace scai::common;

    int nr = 0;

    Settings::getEnvironment( nr, "SCAI_DEVICE" );

    // create 3 different contexts for one and the same device 

    CUDACtx myCuda1( nr );
    CUDACtx myCuda2( nr );
    CUDACtx myCuda3( nr );

    {
        CUDAAccess access1( myCuda1 );

        {
            CUDAAccess access2( myCuda2 );

            {
                CUDAAccess access3( myCuda3 );
            }

            const CUDACtx& current = CUDAAccess::getCurrentCUDACtx();
            BOOST_CHECK_EQUAL( &myCuda2, &current );
        }

        const CUDACtx& current = CUDAAccess::getCurrentCUDACtx();
        BOOST_CHECK_EQUAL( &myCuda1, &current );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

