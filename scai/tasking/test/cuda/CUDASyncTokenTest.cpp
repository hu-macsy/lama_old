/**
 * @file CUDSyncTokenTest.cpp
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
 * @brief Test of class CUDASyncToken
 * @author: Thomas Brandes
 * @date 08.03.2016
 **/

#include <boost/test/unit_test.hpp>

#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/cuda/CUDADevice.hpp>
#include <scai/common/Settings.hpp>

#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>
#include <scai/tasking/cuda/CUDAStreamPool.hpp>

#include <scai/tasking/test/cuda/CUDAKernel.hpp>

#include <iostream>

using namespace scai;
using namespace tasking;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDASyncTokenTest );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    int deviceNr = 0;

    common::Settings::getEnvironment( deviceNr, "SCAI_DEVICE" );

    common::CUDADevice myCuda( deviceNr );

    {
        CUDAStreamSyncToken token( myCuda, CUDAStreamSyncToken::TransferStream );

        // wait is done implicitly at end of this scope
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( destructorTest )
{
    int deviceNr = 0;

    common::Settings::getEnvironment( deviceNr, "SCAI_DEVICE" );

    BOOST_CHECK_THROW(
        {
            common::CUDADevice myCuda( deviceNr );

            {
                new CUDAStreamSyncToken( myCuda, CUDAStreamSyncToken::TransferStream );
        
                // Note: token not freed, no synchronization, no release
            }

            // destructor of device will free the stream pool and detect the unreleased token

        }, common::Exception );

     // Note: there might be an additional ERROR logging message at the end of the test run
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( asyncTest )
{
    int deviceNr = 0;

    common::Settings::getEnvironment( deviceNr, "SCAI_DEVICE" );

    common::CUDADevice myCuda( deviceNr );

    common::CUDAAccess cudaAccess( myCuda );

    const int N = 100000;

    CUdeviceptr pointer = 0;

    size_t size = sizeof( float ) * N;

    SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, size ), "cuMemAlloc( size = " << size << " ) failed." )

    float* fpointer = reinterpret_cast<float*>( pointer );

    {
        CUDAStreamSyncToken token( myCuda, CUDAStreamSyncToken::ComputeStream );

        // launch kernel asynchronously

        init( fpointer, N, 3.0 );

        // wait is done here implicitly
    }

    float s = sum( fpointer, N );

    BOOST_CHECK_CLOSE( 3.0f * N, s, 0.01 );

    SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << pointer << " ) failed" )
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

