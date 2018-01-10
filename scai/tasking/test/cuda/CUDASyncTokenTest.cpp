/**
 * @file tasking/test/cuda/CUDASyncTokenTest.cpp
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
 * @brief Test of class CUDASyncToken
 * @author Thomas Brandes
 * @date 08.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
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
    common::CUDACtx myCuda( deviceNr );
    {
        CUDAStreamSyncToken token( myCuda, StreamType::TransferStream );
        // wait is done implicitly at end of this scope
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( destructorTest )
{
    int deviceNr = 0;

    common::Settings::getEnvironment( deviceNr, "SCAI_DEVICE" );

    common::CUDACtx myCuda( deviceNr );

    SyncToken* token = new CUDAStreamSyncToken( myCuda, StreamType::TransferStream );

    // Note: token not freed, no synchronization, no release

    BOOST_CHECK( ! CUDAStreamPool::getPool( myCuda ).isEmpty() );

    delete token;

    BOOST_CHECK( CUDAStreamPool::getPool( myCuda ).isEmpty() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( asyncTest )
{
    int deviceNr = 0;
    common::Settings::getEnvironment( deviceNr, "SCAI_DEVICE" );
    common::CUDACtx myCuda( deviceNr );
    common::CUDAAccess cudaAccess( myCuda );
    const int N = 100000;
    CUdeviceptr pointer = 0;
    size_t size = sizeof( float ) * N;
    SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, size ), "cuMemAlloc( size = " << size << " ) failed." )
    float* fpointer = reinterpret_cast<float*>( pointer );
    {
        CUDAStreamSyncToken token( myCuda, StreamType::ComputeStream );
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

