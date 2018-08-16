/**
 * @file CUDAAccessTest.cpp
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
 * @brief Test of class CUDAAccess
 * @author Thomas Brandes
 * @date 08.03.2016
 */

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
    scai::common::Exception );
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
    scai::common::Exception );
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

