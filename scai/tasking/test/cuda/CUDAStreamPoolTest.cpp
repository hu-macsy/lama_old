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
 * @brief Test of class CUDAStreamPool
 * @author: Thomas Brandes
 * @date 08.03.2016
 **/

#include <boost/test/unit_test.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/exception/Exception.hpp>
#include <scai/common/Settings.hpp>

#include <scai/tasking/cuda/CUDAStreamPool.hpp>
#include <scai/tasking/test/cuda/CUDAKernel.hpp>

#include <iostream>

using namespace scai;
using namespace tasking;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDAStreamPoolTest );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    int deviceNr = 0;

    common::Settings::getEnvironment( deviceNr, "SCAI_DEVICE" );

    common::CUDACtx myCuda( deviceNr );

    CUDAStreamPool& pool = CUDAStreamPool::getPool( myCuda );

    CUstream stream1 = pool.reserveStream( CUDAStreamPool::ComputeStream );
    CUstream stream2 = pool.reserveStream( CUDAStreamPool::TransferStream );

    pool.releaseStream( stream1 );

    BOOST_CHECK_THROW(
        {
            pool.releaseStream( stream1 );
        }, common::Exception );

    pool.releaseStream( stream2 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

