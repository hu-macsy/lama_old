/**
 * @file CUDACtxTest.cpp
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
 * @brief Test of class CUDACtx
 * @author: Thomas Brandes
 * @date 08.03.2016
 **/

#include <boost/test/unit_test.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>

#include <scai/common/test/cuda/CUDAKernel.hpp>
#include <scai/common/bind.hpp>

#include <iostream>


static void inc( int* val ) 
{
    *val += 1;
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDACtxTest );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    int blas_version = 0;

    BOOST_CHECK_THROW (
    {
        SCAI_CUBLAS_CALL( cublasGetVersion( 0, &blas_version ), "get version" );

    }, scai::common::Exception );

    int nr = 0;

    scai::common::Settings::getEnvironment( nr, "SCAI_DEVICE" );

    scai::common::CUDACtx myCuda( nr );

    BOOST_CHECK_EQUAL( nr, myCuda.getDeviceNr() );

    cublasHandle_t blasHandle = myCuda.getcuBLASHandle();

    blas_version = 0;

    SCAI_CUBLAS_CALL( cublasGetVersion( blasHandle, &blas_version ), "get version" );

    BOOST_CHECK( blas_version > 0 );

    cusparseHandle_t sparseHandle = myCuda.getcuSparseHandle();

    int sparse_version = 0;

    SCAI_CUSPARSE_CALL( cusparseGetVersion( sparseHandle, &sparse_version ), "get version" );

    BOOST_CHECK( sparse_version > 0 );

    cusparseMatDescr_t dA = NULL;
    cusparseDiagType_t dtype = cusparseDiagType_t ( 10 );    //  no legal enum value

    SCAI_CUSPARSE_CALL( cusparseCreateMatDescr( &dA ), "create mat" );

    BOOST_CHECK_THROW (
    {
        SCAI_CUSPARSE_CALL( cusparseSetMatDiagType( dA, dtype ), "set unknown diag type" );

    }, scai::common::Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( shutdownTest )
{
    int nr = 0;

    int val = 0;

    scai::common::Settings::getEnvironment( nr, "SCAI_DEVICE" );

    {
        scai::common::CUDACtx myCuda( nr );

        myCuda.addShutdown( scai::common::bind( &inc, &val ) );

        BOOST_CHECK_EQUAL( 0, val );

        myCuda.addShutdown( scai::common::bind( &inc, &val ) );

        // shutdown routine will only be called with destructor
    }

    // inc should have been called twice by destructor

    BOOST_CHECK_EQUAL( 2, val );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

