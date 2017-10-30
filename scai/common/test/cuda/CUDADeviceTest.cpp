/**
 * @file test/cuda/CUDADeviceTest.cpp
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
 * @brief Test of class CUDACtx
 * @author Thomas Brandes
 * @date 08.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>

#include <scai/common/test/cuda/CUDAKernel.hpp>

#include <iostream>
#include <functional>


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
        myCuda.addShutdown( std::bind( &inc, &val ) );
        BOOST_CHECK_EQUAL( 0, val );
        myCuda.addShutdown( std::bind( &inc, &val ) );
        // shutdown routine will only be called with destructor
    }
    // inc should have been called twice by destructor
    BOOST_CHECK_EQUAL( 2, val );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

