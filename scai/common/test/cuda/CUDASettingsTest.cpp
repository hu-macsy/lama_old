/**
 * @file CUDASettingsTest.cpp
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
 * @brief Test of class CUDASettings
 * @author: Thomas Brandes
 * @date 31.03.2016
 **/

#include <boost/test/unit_test.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDASettings.hpp>

#include <scai/common/Settings.hpp>

#include <scai/common/test/cuda/CUDAKernel.hpp>
#include <scai/common/bind.hpp>

#include <iostream>

using namespace scai;
using namespace common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDASettingsTest );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( useTextureTest )
{

    int nr = 0;

    scai::common::Settings::getEnvironment( nr, "SCAI_DEVICE" );
    scai::common::CUDACtx myCuda( nr );
    scai::common::CUDAAccess tmpAccess( myCuda );

    bool useTexture = CUDASettings::useTexture();
    bool useSharedMem = CUDASettings::useSharedMem();

    // cannot be rewritten this way

    Settings::putEnvironment( "SCAI_CUDA_USE_TEXTURE", not useTexture );

    BOOST_CHECK_EQUAL( useTexture, CUDASettings::useTexture() );

    // but can be rewritten in this way

    CUDASettings::set( useSharedMem, not useTexture );

    BOOST_CHECK_EQUAL( not useTexture, CUDASettings::useTexture() );
    BOOST_CHECK_EQUAL( useSharedMem, CUDASettings::useSharedMem() );

    CUDASettings::set( not useSharedMem, useTexture );
    BOOST_CHECK_EQUAL( not useSharedMem, CUDASettings::useSharedMem() );
}

BOOST_AUTO_TEST_CASE( blockSizeTest )
{
    int nr = 0;

    scai::common::Settings::getEnvironment( nr, "SCAI_DEVICE" );
    scai::common::CUDACtx myCuda( nr );
    scai::common::CUDAAccess tmpAccess( myCuda );

    // This is the default block size for parallelism 

    int bsize = CUDASettings::getBlockSize();

    // For large sized problems we get the full block size
 
    int bsize100 = CUDASettings::getBlockSize( 100 * bsize );
    BOOST_CHECK_EQUAL( bsize, bsize100 );

    // For small sized problems we get the half block size

    int bsize2 = CUDASettings::getBlockSize( bsize );
    BOOST_CHECK_EQUAL( bsize, 2 * bsize2 );
}

BOOST_AUTO_TEST_SUITE_END();

