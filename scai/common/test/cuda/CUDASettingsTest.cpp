/**
 * @file CUDASettingsTest.cpp
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
 * @brief Test of class CUDASettings
 * @author Thomas Brandes
 * @date 31.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDASettings.hpp>

#include <scai/common/Settings.hpp>

#include <scai/common/test/cuda/CUDAKernel.hpp>

#include <iostream>
#include <functional>

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

