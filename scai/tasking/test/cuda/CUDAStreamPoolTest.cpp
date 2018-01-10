/**
 * @file tasking/test/cuda/CUDAStreamPoolTest.cpp
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
 * @brief Test of class CUDAStreamPool
 * @author Thomas Brandes
 * @date 08.03.2016
 */

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
    CUstream stream1 = pool.reserveStream( StreamType::ComputeStream );
    CUstream stream2 = pool.reserveStream( StreamType::TransferStream );
    pool.releaseStream( stream1 );
    BOOST_CHECK_THROW(
    {
        pool.releaseStream( stream1 );
    }, common::Exception );
    pool.releaseStream( stream2 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

