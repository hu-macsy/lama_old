/**
 * @file utilskernel/test/UtilsKernelTest.cpp
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
 * @brief ToDo: Missing description in ./utilskernel/test/UtilsKernelTest.cpp
 * @author Eric Schricker
 * @date 18.02.2016
 */

#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE utilskernelTest
#define BOOST_TEST_NO_MAIN

#include <scai/testsupport/hmemoTestMain.hpp>
#include <scai/logging.hpp>
#include <scai/tracing.hpp>

scai::hmemo::ContextPtr testContext;

int main( int argc, char* argv[] )
{
    testContext = scai::hmemo::Context::getContextPtr();
    SCAI_LOG_THREAD( "main" )
    SCAI_REGION( "Main.UtilsKernelTest" )
    return scai::testsupport::hmemoTestMain(argc, argv);
}
