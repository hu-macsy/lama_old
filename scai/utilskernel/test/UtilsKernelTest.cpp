/**
 * @file utilskernel/test/UtilsKernelTest.cpp
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

int main( int argc, char* argv[] )
{
    SCAI_LOG_THREAD( "main" )
    SCAI_REGION( "Main.UtilsKernelTest" )
    return scai::testsupport::hmemoTestMain(argc, argv);
}
