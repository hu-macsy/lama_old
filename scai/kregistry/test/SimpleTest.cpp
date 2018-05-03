/**
 * @file kregistry/test/SimpleTest.cpp
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
 * @brief ToDo: Missing description in ./kregistry/test/SimpleTest.cpp
 * @author Thomas Brandes
 * @date 16.10.2015
 */
#include <boost/test/unit_test.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace scai::common;
using namespace scai::kregistry;

static void dummyRoutine()
{
}

BOOST_AUTO_TEST_CASE( SimpleTest )
{
    // This simple test registers a function in the kernel registry and uses it later
    KernelRegistry::set( dummyRoutine, "dummy", ContextType::Host, KernelRegistry::KERNEL_ADD );
    KernelContextFunction<void(* )()> f( "dummy" );
    f[ ContextType::Host ]();  // just call it
    // throw exception if called for CUDA, not registered
    BOOST_CHECK_THROW(
    {
        f[ ContextType::CUDA ]();

    }, KernelRegistryException );
    BOOST_CHECK_THROW(
    {
        KernelContextFunction<void(* )()> g( "dummy1" );  // wrong name

    }, KernelRegistryException );
    BOOST_CHECK_THROW(
    {
        KernelContextFunction<int(* )()> g( "dummy" );  // wrong signature

    }, KernelRegistryException );
}
