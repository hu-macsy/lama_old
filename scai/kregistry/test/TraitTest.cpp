/**
 * @file kregistry/test/TraitTest.cpp
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
 * @brief ToDo: Missing description in ./kregistry/test/TraitTest.cpp
 * @author Thomas Brandes
 * @date 16.10.2015
 */
#include <boost/test/unit_test.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace scai::common;
using namespace scai::kregistry;

static int dummyRoutine()
{
    return 15;
}

/** Trait structure for registration of int routine with name "MyDummy" */

struct TraitDummyRoutine
{
    typedef int ( *FuncType ) ();    // signature of the function
    static const char* getId()
    {
        return "MyDummy";
    }
};

BOOST_AUTO_TEST_CASE( TraitTest )
{
    // Same as simple test but uses a Trait for registration
    // The trait avoids misspelling of the routine name and the signature
    KernelRegistry::set<TraitDummyRoutine>( dummyRoutine, ContextType::CUDA, KernelRegistry::KERNEL_ADD );
    KernelTraitContextFunction<TraitDummyRoutine> f;
    int x = f[ ContextType::CUDA ]();  // just call it
    BOOST_CHECK_EQUAL( 15, x );
    // misspelling of name or signature is no more possible here, so no further test for failure
}

