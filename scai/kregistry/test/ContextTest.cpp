/**
 * @file kregistry/test/ContextTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief ToDo: Missing description in ./kregistry/test/ContextTest.cpp
 * @author Thomas Brandes
 * @date 16.10.2015
 */
#include <boost/test/unit_test.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace scai::kregistry;

using scai::common::context;

static double add1( const double x )
{
    return x + 1.0;
}

static double minus1( const double x )
{
    return x - 1.0;
}

/** Trait to handle function double ( fn ) ( double ) in KernelRegistry. */

struct UnaryAddTrait
{
    typedef double ( *FuncType ) ( double );
    static const char* getId()
    {
        return "add";
    }
};

struct UnaryMinusTrait
{
    typedef double ( *FuncType ) ( double );
    static const char* getId()
    {
        return "minus";
    }
};

BOOST_AUTO_TEST_CASE( ContextTest )
{
    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD;
    // register context::Host
    KernelRegistry::set<UnaryAddTrait>( add1, context::Host, flag );
    KernelRegistry::set<UnaryMinusTrait>( minus1, context::Host, flag );
    // register context::CUDA
    KernelRegistry::set<UnaryAddTrait>( add1, context::CUDA, flag );
    KernelRegistry::set<UnaryMinusTrait>( minus1, context::CUDA, flag );
    // register context::MIC, only add1
    KernelRegistry::set<UnaryAddTrait>( add1, context::MIC, flag );
    // register context::UserContext, only minus1
    KernelRegistry::set<UnaryMinusTrait>( minus1, context::UserContext, flag );
    KernelRegistry::printAll();
    KernelTraitContextFunction<UnaryAddTrait> add;
    KernelTraitContextFunction<UnaryMinusTrait> minus;
    // add can be alled at context::MIC
    BOOST_CHECK_EQUAL( context::MIC, add.validContext( context::MIC ) );
    // minus must be called at context::Host
    BOOST_CHECK_EQUAL( context::Host, minus.validContext( context::MIC ) );
    // add, minus can be called together @ CUDA
    BOOST_CHECK_EQUAL( context::CUDA, add.validContext( minus, context::CUDA ) );
    BOOST_CHECK_EQUAL( context::Host, add.validContext( minus, context::MIC ) );
    BOOST_CHECK_EQUAL( context::Host, add.validContext( minus, context::UserContext ) );
    BOOST_CHECK_EQUAL( context::Host, add.validContext( minus, context::Host ) );
}

