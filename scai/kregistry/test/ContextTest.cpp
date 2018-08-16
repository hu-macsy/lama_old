/**
 * @file kregistry/test/ContextTest.cpp
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
 * @brief ToDo: Missing description in ./kregistry/test/ContextTest.cpp
 * @author Thomas Brandes
 * @date 16.10.2015
 */
#include <boost/test/unit_test.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace kregistry;
using common::ContextType;

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
    // register ContextType::Host
    KernelRegistry::set<UnaryAddTrait>( add1, ContextType::Host, flag );
    KernelRegistry::set<UnaryMinusTrait>( minus1, ContextType::Host, flag );
    // register ContextType::CUDA
    KernelRegistry::set<UnaryAddTrait>( add1, ContextType::CUDA, flag );
    KernelRegistry::set<UnaryMinusTrait>( minus1, ContextType::CUDA, flag );
    // register ContextType::UserContext, only minus1
    KernelRegistry::set<UnaryMinusTrait>( minus1, ContextType::UserContext, flag );
    KernelRegistry::printAll();
    KernelTraitContextFunction<UnaryAddTrait> add;
    KernelTraitContextFunction<UnaryMinusTrait> minus;
    // add, minus can be called together @ CUDA
    BOOST_CHECK_EQUAL( ContextType::CUDA, add.validContext( minus, ContextType::CUDA ) );
    BOOST_CHECK_EQUAL( ContextType::Host, add.validContext( minus, ContextType::UserContext ) );
    BOOST_CHECK_EQUAL( ContextType::Host, add.validContext( minus, ContextType::Host ) );
}

