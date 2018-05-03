/**
 * @file kregistry/test/ReplaceTest.cpp
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
 * @brief ToDo: Missing description in ./kregistry/test/ReplaceTest.cpp
 * @author Thomas Brandes
 * @date 16.10.2015
 */
#include <boost/test/unit_test.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace scai::kregistry;
using namespace scai::common;

template<typename ValueType>
static ValueType add1( const ValueType x )
{
    return x + static_cast<ValueType>( 1 );
}

template<typename ValueType>
static ValueType minus1( const ValueType x )
{
    return x - static_cast<ValueType>( 1 );
}

/** Trait to handle function ValueType ( fn ) ( ValueType ) in KernelRegistry. */

template<typename ValueType>
struct UnaryOpTrait
{
    typedef ValueType ( *FuncType ) ( ValueType );
    static const char* getId()
    {
        return "UnaryOp";
    }
};

BOOST_AUTO_TEST_CASE( ReplaceTest )
{
    // register unary operator for double
    KernelRegistry::set<UnaryOpTrait<double> >( add1<double>, ContextType::Host, KernelRegistry::KERNEL_ADD );
    KernelRegistry::set<UnaryOpTrait<double> >( minus1<double>, ContextType::Host, KernelRegistry::KERNEL_ADD );  // does not overwrite add1
    // register unary operator for float
    KernelRegistry::set<UnaryOpTrait<float> >( add1<float>, ContextType::Host, KernelRegistry::KERNEL_ADD );
    KernelRegistry::set<UnaryOpTrait<float> >( minus1<float>, ContextType::Host, KernelRegistry::KERNEL_REPLACE );  // overrides add1
    KernelTraitContextFunction<UnaryOpTrait<float> > opFloat;
    KernelTraitContextFunction<UnaryOpTrait<double> > opDouble;
    double xd = opDouble[ ContextType::Host ]( 1.0 );
    float  xf = opFloat [ ContextType::Host ]( 1.0f );
    BOOST_CHECK_EQUAL( 2.0, xd );
    BOOST_CHECK_EQUAL( 0.0, xf );
}

