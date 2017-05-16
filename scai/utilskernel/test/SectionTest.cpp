/**
 * @file SectionTest.cpp
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
 * @brief Contains tests for all kernels working on sections     
 * @author Thomas Brandes
 * @date 16.05.2017
 */

// boost
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

// others
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/SectionKernelTrait.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/common/TypeTraits.hpp>

// import scai_numeric_test_types, scai_array_test_types

#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace utilskernel;
using namespace hmemo;
using common::binary;

/* --------------------------------------------------------------------- */

extern ContextPtr testContext;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SectionTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SectonTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assign2Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assign<ValueType> > assign;

    ContextPtr loc = testContext;

    assign.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assign2Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 4;

    // work on this array:   0   1   2   3
    //                      10  11  12  13
    //                      20  21  22  23

    ValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2, nValues );

    // array[0:3:2,0:2] = array[1:3,2:4]
    //   gives this array:  12  13   2   3
    //                      10  11  12  13
    //                      22  23  22  23

    ValueType expValues[] = { 12, 13, 2, 3, 10, 11, 12, 13, 22, 23, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2, nValues1 );

    const IndexType sizes[] = { 2, 2 };
    const IndexType sourceOffset      = n2 + 2;
    const IndexType sourceDistances[] = { n2, 1 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2, 1 };

    LArray<ValueType> source( nValues, rawValues );
    LArray<ValueType> target( nValues, rawValues );

    {
        WriteAccess<ValueType> wTarget( target, loc );
        ReadAccess<ValueType> rSource( source, loc );
        assign[loc]( wTarget.get() + targetOffset, 2, sizes, targetDistances, 
                     rSource.get() + sourceOffset, sourceDistances, binary::COPY, false );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */


BOOST_AUTO_TEST_SUITE_END()

