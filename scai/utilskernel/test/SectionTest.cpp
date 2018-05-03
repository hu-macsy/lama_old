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
#include <scai/utilskernel.hpp>
#include <scai/common/TypeTraits.hpp>

// import scai_numeric_test_types, scai_array_test_types

#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace utilskernel;
using namespace hmemo;

using common::BinaryOp;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SectionTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SectonTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assign0Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assign<ValueType> > assign;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    assign.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assign0Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    ValueType rawValues[] = { 2 };
    ValueType expValues[] = { 4 };

    const IndexType sizes[] = { };
    const IndexType sourceOffset      = 0;
    const IndexType sourceDistances[] = { };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { };

    const IndexType nValues = 1;

    HArray<ValueType> source( nValues, rawValues );
    HArray<ValueType> target( nValues, rawValues );

    IndexType nDims = 0;

    { 
        WriteAccess<ValueType> wTarget( target, loc );
        ReadAccess<ValueType> rSource( source, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assign[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                     rSource.get() + sourceOffset, sourceDistances, BinaryOp::ADD, false );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assign1Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assign<ValueType> > assign;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    assign.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assign1Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 8;

    // work on this array:   0   1   2   3   4  5  6  7 

    ValueType rawValues[] = { 0, 1, 2, 3, 4, 5, 6, 7 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1, nValues );

    // array[0:3:2] = array[ 6:8 ]
    //   gives this array:   6  1   7   3    4   5  6  7

    ValueType expValues[] = { 6, 1, 7, 3,  4, 5,  6,  7 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1, nValues1 );

    const IndexType sizes[] = { 2 };
    const IndexType sourceOffset      = 6;
    const IndexType sourceDistances[] = { 1 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 };

    HArray<ValueType> source( nValues, rawValues );
    HArray<ValueType> target( nValues, rawValues );

    IndexType nDims = 1;

    {
        WriteAccess<ValueType> wTarget( target, loc );
        ReadAccess<ValueType> rSource( source, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assign[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                     rSource.get() + sourceOffset, sourceDistances, BinaryOp::COPY, false );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assign2Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assign<ValueType> > assign;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

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

    HArray<ValueType> source( nValues, rawValues );
    HArray<ValueType> target( nValues, rawValues );

    {
        WriteAccess<ValueType> wTarget( target, loc );
        ReadAccess<ValueType> rSource( source, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assign[loc]( wTarget.get() + targetOffset, 2, sizes, targetDistances, 
                     rSource.get() + sourceOffset, sourceDistances, BinaryOp::COPY, false );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assign3Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assign<ValueType> > assign;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    assign.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assign3Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    // take same arrays as in assign2Test but add a new 2nd dimension

    const IndexType n1 = 3;
    const IndexType n2 = 1;
    const IndexType n3 = 4;

    ValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3, nValues );

    ValueType expValues[] = { 12, 13, 2, 3, 10, 11, 12, 13, 22, 23, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3, nValues1 );

    const IndexType sizes[] = { 2, 1, 2 };
    const IndexType sourceOffset      = n2 * n3 + 2;
    const IndexType sourceDistances[] = { n3 * n2, n2, 1 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n3 * n2, 1, 1 };

    HArray<ValueType> source( nValues, rawValues );
    HArray<ValueType> target( nValues, rawValues );

    const IndexType nDims = 3;

    {
        WriteAccess<ValueType> wTarget( target, loc );
        ReadAccess<ValueType> rSource( source, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assign[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                     rSource.get() + sourceOffset, sourceDistances, BinaryOp::COPY, false );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assign4Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assign<ValueType> > assign;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    assign.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assign4Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    // take same arrays as in assign3Test but add a new 4th dimension

    const IndexType n1 = 3;
    const IndexType n2 = 1;
    const IndexType n3 = 4;
    const IndexType n4 = 1;

    ValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3 * n4, nValues );

    ValueType expValues[] = { 12, 13, 2, 3, 10, 11, 12, 13, 22, 23, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3 * n4, nValues1 );

    const IndexType sizes[] = { 2, 1, 2, 1 };
    const IndexType sourceDistances[] = { n2 * n3 * n4, n3 * n4 , n4, 1 };
    const IndexType sourceOffset      = 1 * sourceDistances[0] + 2 * sourceDistances[3];
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2 * n3 * n4, n3 * n4, n4, 1 };

    HArray<ValueType> source( nValues, rawValues );
    HArray<ValueType> target( nValues, rawValues );

    const IndexType nDims = 4;

    {
        WriteAccess<ValueType> wTarget( target, loc );
        ReadAccess<ValueType> rSource( source, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assign[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                     rSource.get() + sourceOffset, sourceDistances, BinaryOp::COPY, false );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignScalar0Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assignScalar<ValueType> > assignScalar;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    assignScalar.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assignScalar0Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    ValueType rawValues[] = { 2 };
    ValueType expValues[] = { 5 };   // scalar - 2

    const IndexType sizes[] = { };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { };

    const IndexType nValues = 1;

    HArray<ValueType> target( nValues, rawValues );

    IndexType nDims = 0;
    ValueType scalar = 7;

    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assignScalar[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                           scalar, BinaryOp::SUB, true );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignScalar1Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assignScalar<ValueType> > assignScalar;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    assignScalar.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assignScalar1Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 8;

    // work on this array:   0   1   2   3   4  5  6  7 

    ValueType rawValues[] = { 0, 1, 2, 3, 4, 5, 6, 7 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1, nValues );

    // array[0:3:2] = 10 - array[0:3:2] 
    //   gives this array:   10  1   8   3    4   5  6  7

    ValueType expValues[] = { 10, 1, 8, 3,  4, 5,  6,  7 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1, nValues1 );

    const IndexType sizes[] = { 2 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 };

    HArray<ValueType> target( nValues, rawValues );

    IndexType nDims = 1;
    ValueType scalar = 10;

    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assignScalar[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                           scalar, BinaryOp::SUB, true );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignScalar2Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assignScalar<ValueType> > assignScalar;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    assignScalar.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assignScalar2Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 4;

    // work on this array:   0   1   2   3
    //                      10  11  12  13
    //                      20  21  22  23

    ValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2, nValues );

    // array[0:3:2,0:2] *= 2
    //   gives this array:   0   2   2   3
    //                      10  11  12  13
    //                      40  42  22  23

    ValueType expValues[] = { 0, 2, 2, 3, 10, 11, 12, 13, 40, 42, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2, nValues1 );

    const IndexType sizes[] = { 2, 2 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2, 1 };

    HArray<ValueType> target( nValues, rawValues );

    const IndexType nDims = 2;
    ValueType scalarVal = 2;
    bool swap = false;

    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assignScalar[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                           scalarVal, BinaryOp::MULT, swap );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignScalar3Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assignScalar<ValueType> > assignScalar;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    assignScalar.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assignScalar3Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 1;
    const IndexType n3 = 4;

    ValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3, nValues );

    ValueType expValues[] = { 0, 2, 2, 3, 10, 11, 12, 13, 40, 42, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3, nValues1 );

    const IndexType sizes[] = { 2, 1, 2 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2 * n3, n3, 1 };

    HArray<ValueType> target( nValues, rawValues );

    const IndexType nDims = 3;
    ValueType scalarVal = 2;
    bool swap = false;

    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assignScalar[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                           scalarVal, BinaryOp::MULT, swap );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignScalar4Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<SectionKernelTrait::assignScalar<ValueType> > assignScalar;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    assignScalar.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "assignScalar4Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 1;
    const IndexType n3 = 4;
    const IndexType n4 = 1;

    ValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3 * n4, nValues );

    ValueType expValues[] = { 0, 2, 2, 3, 10, 11, 12, 13, 40, 42, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3 * n4, nValues1 );

    const IndexType nDims = 4;

    const IndexType sizes[] = { 2, 1, 2, 1 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2 * n3 * n4, n3 * n4, n4, 1 };

    const IndexType nSizes = sizeof( sizes ) / sizeof( IndexType );
    const IndexType nDistances = sizeof( targetDistances ) / sizeof( IndexType );

    SCAI_ASSERT_EQ_ERROR( nDims, nSizes, "serious mismatch of sizes array" )
    SCAI_ASSERT_EQ_ERROR( nDims, nDistances, "serious mismatch of distances array" )
 
    HArray<ValueType> target( nValues, rawValues );

    ValueType scalarVal = 2;
    bool swap = false;

    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc );
        assignScalar[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                           scalarVal, BinaryOp::MULT, swap );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( unaryOp1Test, SourceValueType, scai_array_test_types )
{
    typedef DefaultReal TargetValueType;

    static LAMAKernel<SectionKernelTrait::unaryOp<TargetValueType, SourceValueType> > unaryOp;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    unaryOp.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "unaryOp1Test<" << common::TypeTraits<TargetValueType>::id() << ", " 
                                           << common::TypeTraits<SourceValueType>::id() << "> for " 
                                           << *testContext << ", done on " << *loc )

    // work on this array:   0   1   2   3   4  5  6  7 
    // array[0:3:2] = array[ 6:8 ]
    //   gives this array:   6  1   7   3    4   5  6  7

    HArray<SourceValueType> source(    { 0, 1, 2, 3, 4, 5, 6, 7 } );
    HArray<TargetValueType> expTarget( { 6, 1, 7, 3, 4, 5, 6, 7 } );

    BOOST_REQUIRE_EQUAL( source.size(), expTarget.size() );

    const IndexType sizes[] = { 2 };
    const IndexType sourceOffset      = 6;
    const IndexType sourceDistances[] = { 1 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 };

    auto target = utilskernel::convertHArray<TargetValueType>( source );

    IndexType nDims = 1;

    {
        WriteAccess<TargetValueType> wTarget( target, loc );
        ReadAccess<SourceValueType> rSource( source, loc );
        SCAI_CONTEXT_ACCESS( loc );
        unaryOp[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                      rSource.get() + sourceOffset, sourceDistances, common::UnaryOp::COPY );
    }

    BOOST_TEST( hostReadAccess( target ) == hostReadAccess( expTarget ), per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( unaryOp2Test, SourceValueType, scai_array_test_types )
{
    typedef DefaultReal TargetValueType;

    static LAMAKernel<SectionKernelTrait::unaryOp<TargetValueType, SourceValueType> > unaryOp;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    unaryOp.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "unaryOp2Test<" << common::TypeTraits<TargetValueType>::id() << ", " 
                                           << common::TypeTraits<SourceValueType>::id() << "> for " 
                                           << *testContext << ", done on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 4;

    // work on this array:   0   1   2   3
    //                      10  11  12  13
    //                      20  21  22  23

    SourceValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( SourceValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2, nValues );

    // array[0:3:2,0:2] = array[1:3,2:4]
    //   gives this array:  12  13   2   3
    //                      10  11  12  13
    //                      22  23  22  23

    TargetValueType expValues[] = { 12, 13, 2, 3, 10, 11, 12, 13, 22, 23, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( TargetValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2, nValues1 );

    const IndexType sizes[] = { 2, 2 };
    const IndexType sourceOffset      = n2 + 2;
    const IndexType sourceDistances[] = { n2, 1 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2, 1 };

    HArray<SourceValueType> source( nValues, rawValues );
    auto target = utilskernel::convertHArray<TargetValueType>( source );

    {
        WriteAccess<TargetValueType> wTarget( target, loc );
        ReadAccess<SourceValueType> rSource( source, loc );
        SCAI_CONTEXT_ACCESS( loc );
        unaryOp[loc]( wTarget.get() + targetOffset, 2, sizes, targetDistances, 
                      rSource.get() + sourceOffset, sourceDistances, common::UnaryOp::COPY );
    }

    ReadAccess<TargetValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( unaryOp3Test, SourceValueType, scai_array_test_types )
{
    typedef DefaultReal TargetValueType;

    static LAMAKernel<SectionKernelTrait::unaryOp<TargetValueType, SourceValueType> > unaryOp;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    unaryOp.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "unaryOp3Test<" << common::TypeTraits<TargetValueType>::id() << ", " 
                                           << common::TypeTraits<SourceValueType>::id() << "> for " 
                                           << *testContext << ", done on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 1;
    const IndexType n3 = 4;

    SourceValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( SourceValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3, nValues );

    // array[0:3:2,0:2] = array[1:3,2:4]
    //   gives this array:  12  13   2   3
    //                      10  11  12  13
    //                      22  23  22  23

    TargetValueType expValues[] = { 12, 13, 2, 3, 10, 11, 12, 13, 22, 23, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( TargetValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3, nValues1 );

    const IndexType sizes[] = { 2, 1, 2 };
    const IndexType sourceOffset      = n3 + 2;
    const IndexType sourceDistances[] = { n2 * n3, n3, 1 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2 * n3, n3, 1 };

    HArray<SourceValueType> source( nValues, rawValues );
    auto target = utilskernel::convertHArray<TargetValueType>( source );

    IndexType nDims = 3;

    {
        WriteAccess<TargetValueType> wTarget( target, loc );
        ReadAccess<SourceValueType> rSource( source, loc );
        SCAI_CONTEXT_ACCESS( loc );
        unaryOp[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                      rSource.get() + sourceOffset, sourceDistances, common::UnaryOp::COPY );
    }

    ReadAccess<TargetValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( unaryOp4Test, SourceValueType, scai_array_test_types )
{
    typedef DefaultReal TargetValueType;

    static LAMAKernel<SectionKernelTrait::unaryOp<TargetValueType, SourceValueType> > unaryOp;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    unaryOp.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "unaryOp4Test<" << common::TypeTraits<TargetValueType>::id() << ", " 
                                           << common::TypeTraits<SourceValueType>::id() << "> for " 
                                           << *testContext << ", done on " << *loc )

    IndexType nDims = 4;

    const IndexType n1 = 3;
    const IndexType n2 = 1;
    const IndexType n3 = 4;
    const IndexType n4 = 1;

    SourceValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( SourceValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3 * n4, nValues );

    // array[0:3:2,0:2] = array[1:3,2:4]
    //   gives this array:  12  13   2   3
    //                      10  11  12  13
    //                      22  23  22  23

    TargetValueType expValues[] = { 12, 13, 2, 3, 10, 11, 12, 13, 22, 23, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( TargetValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3 * n4, nValues1 );

    const IndexType sizes[] = { 2, 1, 2, 1 };
    const IndexType sourceOffset      = n3 * n4 + 2;
    const IndexType sourceDistances[] = { n2 * n3 * n4, n3 * n4,  n4, 1 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2 * n3 * n4, n3 * n4, n4, 1 };

    HArray<SourceValueType> source( nValues, rawValues );
    auto target = utilskernel::convertHArray<TargetValueType>( source );

    {
        WriteAccess<TargetValueType> wTarget( target, loc );
        ReadAccess<SourceValueType> rSource( source, loc );
        SCAI_CONTEXT_ACCESS( loc );
        unaryOp[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, 
                      rSource.get() + sourceOffset, sourceDistances, common::UnaryOp::COPY );
    }

    ReadAccess<TargetValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( UnaryOp1Test, ValueType, scai_numeric_test_types )
{
    static LAMAKernel<SectionKernelTrait::UnaryOp<ValueType> > UnaryOp;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    UnaryOp.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "UnaryOp1Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 8;

    ValueType rawValues[] = { 0, 1, 2, 3, 4, 5, 6, 7 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1, nValues );

    // array[1:4:2] = -array[1:4:2]

    ValueType expValues[] = { 0, -1, 2, -3,  4, 5,  6,  7 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1, nValues1 );

    const IndexType sizes[] = { 2 };
    const IndexType targetOffset      = 1;
    const IndexType targetDistances[] = { 2 };

    HArray<ValueType> target( nValues, rawValues );

    IndexType nDims = 1;

    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc );
        UnaryOp[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, common::UnaryOp::MINUS );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( UnaryOp2Test, ValueType, scai_numeric_test_types )
{
    static LAMAKernel<SectionKernelTrait::UnaryOp<ValueType> > UnaryOp;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    UnaryOp.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "UnaryOp2Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 4;

    // work on this array:   0   1   2   3
    //                      10  11  12  13
    //                      20  21  22  23

    ValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2, nValues );

    // array[0:3:2,0:2] -= array[0:3:2,0:2] 
    //   gives this array:   0  -1   2   3
    //                      10  11  12  13
    //                     -20 -21  22  23

    ValueType expValues[] = { 0, -1, 2, 3, 10, 11, 12, 13, -20, -21, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2, nValues1 );

    const IndexType sizes[] = { 2, 2 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2, 1 };

    HArray<ValueType> target( nValues, rawValues );
    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc );
        UnaryOp[loc]( wTarget.get() + targetOffset, 2, sizes, targetDistances, common::UnaryOp::MINUS );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( UnaryOp3Test, ValueType, scai_numeric_test_types )
{
    static LAMAKernel<SectionKernelTrait::UnaryOp<ValueType> > UnaryOp;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    UnaryOp.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "UnaryOp3Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 1;
    const IndexType n3 = 4;

    ValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3, nValues );

    ValueType expValues[] = { 0, -1, 2, 3, 10, 11, 12, 13, -20, -21, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3, nValues1 );

    const IndexType sizes[] = { 2, 1, 2 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2 * n3, n2, 1 };

    const IndexType nDims = 3;

    HArray<ValueType> target( nValues, rawValues );
    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc );
        UnaryOp[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, common::UnaryOp::MINUS );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( UnaryOp4Test, ValueType, scai_numeric_test_types )
{
    static LAMAKernel<SectionKernelTrait::UnaryOp<ValueType> > UnaryOp;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    UnaryOp.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "UnaryOp4Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    const IndexType n1 = 3;
    const IndexType n2 = 1;
    const IndexType n3 = 4;
    const IndexType n4 = 1;

    ValueType rawValues[] = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };

    const IndexType nValues = sizeof( rawValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3 * n4, nValues );

    ValueType expValues[] = { 0, -1, 2, 3, 10, 11, 12, 13, -20, -21, 22, 23 };
    const IndexType nValues1 = sizeof( expValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1 * n2 * n3 * n4, nValues1 );

    const IndexType sizes[] = { 2, 1, 2, 1 };
    const IndexType targetOffset      = 0;
    const IndexType targetDistances[] = { 2 * n2 * n3 * n4, n3 * n4, n4, 1 };

    const IndexType nDims = 4;

    HArray<ValueType> target( nValues, rawValues );
    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc );
        UnaryOp[loc]( wTarget.get() + targetOffset, nDims, sizes, targetDistances, common::UnaryOp::MINUS );
    }

    ReadAccess<ValueType> rTarget( target );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expValues[i], rTarget[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */


BOOST_AUTO_TEST_SUITE_END()

