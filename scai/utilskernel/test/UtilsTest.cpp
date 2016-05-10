/**
 * @file UtilsTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Contains tests for the class CUDAUtils and OpenMPUtils
 * @author Jan Ecker
 * @date 19.11.2012
 */

// boost
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

// others
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/hmemo.hpp>
#include <scai/common/TypeTraits.hpp>

// import scai_arithmetic_test_types, scai_array_test_types

#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace utilskernel;
using namespace hmemo;

/* --------------------------------------------------------------------- */

extern ContextPtr testContext;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HArrayTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HArrayTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleTest, ValueType, scai_array_test_types )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    ContextPtr loc = testContext;

    setVal.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "scaleTest<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    ValueType valuesValues[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
    ValueType expectedValues[] =
    { 0, 2, 4, 6, 8, 0, 2, 4, 6, 8, 0, 2, 4, 6, 8 };
    const ValueType mult = 2;
    LArray<ValueType> values( nValues, valuesValues );
    {
        WriteAccess<ValueType> wValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        setVal[loc]( wValues.get(), nValues, mult, common::reduction::MULT );
    }
    ReadAccess<ValueType> rValues( values );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expectedValues[i], rValues[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( sumTest, ValueType, scai_array_test_types )
{
    static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;

    ContextPtr loc = testContext;

    reduce.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // print warning if not available for test context

    SCAI_LOG_INFO( logger, "sumTest<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    {
        ValueType valuesValues[] =
        { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
        const ValueType expectedSum = 30;
        LArray<ValueType> values( nValues, valuesValues );
        ReadAccess<ValueType> rValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        const ValueType resultSum = reduce[loc]( rValues.get(), nValues, common::reduction::ADD );
        BOOST_CHECK_EQUAL( expectedSum, resultSum );
    }
    {
        const ValueType expectedSum = 0;
        LArray<ValueType> values;
        ReadAccess<ValueType> rValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        const ValueType resultSum = reduce[loc]( rValues.get(), values.size(), common::reduction::ADD );
        BOOST_CHECK_EQUAL( expectedSum, resultSum );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setValTest, ValueType, scai_array_test_types )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    ContextPtr loc = testContext;

    setVal.getSupportedContext( loc );

    SCAI_LOG_INFO( logger, "setValTest<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // print warning if not available for test context

    {
        const IndexType n = 20;
        LArray<ValueType> values;
        {
            WriteOnlyAccess<ValueType> wValues( values, loc, 3 * n );
            SCAI_CONTEXT_ACCESS( loc );
            setVal[loc]( wValues.get(), 3 * n, 0, common::reduction::COPY );
            // overwrite in the middle to check that there is no out-of-range set
            setVal[loc]( wValues.get() + n, n, 10, common::reduction::COPY );
        }
        ReadAccess<ValueType> rValues( values );

        for ( IndexType i = 0; i < n; i++ )
        {
            BOOST_CHECK_EQUAL( 0, rValues.get()[i + 0 * n] );
            BOOST_CHECK_EQUAL( 10, rValues.get()[i + 1 * n] );
            BOOST_CHECK_EQUAL( 0, rValues.get()[i + 2 * n] );
        }
    }
    {
        const IndexType n = 0;
        LArray<ValueType> values;
        {
            WriteOnlyAccess<ValueType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            setVal[loc]( wValues.get(), n, 7, common::reduction::COPY );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( isSortedTest, ValueType, scai_array_test_types )
{
    static LAMAKernel<UtilKernelTrait::isSorted<ValueType> > isSorted;

    ContextPtr loc = testContext;

    isSorted.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // print warning if not available for test context

    SCAI_LOG_INFO( logger, "isSortedTest<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )
 
    {
        ValueType values1[] =
        { 1, 2, 2, 2, 5, 8 };
        ValueType values2[] =
        { 2, 2, 1, 0 };
        ValueType values3[] =
        { 1, 0, 1 };
        const IndexType nValues1 = sizeof( values1 ) / sizeof( ValueType );
        const IndexType nValues2 = sizeof( values2 ) / sizeof( ValueType );
        const IndexType nValues3 = sizeof( values3 ) / sizeof( ValueType );
        LArray<ValueType> valueArray1( nValues1, values1 );
        LArray<ValueType> valueArray2( nValues2, values2 );
        LArray<ValueType> valueArray3( nValues3, values3 );
        ReadAccess<ValueType> rValues1( valueArray1, loc );
        ReadAccess<ValueType> rValues2( valueArray2, loc );
        ReadAccess<ValueType> rValues3( valueArray3, loc );
        SCAI_CONTEXT_ACCESS( loc );
        // values1 are sorted, ascending = true
        BOOST_CHECK( isSorted[loc]( rValues1.get(), nValues1, true ) );
        BOOST_CHECK( ! isSorted[loc]( rValues1.get(), nValues1, false ) );
        // values2 are sorted, ascending = false
        BOOST_CHECK( isSorted[loc]( rValues2.get(), nValues2, false ) );
        BOOST_CHECK( ! isSorted[loc]( rValues2.get(), nValues2, true ) );
        BOOST_CHECK( isSorted[loc]( rValues2.get(), 0, true ) );
        // only first two values are sorted
        BOOST_CHECK( isSorted[loc]( rValues2.get(), 1, true ) );
        // only first two values are sorted
        BOOST_CHECK( isSorted[loc]( rValues2.get(), 2, true ) );
        // only first two values are sorted
        // values3 are not sorted, neither ascending nor descending
        BOOST_CHECK( ! isSorted[loc]( rValues3.get(), nValues3, false ) );
        BOOST_CHECK( ! isSorted[loc]( rValues3.get(), nValues3, true ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setOrderTest )
{
    // setOrder is only for IndexType

    ContextPtr loc = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "setOrderTest on " << *loc )

    static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
    {
        const IndexType n = 20;
        LArray<IndexType> values;
        {
            WriteOnlyAccess<IndexType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            setOrder[loc]( wValues.get(), n );
        }
        ReadAccess<IndexType> rValues( values );

        for ( IndexType i = 0; i < n; i++ )
        {
            BOOST_CHECK_EQUAL( i, rValues.get()[i] );
        }
    }
    {
        const IndexType n = 0;
        LArray<IndexType> values;
        {
            WriteOnlyAccess<IndexType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            setOrder[loc]( wValues.get(), n );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( invertTest, ValueType, scai_arithmetic_test_types )
{
    static LAMAKernel<UtilKernelTrait::invert<ValueType> > invert;

    ContextPtr loc = testContext;

    invert.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // print warning if not available for test context

    SCAI_LOG_INFO( logger, "invertTest<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )

    {
        // TODO: should it be possible to pass 0 elements? What should be the result?
        ValueType valuesValues[] =
        { 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
        LArray<ValueType> values( nValues, valuesValues );
        {
            WriteAccess<ValueType> wValues( values, loc );
            SCAI_CONTEXT_ACCESS( loc );
            invert[loc]( wValues.get(), nValues );
        }
        ReadAccess<ValueType> rValues( values );

        for ( IndexType i = 0; i < nValues; i++ )
        {
            BOOST_CHECK_EQUAL( ValueType( 1 ) / valuesValues[i], rValues[i] );
        }
    }
    {
        const IndexType n = 0;
        LArray<ValueType> values;
        {
            WriteOnlyAccess<ValueType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            invert[loc]( wValues.get(), n );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()

