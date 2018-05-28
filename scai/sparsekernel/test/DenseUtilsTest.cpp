/**
 * @file DenseUtilsTest.cpp
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
 * @brief Contains tests for the class DenseUtils that run on different devices
 * @author Thomas Brandes
 * @date 19.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/sparsekernel/DenseKernelTrait.hpp>
#include <scai/sparsekernel/DenseUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;

using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;

using common::TypeTraits;
using common::BinaryOp;
using common::UnaryOp;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DenseUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( nonZeroValuesTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    static LAMAKernel<DenseKernelTrait::nonZeroValues<ValueType> > nonZeroValues;

    ContextPtr loc = testContext;
    nonZeroValues.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "nonZeroValues test for " << *testContext << " on " << *loc )

    const ValueType dense_values[] = { 1, 0, 2, -3, 1, 0, 2, 4, -5, 6, 3, 1 };

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    const IndexType dense_n = sizeof( dense_values ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( dense_n, numRows * numColumns );

    HArray<ValueType> dense( numRows * numColumns, dense_values, testContext );

    SCAI_CONTEXT_ACCESS( loc );

    ReadAccess<ValueType> rDense( dense, loc );

    IndexType count = nonZeroValues[loc]( rDense.get(), numRows, numColumns, 0 );

    // two values were zero

    BOOST_CHECK_EQUAL( count, dense_n - 2 );

    count = nonZeroValues[loc]( rDense.get(), numRows, numColumns, 10 );

    // all values are less than 10

    BOOST_CHECK_EQUAL( count, IndexType( 0 ) );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getCSRTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<ValueType> dense( { 1, 0, 2,
                              -3, 1, 0,
                               0, 0, -5,
                               6, 0, 1   }, testContext );

    // Note: compress keeps order of values

    HArray<IndexType> expIA( { 2, 2, 1, 2 } );
    HArray<IndexType> expJA( { 0, 2, 0, 1, 2, 0, 2 } );
    HArray<ValueType> expValues( { 1, 2, -3, 1, -5, 6, 1 } );

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    BOOST_REQUIRE_EQUAL( dense.size(), numRows * numColumns );

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    DenseUtils::getSparseRowSizes( csrIA, numRows, numColumns, dense, testContext );
    BOOST_TEST( hostReadAccess( expIA ) == hostReadAccess( csrIA ), per_element() );

    DenseUtils::convertDense2CSR( csrIA, csrJA, csrValues, numRows, numColumns, dense, testContext );

    BOOST_TEST( hostReadAccess( expJA ) == hostReadAccess( csrJA ), per_element() );
    BOOST_TEST( hostReadAccess( expValues ) == hostReadAccess( csrValues ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setCSRTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "setCSRValues test for " << *testContext )

    HArray<IndexType> csrIA( { 0, 2, 4, 5, 7 }, testContext );
    HArray<IndexType> csrJA( { 0, 2, 0, 1, 2, 0, 2 }, testContext );
    HArray<ValueType> csrValues( { 1, 2, -3, 1, -5, 6, 1 }, testContext );

    // expected dense storage data, stored row-wise

    HArray<ValueType> expDense( { 1, 0, 2,
                                 -3, 1, 0,
                                  0, 0, -5,
                                  6, 0, 1   } );

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    BOOST_REQUIRE_EQUAL( csrValues.size(), csrJA.size() );
    BOOST_REQUIRE_EQUAL( csrIA.size(), numRows + 1 );
    BOOST_REQUIRE_EQUAL( expDense.size(), numRows * numColumns );

    HArray<ValueType> dense;

    DenseUtils::convertCSR2Dense( dense, numRows, numColumns, csrIA, csrJA, csrValues, testContext );

    BOOST_TEST( hostReadAccess( expDense ) == hostReadAccess( dense ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DenseKernelTrait::set<ValueType, DefaultReal> > set;

    ContextPtr loc = testContext;

    set.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "set test for " << *testContext << " on " << *loc )

    const ValueType dense_values[] = { 1, 1, 2,
                                       3, 1, 3,
                                       2, 4, 5,
                                       6, 9, 1
                                     };

    const DefaultReal dense_values1[] = { 1, 1, 2,
                                       2, 1, 1,
                                       1, 1, 1,
                                       2, 2, 2
                                     };

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    const IndexType dense_n = sizeof( dense_values ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( dense_n, numRows * numColumns );

    for ( int op = 0; op < static_cast<int>( BinaryOp::MAX_BINARY_OP ); ++op )
    {
        BinaryOp binOp = BinaryOp( op );

        HArray<ValueType> dense( numRows * numColumns, dense_values, testContext );
        HArray<DefaultReal> dense1( numRows * numColumns, dense_values1, testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );

            WriteAccess<ValueType> wDense( dense, loc );
            ReadAccess<DefaultReal> rDense1( dense1, loc );

            set[loc]( wDense.get(), numRows, numColumns, rDense1.get(), binOp );
        }

        {
            ReadAccess<ValueType> rDense( dense, hostContext );

            for ( IndexType i = 0; i < dense_n; ++i )
            {
                ValueType val2 = static_cast<ValueType>( dense_values1[i] );
                ValueType expectedVal = applyBinary( dense_values[i], binOp, val2 );
                BOOST_CHECK_EQUAL( expectedVal, rDense[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DenseKernelTrait::setValue<ValueType> > setValue;

    ContextPtr loc = testContext;

    setValue.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "setValue test for " << *testContext << " on " << *loc )

    const ValueType dense_values[] = { 1, 1, 2,
                                       3, 1, 3,
                                       2, 4, 5,
                                       6, 9, 1
                                     };

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    const IndexType dense_n = sizeof( dense_values ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( dense_n, numRows * numColumns );

    ValueType val = 2;

    for ( int op = 0; op < static_cast<int>( BinaryOp::MAX_BINARY_OP ); ++op )
    {
        BinaryOp binOp = BinaryOp( op );

        HArray<ValueType> dense( numRows * numColumns, dense_values, testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );

            WriteAccess<ValueType> wDense( dense, loc );

            setValue[loc]( wDense.get(), numRows, numColumns, val, binOp );
        }

        {
            ReadAccess<ValueType> rDense( dense, hostContext );

            for ( IndexType i = 0; i < dense_n; ++i )
            {
                ValueType expectedVal = applyBinary( dense_values[i], binOp, val );
                BOOST_CHECK_EQUAL( expectedVal, rDense[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleRowsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DenseKernelTrait::scaleRows<ValueType> > scaleRows;

    ContextPtr loc = testContext;

    scaleRows.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "scaleRows test for " << *testContext << " on " << *loc )

    const ValueType dense_values[] = { 1, 1, 2,
                                       3, 1, 3,
                                       2, 4, 5,
                                       6, 9, 1
                                     };

    const ValueType row_values[] = { 1, 2, 0, 1 };

    const ValueType dense1_values[] = { 1, 1, 2,
                                        6, 2, 6,
                                        0, 0, 0,
                                        6, 9, 1
                                      };

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    const IndexType dense_n = sizeof( dense_values ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( dense_n, numRows * numColumns );

    HArray<ValueType> rows( numRows, row_values, testContext );
    HArray<ValueType> dense( numRows * numColumns, dense_values, testContext );

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<ValueType> wDense( dense, loc );
        ReadAccess<ValueType> rRows( rows, loc );

        scaleRows[loc]( wDense.get(), numRows, numColumns, rRows.get() );
    }

    {
        ReadAccess<ValueType> rDense( dense, hostContext );

        for ( IndexType i = 0; i < dense_n; ++i )
        {
            BOOST_CHECK_EQUAL( dense1_values[i], rDense[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()

