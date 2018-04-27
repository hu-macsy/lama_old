/**
 * @file ELLUtilsTest.cpp
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
 * @brief Contains tests for kernel implementations of ELLKernelTrait routines.
 * @author Thomas Brandes
 * @date 15.10.2012
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/sparsekernel/ELLKernelTrait.hpp>
#include <scai/sparsekernel/openmp/OpenMPELLUtils.hpp>
#include <scai/utilskernel.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/kregistry.hpp>
#include <scai/hmemo.hpp>
#include <scai/common/Math.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>
#include <scai/sparsekernel/test/TestData1.hpp>
#include <scai/sparsekernel/test/TestData2.hpp>

#include <scai/hmemo/test/ContextFix.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;
using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;
using namespace kregistry;
using common::TypeTraits;
using common::Exception;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ELLUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.ELLUtilsTest" )

/* ------------------------------------------------------------------------------------- */

// BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxDiffValTest, ValueType, scai_numeric_test_types )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( hasDiagonalPropertyTest )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::hasDiagonalProperty> hasDiagonalProperty;

    ContextPtr loc = testContext;
    hasDiagonalProperty.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data2::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    // positive test

    {
        const IndexType numDiagonals = common::Math::min( numRows, numColumns );

        ReadAccess<IndexType> rEllJa( ellJA, loc );
        SCAI_CONTEXT_ACCESS( loc );
        bool diagonalProperty = hasDiagonalProperty[loc]( numDiagonals, rEllJa.get() );
        BOOST_CHECK( diagonalProperty );
    }

    // negative test

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    {
        const IndexType numDiagonals = common::Math::min( numRows, numColumns );

        ReadAccess<IndexType> rEllJa( ellJA, loc );
        SCAI_CONTEXT_ACCESS( loc );

        bool diagonalProperty = hasDiagonalProperty[loc]( numDiagonals, rEllJa.get() );

        BOOST_CHECK( !diagonalProperty );
    }

    // test empty array

    ellJA.clear();

    {
        const IndexType numDiagonals = 0;
        HArray<IndexType> ellJa;
        ReadAccess<IndexType> rEllJa( ellJa, loc );
        SCAI_CONTEXT_ACCESS( loc );
        bool diagonalProperty = hasDiagonalProperty[loc]( numDiagonals, rEllJa.get() );

        BOOST_CHECK( !diagonalProperty );
    }
}
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( fillELlValuesTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<ELLKernelTrait::fillELLValues<ValueType> > fillELLValues;

    ContextPtr loc = testContext;
    fillELLValues.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    {
        WriteAccess<IndexType> wJA( ellJA, loc );
        WriteAccess<ValueType> wValues( ellValues, loc );
        ReadAccess<IndexType> rIA( ellIA, loc );

        SCAI_CONTEXT_ACCESS( loc );

        fillELLValues[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, numValuesPerRow );
    }

    HArray<ValueType> xDummy( numColumns, ValueType( 1 ) );

    {
        ReadAccess<IndexType> rIA( ellIA, hostContext );
        ReadAccess<IndexType> rJA( ellJA, hostContext );
        ReadAccess<ValueType> rValues( ellValues, hostContext );
        ReadAccess<ValueType> rX( xDummy, hostContext );

        ValueType testVal = 0;

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType jj = rIA[i]; jj < numValuesPerRow; ++jj )
            {
                IndexType pos = jj * numRows + i;
                IndexType j   = rJA[ pos ];
                SCAI_ASSERT_VALID_INDEX( j, numColumns, "illegal col pos" )
                testVal += rValues[ pos ] * rX[ j ];
            }
        }

        BOOST_CHECK_EQUAL( ValueType( 0 ), testVal );
    }
}


/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( checkTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::check> check;

    ContextPtr loc = testContext;
    check.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    // check with correct values
    {
        const IndexType valuesIa[] = { 4, 3, 5, 2 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof( IndexType );
        const IndexType minusOne = static_cast<IndexType>( -1 );
        const IndexType valuesJa[] = { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, minusOne, 4, minusOne, 4, 0, 0, 0, 5, 0 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
        const IndexType numRows = nIa;
        const IndexType numValuesPerRow = 5;
        const IndexType numColumns = 6;
        HArray<IndexType> ia( nIa, valuesIa, testContext );
        HArray<IndexType> ja( nJa, valuesJa, testContext );
        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );
        SCAI_CONTEXT_ACCESS( loc );
        BOOST_CHECK_NO_THROW( check[loc]( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ) );
    }
    // check with invalid ia
    {
        HArray<IndexType> ia( { 4, 3, 7, 2 }, testContext );   // contains at least one illegal size
        HArray<IndexType> ja( { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 0, 4, 0, 4, 0, 0, 0, 5, 0 }, testContext );
        
        const IndexType numRows = ia.size();
        const IndexType numValuesPerRow = ja.size() / numRows;
        SCAI_ASSERT_EQ_ERROR( numRows * numValuesPerRow, ja.size(), "illegal size of ja" )
        const IndexType numColumns = 5;

        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );
        SCAI_CONTEXT_ACCESS( loc );
        BOOST_CHECK_THROW( check[loc]( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ),
                           Exception );
    }
    // check with invalid ja
    {
        const IndexType valuesIa[] = { 4, 3, 5, 2 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof( IndexType );
        const IndexType valuesJa[] = { 1, 1, 1, 1, 2, 2, 2, static_cast<IndexType>( -1 ), 3, 3, 3, 0, 4, 0, 4, 0, 0, 0, 5, 0 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
        const IndexType numRows = nIa;
        const IndexType numValuesPerRow = 5;
        const IndexType numColumns = 5;
        HArray<IndexType> ia( nIa, valuesIa, testContext );
        HArray<IndexType> ja( nJa, valuesJa, testContext );
        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );
        SCAI_CONTEXT_ACCESS( loc );
        BOOST_CHECK_THROW( check[loc]( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ),
                           Exception );
    }
    // check with valid empty values
    {
        const IndexType numRows = 0;
        const IndexType numValuesPerRow = 0;
        const IndexType numColumns = 0;
        HArray<IndexType> ia;
        HArray<IndexType> ja;
        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );
        SCAI_CONTEXT_ACCESS( loc );
        BOOST_CHECK_NO_THROW( check[loc]( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getRowTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::getRow<ValueType> > getRow;

    ContextPtr loc = testContext;
    getRow.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    // check with valid dense values
    {
        ValueType valuesValues[] = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
        IndexType valuesIa[] = { 5, 5, 5 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof( IndexType );
        IndexType valuesJa[] = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
        const IndexType i = 1;
        const IndexType numRows = nIa;
        const IndexType numValuesPerRow = nJa / nIa;
        const IndexType numColumns = 5;
        HArray<ValueType> values( nValues, valuesValues, testContext );
        HArray<IndexType> ia( nIa, valuesIa, testContext );
        HArray<IndexType> ja( nJa, valuesJa, testContext );
        HArray<ValueType> row( numColumns, ValueType( 0 ) );
        {
            ReadAccess<ValueType> rValues( values, loc );
            ReadAccess<IndexType> rIa( ia, loc );
            ReadAccess<IndexType> rJa( ja, loc );
            WriteOnlyAccess<ValueType> wRow( row, loc, numColumns );
            SCAI_CONTEXT_ACCESS( loc );
            getRow[loc]( wRow.get(), i, numRows, numColumns, numValuesPerRow, rIa.get(), rJa.get(), rValues.get() );
        }

        std::vector<ValueType> expectedValues( { 0, 1, 2, 3, 4 } );
        BOOST_TEST( hostReadAccess( row ) == expectedValues, boost::test_tools::per_element() );
    }
    // check with valid sparse values
    {
        HArray<IndexType> ia( { 5, 5, 5 }, testContext );
        HArray<ValueType> values( { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 }, testContext );
        HArray<IndexType> ja( { 0, 0, 0, 2, 2, 2, 4, 4, 4, 6, 6, 6, 10, 10, 10 }, testContext );

        const IndexType numRows = ia.size();
        const IndexType numColumns = 11;
        const IndexType numValuesPerRow = ja.size() / numRows;

        BOOST_CHECK_EQUAL( numValuesPerRow * numRows, ja.size() );
        BOOST_CHECK_EQUAL( numValuesPerRow * numRows, values.size() );

        const IndexType i = 1;  // row to get

        HArray<ValueType> row( numColumns, ValueType( 0 ) );
        {
            ReadAccess<ValueType> rValues( values, loc );
            ReadAccess<IndexType> rIa( ia, loc );
            ReadAccess<IndexType> rJa( ja, loc );
            WriteOnlyAccess<ValueType> wRow( row, loc, numColumns );
            SCAI_CONTEXT_ACCESS( loc );
            getRow[loc]( wRow.get(), i, numRows, numColumns, numValuesPerRow, rIa.get(), rJa.get(), rValues.get() );
        }

        HArray<ValueType> expectedValues( { 0, 0, 1, 0, 2, 0, 3, 0, 0, 0, 4 } );

        BOOST_TEST( hostReadAccess( row ) == hostReadAccess( expectedValues ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = testContext;
    getValuePos.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    IndexType ia_values[] = { 5, 5, 5 };

    HArray<ValueType> values(  { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 }, testContext );
    HArray<IndexType> ja(      { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 }, testContext );
    HArray<IndexType> ia(  3, ia_values, testContext );

    ValueType expectedValues[] =  { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };

    const IndexType numRows = ia.size();
    const IndexType numValuesPerRow = ja.size() / numRows;

    BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, ja.size() );
    BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, values.size() );

    auto rValues = hostReadAccess( values );  // for check on host

    ReadAccess<IndexType> rIa( ia, loc );
    ReadAccess<IndexType> rJa( ja, loc );

    SCAI_CONTEXT_ACCESS( loc );

    for ( IndexType i = 0; i < numRows; i++ )
    {
        for ( IndexType j = 0; j < ia_values[i]; j++ )
        {
            IndexType pos = getValuePos[loc]( i, j, numRows, numValuesPerRow, rIa.get(), rJa.get() );
            BOOST_CHECK_EQUAL( expectedValues[j * numRows + i], rValues[ pos ] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleRowsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    LAMAKernel<ELLKernelTrait::scaleRows<ValueType> > scaleRows;

    ContextPtr loc = testContext;
    scaleRows.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    ValueType mValues[] =
    { 1, 2, 3, 4, 5, 2, 2, 2, 2, 2, 4, 2, 0, 1, 3, 0, 0, 0, 0, 3 };
    const IndexType nValues = sizeof( mValues ) / sizeof( ValueType );
    const ValueType expectedValues[] =
    { 2, 4, 15, 8, 10, 4, 4, 10, 4, 4, 8, 4, 0, 2, 6, 0, 0, 0, 0, 6 };
    HArray<ValueType> ellValues( nValues, mValues, testContext );
    {
        const IndexType numRows = 5;
        const IndexType numValuesPerRow = nValues / numRows;
        const IndexType ellIaValues[] =
        { 3, 3, 3, 3, 4 };
        const IndexType n = sizeof( ellIaValues ) / sizeof( IndexType );
        const ValueType values[] = { 2, 2, 5, 2, 2 };
        HArray<IndexType> ellIa( n, ellIaValues );
        HArray<ValueType> diagonal( n, values );
        ReadAccess<IndexType> rEllIa( ellIa, loc );
        WriteAccess<ValueType> wEllValues( ellValues, loc );
        ReadAccess<ValueType> rScaleValues( diagonal, loc );
        SCAI_CONTEXT_ACCESS( loc );
        scaleRows[loc]( wEllValues.get(), numRows, numValuesPerRow, rEllIa.get(), rScaleValues.get() );
    }
    ReadAccess<ValueType> rEllValues( ellValues );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expectedValues[i], rEllValues[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getCSRValuesTest, ValueType, scai_numeric_test_types )
{
    typedef DefaultReal OtherValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::getCSRValues<ValueType, OtherValueType> > getCSRValues;

    ContextPtr loc = testContext;
    getCSRValues.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    ValueType valuesELLValues[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
    const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof( ValueType );
    IndexType valuesELLIa[] =
    { 5, 5, 5 };
    const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof( IndexType );
    IndexType valuesELLJa[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
    const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof( IndexType );
    IndexType valuesCSRIa[] =
    { 0, 5, 10, 15 };
    const IndexType nCSRIa = sizeof( valuesCSRIa ) / sizeof( IndexType );
    OtherValueType expectedCSRValues[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    IndexType expectedCSRJa[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    const IndexType numRows = nELLIa;
    const IndexType numValuesPerRow = nELLValues / numRows;
    // make sure that division did fit
    BOOST_REQUIRE_EQUAL( numValuesPerRow * numRows, nELLValues );
    const IndexType nCSRValues = 15;
    HArray<ValueType> ellValues( nELLValues, valuesELLValues, testContext );
    HArray<IndexType> ellIa( nELLIa, valuesELLIa, testContext );
    HArray<IndexType> ellJa( nELLJa, valuesELLJa, testContext );
    HArray<OtherValueType> csrValues( nCSRValues, OtherValueType( 0 ), testContext );
    HArray<IndexType> csrIa( nCSRIa, valuesCSRIa, testContext );
    HArray<IndexType> csrJa( nCSRValues, IndexType( 0 ), testContext );
    {
        ReadAccess<ValueType> rELLValues( ellValues, loc );
        ReadAccess<IndexType> rELLIa( ellIa, loc );
        ReadAccess<IndexType> rELLJa( ellJa, loc );
        ReadAccess<IndexType> rCSRIa( csrIa, loc );
        WriteOnlyAccess<OtherValueType> wCSRValues( csrValues, loc, nCSRValues );
        WriteOnlyAccess<IndexType> wCSRJa( csrJa, loc, nCSRValues );
        SCAI_CONTEXT_ACCESS( loc );
        getCSRValues[loc] ( wCSRJa.get(), wCSRValues.get(), rCSRIa.get(), numRows, numValuesPerRow, rELLIa.get(),
                            rELLJa.get(), rELLValues.get() );
    }
    ReadAccess<IndexType> rCSRJa( csrJa );
    ReadAccess<OtherValueType> rCSRValues( csrValues );

    for ( IndexType i = 0; i < nCSRValues; i++ )
    {
        BOOST_CHECK_EQUAL( expectedCSRJa[i], rCSRJa[i] );
        BOOST_CHECK_EQUAL( expectedCSRValues[i], rCSRValues[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setCSRValuesTest, ValueType, scai_numeric_test_types )

{
    typedef DefaultReal OtherValueType;
    ContextPtr testContext = Context::getContextPtr();
    LAMAKernel<ELLKernelTrait::setCSRValues<OtherValueType, ValueType> > setCSRValues;

    ContextPtr loc = testContext;
    setCSRValues.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    ValueType valuesCSRValues[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    const IndexType nCSRValues = sizeof( valuesCSRValues ) / sizeof( ValueType );
    IndexType valuesCSRIa[] =
    { 0, 5, 10, 15 };
    const IndexType nCSRIa = sizeof( valuesCSRIa ) / sizeof( IndexType );
    IndexType valuesCSRJa[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    const IndexType nCSRJa = sizeof( valuesCSRJa ) / sizeof( IndexType );
    IndexType valuesELLIa[] =
    { 5, 5, 5 };
    const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof( IndexType );
    OtherValueType expectedELLValues[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
    IndexType expectedELLJa[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
    const IndexType numRows = nELLIa;
    const IndexType nELLValues = 15;
    const IndexType numValuesPerRow = 5;
    HArray<ValueType> csrValues( nCSRValues, valuesCSRValues, testContext );
    HArray<IndexType> csrIa( nCSRIa, valuesCSRIa, testContext );
    HArray<IndexType> csrJa( nCSRJa, valuesCSRJa, testContext );
    HArray<IndexType> ellIa( nELLIa, valuesELLIa, testContext );
    // initialization of ellValues and ellJA, even if not mandatory
    HArray<OtherValueType> ellValues( nELLValues, OtherValueType( 0 ), testContext );
    HArray<IndexType> ellJa( nELLValues, IndexType( 0 ), testContext );
    {
        ReadAccess<ValueType> rCSRValues( csrValues, loc );
        ReadAccess<IndexType> rCSRIa( csrIa, loc );
        ReadAccess<IndexType> rCSRJa( csrJa, loc );
        ReadAccess<IndexType> rELLIa( ellIa, loc );
        WriteOnlyAccess<OtherValueType> wELLValues( ellValues, loc, nELLValues );
        WriteOnlyAccess<IndexType> wELLJa( ellJa, loc, nELLValues );
        SCAI_CONTEXT_ACCESS( loc );
        setCSRValues[loc]( wELLJa.get(), wELLValues.get(), rELLIa.get(), numRows, numValuesPerRow, rCSRIa.get(),
                           rCSRJa.get(), rCSRValues.get() );
    }
    ReadAccess<IndexType> rELLJa( ellJa );
    ReadAccess<OtherValueType> rELLValues( ellValues );

    for ( IndexType i = 0; i < nELLValues; i++ )
    {
        BOOST_CHECK_EQUAL( expectedELLJa[i], rELLJa[i] );
        BOOST_CHECK_EQUAL( expectedELLValues[i], rELLValues[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( compressIATest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();
    LAMAKernel<ELLKernelTrait::compressIA<ValueType> > compressIA;

    ContextPtr loc = testContext;
    compressIA.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    {
        HArray<IndexType> ellIa(  { 5, 5, 5 }, testContext );

        HArray<ValueType> ellValues( { 0, 0, 0,
                                       1, 1, 1,
                                       0, 0, 0,
                                       0, 0, 1,
                                       0, 1, 1 }, testContext );

        HArray<IndexType> ellJa(  { 0, 1, 2,
                                    3, 3, 3,
                                    4, 4, 4,
                                    5, 5, 5,
                                    6, 6, 6 }, testContext );

        const IndexType numRows = ellIa.size();
        const IndexType numValuesPerRow = ellValues.size() / numRows;
 
        BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, ellValues.size() );
        BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, ellJa.size() );

        const ValueType eps = 0.0;

        HArray<IndexType> newSizes1; // output array 1 ( keep diagonal )
        HArray<IndexType> newSizes2; // output array 2 ( do not keep diagonal )

        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );
            WriteOnlyAccess<IndexType> wNewSizes1( newSizes1, loc, numRows );
            WriteOnlyAccess<IndexType> wNewSizes2( newSizes2, loc, numRows );
            compressIA[loc]( wNewSizes1.get(), rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, true );
            compressIA[loc]( wNewSizes2.get(), rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, false );
        }

        BOOST_TEST( hostReadAccess( newSizes1 ) == std::vector<IndexType>(  { 2, 3, 4 } ), boost::test_tools::per_element() );
        BOOST_TEST( hostReadAccess( newSizes2 ) == std::vector<IndexType>(  { 1, 2, 3 } ), boost::test_tools::per_element() );
    }

    // Check with epsilon

    {
        HArray<IndexType> ellIa(  { 5, 5, 5 }, testContext );

        HArray<ValueType> ellValues( { 0.01,   0.02,  0.01,
                                       1,      1,     1,
                                       0.01,   0.01,  -0.01,
                                       -0.001, 0.001, 0.02,
                                       0.001,  1,     1         }, testContext );

        HArray<IndexType> ellJa(  { 0, 1, 2,
                                    3, 3, 3,
                                    4, 4, 4,
                                    5, 5, 5,
                                    6, 6, 6 }, testContext );

        const IndexType numRows = ellIa.size();
        const IndexType numValuesPerRow = ellValues.size() / numRows;

        BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, ellValues.size() );
        BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, ellJa.size() );

        const ValueType eps = 0.01;

        HArray<IndexType> newSizes1; // output array 1 ( keep diagonal )
        HArray<IndexType> newSizes2; // output array 2 ( do not keep diagonal )

        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );
            WriteOnlyAccess<IndexType> wNewSizes1( newSizes1, loc, numRows );
            WriteOnlyAccess<IndexType> wNewSizes2( newSizes2, loc, numRows );
            compressIA[loc]( wNewSizes1.get(), rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, true );
            compressIA[loc]( wNewSizes2.get(), rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, false );
        }

        BOOST_TEST( hostReadAccess( newSizes1 ) == std::vector<IndexType>(  { 2, 3, 4 } ), boost::test_tools::per_element() );
        BOOST_TEST( hostReadAccess( newSizes2 ) == std::vector<IndexType>(  { 1, 3, 3 } ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( compressValuesTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::compressValues<ValueType> > compressValues;

    ContextPtr loc = testContext;
    compressValues.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // Check without epsilon

    {
        /* Input Matrix:     1  4  0  0  0     0  3  4  5  6
                             2  5  0  0  8     1  3  4  5  6
                             3  6  0  7  9     2  3  4  5  6

           Output Matrix:    1  4  0  0        0  3  0  0
                             2  5  8  0        1  3  6  0
                             3  6  7  9        2  3  5  6
        */
        HArray<IndexType> ellIa( { 5, 5, 5 }, testContext );

        HArray<ValueType> ellValues( { 1, 2, 3,
                                       4, 5, 6,
                                       0, 0, 0,
                                       0, 0, 7,
                                       0, 8, 9  }, testContext );

        HArray<ValueType> expEllValues( { 1, 2, 3,
                                          4, 5, 6,
                                          0, 8, 7,
                                          0, 0, 9  }, testContext );

        HArray<IndexType> ellJa(    { 0, 1, 2,
                                      3, 3, 3,
                                      4, 4, 4,
                                      5, 5, 5,
                                      6, 6, 6  }, testContext );

        HArray<IndexType> expEllJa( { 0, 1, 2,
                                      3, 3, 3,
                                      0, 6, 5,
                                      0, 0, 6  }, testContext );

        const IndexType numRows = ellIa.size();
        const IndexType numValuesPerRow = ellValues.size() / numRows;
        const IndexType newNumValuesPerRow = expEllValues.size() / numRows;

        BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, ellValues.size() );
        BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, ellJa.size() );

        const ValueType eps = 0.0;

        SCAI_LOG_INFO( logger, "compress ELL, #rows = " << numRows << " #values/row = " << numValuesPerRow )

        HArray<IndexType> newEllJa( testContext );      // output array
        HArray<ValueType> newEllValues( testContext );  // output array

        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );
            WriteOnlyAccess<ValueType> wNewELLValues( newEllValues, loc, numRows * newNumValuesPerRow );
            WriteOnlyAccess<IndexType> wNewELLJa( newEllJa, loc, numRows * newNumValuesPerRow );
            SCAI_CONTEXT_ACCESS( loc );
            compressValues[loc]( wNewELLJa.get(), wNewELLValues.get(), newNumValuesPerRow,
                                 rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, true );
        }

        BOOST_TEST( hostReadAccess( newEllJa ) == hostReadAccess( expEllJa ), boost::test_tools::per_element() );
        BOOST_TEST( hostReadAccess( newEllValues ) == hostReadAccess( expEllValues ), boost::test_tools::per_element() );
    }

    // Check with epsilon
    {
        ValueType valuesELLValues[] = { 0.02,  2,     3,
                                        4,     5,     6,
                                        0.01,  -0.01, 0.002,
                                        -0.002, 0.01, 7,
                                        -0.01,  8,    9
                                      };
        const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof( ValueType );
        IndexType valuesELLIa[] = { 5, 5, 5 };
        const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof( IndexType );
        IndexType valuesELLJa[] = { 0, 1, 2,
                                    3, 3, 3,
                                    4, 4, 4,
                                    5, 5, 5,
                                    6, 6, 6
                                  };
        const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof( IndexType );
        ValueType expectedELLValues[] = { 0.02, 2, 3,
                                          4,    5, 6,
                                          0,    8, 7,
                                          0,    0, 9
                                        };
        IndexType expectedELLJa[] = { 0, 1, 2,
                                      3, 3, 3,
                                      0, 6, 5,
                                      0, 0, 6
                                    };
        const IndexType numRows = nELLIa;
        const IndexType numValuesPerRow = nELLJa / nELLIa;
        const ValueType eps = 0.01;
        const IndexType numValues = 12;
        const IndexType newNumValuesPerRow = numValues / nELLIa;
        HArray<ValueType> ellValues( nELLValues, valuesELLValues, testContext );
        HArray<IndexType> ellIa( nELLIa, valuesELLIa, testContext );
        HArray<IndexType> ellJa( nELLJa, valuesELLJa, testContext );
        HArray<IndexType> newEllJa( testContext );      // output array
        HArray<ValueType> newEllValues( testContext );  // output array
        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );
            WriteOnlyAccess<ValueType> wNewELLValues( newEllValues, loc, numValues );
            WriteOnlyAccess<IndexType> wNewELLJa( newEllJa, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            compressValues[loc]( wNewELLJa.get(), wNewELLValues.get(), newNumValuesPerRow,
                                 rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, true );
        }
        ReadAccess<ValueType> rNewELLValues( newEllValues );
        ReadAccess<IndexType> rNewELLJa( newEllJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            SCAI_LOG_DEBUG( logger, "Entry " << i << ", exp " << expectedELLJa[i] << ":" << expectedELLValues[i]
                            << ", is "  <<  rNewELLJa[i] << ":" << rNewELLValues[i] )
            BOOST_CHECK_EQUAL( expectedELLValues[i], rNewELLValues[i] );
            BOOST_CHECK_EQUAL( expectedELLJa[i], rNewELLJa[i] );
        }
    }
    // Check if compress destroys diagonal property (it shouldn't!)
    {
        ValueType valuesELLValues[] = { 0, 0, 0,
                                        4, 5, 6,
                                        0, 0, 0,
                                        0, 0, 7,
                                        0, 8, 9
                                      };
        const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof( ValueType );
        IndexType valuesELLIa[] = { 5, 5, 5 };
        const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof( IndexType );
        IndexType valuesELLJa[] = { 0, 1, 2,
                                    3, 3, 3,
                                    4, 4, 4,
                                    5, 5, 5,
                                    6, 6, 6
                                  };
        const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof( IndexType );
        ValueType expectedELLValues[] = { 0, 0, 0,
                                          4, 5, 6,
                                          0, 8, 7,
                                          0, 0, 9
                                        };
        IndexType expectedELLJa[] = { 0, 1, 2,
                                      3, 3, 3,
                                      0, 6, 5,
                                      0, 0, 6
                                    };
        const IndexType numRows = nELLIa;
        const ValueType eps = 0.01;
        const IndexType numValues = 12;
        const IndexType numValuesPerRow = nELLJa / nELLIa;
        const IndexType newNumValuesPerRow = numValues / nELLIa;
        HArray<ValueType> ellValues( nELLValues, valuesELLValues, testContext );
        HArray<IndexType> ellIa( nELLIa, valuesELLIa, testContext );
        HArray<IndexType> ellJa( nELLJa, valuesELLJa, testContext );
        HArray<IndexType> newEllJa( testContext );      // output array
        HArray<ValueType> newEllValues( testContext );  // output array
        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );
            WriteOnlyAccess<ValueType> wNewELLValues( newEllValues, loc, numValues );
            WriteOnlyAccess<IndexType> wNewELLJa( newEllJa, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            compressValues[loc]( wNewELLJa.get(), wNewELLValues.get(), newNumValuesPerRow,
                                 rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, true );
        }
        ReadAccess<ValueType> rNewELLValues( newEllValues );
        ReadAccess<IndexType> rNewELLJa( newEllJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            SCAI_LOG_DEBUG( logger, "Entry " << i << ", exp " << expectedELLJa[i] << ":" << expectedELLValues[i]
                            << ", is "  <<  rNewELLJa[i] << ":" << rNewELLValues[i] )
            BOOST_CHECK_EQUAL( expectedELLValues[i], rNewELLValues[i] );
            BOOST_CHECK_EQUAL( expectedELLJa[i], rNewELLJa[i] );
        }
    }
    // special case
    {
        ValueType valuesELLValues[] = { 1, 0, 0, 0, 0,
                                        0, 1, 1, 1, 1
                                      };
        const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof( ValueType ); // 10
        IndexType valuesELLIa[] = { 1, 2, 2, 2, 2 };
        const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof( IndexType ); // 5
        IndexType valuesELLJa[] = { 0, 0, 0, 0, 0,
                                    0, 1, 2, 3, 4
                                  };
        const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof( IndexType ); // 10
        ValueType expectedELLValues[] = { 1, 1, 1, 1, 1 };
        IndexType expectedELLJa[] = { 0, 1, 2, 3, 4 };
        const IndexType numRows = nELLIa;
        const ValueType eps = 0.01;
        const IndexType numValues = 10;
        const IndexType newNumValues = 5;
        const IndexType numValuesPerRow    = numValues    / nELLIa;
        const IndexType newNumValuesPerRow = newNumValues / nELLIa;
        HArray<ValueType> ellValues( nELLValues, valuesELLValues, testContext );
        HArray<IndexType> ellIa( nELLIa, valuesELLIa, testContext );
        HArray<IndexType> ellJa( nELLJa, valuesELLJa, testContext );
        HArray<IndexType> newEllJa( testContext );      // output array
        HArray<ValueType> newEllValues( testContext );  // output array
        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );
            WriteOnlyAccess<ValueType> wNewELLValues( newEllValues, loc, numValues );
            WriteOnlyAccess<IndexType> wNewELLJa( newEllJa, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            compressValues[loc]( wNewELLJa.get(), wNewELLValues.get(), newNumValuesPerRow,
                                 rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, true );
        }
        ReadAccess<ValueType> rNewELLValues( newEllValues );
        ReadAccess<IndexType> rNewELLJa( newEllJa );

        for ( IndexType i = 0; i < newNumValues; i++ )
        {
            SCAI_LOG_DEBUG( logger, "Entry " << i << ", exp " << expectedELLJa[i] << ":" << expectedELLValues[i]
                            << ", is "  <<  rNewELLJa[i] << ":" << rNewELLValues[i] )
            BOOST_CHECK_EQUAL( expectedELLValues[i], rNewELLValues[i] );
            BOOST_CHECK_EQUAL( expectedELLJa[i], rNewELLJa[i] );
        }
    }
}

BOOST_AUTO_TEST_CASE( matrixMultiplySizesTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::matrixMultiplySizes> matrixMultiplySizes;

    ContextPtr loc         = testContext;
    matrixMultiplySizes.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // Check with symmetric matrix
    {
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 3 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType expectedCIa[] =
        { 4, 4, 3, 4, 4 };
        IndexType numValues = 5; // all matrices have shape 5 x 5
        IndexType aNumValuesPerRow = aNumValues / numValues;
        IndexType bNumValuesPerRow = bNumValues / numValues;
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( testContext );
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            WriteOnlyAccess<IndexType> wCIa( CIa, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            matrixMultiplySizes[loc]( wCIa.get(), numValues, numValues, numValues, false, rAIa.get(), rAJa.get(),
                                      aNumValuesPerRow, rBIa.get(), rBJa.get(), bNumValuesPerRow );
        }
        BOOST_CHECK_EQUAL( numValues, CIa.size() );
        ReadAccess<IndexType> rCIa( CIa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCIa[i], rCIa[i] );
        }
    }
    // Check with asymmetric matrix
    {
        //   A       B
        //   x x 0   x  x
        //   0 x x   x
        //   0 0 x   x
        IndexType valuesAIa[] =
        { 2, 2, 2 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 0, 1, 0, 2, 3, 3 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        IndexType valuesBIa[] =
        { 2, 1, 2, 3 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 2 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType expectedCIa[] =
        { 3, 3, 3 };
        IndexType cNumRows = sizeof( expectedCIa ) / sizeof( IndexType );
        BOOST_REQUIRE_EQUAL( aNumRows, cNumRows );
        // a and a * b have same number rows
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( testContext );
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            WriteOnlyAccess<IndexType> wCIa( CIa, loc, cNumRows );
            SCAI_CONTEXT_ACCESS( loc );
            IndexType numColumns = 5; // does not really matter
            matrixMultiplySizes[loc]( wCIa.get(), aNumRows, numColumns, bNumRows, false, rAIa.get(), rAJa.get(),
                                      aNumValuesPerRow, rBIa.get(), rBJa.get(), bNumValuesPerRow );
        }
        BOOST_CHECK_EQUAL( cNumRows, CIa.size() );
        ReadAccess<IndexType> rCIa( CIa );

        for ( IndexType i = 0; i < cNumRows; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCIa[i], rCIa[i] );
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( matrixMultiplyTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::matrixMultiply<ValueType> > matrixMultiply;

    ContextPtr loc         = testContext;
    matrixMultiply.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    // Check with symmetric matrix
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType valuesCIa[] =
        { 4, 4, 3, 4, 4 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 28, 71, 8, 84, 73, 27, 63, 27, 144, 153, 15, 99, 37, 77, 85, 8, 56, 0, 49, 112 };
        IndexType expectedCJa[] =
        { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 2, 2, 3, 3, 0, 3, 3 };
        IndexType numValues = 20;
        ValueType alpha = 1;
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays
        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = numValues / cNumRows;
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, numValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            IndexType numColumns = 5; // not really needed here but internally used
            bool diagonalProperty = false; // do not care about it here
            matrixMultiply[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, bNumRows,
                                 diagonalProperty, alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow,
                                 rBIa.get(), rBJa.get(), rBValues.get(), bNumValuesPerRow );
        }
        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
    // Check with set alpha
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType valuesCIa[] =
        { 4, 4, 3, 4, 4 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 28, 71, 8, 84, 73, 27, 63, 27, 144, 153, 15, 99, 37, 77, 85, 8, 56, 0, 49, 112 };
        IndexType expectedCJa[] =
        { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 2, 2, 3, 3, 0, 3, 3 };
        IndexType numValues = 20;
        ValueType alpha = 2.5;
        // input arrays, directly initialized on context device
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays
        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = numValues / cNumRows;
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, numValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            bool diagonalProperty = false; // do not care about it
            IndexType numColumns = 15; // does not matter here but internally used for optimizations
            matrixMultiply[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, bNumRows,
                                 diagonalProperty, alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow,
                                 rBIa.get(), rBJa.get(), rBValues.get(), bNumValuesPerRow );
        }
        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i]*alpha, rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
    // Check with asymmetric matrix
    {
        ValueType valuesAValues[] =
        { 2, 4, 4, 3, 1, 5 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 2, 2 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 0, 1, 0, 2, 3, 3 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        BOOST_REQUIRE_EQUAL( aNumValues, nAValues );
        ValueType valuesBValues[] =
        { 4, 3, 7, 5, 9, 0, 6, 8, 0, 0, 0, 9 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 1, 2, 3 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 2 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        BOOST_REQUIRE_EQUAL( bNumValues, nBValues );
        IndexType valuesCIa[] =
        { 3, 3, 3 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 29, 5, 41, 18, 20, 40, 18, 9, 81 };
        IndexType expectedCJa[] =
        { 0, 0, 0, 1, 1, 1, 2, 2, 2 };
        IndexType cNumValues = 9;
        ValueType alpha = 1;
        // input arrays, directly initialized on context device
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays
        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = cNumValues / cNumRows;
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, cNumValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, cNumValues );
            SCAI_CONTEXT_ACCESS( loc );
            bool diagonalProperty = false; // do not care about it
            IndexType numColumns = 15; // does not matter here but internally used for optimizations
            matrixMultiply[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, bNumRows,
                                 diagonalProperty, alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow,
                                 rBIa.get(), rBJa.get(), rBValues.get(), bNumValuesPerRow );
        }
        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < cNumValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
}

BOOST_AUTO_TEST_CASE( matrixAddSizesTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::matrixAddSizes > matrixAddSizes;

    ContextPtr loc         = testContext;
    matrixAddSizes.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    IndexType valuesAIa[] =
    { 2, 3, 2, 3, 4 };
    const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
    IndexType valuesAJa[] =
    { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
    const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
    IndexType valuesBIa[] =
    { 2, 2, 2, 3, 3 };
    const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
    IndexType valuesBJa[] =
    { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
    const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
    IndexType expectedCIa[] =
    { 4, 3, 3, 4, 4 };
    const IndexType expectedCNumRows = sizeof( expectedCIa ) / sizeof( IndexType );
    // for matrix add A and B must have same shape
    BOOST_REQUIRE_EQUAL( aNumRows, bNumRows );
    IndexType cNumRows = aNumRows; // C gets same shape as A and B
    BOOST_REQUIRE_EQUAL( cNumRows, expectedCNumRows );
    // values per row needed, verify that numValues is multiple of numRows
    IndexType aNumValuesPerRow = aNumValues / aNumRows;
    BOOST_REQUIRE_EQUAL( aNumRows * aNumValuesPerRow, aNumValues );
    IndexType bNumValuesPerRow = bNumValues / bNumRows;
    BOOST_REQUIRE_EQUAL( bNumRows * bNumValuesPerRow, bNumValues );
    HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
    HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
    HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
    HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
    HArray<IndexType> CIa( testContext );
    {
        ReadAccess<IndexType> rAIa( AIa, loc );
        ReadAccess<IndexType> rAJa( AJa, loc );
        ReadAccess<IndexType> rBIa( BIa, loc );
        ReadAccess<IndexType> rBJa( BJa, loc );
        WriteOnlyAccess<IndexType> wCIa( CIa, loc, cNumRows );
        SCAI_CONTEXT_ACCESS( loc );
        bool diagonalProperty = false;
        IndexType numColumns = aNumRows; // square matrices here
        matrixAddSizes[loc]( wCIa.get(), aNumRows, numColumns, diagonalProperty, rAIa.get(), rAJa.get(), aNumValuesPerRow,
                             rBIa.get(), rBJa.get(), bNumValuesPerRow );
    }
    ReadAccess<IndexType> rCIa( CIa );

    for ( IndexType i = 0; i < cNumRows; i++ )
    {
        BOOST_CHECK_EQUAL( expectedCIa[i], rCIa[i] );
    }
}

BOOST_AUTO_TEST_CASE( matrixAddTest )
{
    typedef DefaultReal ValueType;
    ContextPtr testContext = ContextFix::testContext;
    LAMAKernel<ELLKernelTrait::matrixAdd<ValueType> > matrixAdd;

    ContextPtr loc         = testContext;
    matrixAdd.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // Check with neutral beta
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType valuesCIa[] =
        { 4, 3, 3, 4, 4 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 3, 9, 11, 12, 3, 1, 15, 3, 9, 12, 8, 8, 7, 12, 15, 3, 0, 0, 9, 7 };
        IndexType expectedCJa[] =
        { 0, 0, 1, 0, 0, 1, 3, 2, 1, 2, 2, 4, 3, 2, 3, 3, 0, 0, 3, 4 };
        IndexType cNumValues = 20;
        IndexType numColumns = 5; // for convenience
        ValueType alpha = 1;
        ValueType beta = 1;
        // input arrays, directly initialized on context device
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays, CValues, CJa
        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = cNumValues / cNumRows;
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, cNumValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, cNumValues );
            SCAI_CONTEXT_ACCESS( loc );
            bool diagonalProperty = false; // does not matter here
            matrixAdd[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, diagonalProperty,
                            alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow, beta, rBIa.get(), rBJa.get(),
                            rBValues.get(), bNumValuesPerRow );

        }

        // sort the columns, otherwise comparison might fail

        {
            ReadAccess<IndexType> rCIA( CIa );
            WriteAccess<IndexType> wCJA( CJa );
            WriteAccess<ValueType> wCValues( CValues );

            OpenMPELLUtils::sortRowElements( wCJA.get(), wCValues.get(), rCIA.get(), cNumRows, cNumValuesPerRow, false );
        }

        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < cNumValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
    // Check with set beta
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType valuesCIa[] =
        { 4, 3, 3, 4, 4 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 6, 13, 20, 20, 3, 1, 23, 3, 18, 15, 16, 8, 14, 17, 22, 3, 0, 0, 9, 7 };
        IndexType expectedCJa[] =
        { 0, 0, 1, 0, 0, 1, 3, 2, 1, 2, 2, 4, 3, 2, 3, 3, 0, 0, 3, 4 };
        IndexType cNumValues = 20;
        IndexType numColumns = 5; // for convenience
        ValueType alpha = 1;
        ValueType beta = 2;
        // input arrays, directly initialized on context device
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays

        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );

        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = cNumValues / cNumRows;

        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, cNumValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, cNumValues );
            SCAI_CONTEXT_ACCESS( loc );
            bool diagonalProperty = false; // does not matter here
            matrixAdd[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, diagonalProperty,
                            alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow, beta, rBIa.get(), rBJa.get(),
                            rBValues.get(), bNumValuesPerRow );

        }
        // sort the columns, otherwise comparison might fail

        {
            ReadAccess<IndexType> rCIA( CIa );
            WriteAccess<IndexType> wCJA( CJa );
            WriteAccess<ValueType> wCValues( CValues );

            OpenMPELLUtils::sortRowElements( wCJA.get(), wCValues.get(), rCIA.get(), cNumRows, cNumValuesPerRow, false );
        }

        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < cNumValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getValuePosColTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc         = testContext;
    getValuePosCol.getSupportedContext( loc );


    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    //    1.0   -   2.0
    //    0.5  0.3   -
    //     -    -   3.0

    const IndexType ia[] = { 2, 2, 1 };
    //  not this way: const IndexType ja[] = { 0, 2, 0, 1, 2, invalidIndex };
    const IndexType ja[] = { 0, 0, 2, 2, 1, invalidIndex };

    const IndexType numRows = 3;
    const IndexType numValuesPerRow = 2;

    HArray<IndexType> ellIA( numRows, ia, testContext );
    HArray<IndexType> ellJA( numRows * numValuesPerRow, ja, testContext );

    HArray<IndexType> row;   // result for rowIndexes
    HArray<IndexType> pos;   // result for positions

    IndexType cnt;

    IndexType columnIndex = 1;   // has 1 entry

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( ellIA, loc );
        ReadAccess<IndexType> rJA( ellJA, loc );
        WriteOnlyAccess<IndexType> wRow( row, loc, numRows );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numRows );
        cnt = getValuePosCol[loc]( wRow.get(), wPos.get(), columnIndex, rIA.get(), numRows, rJA.get(), numValuesPerRow );
    }

    BOOST_REQUIRE_EQUAL( cnt, IndexType( 1 ) );   //  only one entry for column 1

    {
        ReadAccess<IndexType> rPos( pos );
        ReadAccess<IndexType> rRow( row );

        BOOST_CHECK_EQUAL( IndexType( 1 ), rRow[0] );   // is in entry row
        BOOST_CHECK_EQUAL( IndexType( 4 ), rPos[0] );   // value of for (1,1) is at pos 4
    }

    columnIndex = 2;
    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( ellIA, loc );
        ReadAccess<IndexType> rJA( ellJA, loc );
        WriteOnlyAccess<IndexType> wRow( row, loc, numRows );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numRows );
        cnt = getValuePosCol[loc]( wRow.get(), wPos.get(), columnIndex, rIA.get(), numRows, rJA.get(), numValuesPerRow );
    }

    BOOST_REQUIRE_EQUAL( cnt, IndexType( 2 ) );   //  two entries for column 2, order might be arbitrary

    {
        ReadAccess<IndexType> rPos( pos );
        ReadAccess<IndexType> rRow( row );
        ReadAccess<IndexType> rJA( ellJA );
        ReadAccess<IndexType> rIA( ellIA );

        for ( IndexType k = 0; k < cnt; ++k )
        {
            IndexType p = rPos[k];
            IndexType i = rRow[k];
            BOOST_CHECK_EQUAL( rJA[ p ], columnIndex );
            BOOST_CHECK_EQUAL( p % numRows, i );
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( compressTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    LAMAKernel<ELLKernelTrait::compressIA<ValueType> > compressIA;
    LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    LAMAKernel<ELLKernelTrait::compressValues<ValueType> > compressValues;

    ContextPtr loc         = testContext;
    compressIA.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // full test (ia and ja/values) for
    // special case
    {
        ValueType valuesELLValues[] = { 1, 0, 0, 0, 0,
                                        0, 1, 1, 1, 1
                                      };
        IndexType valuesELLIa[]     = { 1, 2, 2, 2, 2 };
        IndexType valuesELLJa[]     = { 0, 0, 0, 0, 0,
                                        0, 1, 2, 3, 4
                                      };

        ValueType expectedELLValues[] = { 1, 1, 1, 1, 1 };
        IndexType expectedELLJa[]     = { 0, 1, 2, 3, 4 };
        IndexType expectedELLIa[]     = { 1, 1, 1, 1, 1 };

        const IndexType numRows = sizeof( valuesELLIa ) / sizeof( IndexType ); // 5
        const ValueType eps = 0.01;
        const IndexType numValues = sizeof( valuesELLValues ) / sizeof( ValueType ); // 10
        const IndexType newNumValues = sizeof( expectedELLValues ) / sizeof( ValueType ); // 5
        const IndexType numValuesPerRow    = numValues    / numRows;
        const IndexType newNumValuesPerRow = newNumValues / numRows;

        IndexType newNumValuesPerRow_calc = invalidIndex;
        HArray<ValueType> ellValues( numValues, valuesELLValues, testContext );
        HArray<IndexType> ellIa( numRows, valuesELLIa, testContext );
        HArray<IndexType> ellJa( numValues, valuesELLJa, testContext );
        HArray<IndexType> newEllIa( testContext ); // output array
        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );
            WriteOnlyAccess<IndexType> wNewELLIa( newEllIa, loc, numRows );
            SCAI_CONTEXT_ACCESS( loc );
            compressIA[loc]( wNewELLIa.get(), rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, true );
            newNumValuesPerRow_calc = reduce[loc]( wNewELLIa.get(), numRows, 0, common::BinaryOp::MAX );
        }
        ReadAccess<IndexType> rNewELLIa( newEllIa );

        for ( IndexType i = 0; i < numRows; i++ )
        {
            BOOST_CHECK_EQUAL( expectedELLIa[i], rNewELLIa[i] );
        }

        BOOST_CHECK_EQUAL ( newNumValuesPerRow, newNumValuesPerRow_calc );

        HArray<IndexType> newEllJa( testContext );      // output array
        HArray<ValueType> newEllValues( testContext );  // output array

        if ( newNumValuesPerRow_calc < numValuesPerRow )
        {
            if ( newNumValuesPerRow < numValuesPerRow )
            {
                {
                    ReadAccess<ValueType> rELLValues( ellValues, loc );
                    ReadAccess<IndexType> rELLIa( ellIa, loc );
                    ReadAccess<IndexType> rELLJa( ellJa, loc );
                    WriteOnlyAccess<ValueType> wNewELLValues( newEllValues, loc, numValues );
                    WriteOnlyAccess<IndexType> wNewELLJa( newEllJa, loc, numValues );
                    SCAI_CONTEXT_ACCESS( loc );
                    compressValues[loc]( wNewELLJa.get(), wNewELLValues.get(), newNumValuesPerRow, 
                                         rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, numValuesPerRow, eps, true );
                }
            }
        }

        ReadAccess<ValueType> rNewELLValues( newEllValues );
        ReadAccess<IndexType> rNewELLJa( newEllJa );

        for ( IndexType i = 0; i < newNumValues; i++ )
        {
            SCAI_LOG_DEBUG( logger, "Entry " << i << ", exp " << expectedELLJa[i] << ":" << expectedELLValues[i]
                            << ", is "  <<  rNewELLJa[i] << ":" << rNewELLValues[i] )
            BOOST_CHECK_EQUAL( expectedELLValues[i], rNewELLValues[i] );
            BOOST_CHECK_EQUAL( expectedELLJa[i], rNewELLJa[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();
    ContextPtr loc         = testContext;

    static LAMAKernel<ELLKernelTrait::normalGEMV<ValueType> > normalGEMV;

    normalGEMV.getSupportedContext( loc );
    BOOST_WARN_EQUAL( loc, testContext );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    SCAI_ASSERT_EQ_ERROR( ellIA.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), numRows * numValuesPerRow, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellValues.size(), numRows * numValuesPerRow, "size mismatch" )

    const ValueType y_values[]   = { 1, -1, 2, -2, 1, 1, -1 };
    const ValueType x_values[]   = { 3, -3, 2, -2 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_y   = sizeof( y_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numColumns, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, n_y, "size mismatch" );

    HArray<ValueType> x( numColumns, x_values, testContext );
    HArray<ValueType> y( numRows, y_values, testContext );

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        HArray<ValueType> res( testContext );

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * ELL * x + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", ELL: ia = " << ellIA << ", ja = " << ellJA << ", values = " << ellValues )
        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( ellIA, loc );
            ReadAccess<IndexType> rJA( ellJA, loc );
            ReadAccess<ValueType> rValues( ellValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            ReadAccess<ValueType> rY( y, loc );
            WriteOnlyAccess<ValueType> wResult( res, loc, numRows );

            normalGEMV[loc]( wResult.get(),
                             alpha, rX.get(), beta, rY.get(),
                             numRows, numColumns, numValuesPerRow, rIA.get(), rJA.get(), rValues.get(), common::MatrixOp::NORMAL );
        }

        HArray<ValueType> expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gevmTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();
    ContextPtr loc         = testContext;

    static LAMAKernel<ELLKernelTrait::normalGEMV<ValueType> > normalGEMV;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    SCAI_ASSERT_EQ_ERROR( ellIA.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), numValuesPerRow * numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellValues.size(), numValuesPerRow * numRows, "size mismatch" )

    const ValueType y_values[]   = { 1, -1, 2, -2 };
    const ValueType x_values[]   = { 3, -2, -2, 3, 1, 0, 1 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_y   = sizeof( y_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numRows, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, n_y, "size mismatch" );

    HArray<ValueType> x( numRows, x_values, testContext );
    HArray<ValueType> y( numColumns, y_values, testContext );

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        HArray<ValueType> res( testContext );

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * x * ELL + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", ELL: ia = " << ellIA << ", ja = " << ellJA << ", values = " << ellValues )
        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( ellIA, loc );
            ReadAccess<IndexType> rJA( ellJA, loc );
            ReadAccess<ValueType> rValues( ellValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            ReadAccess<ValueType> rY( y, loc );
            WriteOnlyAccess<ValueType> wResult( res, loc, numColumns );

            normalGEMV[loc]( wResult.get(),
                             alpha, rX.get(), beta, rY.get(),
                             numRows, numColumns, numValuesPerRow, rIA.get(), rJA.get(), rValues.get(), common::MatrixOp::TRANSPOSE );
        }

        HArray<ValueType> expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( spGEMVTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<ELLKernelTrait::sparseGEMV<ValueType> > sparseGEMV;

    ContextPtr loc = testContext;

    sparseGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;
    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    data1::getRowIndexes( rowIndexes );
    IndexType numNonEmptyRows = rowIndexes.size();

    SCAI_ASSERT_EQ_ERROR( numRows, ellIA.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( numRows * numValuesPerRow, ellJA.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), ellValues.size(), "size mismatch" )

    const ValueType res_values[]   = { 1, -1, 2, -2, 1, 1, -1 };
    const ValueType x_values[]     = { 3, -3, 2, -2 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_res = sizeof( res_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numColumns, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, n_res, "size mismatch" );

    HArray<ValueType> x( numColumns, x_values, testContext );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha; ++icase )
    {
        ValueType alpha = alpha_values[icase];

        HArray<ValueType> res( numRows, res_values, testContext );

        SCAI_LOG_INFO( logger, "compute res += " << alpha << " * CSR * x "
                       << ", with x = " << x
                       << ", CSR: ia = " << ellIA << ", ja = " << ellJA << ", values = " << ellValues )
        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( ellIA, loc );
            ReadAccess<IndexType> rJA( ellJA, loc );
            ReadAccess<ValueType> rValues( ellValues, loc );
            ReadAccess<IndexType> rIndexes( rowIndexes, loc );

            ReadAccess<ValueType> rX( x, loc );
            WriteAccess<ValueType> wResult( res, loc );

            sparseGEMV[loc]( wResult.get(),
                             alpha, rX.get(),
                             numRows, numValuesPerRow,
                             numNonEmptyRows, rIndexes.get(),
                             rIA.get(), rJA.get(), rValues.get(), common::MatrixOp::NORMAL );
        }

        HArray<ValueType> expectedRes( numRows, res_values );

        ValueType beta = 1;  // res = alpha * A * x + 1 * res <-> res += alpha * A * x

        expectedRes = data1::getGEMVNormalResult( alpha, x, beta, expectedRes );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( spGEVMTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<ELLKernelTrait::sparseGEMV<ValueType> > sparseGEMV;

    ContextPtr loc = testContext;

    sparseGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;
    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    data1::getRowIndexes( rowIndexes );
    IndexType numNonEmptyRows = rowIndexes.size();

    SCAI_ASSERT_EQ_ERROR( ellIA.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), numRows * numValuesPerRow, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellValues.size(), ellJA.size(), "size mismatch" )

    HArray<ValueType> x( { 3, -2, -2, 3, 1, 0, 1 }, testContext );
    HArray<ValueType> y( { 1, -1, 2, -2 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numRows, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, y.size(), "size mismatch" );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha; ++icase )
    {
        ValueType alpha = alpha_values[icase];

        HArray<ValueType> res( y );

        SCAI_LOG_INFO( logger, "compute res += " << alpha << " * CSR * x "
                       << ", with x = " << x
                       << ", CSR: ia = " << ellIA << ", ja = " << ellJA << ", values = " << ellValues )
        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( ellIA, loc );
            ReadAccess<IndexType> rJA( ellJA, loc );
            ReadAccess<ValueType> rValues( ellValues, loc );
            ReadAccess<IndexType> rIndexes( rowIndexes, loc );

            ReadAccess<ValueType> rX( x, loc );
            WriteAccess<ValueType> wResult( res, loc );

            sparseGEMV[loc]( wResult.get(), alpha, rX.get(),
                             numRows, numValuesPerRow,
                             numNonEmptyRows, rIndexes.get(),
                             rIA.get(), rJA.get(), rValues.get(), common::MatrixOp::TRANSPOSE );
        }

        ValueType beta = 1;  // res = alpha * x * A + 1 * res <-> res += alpha * x * A

        HArray<ValueType> expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<ELLKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = testContext;

    jacobi.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext << " on " << *loc )

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data2::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    const ValueType rhs_values[]   = { 1, -1, 2, -2 };
    const ValueType old_values[]   = { 3, -2, -2, 3 };

    HArray<ValueType> rhs( numRows, rhs_values, testContext );
    HArray<ValueType> oldSolution( numRows, old_values, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> res( testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( ellIA, loc );
            ReadAccess<IndexType> rJA( ellJA, loc );
            ReadAccess<ValueType> rValues( ellValues, loc );

            ReadAccess<ValueType> rOld( oldSolution, loc );
            ReadAccess<ValueType> rRhs( rhs, loc );
            WriteOnlyAccess<ValueType> wSolution( res, loc, numColumns );

            jacobi[loc]( wSolution.get(), numRows,
                         numValuesPerRow, rIA.get(), rJA.get(), rValues.get(),
                         rOld.get(), rRhs.get(), omega );

        }

        HArray<ValueType> expectedRes( testContext );

        data2::getJacobiResult( expectedRes, oldSolution, omega, rhs );

        auto maxDiff = HArrayUtils::maxDiffNorm( expectedRes, res );

        BOOST_CHECK( maxDiff < 0.1 );

        bool mustBeIdentical = false;

        if ( mustBeIdentical )
        {
            ReadAccess<ValueType> rExpected( expectedRes );
            ReadAccess<ValueType> rComputed( res );

            for ( IndexType i = 0; i < numRows; ++i )
            {
                BOOST_CHECK_EQUAL( rExpected[i], rComputed[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiHaloTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<ELLKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = testContext;

    jacobiHalo.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );

    SCAI_LOG_INFO( logger, "jacobiHalo test for " << *testContext << " on " << *loc )

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;
    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    data1::getRowIndexes( rowIndexes );

    IndexType numNonEmptyRows = rowIndexes.size();

    HArray<ValueType> oldSolution( { 3, -2, -2, 3, 1, 0, 2 }, testContext );
    HArray<ValueType> diag( { 9,  8,  7, 6, 7, 8, 9 }, testContext );

    SCAI_ASSERT_EQ_ERROR( oldSolution.size(), numRows, "test size mismatch" )
    SCAI_ASSERT_EQ_ERROR( diag.size(), numRows, "test size mismatch" )

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> solution( numRows, ValueType( 0 ), testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( ellIA, loc );
            ReadAccess<IndexType> rJA( ellJA, loc );
            ReadAccess<ValueType> rValues( ellValues, loc );
            ReadAccess<IndexType> rIndexes( rowIndexes, loc );

            ReadAccess<ValueType> rOld( oldSolution, loc );
            ReadAccess<ValueType> rDiag( diag, loc );
            WriteAccess<ValueType> wSolution( solution, loc, numColumns );

            jacobiHalo[loc]( wSolution.get(), numRows, rDiag.get(),
                             numValuesPerRow, rIA.get(), rJA.get(), rValues.get(),
                             rIndexes.get(), numNonEmptyRows,
                             rOld.get(), omega );
        }

        HArray<ValueType> expectedSol( numRows, ValueType( 0 ), testContext );

        data1::getJacobiHaloResult( expectedSol, oldSolution, diag, omega );

        auto maxDiff = HArrayUtils::maxDiffNorm( expectedSol, solution );

        BOOST_CHECK( maxDiff < 0.1 );

        bool mustBeIdentical = false;

        if ( mustBeIdentical )
        {
            ReadAccess<ValueType> rExpected( expectedSol );
            ReadAccess<ValueType> rComputed( solution );

            for ( IndexType i = 0; i < numRows; ++i )
            {
                BOOST_CHECK_EQUAL( rExpected[i], rComputed[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxValTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    LAMAKernel<ELLKernelTrait::absMaxVal<ValueType> > absMaxVal;

    ContextPtr loc = testContext;
    absMaxVal.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    ValueType maxVal = 0;

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rIA( ellIA, loc );
        ReadAccess<ValueType> rValues( ellValues, loc );

        maxVal = absMaxVal[loc]( numRows, numValuesPerRow, rIA.get(), rValues.get() );
    }

    ValueType expectedMaxVal = data1::getMaxVal<ValueType>();

    BOOST_CHECK_EQUAL( expectedMaxVal, maxVal );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()
