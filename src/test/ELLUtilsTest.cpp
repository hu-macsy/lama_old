/**
 * @file ELLUtilsTest.cpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Contains tests for the class CUDAELLUtils and OpenMPELLUtils
 * @author: Jan Ecker
 * @date 15.10.2012
 * $Id$
 **/

// boost
#include <boost/test/unit_test.hpp>

// others
#include <lama/ContextAccess.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/LAMAArray.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/ReadAccess.hpp>
#include <lama/WriteAccess.hpp>

#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

/* ------------------------------------------------------------------------------------------------------------------ */

// Dummy type, needed to use the lama interface
typedef bool NoType;

/* ------------------------------------------------------------------------------------------------------------------ */

namespace lama
{
namespace ELLUtilsTest
{

template<typename NoType>
void countNonEmptyRowsBySizesTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN( countNonEmptyRowsBySizes, loc, ELLUtils, Operations );

    // count valid array
    {
        const IndexType values[] =
        { 3, 0, 1, 0, 0, 1, 0, 4 };
        const IndexType n = sizeof( values ) / sizeof(IndexType);

        LAMAArray<IndexType> sizes( n, values );

        ReadAccess<IndexType> rSizes( sizes, loc );

        LAMA_CONTEXT_ACCESS( loc );

        IndexType count = countNonEmptyRowsBySizes( rSizes.get(), n );

        BOOST_CHECK_EQUAL( 4, count );
    }

    // count empty array
    {

        LAMAArray<IndexType> sizes;

        ReadAccess<IndexType> rSizes( sizes, loc );

        LAMA_CONTEXT_ACCESS( loc );

        IndexType count = countNonEmptyRowsBySizes( rSizes.get(), sizes.size() );

        BOOST_CHECK_EQUAL( 0, count );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename NoType>
void setNonEmptyRowsBySizesTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN( setNonEmptyRowsBySizes, loc, ELLUtils, Operations );

    const IndexType values[] =
    { 3, 0, 1, 0, 0, 1, 0, 4, 3, 0 };
    const IndexType valuesResult[] =
    { 0, 2, 5, 7, 8 };
    const IndexType n = 10;
    const IndexType numNonEmptyRows = 5;

    LAMAArray<IndexType> sizes( n, values );
    LAMAArray<IndexType> rowIndexes( numNonEmptyRows, 0 );

    {
        ReadAccess<IndexType> rSizes( sizes, loc );
        WriteAccess<IndexType> wRowIndexes( rowIndexes, loc );

        LAMA_CONTEXT_ACCESS( loc );

        setNonEmptyRowsBySizes( wRowIndexes.get(), numNonEmptyRows, rSizes.get(), n );
    }
    {
        HostReadAccess<IndexType> rRowIndexes( rowIndexes );

        for ( int i = 0; i < numNonEmptyRows; ++i )
        {
            BOOST_CHECK_EQUAL( valuesResult[i], rRowIndexes[i] );
        }
    }

}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename NoType>
void hasDiagonalPropertyTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN( hasDiagonalProperty, loc, ELLUtils, Operations );

    // positive test
    {
        const IndexType ellJaValues[] =
        { 0, 1, 2, 3, 4, 5, 6, 7, 5, 3, 9, 10, 7, 8, 9, 10 };
        const IndexType n = sizeof( ellJaValues ) / sizeof(IndexType);
        const IndexType numDiagonals = 8;

        LAMAArray<IndexType> ellJa( n, ellJaValues );
        ReadAccess<IndexType> rEllJa( ellJa, loc );

        LAMA_CONTEXT_ACCESS( loc );

        bool diagonalProperty = hasDiagonalProperty( numDiagonals, rEllJa.get() );

        BOOST_CHECK_EQUAL( true, diagonalProperty );
    }

    // negative test
    {
        const IndexType ellJaValues[] =
        { 0, 1, 2, 3, 7, 5, 6, 7, 5, 3, 9, 10, 7, 8, 9, 10 };
        const IndexType n = sizeof( ellJaValues ) / sizeof(IndexType);
        const IndexType numDiagonals = 8;

        LAMAArray<IndexType> ellJa( n, ellJaValues );
        ReadAccess<IndexType> rEllJa( ellJa, loc );

        LAMA_CONTEXT_ACCESS( loc );

        bool diagonalProperty = hasDiagonalProperty( numDiagonals, rEllJa.get() );

        BOOST_CHECK_EQUAL( false, diagonalProperty );
    }

    // test empty array
    {
        const IndexType numDiagonals = 0;

        LAMAArray<IndexType> ellJa;
        ReadAccess<IndexType> rEllJa( ellJa, loc );

        LAMA_CONTEXT_ACCESS( loc );

        bool diagonalProperty = hasDiagonalProperty( numDiagonals, rEllJa.get() );

        BOOST_CHECK_EQUAL( false, diagonalProperty );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename NoType>
void checkTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN( check, loc, ELLUtils, Operations );

    // check with correct values
    {
        const IndexType valuesIa[] =
        { 4, 3, 5, 2 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof(IndexType);
        const IndexType valuesJa[] =
        { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, -1, 4, -1, 4, 0, 0, 0, 5, 0 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof(IndexType);

        const IndexType numRows = nIa;
        const IndexType numValuesPerRow = 5;
        const IndexType numColumns = 6;

        LAMAArray<IndexType> ia( nIa, valuesIa );
        LAMAArray<IndexType> ja( nJa, valuesJa );

        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );

        LAMA_CONTEXT_ACCESS( loc );
        BOOST_CHECK_NO_THROW( check( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ) );
    }

    // check with invalid ia
    {
        const IndexType valuesIa[] =
        { 4, 3, 7, 2 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof(IndexType);
        const IndexType valuesJa[] =
        { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 0, 4, 0, 4, 0, 0, 0, 5, 0 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof(IndexType);

        const IndexType numRows = nIa;
        const IndexType numValuesPerRow = 5;
        const IndexType numColumns = 5;

        LAMAArray<IndexType> ia( nIa, valuesIa );
        LAMAArray<IndexType> ja( nJa, valuesJa );

        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );

        LAMA_CONTEXT_ACCESS( loc );
        BOOST_CHECK_THROW( check( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ),
                           Exception );
    }

    // check with invalid ja
    {
        const IndexType valuesIa[] =
        { 4, 3, 5, 2 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof(IndexType);
        const IndexType valuesJa[] =
        { 1, 1, 1, 1, 2, 2, 2, -1, 3, 3, 3, 0, 4, 0, 4, 0, 0, 0, 5, 0 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof(IndexType);

        const IndexType numRows = nIa;
        const IndexType numValuesPerRow = 5;
        const IndexType numColumns = 5;

        LAMAArray<IndexType> ia( nIa, valuesIa );
        LAMAArray<IndexType> ja( nJa, valuesJa );

        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );

        LAMA_CONTEXT_ACCESS( loc );
        BOOST_CHECK_THROW( check( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ),
                           Exception );
    }

    // check with valid empty values
    {
        const IndexType numRows = 0;
        const IndexType numValuesPerRow = 0;
        const IndexType numColumns = 0;

        LAMAArray<IndexType> ia;
        LAMAArray<IndexType> ja;

        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );

        LAMA_CONTEXT_ACCESS( loc );
        BOOST_CHECK_NO_THROW( check( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ) );
    }

    // check with invalid empty values
    {
        const IndexType numRows = 0;
        const IndexType numValuesPerRow = 1;
        const IndexType numColumns = 1;

        LAMAArray<IndexType> ia;
        LAMAArray<IndexType> ja;

        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );

        LAMA_CONTEXT_ACCESS( loc );
        BOOST_CHECK_THROW( check( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ),
                           Exception );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void getRowTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_TT( getRow, loc, ELLUtils, Getter, ValueType, OtherValueType );

    // check with valid dense values
    {
        ValueType valuesValues[] =
        { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof(ValueType);
        IndexType valuesIa[] =
        { 5, 5, 5 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof(IndexType);
        IndexType valuesJa[] =
        { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof(IndexType);
        OtherValueType expectedValues[] =
        { 0, 1, 2, 3, 4 };

        const IndexType i = 1;
        const IndexType numRows = nIa;
        const IndexType numColumns = 5;

        LAMAArray<ValueType> values( nValues, valuesValues );
        LAMAArray<IndexType> ia( nIa, valuesIa );
        LAMAArray<IndexType> ja( nJa, valuesJa );
        LAMAArray<OtherValueType> row( numColumns, 0.0 );

        {
            ReadAccess<ValueType> rValues( values, loc );
            ReadAccess<IndexType> rIa( ia, loc );
            ReadAccess<IndexType> rJa( ja, loc );
            WriteOnlyAccess<OtherValueType> wRow( row, loc, numColumns );

            LAMA_CONTEXT_ACCESS( loc );

            getRow( wRow.get(), i, numRows, numColumns, rIa.get(), rJa.get(), rValues.get() );
        }

        HostReadAccess<OtherValueType> rRow( row );

        for ( IndexType i = 0; i < numColumns; i++ )
        {
            BOOST_CHECK_EQUAL( expectedValues[i], rRow[i] );
        }
    }

    // check with valid sparse values
    {
        ValueType valuesValues[] =
        { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof(ValueType);
        IndexType valuesIa[] =
        { 5, 5, 5 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof(IndexType);
        IndexType valuesJa[] =
        { 0, 0, 0, 2, 2, 2, 4, 4, 4, 6, 6, 6, 10, 10, 10 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof(IndexType);
        OtherValueType expectedValues[] =
        { 0, 0, 1, 0, 2, 0, 3, 0, 0, 0, 4 };

        const IndexType i = 1;
        const IndexType numRows = nIa;
        const IndexType numColumns = 11;

        LAMAArray<ValueType> values( nValues, valuesValues );
        LAMAArray<IndexType> ia( nIa, valuesIa );
        LAMAArray<IndexType> ja( nJa, valuesJa );
        LAMAArray<OtherValueType> row( numColumns, 0.0 );

        {
            ReadAccess<ValueType> rValues( values, loc );
            ReadAccess<IndexType> rIa( ia, loc );
            ReadAccess<IndexType> rJa( ja, loc );
            WriteOnlyAccess<OtherValueType> wRow( row, loc, numColumns );

            LAMA_CONTEXT_ACCESS( loc );

            getRow( wRow.get(), i, numRows, numColumns, rIa.get(), rJa.get(), rValues.get() );
        }

        HostReadAccess<OtherValueType> rRow( row );

        for ( IndexType i = 0; i < numColumns; i++ )
        {
            BOOST_CHECK_EQUAL( expectedValues[i], rRow[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void getValueTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_TT( getValue, loc, ELLUtils, Getter, ValueType, OtherValueType );

    ValueType valuesValues[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
    const IndexType nValues = sizeof( valuesValues ) / sizeof(ValueType);
    IndexType valuesIa[] =
    { 5, 5, 5 };
    const IndexType nIa = sizeof( valuesIa ) / sizeof(IndexType);
    IndexType valuesJa[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
    const IndexType nJa = sizeof( valuesJa ) / sizeof(IndexType);
    OtherValueType expectedValues[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };

    const IndexType numRows = nIa;

    LAMAArray<ValueType> values( nValues, valuesValues );
    LAMAArray<IndexType> ia( nIa, valuesIa );
    LAMAArray<IndexType> ja( nJa, valuesJa );

    ReadAccess<ValueType> rValues( values, loc );
    ReadAccess<IndexType> rIa( ia, loc );
    ReadAccess<IndexType> rJa( ja, loc );

    LAMA_CONTEXT_ACCESS( loc );

    for ( IndexType i = 0; i < numRows; i++ )
    {
        for ( IndexType j = 0; j < valuesIa[i]; j++ )
        {
            OtherValueType result = getValue( i, j, numRows, rIa.get(), rJa.get(), rValues.get() );
            BOOST_CHECK_EQUAL( expectedValues[j * numRows + i], result );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void scaleValueTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_TT( scaleValue, loc, ELLUtils, Scale, ValueType, OtherValueType );

    ValueType mValues[] =
    { 1, 2, 3, 4, 5, 2, 2, 2, 2, 2, 4, 2, 0, 1, 3, 0, 0, 0, 0, 3 };
    const IndexType nValues = sizeof( mValues ) / sizeof(ValueType);
    const ValueType expectedValues[] =
    { 2, 4, 15, 8, 10, 4, 4, 10, 4, 4, 8, 4, 0, 2, 6, 0, 0, 0, 0, 6 };

    LAMAArray<ValueType> ellValues( nValues, mValues );

    {
        const IndexType numRows = 5;
        const IndexType ellIaValues[] =
        { 3, 3, 3, 3, 4 };
        const IndexType n = sizeof( ellIaValues ) / sizeof(IndexType);
        const OtherValueType values[] =
        { 2, 2, 5, 2, 2 };

        LAMAArray<IndexType> ellIa( n, ellIaValues );

        LAMAArray<OtherValueType> scaleValues( n, values );

        ReadAccess<IndexType> rEllIa( ellIa, loc );
        WriteAccess<ValueType> wEllValues( ellValues, loc );
        ReadAccess<OtherValueType> rScaleValues( scaleValues, loc );

        LAMA_CONTEXT_ACCESS( loc );

        scaleValue( numRows, rEllIa.get(), wEllValues.get(), rScaleValues.get() );
    }

    HostReadAccess<ValueType> rEllValues( ellValues );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expectedValues[i], rEllValues[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void getCSRValuesTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_TT( getCSRValues, loc, ELLUtils, Conversions, ValueType, OtherValueType );

    ValueType valuesELLValues[] = 
	{ 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };

    const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof(ValueType);
    IndexType valuesELLIa[] =
    { 5, 5, 5 };
    const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof(IndexType);
    IndexType valuesELLJa[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
    const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof(IndexType);
    IndexType valuesCSRIa[] =
    { 0, 5, 10, 15 };
    const IndexType nCSRIa = sizeof( valuesCSRIa ) / sizeof(IndexType);
    OtherValueType expectedCSRValues[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    IndexType expectedCSRJa[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };

    const IndexType numRows = nELLIa;
    const IndexType nCSRValues = 15;

    LAMAArray<ValueType> ellValues( nELLValues, valuesELLValues );
    LAMAArray<IndexType> ellIa( nELLIa, valuesELLIa );
    LAMAArray<IndexType> ellJa( nELLJa, valuesELLJa );
    LAMAArray<IndexType> csrIa( nCSRIa, valuesCSRIa );

    LAMAArray<OtherValueType> csrValues( nCSRValues, 0.0f );
    LAMAArray<IndexType> csrJa( nCSRValues, 0 );

    {
        ReadAccess<ValueType> rELLValues( ellValues, loc );
        ReadAccess<IndexType> rELLIa( ellIa, loc );
        ReadAccess<IndexType> rELLJa( ellJa, loc );
        ReadAccess<IndexType> rCSRIa( csrIa, loc );

        WriteOnlyAccess<OtherValueType> wCSRValues( csrValues, loc, nCSRValues );
        WriteOnlyAccess<IndexType> wCSRJa( csrJa, loc, nCSRValues );

        LAMA_CONTEXT_ACCESS( loc );

        getCSRValues( wCSRJa.get(), wCSRValues.get(), rCSRIa.get(), numRows, rELLIa.get(), rELLJa.get(),
                      rELLValues.get() );
    }

    HostReadAccess<IndexType> rCSRJa( csrJa );
    HostReadAccess<OtherValueType> rCSRValues( csrValues );

    for ( IndexType i = 0; i < nCSRValues; i++ )
    {
        BOOST_CHECK_EQUAL( expectedCSRJa[i], rCSRJa[i] );
        BOOST_CHECK_EQUAL( expectedCSRValues[i], rCSRValues[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void setCSRValuesTest( ContextPtr loc )
{
    // TODO: Change to use both types
    LAMA_INTERFACE_FN_TT( setCSRValues, loc, ELLUtils, Conversions, OtherValueType, ValueType );

    ValueType valuesCSRValues[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    const IndexType nCSRValues = sizeof( valuesCSRValues ) / sizeof(ValueType);
    IndexType valuesCSRIa[] =
    { 0, 5, 10, 15 };
    const IndexType nCSRIa = sizeof( valuesCSRIa ) / sizeof(IndexType);
    IndexType valuesCSRJa[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    const IndexType nCSRJa = sizeof( valuesCSRJa ) / sizeof(IndexType);
    IndexType valuesELLIa[] =
    { 5, 5, 5 };
    const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof(IndexType);
    OtherValueType expectedELLValues[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
    IndexType expectedELLJa[] =
    { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };

    const IndexType numRows = nELLIa;
    const IndexType nELLValues = 15;
    const IndexType numValuesPerRow = 5;

    LAMAArray<ValueType> csrValues( nCSRValues, valuesCSRValues );
    LAMAArray<IndexType> csrIa( nCSRIa, valuesCSRIa );
    LAMAArray<IndexType> csrJa( nCSRJa, valuesCSRJa );
    LAMAArray<IndexType> ellIa( nELLIa, valuesELLIa );

    LAMAArray<OtherValueType> ellValues( nELLValues, 0.0 );
    LAMAArray<IndexType> ellJa( nELLValues, 0.0 );

    {
        ReadAccess<ValueType> rCSRValues( csrValues, loc );
        ReadAccess<IndexType> rCSRIa( csrIa, loc );
        ReadAccess<IndexType> rCSRJa( csrJa, loc );
        ReadAccess<IndexType> rELLIa( ellIa, loc );

        WriteOnlyAccess<OtherValueType> wELLValues( ellValues, loc, nELLValues );
        WriteOnlyAccess<IndexType> wELLJa( ellJa, loc, nELLValues );

        LAMA_CONTEXT_ACCESS( loc );

        setCSRValues( wELLJa.get(), wELLValues.get(), rELLIa.get(), numRows, numValuesPerRow, rCSRIa.get(),
                      rCSRJa.get(), rCSRValues.get() );
    }

    HostReadAccess<IndexType> rELLJa( ellJa );
    HostReadAccess<OtherValueType> rELLValues( ellValues );

    for ( IndexType i = 0; i < nELLValues; i++ )
    {
        BOOST_CHECK_EQUAL( expectedELLJa[i], rELLJa[i] );
        BOOST_CHECK_EQUAL( expectedELLValues[i], rELLValues[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void compressIATest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( compressIA, loc, ELLUtils, Helper, ValueType );

    // Check without epsilon
    {
        ValueType valuesELLValues[] =
        { 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1 };
        const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof(ValueType);
        IndexType valuesELLIa[] =
        { 5, 5, 5 };
        const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof(IndexType);
        IndexType valuesELLJa[] =
        { 0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6 };
        const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof(IndexType);
        IndexType expectedELLIa[] =
        { 2, 3, 4 };

        const IndexType numRows = nELLIa;
        const ValueType eps = 0.0;

        LAMAArray<ValueType> ellValues( nELLValues, valuesELLValues );
        LAMAArray<IndexType> ellIa( nELLIa, valuesELLIa );
        LAMAArray<IndexType> ellJa( nELLJa, valuesELLJa );

        LAMAArray<IndexType> newEllIa( nELLIa, 0.0 );

        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );

            WriteOnlyAccess<IndexType> wNewELLIa( newEllIa, loc, nELLIa );

            LAMA_CONTEXT_ACCESS( loc );

            compressIA( rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, eps, wNewELLIa.get() );
        }

        HostReadAccess<IndexType> rNewELLIa( newEllIa );

        for ( IndexType i = 0; i < nELLIa; i++ )
        {
            BOOST_CHECK_EQUAL( expectedELLIa[i], rNewELLIa[i] );
        }
    }

    // Check with epsilon
    {
        ValueType valuesELLValues[] =
        { 1, 1, 1, 1, 1, 1, 0.01, 0.01, -0.01, -0.001, 0.001, 0.02, 0.001, 1, 1 };
        const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof(ValueType);
        IndexType valuesELLIa[] =
        { 5, 5, 5 };
        const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof(IndexType);
        IndexType valuesELLJa[] =
        { 0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6 };
        const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof(IndexType);
        IndexType expectedELLIa[] =
        { 2, 3, 4 };

        const IndexType numRows = nELLIa;
        const ValueType eps = 0.01;

        LAMAArray<ValueType> ellValues( nELLValues, valuesELLValues );
        LAMAArray<IndexType> ellIa( nELLIa, valuesELLIa );
        LAMAArray<IndexType> ellJa( nELLJa, valuesELLJa );

        LAMAArray<IndexType> newEllIa( nELLIa, 0.0 );

        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );

            WriteOnlyAccess<IndexType> wNewELLIa( newEllIa, loc, nELLIa );

            LAMA_CONTEXT_ACCESS( loc );

            compressIA( rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, eps, wNewELLIa.get() );
        }

        HostReadAccess<IndexType> rNewELLIa( newEllIa );

        for ( IndexType i = 0; i < nELLIa; i++ )
        {
            BOOST_CHECK_EQUAL( expectedELLIa[i], rNewELLIa[i] );
        }
    }

    // Check if compress destroys diagonal property (it shouldn't!)
    {
        ValueType valuesELLValues[] =
        { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof(ValueType);
        IndexType valuesELLIa[] =
        { 5, 5, 5 };
        const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof(IndexType);
        IndexType valuesELLJa[] =
        { 0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6 };
        const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof(IndexType);
        IndexType expectedELLIa[] =
        { 5, 5, 5 };

        const IndexType numRows = nELLIa;
        const ValueType eps = 0.0;

        LAMAArray<ValueType> ellValues( nELLValues, valuesELLValues );
        LAMAArray<IndexType> ellIa( nELLIa, valuesELLIa );
        LAMAArray<IndexType> ellJa( nELLJa, valuesELLJa );

        LAMAArray<IndexType> newEllIa( nELLIa, 0.0 );

        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );

            WriteOnlyAccess<IndexType> wNewELLIa( newEllIa, loc, nELLIa );

            LAMA_CONTEXT_ACCESS( loc );

            compressIA( rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, eps, wNewELLIa.get() );
        }

        HostReadAccess<IndexType> rNewELLIa( newEllIa );

        for ( IndexType i = 0; i < nELLIa; i++ )
        {
            BOOST_CHECK_EQUAL( expectedELLIa[i], rNewELLIa[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void compressValuesTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( compressValues, loc, ELLUtils, Helper, ValueType );

    // Check without epsilon
    {
        ValueType valuesELLValues[] =
        { 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 7, 0, 8, 9 };
        const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof(ValueType);
        IndexType valuesELLIa[] =
        { 5, 5, 5 };
        const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof(IndexType);
        IndexType valuesELLJa[] =
        { 0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6 };
        const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof(IndexType);
        ValueType expectedELLValues[] =
        { 1, 2, 3, 4, 5, 6, 0, 8, 7, 0, 0, 9 };
        IndexType expectedELLJa[] =
        { 0, 1, 2, 3, 3, 3, 0, 6, 5, 0, 0, 6 };

        const IndexType numRows = nELLIa;
        const ValueType eps = 0.0;
        const IndexType numValues = 12;

        LAMAArray<ValueType> ellValues( nELLValues, valuesELLValues );
        LAMAArray<IndexType> ellIa( nELLIa, valuesELLIa );
        LAMAArray<IndexType> ellJa( nELLJa, valuesELLJa );

        LAMAArray<ValueType> newEllValues( numValues, 0.0 );
        LAMAArray<IndexType> newEllJa( numValues, 0.0 );

        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );

            WriteOnlyAccess<ValueType> wNewELLValues( newEllValues, loc, numValues );
            WriteOnlyAccess<IndexType> wNewELLJa( newEllJa, loc, numValues );

            LAMA_CONTEXT_ACCESS( loc );

            compressValues( rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, eps, wNewELLJa.get(),
                            wNewELLValues.get() );
        }

        HostReadAccess<ValueType> rNewELLValues( newEllValues );
        HostReadAccess<IndexType> rNewELLJa( newEllJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedELLValues[i], rNewELLValues[i] );
            BOOST_CHECK_EQUAL( expectedELLJa[i], rNewELLJa[i] );
        }
    }

    // Check with epsilon
    {
        ValueType valuesELLValues[] =
        { 0.02, 2, 3, 4, 5, 6, 0.01, -0.01, 0.002, -0.002, 0.01, 7, -0.01, 8, 9 };
        const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof(ValueType);
        IndexType valuesELLIa[] =
        { 5, 5, 5 };
        const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof(IndexType);
        IndexType valuesELLJa[] =
        { 0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6 };
        const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof(IndexType);
        ValueType expectedELLValues[] =
        { 0.02, 2, 3, 4, 5, 6, 0, 8, 7, 0, 0, 9 };
        IndexType expectedELLJa[] =
        { 0, 1, 2, 3, 3, 3, 0, 6, 5, 0, 0, 6 };

        const IndexType numRows = nELLIa;
        const ValueType eps = 0.01;
        const IndexType numValues = 12;

        LAMAArray<ValueType> ellValues( nELLValues, valuesELLValues );
        LAMAArray<IndexType> ellIa( nELLIa, valuesELLIa );
        LAMAArray<IndexType> ellJa( nELLJa, valuesELLJa );

        LAMAArray<ValueType> newEllValues( numValues, 0.0 );
        LAMAArray<IndexType> newEllJa( numValues, 0.0 );

        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );

            WriteOnlyAccess<ValueType> wNewELLValues( newEllValues, loc, numValues );
            WriteOnlyAccess<IndexType> wNewELLJa( newEllJa, loc, numValues );

            LAMA_CONTEXT_ACCESS( loc );

            compressValues( rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, eps, wNewELLJa.get(),
                            wNewELLValues.get() );
        }

        HostReadAccess<ValueType> rNewELLValues( newEllValues );
        HostReadAccess<IndexType> rNewELLJa( newEllJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedELLValues[i], rNewELLValues[i] );
            BOOST_CHECK_EQUAL( expectedELLJa[i], rNewELLJa[i] );
        }
    }

    // Check if compress destroys diagonal property (it shouldn't!)
    {
        ValueType valuesELLValues[] =
        { 0, 0, 0, 4, 5, 6, 0, 0, 0, 0, 0, 7, 0, 8, 9 };
        const IndexType nELLValues = sizeof( valuesELLValues ) / sizeof(ValueType);
        IndexType valuesELLIa[] =
        { 5, 5, 5 };
        const IndexType nELLIa = sizeof( valuesELLIa ) / sizeof(IndexType);
        IndexType valuesELLJa[] =
        { 0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6 };
        const IndexType nELLJa = sizeof( valuesELLJa ) / sizeof(IndexType);
        ValueType expectedELLValues[] =
        { 0, 0, 0, 4, 5, 6, 0, 8, 7, 0, 0, 9 };
        IndexType expectedELLJa[] =
        { 0, 1, 2, 3, 3, 3, 0, 6, 5, 0, 0, 6 };

        const IndexType numRows = nELLIa;
        const ValueType eps = 0.01;
        const IndexType numValues = 12;

        LAMAArray<ValueType> ellValues( nELLValues, valuesELLValues );
        LAMAArray<IndexType> ellIa( nELLIa, valuesELLIa );
        LAMAArray<IndexType> ellJa( nELLJa, valuesELLJa );

        LAMAArray<ValueType> newEllValues( numValues, 0.0 );
        LAMAArray<IndexType> newEllJa( numValues, 0.0 );

        {
            ReadAccess<ValueType> rELLValues( ellValues, loc );
            ReadAccess<IndexType> rELLIa( ellIa, loc );
            ReadAccess<IndexType> rELLJa( ellJa, loc );

            WriteOnlyAccess<ValueType> wNewELLValues( newEllValues, loc, numValues );
            WriteOnlyAccess<IndexType> wNewELLJa( newEllJa, loc, numValues );

            LAMA_CONTEXT_ACCESS( loc );

            compressValues( rELLIa.get(), rELLJa.get(), rELLValues.get(), numRows, eps, wNewELLJa.get(),
                            wNewELLValues.get() );
        }

        HostReadAccess<ValueType> rNewELLValues( newEllValues );
        HostReadAccess<IndexType> rNewELLJa( newEllJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedELLValues[i], rNewELLValues[i] );
            BOOST_CHECK_EQUAL( expectedELLJa[i], rNewELLJa[i] );
        }
    }
}

template<typename NoType>
void computeIATest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( computeIA, loc, ELLUtils, MatrixTimesMatrix, IndexType );

    // Check with symmetric matrix
    {
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType nAIa = sizeof( valuesAIa ) / sizeof(IndexType);
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType nAJa = sizeof( valuesAJa ) / sizeof(IndexType);
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 3 };
        const IndexType nBIa = sizeof( valuesBIa ) / sizeof(IndexType);
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType nBJa = sizeof( valuesBJa ) / sizeof(IndexType);
        IndexType expectedCIa[] =
        { 4, 4, 3, 4, 4 };

        IndexType numValues = 5;

        LAMAArray<IndexType> AIa( nAIa, valuesAIa );
        LAMAArray<IndexType> AJa( nAJa, valuesAJa );
        LAMAArray<IndexType> BIa( nBIa, valuesBIa );
        LAMAArray<IndexType> BJa( nBJa, valuesBJa );

        LAMAArray<IndexType> CIa( numValues, 0.0 );
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );

            WriteOnlyAccess<IndexType> wCIa( CIa, loc, numValues );

            LAMA_CONTEXT_ACCESS( loc );

            computeIA( rAIa.get(), rAJa.get(), nAIa, rBIa.get(), rBJa.get(), nBIa, wCIa.get() );
        }

        HostReadAccess<IndexType> rCIa( CIa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCIa[i], rCIa[i] );
        }
    }

    // Check with asymmetric matrix
    {
        IndexType valuesAIa[] =
        { 2, 2, 2 };
        const IndexType nAIa = sizeof( valuesAIa ) / sizeof(IndexType);
        IndexType valuesAJa[] =
        { 0, 1, 0, 2, 3, 3 };
        const IndexType nAJa = sizeof( valuesAJa ) / sizeof(IndexType);
        IndexType valuesBIa[] =
        { 2, 1, 2, 3 };
        const IndexType nBIa = sizeof( valuesBIa ) / sizeof(IndexType);
        IndexType valuesBJa[] =
        { 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 2 };
        const IndexType nBJa = sizeof( valuesBJa ) / sizeof(IndexType);
        IndexType expectedCIa[] =
        { 3, 3, 3 };

        IndexType numValues = 3;

        LAMAArray<IndexType> AIa( nAIa, valuesAIa );
        LAMAArray<IndexType> AJa( nAJa, valuesAJa );
        LAMAArray<IndexType> BIa( nBIa, valuesBIa );
        LAMAArray<IndexType> BJa( nBJa, valuesBJa );

        LAMAArray<IndexType> CIa( numValues, 0.0 );
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );

            WriteOnlyAccess<IndexType> wCIa( CIa, loc, numValues );

            LAMA_CONTEXT_ACCESS( loc );

            computeIA( rAIa.get(), rAJa.get(), nAIa, rBIa.get(), rBJa.get(), nBIa, wCIa.get() );
        }

        HostReadAccess<IndexType> rCIa( CIa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCIa[i], rCIa[i] );
        }
    }
}

template<typename ValueType>
void computeValuesTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( computeValues, loc, ELLUtils, MatrixTimesMatrix, ValueType );

    // Check with symmetric matrix
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof(ValueType);
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType nAIa = sizeof( valuesAIa ) / sizeof(IndexType);
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType nAJa = sizeof( valuesAJa ) / sizeof(IndexType);
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof(ValueType);
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType nBIa = sizeof( valuesBIa ) / sizeof(IndexType);
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType nBJa = sizeof( valuesBJa ) / sizeof(IndexType);
        IndexType valuesCIa[] =
        { 4, 4, 3, 4, 4 };
        const IndexType nCIa = sizeof( valuesCIa ) / sizeof(IndexType);
        ValueType expectedCValues[] =
        { 28, 71, 8, 84, 73, 27, 63, 27, 144, 153, 15, 99, 37, 77, 85, 8, 56, 0, 49, 112 };
        IndexType expectedCJa[] =
        { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 2, 2, 3, 3, 0, 3, 3 };

        IndexType numValues = 20;
        ValueType alpha = 1;

        LAMAArray<ValueType> AValues( nAValues, valuesAValues );
        LAMAArray<IndexType> AIa( nAIa, valuesAIa );
        LAMAArray<IndexType> AJa( nAJa, valuesAJa );
        LAMAArray<ValueType> BValues( nBValues, valuesBValues );
        LAMAArray<IndexType> BIa( nBIa, valuesBIa );
        LAMAArray<IndexType> BJa( nBJa, valuesBJa );
        LAMAArray<IndexType> CIa( nCIa, valuesCIa );

        LAMAArray<ValueType> CValues( numValues, 0.0 );
        LAMAArray<IndexType> CJa( numValues, 0.0 );

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

            LAMA_CONTEXT_ACCESS( loc );

            computeValues( rAIa.get(), rAJa.get(), rAValues.get(), nAIa, rBIa.get(), rBJa.get(), rBValues.get(), nBIa,
                           alpha, rCIa.get(), wCJa.get(), wCValues.get() );
        }

        HostReadAccess<ValueType> rCValues( CValues );
        HostReadAccess<IndexType> rCJa( CJa );

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
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof(ValueType);
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType nAIa = sizeof( valuesAIa ) / sizeof(IndexType);
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType nAJa = sizeof( valuesAJa ) / sizeof(IndexType);
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof(ValueType);
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType nBIa = sizeof( valuesBIa ) / sizeof(IndexType);
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType nBJa = sizeof( valuesBJa ) / sizeof(IndexType);
        IndexType valuesCIa[] =
        { 4, 4, 3, 4, 4 };
        const IndexType nCIa = sizeof( valuesCIa ) / sizeof(IndexType);
        ValueType expectedCValues[] =
        { 28, 71, 8, 84, 73, 27, 63, 27, 144, 153, 15, 99, 37, 77, 85, 8, 56, 0, 49, 112 };
        IndexType expectedCJa[] =
        { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 2, 2, 3, 3, 0, 3, 3 };

        IndexType numValues = 20;
        ValueType alpha = 2.5;

        LAMAArray<ValueType> AValues( nAValues, valuesAValues );
        LAMAArray<IndexType> AIa( nAIa, valuesAIa );
        LAMAArray<IndexType> AJa( nAJa, valuesAJa );
        LAMAArray<ValueType> BValues( nBValues, valuesBValues );
        LAMAArray<IndexType> BIa( nBIa, valuesBIa );
        LAMAArray<IndexType> BJa( nBJa, valuesBJa );
        LAMAArray<IndexType> CIa( nCIa, valuesCIa );

        LAMAArray<ValueType> CValues( numValues, 0.0 );
        LAMAArray<IndexType> CJa( numValues, 0.0 );

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

            LAMA_CONTEXT_ACCESS( loc );

            computeValues( rAIa.get(), rAJa.get(), rAValues.get(), nAIa, rBIa.get(), rBJa.get(), rBValues.get(), nBIa,
                           alpha, rCIa.get(), wCJa.get(), wCValues.get() );
        }

        HostReadAccess<ValueType> rCValues( CValues );
        HostReadAccess<IndexType> rCJa( CJa );

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
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof(ValueType);
        IndexType valuesAIa[] =
        { 2, 2, 2 };
        const IndexType nAIa = sizeof( valuesAIa ) / sizeof(IndexType);
        IndexType valuesAJa[] =
        { 0, 1, 0, 2, 3, 3 };
        const IndexType nAJa = sizeof( valuesAJa ) / sizeof(IndexType);
        ValueType valuesBValues[] =
        { 4, 3, 7, 5, 9, 0, 6, 8, 0, 0, 0, 9 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof(ValueType);
        IndexType valuesBIa[] =
        { 2, 1, 2, 3 };
        const IndexType nBIa = sizeof( valuesBIa ) / sizeof(IndexType);
        IndexType valuesBJa[] =
        { 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 2 };
        const IndexType nBJa = sizeof( valuesBJa ) / sizeof(IndexType);
        IndexType valuesCIa[] =
        { 3, 3, 3 };
        const IndexType nCIa = sizeof( valuesCIa ) / sizeof(IndexType);
        ValueType expectedCValues[] =
        { 29, 5, 41, 18, 20, 40, 18, 9, 81 };
        IndexType expectedCJa[] =
        { 0, 0, 0, 1, 1, 1, 2, 2, 2 };

        IndexType numValues = 9;
        ValueType alpha = 1;

        LAMAArray<ValueType> AValues( nAValues, valuesAValues );
        LAMAArray<IndexType> AIa( nAIa, valuesAIa );
        LAMAArray<IndexType> AJa( nAJa, valuesAJa );
        LAMAArray<ValueType> BValues( nBValues, valuesBValues );
        LAMAArray<IndexType> BIa( nBIa, valuesBIa );
        LAMAArray<IndexType> BJa( nBJa, valuesBJa );
        LAMAArray<IndexType> CIa( nCIa, valuesCIa );

        LAMAArray<ValueType> CValues( numValues, 0.0 );
        LAMAArray<IndexType> CJa( numValues, 0.0 );

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

            LAMA_CONTEXT_ACCESS( loc );

            computeValues( rAIa.get(), rAJa.get(), rAValues.get(), nAIa, rBIa.get(), rBJa.get(), rBValues.get(), nBIa,
                           alpha, rCIa.get(), wCJa.get(), wCValues.get() );
        }

        HostReadAccess<ValueType> rCValues( CValues );
        HostReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
}

template<typename NoType>
void addComputeIATest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( addComputeIA, loc, ELLUtils, MatrixTimesMatrix, IndexType );

    {
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType nAIa = sizeof( valuesAIa ) / sizeof(IndexType);
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType nAJa = sizeof( valuesAJa ) / sizeof(IndexType);
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 3 };
        const IndexType nBIa = sizeof( valuesBIa ) / sizeof(IndexType);
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType nBJa = sizeof( valuesBJa ) / sizeof(IndexType);
        IndexType expectedCIa[] =
        { 4, 3, 3, 4, 4 };

        IndexType numValues = 5;

        LAMAArray<IndexType> AIa( nAIa, valuesAIa );
        LAMAArray<IndexType> AJa( nAJa, valuesAJa );
        LAMAArray<IndexType> BIa( nBIa, valuesBIa );
        LAMAArray<IndexType> BJa( nBJa, valuesBJa );

        LAMAArray<IndexType> CIa( numValues, 0.0 );
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );

            WriteOnlyAccess<IndexType> wCIa( CIa, loc, numValues );

            LAMA_CONTEXT_ACCESS( loc );

            addComputeIA( rAIa.get(), rAJa.get(), nAIa, rBIa.get(), rBJa.get(), nBIa, wCIa.get() );
        }

        HostReadAccess<IndexType> rCIa( CIa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCIa[i], rCIa[i] );
        }
    }
}

template<typename ValueType>
void addComputeValuesTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( addComputeValues, loc, ELLUtils, MatrixTimesMatrix, ValueType );

    // Check with neutral beta
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof(ValueType);
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType nAIa = sizeof( valuesAIa ) / sizeof(IndexType);
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType nAJa = sizeof( valuesAJa ) / sizeof(IndexType);
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof(ValueType);
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType nBIa = sizeof( valuesBIa ) / sizeof(IndexType);
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType nBJa = sizeof( valuesBJa ) / sizeof(IndexType);
        IndexType valuesCIa[] =
        { 4, 3, 3, 4, 4 };
        const IndexType nCIa = sizeof( valuesCIa ) / sizeof(IndexType);
        ValueType expectedCValues[] =
        { 3, 9, 11, 12, 3, 1, 15, 3, 9, 12, 8, 8, 7, 12, 15, 3, 0, 0, 9, 7 };
        IndexType expectedCJa[] =
        { 0, 0, 1, 0, 0, 1, 3, 2, 1, 2, 2, 4, 3, 2, 3, 3, 0, 0, 3, 4 };

        IndexType numValues = 20;
        ValueType beta = 1;

        LAMAArray<ValueType> AValues( nAValues, valuesAValues );
        LAMAArray<IndexType> AIa( nAIa, valuesAIa );
        LAMAArray<IndexType> AJa( nAJa, valuesAJa );
        LAMAArray<ValueType> BValues( nBValues, valuesBValues );
        LAMAArray<IndexType> BIa( nBIa, valuesBIa );
        LAMAArray<IndexType> BJa( nBJa, valuesBJa );
        LAMAArray<IndexType> CIa( nCIa, valuesCIa );

        LAMAArray<ValueType> CValues( numValues, 0.0 );
        LAMAArray<IndexType> CJa( numValues, 0.0 );

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

            LAMA_CONTEXT_ACCESS( loc );

            addComputeValues( rAIa.get(), rAJa.get(), rAValues.get(), nAIa, rBIa.get(), rBJa.get(), rBValues.get(),
                              nBIa, beta, rCIa.get(), wCJa.get(), wCValues.get() );
        }

        HostReadAccess<ValueType> rCValues( CValues );
        HostReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }

    // Check with set beta
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof(ValueType);
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType nAIa = sizeof( valuesAIa ) / sizeof(IndexType);
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType nAJa = sizeof( valuesAJa ) / sizeof(IndexType);
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof(ValueType);
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType nBIa = sizeof( valuesBIa ) / sizeof(IndexType);
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType nBJa = sizeof( valuesBJa ) / sizeof(IndexType);
        IndexType valuesCIa[] =
        { 4, 3, 3, 4, 4 };
        const IndexType nCIa = sizeof( valuesCIa ) / sizeof(IndexType);
        ValueType expectedCValues[] =
        { 6, 13, 20, 20, 3, 1, 23, 3, 18, 15, 16, 8, 14, 17, 22, 3, 0, 0, 9, 7 };
        IndexType expectedCJa[] =
        { 0, 0, 1, 0, 0, 1, 3, 2, 1, 2, 2, 4, 3, 2, 3, 3, 0, 0, 3, 4 };

        IndexType numValues = 20;
        ValueType beta = 2;

        LAMAArray<ValueType> AValues( nAValues, valuesAValues );
        LAMAArray<IndexType> AIa( nAIa, valuesAIa );
        LAMAArray<IndexType> AJa( nAJa, valuesAJa );
        LAMAArray<ValueType> BValues( nBValues, valuesBValues );
        LAMAArray<IndexType> BIa( nBIa, valuesBIa );
        LAMAArray<IndexType> BJa( nBJa, valuesBJa );
        LAMAArray<IndexType> CIa( nCIa, valuesCIa );

        LAMAArray<ValueType> CValues( numValues, 0.0 );
        LAMAArray<IndexType> CJa( numValues, 0.0 );

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

            LAMA_CONTEXT_ACCESS( loc );

            addComputeValues( rAIa.get(), rAJa.get(), rAValues.get(), nAIa, rBIa.get(), rBJa.get(), rBValues.get(),
                              nBIa, beta, rCIa.get(), wCJa.get(), wCValues.get() );
        }

        HostReadAccess<ValueType> rCValues( CValues );
        HostReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
}

// TODO: add SPMV tests

}//namespace ELLUtilsTest
} //namespace lama

/* ------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( ELLUtilsTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.ELLUtilsTest" );

LAMA_AUTO_TEST_CASE_TDUMMY( countNonEmptyRowsBySizesTest, ELLUtilsTest );
LAMA_AUTO_TEST_CASE_TDUMMY( setNonEmptyRowsBySizesTest, ELLUtilsTest );
LAMA_AUTO_TEST_CASE_TDUMMY( hasDiagonalPropertyTest, ELLUtilsTest );
LAMA_AUTO_TEST_CASE_TDUMMY( checkTest, ELLUtilsTest );

// TODO: uncomment after implementing this for CUDA!
//LAMA_AUTO_TEST_CASE_TDUMMY( computeIATest, ELLUtilsTest );
//LAMA_AUTO_TEST_CASE_TDUMMY( addComputeIATest, ELLUtilsTest );
//
//LAMA_AUTO_TEST_CASE_T( compressIATest, ELLUtilsTest );
//LAMA_AUTO_TEST_CASE_T( compressValuesTest, ELLUtilsTest );
//LAMA_AUTO_TEST_CASE_T( computeValuesTest, ELLUtilsTest );
//LAMA_AUTO_TEST_CASE_T( addComputeValuesTest, ELLUtilsTest );

LAMA_AUTO_TEST_CASE_TT( getRowTest, ELLUtilsTest );
LAMA_AUTO_TEST_CASE_TT( getValueTest, ELLUtilsTest );
LAMA_AUTO_TEST_CASE_TT( scaleValueTest, ELLUtilsTest );
LAMA_AUTO_TEST_CASE_TT( getCSRValuesTest, ELLUtilsTest );
LAMA_AUTO_TEST_CASE_TT( setCSRValuesTest, ELLUtilsTest );

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
