/**
 * @file CSRUtilsTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Contains tests for the CSRUtils interface to be tested on different devices
 * @author: Thomas Brandes
 * @date 05.07.2013
 * @since 1.0.1
 **/

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/lama/CSRKernelTrait.hpp>
#include <scai/lama/LAMAKernel.hpp>

#include <scai/lama/openmp/OpenMPCSRUtils.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

/* ------------------------------------------------------------------------------------------------------------------ */
namespace scai
{

namespace lama
{

namespace CSRUtilsTest
{

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType, typename OtherValueType>
void absMaxDiffValTest( ContextPtr loc )
{
    LAMAKernel<CSRKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;

    // input arrays
    //    Array1             Array2
    //
    //    1 2 3 0 0          1 2 0 0 0
    //    0 0 1 1 2          1 0 2 2 1
    const IndexType ia1[] =
    { 0, 3, 6 };
    const IndexType ja1[] =
    { 0, 1, 2, 2, 3, 4 };
    const IndexType ia2[] =
    { 0, 2, 6 };
    const IndexType ja2[] =
    { 0, 1, 0, 2, 3, 4 };
    const ValueType values1[] =
    { 1, 2, 3, 1, 1, 2 };
    const ValueType values2[] =
    { 1, 2, 1, 2, 2, 1 };
    const IndexType numRows = 2;
    // const IndexType numColumns = 5;
    const IndexType numValues1 = sizeof( ja1 ) / sizeof( IndexType );
    const IndexType numValues2 = sizeof( ja2 ) / sizeof( IndexType );
    HArray<IndexType> csrIA1( numRows + 1, ia1 );
    HArray<IndexType> csrJA1( numValues1, ja1 );
    HArray<ValueType> csrValues1( numValues1, values1 );
    HArray<IndexType> csrIA2( numRows + 1, ia2 );
    HArray<IndexType> csrJA2( numValues2, ja2 );
    HArray<ValueType> csrValues2( numValues2, values2 );
    ReadAccess<IndexType> rCSRIA1( csrIA1, loc );
    ReadAccess<IndexType> rCSRJA1( csrJA1, loc );
    ReadAccess<ValueType> rCSRValues1( csrValues1, loc );
    ReadAccess<IndexType> rCSRIA2( csrIA2, loc );
    ReadAccess<IndexType> rCSRJA2( csrJA2, loc );
    ReadAccess<ValueType> rCSRValues2( csrValues2, loc );
    SCAI_CONTEXT_ACCESS( loc );
    ValueType maxVal = absMaxDiffVal[loc]( numRows, false, rCSRIA1.get(), rCSRJA1.get(), rCSRValues1.get(), rCSRIA2.get(),
                                           rCSRJA2.get(), rCSRValues2.get() );
    BOOST_CHECK_EQUAL( 3, maxVal );
    // rows are sorted, so we can also apply sortFlag = true
    maxVal = absMaxDiffVal[loc]( numRows, true, rCSRIA1.get(), rCSRJA1.get(), rCSRValues1.get(), rCSRIA2.get(),
                                rCSRJA2.get(), rCSRValues2.get() );
    BOOST_CHECK_EQUAL( 3, maxVal );
}

template<typename ValueType>
void transposeTestSquare( ContextPtr loc )
{
    LAMAKernel<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;

    //  input array           transpose
    //    1.0   -   2.0       1.0  0.5   -
    //    0.5  0.3   -         -   0.3   -
    //     -    -   3.0       2.0   -   3.0
    const IndexType ia1[] =
    { 0, 2, 4, 5 };
    const IndexType ja1[] =
    { 0, 2, 0, 1, 2 };
    const IndexType ia2[] =
    { 0, 2, 3, 5 };
    const IndexType ja2[] =
    { 0, 1, 1, 0, 2 };
    const ValueType values1[] =
    { 1.0, 2.0, 0.5, 0.3, 3.0 };
    const ValueType values2[] =
    { 1.0, 0.5, 0.3, 2.0, 3.0 };
    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType numValues = 5;
    HArray<IndexType> csrIA( numRows + 1, ia1 );
    HArray<IndexType> csrJA( numValues, ja1 );
    HArray<ValueType> csrValues( numValues, values1 );
    HArray<IndexType> cscIA;
    HArray<IndexType> cscJA;
    HArray<ValueType> cscValues;
  
    {
        ReadAccess<IndexType> rCSRIA( csrIA, loc );
        ReadAccess<IndexType> rCSRJA( csrJA, loc );
        ReadAccess<ValueType> rCSRValues( csrValues, loc );
        WriteOnlyAccess<IndexType> wCSCIA( cscIA, loc, numColumns + 1 );
        WriteOnlyAccess<IndexType> wCSCJA( cscJA, loc, numValues );
        WriteOnlyAccess<ValueType> wCSCValues( cscValues, loc, numValues );
        SCAI_CONTEXT_ACCESS( loc );
        convertCSR2CSC[loc]( wCSCIA.get(), wCSCJA.get(), wCSCValues.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(), numRows,
                             numColumns, numValues );
    }
    {
        ReadAccess<IndexType> rCSCIA( cscIA );
        WriteAccess<IndexType> wCSCJA( cscJA );
        WriteAccess<ValueType> wCSCValues( cscValues );

        for ( int j = 0; j <= numColumns; ++j )
        {
            BOOST_CHECK_EQUAL( rCSCIA[j], ia2[j] );
        }

        // For comparison of cscJA and cscValue we need to sort it
        bool diagonalFlag = false;
        OpenMPCSRUtils::sortRowElements( wCSCJA.get(), wCSCValues.get(), rCSCIA.get(), numColumns, diagonalFlag );

        for ( int j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCJA[j], ja2[j] );
        }

        for ( int j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCValues[j], values2[j] );
        }
    }
}

template<typename ValueType>
void transposeTestNonSquare( ContextPtr loc )
{
    LAMAKernel<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;
    //  input array           transpose
    //    1.0   -   2.0       1.0  0.5   -    4.0
    //    0.5  0.3   -         -   0.3   -    1.5
    //     -    -   3.0       2.0   -   3.0    -
    //    4.0  1.5   -
    const IndexType ia1[] =
    { 0, 2, 4, 5, 7 };
    const IndexType ja1[] =
    { 0, 2, 0, 1, 2, 0, 1 };
    const IndexType ia2[] =
    { 0, 3, 5, 7 };
    const IndexType ja2[] =
    { 0, 1, 3, 1, 3, 0, 2 };
    const ValueType values1[] =
    { 1.0, 2.0, 0.5, 0.3, 3.0, 4.0, 1.5 };
    const ValueType values2[] =
    { 1.0, 0.5, 4.0, 0.3, 1.5, 2.0, 3.0 };
    const IndexType numRows = 4;
    const IndexType numColumns = 3;
    const IndexType numValues = 7;
    HArray<IndexType> csrIA( numRows + 1, ia1 );
    HArray<IndexType> csrJA( numValues, ja1 );
    HArray<ValueType> csrValues( numValues, values1 );
    HArray<IndexType> cscIA;
    HArray<IndexType> cscJA;
    HArray<ValueType> cscValues;

    // CSC <- transpose CSR

    {
        ReadAccess<IndexType> rCSRIA( csrIA, loc );
        ReadAccess<IndexType> rCSRJA( csrJA, loc );
        ReadAccess<ValueType> rCSRValues( csrValues, loc );
        WriteOnlyAccess<IndexType> wCSCIA( cscIA, loc, numColumns + 1 );
        WriteOnlyAccess<IndexType> wCSCJA( cscJA, loc, numValues );
        WriteOnlyAccess<ValueType> wCSCValues( cscValues, loc, numValues );
        SCAI_CONTEXT_ACCESS( loc );
        convertCSR2CSC[loc]( wCSCIA.get(), wCSCJA.get(), wCSCValues.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(), numRows,
                             numColumns, numValues );
    }
  
    // check CSC for correctness

    {
        ContextPtr host = Context::getHostPtr();

        ReadAccess<IndexType> rCSCIA( cscIA, host );
        WriteAccess<IndexType> wCSCJA( cscJA, host );
        WriteAccess<ValueType> wCSCValues( cscValues, host );

        for ( int j = 0; j <= numColumns; ++j )
        {
            BOOST_CHECK_EQUAL( rCSCIA[j], ia2[j] );
        }

        // For comparison of cscJA and cscValue we need to sort it
        bool diagonalFlag = false;
        OpenMPCSRUtils::sortRowElements( wCSCJA.get(), wCSCValues.get(), rCSCIA.get(), numColumns, diagonalFlag );

        for ( int j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCJA[j], ja2[j] );
        }

        for ( int j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCValues[j], values2[j] );
        }
    }
}

} /* end namespace CSRUtilsTest */

} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( CSRUtilsTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRUtilsTest" )

LAMA_AUTO_TEST_CASE_CTT( absMaxDiffValTest, CSRUtilsTest )
LAMA_AUTO_TEST_CASE_CT( transposeTestSquare, CSRUtilsTest, scai::lama )
LAMA_AUTO_TEST_CASE_CT( transposeTestNonSquare, CSRUtilsTest, scai::lama )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
