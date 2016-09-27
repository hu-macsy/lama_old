/**
 * @file CSRUtilsTest.cpp
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
 * @brief Contains tests for the CSRUtils interface to be tested on different devices
 * @author Thomas Brandes
 * @date 05.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/kregistry.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/common/test/TestMacros.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;
using namespace hmemo;
using namespace sparsekernel;
using common::TypeTraits;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CSRUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxDiffValTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();
    kregistry::KernelTraitContextFunction<CSRKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;
    ContextPtr loc = Context::getContextPtr( absMaxDiffVal.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "absMaxDiffVal< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
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
    HArray<IndexType> csrIA1( numRows + 1, ia1, testContext );
    HArray<IndexType> csrJA1( numValues1, ja1, testContext );
    HArray<ValueType> csrValues1( numValues1, values1, testContext );
    HArray<IndexType> csrIA2( numRows + 1, ia2, testContext );
    HArray<IndexType> csrJA2( numValues2, ja2, testContext );
    HArray<ValueType> csrValues2( numValues2, values2, testContext );
    ReadAccess<IndexType> rCSRIA1( csrIA1, loc );
    ReadAccess<IndexType> rCSRJA1( csrJA1, loc );
    ReadAccess<ValueType> rCSRValues1( csrValues1, loc );
    ReadAccess<IndexType> rCSRIA2( csrIA2, loc );
    ReadAccess<IndexType> rCSRJA2( csrJA2, loc );
    ReadAccess<ValueType> rCSRValues2( csrValues2, loc );
    SCAI_CONTEXT_ACCESS( loc );
    ValueType maxVal = absMaxDiffVal[loc->getType()]( numRows, false, rCSRIA1.get(), rCSRJA1.get(), rCSRValues1.get(), rCSRIA2.get(),
                       rCSRJA2.get(), rCSRValues2.get() );
    BOOST_CHECK_EQUAL( 3, maxVal );
    // rows are sorted, so we can also apply sortFlag = true
    maxVal = absMaxDiffVal[loc->getType()]( numRows, true, rCSRIA1.get(), rCSRJA1.get(), rCSRValues1.get(), rCSRIA2.get(),
                                            rCSRJA2.get(), rCSRValues2.get() );
    BOOST_CHECK_EQUAL( 3, maxVal );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeSquareTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();
    kregistry::KernelTraitContextFunction<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;
    ContextPtr loc = Context::getContextPtr( convertCSR2CSC.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected
    SCAI_LOG_INFO( logger, "transpose< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
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
    HArray<IndexType> csrIA( numRows + 1, ia1, testContext );
    HArray<IndexType> csrJA( numValues, ja1, testContext );
    HArray<ValueType> csrValues( numValues, values1, testContext );
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
        convertCSR2CSC[loc->getType()]( wCSCIA.get(), wCSCJA.get(), wCSCValues.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(), numRows,
                                        numColumns, numValues );
    }
    {
        ReadAccess<IndexType> rCSCIA( cscIA );
        WriteAccess<IndexType> wCSCJA( cscJA );
        WriteAccess<ValueType> wCSCValues( cscValues );

        for ( IndexType j = 0; j <= numColumns; ++j )
        {
            BOOST_CHECK_EQUAL( rCSCIA[j], ia2[j] );
        }

        // For comparison of cscJA and cscValue we need to sort it
        bool diagonalFlag = false;
        OpenMPCSRUtils::sortRowElements( wCSCJA.get(), wCSCValues.get(), rCSCIA.get(), numColumns, diagonalFlag );

        for ( IndexType j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCJA[j], ja2[j] );
        }

        for ( IndexType j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCValues[j], values2[j] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeNonSquareTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();
    kregistry::KernelTraitContextFunction<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;
    ContextPtr loc = Context::getContextPtr( convertCSR2CSC.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected
    SCAI_LOG_INFO( logger, "transpose< " << TypeTraits<ValueType>::id() << "> non-square test for " << *testContext << " on " << *loc )
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
    HArray<IndexType> csrIA( numRows + 1, ia1, testContext );
    HArray<IndexType> csrJA( numValues, ja1, testContext );
    HArray<ValueType> csrValues( numValues, values1, testContext );
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
        convertCSR2CSC[loc->getType()]( wCSCIA.get(), wCSCJA.get(), wCSCValues.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(), numRows,
                                        numColumns, numValues );
    }
    //  For comparison later we sort cscJA and cscValue
    kregistry::KernelTraitContextFunction<CSRKernelTrait::sortRowElements<ValueType> > sortRowElements;
    loc = Context::getContextPtr( sortRowElements.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected
    SCAI_LOG_INFO( logger, "sortRowElements< " << TypeTraits<ValueType>::id() << "> for " << *testContext << " on " << *loc )
    {
        ReadAccess<IndexType> rCSCIA( cscIA, loc );
        WriteAccess<IndexType> wCSCJA( cscJA, loc );
        WriteAccess<ValueType> wCSCValues( cscValues, loc );

        for ( IndexType j = 0; j <= numColumns; ++j )
        {
            BOOST_CHECK_EQUAL( rCSCIA[j], ia2[j] );
        }

        bool diagonalFlag = false;
        // For comparison of cscJA and cscValue we need to sort it
        sortRowElements[loc->getType()]( wCSCJA.get(), wCSCValues.get(), rCSCIA.get(), numColumns, diagonalFlag );
    }
    // check CSC for correctness, done on host
    {
        ContextPtr host = Context::getHostPtr();
        ReadAccess<IndexType> rCSCIA( cscIA, host );
        ReadAccess<IndexType> wCSCJA( cscJA, host );
        ReadAccess<ValueType> wCSCValues( cscValues, host );

        for ( IndexType j = 0; j <= numColumns; ++j )
        {
            BOOST_CHECK_EQUAL( rCSCIA[j], ia2[j] );
        }

        for ( IndexType j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCJA[j], ja2[j] );
        }

        for ( IndexType j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCValues[j], values2[j] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_ARITHMETIC_EXT_HOST> scai_ext_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( decompositionTest, ValueType, scai_ext_test_types )
{
    ContextPtr testContext = Context::getContextPtr();
    kregistry::KernelTraitContextFunction<CSRKernelTrait::decomposition<ValueType> > decomposition;
    ContextPtr loc = Context::getContextPtr( decomposition.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected
    SCAI_LOG_INFO( logger, "decomposition< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    const IndexType ia[] = { 0, 4, 8, 12, 15 };
    const IndexType ja[] = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3 };
    const ValueType values[] = { 3.0,  4.0, -5.0,  6.0,
                                 6.0,  5.0, -6.0, 5.0,
                                 9.0, -4.0,  2.0, 3.0,
                                       2.0, -3.0, 1.0 };
    const ValueType rhsValues[] = { 39.0, 43.0, 6.0, 13.0 };
    const ValueType solValues[] = { 1.0, 2.0, -2.0, 3.0 };
    const IndexType numRows = 4;
    const IndexType nnz = 15;

    HArray<IndexType> csrIA( numRows + 1, ia, testContext );
    HArray<IndexType> csrJA( nnz, ja, testContext );
    HArray<ValueType> csrValues( nnz, values, testContext );

    HArray<ValueType> rhs( numRows, rhsValues, testContext );
    HArray<ValueType> solution;

    {
        ReadAccess<IndexType> rCSRIA( csrIA, loc );
        ReadAccess<IndexType> rCSRJA( csrJA, loc );
        ReadAccess<ValueType> rCSRValues( csrValues, loc );
        ReadAccess<ValueType> rRHS( rhs, loc );
        WriteOnlyAccess<ValueType> wSol( solution, loc, numRows );
        SCAI_CONTEXT_ACCESS( loc );
        decomposition[loc->getType()]( wSol.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(),
            rRHS.get(), numRows, nnz, false );
    }

    {
        ContextPtr host = Context::getHostPtr();
        ReadAccess<ValueType> rSol( solution, host );

        for ( int i = 0; i < numRows; ++i )
        {
            ValueType x = rSol[i] - solValues[i];
            BOOST_CHECK_SMALL( common::Math::real( x ), common::TypeTraits<ValueType>::small() );
            BOOST_CHECK_SMALL( common::Math::imag( x ), common::TypeTraits<ValueType>::small() );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
