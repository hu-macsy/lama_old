/**
 * @file COOUtilsTest.cpp
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
 * @brief Contains tests for the COOUtils interface to be tested on different devices
 * @author Thomas Brandes
 * @date 19.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/kregistry/KernelContextFunction.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;

using namespace hmemo;
using namespace sparsekernel;
using namespace kregistry;
using common::TypeTraits;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( COOUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.COOUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( offsets2iaTest )
{
    ContextPtr testContext = ContextFix::testContext;
    KernelTraitContextFunction<COOKernelTrait::offsets2ia > offsets2ia;
    ContextPtr loc = Context::getContextPtr( offsets2ia.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "offsets2ia test for " << *testContext << " on " << *loc )
    // Test without diagonal property
    {
        const IndexType offsets_values[] =
        { 0, 2, 2, 3, 5 };
        const IndexType ia_values[] =
        { 0, 0, 2, 3, 3 };
        const IndexType numOffsets = sizeof( offsets_values ) / sizeof( IndexType );
        const IndexType numRows = numOffsets - 1;
        const IndexType numValues = sizeof( ia_values ) / sizeof( IndexType );
        // verify that offsets and ia fit
        BOOST_REQUIRE_EQUAL( numValues, offsets_values[numRows] );
        HArray<IndexType> offsets( numOffsets, offsets_values, testContext );
        HArray<IndexType> ia;
        ReadAccess<IndexType> rOffsets( offsets, loc );
        const IndexType numDiagonals = 0;
        {
            WriteOnlyAccess<IndexType> wIA( ia, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            offsets2ia[loc->getType()]( wIA.get(), numValues, rOffsets.get(), numRows, numDiagonals );
        }
        ReadAccess<IndexType> rIA( ia );

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( rIA[i], ia_values[i] );
        }
    }
    // Test with diagonal property
    {
        const IndexType offsets_values[] =
        { 0, 2, 5, 7, 9 };
        // Result with diagonals = 0: ia_values[] = { 0, 0, 1, 1, 1, 2, 2, 3, 3 };
        // But with 3 diagonals we get this:
        const IndexType ia_values[] =
        { 0, 1, 2, 0, 1, 1, 2, 3, 3 };
        const IndexType numOffsets = sizeof( offsets_values ) / sizeof( IndexType );
        const IndexType numRows = numOffsets - 1;
        const IndexType numValues = sizeof( ia_values ) / sizeof( IndexType );
        // verify that offsets and ia fit
        BOOST_REQUIRE_EQUAL( numValues, offsets_values[numRows] );
        HArray<IndexType> offsets( numOffsets, offsets_values, testContext );
        HArray<IndexType> ia;
        ReadAccess<IndexType> rOffsets( offsets, loc );
        const IndexType numDiagonals = 3;
        {
            WriteOnlyAccess<IndexType> wIA( ia, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            offsets2ia[loc->getType()]( wIA.get(), numValues, rOffsets.get(), numRows, numDiagonals );
        }
        ReadAccess<IndexType> rIA( ia );

        for ( IndexType i = 0; i < numValues; ++i )
        {
            // SCAI_LOG_TRACE( logger,  "rIA[" << i << "] = " << rIA[i] << ", expects " << ia_values[i] )
            BOOST_CHECK_EQUAL( rIA[i], ia_values[i] );
        }
    }
} // offsets2iaTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setCSRDataTest, ValueType, scai_numeric_test_types )
{
    typedef float CSRValueType;
    ContextPtr testContext = ContextFix::testContext;
    KernelTraitContextFunction<COOKernelTrait::setCSRData<ValueType, CSRValueType> > setCSRData;
    ContextPtr loc = Context::getContextPtr( setCSRData.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "setCSRData< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    // setCSRData is for conversion of CSR storage to COO storage
    // is usually just a copy but has some reordering if diagonal property is required
    // here we test only for csrJA
    {
        const IndexType offsets_values[] = { 0, 2, 5, 7, 9 };
        const CSRValueType csrja_values[]   = { 0, 5, 1, 4, 5, 2, 0, 4, 3 };
        const ValueType cooja_values[]   = { 0, 1, 2, 5, 4, 5, 0, 4, 3 };
        const IndexType numOffsets = sizeof( offsets_values ) / sizeof( IndexType );
        const IndexType numRows = numOffsets - 1;
        const IndexType numDiagonals = 3;
        const IndexType numValues = sizeof( csrja_values ) / sizeof( CSRValueType );
        // verify that offsets and ia fit
        BOOST_REQUIRE_EQUAL( numValues, offsets_values[numRows] );
        BOOST_REQUIRE( numDiagonals <= numRows );
        HArray<IndexType> offsets( numOffsets, offsets_values, testContext );
        HArray<CSRValueType> csrJA( numValues, csrja_values, testContext );
        HArray<ValueType> cooJA;
        ReadAccess<IndexType> rOffsets( offsets, loc );
        ReadAccess<CSRValueType> rCSRJA( csrJA, loc );
        {
            WriteOnlyAccess<ValueType> wCOOJA( cooJA, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            setCSRData[loc->getType()]( wCOOJA.get(), rCSRJA.get(), numValues, rOffsets.get(), numRows, numDiagonals );
        }
        ReadAccess<ValueType> rCOOJA( cooJA );

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( rCOOJA[i], cooja_values[i] );
        }
    }
} // setCSRData


/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getValuePosColTest )
{
    ContextPtr testContext = Context::getContextPtr();

    kregistry::KernelTraitContextFunction<COOKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = Context::getContextPtr( getValuePosCol.validContext( testContext->getType() ) );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    //    1.0   -   2.0      
    //    0.5  0.3   -       
    //     -    -   3.0      

    const IndexType ia[] = { 0, 0, 1, 1, 2 };
    const IndexType ja[] = { 0, 2, 0, 1, 2 };

    const IndexType numRows = 3;
    const IndexType numValues = 5;

    HArray<IndexType> cooIA( numValues, ia, testContext );
    HArray<IndexType> cooJA( numValues, ja, testContext );

    HArray<IndexType> row;   // result for rowIndexes 
    HArray<IndexType> pos;   // result for positions

    IndexType cnt;

    IndexType columnIndex = 1;   // has 1 entry

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );
        WriteOnlyAccess<IndexType> wRow( row, loc, numRows );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numRows );
        cnt = getValuePosCol[loc->getType()]( wRow.get(), wPos.get(), columnIndex, rIA.get(), numRows, rJA.get(), numValues );
    }

    BOOST_REQUIRE_EQUAL( cnt, 1 );   //  only one entry for column 1

    {
        ReadAccess<IndexType> rPos( pos );
        ReadAccess<IndexType> rRow( row );

        BOOST_CHECK_EQUAL( 1, rRow[0] );   // is in entry row
        BOOST_CHECK_EQUAL( 3, rPos[0] );   // value of for (1,1) is at pos 3
    }

    columnIndex = 2;
    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );
        WriteOnlyAccess<IndexType> wRow( row, loc, numRows );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numRows );
        cnt = getValuePosCol[loc->getType()]( wRow.get(), wPos.get(), columnIndex, rIA.get(), numRows, rJA.get(), numValues );
    }

    BOOST_REQUIRE_EQUAL( cnt, 2 );   //  two entries for column 2, order might be arbitrary

    {
        ReadAccess<IndexType> rPos( pos );
        ReadAccess<IndexType> rRow( row );
        ReadAccess<IndexType> rJA( cooJA );
        ReadAccess<IndexType> rIA( cooIA );

        for ( IndexType k = 0; k < cnt; ++k )
        {
            IndexType p = rPos[k];
            IndexType i = rRow[k];
            BOOST_CHECK_EQUAL( rJA[ p ], columnIndex );
            BOOST_CHECK( rIA[p] == i );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()

