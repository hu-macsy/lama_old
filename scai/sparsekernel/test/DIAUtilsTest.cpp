/**
 * @file DIAUtilsTest.cpp
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
 * @brief Contains tests for the DIAUtils interface to be tested on different devices
 * @author Thomas Brandes
 * @date 15.12.2016 
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/sparsekernel/DIAKernelTrait.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;

using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;

using common::TypeTraits;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DIAUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DIAUtilsTest" )

/* ------------------------------------------------------------------------------------- */

template<typename ValueType> 
static void getDIATestData( 
    IndexType& numRows, 
    IndexType& numColumns, 
    IndexType& numDiagonals,
    HArray<IndexType>& diaOffsets,
    HArray<ValueType>& diaValues )
{
    /*                -5 -4 -3 -2 -1  0  1  2  3 

        Matrix:     x  x  x  x  x  x  6  0  0  4 
                       x  x  x  x  x  7  0  0  0  x  
                          x  x  x  x  0  0  9  4  x  x  
                             x  x  x  2  5  0  3  x  x  x 
                                x  x  2  0  0  1  x  x  x  x 
                                   x  0  0  0  0  x  x  x  x  x 
                                      0  1  0  2  x  x  x  x  x  x  */

    const IndexType diag_offsets[] = { 0, 
                                       static_cast<IndexType>( -5 ), 
                                       static_cast<IndexType>( -4 ), 
                                       static_cast<IndexType>( -3 ),
                                       static_cast<IndexType>( -2 ), 
                                       static_cast<IndexType>( -1 ), 
                                       1, 3 };

    const ValueType x = 0;  // just a stupid value as these entries should never be used

    const ValueType diag_values[]  = { 6, 0, 9, 3, x, x, x,
                                       x, x, x, x, x, 0, 1,
                                       x, x, x, x, 2, 0, 0,
                                       x, x, x, 2, 0, 0, 2,
                                       x, x, 0, 5, 0, 0, x,
                                       x, 7, 0, 0, 1, x, x,
                                       0, 0, 4, x, x, x, x,
                                       4, x, x, x, x, x, x };

    numRows      = 7;
    numColumns   = 4;
    numDiagonals = sizeof( diag_offsets ) / sizeof( IndexType );

    const IndexType diag_nvalues = sizeof( diag_values ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( diag_nvalues, numRows * numDiagonals );

    diaOffsets.init( diag_offsets, numDiagonals );
    diaValues.init( diag_values, diag_nvalues );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getCSRTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DIAKernelTrait::getCSRSizes<ValueType> > getCSRSizes;
    static LAMAKernel<DIAKernelTrait::getCSRValues<ValueType, ValueType> > getCSRValues;
    static LAMAKernel<UtilKernelTrait::scan<IndexType> > scan;

    ContextPtr loc = testContext;

    getCSRSizes.getSupportedContext( loc, getCSRValues, scan );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "getCSRSizes/getCSRValues test for " << *testContext << " on " << *loc )

    /*                -5 -4 -3 -2 -1  0  1  2  3 

        Matrix:     x  x  x  x  x  x  6  0  0  4 
                       x  x  x  x  x  7  0  0  0  x  
                          x  x  x  x  0  0  9  4  x  x  
                             x  x  x  2  5  0  3  x  x  x 
                                x  x  2  0  0  1  x  x  x  x 
                                   x  0  0  0  0  x  x  x  x  x 
                                      0  1  0  2  x  x  x  x  x  x  */

    const IndexType diag_offsets[] = { 0, 
                                       static_cast<IndexType>( -5 ), 
                                       static_cast<IndexType>( -4 ), 
                                       static_cast<IndexType>( -3 ),
                                       static_cast<IndexType>( -2 ), 
                                       static_cast<IndexType>( -1 ), 
                                       1, 3 };

    const ValueType x = 0;

    const ValueType diag_values[]  = { 6, 0, 9, 3, x, x, x,
                                       x, x, x, x, x, 0, 1,
                                       x, x, x, x, 2, 0, 0,
                                       x, x, x, 2, 0, 0, 2,
                                       0, 0, 0, 5, 0, 0, 0,
                                       0, 7, 0, 0, 1, 0, 0,
                                       0, 0, 4, 0, 0, 0, 0,
                                       4, 0, 0, 0, 0, 0, 0 };

    const IndexType ia_values[]  = { 2,    1, 2,    3,       2,    0, 2 };
    const IndexType ja_values[]  = { 0, 3, 0, 2, 3, 3, 0, 1, 0, 3,    1, 3 };
    const ValueType csr_values[] = { 6, 4, 7, 9, 4, 3, 2, 5, 2, 1,    1, 2 };

    const IndexType numRows      = 7;
    const IndexType numColumns   = 4;
    const IndexType numDiagonals = sizeof( diag_offsets ) / sizeof( IndexType );

    const IndexType diag_nvalues = sizeof( diag_values ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( diag_nvalues, numRows * numDiagonals );

    HArray<ValueType> diaValues( numRows * numDiagonals, diag_values, testContext );
    HArray<IndexType> diaOffsets( numDiagonals, diag_offsets, testContext );

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    ValueType eps = 0;
    bool diagonalProperty = false;

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rOffsets( diaOffsets, loc );
        ReadAccess<ValueType> rValues( diaValues, loc );
        WriteOnlyAccess<IndexType> wIA( csrIA, loc, numRows );

        getCSRSizes[loc]( wIA.get(), diagonalProperty, numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get(), eps );
    }

    {
        ReadAccess<IndexType> rIA( csrIA, hostContext );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_CHECK_EQUAL( rIA[i], ia_values[i] );
        }
    }

    IndexType numValues = 0;

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<IndexType> wIA( csrIA, loc );
        wIA.resize( numRows + 1 );
        numValues = scan[loc]( wIA.get(), numRows );
    }

    BOOST_REQUIRE_EQUAL( 12, numValues );

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( csrIA, loc );
        ReadAccess<ValueType> rValues( diaValues, loc );
        ReadAccess<IndexType> rOffsets( diaOffsets, loc );
        WriteOnlyAccess<IndexType> wJA( csrJA, loc, numValues );
        WriteOnlyAccess<ValueType> wValues( csrValues, loc, numValues );

        getCSRValues[loc]( wJA.get(), wValues.get(), rIA.get(), diagonalProperty, 
                           numRows, numColumns, numDiagonals,
                           rOffsets.get(), rValues.get(), eps );
    }

    {
        ReadAccess<IndexType> rJA( csrJA, hostContext );
        ReadAccess<ValueType> rValues( csrValues, hostContext );

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( rJA[i], ja_values[i] );
            BOOST_CHECK_EQUAL( rValues[i], csr_values[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = testContext;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "normalGEMV test for " << *testContext << " on " << *loc )

    HArray<ValueType> diaValues( testContext );
    HArray<IndexType> diaOffsets( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

    ValueType alpha = 1;
    ValueType beta  = -1;

    const ValueType y_values[]  = { 1, -1, 2, -2, 1, 1, -1 };
    const ValueType x_values[]  = { 3, -3, 2, -2, 3, 1, 2 };
    const ValueType res_values[]  = { 9, 22, 8, -13, 3, -1, -6 };

    HArray<ValueType> x( numRows, x_values, testContext );
    HArray<ValueType> y( numRows, y_values, testContext );

    HArray<ValueType> res( testContext );

    SCAI_LOG_INFO( logger, "compute res = " << alpha << " * x + " << beta << " * y " 
                            << ", with x = " << x << ", y = " << y )
    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rOffsets( diaOffsets, loc );
        ReadAccess<ValueType> rValues( diaValues, loc );

        ReadAccess<ValueType> rX( x, loc );
        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( res, loc, numRows );

        normalGEMV[loc]( wResult.get(), 
                         alpha, rX.get(), beta, rY.get(), 
                         numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get() );

    }

    {
        ReadAccess<ValueType> rResult( res, hostContext );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_CHECK_EQUAL( rResult[i], res_values[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()

