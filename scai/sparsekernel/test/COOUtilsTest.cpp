/**
 * @file COOUtilsTest.cpp
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
 * @brief Contains tests for the COOUtils interface to be tested on different devices
 * @author: Thomas Brandes
 * @date 19.07.2013
 * @since 1.1.0
 **/

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/kregistry/KernelContextFunction.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>

using namespace scai::hmemo;

using namespace scai::kregistry;

using namespace scai::sparsekernel;

/* ------------------------------------------------------------------------------------------------------------------ */

// Dummy type, needed to use the lama interface
typedef bool NoType;

/* ------------------------------------------------------------------------------------------------------------------ */

namespace scai
{
namespace lama
{
namespace COOUtilsTest
{

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename NoType>
void offsets2iaTest( ContextPtr loc )
{
    KernelTraitContextFunction<COOKernelTrait::offsets2ia > offsets2ia;

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
        HArray<IndexType> offsets( numOffsets );

        initArray( offsets, offsets_values, numOffsets );

        HArray<IndexType> ia;
        ReadAccess<IndexType> rOffsets( offsets, loc );
        const IndexType numDiagonals = 0;
        {
            WriteOnlyAccess<IndexType> wIA( ia, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            offsets2ia[loc->getType()]( wIA.get(), numValues, rOffsets.get(), numRows, numDiagonals );
        }
        ReadAccess<IndexType> rIA( ia );

        for ( int i = 0; i < numValues; ++i )
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
        HArray<IndexType> offsets( numOffsets );

        initArray( offsets, offsets_values, numOffsets );

        HArray<IndexType> ia;
        ReadAccess<IndexType> rOffsets( offsets, loc );
        const IndexType numDiagonals = 3;
        {
            WriteOnlyAccess<IndexType> wIA( ia, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            offsets2ia[loc->getType()]( wIA.get(), numValues, rOffsets.get(), numRows, numDiagonals );
        }
        ReadAccess<IndexType> rIA( ia );

        for ( int i = 0; i < numValues; ++i )
        {
            // SCAI_LOG_TRACE( logger,  "rIA[" << i << "] = " << rIA[i] << ", expects " << ia_values[i] )
            BOOST_CHECK_EQUAL( rIA[i], ia_values[i] );
        }
    }
} // offsets2iaTest

template<typename NoType>
void setCSRDataTest( ContextPtr loc )
{
    KernelTraitContextFunction<COOKernelTrait::setCSRData<IndexType, IndexType> > setCSRData;

    // setCSRData is for conversion of CSR storage to COO storage
    // is usually just a copy but has some reordering if diagonal property is required
    // here we test only for csrJA
    {
        const IndexType offsets_values[] =
        { 0, 2, 5, 7, 9 };
        const IndexType csrja_values[] =
        { 0, 5, 1, 4, 5, 2, 0, 4, 3 };
        const IndexType cooja_values[] =
        { 0, 1, 2, 5, 4, 5, 0, 4, 3 };
        const IndexType numOffsets = sizeof( offsets_values ) / sizeof( IndexType );
        const IndexType numRows = numOffsets - 1;
        const IndexType numDiagonals = 3;
        const IndexType numValues = sizeof( csrja_values ) / sizeof( IndexType );
        // verify that offsets and ia fit
        BOOST_REQUIRE_EQUAL( numValues, offsets_values[numRows] );
        BOOST_REQUIRE( numDiagonals <= numRows );
        HArray<IndexType> offsets( numOffsets );

        initArray( offsets, offsets_values, numOffsets );

        HArray<IndexType> csrJA( numValues );

        initArray( csrJA, csrja_values, numValues );

        HArray<IndexType> cooJA;
        ReadAccess<IndexType> rOffsets( offsets, loc );
        ReadAccess<IndexType> rCSRJA( csrJA, loc );
        {
            WriteOnlyAccess<IndexType> wCOOJA( cooJA, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            setCSRData[loc->getType()]( wCOOJA.get(), rCSRJA.get(), numValues, rOffsets.get(), numRows, numDiagonals );
        }
        ReadAccess<IndexType> rCOOJA( cooJA );

        for ( int i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( rCOOJA[i], cooja_values[i] );
        }
    }
} // setCSRData

} /* end namespace COOUtilsTest */

} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( COOUtilsTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.COOUtilsTest" )

LAMA_AUTO_TEST_CASE_CTDUMMY( offsets2iaTest, COOUtilsTest )
LAMA_AUTO_TEST_CASE_CTDUMMY( setCSRDataTest, COOUtilsTest )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
