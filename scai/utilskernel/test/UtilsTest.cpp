/**
 * @file UtilsTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
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
#include <scai/utilskernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/SparseKernelTrait.hpp>
#include <scai/hmemo.hpp>
#include <scai/common/TypeTraits.hpp>

// import scai_numeric_test_types, scai_array_test_types

#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace utilskernel;
using namespace hmemo;
using common::BinaryOp;
using common::CompareOp;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( UtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.UtilsTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleTest, ValueType, scai_array_test_types )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    setVal.getSupportedContext( loc );
    BOOST_WARN_EQUAL( loc.get(), testContext.get() );

    SCAI_LOG_INFO( logger, "scaleTest<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )
    ValueType valuesValues[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
    ValueType expectedValues[] =
    { 0, 2, 4, 6, 8, 0, 2, 4, 6, 8, 0, 2, 4, 6, 8 };
    const ValueType mult = 2;
    HArray<ValueType> values( nValues, valuesValues );
    {
        WriteAccess<ValueType> wValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        setVal[loc]( wValues.get(), nValues, mult, BinaryOp::MULT );
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

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    reduce.getSupportedContext( loc );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // print warning if not available for test context
    SCAI_LOG_INFO( logger, "sumTest<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )
    {
        ValueType valuesValues[] =
        { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
        const ValueType expectedSum = 30;
        HArray<ValueType> values( nValues, valuesValues );
        ReadAccess<ValueType> rValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        const ValueType resultSum = reduce[loc]( rValues.get(), nValues, ValueType( 0 ), BinaryOp::ADD );
        BOOST_CHECK_EQUAL( expectedSum, resultSum );
    }
    {
        const ValueType expectedSum = 0;
        HArray<ValueType> values;
        ReadAccess<ValueType> rValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        const ValueType resultSum = reduce[loc]( rValues.get(), values.size(), ValueType( 0 ), BinaryOp::ADD );
        BOOST_CHECK_EQUAL( expectedSum, resultSum );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( reduce2Test, ValueType, scai_array_test_types )
{
    static LAMAKernel<UtilKernelTrait::reduce2<ValueType> > reduce2;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    reduce2.getSupportedContext( loc );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // print warning if not available for test context
    SCAI_LOG_INFO( logger, "reduce2Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )
    {
        ValueType valuesValues1[] = { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
        ValueType valuesValues2[] = { 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0 };

        ValueType expectedVal = 3 * 4 ; // maximal product

        const IndexType nValues1 = sizeof( valuesValues1 ) / sizeof( ValueType );
        const IndexType nValues2 = sizeof( valuesValues2 ) / sizeof( ValueType );

        SCAI_ASSERT_EQ_ERROR( nValues1, nValues2, "size mismatch" )

        HArray<ValueType> values1( nValues1, valuesValues1 );
        HArray<ValueType> values2( nValues2, valuesValues2 );

        ReadAccess<ValueType> rValues1( values1, loc );
        ReadAccess<ValueType> rValues2( values2, loc );

        SCAI_CONTEXT_ACCESS( loc );

        BinaryOp binop = BinaryOp::MULT;
        BinaryOp redop = BinaryOp::MAX;

        ValueType zero = 0;

        const ValueType result = reduce2[loc]( rValues1.get(), rValues2.get(), nValues1, binop, zero, redop );

        BOOST_CHECK_EQUAL( expectedVal, result );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setValTest, ValueType, scai_array_test_types )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    setVal.getSupportedContext( loc );
    SCAI_LOG_INFO( logger, "setValTest<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // print warning if not available for test context
    {
        const IndexType n = 20;
        HArray<ValueType> values;
        {
            WriteOnlyAccess<ValueType> wValues( values, loc, 3 * n );
            SCAI_CONTEXT_ACCESS( loc );
            setVal[loc]( wValues.get(), 3 * n, 0, BinaryOp::COPY );
            // overwrite in the middle to check that there is no out-of-range set
            setVal[loc]( wValues.get() + n, n, 10, BinaryOp::COPY );
        }
        ReadAccess<ValueType> rValues( values );

        for ( IndexType i = 0; i < n; i++ )
        {
            BOOST_CHECK_EQUAL( ValueType( 0 ), rValues.get()[i + 0 * n] );
            BOOST_CHECK_EQUAL( ValueType( 10 ), rValues.get()[i + 1 * n] );
            BOOST_CHECK_EQUAL( ValueType( 0 ), rValues.get()[i + 2 * n] );
        }
    }
    {
        const IndexType n = 0;
        HArray<ValueType> values;
        {
            WriteOnlyAccess<ValueType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            setVal[loc]( wValues.get(), n, 7, BinaryOp::COPY );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( isSortedTest, ValueType, scai_array_test_types )
{
    static LAMAKernel<UtilKernelTrait::isSorted<ValueType> > isSorted;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

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
        HArray<ValueType> valueArray1( nValues1, values1 );
        HArray<ValueType> valueArray2( nValues2, values2 );
        HArray<ValueType> valueArray3( nValues3, values3 );
        ReadAccess<ValueType> rValues1( valueArray1, loc );
        ReadAccess<ValueType> rValues2( valueArray2, loc );
        ReadAccess<ValueType> rValues3( valueArray3, loc );
        SCAI_CONTEXT_ACCESS( loc );
        // values1 are sorted, operator = LE
        BOOST_CHECK( isSorted[loc]( rValues1.get(), nValues1, CompareOp::LE ) );
        BOOST_CHECK( ! isSorted[loc]( rValues1.get(), nValues1, CompareOp::GE ) );
        // values2 are sorted, ascending = false
        BOOST_CHECK( isSorted[loc]( rValues2.get(), nValues2, CompareOp::GE ) );
        BOOST_CHECK( ! isSorted[loc]( rValues2.get(), nValues2, CompareOp::LE ) );
        BOOST_CHECK( isSorted[loc]( rValues2.get(), 0, CompareOp::LE ) );
        // only first two values are sorted
        BOOST_CHECK( isSorted[loc]( rValues2.get(), 1, CompareOp::LE ) );
        // only first two values are sorted
        BOOST_CHECK( isSorted[loc]( rValues2.get(), 2, CompareOp::LE ) );
        // only first two values are sorted
        // values3 are not sorted, neither ascending nor descending
        BOOST_CHECK( ! isSorted[loc]( rValues3.get(), nValues3, CompareOp::GE ) );
        BOOST_CHECK( ! isSorted[loc]( rValues3.get(), nValues3, CompareOp::LE ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setOrderTest )
{
    // setOrder is only for IndexType

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    SCAI_LOG_INFO( logger, "setOrderTest on " << *loc )
    static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
    {
        const IndexType n = 20;
        HArray<IndexType> values;
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
        HArray<IndexType> values;
        {
            WriteOnlyAccess<IndexType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            setOrder[loc]( wValues.get(), n );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( binaryOpScalar1Test, ValueType, scai_numeric_test_types )
{
    static LAMAKernel<UtilKernelTrait::binaryOpScalar<ValueType> > binop;

    auto testContext = Context::getContextPtr();
    auto loc         = testContext;

    binop.getSupportedContext( loc );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // print warning if not available for test context

    SCAI_LOG_INFO( logger, "binaryOpScalar1Test<" << common::TypeTraits<ValueType>::id() << "> for " << *testContext << ", done on " << *loc )
    {
        // TODO: should it be possible to pass 0 elements? What should be the result?
        ValueType valuesValues[] =
        { 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
        HArray<ValueType> values( nValues, valuesValues );
        {
            WriteAccess<ValueType> wValues( values, loc );
            SCAI_CONTEXT_ACCESS( loc );
            binop[loc]( wValues.get(), wValues.get(), ValueType( 1 ), nValues, BinaryOp::DIVIDE, true );
        }

        ReadAccess<ValueType> rValues( values );

        for ( IndexType i = 0; i < nValues; i++ )
        {
            BOOST_CHECK_EQUAL( ValueType( 1 ) / valuesValues[i], rValues[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( countNonZerosTest )
{
    ContextPtr testContext = Context::getContextPtr();

    static LAMAKernel<SparseKernelTrait::countNonZeros<IndexType> > countNonZeros;

    ContextPtr loc = Context::getContextPtr( countNonZeros.validContext( testContext->getType() ) );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "countNonZeros for " << *testContext << " on " << *loc )

    // count valid array
    {
        const IndexType values[] = { 3, 0, 2, 0, 0, 2, 0, 4 };
        const IndexType n = sizeof( values ) / sizeof( IndexType );
        HArray<IndexType> sizes( n, values, testContext );
        ReadAccess<IndexType> rSizes( sizes, loc );
        SCAI_CONTEXT_ACCESS( loc );

        IndexType zero = 0;
        IndexType eps  = 0;
        BOOST_CHECK_EQUAL( IndexType( 4 ), countNonZeros[loc]( rSizes.get(), n, zero, eps ) );

        zero = 2;
        BOOST_CHECK_EQUAL( IndexType( 6 ), countNonZeros[loc]( rSizes.get(), n, zero, eps ) );
        eps  = 1;  //  1 - 3 are also considered to be zero
        BOOST_CHECK_EQUAL( IndexType( 5 ), countNonZeros[loc]( rSizes.get(), n, zero, eps ) );
    }

    // count empty array
    {
        HArray<IndexType> sizes;
        ReadAccess<IndexType> rSizes( sizes, loc );
        SCAI_CONTEXT_ACCESS( loc );
        IndexType zero = 1;
        IndexType eps  = 0;
        BOOST_CHECK_EQUAL( IndexType( 0 ), countNonZeros[loc]( rSizes.get(), 0, zero, eps ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( compressTest )
{
    ContextPtr testContext = Context::getContextPtr();

    static LAMAKernel<SparseKernelTrait::compress<IndexType, IndexType> > compress;

    ContextPtr loc = Context::getContextPtr( compress.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    const IndexType theValues[] = { 3, 0, 1, 0, 0, 1, 0, 4, 3, 0 };
    const IndexType theSparseIndexes[] = { 0, 2, 5, 7, 8 };
    const IndexType nDense = 10;
    const IndexType nSparse = 5;

    HArray<IndexType> denseArray( nDense, theValues, testContext );
    HArray<IndexType> sparseIndexes( nDense, IndexType( 0 ), testContext );

    {
        ReadAccess<IndexType> rArray( denseArray, loc );
        WriteAccess<IndexType> wSparseIndexes( sparseIndexes, loc );
        SCAI_CONTEXT_ACCESS( loc );
        IndexType cnt = compress[loc]( NULL, wSparseIndexes.get(), rArray.get(), nDense, 0, 0 );

        BOOST_REQUIRE_EQUAL( nSparse, cnt );
    }

    // sparseIndexes are sorted

    {
        ReadAccess<IndexType> rSparseIndexes( sparseIndexes );

        for ( IndexType i = 0; i < nSparse; ++i )
        {
            BOOST_CHECK_EQUAL( theSparseIndexes[i], rSparseIndexes[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setInversePermTest )
{
    ContextPtr testContext = Context::getContextPtr();

    static LAMAKernel<UtilKernelTrait::setInversePerm> setInversePerm;

    ContextPtr loc = Context::getContextPtr( setInversePerm.validContext( testContext->getType() ) );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    {
        IndexType valuesPerm[] = { 5, 0, 2, 3, 1, 4 };
        const IndexType nPerm = sizeof( valuesPerm ) / sizeof( IndexType );
        IndexType expectedPerm[] = { 1, 4, 2, 3, 5, 0 };
        const IndexType numRows = 6;
        HArray<IndexType> perm( nPerm, valuesPerm, testContext );
        HArray<IndexType> inversePerm;  // will be allocated/used on loc
        {
            ReadAccess<IndexType> rPerm( perm, loc );
            WriteOnlyAccess<IndexType> wInversePerm( inversePerm, loc, numRows );
            SCAI_CONTEXT_ACCESS( loc );
            setInversePerm[loc]( wInversePerm.get(), rPerm.get(), numRows );
        }
        ReadAccess<IndexType> rInversePerm( inversePerm );

        for ( IndexType i = 0; i < numRows; i++ )
        {
            BOOST_CHECK_EQUAL( expectedPerm[i], rInversePerm.get()[i] );
        }
    }
    {
        const IndexType numRows = 0;
        HArray<IndexType> perm( numRows );
        HArray<IndexType> inversePerm( numRows );
        {
            ReadAccess<IndexType> rPerm( perm, loc );
            WriteOnlyAccess<IndexType> wInversePerm( inversePerm, loc, numRows );
            SCAI_CONTEXT_ACCESS( loc );
            setInversePerm[loc]( wInversePerm.get(), rPerm.get(), numRows );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( sortIndexesTest )
{
    ContextPtr testContext = Context::getContextPtr();
    static LAMAKernel<UtilKernelTrait::sort<IndexType> > sort;
    ContextPtr loc = Context::getContextPtr( sort.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // just one simple example of a typical sort with permutation

    {
        bool ascending = false;

        IndexType valuesArr[]    = { 5, 2, 4, 4, 2, 7 };
        IndexType expectedArr[]  = { 7, 5, 4, 4, 2, 2 };
        IndexType expectedPerm[] = { 5, 0, 2, 3, 1, 4 };

        const IndexType n = sizeof( valuesArr ) / sizeof( IndexType );
        HArray<IndexType> array( n, valuesArr, testContext );

        HArray<IndexType> perm;

        {
            WriteOnlyAccess<IndexType> wPerm( perm, loc, n );
            WriteAccess<IndexType> wArr( array, loc );
            SCAI_CONTEXT_ACCESS( loc );
            sort[loc]( wPerm.get(), wArr.get(), wArr.get(), n, ascending );
        }

        ReadAccess<IndexType> rArr( array );
        ReadAccess<IndexType> rPerm( perm );

        for ( IndexType i = 0; i < n; i++ )
        {
            BOOST_CHECK_EQUAL( expectedArr[i], rArr.get()[i] );
            BOOST_CHECK_EQUAL( expectedPerm[i], rPerm.get()[i] );
        }
    }
    {
        bool ascending = true;

        IndexType valuesArr[]    = { 5, 2, 4, 4, 2, 7 };
        IndexType expectedArr[]  = { 2, 2, 4, 4, 5, 7 };

        const IndexType n = sizeof( valuesArr ) / sizeof( IndexType );

        HArray<IndexType> unsortedArray( n, valuesArr, testContext );
        HArray<IndexType> sortedArray;

        {
            WriteOnlyAccess<IndexType> wOut( sortedArray, loc, n );
            ReadAccess<IndexType> rIn( unsortedArray, loc );
            SCAI_CONTEXT_ACCESS( loc );
            sort[loc]( NULL, wOut.get(), rIn.get(), n, ascending );
        }

        ReadAccess<IndexType> rArr( sortedArray );

        for ( IndexType i = 0; i < n; i++ )
        {
            BOOST_CHECK_EQUAL( expectedArr[i], rArr[i] );
        }
    }
    {
        IndexType valuesArr[]    = { 5, 2, 4, 4, 2, 7 };
        IndexType valuesPerm[]   = { 0, 1, 2, 3, 4, 5 };
        IndexType expectedArr[]  = { 2, 2, 4, 4, 5, 7 };
        IndexType expectedPerm[] = { 1, 4, 2, 3, 0, 5 };

        const IndexType nPerm = sizeof( valuesPerm ) / sizeof( IndexType );
        const IndexType n = sizeof( valuesArr ) / sizeof( IndexType );
        const IndexType numRows = 6;
        HArray<IndexType> perm( nPerm, valuesPerm, testContext );
        HArray<IndexType> ilg( n, valuesArr, testContext );
        {
            WriteAccess<IndexType> wPerm( perm, loc );
            WriteAccess<IndexType> wArr( ilg, loc );
            SCAI_CONTEXT_ACCESS( loc );
            sort[loc]( wPerm.get(), wArr.get(), wArr.get(), numRows, true );
        }
        ReadAccess<IndexType> rArr( ilg );
        ReadAccess<IndexType> rPerm( perm );

        for ( IndexType i = 0; i < numRows; i++ )
        {
            BOOST_CHECK_EQUAL( expectedArr[i], rArr.get()[i] );
            BOOST_CHECK_EQUAL( expectedPerm[i], rPerm.get()[i] );
        }
    }
    {
        const IndexType n = 0;
        HArray<IndexType> perm( n );
        HArray<IndexType> ilg( n );
        {
            WriteOnlyAccess<IndexType> wPerm( perm, loc, n );
            WriteOnlyAccess<IndexType> wArr( ilg, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            sort[loc]( wPerm.get(), wArr.get(), wArr.get(), n, false );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( fullSortTest, ValueType, scai_array_test_types )
{
    typedef typename common::TypeTraits<ValueType>::RealType DefaultReal;

    ContextPtr testContext = Context::getContextPtr();

    static LAMAKernel<UtilKernelTrait::sort<DefaultReal> > sort;
    static LAMAKernel<UtilKernelTrait::isSorted<DefaultReal> > isSorted;

    ContextPtr loc = testContext;
    sort.getSupportedContext( loc, isSorted );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    const IndexType n = 1000;
    const bool ascending = true;

    HArray<DefaultReal> values( loc );
    {
        ContextPtr host = Context::getHostPtr();   // currently only available on host

        WriteOnlyAccess<DefaultReal> wValues( values, host, n );

        for ( IndexType i = 0; i < n; ++i )
        {
            wValues[i] = common::Math::random<ValueType>( 2 * n );
        }
    }

    HArray<IndexType> perm;

    WriteOnlyAccess<IndexType> wPerm( perm, loc, n );
    WriteAccess<DefaultReal> wValues( values, loc );

    SCAI_CONTEXT_ACCESS( loc );

    sort[loc]( wPerm.get(), wValues.get(), wValues.get(), n, ascending );
    bool valuesSorted = isSorted[loc]( wValues.get(), n, ascending ? CompareOp::LE : CompareOp::GE );

    BOOST_CHECK( valuesSorted );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( sortInPlaceTest )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = Context::getContextPtr();
    static LAMAKernel<UtilKernelTrait::sortInPlace<ValueType> > sortInPlace;
    ContextPtr loc = Context::getContextPtr( sortInPlace.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // just one simple example of a typical sort with permutation

    {
        bool ascending = true;

        IndexType indexesArr[]  = { 10, 5, 13, 1, 3, 0 };
        ValueType valuesArr[]   = {  1, 2, 3, 4, 5, 6  };

        IndexType expectedIndexesArr[]  = { 0, 1, 3, 5, 10, 13 };
        ValueType expectedValuesArr[]   = { 6, 4, 5, 2,  1,  3  };

        const IndexType n = sizeof( indexesArr ) / sizeof( IndexType );

        HArray<IndexType> indexes( n, indexesArr, testContext );
        HArray<ValueType> values( n, valuesArr, testContext );

        {
            WriteAccess<IndexType> wIndexes( indexes, loc );
            WriteAccess<ValueType> wValues( values, loc );
            SCAI_CONTEXT_ACCESS( loc );
            sortInPlace[loc]( wIndexes.get(), wValues.get(), n, ascending );
        }

        {
            ReadAccess<IndexType> rIndexes( indexes );
            ReadAccess<ValueType> rValues( values );

            for ( IndexType i = 0; i < n; i++ )
            {
                BOOST_CHECK_EQUAL( expectedIndexesArr[i], rIndexes[i] );
                BOOST_CHECK_EQUAL( expectedValuesArr[i], rValues[i] );
            }
        }

        {
            WriteAccess<IndexType> wIndexes( indexes, loc );
            WriteAccess<ValueType> wValues( values, loc );
            SCAI_CONTEXT_ACCESS( loc );
            sortInPlace[loc]( wIndexes.get(), wValues.get(), n, ascending );
        }

        {
            ReadAccess<IndexType> rIndexes( indexes );
            ReadAccess<ValueType> rValues( values );

            for ( IndexType i = 0; i < n; i++ )
            {
                BOOST_CHECK_EQUAL( expectedIndexesArr[i], rIndexes[i] );
                BOOST_CHECK_EQUAL( expectedValuesArr[i], rValues[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( sparseAddTest, ValueType, scai_array_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    static LAMAKernel<SparseKernelTrait::countAddSparse > countAddSparse;
    static LAMAKernel<SparseKernelTrait::addSparse<ValueType> > addSparse;

    ContextPtr loc = testContext;
    addSparse.getSupportedContext( loc, countAddSparse );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    const IndexType indexes1_values[]   = { 0,    2,   5,  7     };
    const IndexType indexes2_values[]   = {    1, 2,   5,      8 };
    const IndexType indexes_values[]    = { 0, 1, 2,   5,  7,  8 };

    const ValueType values1_values[]    = { 1,       2,     3,    4 };
    const ValueType values2_values[]    = {      5,  6,     7,         8 };
    const ValueType values_values[]     = { 1,  10, 14,    17,    4,  16 };

    IndexType nR = sizeof( indexes_values ) / sizeof( IndexType );
    IndexType n1 = sizeof( indexes1_values ) / sizeof( IndexType );
    IndexType n2 = sizeof( indexes2_values ) / sizeof( IndexType );

    SCAI_ASSERT_EQ_ERROR( n1, sizeof( values1_values ) / sizeof( ValueType ), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( n2, sizeof( values2_values ) / sizeof( ValueType ), "size mismatch" )

    HArray<IndexType> indexes1( n1, indexes1_values, testContext );
    HArray<IndexType> indexes2( n2, indexes2_values, testContext );

    HArray<ValueType> values1( n1, values1_values, testContext );
    HArray<ValueType> values2( n2, values2_values, testContext );

    HArray<IndexType> indexes;
    HArray<ValueType> values;

    ValueType alpha = 1;
    ValueType beta  = 2;

    ValueType zero = 0;

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIndexes1( indexes1, loc );
        ReadAccess<IndexType> rIndexes2( indexes2, loc );

        IndexType n = countAddSparse[loc]( rIndexes1.get(), n1, rIndexes2.get(), n2 );

        BOOST_CHECK_EQUAL( n, nR );

        ReadAccess<ValueType> rValues1( values1, loc );
        ReadAccess<ValueType> rValues2( values2, loc );

        WriteOnlyAccess<IndexType> wIndexes( indexes, loc, n );
        WriteOnlyAccess<ValueType> wValues( values, loc, n );

        IndexType nc = addSparse[loc]( wIndexes.get(), wValues.get(), 
                                       rIndexes1.get(), rValues1.get(), zero, n1, alpha,
                                       rIndexes2.get(), rValues2.get(), zero, n2, beta );

        BOOST_CHECK_EQUAL( n, nc );
    }

    BOOST_CHECK_EQUAL( values.size(), indexes.size() );

    {
        ReadAccess<ValueType> rValues( values );
        ReadAccess<IndexType> rIndexes( indexes );

        for ( IndexType i = 0; i < values.size(); ++i )
        {
            BOOST_CHECK_EQUAL( indexes_values[i], rIndexes[i] );
            BOOST_CHECK_EQUAL( values_values[i], rValues[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */


BOOST_AUTO_TEST_SUITE_END()

