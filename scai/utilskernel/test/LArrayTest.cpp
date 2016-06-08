/**
 * @file LArrayTest.cpp
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
 * @endlicense
 *
 * @brief Tests for the class LArray
 * @author Thomas Brandes
 * @date 22.02.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/utilskernel/LArray.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace utilskernel;
using namespace hmemo;
using namespace common;

extern ContextPtr testContext;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LArrayTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.LArrayTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( indexTest )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "indexTest on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    LArray<IndexType> array( testContext );

    const IndexType N = 10;

    array.resize( N );

    for ( IndexType i = 0; i < N; ++i )
    {
        array[i] = i;
    }

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_CHECK_EQUAL( i, array[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assignTest )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "assignTest on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    LArray<IndexType> array( testContext );

    const IndexType N = 10;
    const IndexType val = 3;

    array.resize( N );

    array = val;

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_CHECK_EQUAL( val, array[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "constructorTest on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    const IndexType N = 10;
    const IndexType val = 3;

    LArray<IndexType> array( N, val, testContext );

    BOOST_CHECK( array.isValid( testContext ) );

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_CHECK_EQUAL( val, array[i] );
    }

    const IndexType myVals[N] = { 1, 5, 9, 4, 6, 3, 7, 8, 0, 2 };

    LArray<IndexType> array1( N, myVals, testContext );

    BOOST_CHECK( array.isValid( testContext ) );

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_CHECK_EQUAL( myVals[i], array1[i] );
    }
}

/* --------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_ARITHMETIC_ARRAY_CUDA> ArrayRedTypes;

// ToDo: introduce a predicate in COMMON to check if a certain type is supported on a context

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( reductionTest, ValueType, ArrayRedTypes )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "reductionTest on " << *testContext )

    // ToDo: example with complex numbers

    const ValueType myVals[] = { 9, 5, 1, 4, 6, 3, 7, 8, 2, 0 };
    const IndexType N = sizeof( myVals ) / sizeof( ValueType );

    ValueType expectedMin = TypeTraits<ValueType>::getMax();
    ValueType expectedMax = TypeTraits<ValueType>::getMin();;
    ValueType expectedSum = 0;

    for ( IndexType i = 0; i < N; ++i )
    {
        expectedMin = Math::min( expectedMin, myVals[i] );
        expectedMax = Math::max( expectedMax, myVals[i] );
        expectedSum += myVals[i];
    }

    LArray<ValueType> array( N, myVals, testContext );

    // Constructor should have provided a valid copy on testContext

    BOOST_CHECK( array.isValid( testContext ) );

    // reduction ops will be executed on testContext

    BOOST_CHECK_EQUAL( expectedMin, array.min() );
    BOOST_CHECK_EQUAL( expectedMax, array.max() );
    BOOST_CHECK_EQUAL( expectedSum, array.sum() );
}

/* --------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_ARITHMETIC_CUDA> ArithmeticRedTypes;

// ToDo: introduce a predicate in COMMON to check if a certain type is supported on a context

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( normTest, ValueType, ArithmeticRedTypes )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "normTest on " << *testContext )

    // ToDo: example with complex numbers

    const IndexType N = 13;

    scoped_array<ValueType> myVals( new ValueType[N] );

    for ( IndexType i = 0; i < N; ++i )
    {
        // random numbers between -1.0 and 1.0, real and imag part for complex

        Math::random( myVals[i] );
    }

    typedef typename TypeTraits<ValueType>::AbsType AbsType;

    AbsType expectedL1Norm = 0;
    AbsType expectedL2Norm = 0;
    AbsType expectedMaxNorm = 0;

    for ( IndexType i = 0; i < N; ++i )
    {
        expectedL1Norm += Math::abs( Math::real( myVals[i] ) );
        expectedL1Norm += Math::abs( Math::imag( myVals[i] ) );
        expectedL2Norm += Math::real( myVals[i] * Math::conj( myVals[i] ) );
        expectedMaxNorm = Math::max( expectedMaxNorm, Math::abs( myVals[i] ) );
    }

    expectedL2Norm = Math::sqrt( expectedL2Norm );

    LArray<ValueType> array( N, myVals.get(), testContext );

    // Constructor should have provided a valid copy on testContext

    BOOST_CHECK( array.isValid( testContext ) );

    // reduction ops will be executed on testContext

    BOOST_CHECK_CLOSE( expectedL1Norm, AbsType( array.l1Norm() ), 0.1  );
    BOOST_CHECK_CLOSE( expectedL2Norm, AbsType( array.l2Norm() ), 0.1 );
    BOOST_CHECK_CLOSE( expectedMaxNorm, AbsType( array.maxNorm() ), 0.1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( binReductionTest, ValueType, ArithmeticRedTypes )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "binReductionTest<" << TypeTraits<ValueType>::id() << " on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    const ValueType myVals1[] = { 9, 5, 1, 4, 6, 3, 7, 8, 2, 0 };
    const ValueType myVals2[] = { 9, 5, 1, 3, 6, 3, 7, 8, 2, 0 };

    const IndexType N = sizeof( myVals1 ) / sizeof( ValueType );

    ValueType expectedDotProduct = 0;
    ValueType expectedMaxDiffNorm = 0;

    for ( IndexType i = 0; i < N; ++i )
    {
        expectedDotProduct += myVals1[i] * Math::conj( myVals2[i] );
        ValueType absDiff = Math::abs( myVals2[i] - myVals1[i] );
        expectedMaxDiffNorm = Math::max( expectedMaxDiffNorm, absDiff );
    }

    LArray<ValueType> array1( N, myVals1, testContext );
    LArray<ValueType> array2( N, myVals2, testContext );

    BOOST_CHECK( array1.isValid( testContext ) );
    BOOST_CHECK( array2.isValid( testContext ) );

    BOOST_CHECK_EQUAL( expectedMaxDiffNorm, array1.maxDiffNorm( array2 ) );
    BOOST_CHECK_EQUAL( expectedDotProduct, array1.dotProduct( array2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( invertTest, ValueType, ArithmeticRedTypes )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "invertTest<" << TypeTraits<ValueType>::id() << " on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    const ValueType myVals[] = { 9, 5, 1, 4, 6, 3, 7, 8, 2, 3 };

    const IndexType N = sizeof( myVals ) / sizeof( ValueType );

    LArray<ValueType> array( N, myVals, testContext );

    // make sure not to divide by zero

    BOOST_CHECK( Math::real( array.maxNorm() ) > 0 );

    array.invert();

    for ( IndexType i = 0; i < N; ++i )
    {
        typedef typename TypeTraits<ValueType>::AbsType AbsType;

        ValueType x1 = 1 / myVals[i];
        ValueType x2 = array[i];

        BOOST_CHECK_CLOSE( AbsType( x1 ), AbsType( x2 ), 0.01 );
    }

    array.invert();

    for ( IndexType i = 0; i < N; ++i )
    {
        typedef typename TypeTraits<ValueType>::AbsType AbsType;

        AbsType x1 =  myVals[i];
        ValueType x2 = array[i];
        BOOST_CHECK_CLOSE( x1, AbsType( x2 ), 0.01 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( conjTest, ValueType, ArithmeticRedTypes )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "conjTest<" << TypeTraits<ValueType>::id() << " on " << *testContext )

    const IndexType N = 16;

    LArray<ValueType> array;
    HArrayUtils::setRandom( array, N, 1.0f, testContext );

    LArray<ValueType> conjArray( array );
    conjArray.conj();   // build in place


    if ( isComplex( TypeTraits<ValueType>::stype )  )
    {
        BOOST_CHECK( Math::real( array.maxDiffNorm( conjArray ) ) > 0 );
        conjArray.conj();
        BOOST_CHECK_EQUAL( 0, Math::real( array.maxDiffNorm( conjArray ) ) );
    }
    else
    {
        // not complex: both arrays should be the same

        BOOST_CHECK_EQUAL( 0, Math::real( array.maxDiffNorm( conjArray ) ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignOperatorTest, ValueType, ArithmeticRedTypes )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "assignOperatorTest<" << TypeTraits<ValueType>::id() << " on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    const ValueType myVals[]  = { 9, 5, 1, 4, 6, 3, 7, 8, 2, 0 };
    const ValueType myVals1[] = { 1, -1, 2, -2, 3, -3, 2, -1, -2, 1 };

    const IndexType N = sizeof( myVals ) / sizeof( ValueType );

    LArray<ValueType> array( N, myVals, testContext );
    LArray<ValueType> other( N, myVals1, testContext );

    array /= ValueType( 2 );
    array *= ValueType( 2 );
    array += ValueType( 2 );
    array -= ValueType( 2 );
    array /= other;
    array *= other;
    array += other;
    array -= other;

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_CHECK_CLOSE( Math::real( myVals[i] ), Math::real( array[i] ), 0.01 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
