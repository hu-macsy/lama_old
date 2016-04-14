/**
 * @file HArrayUtilsTest.cpp
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
 * @brief Tests for the class HArrayUtils
 * @author: Thomas Brandes
 * @date 22.01.2016
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/utilskernel/test/TestMacros.hpp>
#include <scai/utilskernel/test/HArrays.hpp>
#include <scai/common/ReductionOp.hpp>
#include <scai/common/exception/Exception.hpp>

using namespace scai;
using namespace scai::utilskernel;
using namespace scai::hmemo;
using namespace scai::common;

typedef boost::mpl::list<SCAI_ARITHMETIC_HOST> test_types;

// typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HArrayUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HArrayUtilsTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( FactoryTest )
{
    HArrays allArrays;    // is created by factory

    size_t nTypes = SCAI_ARRAY_TYPE_CNT;

    SCAI_LOG_INFO( logger, "Test all arrys of factory to be zero, #arrays = " << allArrays.size() )

    BOOST_CHECK_EQUAL( nTypes, allArrays.size() );

    for ( size_t i = 0; i < allArrays.size(); ++i )
    {
        _HArray& array = *allArrays[i];

        BOOST_CHECK_EQUAL( 0, array.size() ); 
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetScalarTest, ValueType, test_types )
{
    const IndexType N = 10;

    const ValueType a = 1;
    const ValueType b = 2;
    const ValueType c = 3;

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    HArray<ValueType> array( N );

    HArrayUtils::setScalar( array, a, reduction::COPY, ctx );  // array = a
    HArrayUtils::setScalar( array, b, reduction::ADD, ctx  );  // array += b
    HArrayUtils::setScalar( array, c, reduction::MULT, ctx  ); // array *= c

    ValueType expectedValue = ( a + b ) * c;

    {
        ReadAccess<ValueType> read( array, host );

        for ( IndexType i = 0; i < N; ++i )
        {
             BOOST_CHECK_EQUAL( expectedValue, read[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetValueTest, ValueType, test_types )
{
    const IndexType N = 10;
    const IndexType k = 3;    // one value in range 0 .. N-1

    const ValueType a = 1;
    const ValueType b = 2;

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    HArray<ValueType> array( N );

    HArrayUtils::setScalar( array, a, reduction::COPY, ctx );  // array = a

    HArrayUtils::setVal( array, k, b );  // array[k] = b

    {
        ReadAccess<ValueType> read( array, host );

        for ( IndexType i = 0; i < N; ++i )
        {
            if ( i == k )
            {
                BOOST_CHECK_EQUAL( b, read[i] );
            }
            else
            {
                BOOST_CHECK_EQUAL( a, read[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( GatherTest, ValueType, test_types )
{
    ValueType sourceVals[] = { 3, 1, 4, 2 };
    IndexType indexVals[]  = { 0, 2, 1, 2, 1, 3 };

    const IndexType M = sizeof( sourceVals ) / sizeof( ValueType );
    const IndexType N = sizeof( indexVals ) / sizeof( IndexType );

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_REQUIRE( indexVals[i] < M );
    }

    // target = source[ indexes ]

    LArray<ValueType> source( M, sourceVals );
    LArray<IndexType> indexes( N, indexVals );
    LArray<ValueType> target;

    BOOST_CHECK( HArrayUtils::validIndexes( indexes, M ) );
    BOOST_CHECK( !HArrayUtils::validIndexes( indexes, 1 ) );

    HArrayUtils::gather( target, source, indexes );

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType x1 = source[indexes[i]];
        ValueType x2 = target[i];
        BOOST_CHECK_EQUAL( x1, x2 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ScatterTest, ValueType, test_types )
{
    ValueType sourceVals[] = { 3, 1, 4, 2 };
    IndexType indexVals[]  = { 0, 2, 1, 3 };

    const IndexType M = sizeof( sourceVals ) / sizeof( ValueType );
    const IndexType N = sizeof( indexVals ) / sizeof( IndexType );

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_REQUIRE( indexVals[i] < M );
    }

    // target = source[ indexes ]

    LArray<ValueType> source( M, sourceVals );
    LArray<IndexType> indexes( N, indexVals );
    LArray<ValueType> target;

    BOOST_CHECK( HArrayUtils::validIndexes( indexes, M ) );

    HArrayUtils::scatter( target, indexes, source );

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType x1 = target[indexes[i]];
        ValueType x2 = source[i];

        BOOST_CHECK_EQUAL( x1, x2 );
    }
}

/* --------------------------------------------------------------------- */

typedef boost::mpl::list<IndexType, SCAI_ARITHMETIC_HOST> array_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( scanTest, ValueType, array_types )
{
    ContextPtr loc = Context::getContextPtr();

    ValueType vals[]     = { 3, 1, 4, 2 };

    const IndexType n = sizeof( vals ) / sizeof( ValueType );

    scoped_array<ValueType> scans( new ValueType[n+1] );

    scans[0] = 0;

    for ( IndexType i = 0; i < n; ++i )
    {
        scans[i+1] = scans[i] + vals[i];
    }

    LArray<ValueType> array;
    array.reserve( loc, n+1 );
    array.init( vals, n );
    LArray<ValueType> correct( n+1, scans.get(), loc );

    ValueType total = HArrayUtils::scan( array );
    ValueType lastVal = array[n];

    BOOST_CHECK_EQUAL( array.size(), n + 1 );
    BOOST_CHECK_EQUAL( array.maxDiffNorm( correct ), 0 );
    BOOST_CHECK_EQUAL( total, lastVal );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sortTest, ValueType, array_types )
{
    ContextPtr loc = Context::getContextPtr();

    ValueType vals[] = { 13, 1, 14, -2 };

    const IndexType n = sizeof( vals ) / sizeof( ValueType );

    LArray<ValueType> array( n, vals, loc );
    LArray<IndexType> perm;

    HArrayUtils::sort( array, perm );

    BOOST_CHECK_EQUAL( perm.size(), n );

    LArray<ValueType> origArray( n, vals, loc );
    LArray<ValueType> array1;
  
    HArrayUtils::gather( array1, origArray, perm );

    BOOST_CHECK_EQUAL( array.maxDiffNorm( array1 ), 0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( arrayPlusArrayTest, ValueType, test_types )
{
    ContextPtr loc = Context::getContextPtr();

    ValueType sourceVals1[] = { 3, 1, 4, 2 };
    ValueType sourceVals2[] = { 2, -1, -1, -5 };

    const IndexType n1 = sizeof( sourceVals1 ) / sizeof( ValueType );
    const IndexType n2 = sizeof( sourceVals2 ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n1, n2 );

    LArray<ValueType> x1( n1, sourceVals1, loc );
    LArray<ValueType> x2( n2, sourceVals2, loc );
    LArray<ValueType> xf( n2-1, sourceVals2, loc );  // wrong sized array

    LArray<ValueType> target( loc );

    ValueType factors[] = { -1, 0, 1, 2 };
    int NCASE = sizeof( factors ) / sizeof( ValueType );

    for ( int k1 = 0; k1 < NCASE; ++k1 )
    {
        ValueType alpha = factors[k1];

        for ( int k2 = 0; k2 < NCASE; ++k2 )
        {
            ValueType beta = factors[k2];

            SCAI_LOG_DEBUG( logger, "target = " << alpha << " * x1 + " << beta << " * x2" )

            target.purge();

            if ( beta == ValueType( 0 ) || alpha == ValueType( 0 ) )
            {
                // as one factor is zero, sizes must not match 

                HArrayUtils::arrayPlusArray( target, alpha , x1, beta, xf, loc );

                if ( alpha != ValueType( 0 ) )
                {
                    BOOST_CHECK_EQUAL( x1.size(), target.size() );
                }
                else if ( beta != ValueType( 0 ) )
                {
                    BOOST_CHECK_EQUAL( xf.size(), target.size() );
                }
            }
            else
            {
                // alpha, beta != 0, so sizes must match

                BOOST_CHECK_THROW (
                {
                     HArrayUtils::arrayPlusArray( target, alpha , x1, beta, xf, loc );
                }, Exception );
            }

            HArrayUtils::arrayPlusArray( target, alpha, x1, beta, x2, loc );

            if ( alpha == ValueType( 0 ) && beta == ValueType( 0 ) )
            {
                // target should be unchanged, size was 0 due to purge

                BOOST_CHECK_EQUAL( 0, target.size() );
                continue;
            }

            BOOST_CHECK_EQUAL( x1.size(), target.size() );

            for ( IndexType i = 0; i < n1; ++i )
            {
                ValueType v = target[i];
                BOOST_CHECK_EQUAL( v, alpha * sourceVals1[i] + beta * sourceVals2[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( arrayPlusArrayAliasTest, ValueType, test_types )
{
    ContextPtr loc = Context::getContextPtr();

    ValueType sourceVals1[] = { 3, 1, 4, 2 };
    ValueType sourceVals2[] = { 2, -1, -1, -5 };

    const IndexType n = sizeof( sourceVals1 ) / sizeof( ValueType );
    const IndexType n2 = sizeof( sourceVals2 ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n, n2 );

    ValueType factors[] = { -1, 0, 1, 2 };
    int NCASE = sizeof( factors ) / sizeof( ValueType );

    for ( int k1 = 0; k1 < NCASE; ++k1 )
    {
        ValueType alpha = factors[k1];

        for ( int k2 = 0; k2 < NCASE; ++k2 )
        {
            ValueType beta = factors[k2];

            {
                LArray<ValueType> x1( n, sourceVals1, loc );
                LArray<ValueType> x2( n, sourceVals2, loc );

                SCAI_LOG_DEBUG( logger, "x1 = " << alpha << " * x1 + " << beta << " * x2" )

                // target array aliased to x1

                HArrayUtils::arrayPlusArray( x1, alpha, x1, beta, x2, loc );

                for ( IndexType i = 0; i < n; ++i )
                {
                    ValueType v = x1[i];
                    BOOST_CHECK_EQUAL( v, alpha * sourceVals1[i] + beta * sourceVals2[i] );
                }
            }

            {
                LArray<ValueType> x1( n, sourceVals1, loc );
                LArray<ValueType> x2( n, sourceVals2, loc );

                SCAI_LOG_DEBUG( logger, "x2 = " << alpha << " * x1 + " << beta << " * x2" )

                // target array aliased to x2

                HArrayUtils::arrayPlusArray( x2, alpha, x1, beta, x2, loc );

                for ( IndexType i = 0; i < n; ++i )
                {
                    ValueType v = x2[i];
                    BOOST_CHECK_EQUAL( v, alpha * sourceVals1[i] + beta * sourceVals2[i] );
                }
            }

            {
                LArray<ValueType> x1( n, sourceVals1, loc );

                SCAI_LOG_DEBUG( logger, "x1 = " << alpha << " * x1 + " << beta << " * x1" )

                // target array aliased to x1 and x2

                HArrayUtils::arrayPlusArray( x1, alpha, x1, beta, x1, loc );

                for ( IndexType i = 0; i < n; ++i )
                {
                    ValueType v = x1[i];
                    BOOST_CHECK_EQUAL( v, alpha * sourceVals1[i] + beta * sourceVals1[i] );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
