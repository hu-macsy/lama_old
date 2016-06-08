/**
 * @file HArrayUtilsTest.cpp
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
 * @brief Tests for the class HArrayUtils
 * @author Thomas Brandes
 * @date 22.01.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/utilskernel/ReductionOp.hpp>

#include <scai/utilskernel/test/TestMacros.hpp>
#include <scai/utilskernel/test/HArrays.hpp>

#include <scai/common/Math.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/exception/Exception.hpp>

#include <typeinfo>

using namespace scai;
using namespace scai::utilskernel;
using namespace scai::hmemo;
using namespace scai::common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HArrayUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HArrayUtilsTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( FactoryTest )
{
    HArrays allArrays;    // is created by factory

    size_t nTypes = SCAI_COMMON_COUNT_NARG( SCAI_ARITHMETIC_ARRAY_HOST );

    SCAI_LOG_INFO( logger, "Test all arrys of factory to be zero, #arrays = " << allArrays.size() )

    BOOST_CHECK_EQUAL( nTypes, allArrays.size() );

    for ( size_t i = 0; i < allArrays.size(); ++i )
    {
        _HArray& array = *allArrays[i];

        BOOST_CHECK_EQUAL( 0, array.size() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( UntypedTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr ctx  = Context::getContextPtr();

    HArrays allArrays1;    // vector with HArray of each type
    HArrays allArrays2;    // another vector to get each pair combination

    const IndexType permVals[] = { 1, 0, 3, 2 };

    const IndexType n = sizeof( permVals ) / sizeof( IndexType );

    LArray<IndexType> perm( n, permVals );

    LArray<IndexType> order;
    HArrayUtils::setOrder( order, n );    // 0, 1, 2, .., n-1

    for ( size_t i1 = 0; i1 < allArrays1.size(); ++i1 )
    {
        _HArray& array1 = *allArrays1[i1];

        HArrayUtils::assign( array1, order );

        for ( size_t i2 = 0; i2 < allArrays2.size(); ++i2 )
        {
            _HArray& array2 = *allArrays2[i2];

            // create array with same type, context

            common::unique_ptr<_HArray> tmp( array2.copy() );

            // array2 = array1[ perm ], tmp[ perm ] = array1

            HArrayUtils::assignGather( array2, array1, perm, ctx );

            BOOST_CHECK_EQUAL( array2.size(), perm.size() );

            tmp->resize( n ); // no init required as all values are set

            HArrayUtils::assignScatter( *tmp, perm, array1, ctx );

            // as perm is its inverse, tmp and array2 should be the same

            for ( IndexType k = 0; k < n; ++k )
            {
                ValueType val1 = HArrayUtils::getVal<ValueType>( array2, k );
                ValueType val2 = HArrayUtils::getVal<ValueType>( *tmp, k );
                BOOST_CHECK_EQUAL( val1, val2 );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetScalarTest, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( SetValueTest, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( GatherTest, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( ScatterTest, ValueType, scai_arithmetic_test_types )
{
    ValueType sourceVals[] = { 3, 1, 4, 2 };
    IndexType indexVals[]  = { 0, 2, 1, 3 };

    const IndexType N = sizeof( sourceVals ) / sizeof( ValueType );
    const IndexType N1 = sizeof( indexVals ) / sizeof( IndexType );

    BOOST_REQUIRE_EQUAL( N, N1 );

    const IndexType M = 5; // size of target arry, can be larger

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_REQUIRE( indexVals[i] < M );
    }

    // target = source[ indexes ]

    LArray<ValueType> source( N, sourceVals );
    LArray<IndexType> indexes( N, indexVals );

    BOOST_CHECK_THROW (
    {
        LArray<ValueType> target;
        HArrayUtils::scatter( target, indexes, source );

    }, Exception );

    LArray<ValueType> target( M );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( randomTest, ValueType, array_types )
{
    ContextPtr loc = Context::getContextPtr();

    LArray<ValueType> array( loc );

    const IndexType n = 100;

    float fillRate = 1.0f;

    HArrayUtils::setRandom( array, n, fillRate, loc );

    BOOST_CHECK_EQUAL( array.size(), n );

    ValueType sum = array.sum();

    if ( typeid( ValueType ).name() == typeid( IndexType ).name() )
    {
        // no check for integer types
    }
    else
    {
        typedef typename TypeTraits<ValueType>::AbsType AbsType;

        AbsType asum = Math::abs( sum );

        // random numbers are between -1.0 and 1.0, so should sum up approximately to 0

        BOOST_CHECK( asum / AbsType( n ) < AbsType( 0.2 ) );
    }
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

BOOST_AUTO_TEST_CASE( setOrderTest )
{
    ContextPtr loc = Context::getContextPtr();

    const IndexType n = 10;

    LArray<IndexType> array;

    HArrayUtils::setOrder( array, n, loc );

    BOOST_CHECK_EQUAL( array.size(), n );

    for ( IndexType i = 0; i < n; ++i )
    {
        IndexType elem = array[i];
        BOOST_CHECK_EQUAL( i, elem );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( axpyTest, ValueType, scai_arithmetic_test_types )
{
    ContextPtr loc = Context::getContextPtr();

    ValueType xVals[] = { 3, 1, 4, 2 };
    ValueType yVals[] = { 2, -1, -1, -5 };

    const IndexType nx = sizeof( xVals ) / sizeof( ValueType );
    const IndexType ny = sizeof( yVals ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( nx, ny );

    LArray<ValueType> x( nx, xVals, loc );
    LArray<ValueType> y( ny, yVals, loc );

    ValueType alpha = 2.0;

    HArrayUtils::axpy( y, alpha, x, loc ); // y += alpha * x

    for ( IndexType i = 0; i < ny; ++i )
    {
        ValueType elem = y[i];
        BOOST_CHECK_EQUAL( yVals[i] + alpha * xVals[i], elem );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( arrayPlusArrayTest, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( arrayPlusArrayAliasTest, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( sparseTest, ValueType, scai_arithmetic_test_types )
{
    ContextPtr loc = Context::getContextPtr();

    ValueType denseValues[] = { 3, 0, 0, 1, 4, 0, 2, 7, 0, -1, -3 };
    const IndexType n = sizeof( denseValues ) / sizeof( ValueType );

    LArray<ValueType> denseArray( n, denseValues, loc );

    /* Idea: convert to sparse and back to dense, must contain original data */

    LArray<ValueType> sparseArray( loc );
    LArray<IndexType> sparseIndexes( loc );

    HArrayUtils::buildSparseArray( sparseArray, sparseIndexes, denseArray, loc );

    BOOST_CHECK_EQUAL( sparseArray.size(), 7 );
    BOOST_REQUIRE_EQUAL( sparseArray.size(), sparseIndexes.size() );

    // The spare indexes must be sorted, ascending = true

    BOOST_CHECK( HArrayUtils::isSorted( sparseIndexes, true, loc ) );

    denseArray.purge();  // will also reset size

    HArrayUtils::buildDenseArray( denseArray, n, sparseArray, sparseIndexes, loc );

    BOOST_REQUIRE_EQUAL( denseArray.size(), n );

    for ( IndexType i = 0; i < n; ++i )
    {
        ValueType v = denseArray[i];
        BOOST_CHECK_EQUAL( v, denseValues[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
