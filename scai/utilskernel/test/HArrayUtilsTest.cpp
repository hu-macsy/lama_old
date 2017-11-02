/**
 * @file HArrayUtilsTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Tests for the class HArrayUtils
 * @author Thomas Brandes
 * @date 22.01.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/utilskernel/test/TestMacros.hpp>
#include <scai/utilskernel/test/HArrays.hpp>

#include <scai/common/Math.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/exception/Exception.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/mepr/ScalarTypeHelper.hpp>

#include <typeinfo>
#include <memory>

using std::unique_ptr;

using namespace scai;
using namespace scai::utilskernel;
using namespace scai::hmemo;
using namespace scai::common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HArrayUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HArrayUtilsTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( factoryTest )
{
    HArrays allArrays;    // is created by factory

    size_t nTypes = SCAI_COMMON_COUNT_NARG( SCAI_ALL_TYPES );

    SCAI_LOG_INFO( logger, "Test all arrys of factory to be zero, #arrays = " << allArrays.size() )
    BOOST_CHECK_EQUAL( nTypes, allArrays.size() );

    for ( size_t i = 0; i < allArrays.size(); ++i )
    {
        _HArray& array = *allArrays[i];
        BOOST_CHECK_EQUAL( IndexType( 0 ), array.size() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( untypedTest )
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

        // Note: HArrayUtils are only available for SCAI_ARRAY_TYPES_HOST

        if ( !common::mepr::ScalarTypeHelper<SCAI_ARRAY_TYPES_HOST_LIST>::contains( array1.getValueType() ) )
        {
            BOOST_CHECK_THROW(
            {
                HArrayUtils::assign( array1, order );
            }, common::Exception );

            continue;
        }
        else
        {
            HArrayUtils::assign( array1, order );
        }

        for ( size_t i2 = 0; i2 < allArrays2.size(); ++i2 )
        {
            _HArray& array2 = *allArrays2[i2];

            if ( !common::mepr::ScalarTypeHelper<SCAI_ARRAY_TYPES_HOST_LIST>::contains( array2.getValueType() ) )
            {
                BOOST_CHECK_THROW(
                {
                    HArrayUtils::gather( array2, array1, perm, BinaryOp::COPY, ctx );
                }, common::Exception );

                continue;
            }

            // create array with same type, context
            std::unique_ptr<_HArray> tmp( array2.copy() );
            // array2 = array1[ perm ], tmp[ perm ] = array1
            HArrayUtils::gather( array2, array1, perm, BinaryOp::COPY, ctx );
            BOOST_CHECK_EQUAL( array2.size(), perm.size() );
            tmp->resize( n ); // no init required as all values are set
            bool unique = true;  // perm has unique indexes
            HArrayUtils::scatter( *tmp, perm, unique, array1, BinaryOp::COPY, ctx );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( setScalarTest, ValueType, scai_numeric_test_types )
{
    const IndexType N = 10;
    const ValueType a = 1;
    const ValueType b = 2;
    const ValueType c = 3;
    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();
    HArray<ValueType> array( N );
    HArrayUtils::setScalar( array, a, BinaryOp::COPY, ctx );  // array = a
    HArrayUtils::setScalar( array, b, BinaryOp::ADD, ctx  );  // array += b
    HArrayUtils::setScalar( array, c, BinaryOp::MULT, ctx  ); // array *= c
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

BOOST_AUTO_TEST_CASE_TEMPLATE( setValueTest, ValueType, scai_numeric_test_types )
{
    const IndexType N = 10;
    const IndexType k = 3;    // one value in range 0 .. N-1
    const ValueType a = 1;
    const ValueType b = 2;
    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();
    HArray<ValueType> array( N );
    HArrayUtils::setScalar( array, a, BinaryOp::COPY, ctx );  // array = a
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

BOOST_AUTO_TEST_CASE_TEMPLATE( UnaryOpOpTest, ValueType, scai_array_test_types )
{
    // check of all UnaryOp array operations

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    const ValueType values[] = { 1, 2, 3, 4  };
    const IndexType n = sizeof( values ) / sizeof( ValueType );

    for ( IndexType i = 0; i < UnaryOp::MAX_UNARY_OP; ++i )
    {
        UnaryOp op = UnaryOp( i );

        HArray<ValueType> array( ctx );

        array.setRawData( n, values );

        SCAI_LOG_DEBUG( logger, "test UnaryOp op " << op << " for " << array )

        if ( ! isUnarySupported<ValueType>( op ) )
        {
            BOOST_CHECK_THROW(
            {
                HArrayUtils::UnaryOpOp( array, array, op, ctx );
            }, Exception );

            continue;  // not all operations are supported for IndexType
        }

        HArrayUtils::UnaryOpOp( array, array, op, ctx );

        ReadAccess<ValueType> read( array, host );

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType res = applyUnary( op, values[i] );

            typedef typename TypeTraits<ValueType>::AbsType AbsType;

            AbsType diff = common::Math::abs( read[i] - res  );

            // might happen that result on other devices are not exactly the same, give message

            if ( diff != AbsType( 0 ) )
            {
                BOOST_TEST_MESSAGE( "Result " << read[i] << " on " << *ctx
                                    << " is different from expected result " << res
                                    << " = " << op << " " << values[i]
                                    << ", diff = " << diff );
            }

            // but they must be close, otherwise fail

            if ( common::isNumeric( TypeTraits<ValueType>::stype ) && AbsType( res ) > 1 )
            {
                // large numbers due to EXP function, so take relative error

                diff /= AbsType( res );
            }

            BOOST_CHECK( diff <= TypeTraits<AbsType>::small() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( binaryOpTest, ValueType, scai_numeric_test_types )
{
    // check of all UnaryOp array operations

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    const ValueType values1[] = { 1.0, 1.2, -1.3, -1.0 };
    const ValueType values2[] = { 0.5, -0.7, 0.3, -1.3 };

    const IndexType n = sizeof( values1 ) / sizeof( ValueType );

    for ( IndexType i = 0; i < BinaryOp::MAX_BINARY_OP; ++i )
    {
        BinaryOp op = BinaryOp( i );

        if ( op == BinaryOp::COPY )
        {
            continue;    // not implemented
        }

        HArray<ValueType> array1( ctx );
        HArray<ValueType> array2( ctx );
        HArray<ValueType> array3( ctx );

        array1.setRawData( n, values1 );
        array2.setRawData( n, values2 );

        if ( op == BinaryOp::POW )
        {
            // first argument must not be negative, build ABS

            HArrayUtils::UnaryOpOp( array1, array1, UnaryOp::ABS, ctx );
        }

        HArrayUtils::binaryOp( array3, array1, array2, op, ctx );

        BOOST_REQUIRE_EQUAL( n, array3.size() );

        ReadAccess<ValueType> read( array3, host );  // read result array

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType x1  = values1[i];

            if ( op == BinaryOp::POW )
            {
                x1 = common::Math::abs( x1 );
            }

            ValueType res = applyBinary( x1, op, values2[i] );

            typedef typename TypeTraits<ValueType>::AbsType AbsType;

            AbsType diff = common::Math::abs( read[i] - res  );

            // might happen that result on other devices are not exactly the same, give message

            if ( diff <= TypeTraits<AbsType>::small() )
            {
                BOOST_TEST_MESSAGE( "Result " << read[i] << " on " << *ctx
                                    << " is different from expected result " << res
                                    << " = " << x1 << " " << op << " " << values2[i]
                                    << ", diff = " << diff );
            }

            // but they must be close, otherwise fail

            BOOST_CHECK( diff <= TypeTraits<AbsType>::small() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( binaryOpIndexTypeTest )
{
    // check of all UnaryOp array operations

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    const IndexType values1[] = { 1, 3, 7, 13 };
    const IndexType values2[] = { 7, 3, 14, 2 };

    const IndexType n = sizeof( values1 ) / sizeof( IndexType );

    for ( IndexType i = 0; i < BinaryOp::MAX_BINARY_OP; ++i )
    {
        BinaryOp op = BinaryOp( i );

        HArray<IndexType> array1( ctx );
        HArray<IndexType> array2( ctx );
        HArray<IndexType> array3( ctx );

        array1.setRawData( n, values1 );
        array2.setRawData( n, values2 );

        if ( ! isBinarySupported<IndexType>( op ) )
        {
            BOOST_CHECK_THROW(
            {
                HArrayUtils::binaryOp( array3, array1, array2, op, ctx );
            }, Exception );

            continue;  // not all operations are supported for IndexType
        }

        HArrayUtils::binaryOp( array3, array1, array2, op, ctx );

        BOOST_REQUIRE_EQUAL( n, array3.size() );

        ReadAccess<IndexType> read( array3, host );  // read result array

        for ( IndexType i = 0; i < n; ++i )
        {
            IndexType res = applyBinary( values1[i], op, values2[i] );

            BOOST_CHECK_EQUAL( res, read[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( binaryOpTestScalar1, ValueType, scai_numeric_test_types )
{
    // test array operations: array = scalar <op> array

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    // be careful about values, e.g. pow( x, y ), x >= 0 

    const ValueType scalar = 3.5;
    const ValueType values[] = { 1.0, 1.2, 2.0, 1.3 };

    const IndexType n = sizeof( values ) / sizeof( ValueType );

    for ( IndexType i = 0; i < BinaryOp::MAX_BINARY_OP; ++i )
    {
        BinaryOp op = BinaryOp( i );

        HArray<ValueType> array1( ctx );
        HArray<ValueType> array2( ctx );

        array1.setRawData( n, values );

        // array2 = scalar <op> array1

        HArrayUtils::compute( array2, scalar, op, array1, ctx );

        BOOST_REQUIRE_EQUAL( n, array2.size() );

        ReadAccess<ValueType> read( array2, host );  // read result array

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType res = applyBinary( scalar, op, values[i] );

            typedef typename TypeTraits<ValueType>::AbsType AbsType;

            // might happen that result on other devices are not exactly the same

            AbsType diff = common::Math::abs( read[i] - res  );

            BOOST_CHECK( diff <= TypeTraits<AbsType>::small() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( binaryOpTestScalar2, ValueType, scai_numeric_test_types )
{
    // test array operations: array2 = array1 <op> scalar

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    const ValueType scalar = 3.5;
    const ValueType values[] = { 1.0, 1.2, 2.0, 1.3 };

    const IndexType n = sizeof( values ) / sizeof( ValueType );

    for ( IndexType i = 0; i < BinaryOp::MAX_BINARY_OP; ++i )
    {
        BinaryOp op = BinaryOp( i );

        HArray<ValueType> array1( ctx );
        HArray<ValueType> array2( ctx );

        array1.setRawData( n, values );

        HArrayUtils::compute( array2, array1, op, scalar, ctx );

        BOOST_REQUIRE_EQUAL( n, array2.size() );

        ReadAccess<ValueType> read( array2, host );  // read result array

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType res = applyBinary( values[i], op, scalar );

            typedef typename TypeTraits<ValueType>::AbsType AbsType;

            // might happen that result on other devices are not exactly the same

            AbsType diff = common::Math::abs( read[i] - res  );

            SCAI_LOG_TRACE( logger, values[i] << " " << op << " " << scalar << ", array[" << i << "] = " << read[i] << ", res = " << res )

            BOOST_CHECK( diff <= TypeTraits<AbsType>::small() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( binaryOpSparseNewTest, ValueType, scai_numeric_test_types )
{
    // check of all UnaryOp array operations

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    // be careful about values

    const ValueType values1[]  = { 1.0, 1.2, 2.0, 1.3 };
    const IndexType indexes1[] = { 2, 4, 6, 7 };
    const ValueType values2[]  = { 0.5, 0.7, 0.3, 1.3 };
    const IndexType indexes2[] = { 2, 3, 6, 8 };

    const IndexType n = 10;

    const IndexType nnz1 = sizeof( values1 ) / sizeof( ValueType );
    const IndexType nnz2 = sizeof( values2 ) / sizeof( ValueType );

    for ( IndexType i = 0; i < BinaryOp::MAX_BINARY_OP; ++i )
    {
        BinaryOp op = BinaryOp( i );

        HArray<IndexType> ia1( ctx );
        HArray<IndexType> ia2( ctx );
        HArray<IndexType> ia3( ctx );

        HArray<ValueType> array1( ctx );
        HArray<ValueType> array2( ctx );
        HArray<ValueType> array3( ctx );

        ia1.setRawData( nnz1, indexes1 );
        array1.setRawData( nnz1, values1 );
        ia2.setRawData( nnz1, indexes2 );
        array2.setRawData( nnz2, values2 );

        ValueType zero1 = 1;
        ValueType zero2 = 2;

        HArrayUtils::binaryOpSparse( ia3, array3, ia1, array1, zero1, ia2, array2, zero2, op, ctx );

        BOOST_REQUIRE_EQUAL( ia3.size(), array3.size() ); 
        BOOST_REQUIRE( ia3.size() >= ia2.size() ); 
        BOOST_REQUIRE( ia3.size() >= ia1.size() ); 

        IndexType pos1 = 0;
        IndexType pos2 = 0;
        IndexType pos3 = 0;

        ReadAccess<IndexType> rIA1( ia1, host );  // read result array
        ReadAccess<IndexType> rIA2( ia2, host );  // read result array
        ReadAccess<IndexType> rIA3( ia3, host );  // read result array

        ReadAccess<ValueType> rValues1( array1, host );  // read result array
        ReadAccess<ValueType> rValues2( array2, host );  // read result array
        ReadAccess<ValueType> rValues3( array3, host );  // read result array

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType res = 0;

            if ( rIA3[pos3] != i )
            {
                continue;
            }

            if ( rIA2[pos2] == i && rIA1[pos1] == i )
            {
                res = applyBinary( rValues1[pos1], op, rValues2[pos2] );
                pos1++;
                pos2++;
            }
            else if ( rIA1[pos1] == i )
            {
                res = applyBinary( rValues1[pos1], op, zero2 );
                pos1++;
            }
            else if ( rIA2[pos2] == i )
            {
                res = applyBinary( zero1, op, rValues2[pos2] );
                pos2++;
            }
            else
            {
                BOOST_CHECK( false );  // fail here
            }

            typedef typename TypeTraits<ValueType>::AbsType AbsType;

            AbsType diff = common::Math::abs( rValues3[pos3] - res  );

            SCAI_LOG_TRACE( logger, "res = " << res << ", computed " << rValues3[pos3] )

            BOOST_CHECK( diff <= TypeTraits<AbsType>::small() );
 
            pos3++;
        }

        BOOST_CHECK_EQUAL( pos1, array1.size() );
        BOOST_CHECK_EQUAL( pos2, array2.size() );
        BOOST_CHECK_EQUAL( pos3, array3.size() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( binaryOpSparseSameTest, ValueType, scai_numeric_test_types )
{
    // check of all binary ops on sparse arrays with same pattern

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    // be careful about values, must be valid arguments for all binary ops

    const IndexType indexes[] = { 2, 3, 6, 8 };
    const ValueType values1[] = { 1.0, 1.2, 2.0, 1.3 };
    const ValueType values2[] = { 0.5, 0.7, 0.3, 1.3 };

    const IndexType nnz  = sizeof( indexes ) / sizeof( IndexType );

    const IndexType nnz1 = sizeof( values1 ) / sizeof( ValueType );
    const IndexType nnz2 = sizeof( values2 ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( nnz, nnz1, "serious size mismatch" )
    SCAI_ASSERT_EQ_ERROR( nnz, nnz2, "serious size mismatch" )

    for ( IndexType i = 0; i < BinaryOp::MAX_BINARY_OP; ++i )
    {
        BinaryOp op = BinaryOp( i );

        LArray<IndexType> ia1( ctx );
        LArray<IndexType> ia2( ctx );
        LArray<IndexType> ia3( ctx );

        LArray<ValueType> array1( ctx );
        LArray<ValueType> array2( ctx );
        LArray<ValueType> array3( ctx );

        ia1.setRawData( nnz, indexes );
        array1.setRawData( nnz, values1 );
        ia2.setRawData( nnz, indexes );
        array2.setRawData( nnz, values2 );

        // zero values are not really needed if indexes are all same

        ValueType zero1 = 1;
        ValueType zero2 = 2;

        HArrayUtils::binaryOpSparse( ia3, array3, ia1, array1, zero1, ia2, array2, zero2, op, ctx );

        BOOST_REQUIRE_EQUAL( ia3.size(), array3.size() ); 
        BOOST_REQUIRE_EQUAL( nnz, ia3.size() );

        // ia3 must be equal to ia1, ia2

        BOOST_REQUIRE_EQUAL( IndexType( 0 ), ia3.maxDiffNorm( ia1 ) );

        // array3 must be array1 op array2 for all elements

        ReadAccess<ValueType> rValues1( array1, host );  // read result array
        ReadAccess<ValueType> rValues2( array2, host );  // read result array
        ReadAccess<ValueType> rValues3( array3, host );  // read result array

        typedef typename TypeTraits<ValueType>::AbsType AbsType;

        for ( IndexType i = 0; i < nnz; ++i )
        {
            ValueType expectedValue = applyBinary( rValues1[i], op, rValues2[i] );

            AbsType diff = common::Math::abs( expectedValue - rValues3[i] );

            BOOST_CHECK( diff <= TypeTraits<AbsType>::small() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( allSparseTest, ValueType, scai_numeric_test_types )
{
    // check of all binary ops on sparse arrays with same pattern

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    // SparseArray1 = { 1, 2, 3, 2 }

    ValueType zero1 = 2;
    const IndexType indexes1[] = { 0, 2 };
    const ValueType values1[]  = { 1, 3 };

    // SparseArray2 = { 1, 2, 3, 2 }

    ValueType zero2 = 1;
    const IndexType indexes2[] = { 1, 2, 3 };
    const ValueType values2[]  = { 2, 3, 2 };

    const IndexType nnz1 = sizeof( values1 ) / sizeof( ValueType );
    const IndexType nnz2 = sizeof( values2 ) / sizeof( ValueType );

    LArray<IndexType> ia1( ctx );
    LArray<IndexType> ia2( ctx );

    LArray<ValueType> array1( ctx );
    LArray<ValueType> array2( ctx );

    ia1.setRawData( nnz1, indexes1 );
    array1.setRawData( nnz1, values1 );
    ia2.setRawData( nnz2, indexes2 );
    array2.setRawData( nnz2, values2 );

    // zero values are not really needed if indexes are all same

    bool allFlag;

    IndexType n = HArrayUtils::allSparse( allFlag, ia1, array1, zero1, ia2, array2, zero2, CompareOp::EQ, ctx );

    BOOST_CHECK_EQUAL( n, 4 );  
    BOOST_CHECK( allFlag );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( copysignTest, ValueType, scai_numeric_test_types )
{
    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();
    const ValueType magnitude[] = {  1.0, 1.2,  1.3, 1.0 };
    const ValueType sign[]      = { -1.0, 1.0, -2.0, 2.0 };
    const ValueType result[]    = { -1.0, 1.2, -1.3, 1.0 };
    const IndexType n = sizeof( magnitude ) / sizeof( ValueType );
    HArray<ValueType> magArray( ctx );
    HArray<ValueType> signArray( ctx );
    HArray<ValueType> resultArray( ctx );
    magArray.setRawData( n, magnitude );
    signArray.setRawData( n, sign );
    resultArray.setRawData( n, result );
    HArrayUtils::binaryOp( resultArray, magArray, signArray, BinaryOp::COPY_SIGN, ctx );
    {
        ReadAccess<ValueType> readMag( magArray, host );
        ReadAccess<ValueType> readSign( signArray, host );
        ReadAccess<ValueType> readResult( resultArray, host );

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType x = readResult[i] - common::Math::copysign( readMag[i], readSign[i] );
            BOOST_CHECK_SMALL( common::Math::real( x ), common::TypeTraits<ValueType>::small() );
            BOOST_CHECK_SMALL( common::Math::imag( x ), common::TypeTraits<ValueType>::small() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( GatherTest, ValueType, scai_numeric_test_types )
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
    HArrayUtils::gatherImpl( target, source, indexes, BinaryOp::COPY );

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType x1 = source[indexes[i]];
        ValueType x2 = target[i];
        BOOST_CHECK_EQUAL( x1, x2 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SparseGatherTest, ValueType, scai_numeric_test_types )
{
    // Define the sparse array, IA, VA, nnz 

    ValueType sourceVA_vals[] = { 3, 1, 4, 2 };
    IndexType sourceIA_vals[] = { 0, 1, 5, 7 };

    const IndexType M = 9;

    const IndexType nnz = sizeof( sourceVA_vals ) / sizeof( ValueType );
    const IndexType nnz1 = sizeof( sourceIA_vals ) / sizeof( IndexType );

    LArray<ValueType> sourceVA( nnz, sourceVA_vals );
    LArray<IndexType> sourceIA( nnz1, sourceIA_vals );

    BOOST_CHECK_EQUAL( sourceVA.size(), sourceIA.size() );

    IndexType indexes_vals[]  = { 0, 2, 1, 2, 7, 5 };
    ValueType target_vals[]   = { 3, 0, 1, 0, 2, 4 };   // expected

    const IndexType N = sizeof( indexes_vals ) / sizeof( IndexType );

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_REQUIRE( indexes_vals[i] < M );
    }

    // target = sourceSparse[ indexes ]

    LArray<IndexType> indexes( N, indexes_vals );

    LArray<ValueType> target;

    HArrayUtils::sparseGather( target, sourceVA, sourceIA, indexes, BinaryOp::COPY );

    BOOST_REQUIRE_EQUAL( target.size(), indexes.size() );

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType t = target[i];
        BOOST_CHECK_EQUAL( t, target_vals[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ScatterTest, ValueType, scai_numeric_test_types )
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

    bool uniqueIndexes = true;

    // target = source[ indexes ]
    LArray<ValueType> source( N, sourceVals );
    LArray<IndexType> indexes( N, indexVals );
    BOOST_CHECK_THROW (
    {
        // scatter on an empty array with more than one index should always fail

        LArray<ValueType> target;
        HArrayUtils::scatterImpl( target, indexes, uniqueIndexes, source, BinaryOp::COPY );

    }, Exception );
    LArray<ValueType> target( M );
    HArrayUtils::scatterImpl( target, indexes, uniqueIndexes, source, BinaryOp::COPY );

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType x1 = target[indexes[i]];
        ValueType x2 = source[i];
        BOOST_CHECK_EQUAL( x1, x2 );
    }
}

/* --------------------------------------------------------------------- */

typedef boost::mpl::list<IndexType, SCAI_NUMERIC_TYPES_HOST> array_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( scan1Test, ValueType, array_types )
{
    ContextPtr loc = Context::getContextPtr();
    ValueType vals[]     = { 3, 1, 4, 2 };
    const IndexType n = sizeof( vals ) / sizeof( ValueType );
    unique_ptr<ValueType[]> scans( new ValueType[n + 1] );
    scans[0] = 0;

    for ( IndexType i = 0; i < n; ++i )
    {
        scans[i + 1] = scans[i] + vals[i];
    }

    LArray<ValueType> array;
    array.reserve( loc, n + 1 );
    array.setRawData( n, vals );
    LArray<ValueType> correct( n + 1, scans.get(), loc );
    SCAI_LOG_DEBUG( logger, "scan1( " << array << " )" )
    ValueType total = HArrayUtils::scan1( array );
    ValueType lastVal = array[n];
    BOOST_CHECK_EQUAL( array.size(), n + 1 );
    BOOST_CHECK_EQUAL( array.maxDiffNorm( correct ), ValueType( 0 ) );
    BOOST_CHECK_EQUAL( total, lastVal );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scanTest, ValueType, array_types )
{
    ContextPtr loc = Context::getContextPtr();

    ValueType first      = 1;
    ValueType vals[]     = { 3, 1, 4, 2 };
    const IndexType n = sizeof( vals ) / sizeof( ValueType );
    unique_ptr<ValueType[]> scans( new ValueType[n] );

    scans[0] = first + vals[0];

    for ( IndexType i = 1; i < n; ++i )
    {
        scans[i] = scans[i-1] + vals[i];
    }

    ValueType last = scans[n-1];

    LArray<ValueType> array( n, vals, loc );
    LArray<ValueType> correct( n, scans.get(), loc );

    bool exclusive = false;

    ValueType total = HArrayUtils::scan( array, first, exclusive, loc );

    BOOST_CHECK_EQUAL( array.size(), n );
    BOOST_CHECK_EQUAL( array.maxDiffNorm( correct ), ValueType( 0 ) );
    BOOST_CHECK_EQUAL( total, last );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( unscanTest, ValueType, array_types )
{
    ContextPtr loc = Context::getContextPtr();
    ValueType vals[]     = { 0, 4, 7, 11, 15 };
    const IndexType n = sizeof( vals ) / sizeof( ValueType );
    unique_ptr<ValueType[]> unscans( new ValueType[n - 1] );
    unscans[0] = 0;

    for ( IndexType i = 0; i < n - 1; ++i )
    {
        unscans[i] = vals[i + 1] - vals[i];
    }

    LArray<ValueType> array( n, vals, loc );
    LArray<ValueType> correct( n - 1, unscans.get(), loc );

    ValueType firstVal  = array[0];
    ValueType returnVal = HArrayUtils::unscan( array );

    BOOST_REQUIRE_EQUAL( array.size(), n - 1 );
    BOOST_CHECK_EQUAL( array.maxDiffNorm( correct ), ValueType( 0 ) );
    BOOST_CHECK_EQUAL( firstVal, returnVal );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( randomTest, ValueType, array_types )
{
    ContextPtr loc = Context::getContextPtr();
    const IndexType n = 100;
    LArray<ValueType> array( n, loc );
    array.setRandom( 1 );
    ValueType sum = array.sum();

    SCAI_LOG_INFO( logger, "Draw " << n << " random values of type " << common::TypeTraits<ValueType>::id()
                   << ", sum up to " << sum )

    if ( typeid( ValueType ).name() == typeid( IndexType ).name() )
    {
        // no check for integer types
    }
    else
    {
        typedef typename TypeTraits<ValueType>::AbsType AbsType;
        // AbsType asum = Math::abs( sum );
        AbsType asum = sum;   // only real part for complex numbers
        // random numbers are between 0 and 1.0, so should sum up approximately to n/2
        BOOST_CHECK( AbsType( asum - n / 2 )  < AbsType( n / 5 ) );
    }
  
    float fillRate = 0.1f;

    HArray<IndexType> nonZeroIndexes;
    HArrayUtils::randomSparseIndexes( nonZeroIndexes, 10 * n, fillRate );

    // nonZeroIndexes should have size around n 

    SCAI_LOG_DEBUG( logger, "nonZeroIndexes.size() = " << nonZeroIndexes.size() << ", should be close to " << n )

    // not possible: BOOST_CHECK_CLOSE( n, nonZeroIndexes.size(), 20  );
    // not possible for unsigned int abs( n - nonZeroIndexes.size() )

    BOOST_CHECK( common::applyBinary( n, common::BinaryOp::ABS_DIFF, nonZeroIndexes.size() ) < ( n / 5 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sortPermTest, ValueType, array_types )
{
    typedef typename common::TypeTraits<ValueType>::AbsType RealValueType;

    ContextPtr loc = Context::getContextPtr();

    RealValueType vals[] = { 13, 5, 14, 2 };
    const IndexType n = sizeof( vals ) / sizeof( RealValueType );
    LArray<RealValueType> array( n, vals, loc );

    LArray<IndexType> perm;
    LArray<RealValueType> array1;
    HArrayUtils::sort( &perm, &array1, array, true );

    BOOST_REQUIRE_EQUAL( perm.size(), n );
    BOOST_REQUIRE_EQUAL( array1.size(), n );

    LArray<RealValueType> array2;   // = array[perm]
    HArrayUtils::gatherImpl( array2, array, perm, BinaryOp::COPY );

    BOOST_CHECK_EQUAL( array2.maxDiffNorm( array1 ), RealValueType( 0 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sortSparseTest, ValueType, array_types )
{
    typedef typename common::TypeTraits<ValueType>::AbsType RealValueType;

    ContextPtr loc = Context::getContextPtr();

    IndexType raw_indexes[] = { 13, 5, 14, 2, 1, 9, 16, 31, 32, 17 };
    RealValueType raw_vals[] = { 13, 5, 14, 2, 1, 9, 16, 31, 32, 17 };

    const IndexType n = sizeof( raw_vals ) / sizeof( RealValueType );

    LArray<IndexType> indexes( n, raw_indexes, loc );
    LArray<RealValueType> values( n, raw_vals, loc );

    for ( IndexType i = 0; i < 2; ++i )
    { 
        bool ascending = i == 0 ? false : true;
        CompareOp op = ascending ? CompareOp::LE : CompareOp::GE;

        HArrayUtils::sortSparseEntries( indexes, values, ascending );

        BOOST_CHECK( HArrayUtils::isSorted( indexes, op ) );
        BOOST_CHECK( HArrayUtils::isSorted( values, op ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( inversePermTest )
{
    IndexType valuesPerm[]   = { 1, 0, 3, 2, 5, 6, 7, 4 };
    IndexType expectedPerm[] = { 1, 0, 3, 2, 7, 4, 5, 6 };

    ContextPtr loc = Context::getContextPtr();

    const IndexType n = sizeof( valuesPerm ) / sizeof( IndexType );
    LArray<IndexType> perm( n, valuesPerm, loc );
    LArray<IndexType> expectedInvPerm( n, expectedPerm, loc );
    
    LArray<IndexType> invPerm;

    HArrayUtils::inversePerm( invPerm, perm, loc );
    BOOST_CHECK_EQUAL( expectedInvPerm.maxDiffNorm( invPerm ), IndexType( 0 ) );

    perm[0] = 0;

    // now is perm no more a permutation

    BOOST_CHECK_THROW (
    {
        HArrayUtils::inversePerm( invPerm, perm, loc );
    }, 
    common::Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sortValuesTest, ValueType, array_types )
{
    typedef typename common::TypeTraits<ValueType>::AbsType RealValueType;

    ContextPtr loc = Context::getContextPtr();

    bool descending = false;
    bool ascending = true;

    RealValueType vals[] = { 13, 5, 14, 2, 8, 2, 1 };

    const IndexType n = sizeof( vals ) / sizeof( RealValueType );
    LArray<RealValueType> array( n, vals, loc );

    HArrayUtils::sort( NULL, &array, array, descending );
    BOOST_CHECK( HArrayUtils::isSorted( array, CompareOp::GE ) );

    RealValueType expectedVals[] = { 1, 2,  2, 5, 8, 13, 14 };
    LArray<RealValueType> array1( n, expectedVals, loc );
    LArray<RealValueType> array2;
    HArrayUtils::sort( NULL, &array2, array, ascending );

    BOOST_CHECK_EQUAL( array2.maxDiffNorm( array1 ), RealValueType( 0 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( bucketSortTest )
{
    ContextPtr loc = Context::getContextPtr();

    LArray<IndexType> emptyArray;

    LArray<IndexType> perm;
    LArray<IndexType> offsets;

    IndexType numBuckets = 5;

    HArrayUtils::bucketSort( offsets, perm, emptyArray, numBuckets, loc );

    BOOST_CHECK_EQUAL( perm.size(), IndexType( 0 ) );
    BOOST_CHECK_EQUAL( offsets.size(), numBuckets + 1 );

    IndexType vals[] = { 1, 3, 2, 0, 1, 3, 1 , 2, 0, 1 };
    const IndexType n = sizeof( vals ) / sizeof( IndexType );
    LArray<IndexType> array( n, vals, loc );
    numBuckets = array.max() + 1;

    HArrayUtils::bucketSort( offsets, perm, array, numBuckets, loc );

    BOOST_CHECK_EQUAL( offsets.size(), numBuckets + 1 );
    BOOST_CHECK_EQUAL( perm.size(), n );

    LArray<IndexType> sortedArray;
    HArrayUtils::gatherImpl( sortedArray, array, perm, BinaryOp::COPY );
    BOOST_CHECK( HArrayUtils::isSorted( sortedArray, CompareOp::LE, loc ) );

    // number of buckets = 1, so only two values array[i] == 0 are taken

    numBuckets = 1;
    HArrayUtils::bucketSort( offsets, perm, array, numBuckets, loc );
    BOOST_CHECK_EQUAL( perm.size(), IndexType( 2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( bucketCountTest )
{
    ContextPtr loc = Context::getContextPtr();

    LArray<IndexType> emptyArray;

    LArray<IndexType> sizes;

    IndexType numBuckets = 5;

    HArrayUtils::bucketCount( sizes, emptyArray, numBuckets, loc );

    BOOST_CHECK_EQUAL( sizes.size(), numBuckets );

    IndexType vals[] = { 1, 3, 2, 0, 1, 3, 1 , 2, 0, 1 };
    const IndexType n = sizeof( vals ) / sizeof( IndexType );
    LArray<IndexType> array( n, vals, loc );

    numBuckets = array.max() + 1;
    HArrayUtils::bucketCount( sizes, array, numBuckets, loc );
    BOOST_CHECK_EQUAL( sizes.size(), numBuckets );
    BOOST_CHECK_EQUAL( sizes.sum(), n );

    BOOST_CHECK_EQUAL( IndexType( 4 ), numBuckets );
    BOOST_CHECK_EQUAL( IndexType( 2 ), sizes[0] );
    BOOST_CHECK_EQUAL( IndexType( 4 ), sizes[1] );
    BOOST_CHECK_EQUAL( IndexType( 2 ), sizes[2] );
    BOOST_CHECK_EQUAL( IndexType( 2 ), sizes[3] );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( mergeSortTest )
{
    typedef RealType ValueType;

    ContextPtr loc = Context::getContextPtr();

    ValueType vals[]    = { 0, 5, 18, 3, 6, 9, 10, 11, 12 };
    IndexType offsets[] = { 0,        3, 3,        6,        9 };

    // Note: one empty subarray

    const IndexType n  = sizeof( vals ) / sizeof( ValueType );        // number values to sort
    const IndexType nb = sizeof( offsets ) / sizeof( IndexType ) - 1; // number of sorted subarray

    BOOST_REQUIRE_EQUAL( IndexType( 0 ), offsets[0] );
    BOOST_REQUIRE_EQUAL( IndexType( n ), offsets[nb] );

    // make LAMA Arrays of the data to use it in different context

    LArray<ValueType> array( n, vals, loc );
    LArray<IndexType> sortOffsets( nb + 1, offsets, loc );

    bool ascending = true;

    HArrayUtils::mergeSort( array, sortOffsets, ascending );

    BOOST_CHECK( HArrayUtils::isSorted( array, CompareOp::LE, loc ) );

    array.setRawData( n, vals );

    LArray<IndexType> perm;                   // will contain the sort permutation
    HArrayUtils::setOrder( perm, n, loc );

    HArrayUtils::mergeSort( array, perm, sortOffsets, ascending );

    // perm = [ 0, 3, 1, 4, 5, 6, 7, 8, 2 ]

    LArray<ValueType> arrayUnsorted( n, vals, loc );

    for ( IndexType i = 0; i < n; ++i )
    {
        SCAI_LOG_DEBUG( logger, "array_sorted[" << i << "] = " << array[i]
                        << ", perm[" << i << "] = " << perm[i]
                        << ", array_unsorted[" << i << "] = " << arrayUnsorted[i] )
    }

    // check that array_sorted  ==  array_unsorted[perm]

    HArrayUtils::gather( array, arrayUnsorted, perm, BinaryOp::SUB, loc );

    BOOST_CHECK_EQUAL( 0, array.maxNorm() );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( setSequenceTest, ValueType, array_types )
{
    ContextPtr loc = Context::getContextPtr();
    const IndexType n = 10;
    const ValueType start = 5;
    const ValueType inc = 10;
    LArray<ValueType> array;
    HArrayUtils::setSequence( array, start, inc, n, loc );
    BOOST_CHECK_EQUAL( array.size(), n );

    for ( IndexType i = 0; i < n; ++i )
    {
        ValueType elem = array[i];
        BOOST_CHECK_EQUAL( start + static_cast<ValueType>( i ) * inc, elem );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setArraySectionTest, ValueType, array_types )
{
    ContextPtr loc = Context::getContextPtr();

    ValueType xVals[] = { 3, 1, 4, 2, 2, 1, 4, 3 };
    ValueType yVals[] = { 2, 1, 1, 5, 5, 3, 9 };

    const IndexType nx = sizeof( xVals ) / sizeof( ValueType );
    const IndexType ny = sizeof( yVals ) / sizeof( ValueType );

    LArray<ValueType> x( nx, xVals, loc );
    LArray<ValueType> y( ny, yVals, loc );

    const IndexType n = 3;
    const IndexType ofs_x = 1;
    const IndexType inc_x = 2;
    const IndexType ofs_y = 1;
    const IndexType inc_y = 2;

    HArrayUtils::setArraySection( x, ofs_x, inc_x, y, ofs_y, inc_y, n, BinaryOp::COPY, loc );

    for ( IndexType i = 0; i < n; ++i )
    {
        ValueType elem = x[ ofs_x + i * inc_x];
        BOOST_CHECK_EQUAL( yVals[ ofs_y + i * inc_y ], elem );
    }
}

/* --------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_ALL_TYPES> all_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( setArrayFailTest, ValueType, all_types )
{
    ContextPtr loc = Context::getContextPtr();

    ValueType xVals[] = { 3, 1, 4, 2, 2, 1, 4 };

    const IndexType nx = sizeof( xVals ) / sizeof( ValueType );

    HArray<ValueType> x( nx, xVals, loc );
    HArray<ValueType> y;

    // verify that unsupported array types throw an exception

    if ( common::mepr::TypeListUtilsV<ValueType, SCAI_ARRAY_TYPES_HOST_LIST>::contains )
    {
        SCAI_LOG_INFO( logger, "supported value type " << TypeTraits<ValueType>::id() << " can be used in HArray utilities" )

        HArrayUtils::setArray( y, x, BinaryOp::COPY, loc );

        BOOST_REQUIRE_EQUAL( y.size(), x.size() );

        {
            ReadAccess<ValueType> rY( y );

            for ( IndexType i = 0; i < nx; ++i )
            {
                BOOST_CHECK_EQUAL( rY[i], xVals[i] );
            }
        }
    }
    else
    {
        SCAI_LOG_INFO( logger, "unsupported value type " << TypeTraits<ValueType>::id() << " cannot be used in HArray utilities" )

        BOOST_CHECK_THROW(
        {
            HArrayUtils::setArray( y, x, BinaryOp::COPY, loc );
        },
        common::Exception );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( axpyTest, ValueType, scai_numeric_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( arrayPlusArrayTest, ValueType, scai_numeric_test_types )
{
    ContextPtr loc = Context::getContextPtr();
    ValueType sourceVals1[] = { 3, 1, 4, 2 };
    ValueType sourceVals2[] = { 2, -1, -1, -5 };
    const IndexType n1 = sizeof( sourceVals1 ) / sizeof( ValueType );
    const IndexType n2 = sizeof( sourceVals2 ) / sizeof( ValueType );
    BOOST_REQUIRE_EQUAL( n1, n2 );
    LArray<ValueType> x1( n1, sourceVals1, loc );
    LArray<ValueType> x2( n2, sourceVals2, loc );
    LArray<ValueType> xf( n2 - 1, sourceVals2, loc ); // wrong sized array
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

BOOST_AUTO_TEST_CASE_TEMPLATE( arrayPlusScalarTest, ValueType, scai_numeric_test_types )
{
    ContextPtr loc = Context::getContextPtr();
    ValueType sourceVals1[] = { 3, 1, 4, 2 };
    const IndexType n1 = sizeof( sourceVals1 ) / sizeof( ValueType );
    LArray<ValueType> x1( n1, sourceVals1, loc );
    LArray<ValueType> target( loc );
    ValueType factors[] = { -1, 0, 1, 2 };
    int NCASE = sizeof( factors ) / sizeof( ValueType );

    for ( int k1 = 0; k1 < NCASE; ++k1 )
    {
        ValueType alpha = factors[k1];

        for ( int k2 = 0; k2 < NCASE; ++k2 )
        {
            ValueType beta = factors[k2];
            SCAI_LOG_DEBUG( logger, "target = " << alpha << " * x1 + " << beta )
            target.purge();

            HArrayUtils::arrayPlusScalar( target, alpha, x1, beta, loc );

            BOOST_CHECK_EQUAL( x1.size(), target.size() );

            for ( IndexType i = 0; i < n1; ++i )
            {
                ValueType v = target[i];
                BOOST_CHECK_EQUAL( v, alpha * sourceVals1[i] + beta );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( arrayTimesArrayTest, ValueType, scai_numeric_test_types )
{
    ContextPtr loc = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();
    ValueType sourceVals1[] = { 3,  1,  4,  2 };
    ValueType sourceVals2[] = { 2, -1, -1, -5 };
    const IndexType n1 = sizeof( sourceVals1 ) / sizeof( ValueType );
    const IndexType n2 = sizeof( sourceVals2 ) / sizeof( ValueType );
    BOOST_REQUIRE_EQUAL( n1, n2 );
    LArray<ValueType> x1( n1, sourceVals1, loc );
    LArray<ValueType> x2( n2, sourceVals2, loc );
    LArray<ValueType> xf( n2 - 1, sourceVals2, loc ); // wrong sized array
    LArray<ValueType> target( loc );
    ValueType alpha = 1.0;
    ValueType beta = 2.0;
    ValueType factors[] = { 6, -1, -4, -10 };
    IndexType NCASE = sizeof( factors ) / sizeof( ValueType );

    HArrayUtils::arrayTimesArray( target, alpha, x1, x2, loc );

    {
        BOOST_CHECK_EQUAL( NCASE, target.size() );
        ReadAccess<ValueType> rTarget( target, host );

        for ( IndexType i = 0; i < NCASE; ++i )
        {
            BOOST_CHECK_EQUAL( factors[i], rTarget[i] );
        }
    }

    HArrayUtils::arrayTimesArray( target, beta, x1, x2, loc );

    {
        BOOST_CHECK_EQUAL( NCASE, target.size() );
        ReadAccess<ValueType> rTarget( target, host );

        for ( IndexType i = 0; i < NCASE; ++i )
        {
            BOOST_CHECK_EQUAL( beta * factors[i], rTarget[i] );
        }
    }

    BOOST_CHECK_THROW (
    {
        HArrayUtils::arrayTimesArray( target, alpha, x1, xf, loc );
    }, Exception );
}

/* --------------------------------------------------------------------- */

// BOOST_AUTO_TEST_CASE_TEMPLATE( arrayPlusArrayAliasTest, ValueType, scai_numeric_test_types )

BOOST_AUTO_TEST_CASE( arrayPlusArrayAliasTest )
{
    typedef float ValueType;
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

BOOST_AUTO_TEST_CASE_TEMPLATE( sparseTest, ValueType, scai_numeric_test_types )
{
    ContextPtr loc = Context::getContextPtr();
    ValueType denseValues[] = { 3, 0, 0, 1, 4, 0, 2, 7, 0, -1, -3 };
    const IndexType n = sizeof( denseValues ) / sizeof( ValueType );
    LArray<ValueType> denseArray( n, denseValues, loc );
    /* Idea: convert to sparse and back to dense, must contain original data */
    LArray<ValueType> sparseArray( loc );
    LArray<IndexType> sparseIndexes( loc );
    HArrayUtils::buildSparseArray( sparseArray, sparseIndexes, denseArray, loc );
    BOOST_CHECK_EQUAL( sparseArray.size(), IndexType( 7 ) );
    BOOST_REQUIRE_EQUAL( sparseArray.size(), sparseIndexes.size() );
    // The sparse indexes must be sorted, no doubles
    BOOST_CHECK( HArrayUtils::isSorted( sparseIndexes, CompareOp::LT, loc ) );
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

BOOST_AUTO_TEST_CASE( insertTest )
{
    ContextPtr loc = Context::getContextPtr();

    typedef RealType ValueType;

    IndexType indexes[] = { 5, 7, 9, 1, 8 };
    IndexType pos[]     = { 0, 1, 2, 0, 3 };

    ValueType values[]   = { 5, 7, 9, 1, 8 };
    ValueType values1[]  = { 1, 5, 7, 8, 9 };

    const IndexType n = sizeof( indexes ) / sizeof( IndexType );

    HArray<IndexType> indexArray;

    for ( IndexType i = 0; i < n; ++i )
    {
        IndexType p = HArrayUtils::insertSorted( indexArray, indexes[i], loc );
        SCAI_LOG_TRACE( logger, i << " of " << n << ": insertIndex " << indexes[i] << " at pos = " << p << ", expected " << pos[i] )
        BOOST_CHECK_EQUAL( pos[i], p );
    }

    BOOST_CHECK_EQUAL( indexArray.size(), n );
    BOOST_CHECK( HArrayUtils::isSorted( indexArray, CompareOp::LT, loc ) );

    LArray<ValueType> valueArray;

    for ( IndexType i = 0; i < n; ++i )
    {
         HArrayUtils::insertAtPos( valueArray, pos[i], values[i] );
    }

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( valueArray[i], values1[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( findPosTest )
{
    ContextPtr loc = Context::getContextPtr();
    IndexType index_values[] = { 0, 3, 7, 11, 16, 19 };
    const IndexType n = sizeof( index_values ) / sizeof( IndexType );
    HArray<IndexType> indexArray( n, index_values, loc );

    IndexType pos = HArrayUtils::findPosInSortedIndexes( indexArray, IndexType( 5 ) );
    BOOST_CHECK_EQUAL( pos, nIndex );

    pos = HArrayUtils::findPosInSortedIndexes( indexArray, 12 );
    BOOST_CHECK_EQUAL( pos, nIndex );

    pos = HArrayUtils::findPosInSortedIndexes( indexArray, 23 );
    BOOST_CHECK_EQUAL( pos, nIndex );

    pos = HArrayUtils::findPosInSortedIndexes( indexArray, 0 );
    BOOST_CHECK_EQUAL( pos, IndexType( 0 ) );

    pos = HArrayUtils::findPosInSortedIndexes( indexArray, 11 );
    BOOST_CHECK_EQUAL( pos, IndexType( 3 ) );

    pos = HArrayUtils::findPosInSortedIndexes( indexArray, 19 );
    BOOST_CHECK_EQUAL( pos, IndexType( 5 ) );

    indexArray.resize( 1 );

    pos = HArrayUtils::findPosInSortedIndexes( indexArray, 0 );
    BOOST_CHECK_EQUAL( pos, IndexType( 0 ) );

    pos = HArrayUtils::findPosInSortedIndexes( indexArray, 19 );
    BOOST_CHECK_EQUAL( pos, nIndex );

    indexArray.clear();

    pos = HArrayUtils::findPosInSortedIndexes( indexArray, 0 );
    BOOST_CHECK_EQUAL( pos, nIndex );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sparseAddTestNew, ValueType, scai_array_test_types )
{
    // sparseArray = alpha * sparseArray1 + beta * sparseArray2, result gets new pattern

    ContextPtr testContext = Context::getContextPtr();

    const IndexType indexes1_values[]   = { 0, 2, 5, 7 };
    const IndexType indexes2_values[]   = { 1, 2, 5, 8 };
    const IndexType indexes_values[]    = { 0, 1, 2, 5, 7, 8 };

    const ValueType values1_values[] = { 1, 2, 3, 4 };
    const ValueType values2_values[] = { 5, 6, 7, 8 };
    const ValueType values_values[]  = { 1, 5, 8, 10, 4, 8 };

    IndexType n1 = sizeof( indexes1_values ) / sizeof( IndexType );
    IndexType n2 = sizeof( indexes2_values ) / sizeof( IndexType );
    IndexType n  = sizeof( indexes_values ) / sizeof( IndexType );

    SCAI_ASSERT_EQ_ERROR( n1, sizeof( values1_values ) / sizeof( ValueType ), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( n2, sizeof( values2_values ) / sizeof( ValueType ), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( n, sizeof( values_values ) / sizeof( ValueType ), "size mismatch" )

    HArray<IndexType> indexes1( n1, indexes1_values, testContext );
    HArray<IndexType> indexes2( n2, indexes2_values, testContext );

    HArray<ValueType> values1( n1, values1_values, testContext );
    HArray<ValueType> values2( n2, values2_values, testContext );

    HArray<IndexType> indexes;
    HArray<ValueType> values;

    ValueType one = 1;
    ValueType zero = 0;
    HArrayUtils::addSparse( indexes, values, indexes1, values1, zero, one, indexes2, values2, zero, one );

    BOOST_REQUIRE_EQUAL( n, indexes.size() );
    BOOST_REQUIRE_EQUAL( n, values.size() );

    LArray<IndexType> okayIndexes( n, indexes_values, testContext );
    LArray<ValueType> okayValues( n, values_values, testContext );

    BOOST_CHECK_EQUAL( IndexType( 0 ), okayIndexes.maxDiffNorm( indexes ) );
    BOOST_CHECK_EQUAL( ValueType( 0 ), okayValues.maxDiffNorm( values ) );

    // just switch the arguments

    HArrayUtils::addSparse( indexes, values, indexes1, values1, zero, one, indexes2, values2, zero, one );

    BOOST_CHECK_EQUAL( IndexType( 0 ), okayIndexes.maxDiffNorm( indexes ) );
    BOOST_CHECK_EQUAL( ValueType( 0 ), okayValues.maxDiffNorm( values ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sparseAddTestSame, ValueType, scai_array_test_types )
{
    // sparseArray = alpha * sparseArray1 + beta * sparseArray2,  both sparse arrays have same pattern

    ContextPtr testContext = Context::getContextPtr();

    const IndexType indexes_values[]    = { 0, 2, 5, 7, 8 };

    const ValueType values1_values[] = { 1, 2, 3, 4, 5 };
    const ValueType values2_values[] = { 5, 6, 7, 8, 9 };
    const ValueType values_values[]  = { 6, 8, 10, 12, 14 };

    IndexType n  = sizeof( indexes_values ) / sizeof( IndexType );

    HArray<IndexType> indexes1( n, indexes_values, testContext );
    HArray<IndexType> indexes2( n, indexes_values, testContext );

    HArray<ValueType> values1( n, values1_values, testContext );
    HArray<ValueType> values2( n, values2_values, testContext );

    HArray<IndexType> indexes;
    HArray<ValueType> values;

    ValueType one = 1;
    ValueType zero = 0;
    HArrayUtils::addSparse( indexes, values, indexes1, values1, zero, one, indexes2, values2, zero, one );

    BOOST_REQUIRE_EQUAL( n, indexes.size() );
    BOOST_REQUIRE_EQUAL( n, values.size() );

    LArray<IndexType> okayIndexes( n, indexes_values, testContext );
    LArray<ValueType> okayValues( n, values_values, testContext );

    BOOST_CHECK_EQUAL( IndexType( 0 ), okayIndexes.maxDiffNorm( indexes ) );
    BOOST_CHECK_EQUAL( ValueType( 0 ), okayValues.maxDiffNorm( values ) );

    // just switch the arguments

    HArrayUtils::addSparse( indexes, values, indexes1, values1, zero, one, indexes2, values2, zero, one );

    BOOST_CHECK_EQUAL( IndexType( 0 ), okayIndexes.maxDiffNorm( indexes ) );
    BOOST_CHECK_EQUAL( ValueType( 0 ), okayValues.maxDiffNorm( values ) );

    // alias on result and operand

    HArrayUtils::addSparse( indexes1, values1, indexes1, values1, zero, one, indexes2, values2, zero, one );

    BOOST_CHECK_EQUAL( IndexType( 0 ), okayIndexes.maxDiffNorm( indexes1 ) );
    BOOST_CHECK_EQUAL( ValueType( 0 ), okayValues.maxDiffNorm( values1 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( mergeSparseTest, ValueType, scai_array_test_types )
{
    // sparseArray = alpha * sparseArray1 + beta * sparseArray2, result gets new pattern

    ContextPtr testContext = Context::getContextPtr();

    const IndexType indexes1_values[]   = { 0,    2, 5,    7 };
    const IndexType indexes2_values[]   = {    1, 2, 5, 5,    8 };
    const IndexType indexes_values[]    = { 0, 1, 2, 5,    7, 8 };

    const ValueType values1_values[] =    { 1,    2, 3,    4 };
    const ValueType values2_values[] =    {    5, 2, 3, 3,    1 };
    const ValueType values_values[]  =    { 1, 5, 2, 3,    4, 1 };

    IndexType n1 = sizeof( indexes1_values ) / sizeof( IndexType );
    IndexType n2 = sizeof( indexes2_values ) / sizeof( IndexType );
    IndexType n  = sizeof( indexes_values ) / sizeof( IndexType );

    SCAI_ASSERT_EQ_ERROR( n1, sizeof( values1_values ) / sizeof( ValueType ), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( n2, sizeof( values2_values ) / sizeof( ValueType ), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( n, sizeof( values_values ) / sizeof( ValueType ), "size mismatch" )

    HArray<IndexType> indexes1( n1, indexes1_values, testContext );
    HArray<IndexType> indexes2( n2, indexes2_values, testContext );

    HArray<ValueType> values1( n1, values1_values, testContext );
    HArray<ValueType> values2( n2, values2_values, testContext );

    HArray<IndexType> indexes;
    HArray<ValueType> values;

    common::BinaryOp op = common::BinaryOp::COPY;    // replace enries at same place

    HArrayUtils::mergeSparse( indexes, values, indexes1, values1, indexes2, values2, op );

    BOOST_REQUIRE_EQUAL( n, indexes.size() );
    BOOST_REQUIRE_EQUAL( n, values.size() );

    LArray<IndexType> okayIndexes( n, indexes_values, testContext );
    LArray<ValueType> okayValues( n, values_values, testContext );

    BOOST_CHECK_EQUAL( IndexType( 0 ), okayIndexes.maxDiffNorm( indexes ) );
    BOOST_CHECK_EQUAL( ValueType( 0 ), okayValues.maxDiffNorm( values ) );

    // just switch the arguments

    HArrayUtils::mergeSparse( indexes, values, indexes1, values1, indexes2, values2, op );

    BOOST_CHECK_EQUAL( IndexType( 0 ), okayIndexes.maxDiffNorm( indexes ) );
    BOOST_CHECK_EQUAL( ValueType( 0 ), okayValues.maxDiffNorm( values ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */

