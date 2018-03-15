/**
 * @file HArrayTest.cpp
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
 * @brief Basic tests for Heterogeneous arrays.
 * @author Thomas Brandes
 * @date 08.07.2015
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/assert.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/HArrayRef.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/test/TestMacros.hpp>

using namespace boost;
using namespace scai;
using namespace scai::hmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HArrayTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HArrayTest" )

/* --------------------------------------------------------------------- */

template<typename ValueType>
void readTest( const ValueType values[], const IndexType N, const ValueType sum )
{
    ValueType mySum = 0;

    for ( IndexType i = 0; i < N; ++i )
    {
        mySum += values[i];
    }

    BOOST_CHECK_EQUAL( sum, mySum );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    ContextPtr contextPtr = Context::getContextPtr();
    SCAI_LOG_INFO( logger, "constructorTest on " << *contextPtr );
    const IndexType N = 100;
    const IndexType zero = 0;
    HArray<float> array;
    BOOST_CHECK_EQUAL( array.size(), zero );
    BOOST_CHECK_EQUAL( array.isValid( contextPtr ), false );
    array.resize( N );
    BOOST_CHECK_EQUAL( array.size(), N );
    BOOST_CHECK_EQUAL( array.capacity( contextPtr ), zero );
    BOOST_CHECK_EQUAL( array.isValid( contextPtr ), false );
    {
        WriteAccess<float> write( array, contextPtr );
        // even if no values are set here the HArray assumes that it is written
    }
    BOOST_CHECK_EQUAL( array.isValid( contextPtr ), true );
    BOOST_CHECK_EQUAL( array.capacity( contextPtr ), N );

    if ( contextPtr->getType() == common::ContextType::Host )
    {
        SCAI_LOG_INFO( logger, "constructorTest +++ on " << *contextPtr );
        // test the other constructors just on the Host
        HArray<float> array1( N );
        HArray<float> array2( N, 5.0f );
        BOOST_CHECK_EQUAL( array1.size(), N );
        BOOST_CHECK_EQUAL( array2.size(), N );
        BOOST_CHECK_EQUAL( array1.capacity( contextPtr ), N );
        BOOST_CHECK_EQUAL( array2.capacity( contextPtr ), N );
        BOOST_CHECK_EQUAL( array1.isValid( contextPtr ), false );
        BOOST_CHECK_EQUAL( array2.isValid( contextPtr ), true );
        // read access on uninitialized non-empty array
        // gives a warning as there is no valid data
        {
            ReadAccess<float>read( array1, contextPtr );
        }
        {
            ReadAccess<float> read( array2, contextPtr );
            readTest<float>( read.get(), N, N * 5 );
        }
    }
}

BOOST_AUTO_TEST_CASE( initializerListConstructorTest )
{
    const auto context = Context::getContextPtr();
    const auto array = HArray<int> ( { 5, 6, 2, 3, -5 }, context );

    BOOST_CHECK_EQUAL( array.size(), 5 );
    BOOST_CHECK_EQUAL( array.capacity( context ), 5 );

    ReadAccess<int> rArray( array );
    BOOST_CHECK_EQUAL( rArray[0], 5 );
    BOOST_CHECK_EQUAL( rArray[1], 6 );
    BOOST_CHECK_EQUAL( rArray[2], 2 );
    BOOST_CHECK_EQUAL( rArray[3], 3 );
    BOOST_CHECK_EQUAL( rArray[4], -5 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( emptyInitializerListConstructorTest )
{
    const auto context = Context::getContextPtr();
    const auto array = HArray<IndexType> ( {  }, context );

    BOOST_CHECK_EQUAL( array.size(), 0 );
    BOOST_CHECK_EQUAL( array.capacity( context ), 0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getTest )
{
    const auto context = Context::getContextPtr(); 

    const IndexType pos = 3;

    const auto array = HArray<int> ( { 5, 6, 2, 3, -5 }, context );

    int val = array[pos];

    ReadAccess<int> rArray( array );
    BOOST_CHECK_EQUAL( rArray[pos], val );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( setTest )
{
    const auto context = Context::getContextPtr();

    const IndexType pos = 3;

    auto array = HArray<int> ( 5, 0, context );

    int val = 2;

    array[pos] = val;

    ReadAccess<int> rArray( array );
    BOOST_CHECK_EQUAL( rArray[pos], val );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( moveConstructorTest )
{
    ContextPtr contextPtr = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "moveConstructorTest on " << *contextPtr );

    typedef SCAI_TEST_TYPE ValueType;

    const IndexType N = 10;
    ValueType val = 1;
    HArray<ValueType> array( N, val, contextPtr );

    // save the pointer to the allocated data @ context

    const ValueType* ptr = NULL;
    {
        ReadAccess<ValueType> ra( array, contextPtr );
        ptr = ra.get();
    }

    HArray<ValueType> marray( std::move( array ) );

    BOOST_CHECK_EQUAL( marray.size(), N );

    BOOST_CHECK( !array.isValid( contextPtr ) );
    BOOST_CHECK( marray.isValid( contextPtr ) );

    // verify that the moved array uses still the same pointer

    {
        ReadAccess<ValueType> ra( marray, contextPtr );
        BOOST_CHECK_EQUAL( ra.get(), ptr );
    }
}

/* --------------------------------------------------------------------- */


BOOST_AUTO_TEST_CASE( moveAssignmentTest )
{
    ContextPtr contextPtr = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "moveConstructorTest on " << *contextPtr );

    typedef SCAI_TEST_TYPE ValueType;

    const IndexType N = 10;
    ValueType val = 1;
    HArray<ValueType> array( N, val, contextPtr );

    // save the pointer to the allocated dat @ context

    const ValueType* ptr = NULL;
    {
        ReadAccess<ValueType> ra( array, contextPtr );
        ptr = ra.get();
    }

    HArray<ValueType> marray( N + 2, val, contextPtr );

    // in contrary to the move constructor the data of marray must be freed

    marray = std::move( array );

    BOOST_CHECK_EQUAL( marray.size(), N );
    BOOST_CHECK_EQUAL( array.size(), 0 );

    BOOST_CHECK( !array.isValid( contextPtr ) );
    BOOST_CHECK( marray.isValid( contextPtr ) );

    // verify that the moved array uses still the same pointer

    {
        ReadAccess<ValueType> ra( marray, contextPtr );
        BOOST_CHECK_EQUAL( ra.get(), ptr );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyConstructorTest )
{
    ContextPtr contextPtr = Context::getContextPtr();
    ContextPtr hostPtr = Context::getHostPtr();

    SCAI_LOG_INFO( logger, "copy test on " << *contextPtr );

    typedef SCAI_TEST_TYPE ValueType;

    const IndexType N = 10;
    ValueType val = 1;
    HArray<ValueType> array( N, val, contextPtr );

    HArray<ValueType> marray( array );

    BOOST_CHECK_EQUAL( marray.size(), N );
    BOOST_CHECK( array.isValid( contextPtr ) );
    BOOST_CHECK( marray.isValid( contextPtr ) );

    BOOST_CHECK_EQUAL( marray.isValid( hostPtr ), array.isValid( hostPtr) );

    // verify correct data after validation checks as WriteAccess changes validity

    {
        WriteAccess<ValueType> wA( marray, contextPtr );
        ValueType singleVal;
        wA.getValue( singleVal, 0 );
        BOOST_CHECK_EQUAL( singleVal, val );
    }
}

/* --------------------------------------------------------------------- */


BOOST_AUTO_TEST_CASE( assignmentTest )
{
    ContextPtr contextPtr = Context::getContextPtr();
    ContextPtr hostPtr = Context::getHostPtr();

    SCAI_LOG_INFO( logger, "assignment test on " << *contextPtr );

    typedef SCAI_TEST_TYPE ValueType;

    const IndexType N = 10;
    ValueType val = 1;
    HArray<ValueType> array( N, val, contextPtr );

    HArray<ValueType> marray( N + 2, val + 1, contextPtr );

    {
        WriteAccess<ValueType> wA( marray, contextPtr );
        ValueType singleVal;
        wA.getValue( singleVal, 0 );
        BOOST_CHECK_EQUAL( singleVal, val + 1 );
    }

    marray = array;
 
    BOOST_CHECK_EQUAL( marray.size(), N );
    BOOST_CHECK( array.isValid( contextPtr ) );
    BOOST_CHECK( marray.isValid( contextPtr ) );

    BOOST_CHECK_EQUAL( marray.isValid( hostPtr ), array.isValid( hostPtr) );

    {
        WriteAccess<ValueType> wA( marray, contextPtr );
        ValueType singleVal;
        wA.getValue( singleVal, 0 );
        BOOST_CHECK_EQUAL( singleVal, val );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assignTest )
{
    ContextPtr contextPtr = Context::getContextPtr();
    ContextPtr hostPtr = Context::getHostPtr();

    SCAI_LOG_INFO( logger, "assign test on " << *contextPtr );

    typedef SCAI_TEST_TYPE ValueType;

    const IndexType N = 10;
    ValueType val = 1;
    HArray<ValueType> array( N, val, contextPtr );

    HArray<ValueType> marray( N + 2, val + 1, contextPtr );

    marray.assign( array, hostPtr );

    BOOST_CHECK_EQUAL( marray.size(), N );
    BOOST_CHECK( array.isValid( contextPtr ) );
    BOOST_CHECK( marray.isValid( hostPtr ) );

    // self assignment also makes sure that data is valid at specified context

    array.assign( array, hostPtr );
    BOOST_CHECK( array.isValid( hostPtr ) );

}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( initTest )
{
    ContextPtr hostContext = Context::getHostPtr();
    ContextPtr testContext = Context::getContextPtr();
    const double values[] = { 1, 3, 5, 7, 4 };
    const IndexType N = sizeof( values ) / sizeof( double );
    SCAI_LOG_INFO( logger, "initTest ( " << N << " values ) on " << *testContext );
    HArray<double> array( testContext );
    array.setRawData( N, values );
    BOOST_CHECK_EQUAL( array.size(), N );
    // init values must be valid at chosen context
    BOOST_CHECK( array.isValid( testContext ) );
    // just check for correct values on the host
    {
        ReadAccess<double>read( array, hostContext );

        for ( IndexType i = 0; i < N; ++i )
        {
            BOOST_CHECK_EQUAL( values[i], read[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( init1Test )
{
    const IndexType N = 4;
    const IndexType val = 1;

    HArray<IndexType> array;
    array.setSameValue( N, val );

    BOOST_REQUIRE_EQUAL( array.size(), N );

    ContextPtr hostContext = Context::getHostPtr();

    {
        ReadAccess<IndexType>read( array, hostContext );

        for ( IndexType i = 0; i < N; ++i )
        {
            BOOST_CHECK_EQUAL( val, read[i] );
        }
    }

    IndexType N2 = N / 2;
    array.setSameValue( N2, val + 1 );

    ReadAccess<IndexType>read( array, hostContext );

    for ( IndexType i = 0; i < N2; ++i )
    {
        BOOST_CHECK_EQUAL( val + 1, read[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( resizeTest )
{
    ContextPtr contextPtr = Context::getContextPtr();
    SCAI_LOG_INFO( logger, "resizeTest with context = " << *contextPtr )
    HArray<IndexType> hArray; // default, not allocated at all
    const IndexType N = 10;
    {
        WriteAccess<IndexType> writeAccess( hArray, contextPtr );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( N );
        IndexType* data = writeAccess.get();
        const Memory& mem = writeAccess.getMemory();
        mem.memset( data, 0, N * sizeof( IndexType ) );
    }
    hArray.purge();
    {
        WriteAccess<IndexType> writeAccess( hArray, contextPtr );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( N );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( resize1Test )
{
    ContextPtr contextPtr = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "resize1Test with context = " << *contextPtr )

    const IndexType N = 5;

    HArray<IndexType> hArray; // default, not allocated at all

    {
        WriteOnlyAccess<IndexType> writeAccess( hArray, contextPtr, 0 );
    }
    {
        WriteOnlyAccess<IndexType> writeAccess( hArray, contextPtr, N );

        if ( contextPtr->getType() == common::ContextType::Host )
        {
            for ( IndexType i = 0; i < N; ++i )
            {
                writeAccess[i] = 0;
            }
        }
    }

    BOOST_CHECK_EQUAL( hArray.size(), N );
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), N );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( capacityTest )
{
    ContextPtr contextPtr = Context::getContextPtr();  // test context
    ContextPtr hostPtr    = Context::getHostPtr();
    SCAI_LOG_INFO( logger, "capacityTest on " << *contextPtr );
    const IndexType N = 10;
    const IndexType zero = 0;
    HArray<IndexType> hArray; // default, not allocated at all
    hArray.reserve( contextPtr, N );
    BOOST_CHECK_EQUAL( hArray.size(), zero );
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), N );
    {
        WriteAccess<IndexType> access( hArray, contextPtr );
        BOOST_CHECK_EQUAL( access.capacity(), N );
        access.reserve( 2 * N );
    }
    BOOST_CHECK_EQUAL( hArray.size(), zero );
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), 2 * N );
    hArray.clear();
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), 2 * N );
    hArray.purge();
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), zero );
    {
        WriteOnlyAccess<IndexType> access( hArray, hostPtr, N );
    }
    BOOST_CHECK_EQUAL( hArray.size(), N );
    BOOST_CHECK_EQUAL( hArray.capacity( hostPtr ), N );
    {
        WriteAccess<IndexType> access( hArray, contextPtr );
    }
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), N );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( swapTest )
{
    const IndexType n1 = 10;
    const IndexType n2 = 5;
    ContextPtr contextPtr = Context::getContextPtr();  // test context
    SCAI_LOG_INFO( logger, "swapTest with valid copies on " << *contextPtr );

    HArray<double> arr1( n1, 1 );
    HArray<double> arr2( n2, 2 );
    HArray<float>  arr3( n1, 3 );

    // now make them valid on test device
    {
        WriteAccess<double> write1( arr1, contextPtr );
        WriteAccess<double> write2( arr2, contextPtr );
    }
    arr1.swap( arr2 );
    BOOST_CHECK_EQUAL( arr2.size(), n1 );
    BOOST_CHECK_EQUAL( arr1.size(), n2 );
    ContextPtr hostPtr = Context::getHostPtr();
    {
        ReadAccess<double> read( arr1, hostPtr );

        for ( IndexType i = 0; i < arr1.size(); ++i )
        {
            BOOST_CHECK_EQUAL( 2, read[i] );
        }
    }
    {
        ReadAccess<double> read( arr2, hostPtr );

        for ( IndexType i = 0; i < arr2.size(); ++i )
        {
            BOOST_CHECK_EQUAL( 1, read[i] );
        }
    }

    BOOST_CHECK_THROW(
    {
        arr1._HArray::swap( arr3 );
    }, common::Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( createTest )
{
    using namespace common;
    HArray<float> A( 10, 1.0f );
    HArray<double> B( 10, 3.1415 );
    HArray<IndexType> C( 10, 5 );
    std::vector<scai::common::ScalarType> values;
    _HArray::getCreateValues( values );

    BOOST_CHECK( _HArray::canCreate( ScalarType::FLOAT ) );
    BOOST_CHECK( _HArray::canCreate( ScalarType::DOUBLE ) );
    BOOST_CHECK( _HArray::canCreate( TypeTraits<IndexType>::stype ) );
    BOOST_CHECK( !_HArray::canCreate( ScalarType::INTERNAL ) );
    _HArray* ca1 = _HArray::create( ScalarType::FLOAT );
    BOOST_REQUIRE( ca1 );
    HArray<double>* da1 = dynamic_cast<HArray<double>*>( ca1 );
    HArray<float>* fa1 = dynamic_cast<HArray<float>*>( ca1 );
    BOOST_CHECK( da1 == NULL );
    BOOST_CHECK( fa1 != NULL );
    BOOST_CHECK_THROW(
    {
        _HArray* ca1 = _HArray::create( ScalarType::INTERNAL );
        ca1->clear();
    }, Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( newArrayTest )
{
    using namespace common;
    HArray<float> A( 10, 1.0f );
  
    std::unique_ptr<_HArray> tmpA( A.newArray() );

    BOOST_CHECK_EQUAL( IndexType( 0 ), tmpA->size() );
    BOOST_CHECK_EQUAL( tmpA->getValueType(), A.getValueType() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( dynCopyTest )
{
    const IndexType N = 10;

    using namespace common;
    HArray<float> A( N, 1.0f );
  
    std::unique_ptr<_HArray> tmpA( A.copy() );

    BOOST_CHECK_EQUAL( N, tmpA->size() );
    BOOST_CHECK_EQUAL( tmpA->getValueType(), A.getValueType() );

    HArray<float>& A2 = reinterpret_cast<HArray<float>&>( *tmpA );

    WriteAccess<float> r1( A );
    WriteAccess<float> r2( A2 );

    for ( IndexType i = 0; i < 10; ++i )
    {
        BOOST_CHECK_EQUAL( r1[i], r2[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( validTest )
{
    HArray<float> A( 10 );
    ContextPtr hostContext = Context::getHostPtr();
    ContextPtr testContext = Context::getContextPtr();
    // Array not allocated at all, should also give some default for validContext
    ContextPtr validContext = A.getValidContext();
    BOOST_CHECK( validContext.get() );
    HArray<float> B;
    {
        // read access on zero sized array, should be okay
        ReadAccess<float> read( B, hostContext );
    }
    HArray<float> C( 10 );
    {
        // read access on undefined array, might give warning
        ReadAccess<float> read( C, hostContext );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( moveTest )
{
    // comparison operator used to sort array by their sizes

    struct less_than_key
    {
        inline bool operator() ( const _HArray& array1, const _HArray& array2 )
        {   
            return array1.size() < array2.size();
        }
    };

    typedef scai::DefaultReal ValueType;

    const IndexType N = 20;   // number of arrays

    const IndexType N_VALUES = 10000;

    const ValueType initVal = 17;

    const ValueType* ptrArray[N];  //  saves pointer data to arrays to check for correct moves

    std::vector<HArray<ValueType> > myData;

    // construct a certain number of different HArrays to store in a vector

    for ( IndexType i = 0; i < N; ++i )
    {
        myData.push_back( HArray<ValueType>( IndexType( N_VALUES - i ), initVal + i ) );
        // save the pointer to the host data 
        ReadAccess<ValueType> rA( myData[i] );
        ptrArray[i] = rA.get();
    }

    // sort the array by sizes, here it gives exactly the reverse vector

    std::sort( myData.begin(), myData.end(), less_than_key() );

    for ( IndexType i = 0; i < N; ++i )
    {
        ReadAccess<ValueType> rA( myData[i] );

        BOOST_CHECK_EQUAL( rA.get(), ptrArray[ N - 1 - i ] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( remoteTest )
{
    ContextPtr hostContext = Context::getHostPtr();
    ContextPtr remoteContext = Context::getContextPtr();

    // Note: in hmemo we can use any value type for HArray

    typedef double ValueType;

    const IndexType Nh = 50;
    const IndexType N = 2 * Nh;

    HArray<ValueType> hostA( N, 5, hostContext );
    HArray<ValueType> remA( N, 2, remoteContext );

    // put single values on remote context

    for ( IndexType i = 0; i < N; i += 2 )
    {
        ReadAccess<ValueType> readA( hostA, hostContext );
        WriteAccess<ValueType> writeA( remA, remoteContext );
        ValueType elem = readA[i];
        writeA.setValue( elem, i );
    }

    ValueType sum = 0;

    {
        ReadAccess<double> readA( remA, hostContext );

        for ( IndexType i = 0; i < N; ++i )
        {
            sum += readA[i];
        }
    }

    {
        // make incarnation of remA invalid
        WriteAccess<ValueType> write( remA, remoteContext );
    }

    BOOST_CHECK_EQUAL( 2 * Nh + 5 * Nh, sum );

    // now we read value from remote context

    for ( IndexType i = 1; i < N; i += 2 )
    {
        ReadAccess<ValueType> readA( remA, remoteContext );
        WriteAccess<ValueType> writeA( hostA, hostContext );
        ValueType elem;
        readA.getValue( elem, i );
        writeA[i] = elem;
    }

    sum = 0;

    {
        ReadAccess<double> readA( hostA, hostContext );

        for ( IndexType i = 0; i < N; ++i )
        {
            sum += readA[i];
        }
    }

    BOOST_CHECK_EQUAL( 2 * Nh + 5 * Nh, sum );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
