/**
 * @file LAMAArrayTest.cpp
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
 * @brief Basic tests for LAMA arrays.
 * @author: Thomas Brandes
 * @date 08.07.2015
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/assert.hpp>

#include <scai/hmemo/LAMAArray.hpp>
#include <scai/hmemo/LAMAArrayRef.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

using namespace boost;
using namespace scai::common;
using namespace scai::hmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LAMAArrayTest )

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


BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    ContextPtr host = Context::getHostPtr();

    const IndexType N = 100;

    LAMAArray<float> array1( N );
    LAMAArray<float> array2( N, 5.0f );
   
    BOOST_CHECK_EQUAL( array1.size(), N );
    BOOST_CHECK_EQUAL( array2.size(), N );

    BOOST_CHECK_EQUAL( array1.capacity( host ), N );
    BOOST_CHECK_EQUAL( array2.capacity( host ), N );

    BOOST_CHECK_EQUAL( array1.isValid( host ), false );
    BOOST_CHECK_EQUAL( array2.isValid( host ), true );

    // read access on uninitialized non-empty array
    // gives a warning as there is no valid data

    {
        ReadAccess<float>read( array1, host );
    }

    {
        ReadAccess<float> read( array2, host );
        readTest<float>( read.get(), N, N * 5 );
    }
}

BOOST_AUTO_TEST_CASE( releaseTest )
{
    ContextPtr hostContext = Context::getHostPtr();

    LAMAArray<IndexType> lamaArray; // default, not allocated at all
    ReadAccess<IndexType> readTestAccess( lamaArray, hostContext );
    readTestAccess.release();
    WriteAccess<IndexType> writeAccess( lamaArray, hostContext );
    writeAccess.resize( 10 );

    IndexType* data = writeAccess.get();

    for ( IndexType i = 0; i < 10; i++ )
    {
        data[i] = 3;
    }

    writeAccess.release();
    BOOST_CHECK_THROW( { writeAccess.resize( 20 ); }, Exception );
    BOOST_CHECK_THROW( { writeAccess.get()[0] = static_cast<IndexType> ( 5.0 ); }, Exception );
    ReadAccess<IndexType> readAccess( lamaArray, hostContext );
    const IndexType* readData = readAccess.get();

    for ( IndexType i = 0; i < 5; i++ )
    {
        BOOST_CHECK_EQUAL( 3, readData[i] );
    }

    readAccess.release();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( resizeTest )
{
    ContextPtr contextPtr = Context::getHostPtr();

    LAMAArray<IndexType> lamaArray; // default, not allocated at all
    {
        WriteAccess<IndexType> writeAccess( lamaArray, contextPtr );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( 10 );
        IndexType* data = writeAccess.get();

        for ( IndexType i = 0; i < 10; i++ )
        {
            data[i] = static_cast<IndexType>( 3.0 );
        }
    }
    lamaArray.purge();
    {
        WriteAccess<IndexType> writeAccess( lamaArray, contextPtr );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( 10 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( swapTest )
{
    const IndexType n1 = 10;
    const IndexType n2 = 5;

    LAMAArray<double> arr1( n1, 1 ); 
    LAMAArray<double> arr2( n2, 2 ); 

    arr1.swap( arr2 );

    BOOST_CHECK_EQUAL( arr2.size(), n1 );
    BOOST_CHECK_EQUAL( arr1.size(), n2 );

    ContextPtr contextPtr = Context::getHostPtr();

    {
        ReadAccess<double> read( arr1, contextPtr );
        for ( IndexType i = 0; i < arr1.size(); ++i )
        {
            BOOST_CHECK_EQUAL( 2, read[i] );
        }
    }

    {
        ReadAccess<double> read( arr2, contextPtr );
        for ( IndexType i = 0; i < arr2.size(); ++i )
        {
            BOOST_CHECK_EQUAL( 1, read[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( accessTest )
{
    const IndexType n = 10;
    const double value = 1.0;
    const double value2 = 2.0;

    ContextPtr contextPtr = Context::getHostPtr();

    LAMAArray<double> lamaArray( n, value );
    {
        ReadAccess<double> lamaArrayRAccess( lamaArray, contextPtr );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value, lamaArrayRAccess.get()[i] );
        }

        WriteAccess<double> tmpWriteAccess( lamaArray, contextPtr );
    }
    {
        WriteAccess<double> lamaArrayWAccess( lamaArray, contextPtr );

        for ( IndexType i = 0; i < n; ++i )
        {
            lamaArrayWAccess.get()[i] = value2;
        }

        ReadAccess<double> tmpReadAccess( lamaArray, contextPtr );

        lamaArrayWAccess.release();
        ReadAccess<double> lamaArrayRAccess( lamaArray, contextPtr );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value2, lamaArrayRAccess.get()[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( aliasTest )
{
    const IndexType N = 10;
    const double value = 1.0;
    LAMAArray<double> lamaArray( N, value );

    ContextPtr contextPtr = Context::getHostPtr();

    {
        // read and write access at same time by same thread

        ReadAccess<double> read( lamaArray, contextPtr );
        WriteAccess<double> write( lamaArray, contextPtr );

        for ( IndexType i = 0; i < N; ++i )
        {
            write.get()[i] = 2.0 * read.get()[i];
        }
    }
    {
        // verify that operation was really on the same array

        ReadAccess<double> read( lamaArray, contextPtr );
        readTest<double>( read.get(), N, N * value * 2.0 );
    }
    {
        // with a single write access resize is possilbe

        WriteAccess<double> write( lamaArray, contextPtr );
        write.resize( 2 * N );
    }
    {
        // with read and write at the same time resize throws Exception

        WriteAccess<double> write( lamaArray, contextPtr );
        ReadAccess<double> read( lamaArray, contextPtr );
        BOOST_CHECK_THROW(
        { 
            write.resize( 3 * N );
        }, Exception );
    }
    {
        // read and write access at same time by same thread

        WriteAccess<double> write( lamaArray, contextPtr );
        ReadAccess<double> read( lamaArray, contextPtr );

        // a clear is not possible as it affects the other access
        // Note: clear is the same as resize( 0 )

        BOOST_CHECK_THROW(
        {
            write.clear();
        }, Exception );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( createTest )
{
    LAMAArray<float> A( 10, 1.0f );
    LAMAArray<double> B( 10, 3.1415 );
    LAMAArray<IndexType> C( 10, 5 );

    std::vector<scai::common::scalar::ScalarType> values;

    ContextArray::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); ++i )
    {
        std::cout << "Registered values[" << i << "] = " << values[i] << std::endl;
    }

    BOOST_CHECK( ContextArray::canCreate( scalar::FLOAT ) );
    BOOST_CHECK( ContextArray::canCreate( scalar::DOUBLE ) );
    BOOST_CHECK( ContextArray::canCreate( scalar::INDEX_TYPE ) );
    BOOST_CHECK( !ContextArray::canCreate( scalar::INTERNAL ) );

    ContextArray* ca1 = ContextArray::create( scalar::FLOAT );

    BOOST_REQUIRE( ca1 );

    LAMAArray<double>* da1 = dynamic_cast<LAMAArray<double>*>( ca1 );
    LAMAArray<float>* fa1 = dynamic_cast<LAMAArray<float>*>( ca1 );

    BOOST_CHECK( da1 == NULL );
    BOOST_CHECK( fa1 != NULL );

    BOOST_CHECK_THROW(
    {
        ContextArray* ca1 = ContextArray::create( scalar::INTERNAL );
        ca1->clear();
    }, Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( validTest )
{
    LAMAArray<float> A( 10 );

    ContextPtr hostContext = Context::getHostPtr();

    // Array not allocated at all, should also give some default for validContext

    ContextPtr validContext = A.getValidContext();

    BOOST_CHECK( validContext.get() );

    LAMAArray<float> B;
    {
        // read access on zero sized array, should be okay
        ReadAccess<float> read( B, hostContext );
    }
    LAMAArray<float> C( 10 );
    {
        // read access on undefined array, might give warning
        ReadAccess<float> read( C, hostContext );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
