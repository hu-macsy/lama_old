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

#include <memory/LAMAArray.hpp>
#include <memory/LAMAArrayRef.hpp>
#include <memory/HostWriteAccess.hpp>
#include <memory/HostReadAccess.hpp>

using namespace boost;
using namespace common;
using namespace memory;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LAMAArrayTest )

/* --------------------------------------------------------------------- */

template<typename ValueType>
void readTest( const ValueType values[], const size_t N, const ValueType sum )
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
    ContextPtr host = Context::getContext( context::Host );

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
    // should throw exception as there is no valid data

    BOOST_CHECK_THROW(
    {
        HostReadAccess<float>read( array1 );

    }, Exception )

    {
        HostReadAccess<float> read( array2 );
        readTest<float>( read.get(), N, N * 5 );
    }
}

BOOST_AUTO_TEST_CASE( releaseTest )
{
    LAMAArray<IndexType> lamaArray; // default, not allocated at all
    HostReadAccess<IndexType> readTestAccess( lamaArray );
    readTestAccess.release();
    HostWriteAccess<IndexType> writeAccess( lamaArray );
    writeAccess.resize( 10 );

    for ( IndexType i = 0; i < 10; i++ )
    {
        writeAccess[i] = 3;
    }

    writeAccess.release();
    BOOST_CHECK_THROW( { writeAccess.resize( 20 ); }, Exception );
    BOOST_CHECK_THROW( { writeAccess[0] = static_cast<IndexType> ( 5.0 ); }, Exception );
    HostReadAccess<IndexType> readAccess( lamaArray );

    for ( IndexType i = 0; i < 5; i++ )
    {
        BOOST_CHECK_EQUAL( 3, readAccess[i] );
    }

    readAccess.release();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( resizeTest )
{
    LAMAArray<IndexType> lamaArray; // default, not allocated at all
    {
        HostWriteAccess<IndexType> writeAccess( lamaArray );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( 10 );

        for ( IndexType i = 0; i < 10; i++ )
        {
            writeAccess[i] = static_cast<IndexType>( 3.0 );
        }
    }
    lamaArray.purge();
    {
        HostWriteAccess<IndexType> writeAccess( lamaArray );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( 10 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( accessTest )
{
    const IndexType n = 10;
    const double value = 1.0;
    const double value2 = 2.0;
    LAMAArray<double> lamaArray( n, value );
    {
        HostReadAccess<double> lamaArrayRAccess( lamaArray );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value, lamaArrayRAccess[i] );
        }

        HostWriteAccess<double> tmpWriteAccess( lamaArray );
    }
    {
        HostWriteAccess<double> lamaArrayWAccess( lamaArray );

        for ( IndexType i = 0; i < n; ++i )
        {
            lamaArrayWAccess[i] = value2;
        }

        HostReadAccess<double> tmpReadAccess( lamaArray );

        lamaArrayWAccess.release();
        HostReadAccess<double> lamaArrayRAccess( lamaArray );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value2, lamaArrayRAccess[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( createTest )
{
    LAMAArray<float> A( 10, 1.0f );
    LAMAArray<double> B( 10, 3.1415 );
    LAMAArray<IndexType> C( 10, 5 );

    std::vector<ScalarType> values;

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

BOOST_AUTO_TEST_SUITE_END();
