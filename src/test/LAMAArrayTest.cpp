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
 * @brief Contains the implementation of the class LAMAArrayTest.
 * @author: Alexander Büchel, Lauretta Schubert
 * @date 18.04.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <test/TestMacros.hpp>

#include <lama/LAMAArray.hpp>
#include <lama/HostWriteAccess.hpp>
#include <lama/HostReadAccess.hpp>

using namespace boost;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LAMAArrayTest )

/* --------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( logger, "Test.LAMAArrayTest" )

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( releaseTest )
{
    LAMAArray<IndexType> lamaArray; // default, not allocated at all
    LAMA_LOG_INFO( logger, "make read test access on empty array\n" );
    HostReadAccess<IndexType> readTestAccess( lamaArray );
    LAMA_LOG_INFO( logger, "release read on empty array\n" );
    readTestAccess.release();
    LAMA_LOG_INFO( logger, "make write test access on empty array\n" );
    HostWriteAccess<IndexType> writeAccess( lamaArray );
    writeAccess.resize( 10 );

    for ( IndexType i = 0; i < 10; i++ )
    {
        writeAccess[i] = 3;
    }

    writeAccess.release();
    LAMA_CHECK_THROW( { writeAccess.resize( 20 ); }, Exception );
    LAMA_CHECK_THROW( { writeAccess[0] = static_cast<IndexType> ( 5.0 ); }, Exception );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( accessTest , ValueType, test_types )
{
    const IndexType n = 10;
    const ValueType value = 1.0;
    const ValueType value2 = 2.0;
    LAMAArray<ValueType> lamaArray( n, value );
    {
        HostReadAccess<ValueType> lamaArrayRAccess( lamaArray );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value, lamaArrayRAccess[i] );
        }

        LAMA_CHECK_THROW(
        {   HostWriteAccess<ValueType> tmpWriteAccess( lamaArray );}, Exception );
    }
    {
        HostWriteAccess<ValueType> lamaArrayWAccess( lamaArray );

        for ( IndexType i = 0; i < n; ++i )
        {
            lamaArrayWAccess[i] = value2;
        }

        LAMA_CHECK_THROW(
        {   HostReadAccess<ValueType> tmpReadAccess( lamaArray );}, Exception );
        lamaArrayWAccess.release();
        HostReadAccess<ValueType> lamaArrayRAccess( lamaArray );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value2, lamaArrayRAccess[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( refTest, ValueType, test_types )
{
    const IndexType n = 10;
    const ValueType value = 3.5;
    ValueType myData[10] =
    {   1, 2, 3, 4, 5, 5, 4, 3, 2, 1};
    const ValueType* myData1 = myData;
    {
        LAMAArrayRef<ValueType> lamaArray( myData, 10 );
        HostWriteAccess<ValueType> lamaArrayWAccess( lamaArray );
        // resize of a LAMA array with referenced data is not possible
        LAMA_CHECK_THROW(
        {   lamaArrayWAccess.resize( 20 );}, Exception );

        for ( IndexType i = 0; i < n; ++i )
        {
            lamaArrayWAccess[i] = value;
        }

        // destructor of lamaArray should write data back
    }

// we used the LAMA array for computations on myData, verify it

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value, myData[i] );
    }

    {
        // this will create a LAMA array with a const reference,
        LAMAArrayRef<ValueType> lamaArray( myData1, 10 );
        BOOST_CHECK_EQUAL( 10, lamaArray.size() );
        // Write access should not be allowed
        LAMA_CHECK_THROW(
        {
            HostWriteAccess<ValueType> lamaArrayWAccess( lamaArray );
        }
        , Exception );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    LAMAArray<IndexType> lamaArray;
    LAMA_WRITEAT_TEST( lamaArray );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
