/**
 * @file HArrayTest.cpp
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

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/HArrayRef.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

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

    HArray<float> array;

    BOOST_CHECK_EQUAL( array.size(), 0 );
    BOOST_CHECK_EQUAL( array.isValid( contextPtr ), false );

    array.resize( N );

    BOOST_CHECK_EQUAL( array.size(), N );
    BOOST_CHECK_EQUAL( array.capacity( contextPtr ), 0 );
    BOOST_CHECK_EQUAL( array.isValid( contextPtr ), false );

    {
        WriteAccess<float> write( array, contextPtr );
        // even if no values are set here the HArray assumes that it is written
    }

    BOOST_CHECK_EQUAL( array.isValid( contextPtr ), true );
    BOOST_CHECK_EQUAL( array.capacity( contextPtr ), N );

    if ( contextPtr->getType() == common::context::Host )
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

BOOST_AUTO_TEST_CASE( capacityTest )
{
    ContextPtr contextPtr = Context::getContextPtr();  // test context
    ContextPtr hostPtr    = Context::getHostPtr();

    SCAI_LOG_INFO( logger, "capacityTest on " << *contextPtr );

    static IndexType N = 10;

    HArray<IndexType> hArray; // default, not allocated at all

    hArray.reserve( contextPtr, N );

    BOOST_CHECK_EQUAL( hArray.size(), 0 );
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), N );

    {
        WriteAccess<IndexType> access( hArray, contextPtr );
        BOOST_CHECK_EQUAL( access.capacity(), N );
        access.reserve( 2 * N );
    }

    BOOST_CHECK_EQUAL( hArray.size(), 0 );
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), 2 * N );

    hArray.clear();  
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), 2 * N );

    hArray.purge();
    BOOST_CHECK_EQUAL( hArray.capacity( contextPtr ), 0 );

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
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( createTest )
{
    using namespace common;

    HArray<float> A( 10, 1.0f );
    HArray<double> B( 10, 3.1415 );
    HArray<IndexType> C( 10, 5 );

    std::vector<scai::common::scalar::ScalarType> values;

    _HArray::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); ++i )
    {
        //std::cout << "Registered values[" << i << "] = " << values[i] << std::endl;
    }

    BOOST_CHECK( _HArray::canCreate( scalar::FLOAT ) );
    BOOST_CHECK( _HArray::canCreate( scalar::DOUBLE ) );
    BOOST_CHECK( _HArray::canCreate( scalar::INDEX_TYPE ) );
    BOOST_CHECK( !_HArray::canCreate( scalar::INTERNAL ) );

    _HArray* ca1 = _HArray::create( scalar::FLOAT );

    BOOST_REQUIRE( ca1 );

    HArray<double>* da1 = dynamic_cast<HArray<double>*>( ca1 );
    HArray<float>* fa1 = dynamic_cast<HArray<float>*>( ca1 );

    BOOST_CHECK( da1 == NULL );
    BOOST_CHECK( fa1 != NULL );

    BOOST_CHECK_THROW(
    {
        _HArray* ca1 = _HArray::create( scalar::INTERNAL );
        ca1->clear();
    }, Exception );
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

BOOST_AUTO_TEST_SUITE_END();
