/**
 * @file MICContextTest.cpp
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
 * @brief Extensive testing of MICContext.
 * @author: Thomas Brandes
 * @date 13.09.2015
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/hmemo.hpp>
#include <scai/hmemo/mic/MICContext.hpp>

using namespace scai;

using namespace hmemo;

typedef boost::mpl::list<double, float> test_types;

BOOST_AUTO_TEST_SUITE( MICContextTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.MICContextTest" );

template<typename ValueType>
static void scal( int n, ValueType alpha, ValueType* array )
{
    void* arrayPtr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( arrayPtr, n, alpha )
    {
        ValueType* array = static_cast<ValueType*>( arrayPtr );

        #pragma omp parallel for 

        for ( IndexType i = 0; i < n; ++i )
        {
            array[i] *= alpha;
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getContextTest )
{
    // Test will take the default MIC device
    ContextPtr micContext1 = Context::getContextPtr( context::MIC );
    ContextPtr micContext2 = Context::getContextPtr( context::MIC );
    // Two queries for the same context should deliver same pointer
    BOOST_CHECK( micContext1.get() == micContext2.get() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getMemoryTest )
{
    // Test will take the default MIC device
    ContextPtr micContext = Context::getContextPtr( context::MIC );
    MemoryPtr micMemory = micContext->getMemoryPtr();
    BOOST_CHECK( micMemory.get() != NULL );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocateTest )
{
    LAMAArray<int> ctxArray; // default, not allocated at all
    ContextPtr micContext = Context::getContextPtr( context::MIC );
    {
        WriteAccess<int> array( ctxArray, micContext );
        array.resize( 10 );
        // Destructors will be called for array, ctxArray, micContext
    }
    BOOST_CHECK_EQUAL( 10, ctxArray.size() );
    BOOST_CHECK( ctxArray.capacity( micContext ) >= 10 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( useTest )
{
    ContextPtr micContext = Context::getContextPtr( context::MIC );

    const IndexType n = 100;
    const float value = 1.4;
    const float alpha = 0.5;
    LAMAArray<float> vector( n, value );
    {
        WriteAccess<float> micV( vector, micContext );
        SCAI_CONTEXT_ACCESS( micContext );
        scal( n, alpha, micV.get() );
    }
    ReadAccess<float> hostV( vector );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value * alpha, hostV.get()[i] );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

