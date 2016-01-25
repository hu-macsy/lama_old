/**
 * @file SolutionProxyTest.cpp
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
 * @brief Contains the implementation of the class SolutionProxyTest
 * @author Alexander Büchel, Matthias Makulla
 * @date 03.02.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lamasolver/SolutionProxy.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

typedef Vector VectorType;
typedef SolutionProxy ProxyType;

/* --------------------------------------------------------------------- */

struct SolutionProxyTestConfig
{
    SolutionProxyTestConfig()
    {
        VectorType* vec = new DenseVector<double>( 3, -5.0 );
        mProxy = vec;
    }

    ~SolutionProxyTestConfig()
    {
        VectorType& vec = mProxy.getReference();
        delete &vec;
    }

    ProxyType mProxy;
};

BOOST_FIXTURE_TEST_SUITE( SolutionProxyTest , SolutionProxyTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.SolutionProxyTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testOperators )
{
    BOOST_CHECK_EQUAL( true, mProxy.isDirty() );
    mProxy.setDirty( false );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), ( *mProxy )( 0 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), ( *mProxy )( 1 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), ( *mProxy )( 2 ) );
    BOOST_CHECK_EQUAL( true, mProxy.isDirty() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testSetAndIsDirty )
{
    mProxy.setDirty( true );
    BOOST_CHECK_EQUAL( true, mProxy.isDirty() );
    mProxy.setDirty( false );
    BOOST_CHECK_EQUAL( false, mProxy.isDirty() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testGetConstReference )
{
    mProxy.setDirty( false );
    const VectorType& vec = mProxy.getConstReference();
    BOOST_CHECK_EQUAL( false, mProxy.isDirty() );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec( 0 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec( 1 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec( 2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testGetReference )
{
    mProxy.setDirty( false );
    VectorType& vec = mProxy.getReference();
    BOOST_CHECK_EQUAL( true, mProxy.isDirty() );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec( 0 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec( 1 ) );
    BOOST_CHECK_EQUAL( Scalar( -5.0 ), vec( 2 ) );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
