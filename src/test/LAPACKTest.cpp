/**
 * @file LAPACKTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Test class for calling LAPACK routines via the LAMAInterface
 * @author Thomas Brandes
 * @date 10.04.2013
 * $Id$
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/ContextFactory.hpp>
#include <lama/LAMAArray.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/WriteAccess.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/ContextAccess.hpp>

#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

// extern bool base_test_case;
// extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LAPACKTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.LAPACKTest" )

typedef boost::mpl::list<float,double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( inverseTest, T, test_types ) {

    typedef T ValueType;

    const IndexType n = 3;

    // set up values for A and B with A * B = identiy

    static ValueType avalues[] = { 2.0, 0.0, -1.0, -3.0, 0.0, 2.0, -2.0, -1.0, 0.0 };
    static ValueType bvalues[] = { 2.0, 1.0,  0.0, -4.0, -2.0, -1.0, 3.0, 2.0, 0.0 };

    LAMAArray<ValueType> a( n * n, avalues );
    LAMAArray<IndexType> permutation( n );

    ContextPtr loc = ContextFactory::getContext( Context::Host );

    {
        WriteAccess<ValueType> wA( a, loc );

        LAMA_INTERFACE_FN_T( getinv, loc, BLAS, LAPACK, ValueType )

        getinv( n, wA.get(), n );
    }
    
    {
        HostReadAccess<ValueType> rA( a );

        for ( int i = 0; i < n * n ; ++i )
        {
            BOOST_CHECK_CLOSE( rA[i], bvalues[i], 1 );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getrifTest, T, test_types ) {

    typedef T ValueType;

    const IndexType n = 3;

    // set up values for A and B with A * B = identiy

    static ValueType avalues[] = { 2.0, 0.0, -1.0, -3.0, 0.0, 2.0, -2.0, -1.0, 0.0 };
    static ValueType bvalues[] = { 2.0, 1.0,  0.0, -4.0, -2.0, -1.0, 3.0, 2.0, 0.0 };

    LAMAArray<ValueType> a( n * n, avalues );
    LAMAArray<IndexType> permutation( n );

    ContextPtr loc = ContextFactory::getContext( Context::Host );

    {
        WriteAccess<ValueType> wA( a, loc );
        WriteAccess<IndexType> wPermutation( permutation, loc );

        LAMA_INTERFACE_FN_T( getrf, loc, BLAS, LAPACK, ValueType )

        int error = getrf( CblasRowMajor, n, n, wA.get(), n, wPermutation.get() );
    
        BOOST_CHECK_EQUAL( 0, error );

        LAMA_INTERFACE_FN_T( getri, loc, BLAS, LAPACK, ValueType )

        error = getri( CblasRowMajor, n, wA.get(), n, wPermutation.get() );
    
        BOOST_CHECK_EQUAL( 0, error );
    }
    
    {
        HostReadAccess<ValueType> rA( a );

        for ( int i = 0; i < n * n ; ++i )
        {
            BOOST_CHECK_CLOSE( rA[i], bvalues[i], 1 );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( tptrsTest, T, test_types ) {

    typedef T ValueType;

    const IndexType n = 3;
    const IndexType ntri = n * ( n + 1 ) / 2;

    // set up values for A, X and B with A * X = B

    static ValueType avalues[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    static ValueType bvalues1[] = { 1.0, 2.0, 3.0 };
    static ValueType bvalues2[] = { 3.0, 2.0, 1.0 };
    static ValueType xvalues[] = { 1.0, 1.0, 1.0 };

    LAMAArray<ValueType> a( ntri, avalues );
    LAMAArray<ValueType> b1( n, bvalues1 );
    LAMAArray<ValueType> b2( n, bvalues2 );

    ContextPtr loc = ContextFactory::getContext( Context::Host );

    {
        ReadAccess<ValueType> rA( a, loc );
        WriteAccess<ValueType> wB1( b1, loc );
        WriteAccess<ValueType> wB2( b2, loc );

        LAMA_INTERFACE_FN_T( tptrs, loc, BLAS, LAPACK, ValueType )

        //  A            X    B
        //  1  0   0     1    1
        //  1  1   0     1    2
        //  1  1   1     1    3

        int error = tptrs( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                           n, 1, rA.get(), wB1.get(), n );
    
        BOOST_CHECK_EQUAL( 0, error );

        //  A            X    B
        //  1  1   1     1    3
        //  0  1   1     1    2
        //  0  0   1     1    1

        error = tptrs( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                       n, 1, rA.get(), wB2.get(), n );
    
        BOOST_CHECK_EQUAL( 0, error );
    }
    
    {
        HostReadAccess<ValueType> rX1( b1 );
        HostReadAccess<ValueType> rX2( b2 );

        for ( int i = 0; i < n ; ++i )
        {
            BOOST_CHECK_CLOSE( rX1[i], xvalues[i], 1 );
            BOOST_CHECK_CLOSE( rX2[i], xvalues[i], 1 );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
