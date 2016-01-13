/**
 * @file LAPACKTest.cpp
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
 * @brief Test class for calling LAPACK routines via the LAMAInterface
 * @author Thomas Brandes
 * @date 10.04.2013
 * @since 1.0.0
 */

#include <scai/hmemo.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/LAMAKernel.hpp>
#include <scai/lama/LArray.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;


// extern bool base_test_case;
// extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LAPACKTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.LAPACKTest" )

typedef boost::mpl::list<float, double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( inverseTest, ValueType, test_types )
{
    const IndexType n = 3;
// set up values for A and B with A * B = identiy
    static ValueType avalues[] =
    {   2.0, 0.0, -1.0, -3.0, 0.0, 2.0, -2.0, -1.0, 0.0};
    static ValueType bvalues[] =
    {   2.0, 1.0, 0.0, -4.0, -2.0, -1.0, 3.0, 2.0, 0.0};
    LArray<ValueType> a( n * n, avalues );
    LArray<IndexType> permutation( n );
    ContextPtr loc = Context::getHostPtr();
    {
        WriteAccess<ValueType> wA( a, loc );
        LAMAKernel<scai::blaskernel::BLASKernelTrait::getinv<ValueType> > getinv;
        getinv[loc]( n, wA.get(), n );
    }
    {
        ReadAccess<ValueType> rA( a );

        for ( int i = 0; i < n * n; ++i )
        {
            BOOST_CHECK_CLOSE( rA[i], bvalues[i], 1 );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getrifTest, ValueType, test_types )
{
    static LAMAKernel<scai::blaskernel::BLASKernelTrait::getrf<ValueType> > getrf;
    static LAMAKernel<scai::blaskernel::BLASKernelTrait::getri<ValueType> > getri;

    const IndexType n = 3;
// set up values for A and B with A * B = identiy
    {
        //CblasRowMajor
        static ValueType avalues[] =
        {   2.0, 0.0, -1.0, -3.0, 0.0, 2.0, -2.0, -1.0, 0.0};
        static ValueType bvalues[] =
        {   2.0, 1.0, 0.0, -4.0, -2.0, -1.0, 3.0, 2.0, 0.0};
        LArray<ValueType> a( n * n, avalues );
        LArray<IndexType> permutation( n );
        ContextPtr loc = Context::getHostPtr();
        {
            WriteAccess<ValueType> wA( a, loc );
            WriteAccess<IndexType> wPermutation( permutation, loc );
            int error = getrf[loc]( CblasRowMajor, n, n, wA.get(), n, wPermutation.get() );
            BOOST_CHECK_EQUAL( 0, error );
            error = getri[loc]( CblasRowMajor, n, wA.get(), n, wPermutation.get() );
            BOOST_CHECK_EQUAL( 0, error );
        }
        {
            ReadAccess<ValueType> rA( a );

            for ( int i = 0; i < n * n; ++i )
            {
                BOOST_CHECK_CLOSE( rA[i], bvalues[i], 1 );
            }
        }
    }
    {
        //CblasColumnMajor
        static ValueType avalues[] =
        {   2.0, -3.0, -2.0, 0.0, 0.0, -1.0, -1.0, 2.0, 0.0};
        static ValueType bvalues[] =
        {   2.0, -4.0, 3.0, 1.0, -2.0, 2.0, 0.0, -1.0, 0.0};
        LArray<ValueType> a( n * n, avalues );
        LArray<IndexType> permutation( n );
        ContextPtr loc = Context::getHostPtr();
        {
            WriteAccess<ValueType> wA( a, loc );
            WriteAccess<IndexType> wPermutation( permutation, loc );
            int error = getrf[loc]( CblasRowMajor, n, n, wA.get(), n, wPermutation.get() );
            BOOST_CHECK_EQUAL( 0, error );
            error = getri[loc]( CblasRowMajor, n, wA.get(), n, wPermutation.get() );
            BOOST_CHECK_EQUAL( 0, error );
        }
        {
            ReadAccess<ValueType> rA( a );

            for ( int i = 0; i < n * n; ++i )
            {
                BOOST_CHECK_CLOSE( rA[i], bvalues[i], 1 );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( tptrsTest, ValueType, test_types )
{
    static LAMAKernel<scai::blaskernel::BLASKernelTrait::tptrs<ValueType> > tptrs;

    {
        const IndexType n = 3;
        const IndexType ntri = n * ( n + 1 ) / 2;
// set up values for A, X and B with A * X = B
        static ValueType avalues[] =
        {   1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        static ValueType bvalues1[] =
        {   1.0, 2.0, 3.0};
        static ValueType bvalues2[] =
        {   3.0, 2.0, 1.0};
        static ValueType xvalues[] =
        {   1.0, 1.0, 1.0};
        LArray<ValueType> a( ntri, avalues );
        LArray<ValueType> b1( n, bvalues1 );
        LArray<ValueType> b2( n, bvalues2 );
        ContextPtr loc = Context::getHostPtr();
        {
            ReadAccess<ValueType> rA( a, loc );
            WriteAccess<ValueType> wB1( b1, loc );
            WriteAccess<ValueType> wB2( b2, loc );
            //  A            X    B
            //  1  0   0     1    1
            //  1  1   0     1    2
            //  1  1   1     1    3
            int error = tptrs[loc]( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                               n, 1, rA.get(), wB1.get(), n );
            BOOST_CHECK_EQUAL( 0, error );
            //  A            X    B
            //  1  1   1     1    3
            //  0  1   1     1    2
            //  0  0   1     1    1
            error = tptrs[loc]( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                           n, 1, rA.get(), wB2.get(), n );
            BOOST_CHECK_EQUAL( 0, error );
        }
        {
            ReadAccess<ValueType> rX1( b1 );
            ReadAccess<ValueType> rX2( b2 );

            for ( int i = 0; i < n; ++i )
            {
//            printf("Lower:");
                BOOST_CHECK_CLOSE( rX1[i], xvalues[i], 1 );
//            printf("\nUpper:");
                BOOST_CHECK_CLOSE( rX2[i], xvalues[i], 1 );
//            printf("\n");
            }
        }
    }
    {
        //Bigger Matrix with a different kind of elements, testing upper/lower for column major
        const IndexType n = 4;
        const IndexType ntri = n * ( n + 1 ) / 2;
// set up values for A, X and B with A * X = B
        //Matrix A -- lower non-unit column major
        static ValueType avalues1[] =
        {   1.0, 2.0, 4.0, 7.0, 3.0, 5.0, 8.0, 6.0, 9.0, 10.0};
        //Matrix A -- upper non-unit column major
        static ValueType avalues2[] =
        {   1.0, 2.0, 5.0, 3.0, 6.0, 8.0, 4.0, 7.0, 9.0, 10.0};
        //Vector B -- lower non-unit column major
        static ValueType bvalues1[] =
        {   1.0, 11.0, 49.0, 146.0};
        //Vector B -- upper non-unit column major
        static ValueType bvalues2[] =
        {   50.0, 94.0, 103.0, 70.0};
        //Vector X -- for all the same
        static ValueType xvalues[] =
        {   1.0, 3.0, 5.0, 7.0};
        LArray<ValueType> a1( ntri, avalues1 );
        LArray<ValueType> a2( ntri, avalues2 );
        LArray<ValueType> b1( n, bvalues1 );
        LArray<ValueType> b2( n, bvalues2 );
        ContextPtr loc = Context::getHostPtr();
        {
            ReadAccess<ValueType> rA1( a1, loc );
            WriteAccess<ValueType> wB1( b1, loc );
            ReadAccess<ValueType> rA2( a2, loc );
            WriteAccess<ValueType> wB2( b2, loc );
            //  A            X    B
            //  1  0  0  0   1    1
            //  2  3  0  0   3    11
            //  4  5  6  0   5    49
            //  7  8  9  10  7    146
            int error1 = tptrs[loc]( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                                     n, 1, rA1.get(), wB1.get(), n );
            //  A            X    B
            //  1  2  3  4   1    50
            //  0  5  6  7   3    94
            //  0  0  8  9   5    103
            //  0  0  0  10  7    70
            int error2 = tptrs[loc]( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                                     n, 1, rA2.get(), wB2.get(), n );
            BOOST_CHECK_EQUAL( 0, error1 );
            BOOST_CHECK_EQUAL( 0, error2 );
        }
        {
            ReadAccess<ValueType> rX1( b1 );
            ReadAccess<ValueType> rX2( b2 );

            for ( int i = 0; i < n; ++i )
            {
                BOOST_CHECK_CLOSE( rX1[i], xvalues[i], 1 );
                BOOST_CHECK_CLOSE( rX2[i], xvalues[i], 1 );
            }
        }
    }
//Test not working for MKL
    /*{ //Bigger Matrix with a different kind of elements, testing upper/lower for row major
        const IndexType n = 4;
        const IndexType ntri = n * ( n + 1 ) / 2;

    // set up values for A, X and B with A * X = B

        //Matrix A -- lower/upper non-unit row major
        static ValueType avalues3[] =
        {   1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
        //Vector B -- lower non-unit row major
        static ValueType bvalues3[] =
        {   1.0, 11.0, 49.0, 146.0};
        //Vector B -- upper non-unit row major
        static ValueType bvalues4[] =
        {   50.0, 94.0, 103.0, 70.0};
        //Vector X -- for all the same
        static ValueType xvalues[] =
        {   1.0, 3.0, 5.0, 7.0};
        printf("######Test1\n");
        LArray<ValueType> a3( ntri, avalues3 );
        LArray<ValueType> b3( n, bvalues3 );
        LArray<ValueType> b4( n, bvalues4 );

        ContextPtr loc = Context::getContextPtr( context::Host );

        {
            ReadAccess<ValueType> rA3( a3, loc );
            WriteAccess<ValueType> wB3( b3, loc );
            WriteAccess<ValueType> wB4( b4, loc );

            LAMA_INTERFACE_FN_T( tptrs, loc, BLAS, LAPACK, ValueType )

            //  A            X    B
            //  1  0  0  0   1    1
            //  2  3  0  0   3    11
            //  4  5  6  0   5    49
            //  7  8  9  10  7    146

            int error3 = tptrs( CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                            n, 1, rA3.get(), wB3.get(), n );
            printf("######Test2\n");
            //  A            X    B
            //  1  2  3  4   1    50
            //  0  5  6  7   3    94
            //  0  0  8  9   5    103
            //  0  0  0  10  7    70

            int error4 = tptrs( CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                            n, 1, rA3.get(), wB4.get(), n );
            printf("######Test3\n");
            BOOST_CHECK_EQUAL( 0, error3 );
            BOOST_CHECK_EQUAL( 0, error4 );
            printf("######Test4\n");
        }

        {
            ReadAccess<ValueType> rX3( b3 );
            ReadAccess<ValueType> rX4( b4 );

            for ( int i = 0; i < n; ++i )
            {
                BOOST_CHECK_CLOSE( rX3[i], xvalues[i], 1 );
                BOOST_CHECK_CLOSE( rX4[i], xvalues[i], 1 );
            }
        }
    }*/
    {
        //Test with decimal marked numbers
        const IndexType n = 3;
        const IndexType ntri = n * ( n + 1 ) / 2;
// set up values for A, X and B with A * X = B
        static ValueType avalues[] =
        {   1.2, 2.3, 6.7, 4.5, 8.11, 10.13};
        static ValueType bvalues1[] =
        {   2.4, 18.1, 88.38};
        static ValueType bvalues2[] =
        {   31.8, 60.65, 50.65};
        static ValueType xvalues[] =
        {   2.0, 3.0, 5.0};
        LArray<ValueType> a( ntri, avalues );
        LArray<ValueType> b1( n, bvalues1 );
        LArray<ValueType> b2( n, bvalues2 );
        ContextPtr loc = Context::getHostPtr();
        {
            ReadAccess<ValueType> rA( a, loc );
            WriteAccess<ValueType> wB1( b1, loc );
            WriteAccess<ValueType> wB2( b2, loc );
            //  A                  X    B
            //  1.2  0      0      2    2.4
            //  2.3  4.5    0      3    18.1
            //  6.7  8.11   10.13  5    88.38
            int error = tptrs[loc]( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                               n, 1, rA.get(), wB1.get(), n );
            BOOST_CHECK_EQUAL( 0, error );
            //  A                  X    B
            //  1.2  2.3   4.5     2    31.8
            //  0    6.7   8.11    3    60.65
            //  0    0     10.13   5    50.65
            error = tptrs[loc]( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                           n, 1, rA.get(), wB2.get(), n );
            BOOST_CHECK_EQUAL( 0, error );
        }
        {
            ReadAccess<ValueType> rX1( b1 );
            ReadAccess<ValueType> rX2( b2 );

            for ( int i = 0; i < n; ++i )
            {
//            printf("Lower:");
                BOOST_CHECK_CLOSE( rX1[i], xvalues[i], 1 );
//            printf("\nUpper:");
                BOOST_CHECK_CLOSE( rX2[i], xvalues[i], 1 );
//            printf("\n");
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
