/**
 * @file LAPACKTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Test class for calling LAPACK routines via the LAMAInterface
 * @author Thomas Brandes
 * @date 10.04.2013
 */

#include <scai/hmemo.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace hmemo;
using namespace kregistry;
using common::TypeTraits;

/* --------------------------------------------------------------------- */

// use of Fixture ContextFix provides the testContext

BOOST_GLOBAL_FIXTURE( ContextFix );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LAPACKTest )

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list<float, double> test_types;

SCAI_LOG_DEF_LOGGER( logger, "Test.LAPACKTest" )

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( inverseTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    const IndexType n = 3;
    // set up values for A and B with A * B = identiy
    static ValueType avalues[] = { 2.0, 0.0, -1.0, -3.0, 0.0, 2.0, -2.0, -1.0, 0.0};
    static ValueType bvalues[] = { 2.0, 1.0, 0.0, -4.0, -2.0, -1.0, 3.0, 2.0, 0.0};
    HArray<ValueType> a( n * n, avalues, testContext );
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::getinv<ValueType> > getinv;
    ContextPtr loc = Context::getContextPtr( getinv.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "getinv<" << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    {
        SCAI_CONTEXT_ACCESS( loc );
        WriteAccess<ValueType> wA( a, loc );
        getinv[loc->getType()]( n, wA.get(), n );
    }
    {
        // comparison only possible on Host
        ContextPtr host = Context::getHostPtr();
        ReadAccess<ValueType> rA( a, host );

        for ( IndexType i = 0; i < n * n; ++i )
        {
            typedef typename TypeTraits<ValueType>::AbsType CompareType;
            CompareType x1 = common::Math::abs( rA[i] );
            CompareType x2 = common::Math::abs( bvalues[i] );
            BOOST_CHECK_CLOSE( x1, x2, 1 );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getrifTest, ValueType, test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    static kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::getrf<ValueType> > getrf;
    static kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::getri<ValueType> > getri;
    SCAI_LOG_INFO( logger, "getrif<" << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << "host" )
    const IndexType n = 3;
// set up values for A and B with A * B = identiy
    {
        //CblasRowMajor
        static ValueType avalues[] = { 2.0, 0.0, -1.0, -3.0, 0.0, 2.0, -2.0, -1.0, 0.0};
        static ValueType bvalues[] = { 2.0, 1.0, 0.0, -4.0, -2.0, -1.0, 3.0, 2.0, 0.0};
        HArray<ValueType> a( n * n, avalues, testContext );
        HArray<IndexType> permutation( n );
        ContextPtr loc = Context::getHostPtr();
        {
            WriteAccess<ValueType> wA( a, loc );
            WriteAccess<IndexType> wPermutation( permutation, loc );
            getrf[loc->getType()]( CblasRowMajor, n, n, wA.get(), n, wPermutation.get() );
            getri[loc->getType()]( CblasRowMajor, n, wA.get(), n, wPermutation.get() );
        }
        {
            ReadAccess<ValueType> rA( a );

            for ( IndexType i = 0; i < n * n; ++i )
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
        HArray<ValueType> a( n * n, avalues, testContext );
        HArray<IndexType> permutation( n );
        ContextPtr loc = Context::getHostPtr();
        {
            WriteAccess<ValueType> wA( a, loc );
            WriteAccess<IndexType> wPermutation( permutation, loc );
            getrf[loc->getType()]( CblasRowMajor, n, n, wA.get(), n, wPermutation.get() );
            getri[loc->getType()]( CblasRowMajor, n, wA.get(), n, wPermutation.get() );
        }
        {
            ReadAccess<ValueType> rA( a );

            for ( IndexType i = 0; i < n * n; ++i )
            {
                BOOST_CHECK_CLOSE( rA[i], bvalues[i], 1 );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( tptrsTest, ValueType, test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    static kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::tptrs<ValueType> > tptrs;
    {
        const IndexType n = 3;
        const IndexType ntri = n * ( n + 1 ) / 2;
// set up values for A, X and B with A * X = B
        static ValueType avalues[] =
        {   1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        static ValueType bvalues1[] = {   1.0, 2.0, 3.0};
        static ValueType bvalues2[] = {   3.0, 2.0, 1.0};
        static ValueType xvalues[]  = {   1.0, 1.0, 1.0};
        HArray<ValueType> a( ntri, avalues, testContext );
        HArray<ValueType> b1( n, bvalues1, testContext );
        HArray<ValueType> b2( n, bvalues2, testContext );
        ContextPtr loc = Context::getHostPtr();
        {
            ReadAccess<ValueType> rA( a, loc );
            WriteAccess<ValueType> wB1( b1, loc );
            WriteAccess<ValueType> wB2( b2, loc );
            //  A            X    B
            //  1  0   0     1    1
            //  1  1   0     1    2
            //  1  1   1     1    3
            tptrs[loc->getType()]( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                                   n, 1, rA.get(), wB1.get(), n );

            //  A            X    B
            //  1  1   1     1    3
            //  0  1   1     1    2
            //  0  0   1     1    1
            tptrs[loc->getType()]( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                                   n, 1, rA.get(), wB2.get(), n );
        }
        {
            ReadAccess<ValueType> rX1( b1 );
            ReadAccess<ValueType> rX2( b2 );

            for ( IndexType i = 0; i < n; ++i )
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
        HArray<ValueType> a1( ntri, avalues1, testContext );
        HArray<ValueType> a2( ntri, avalues2, testContext );
        HArray<ValueType> b1( n, bvalues1, testContext );
        HArray<ValueType> b2( n, bvalues2, testContext );
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
            tptrs[loc->getType()]( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                                   n, 1, rA1.get(), wB1.get(), n );
            //  A            X    B
            //  1  2  3  4   1    50
            //  0  5  6  7   3    94
            //  0  0  8  9   5    103
            //  0  0  0  10  7    70
            tptrs[loc->getType()]( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                                   n, 1, rA2.get(), wB2.get(), n );
        }
        {
            ReadAccess<ValueType> rX1( b1 );
            ReadAccess<ValueType> rX2( b2 );

            for ( IndexType i = 0; i < n; ++i )
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
        HArray<ValueType> a( ntri, avalues, testContext );
        HArray<ValueType> b1( n, bvalues1, testContext );
        HArray<ValueType> b2( n, bvalues2, testContext );
        ContextPtr loc = Context::getHostPtr();
        {
            ReadAccess<ValueType> rA( a, loc );
            WriteAccess<ValueType> wB1( b1, loc );
            WriteAccess<ValueType> wB2( b2, loc );
            //  A                  X    B
            //  1.2  0      0      2    2.4
            //  2.3  4.5    0      3    18.1
            //  6.7  8.11   10.13  5    88.38
            tptrs[loc->getType()]( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                                   n, 1, rA.get(), wB1.get(), n );

            //  A                  X    B
            //  1.2  2.3   4.5     2    31.8
            //  0    6.7   8.11    3    60.65
            //  0    0     10.13   5    50.65
            tptrs[loc->getType()]( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                                   n, 1, rA.get(), wB2.get(), n );
        }
        {
            ReadAccess<ValueType> rX1( b1 );
            ReadAccess<ValueType> rX2( b2 );

            for ( IndexType i = 0; i < n; ++i )
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
