/**
 * @file LAPACKTest.cpp
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
        auto eps = common::TypeTraits<ValueType>::small();

        // comparison only possible on Host
        ContextPtr host = Context::getHostPtr();
        ReadAccess<ValueType> rA( a, host );

        for ( IndexType i = 0; i < n * n; ++i )
        {
            const auto x1 = common::Math::abs( rA[i] );
            const auto x2 = common::Math::abs( bvalues[i] );
            BOOST_CHECK( common::Math::abs( x1 - x2 ) < eps );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getrif2Test )
{
    typedef DefaultReal ValueType;    // this routine checks correct pivoting and use of border

    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::getrf<ValueType> > getrf;
    static kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::getri<ValueType> > getri;

    SCAI_LOG_INFO( logger, "getrif<" << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << "host" )

    //  matrices a, b used with a * b = identity
    //
    //        x   x   x  x  x  x  x        x   x   x  x   x   x   x
    //        x   0   1  0  0  0  x        x   0   0  0   0   1   x
    //        x   0   0  0  1  0  x        x   1   0  0   0   0   x
    //        x   0   0  0  0  1  x        x   0   0  0   1   0   x
    //        x   0   0  1  0  0  x        x   0   1  0   0   0   x
    //        x   1   0  0  0  0  x        x   0   0  1   0   0   x
    //        x   x   x  x  x  x  x        x   x   x  x   x   x   x

    static IndexType iperm[] = {  1, 3, 4, 2, 0 };  // used to set up a permutation matrix

    const IndexType n = sizeof( iperm ) / sizeof( IndexType );

    const IndexType border = 1;
    const IndexType lda = n + 2 * border;
    const ValueType x   = 19;

    for ( IndexType iorder = 1; iorder < 2; ++iorder )
    {
        // const CBLAS_ORDER order = iorder == 0 ? CblasRowMajor : CblasColMajor;

        // matrices get a border to check for working lda

        HArray<ValueType> a( lda * lda, x, testContext );
        HArray<ValueType> b( lda * lda, x, testContext );

        {
            WriteAccess<ValueType> wA( a, hostContext );
            WriteAccess<ValueType> wB( b, hostContext );

            for ( IndexType i = 0; i < n; ++i )
            {
                for ( IndexType j = 0; j < n; ++j )
                {
                    wA[ lda * ( i + border ) + j + border ] = ValueType( 0 );
                    wB[ lda * ( j + border ) + i + border ] = ValueType( 0 );
                }
            }

            for ( IndexType i = 0; i < n; ++i )
            {
                wA[ lda * ( i + border ) + iperm[i] + border ] = ValueType( 1 );
                wB[ lda * ( iperm[i] + border ) + i + border ] = ValueType( 1 );
            }
        }

        SCAI_LOG_INFO( logger, "Compute inverse of a permutation matrix, size = " << n  )

        // now compute the inverse of the permutation matrix as a submatrix

        HArray<IndexType> permutation( n );
        ContextPtr loc = hostContext;
        {
            WriteAccess<ValueType> wA( a, loc );
            WriteAccess<IndexType> wPermutation( permutation, loc );
            ValueType* aData = wA.get() + border * lda + border;
            getrf[loc->getType()]( n, n, aData, lda, wPermutation.get() );
            getri[loc->getType()]( n, aData, lda, wPermutation.get() );
        }

        // now check for correct results

        {
            ReadAccess<ValueType> rA( a, hostContext );
            ReadAccess<ValueType> rB( b, hostContext );

            for ( IndexType k = 0; k <  lda * lda; ++k )
            {
                BOOST_CHECK_EQUAL( rA[k], rB[k] );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getrifTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    static kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::getrf<ValueType> > getrf;
    static kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::getri<ValueType> > getri;
    SCAI_LOG_INFO( logger, "getrif<" << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << "host" )
    const IndexType n = 3;

    //  matrices a, b used with a * b = identity
    //
    //        2   0  -1       2   1   0
    //       -3   0   2      -4  -2  -1
    //       -2  -1   0       3   2   0

    {
        // CblasRowMajor
        static ValueType avalues[] = { 2.0, 0.0, -1.0, -3.0, 0.0, 2.0, -2.0, -1.0, 0.0};
        static ValueType bvalues[] = { 2.0, 1.0, 0.0, -4.0, -2.0, -1.0, 3.0, 2.0, 0.0};
        HArray<ValueType> a( n * n, avalues, testContext );
        HArray<IndexType> permutation( n );
        ContextPtr loc = Context::getHostPtr();
        {
            WriteAccess<ValueType> wA( a, loc );
            WriteAccess<IndexType> wPermutation( permutation, loc );
            getrf[loc->getType()]( n, n, wA.get(), n, wPermutation.get() );
            getri[loc->getType()]( n, wA.get(), n, wPermutation.get() );
        }

        RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

        {
            ReadAccess<ValueType> rA( a );

            for ( IndexType i = 0; i < n * n; ++i )
            {
                BOOST_CHECK( common::Math::abs( rA[i] - bvalues[i] ) < eps );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( tptrsTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    static kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::tptrs<ValueType> > tptrs;

    auto eps = common::TypeTraits<ValueType>::small();
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
            tptrs[loc->getType()]( CblasLower, common::MatrixOp::NORMAL, CblasNonUnit,
                                   n, 1, rA.get(), wB1.get(), n );

            //  A            X    B
            //  1  1   1     1    3
            //  0  1   1     1    2
            //  0  0   1     1    1
            tptrs[loc->getType()]( CblasUpper, common::MatrixOp::NORMAL, CblasNonUnit,
                                   n, 1, rA.get(), wB2.get(), n );
        }
        {
            ReadAccess<ValueType> rX1( b1 );
            ReadAccess<ValueType> rX2( b2 );

            for ( IndexType i = 0; i < n; ++i )
            {
                BOOST_CHECK( common::Math::abs( rX1[i] - xvalues[i] ) < eps );
                BOOST_CHECK( common::Math::abs( rX2[i] - xvalues[i] ) < eps );
            }
        }
    }
    {
        //Bigger Matrix with a different kind of elements, testing upper/lower for column major
        const IndexType n = 4;
        const IndexType ntri = n * ( n + 1 ) / 2;

        // set up values for A, X and B with A * X = B

        //Matrix A -- lower non-unit column major

        static ValueType avalues1[] =  { 1.0, 2.0, 4.0, 7.0, 3.0, 5.0, 8.0, 6.0, 9.0, 10.0 };
        //Matrix A -- upper non-unit column major
        static ValueType avalues2[] =  { 1.0, 2.0, 5.0, 3.0, 6.0, 8.0, 4.0, 7.0, 9.0, 10.0 };
        //Vector B -- lower non-unit column major
        static ValueType bvalues1[] =  { 1.0, 11.0, 49.0, 146.0 };
        //Vector B -- upper non-unit column major
        static ValueType bvalues2[] =  { 50.0, 94.0, 103.0, 70.0 };
        //Vector X -- for all the same
        static ValueType xvalues[] =   { 1.0, 3.0, 5.0, 7.0 };

        HArray<ValueType> a1( ntri, avalues1, testContext );
        HArray<ValueType> a2( ntri, avalues2, testContext );
        HArray<ValueType> b1( n, bvalues1, testContext );
        HArray<ValueType> b2( n, bvalues2, testContext );
        ContextPtr loc = Context::getHostPtr();
        {
            const IndexType ldb = n;
            const IndexType nrhs = 1;
            ReadAccess<ValueType> rA1( a1, loc );
            WriteAccess<ValueType> wB1( b1, loc );
            ReadAccess<ValueType> rA2( a2, loc );
            WriteAccess<ValueType> wB2( b2, loc );
            //  A            X    B
            //  1  0  0  0   1    1
            //  2  3  0  0   3    11
            //  4  5  6  0   5    49
            //  7  8  9  10  7    146
            tptrs[loc->getType()]( CblasLower, common::MatrixOp::NORMAL, CblasNonUnit,
                                   n, nrhs, rA1.get(), wB1.get(), ldb );
            //  A            X    B
            //  1  2  3  4   1    50
            //  0  5  6  7   3    94
            //  0  0  8  9   5    103
            //  0  0  0  10  7    70
            tptrs[loc->getType()]( CblasUpper, common::MatrixOp::NORMAL, CblasNonUnit,
                                   n, nrhs, rA2.get(), wB2.get(), ldb );
        }
        {
            ReadAccess<ValueType> rX1( b1 );
            ReadAccess<ValueType> rX2( b2 );

            for ( IndexType i = 0; i < n; ++i )
            {
                BOOST_CHECK( common::Math::abs( rX1[i] - xvalues[i] ) < eps );
                BOOST_CHECK( common::Math::abs( rX2[i] - xvalues[i] ) < eps );
            }
        }
    }

    bool isColMajor = true;

    if ( !isColMajor )
    {
        // Test with decimal marked numbers and row-major

        const IndexType n = 3;
        const IndexType ntri = n * ( n + 1 ) / 2;

        // set up values for A, X and B with A * X = B
        //  1.2
        //  2.3  4.5
        //  6.7  8.11  10.13

        static ValueType avalues[]  =  { 1.2, 2.3, 4.5, 6.7, 8.11, 10.13 };
        static ValueType bvalues1[] =  { 2.4, 18.1, 88.38 };
        static ValueType bvalues2[] =  { 31.8, 60.65, 50.65 };
        static ValueType xvalues[]  =  { 2.0, 3.0, 5.0};

        HArray<ValueType> a( ntri, avalues, testContext );
        HArray<ValueType> b1( n, bvalues1, testContext );
        HArray<ValueType> b2( n, bvalues2, testContext );

        ContextPtr loc = Context::getHostPtr();
        {
            IndexType ldb = 1;
            ReadAccess<ValueType> rA( a, loc );
            WriteAccess<ValueType> wB1( b1, loc );
            WriteAccess<ValueType> wB2( b2, loc );
            //  A                  X    B
            //  1.2  0      0      2    2.4
            //  2.3  4.5    0      3    18.1
            //  6.7  8.11   10.13  5    88.38
            tptrs[loc->getType()]( CblasLower, common::MatrixOp::NORMAL, CblasNonUnit,
                                   n, 1, rA.get(), wB1.get(), ldb );

            //  A                  X    B
            //  1.2  2.3   4.5     2    31.8
            //  0    6.7   8.11    3    60.65
            //  0    0     10.13   5    50.65
            tptrs[loc->getType()]( CblasUpper, common::MatrixOp::NORMAL, CblasNonUnit,
                                   n, 1, rA.get(), wB2.get(), ldb );
        }
        {

            ReadAccess<ValueType> rX1( b1 );
            ReadAccess<ValueType> rX2( b2 );

            for ( IndexType i = 0; i < n; ++i )
            {
                BOOST_CHECK( common::Math::abs( rX1[i] - xvalues[i] ) < eps );
                BOOST_CHECK( common::Math::abs( rX2[i] - xvalues[i] ) < eps );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
