/**
 * @file BLAS1Test.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Contains tests for the blas1 methods.
 * @author Bea Hornef
 * @date 16.03.2015
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/hmemo.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

#include <scai/blaskernel/test/TestMacros.hpp>
#include <scai/blaskernel/openmp/OpenMPBLAS1.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace blaskernel;
using namespace hmemo;
using common::TypeTraits;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( BLAS1Test )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.BLAS1Test" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( asumOpenMPTest )
{
    // direct use of OpenMP::asum without registry

    typedef SCAI_TEST_TYPE ValueType;

    ValueType values[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
    const IndexType nValues = sizeof( values ) / sizeof( ValueType );
    const IndexType incX1 = 1;
    const IndexType incX2 = 2;
    const ValueType result1 = 21.0;
    const ValueType result2 = 9.0;
    HArray<ValueType> AValues( nValues, values );
    {
        ReadAccess<ValueType> rAValues( AValues );

        ValueType sum = OpenMPBLAS1::asum( 0, rAValues.get(), incX1 );
        BOOST_CHECK_EQUAL( sum, ValueType( 0 ) );

        sum = OpenMPBLAS1::asum( 3, rAValues.get(), 0 );
        BOOST_CHECK_EQUAL( sum, ValueType( 0 ) );

        sum = OpenMPBLAS1::asum( 6, rAValues.get(), incX1 );
        BOOST_CHECK_EQUAL( sum, result1 );

        sum = OpenMPBLAS1::asum( 3, rAValues.get(), incX2 );
        BOOST_CHECK_EQUAL( sum, result2 );
    }

}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( asumTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;   // initialization/free  of context not here
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::asum<ValueType> > asum;
    ContextPtr loc = Context::getContextPtr( asum.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "asum< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    {
        ValueType values[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        const IndexType nValues = sizeof( values ) / sizeof( ValueType );
        const IndexType incX1 = 1;
        const IndexType incX2 = 2;
        const ValueType result1 = 21.0;
        const ValueType result2 = 9.0;
        HArray<ValueType> AValues( nValues, values, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            // n <= 0
            ValueType sum = asum[loc->getType()]( 0, rAValues.get(), incX1 );
            BOOST_CHECK_EQUAL( sum, 0.0 );
            // incX <= 0
            sum = asum[loc->getType()]( 3, rAValues.get(), 0 );
            BOOST_CHECK_EQUAL( sum, 0.0 );
            // std::cout << "test 1 (incX = 1)" << std::endl;
            sum = asum[loc->getType()]( 6, rAValues.get(), incX1 );
            BOOST_CHECK_EQUAL( sum, result1 );
            // std::cout << "test 2 (incX = 2)" << std::endl;
            sum = asum[loc->getType()]( 3, rAValues.get(), incX2 );
            BOOST_CHECK_EQUAL( sum, result2 );
        }
    }
} // asumTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( axpyTest, ValueType, scai_array_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::axpy<ValueType> > axpy;
    ContextPtr loc = Context::getContextPtr( axpy.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "axpy< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    IndexType lenx = 6;
    IndexType leny = 9;

    // check with n <= 0
    {
        int x_int[] = { 1, 2, -3, 4, 5, -6 };
        int y_int[] = { 1, 2, -3, 4, 5, -6, 7, 8, -9 };

        ValueType* x = new ValueType[lenx];
        ValueType* y = new ValueType[leny];

        std::copy( x_int, x_int + lenx, x );
        std::copy( y_int, y_int + leny, y );

        const IndexType incX = 2;
        const IndexType incY = 3;
        HArray<ValueType> Ax( lenx, x, testContext );
        HArray<ValueType> Ay( leny, y, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            axpy[loc->getType()]( 0, 5.0, wAx.get(), incX, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );

            for ( IndexType i = 0; i < leny; ++i )
            {
                BOOST_CHECK_EQUAL( y[i], rAy[i] );
            }
        }
    }
    // check with incX <= 0 and incY <= 0
    {
        int x_int[] = { 1, 2, -3, 4, 5, -6 };
        int y_int[] = { 1, 2, -3, 4, 5, -6, 7, 8, -9 };

        ValueType* x = new ValueType[lenx];
        ValueType* y = new ValueType[leny];

        std::copy( x_int, x_int + lenx, x );
        std::copy( y_int, y_int + leny, y );

        HArray<ValueType> Ax( lenx, x, testContext );
        HArray<ValueType> Ay( leny, y, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            axpy[loc->getType()]( 3, 5.0, wAx.get(), 0, wAy.get(), 0 );
        }
        {
            ReadAccess<ValueType> rAy( Ay );

            for ( IndexType i = 0; i < leny; ++i )
            {
                BOOST_CHECK_EQUAL( y[i], rAy[i] );
            }
        }
    }
    // check with n > 0 and incX > 0 and incY > 0
    {
        int x_int[] = { 1, 2, -3, 4, 5, -6 };
        int y_int[] = { 1, 2, -3, 4, 5, -6, 7, 8, -9 };
        int yResult_int[] = { 6, 2, -3, -11, 5, -6, 32, 8, -9 };

        ValueType* x = new ValueType[lenx];
        ValueType* y = new ValueType[leny];
        ValueType* yResult = new ValueType[leny];

        std::copy( x_int, x_int + lenx, x );
        std::copy( y_int, y_int + leny, y );
        std::copy( yResult_int, yResult_int + leny, yResult );

        const IndexType incX = 2;
        const IndexType incY = 3;
        HArray<ValueType> Ax( testContext );
        Ax.setRawData( lenx, x );
        HArray<ValueType> Ay( testContext );
        Ay.setRawData( leny, y );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            axpy[loc->getType()]( 3, 5.0, wAx.get(), incX, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );

            for ( IndexType i = 0; i < leny; ++i )
            {
                BOOST_CHECK_EQUAL( yResult[i], rAy[i] );
            }
        }
    }
} // axpyTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( copyTest, ValueType, scai_array_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::copy<ValueType> > copy;
    ContextPtr loc = Context::getContextPtr( copy.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "copy< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    int lenx = 6;
    int leny = 9;

    // check with n <= 0
    {
        int x_int[] = { 1, 2, -3, 4, 5, -6 };
        int y_int[] = { -9, 8, -7, 6, 5, -4, 3, 2, -1 };

        ValueType* x = new ValueType[lenx];
        ValueType* y = new ValueType[leny];

        std::copy( x_int, x_int + lenx, x );
        std::copy( y_int, y_int + leny, y );

        const IndexType incX = 2;
        const IndexType incY = 3;
        HArray<ValueType> Ax( testContext );
        Ax.setRawData( 6, x );
        HArray<ValueType> Ay( testContext );
        Ay.setRawData( 9, y );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            copy[loc->getType()]( 0, wAx.get(), incX, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );

            for ( int i = 0; i < 9; ++i )
            {
                BOOST_CHECK_EQUAL( y[i], rAy[i] );
            }
        }
    }
    // check with incX <= 0 and incY <= 0
    {
        int x_int[] = { 1, 2, -3, 4, 5, -6 };
        int y_int[] = { -9, 8, -7, 6, 5, -4, 3, 2, -1 };

        ValueType* x = new ValueType[lenx];
        ValueType* y = new ValueType[leny];

        std::copy( x_int, x_int + lenx, x );
        std::copy( y_int, y_int + leny, y );

        const IndexType incX = 0;
        const IndexType incY = 0;

        HArray<ValueType> Ax( 6, x, testContext );
        HArray<ValueType> Ay( 9, y, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            copy[loc->getType()]( 3, wAx.get(), incX, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );

            for ( int i = 0; i < 9; ++i )
            {
                BOOST_CHECK_EQUAL( y[i], rAy[i] );
            }
        }
    }
    // check with n > 0 and incX > 0 and incY > 0
    {
        int x_int[] = { 1, 2, -3, 4, 5, -6 };
        int y_int[] = { -9, 8, -7, 6, 5, -4, 3, 2, -1 };
        int yResult_int[] = { 1, 8, -7, -3, 5, -4, 5, 2, -1 };

        ValueType* x = new ValueType[lenx];
        ValueType* y = new ValueType[leny];
        ValueType* yResult = new ValueType[leny];

        std::copy( x_int, x_int + lenx, x );
        std::copy( y_int, y_int + leny, y );
        std::copy( yResult_int, yResult_int + leny, yResult );

        const IndexType incX = 2;
        const IndexType incY = 3;
        HArray<ValueType> Ax( 6, x, testContext );
        HArray<ValueType> Ay( 9, y, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            copy[loc->getType()]( 3, wAx.get(), incX, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );

            for ( int i = 0; i < 9; ++i )
            {
                BOOST_CHECK_EQUAL( yResult[i], rAy[i] );
            }
        }
    }
} // copyTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( dotTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::dot<ValueType> > dot;
    ContextPtr loc = Context::getContextPtr( dot.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "dot< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    {
        ValueType x[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0, 7.0, 8.0, -9.0 };
        const IndexType incX = 2;
        const IndexType incY = 3;
        HArray<ValueType> Ax( 6, x, testContext );
        HArray<ValueType> Ay( 9, y, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            // n <= 0
            ValueType result = dot[loc->getType()]( 0, wAx.get(), incX, wAy.get(), incY );
            BOOST_CHECK_EQUAL( result, 0.0 );
            // incX <= 0 and incY <= 0
            result = dot[loc->getType()]( 3, wAx.get(), 0, wAy.get(), 0 );
            BOOST_CHECK_EQUAL( result, 0.0 );
            // n > 0 and incX > 0 and incY > 0
            result = dot[loc->getType()]( 3, wAx.get(), incX, wAy.get(), incY );
            BOOST_CHECK_EQUAL( result, 24.0 );
        }
    }
} // dotTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( iamaxTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::iamax<ValueType> > iamax;
    ContextPtr loc = Context::getContextPtr( iamax.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "iamax< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    {
        ValueType values[] =  { 1.0, 2.0, 3.0, 6.0, 5.0, 6.0 };
        const IndexType nValues = sizeof( values ) / sizeof( ValueType );
        const IndexType incX1 = 1;
        const IndexType incX2 = 2; // { 1, 3, 5}
        const IndexType result1 = 3;
        const IndexType result2 = 2;
        const IndexType zero = 0;
        HArray<ValueType> AValues( nValues, values, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            // n <= 0
            IndexType smallestIndexOfMax = iamax[loc->getType()]( 0, rAValues.get(), incX1 );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, zero );
            // incX <= 0
            smallestIndexOfMax = iamax[loc->getType()]( nValues / incX1, rAValues.get(), 0 );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, zero );
            // n > 0 and incX > 0
            smallestIndexOfMax = iamax[loc->getType()]( nValues / incX1, rAValues.get(), incX1 );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, result1 );
            smallestIndexOfMax = iamax[loc->getType()]( nValues / incX2, rAValues.get(), incX2 );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, result2 );
        }
    }
} // iamaxTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( nrm2Test, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::nrm2<ValueType> > nrm2;
    ContextPtr loc = Context::getContextPtr( nrm2.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "nrm2< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    {
        ValueType values[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        const IndexType nValues = sizeof( values ) / sizeof( ValueType );
        const IndexType incX1 = 1;
        const IndexType incX2 = 2;
        const ValueType result1 = 91.0;
        const ValueType result2 = 35.0;
        HArray<ValueType> AValues( nValues, values, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            // n <= 0, but be careful if IndexType is unsigned
            ValueType euclideanNorm = nrm2[loc->getType()]( 0, rAValues.get(), incX1 );
            BOOST_CHECK_EQUAL( euclideanNorm, ValueType( 0 ) );
            // incX <= 0
            euclideanNorm = nrm2[loc->getType()]( 5, rAValues.get(), 0 );
            BOOST_CHECK_EQUAL( euclideanNorm, ValueType( 0 ) );
            // n > 0 and incX > 0
            euclideanNorm = nrm2[loc->getType()]( nValues / incX1, rAValues.get(), incX1 );
            SCAI_CHECK_CLOSE( euclideanNorm, common::Math::sqrt( result1 ), 1e-4 );
            euclideanNorm = nrm2[loc->getType()]( nValues / incX2, rAValues.get(), incX2 );
            SCAI_CHECK_CLOSE( euclideanNorm, common::Math::sqrt( result2 ), 1e-4 );
        }
    }
} // nrm2Test

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( scalTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::scal<ValueType> > scal;
    ContextPtr loc = Context::getContextPtr( scal.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "scal< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    // same values for all subcases
    ValueType values[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
    const IndexType nValues = sizeof( values ) / sizeof( ValueType );
    // check with n <= 0
    {
        HArray<ValueType> AValues( nValues, values, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            WriteAccess<ValueType> rAValues( AValues, loc );
            scal[loc->getType()]( 0, 2.0, rAValues.get(), 2 );
        }
        {
            ReadAccess<ValueType> rAValues( AValues );

            for ( int i = 0; i < 8; ++i )
            {
                BOOST_CHECK_EQUAL( rAValues[i], values[i] );
            }
        }
    }
    // check with incX <= 0
    {
        HArray<ValueType> AValues( nValues, values, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            WriteAccess<ValueType> rAValues( AValues, loc );
            scal[loc->getType()]( 3, 2.0, rAValues.get(), 0 );
        }
        {
            ReadAccess<ValueType> rAValues( AValues );

            for ( int i = 0; i < 8; ++i )
            {
                BOOST_CHECK_EQUAL( rAValues[i], values[i] );
            }
        }
    }
    // check with n > 0 and incX > 0
    {
        HArray<ValueType> AValues( nValues, values, testContext );
        IndexType incX = 3;
        {
            SCAI_CONTEXT_ACCESS( loc );
            WriteAccess<ValueType> rAValues( AValues, loc );
            scal[loc->getType()]( 3, 2.4, rAValues.get(), incX );
        }
        {
            ReadAccess<ValueType> rAValues( AValues );
            SCAI_CHECK_CLOSE( 2.4, rAValues[0], 1e-5 );
            SCAI_CHECK_CLOSE( 9.6, rAValues[3], 1e-5 );
            SCAI_CHECK_CLOSE( 16.8, rAValues[6], 1e-5 );
        }
    }
} // scalTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( sumTest, ValueType, scai_array_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::sum<ValueType> > sum;
    ContextPtr loc = Context::getContextPtr( sum.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "sum< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    int len = 5;

    // check with n == 0
    {

        int x_int[] = { 1, 2, 3, 4, 5 };
        int y_int[] = { 7, 6, 5, 4, 3 };
        int z_int[] = { 4, 3, -2, 0, -17 };

        ValueType* x = new ValueType[len];
        ValueType* y = new ValueType[len];
        ValueType* z = new ValueType[len];

        std::copy( x_int, x_int + len, x );
        std::copy( y_int, y_int + len, y );
        std::copy( z_int, z_int + len, z );

        HArray<ValueType> Ax( len, x, testContext );
        HArray<ValueType> Ay( len, y, testContext );
        HArray<ValueType> Az( len, z, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            ReadAccess<ValueType> rAy( Ay, loc );
            WriteAccess<ValueType> wAz( Az, loc );
            sum[loc->getType()]( 0, 3.0, rAx.get(), 4.0, rAy.get(), wAz.get() );
        }
        {
            ReadAccess<ValueType> rAz( Az );

            for ( int i = 0; i < len; ++i )
            {
                BOOST_CHECK_EQUAL( rAz[i], z[i] );
            }
        }
    }
    // check with n > 0 and incX > 0
    {
        ValueType x_int[] = { 1, 2, 3, 4, 5 };
        ValueType y_int[] = { 7, 6, 5, 4, 3 };

        ValueType* x = new ValueType[len];
        ValueType* y = new ValueType[len];

        std::copy( x_int, x_int + len, x );
        std::copy( y_int, y_int + len, y );

        HArray<ValueType> Ax( len, x, testContext );
        HArray<ValueType> Ay( len, y, testContext );
        HArray<ValueType> Az( testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            ReadAccess<ValueType> rAy( Ay, loc );
            WriteOnlyAccess<ValueType> wAz( Az, loc, len );
            sum[loc->getType()]( len, 3.0, rAx.get(), 4.0, rAy.get(), wAz.get() );
        }
        {
            ReadAccess<ValueType> rAz( Az );
            BOOST_CHECK_EQUAL( 31.0, rAz[0] );
            BOOST_CHECK_EQUAL( 30.0, rAz[1] );
            BOOST_CHECK_EQUAL( 29.0, rAz[2] );
            BOOST_CHECK_EQUAL( 28.0, rAz[3] );
            BOOST_CHECK_EQUAL( 27.0, rAz[4] );
        }
    }
} // sumTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, ValueType, scai_array_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::swap<ValueType> > swap;
    ContextPtr loc = Context::getContextPtr( swap.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "swap< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    int lenx = 5;
    int leny = 7;

    // check with n <= 0
    {
        ValueType x_int[] = {   1, 2, 3, 4, 5 };
        ValueType y_int[] = {   7, 6, 5, 4, 3, 2, 1 };

        ValueType* x = new ValueType[lenx];
        ValueType* y = new ValueType[leny];

        std::copy( x_int, x_int + lenx, x );
        std::copy( y_int, y_int + leny, y );

        const IndexType incX = 2;
        const IndexType incY = 3;
        HArray<ValueType> Ax( lenx, x, testContext );
        HArray<ValueType> Ay( leny, y, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            WriteAccess<ValueType> wAValues1( Ax, loc );
            WriteAccess<ValueType> wAValues2( Ay, loc );
            swap[loc->getType()]( 0, wAValues1.get(), incX, wAValues2.get(), incY );
        }
        {
            ReadAccess<ValueType> rAx( Ax );
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( 1.0, rAx[0] );
            BOOST_CHECK_EQUAL( 7.0, rAy[0] );
            BOOST_CHECK_EQUAL( 3.0, rAx[2] );
            BOOST_CHECK_EQUAL( 4.0, rAy[3] );
            BOOST_CHECK_EQUAL( 5.0, rAx[4] );
            BOOST_CHECK_EQUAL( 1.0, rAy[6] );
        }
    }
    // check with incX == 0 and inc == 0
    {
        ValueType x_int[] = {   1, 2, 3, 4, 5 };
        ValueType y_int[] = {   7, 6, 5, 4, 3, 2, 1 };

        ValueType* x = new ValueType[lenx];
        ValueType* y = new ValueType[leny];

        std::copy( x_int, x_int + lenx, x );
        std::copy( y_int, y_int + leny, y );

        const IndexType nValues = 3;
        HArray<ValueType> Ax( lenx, x, testContext );
        HArray<ValueType> Ay( leny, y, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            WriteAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            const IndexType incX = 0;
            const IndexType incY = 0;
            swap[loc->getType()]( nValues, wAx.get(), incX, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAx( Ax );
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( 1.0, rAx[0] );
            BOOST_CHECK_EQUAL( 7.0, rAy[0] );
            BOOST_CHECK_EQUAL( 3.0, rAx[2] );
            BOOST_CHECK_EQUAL( 4.0, rAy[3] );
            BOOST_CHECK_EQUAL( 5.0, rAx[4] );
            BOOST_CHECK_EQUAL( 1.0, rAy[6] );
        }
    }
    // check with n > 0, incX > 0 and incY > 0
    {
        ValueType x_int[] = {   1, 2, 3, 4, 5 };
        ValueType y_int[] = {   7, 6, 5, 4, 3, 2, 1 };

        ValueType* x = new ValueType[lenx];
        ValueType* y = new ValueType[leny];

        std::copy( x_int, x_int + lenx, x );
        std::copy( y_int, y_int + leny, y );

        const IndexType nValues = 3;
        const IndexType incX = 2;
        const IndexType incY = 3;
        HArray<ValueType> Ax( lenx, x, testContext );
        HArray<ValueType> Ay( leny, y, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            WriteAccess<ValueType> wAValues1( Ax, loc );
            WriteAccess<ValueType> wAValues2( Ay, loc );
            swap[loc->getType()]( nValues, wAValues1.get(), incX, wAValues2.get(), incY );
        }
        {
            ReadAccess<ValueType> rAValues1( Ax );
            ReadAccess<ValueType> rAValues2( Ay );
            BOOST_CHECK_EQUAL( 7.0, rAValues1[0] );
            BOOST_CHECK_EQUAL( 1.0, rAValues2[0] );
            BOOST_CHECK_EQUAL( 4.0, rAValues1[2] );
            BOOST_CHECK_EQUAL( 3.0, rAValues2[3] );
            BOOST_CHECK_EQUAL( 1.0, rAValues1[4] );
            BOOST_CHECK_EQUAL( 5.0, rAValues2[6] );
        }
    }
} // swapTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
