/**
 * @file BLAS1Test.cpp
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
 * @brief Contains tests for the blas1 methods.
 * @author: Bea Hornef
 * @date 5.7.2013
 * @since 1.0.0
 **/

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/hmemo.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

#include <scai/blaskernel/test/TestMacros.hpp>
#include <scai/common/Math.hpp>

using namespace scai::hmemo;

namespace scai
{
namespace blaskernel
{
namespace BLAS1Test
{

template<typename ValueType>
void asumTest( ContextPtr loc )
{
    scai::kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::asum<ValueType> > asum;
    {
        ValueType values[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        const IndexType nValues = sizeof( values ) / sizeof( ValueType );
        const IndexType incX1 = 1;
        const IndexType incX2 = 2;
        const ValueType result1 = 21.0;
        const ValueType result2 = 9.0;
//        HArray<ValueType> AValues( nValues, values );
        HArray<ValueType> AValues( nValues );
        initArray( AValues, values, nValues );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            // n <= 0
            ValueType sum = asum[loc->getType()]( -1, rAValues.get(), incX1 );
            BOOST_CHECK_EQUAL( sum, 0.0 );
            // incX <= 0
            sum = asum[loc->getType()]( 3, rAValues.get(), -incX2 );
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

template<typename ValueType>
void axpyTest( ContextPtr loc )
{
    scai::kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::axpy<ValueType> > axpy;
    // check with n <= 0
    {
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0, 7.0, 8.0, -9.0 };
        const IndexType incX = 2;
        const IndexType incY = 3;
//        HArray<ValueType> Ax( 6, x );
//        HArray<ValueType> Ay( 9, y );
        HArray<ValueType> Ax( 6 );
        initArray( Ax, x, 6 );
        HArray<ValueType> Ay( 9 );
        initArray( Ay, y, 9 );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            axpy[loc->getType()]( -2, 5.0, wAx.get(), incX, wAy.get(), incY );
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
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0, 7.0, 8.0, -9.0 };
//        HArray<ValueType> Ax( 6, x );
//        HArray<ValueType> Ay( 9, y );
        HArray<ValueType> Ax( 6 );
        initArray( Ax, x, 6 );
        HArray<ValueType> Ay( 9 );
        initArray( Ay, y, 9 );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            axpy[loc->getType()]( 3, 5.0, wAx.get(), 0, wAy.get(), 0 );
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
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0, 7.0, 8.0, -9.0 };
        ValueType yResult[] =
        { 6.0, 2.0, -3.0, -11.0, 5.0, -6.0, 32.0, 8.0, -9.0 };
        const IndexType incX = 2;
        const IndexType incY = 3;
//        HArray<ValueType> Ax( 6, x );
//        HArray<ValueType> Ay( 9, y );
        HArray<ValueType> Ax( 6 );
        initArray( Ax, x, 6 );
        HArray<ValueType> Ay( 9 );
        initArray( Ay, y, 9 );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            axpy[loc->getType()]( 3, 5.0, wAx.get(), incX, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );

            for ( int i = 0; i < 9; ++i )
            {
                BOOST_CHECK_EQUAL( yResult[i], rAy[i] );
            }
        }
    }
} // axpyTest

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void copyTest( ContextPtr loc )
{
    scai::kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::copy<ValueType> > copy;

    // check with n <= 0
    {
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { -9.0, 8.0, -7.0, 6.0, 5.0, -4.0, 3.0, 2.0, -1.0 };
        const IndexType incX = 2;
        const IndexType incY = 3;
//        HArray<ValueType> Ax( 6, x );
//        HArray<ValueType> Ay( 9, y );
        HArray<ValueType> Ax( 6 );
        initArray( Ax, x, 6 );
        HArray<ValueType> Ay( 9 );
        initArray( Ay, y, 9 );
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
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { -9.0, 8.0, -7.0, 6.0, 5.0, -4.0, 3.0, 2.0, -1.0 };
        const IndexType incX = 2;
        const IndexType incY = 3;
//        HArray<ValueType> Ax( 6, x );
//        HArray<ValueType> Ay( 9, y );
        HArray<ValueType> Ax( 6 );
        initArray( Ax, x, 6 );
        HArray<ValueType> Ay( 9 );
        initArray( Ay, y, 9 );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            copy[loc->getType()]( 3, wAx.get(), -incX, wAy.get(), -incY );
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
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { -9.0, 8.0, -7.0, 6.0, 5.0, -4.0, 3.0, 2.0, -1.0 };
        ValueType yResult[] =
        { 1.0, 8.0, -7.0, -3.0, 5.0, -4.0, 5.0, 2.0, -1.0 };
        const IndexType incX = 2;
        const IndexType incY = 3;
//        HArray<ValueType> Ax( 6, x );
//        HArray<ValueType> Ay( 9, y );
        HArray<ValueType> Ax( 6 );
        initArray( Ax, x, 6 );
        HArray<ValueType> Ay( 9 );
        initArray( Ay, y, 9 );
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
template<typename ValueType>
void dotTest( ContextPtr loc )
{
    scai::kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::dot<ValueType> > dot;
    {
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0, 7.0, 8.0, -9.0 };
        const IndexType incX = 2;
        const IndexType incY = 3;
//        HArray<ValueType> Ax( 6, x );
//        HArray<ValueType> Ay( 9, y );
        HArray<ValueType> Ax( 6 );
        initArray( Ax, x, 6 );
        HArray<ValueType> Ay( 9 );
        initArray( Ay, y, 9 );
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

template<typename ValueType>
void iamaxTest( ContextPtr loc )
{
    scai::kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::iamax<ValueType> > iamax;

    {
        ValueType values[] =
        { 1.0, 2.0, 3.0, 6.0, 5.0, 6.0 };
        const IndexType nValues = sizeof( values ) / sizeof( ValueType );
        const IndexType incX1 = 1;
        const IndexType incX2 = 2; // { 1, 3, 5}
        const IndexType result1 = 3;
        const IndexType result2 = 2;
//        HArray<ValueType> AValues( nValues, values );
        HArray<ValueType> AValues( nValues );
        initArray( AValues, values, nValues );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            // n <= 0
            IndexType smallestIndexOfMax = iamax[loc->getType()]( 0, rAValues.get(), incX1 );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, 0 );
            // incX <= 0
            smallestIndexOfMax = iamax[loc->getType()]( nValues / incX1, rAValues.get(), -incX2 );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, 0 );
            // n > 0 and incX > 0
            smallestIndexOfMax = iamax[loc->getType()]( nValues / incX1, rAValues.get(), incX1 );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, result1 );
            smallestIndexOfMax = iamax[loc->getType()]( nValues / incX2, rAValues.get(), incX2 );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, result2 );
        }
    }
} // iamaxTest

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void nrm2Test( ContextPtr loc )
{
    scai::kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::nrm2<ValueType> > nrm2;

    {
        ValueType values[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        const IndexType nValues = sizeof( values ) / sizeof( ValueType );
        const IndexType incX1 = 1;
        const IndexType incX2 = 2;
        const ValueType result1 = 91.0;
        const ValueType result2 = 35.0;
//        HArray<ValueType> AValues( nValues, values );
        HArray<ValueType> AValues( nValues );
        initArray( AValues, values, nValues );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            // n <= 0
            ValueType euclideanNorm = nrm2[loc->getType()]( 0, rAValues.get(), incX1 );
            BOOST_CHECK_EQUAL( euclideanNorm, 0.0 );
            // incX <= 0
            euclideanNorm = nrm2[loc->getType()]( -1, rAValues.get(), 0 );
            BOOST_CHECK_EQUAL( euclideanNorm, 0.0 );
            // n > 0 and incX > 0
            euclideanNorm = nrm2[loc->getType()]( nValues / incX1, rAValues.get(), incX1 );
            SCAI_CHECK_CLOSE( euclideanNorm, scai::common::Math::sqrt( result1 ), 1e-4 );
            euclideanNorm = nrm2[loc->getType()]( nValues / incX2, rAValues.get(), incX2 );
            SCAI_CHECK_CLOSE( euclideanNorm, scai::common::Math::sqrt( result2 ), 1e-4 );
        }
    }
} // nrm2Test

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void scalTest( ContextPtr loc )
{
    scai::kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::scal<ValueType> > scal;

    // check with n <= 0
    {
        ValueType values[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
//        HArray<ValueType> AValues( 8, values );
        HArray<ValueType> AValues( 8 );
        initArray( AValues, values, 8 );
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
        ValueType values[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
//        HArray<ValueType> AValues( 8, values );
        HArray<ValueType> AValues( 8 );
        initArray( AValues, values, 8 );
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
        ValueType values[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
        const IndexType incX = 3;
//        HArray<ValueType> AValues( 8, values );
        HArray<ValueType> AValues( 8 );
        initArray( AValues, values, 8 );
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

template<typename ValueType>
void sumTest( ContextPtr loc )
{
    scai::kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::sum<ValueType> > sum;

    // check with n <= 0
    {
        ValueType x[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0 };
        ValueType y[] =
        { 7.0, 6.0, 5.0, 4.0, 3.0 };
        ValueType z[] =
        { 4.0, 3.0, -2.0, 0.0, -17.0 };
//        HArray<ValueType> Ax( 5, x );
//        HArray<ValueType> Ay( 5, y );
//        HArray<ValueType> Az( 5, z );
        HArray<ValueType> Ax( 5 );
        initArray( Ax, x, 5 );
        HArray<ValueType> Ay( 5 );
        initArray( Ay, y, 5 );
        HArray<ValueType> Az( 5 );
        initArray( Az, z, 5 );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            ReadAccess<ValueType> rAy( Ay, loc );
            WriteAccess<ValueType> wAz( Az, loc );
            sum[loc->getType()]( -1, 3.0, rAx.get(), 4.0, rAy.get(), wAz.get() );
        }
        {
            ReadAccess<ValueType> rAz( Az );

            for ( int i = 0; i < 5; ++i )
            {
                BOOST_CHECK_EQUAL( rAz[i], z[i] );
            }
        }
    }
    // check with n > 0 and incX > 0
    {
        ValueType x[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0 };
        ValueType y[] =
        { 7.0, 6.0, 5.0, 4.0, 3.0 };
//        HArray<ValueType> Ax( 5, x );
//        HArray<ValueType> Ay( 5, y );
//        HArray<ValueType> Az( 5 );
        HArray<ValueType> Ax( 5 );
        initArray( Ax, x, 5 );
        HArray<ValueType> Ay( 5 );
        initArray( Ay, y, 5 );
        HArray<ValueType> Az( 5 );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            ReadAccess<ValueType> rAy( Ay, loc );
            WriteAccess<ValueType> wAz( Az, loc );
            sum[loc->getType()]( 5, 3.0, rAx.get(), 4.0, rAy.get(), wAz.get() );
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

template<typename ValueType>
void swapTest( ContextPtr loc )
{
    scai::kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::swap<ValueType> > swap;

    // check with n <= 0
    {
        ValueType x[] =
        {   1.0, 2.0, 3.0, 4.0, 5.0};
        ValueType y[] =
        {   7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
        const IndexType incX = 2;
        const IndexType incY = 3;
//        HArray<ValueType> Ax( 5, x );
//        HArray<ValueType> Ay( 7, y );
        HArray<ValueType> Ax( 5 );
        initArray( Ax, x, 5 );
        HArray<ValueType> Ay( 7 );
        initArray( Ay, y, 7 );
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
    // check with incX <= 0 and incY <= 0
    {
        ValueType x[] =
        {   1.0, 2.0, 3.0, 4.0, 5.0};
        ValueType y[] =
        {   7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
        const IndexType nValues = 3;
//        HArray<ValueType> Ax( 5, x );
//        HArray<ValueType> Ay( 7, y );
        HArray<ValueType> Ax( 5 );
        initArray( Ax, x, 5 );
        HArray<ValueType> Ay( 7 );
        initArray( Ay, y, 7 );
        {
            SCAI_CONTEXT_ACCESS( loc );
            WriteAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            swap[loc->getType()]( nValues, wAx.get(), 0, wAy.get(), -1 );
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
        ValueType x[] =
        {   1.0, 2.0, 3.0, 4.0, 5.0};
        ValueType y[] =
        {   7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
        const IndexType nValues = 3;
        const IndexType incX = 2;
        const IndexType incY = 3;
//        HArray<ValueType> AValues1( 5, values1 );
//        HArray<ValueType> AValues2( 7, values2 );
        HArray<ValueType> Ax( 5 );
        initArray( Ax, x, 5 );
        HArray<ValueType> Ay( 7 );
        initArray( Ay, y, 7 );
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

} // namespace BLAS1Test

} /* end namespace blaskernel */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( BLAS1Test )

SCAI_LOG_DEF_LOGGER( logger, "Test.BLAS1Test" )

LAMA_AUTO_TEST_CASE_CT( asumTest, BLAS1Test, scai::blaskernel )
LAMA_AUTO_TEST_CASE_CT( axpyTest, BLAS1Test, scai::blaskernel )
LAMA_AUTO_TEST_CASE_CT( copyTest, BLAS1Test, scai::blaskernel )
LAMA_AUTO_TEST_CASE_CT( dotTest, BLAS1Test, scai::blaskernel )
LAMA_AUTO_TEST_CASE_CT( iamaxTest, BLAS1Test, scai::blaskernel )
LAMA_AUTO_TEST_CASE_CT( nrm2Test, BLAS1Test, scai::blaskernel )
LAMA_AUTO_TEST_CASE_CT( scalTest, BLAS1Test, scai::blaskernel )
LAMA_AUTO_TEST_CASE_CT( sumTest, BLAS1Test, scai::blaskernel )
LAMA_AUTO_TEST_CASE_CT( swapTest, BLAS1Test, scai::blaskernel )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
