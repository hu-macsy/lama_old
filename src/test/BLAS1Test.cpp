/**
 * @file BLAS1Test.cpp
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
 * @brief Contains tests for the blas1 methods.
 * @author: Bea Hornef
 * @date 5.7.2013
 * @since 1.0.0
 **/

// math for sqrt
#include <cmath>

// boost
#include <boost/test/unit_test.hpp>

// others
#include <lama/ContextAccess.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/LAMAArray.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/ReadAccess.hpp>
#include <lama/Scalar.hpp>
#include <lama/WriteAccess.hpp>

#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

namespace lama
{
namespace BLAS1Test
{

template<typename ValueType>
void asumTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( asum, loc, BLAS, BLAS1, ValueType );

    {
        ValueType values[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        const IndexType nValues = sizeof( values ) / sizeof(ValueType);
        const IndexType incX1 = 1;
        const IndexType incX2 = 2;
        const ValueType result1 = 21.0;
        const ValueType result2 = 9.0;

        LAMAArray<ValueType> AValues( nValues, values );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> rAValues( AValues, loc );

            // n <= 0
            ValueType sum = asum( -1, rAValues.get(), incX1, NULL );
            BOOST_CHECK_EQUAL( sum, 0.0 );

            // incX <= 0
            sum = asum( 3, rAValues.get(), -incX2, NULL );
            BOOST_CHECK_EQUAL( sum, 0.0 );

            // std::cout << "test 1 (incX = 1)" << std::endl;
            sum = asum( 6, rAValues.get(), incX1, NULL );
            BOOST_CHECK_EQUAL( sum, result1 );

            // std::cout << "test 2 (incX = 2)" << std::endl;
            sum = asum( 3, rAValues.get(), incX2, NULL );
            BOOST_CHECK_EQUAL( sum, result2 );
        }
    }
} // asumTest

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void axpyTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( axpy, loc, BLAS, BLAS1, ValueType );

    // check with n <= 0
    {
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0, 7.0, 8.0, -9.0 };

        const IndexType incX = 2;
        const IndexType incY = 3;

        LAMAArray<ValueType> Ax( 6, x );
        LAMAArray<ValueType> Ay( 9, y );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );

            axpy( -2, 5.0, wAx.get(), incX, wAy.get(), incY, NULL );
        }

        {
            HostReadAccess<ValueType> rAy( Ay );
            for( int i = 0; i < 9; ++i )
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

        LAMAArray<ValueType> Ax( 6, x );
        LAMAArray<ValueType> Ay( 9, y );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );

            axpy( 3, 5.0, wAx.get(), 0, wAy.get(), 0, NULL );
        }

        {
            HostReadAccess<ValueType> rAy( Ay );
            for( int i = 0; i < 9; ++i )
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

        LAMAArray<ValueType> Ax( 6, x );
        LAMAArray<ValueType> Ay( 9, y );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );

            axpy( 3, 5.0, wAx.get(), incX, wAy.get(), incY, NULL );
        }

        {
            HostReadAccess<ValueType> rAy( Ay );
            for( int i = 0; i < 9; ++i )
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
    LAMA_INTERFACE_FN_T( copy, loc, BLAS, BLAS1, ValueType );

    // check with n <= 0
    {
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { -9.0, 8.0, -7.0, 6.0, 5.0, -4.0, 3.0, 2.0, -1.0 };

        const IndexType incX = 2;
        const IndexType incY = 3;

        LAMAArray<ValueType> Ax( 6, x );
        LAMAArray<ValueType> Ay( 9, y );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );

            copy( 0, wAx.get(), incX, wAy.get(), incY, NULL );
        }

        {
            HostReadAccess<ValueType> rAy( Ay );
            for( int i = 0; i < 9; ++i )
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

        LAMAArray<ValueType> Ax( 6, x );
        LAMAArray<ValueType> Ay( 9, y );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );

            copy( 3, wAx.get(), -incX, wAy.get(), -incY, NULL );
        }

        {
            HostReadAccess<ValueType> rAy( Ay );
            for( int i = 0; i < 9; ++i )
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

        LAMAArray<ValueType> Ax( 6, x );
        LAMAArray<ValueType> Ay( 9, y );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );

            copy( 3, wAx.get(), incX, wAy.get(), incY, NULL );
        }

        {
            HostReadAccess<ValueType> rAy( Ay );
            for( int i = 0; i < 9; ++i )
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
    LAMA_INTERFACE_FN_T( dot, loc, BLAS, BLAS1, ValueType );

    {
        ValueType x[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType y[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0, 7.0, 8.0, -9.0 };

        const IndexType incX = 2;
        const IndexType incY = 3;

        LAMAArray<ValueType> Ax( 6, x );
        LAMAArray<ValueType> Ay( 9, y );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );

            // n <= 0
            ValueType result = dot( 0, wAx.get(), incX, wAy.get(), incY, NULL );
            BOOST_CHECK_EQUAL( result, 0.0 );

            // incX <= 0 and incY <= 0
            result = dot( 3, wAx.get(), 0, wAy.get(), 0, NULL );
            BOOST_CHECK_EQUAL( result, 0.0 );

            // n > 0 and incX > 0 and incY > 0
            result = dot( 3, wAx.get(), incX, wAy.get(), incY, NULL );
            BOOST_CHECK_EQUAL( result, 24.0 );
        }
    }
} // dotTest

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void iamaxTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( iamax, loc, BLAS, BLAS1, ValueType );

    {
        ValueType values[] =
        { 1.0, 2.0, 3.0, 6.0, 5.0, 6.0 };
        const IndexType nValues = sizeof( values ) / sizeof(ValueType);
        const IndexType incX1 = 1;
        const IndexType incX2 = 2; // { 1, 3, 5}
        const IndexType result1 = 3;
        const IndexType result2 = 2;

        LAMAArray<ValueType> AValues( nValues, values );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> rAValues( AValues, loc );

            // n <= 0
            IndexType smallestIndexOfMax = iamax( 0, rAValues.get(), incX1, NULL );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, 0 );

            // incX <= 0
            smallestIndexOfMax = iamax( nValues / incX1, rAValues.get(), -incX2, NULL );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, 0 );

            // n > 0 and incX > 0
            smallestIndexOfMax = iamax( nValues / incX1, rAValues.get(), incX1, NULL );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, result1 );

            smallestIndexOfMax = iamax( nValues / incX2, rAValues.get(), incX2, NULL );
            BOOST_CHECK_EQUAL( smallestIndexOfMax, result2 );

        }
    }
} // iamaxTest

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void nrm2Test( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( nrm2, loc, BLAS, BLAS1, ValueType );

    {
        ValueType values[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        const IndexType nValues = sizeof( values ) / sizeof(ValueType);
        const IndexType incX1 = 1;
        const IndexType incX2 = 2;
        const ValueType result1 = 91.0;
        const ValueType result2 = 35.0;

        LAMAArray<ValueType> AValues( nValues, values );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> rAValues( AValues, loc );

            // n <= 0
            ValueType euclideanNorm = nrm2( 0, rAValues.get(), incX1, NULL );
            BOOST_CHECK_EQUAL( euclideanNorm, 0.0 );

            // incX <= 0
            euclideanNorm = nrm2( -1, rAValues.get(), 0, NULL );
            BOOST_CHECK_EQUAL( euclideanNorm, 0.0 );

            // n > 0 and incX > 0
            euclideanNorm = nrm2( nValues / incX1, rAValues.get(), incX1, NULL );
            BOOST_CHECK_CLOSE( euclideanNorm, ::sqrt(result1), 1e-4 );

            euclideanNorm = nrm2( nValues / incX2, rAValues.get(), incX2, NULL );
            BOOST_CHECK_CLOSE( euclideanNorm, ::sqrt(result2), 1e-4 );
        }
    }
} // nrm2Test

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void scalTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( scal, loc, BLAS, BLAS1, ValueType );

    // check with n <= 0
    {
        ValueType values[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
        LAMAArray<ValueType> AValues( 8, values );

        {
            LAMA_CONTEXT_ACCESS( loc );

            WriteAccess<ValueType> rAValues( AValues, loc );
            scal( 0, 2.0, rAValues.get(), 2, NULL );
        }

        {
            HostReadAccess<ValueType> rAValues( AValues );
            for( int i = 0; i < 8; ++i )
            {
                BOOST_CHECK_EQUAL( rAValues[i], values[i] );
            }
        }
    }

    // check with incX <= 0
    {
        ValueType values[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
        LAMAArray<ValueType> AValues( 8, values );

        {
            LAMA_CONTEXT_ACCESS( loc );

            WriteAccess<ValueType> rAValues( AValues, loc );
            scal( 3, 2.0, rAValues.get(), 0, NULL );
        }

        {
            HostReadAccess<ValueType> rAValues( AValues );
            for( int i = 0; i < 8; ++i )
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

        LAMAArray<ValueType> AValues( 8, values );

        {
            LAMA_CONTEXT_ACCESS( loc );

            WriteAccess<ValueType> rAValues( AValues, loc );
            scal( 3, 2.4, rAValues.get(), incX, NULL );
        }

        {
            HostReadAccess<ValueType> rAValues( AValues );
            BOOST_CHECK_CLOSE( 2.4, rAValues[0], 1e-5 );
            BOOST_CHECK_CLOSE( 9.6, rAValues[3], 1e-5 );
            BOOST_CHECK_CLOSE( 16.8, rAValues[6], 1e-5 );
        }
    }
} // scalTest

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void sumTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( sum, loc, BLAS, BLAS1, ValueType );

    // check with n <= 0
    {
        ValueType x[] =
        { 1.0, 2.0, 3.0, 4.0, 5.0 };
        ValueType y[] =
        { 7.0, 6.0, 5.0, 4.0, 3.0 };
        ValueType z[] =
        { 4.0, 3.0, -2.0, 0.0, -17.0 };

        LAMAArray<ValueType> Ax( 5, x );
        LAMAArray<ValueType> Ay( 5, y );
        LAMAArray<ValueType> Az( 5, z );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> rAx( Ax, loc );
            ReadAccess<ValueType> rAy( Ay, loc );
            WriteAccess<ValueType> wAz( Az, loc );

            sum( -1, 3.0, rAx.get(), 4.0, rAy.get(), wAz.get(), NULL );
        }

        {
            HostReadAccess<ValueType> rAz( Az );

            for( int i = 0; i < 5; ++i )
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

        LAMAArray<ValueType> Ax( 5, x );
        LAMAArray<ValueType> Ay( 5, y );
        LAMAArray<ValueType> Az( 5 );

        {
            LAMA_CONTEXT_ACCESS( loc );

            ReadAccess<ValueType> rAx( Ax, loc );
            ReadAccess<ValueType> rAy( Ay, loc );
            WriteAccess<ValueType> wAz( Az, loc );

            sum( 5, 3.0, rAx.get(), 4.0, rAy.get(), wAz.get(), NULL );
        }

        {
            HostReadAccess<ValueType> rAz( Az );
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
    LAMA_INTERFACE_FN_T( swap, loc, BLAS, BLAS1, ValueType );

    // check with n <= 0
    {
        ValueType x[] =
        {   1.0, 2.0, 3.0, 4.0, 5.0};
        ValueType y[] =
        {   7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};

        const IndexType incX = 2;
        const IndexType incY = 3;

        LAMAArray<ValueType> Ax( 5, x );
        LAMAArray<ValueType> Ay( 7, y );

        {
            LAMA_CONTEXT_ACCESS( loc );

            WriteAccess<ValueType> wAValues1( Ax, loc );
            WriteAccess<ValueType> wAValues2( Ay, loc );

            swap( 0, wAValues1.get(), incX, wAValues2.get(), incY, NULL );
        }

        {
            HostReadAccess<ValueType> rAx( Ax );
            HostReadAccess<ValueType> rAy( Ay );
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

        LAMAArray<ValueType> Ax( 5, x );
        LAMAArray<ValueType> Ay( 7, y );

        {
            LAMA_CONTEXT_ACCESS( loc );

            WriteAccess<ValueType> wAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );

            swap( nValues, wAx.get(), 0, wAy.get(), -1, NULL );
        }

        {
            HostReadAccess<ValueType> rAx( Ax );
            HostReadAccess<ValueType> rAy( Ay );
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
        ValueType values1[] =
        {   1.0, 2.0, 3.0, 4.0, 5.0};
        ValueType values2[] =
        {   7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};

        const IndexType nValues = 3;
        const IndexType incX = 2;
        const IndexType incY = 3;

        LAMAArray<ValueType> AValues1( 5, values1 );
        LAMAArray<ValueType> AValues2( 7, values2 );

        {
            LAMA_CONTEXT_ACCESS( loc );

            WriteAccess<ValueType> wAValues1( AValues1, loc );
            WriteAccess<ValueType> wAValues2( AValues2, loc );

            swap( nValues, wAValues1.get(), incX, wAValues2.get(), incY, NULL );
        }

        {
            HostReadAccess<ValueType> rAValues1( AValues1 );
            HostReadAccess<ValueType> rAValues2( AValues2 );
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

template<typename ValueType>
void viamaxTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( viamax, loc, BLAS, BLAS1, ValueType );

    ValueType values[] =
    {   1.0, 2.0, 3.0, 6.0, 5.0, 6.0, 10.0, 10.0, -10.0, -10.0};
    const IndexType incX1 = 1;
    const IndexType incX2 = 2;
    const ValueType result1 = 6.0;
    const ValueType result2 = 5.0;

    LAMAArray<ValueType> AValues( 10, values );

    {
        LAMA_CONTEXT_ACCESS( loc );

        ReadAccess<ValueType> rAValues( AValues, loc );

        // check with incX <= 0
        ValueType max = viamax( 6, rAValues.get(), -incX1, NULL );
        BOOST_CHECK_EQUAL( max, 1.0 );

        // check with n <= 0
        max = viamax( -1, rAValues.get(), incX1, NULL );
        BOOST_CHECK_EQUAL( max, 1.0 );

        // check with n > 0, incX > 0
        max = viamax( 6, rAValues.get(), incX1, NULL );
        BOOST_CHECK_EQUAL( max, result1 );

        max = viamax( 3, rAValues.get(), incX2, NULL );
        BOOST_CHECK_EQUAL( max, result2 );
    }
} // viamaxTest

} // namespace BLAS1Test
} // namespace lama

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( BLAS1Test )

LAMA_LOG_DEF_LOGGER( logger, "Test.BLAS1Test" )

LAMA_AUTO_TEST_CASE_CT( asumTest, BLAS1Test )
LAMA_AUTO_TEST_CASE_CT( axpyTest, BLAS1Test )
LAMA_AUTO_TEST_CASE_CT( copyTest, BLAS1Test )
LAMA_AUTO_TEST_CASE_CT( dotTest, BLAS1Test )
LAMA_AUTO_TEST_CASE_CT( iamaxTest, BLAS1Test )
LAMA_AUTO_TEST_CASE_CT( nrm2Test, BLAS1Test )
LAMA_AUTO_TEST_CASE_CT( scalTest, BLAS1Test )
LAMA_AUTO_TEST_CASE_CT( sumTest, BLAS1Test )
LAMA_AUTO_TEST_CASE_CT( swapTest, BLAS1Test )
LAMA_AUTO_TEST_CASE_CT( viamaxTest, BLAS1Test )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()