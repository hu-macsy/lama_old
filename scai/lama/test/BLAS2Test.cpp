/**
 * @file BLAS2Test.cpp
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
 * @brief Contains tests for the blas2 methods.
 * @author: Bea Hornef
 * @date 15.7.2013
 * @since 1.0.0
 **/

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/memory.hpp>
#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/Scalar.hpp>

#include <test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::memory;

namespace scai
{
namespace lama
{
namespace BLAS2Test
{

template<typename ValueType>
void gemvTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( gemv, loc, BLAS, BLAS2, ValueType );
    // CblasRowMajor and CblasNoTrans
    {
        ValueType matrix[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType x[] =
        { 2.0, -1.0, 4.0 };
        ValueType y[] =
        { 10.0, -20.0, 30.0 };
        const IndexType m = 2;
        const IndexType n = 3;
        const ValueType alpha = 17.0;
        const IndexType lda = 3;
        const IndexType incX = 1;
        const ValueType beta = 13.0;
        const IndexType incY = 2;
        const ValueType result[] =
        { -74.0, 33.0 };
        LAMAArray<ValueType> Am( 6, matrix );
        LAMAArray<ValueType> Ax( 3, x );
        LAMAArray<ValueType> Ay( 3, y );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAm( Am, loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            gemv( CblasRowMajor, CblasNoTrans, m, n, alpha, rAm.get(), lda, rAx.get(), incX, beta, wAy.get(), incY,
                  NULL );
        }
        {
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( result[0], rAy[0] );
            BOOST_CHECK_EQUAL( result[1], rAy[2] );
        }
    }
    // CblasColMajor and CblasNoTrans
    {
        ValueType matrix[] =
        { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };
        ValueType x[] =
        { 2.0, -1.0, 4.0 };
        ValueType y[] =
        { 10.0, -20.0, 30.0 };
        const IndexType m = 2;
        const IndexType n = 3;
        const ValueType alpha = 17.0;
        const IndexType lda = 2;
        const IndexType incX = 1;
        const ValueType beta = 13.0;
        const IndexType incY = 2;
        const ValueType result[] =
        { -74.0, 33.0 };
        LAMAArray<ValueType> Am( 6, matrix );
        LAMAArray<ValueType> Ax( 3, x );
        LAMAArray<ValueType> Ay( 3, y );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAm( Am, loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            gemv( CblasColMajor, CblasNoTrans, m, n, alpha, rAm.get(), lda, rAx.get(), incX, beta, wAy.get(), incY,
                  NULL );
        }
        {
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( result[0], rAy[0] );
            BOOST_CHECK_EQUAL( result[1], rAy[2] );
        }
    }
    // CblasRowMajor and CblasTrans
    {
        ValueType matrix[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType x[] =
        { 2.0, -1.0, 4.0 };
        ValueType y[] =
        { 10.0, -20.0, 30.0 };
        const IndexType m = 2;
        const IndexType n = 3;
        const ValueType alpha = 17.0;
        const IndexType lda = 3;
        const IndexType incX = 2;
        const ValueType beta = 13.0;
        const IndexType incY = 1;
        const ValueType result[] =
        { 436.0, 148.0, -120.0 };
        LAMAArray<ValueType> Am( 6, matrix );
        LAMAArray<ValueType> Ax( 3, x );
        LAMAArray<ValueType> Ay( 3, y );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAm( Am, loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            gemv( CblasRowMajor, CblasTrans, m, n, alpha, rAm.get(), lda, rAx.get(), incX, beta, wAy.get(), incY,
                  NULL );
        }
        {
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( result[0], rAy[0] );
            BOOST_CHECK_EQUAL( result[1], rAy[1] );
            BOOST_CHECK_EQUAL( result[2], rAy[2] );
        }
    }
    // CblasColMajor and CblasTrans
    {
        ValueType matrix[] =
        { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };
        ValueType x[] =
        { 2.0, -1.0, 4.0 };
        ValueType y[] =
        { 10.0, -20.0, 30.0 };
        const IndexType m = 2;
        const IndexType n = 3;
        const ValueType alpha = 17.0;
        const IndexType lda = 2;
        const IndexType incX = 2;
        const ValueType beta = 13.0;
        const IndexType incY = 1;
        const ValueType result[] =
        { 436.0, 148.0, -120.0 };
        LAMAArray<ValueType> Am( 6, matrix );
        LAMAArray<ValueType> Ax( 3, x );
        LAMAArray<ValueType> Ay( 3, y );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAm( Am, loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            gemv( CblasColMajor, CblasTrans, m, n, alpha, rAm.get(), lda, rAx.get(), incX, beta, wAy.get(), incY,
                  NULL );
        }
        {
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( result[0], rAy[0] );
            BOOST_CHECK_EQUAL( result[1], rAy[1] );
            BOOST_CHECK_EQUAL( result[2], rAy[2] );
        }
    }
} // gemvTest

} // namespace BLAS2Test
} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( BLAS2Test )

SCAI_LOG_DEF_LOGGER( logger, "Test.BLAS2Test" )

LAMA_AUTO_TEST_CASE_CT( gemvTest, BLAS2Test )
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
