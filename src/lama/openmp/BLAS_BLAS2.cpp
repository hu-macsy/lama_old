/**
 * @file BLAS_BLAS2.cpp
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
 * @brief BLAS_BLAS2.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/openmp/BLAS_BLAS2.hpp>
#include <lama/BLASInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/Settings.hpp>

// macros
#include <lama/macros/unused.hpp>

// tracing with LAMA_REGION
#include <lama/tracing.hpp>

#include <lama/cblas.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( BLAS_BLAS2::logger, "BLAS.BLAS2" )

/* ---------------------------------------------------------------------------------------*/
/*    gemv                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline
void wrapperGemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const int m,
    const int n,
    ValueType alpha,
    const ValueType* a,
    const int lda,
    const ValueType* x,
    const int incX,
    ValueType beta,
    ValueType* y,
    const int incY );

template<>
void wrapperGemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const int m,
    const int n,
    float alpha,
    const float* a,
    const int lda,
    const float* x,
    const int incX,
    float beta,
    float* y,
    const int incY )
{
    cblas_sgemv( order, transA, m, n, alpha, a, lda, x, incX, beta, y, incY );
}

template<>
void wrapperGemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const int m,
    const int n,
    double alpha,
    const double* a,
    const int lda,
    const double* x,
    const int incX,
    double beta,
    double* y,
    const int incY )
{
    cblas_dgemv( order, transA, m, n, alpha, a, lda, x, incX, beta, y, incY );
}

template<>
void wrapperGemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const int m,
    const int n,
    ComplexFloat alpha,
    const ComplexFloat* a,
    const int lda,
    const ComplexFloat* x,
    const int incX,
    ComplexFloat beta,
    ComplexFloat* y,
    const int incY )
{
    // Attention: alpha, beta must be passed here as a pointer
    cblas_cgemv( order, transA, m, n, &alpha, a, lda, x, incX, &beta, y, incY );
}

template<>
void wrapperGemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const int m,
    const int n,
    ComplexDouble alpha,
    const ComplexDouble* a,
    const int lda,
    const ComplexDouble* x,
    const int incX,
    ComplexDouble beta,
    ComplexDouble* y,
    const int incY )
{
    // Attention: alpha, beta must be passed here as a pointer
    cblas_zgemv( order, transA, m, n, &alpha, a, lda, x, incX, &beta, y, incY );
}

template<>
void wrapperGemv(
    const CBLAS_ORDER,
    const CBLAS_TRANSPOSE,
    const int,
    const int,
    LongDouble,
    const LongDouble*,
    const int,
    const LongDouble*,
    const int,
    LongDouble,
    LongDouble*,
    const int )
{
    LAMA_THROWEXCEPTION( "LongDouble not supported by BLAS, please set LAMA_USE_BLAS=0" )
}

template<typename ValueType>
void BLAS_BLAS2::gemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const IndexType m,
    const IndexType n,
    const ValueType alpha,
    const ValueType* a,
    const IndexType lda,
    const ValueType* x,
    const IndexType incX,
    const ValueType beta,
    ValueType* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "BLAS.BLAS2.gemv" )

    LAMA_LOG_INFO( logger,
                   "gemv<" << Scalar::getType<ValueType>() << ">: " << "m = " << m << ", n = " << n << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

    if( m == 0 )
    {
        return; // empty X, Y, A
    }

    // N == 0: empty A, but deal with X, Y, we can handle this here

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    // ToDo: error handling

    wrapperGemv( order, transA, m, n, alpha, a, lda, x, incX, beta, y, incY );

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void BLAS_BLAS2::setInterface( BLASInterface& BLAS )
{
    // using BLAS wrappers might be disabled explicitly by environment variable

    int level = 0;

    bool useBLAS = Settings::getEnvironment( level, "LAMA_USE_BLAS" );

    if( !useBLAS || ( level <= 0 ) )
    {
        LAMA_LOG_INFO( logger, "BLAS2 wrapper routines for Host Interface are disabled (LAMA_USE_BLAS not set or 0)" )
        return;
    }
    else if( level > 2 )
    {
        LAMA_LOG_INFO( logger,
                       "BLAS2 wrapper routines for Host Interface are disabled (LAMA_USE_BLAS = " << level << ")" )
        return;
    }

    LAMA_LOG_INFO( logger, "set BLAS2 wrapper routines for Host Context in Interface" )

    // REGISTER1: give these routines priority in case of overriding

#define LAMA_BLAS2_REGISTER(z, I, _)                                             \
    LAMA_INTERFACE_REGISTER1_T( BLAS, gemv, ARITHMETIC_TYPE##I )                 \

    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_BLAS2_REGISTER, _ )

#undef LAMA_BLAS2_REGISTER

    // all other routines are not used in LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS2 routines                                */
/* --------------------------------------------------------------------------- */

bool BLAS_BLAS2::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool BLAS_BLAS2::initialized = registerInterface();

} /* namespace lama */

