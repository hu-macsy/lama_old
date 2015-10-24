/**
 * @file MICBLAS2.cpp
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
 * @brief MICBLAS2.cpp
 * @author Thomas Brandes
 * @date 05.07.2013
 * @since 1.1.0
 */

// hpp
#include <scai/lama/mic/MICBLAS2.hpp>

// local library
#include <scai/lama/BLASKernelTrait.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>

// external
#include <mkl.h>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( MICBLAS2::logger, "MIC.BLAS2" )

inline static char trans2C( CBLAS_TRANSPOSE trans )
{
    // Code-Style C Dehning

    switch( trans )
    {
        case CblasNoTrans:
            return 'N';

        case CblasTrans:
            return 'T';

        case CblasConjTrans:
            return 'C';

        default:
            return ' ';
    }
}

inline static char trans2CT( CBLAS_TRANSPOSE trans )
{
    // Code-Style C Dehning

    switch( trans )
    {
        case CblasNoTrans:
            return 'T';

        case CblasTrans:
            return 'N';

        case CblasConjTrans:
            return 'N';

        default:
            return ' ';
    }
}

/** gemv */

template<>
void MICBLAS2::gemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const IndexType m,
    const IndexType n,
    const float alpha,
    const float* a,
    const IndexType lda,
    const float* x,
    const IndexType incX,
    const float beta,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    SCAI_LOG_INFO( logger,
                   "gemv<float>: m = " << m << ", n = " << n << ", lda = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

    if( m == 0 )
    {
        return; // empty X, Y, A
    }

    // n == 0: empty A, but deal with X, Y, we can handle this here

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }

    char ta = ' ';

    switch( order )
    {
        case CblasColMajor:
            ta = trans2C( transA );
            break;

        case CblasRowMajor:
            ta = trans2CT( transA );
            std::swap( m, n );
            break;

        default:
            COMMON_THROWEXCEPTION( "Illegal order setting " << order )
    }

    void* yPtr = y;
    const void* xPtr = x;
    const void* aPtr = a;

#pragma offload target( mic ), in( ta, m, n, alpha, aPtr, lda, xPtr, incX, beta, yPtr, incY )
    {
        const float* x = static_cast<const float*>( xPtr );
        const float* a = static_cast<const float*>( aPtr );
        float* y = static_cast<float*>( yPtr );

        sgemv( &ta, &m, &n, &alpha, a, &lda, x, &incX, &beta, y, &incY );
    }
}

template<>
void MICBLAS2::gemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const IndexType m,
    const IndexType n,
    const double alpha,
    const double* a,
    const IndexType lda,
    const double* x,
    const IndexType incX,
    const double beta,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    SCAI_LOG_INFO( logger,
                   "gemv<double>: m = " << m << ", n = " << n << ", lda = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

    if( m == 0 )
    {
        return; // empty X, Y, A
    }

    // n == 0: empty A, but deal with X, Y, we can handle this here

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }

    char ta = ' ';

    switch( order )
    {
        case CblasColMajor:
            ta = trans2C( transA );
            break;

        case CblasRowMajor:
            ta = trans2CT( transA );
            std::swap( m, n );
            break;

        default:
            COMMON_THROWEXCEPTION( "Illegal order setting " << order )
    }

    void* yPtr = y;
    const void* xPtr = x;
    const void* aPtr = a;

#pragma offload target( mic ), in( ta, m, n, alpha, aPtr, lda, xPtr, incX, beta, yPtr, incY )
    {
        const double* x = static_cast<const double*>( xPtr );
        const double* a = static_cast<const double*>( aPtr );
        double* y = static_cast<double*>( yPtr );

        dgemv( &ta, &m, &n, &alpha, a, &lda, x, &incX, &beta, y, &incY );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MICBLAS2::setInterface( BLASKernelTrait& BLAS )
{
    SCAI_LOG_INFO( logger, "set BLAS2 routines for MIC in Interface" )

    // Note: macro takes advantage of same name for routines and type definitions
    //       ( e.g. routine CUDABLAS1::sum<ValueType> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, gemv, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, gemv, double )

    // all other routines are not used in LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS2 routines                                */
/* --------------------------------------------------------------------------- */

bool MICBLAS2::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( context::MIC );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool MICBLAS2::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
