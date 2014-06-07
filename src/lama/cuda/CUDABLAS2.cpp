/**
 * @file CUDABLAS2.cpp
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
 * @brief CUDABLAS2.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/cuda/CUDABLAS2.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>

// macros
#include <lama/macros/unused.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDABLAS2::logger, "CUDA.BLAS2" )

/* ---------------------------------------------------------------------------------------*/

/** cublasCast converts pointers to LAMA complex numbers to 
 *  cuBlas pointers for complex numbers. This is safe as both
 *  are internally represented in the same way.
 */

static inline cuFloatComplex* cublasCast( ComplexFloat* x )
{
    return reinterpret_cast<cuFloatComplex*>( x );
}

static inline cuDoubleComplex* cublasCast( ComplexDouble* x )
{
    return reinterpret_cast<cuDoubleComplex*>( x );
}

static inline const cuFloatComplex* cublasCast( const ComplexFloat* x )
{
    return reinterpret_cast<const cuFloatComplex*>( x );
}

static inline const cuDoubleComplex* cublasCast( const ComplexDouble* x )
{
    return reinterpret_cast<const cuDoubleComplex*>( x );
}

/* ---------------------------------------------------------------------------------------*/
/*    gemv                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename T>
static inline
void wrapperGemv( char trans_char, IndexType m, IndexType n,
                  T alpha, const T* A, IndexType lda, const T* x, IndexType incX, 
                  T beta, T* y, IndexType incY );

template<>
void wrapperGemv( char trans, IndexType m, IndexType n,
                  float alpha, const float* A, IndexType lda, const float* x, IndexType incX, 
                  float beta, float* y, IndexType incY )
{
    cublasSgemv( trans, m, n, alpha, A, lda, x, incX, beta, y, incY );
}

template<>
void wrapperGemv( char trans, IndexType m, IndexType n,
                  double alpha, const double* A, IndexType lda, const double* x, IndexType incX, 
                  double beta, double* y, IndexType incY )
{
    cublasDgemv( trans, m, n, alpha, A, lda, x, incX, beta, y, incY );
}

template<>
void wrapperGemv( char trans, IndexType m, IndexType n,
                  ComplexFloat alpha, const ComplexFloat* A, IndexType lda, const ComplexFloat* x, IndexType incX, 
                  ComplexFloat beta, ComplexFloat* y, IndexType incY )
{
    cublasCgemv( trans, m, n, *cublasCast( &alpha ), cublasCast( A ), lda, 
                 cublasCast( x ), incX, *cublasCast( &beta ), cublasCast( y ), incY );
}

template<>
void wrapperGemv( char trans, IndexType m, IndexType n, ComplexDouble alpha, 
                  const ComplexDouble* A, IndexType lda, 
                  const ComplexDouble* x, IndexType incX, 
                  ComplexDouble beta, ComplexDouble* y, IndexType incY )
{
    cublasZgemv( trans, m, n, *cublasCast( &alpha ), cublasCast( A ), lda, 
                 cublasCast( x ), incX, *cublasCast( &beta ), cublasCast( y ), incY );
}

/** gemv */

template<typename T>
void CUDABLAS2::gemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE trans,
    const IndexType m,
    const IndexType n,
    const T alpha,
    const T* const A,
    const IndexType lda,
    const T* const x,
    const IndexType incx,
    const T beta,
    T* const y,
    const IndexType incy,
    SyncToken* syncToken )
{
    IndexType order_m = m;
    IndexType order_n = n;
    char trans_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }

        order_m = n;
        order_n = m;
    }
    else
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    LAMA_LOG_INFO( logger, "gemv<" << Scalar::getType<T>() << "> with cuBLAS: m = " << order_m << " x " << order_n )

    wrapperGemv( trans_char, order_m, order_n, alpha, A, lda, x, incx, beta, y, incy );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDABLAS2::setInterface( BLASInterface& BLAS )
{
    LAMA_LOG_INFO( logger, "set BLAS2 routines for CUDA in Interface" )

#define LAMA_BLAS2_REGISTER(z, I, _)                                            \
    LAMA_INTERFACE_REGISTER_T( BLAS, gemv, ARITHMETIC_TYPE##I )                 \
     
    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_BLAS2_REGISTER, _ )

#undef LAMA_BLAS2_REGISTER

    // other routines are not used by LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDABLAS2::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDABLAS2::initialized = registerInterface();

} /* namespace lama */
