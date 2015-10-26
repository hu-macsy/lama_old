/**
 * @file MICBLAS3.cpp
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
 * @brief MICBLAS3.cpp
 * @author Thomas Brandes
 * @date 05.07.2013
 * @since 1.1.0
 */

// hpp
#include <scai/lama/mic/MICBLAS3.hpp>

// local library
#include <scai/lama/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

// external
#include <mkl.h>

namespace scai
{

using tasking::SyncToken;

namespace lama
{

SCAI_LOG_DEF_LOGGER( MICBLAS3::logger, "MIC.BLAS3" )

inline static char trans2C( CBLAS_TRANSPOSE trans )
{
    // Code-Style C Dehning

    switch( trans )
    {
        case CblasTrans:
            return 'T';

        case CblasConjTrans:
            return 'C';

        case CblasNoTrans:
            return 'N';

        default:
            return ' ';
    }
}

template<>
void MICBLAS3::gemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const float alpha,
    const float* a,
    const IndexType lda,
    const float* b,
    const IndexType ldb,
    const float beta,
    float* c,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char ta = trans2C( transA );
    char tb = trans2C( transB );

    const void* aPtr = a;
    const void* bPtr = b;
    void* cPtr = c;

    switch( order )
    {
        case CblasColMajor:
            break;

        case CblasRowMajor:
            std::swap( ta, tb );
            std::swap( aPtr, bPtr );
            std::swap( lda, ldb );
            std::swap( m, n );
            break;

        default:
            COMMON_THROWEXCEPTION( "Illegal order setting " << order )
    }

    SCAI_LOG_INFO( logger, "gemm, ta = " << ta << ", tb = " << tb << ", a has shape " << m << " x " << n )

#pragma offload target( mic ), in( ta, tb, m, n, k, alpha, aPtr, lda, bPtr, ldb, beta, cPtr, ldc )
    {
        const float* a = static_cast<const float*>( aPtr );
        const float* b = static_cast<const float*>( bPtr );
        float* c = static_cast<float*>( cPtr );

        sgemm( &ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc );
    }
}

template<>
void MICBLAS3::gemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const double alpha,
    const double* a,
    const IndexType lda,
    const double* b,
    const IndexType ldb,
    const double beta,
    double* c,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char ta = trans2C( transA );
    char tb = trans2C( transB );

    const void* aPtr = a;
    const void* bPtr = b;
    void* cPtr = c;

    switch( order )
    {
        case CblasColMajor:
            break;

        case CblasRowMajor:
            std::swap( ta, tb );
            std::swap( aPtr, bPtr );
            std::swap( lda, ldb );
            std::swap( m, n );
            break;

        default:
            COMMON_THROWEXCEPTION( "Illegal order setting " << order )
    }

    SCAI_LOG_INFO( logger, "gemm, ta = " << ta << ", tb = " << tb << ", a has shape " << m << " x " << n )

#pragma offload target( mic ), in( ta, tb, m, n, k, alpha, aPtr, lda, bPtr, ldb, beta, cPtr, ldc )
    {
        const double* a = static_cast<const double*>( aPtr );
        const double* b = static_cast<const double*>( bPtr );
        double* c = static_cast<double*>( cPtr );

        dgemm( &ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MICBLAS3::registerKernels()
{
    SCAI_LOG_INFO( logger, "register BLAS3 kernels for MIC in Kernel Registry" )

    using namespace scai::kregistry;

    // ctx will contain the context for which registration is done, here MIC

    common::ContextType ctx = common::context::MIC;

    KernelRegistry::set<BLASKernelTrait::gemm<float> >( gemm, ctx );
    KernelRegistry::set<BLASKernelTrait::gemm<double> >( gemm, ctx );
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS3 routines                                */
/* --------------------------------------------------------------------------- */

bool MICBLAS3::registerInterface()
{
    registerKernels();
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool MICBLAS3::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
