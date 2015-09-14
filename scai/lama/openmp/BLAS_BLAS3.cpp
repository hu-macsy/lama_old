/**
 * @file BLAS_BLAS3.cpp
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
 * @brief BLAS_BLAS3.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/openmp/BLAS_BLAS3.hpp>

// others
#include <scai/lama/BLASInterface.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>
#include <scai/common/Settings.hpp>

// macros
#include <scai/lama/macros/unused.hpp>

// tracing with SCAI_REGION
#include <scai/tracing.hpp>
#include <scai/lama/cblas.hpp>

#include <boost/preprocessor.hpp>

namespace scai
{

namespace lama
{

using scai::tasking::SyncToken;
using scai::common::getScalarType;

SCAI_LOG_DEF_LOGGER( BLAS_BLAS3::logger, "BLAS.BLAS3" )

template<typename ValueType>
static inline
void wrapperGemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const int m,
    const int n,
    const int k,
    const ValueType alpha,
    const ValueType* a,
    const int lda,
    const ValueType* b,
    const int ldb,
    const ValueType beta,
    ValueType* c,
    const int ldc );

template<>
void wrapperGemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const int m,
    const int n,
    const int k,
    const float alpha,
    const float* a,
    const int lda,
    const float* b,
    const int ldb,
    const float beta,
    float* c,
    const int ldc )
{
    cblas_sgemm( order, transA, transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );
}

template<>
void wrapperGemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const int m,
    const int n,
    const int k,
    const double alpha,
    const double* a,
    const int lda,
    const double* b,
    const int ldb,
    const double beta,
    double* c,
    const int ldc )
{
    cblas_dgemm( order, transA, transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );
}

template<>
void wrapperGemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const int m,
    const int n,
    const int k,
    const ComplexFloat alpha,
    const ComplexFloat* a,
    const int lda,
    const ComplexFloat* b,
    const int ldb,
    const ComplexFloat beta,
    ComplexFloat* c,
    const int ldc )
{
    // Attention: alpha and beta are passed by a pointer
    cblas_cgemm( order, transA, transB, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc );
}

template<>
void wrapperGemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const int m,
    const int n,
    const int k,
    const ComplexDouble alpha,
    const ComplexDouble* a,
    const int lda,
    const ComplexDouble* b,
    const int ldb,
    const ComplexDouble beta,
    ComplexDouble* c,
    const int ldc )
{
    // Attention: alpha and beta are passed by a pointer
    cblas_zgemm( order, transA, transB, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc );
}

template<>
void wrapperGemm(
    const CBLAS_ORDER,
    const CBLAS_TRANSPOSE,
    const CBLAS_TRANSPOSE,
    const int,
    const int,
    const int,
    const LongDouble,
    const LongDouble*,
    const int,
    const LongDouble*,
    const int,
    const LongDouble,
    LongDouble*,
    const int )
{
    COMMON_THROWEXCEPTION( "LongDouble not supported by BLAS, please set USE_BLAS=0" )
}

template<>
void wrapperGemm(
    const CBLAS_ORDER,
    const CBLAS_TRANSPOSE,
    const CBLAS_TRANSPOSE,
    const int,
    const int,
    const int,
    const ComplexLongDouble,
    const ComplexLongDouble*,
    const int,
    const ComplexLongDouble*,
    const int,
    const ComplexLongDouble,
    ComplexLongDouble*,
    const int )
{
    COMMON_THROWEXCEPTION( "ComplexLongDouble not supported by BLAS, please set USE_BLAS=0" )
}

template<typename ValueType>
void BLAS_BLAS3::gemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const ValueType alpha,
    const ValueType* A,
    const IndexType lda,
    const ValueType* B,
    const IndexType ldb,
    const ValueType beta,
    ValueType* C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    SCAI_REGION( "BLAS.BLAS3.gemm" )

    SCAI_LOG_INFO( logger,
                   "gemm<" << getScalarType<ValueType>() << ">: " << "m = " << m << ", n = " << n << ", k = " << k << ", lda = " << lda << ", ldb = " << ldb << ", ldc = " << ldc << ", alpha = " << alpha << ", beta = " << beta )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    wrapperGemm( order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void BLAS_BLAS3::setInterface( BLASInterface& BLAS )
{
    // using BLAS wrappers might be disabled explicitly by environment variable

    int level = 0;

    bool useBLAS = common::Settings::getEnvironment( level, "USE_BLAS" );

    if( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS3 wrapper routines for Host Interface are disabled (USE_BLAS not set or 0)" )
        return;
    }
    else if( level > 2 )
    {
        SCAI_LOG_INFO( logger,
                       "BLAS3 wrapper routines for Host Interface are disabled (USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_INFO( logger, "set BLAS3 wrapper routines for Host Context in Interface" )

    // Note: macro takes advantage of same name for routines and type definitions
    //       ( e.g. routine CUDABLAS1::sum<ValueType> is set for BLAS::BLAS1::sum variable

#define LAMA_BLAS3_REGISTER(z, I, _)                                                  \
    LAMA_INTERFACE_REGISTER1_T( BLAS, gemm, ARITHMETIC_HOST_TYPE_##I )                \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_BLAS3_REGISTER, _ )

#undef LAMA_BLAS3_REGISTER
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS3 routines                                */
/* --------------------------------------------------------------------------- */

bool BLAS_BLAS3::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( hmemo::context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool BLAS_BLAS3::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
