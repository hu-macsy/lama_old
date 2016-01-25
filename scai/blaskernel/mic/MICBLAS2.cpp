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
#include <scai/blaskernel/mic/MICBLAS2.hpp>

// local library
#include <scai/blaskernel/mic/MICBLASWrapper.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

// scai library
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/hmemo/mic/MICSyncToken.hpp>
#include <scai/hmemo/mic/MICContext.hpp>

#include <scai/common/TypeTraits.hpp>

// external
#include <mkl.h>

// booost
#include <boost/preprocessor.hpp>

namespace scai
{

using tasking::MICSyncToken;
using hmemo::MICContext;

namespace blaskernel
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

template<typename ValueType>
void MICBLAS2::gemv(
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
    const IndexType incY )
{
    SCAI_LOG_INFO( logger,
                   "gemv<float>: m = " << m << ", n = " << n << ", lda = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

    if( m == 0 )
    {
        return; // empty X, Y, A
    }

    // n == 0: empty A, but deal with X, Y, we can handle this here

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
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

    const ValueType* alphaPtr = &alpha;
    const ValueType* betaPtr = &beta;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( ta, m, n, alphaPtr[0:1], aPtr, lda, xPtr, incX, betaPtr[0:1], yPtr, incY )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* a = static_cast<const ValueType*>( aPtr );
        ValueType* y = static_cast<ValueType*>( yPtr );

//        sgemv( &ta, &m, &n, &alpha, a, &lda, x, &incX, &beta, y, &incY );
        MICBLASWrapper::gemv( ta, m, n, *alphaPtr, a, lda, x, incX, *betaPtr, y, incY );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MICBLAS2::registerKernels( bool deleteFlag )
{
    SCAI_LOG_INFO( logger, "register BLAS2 kernels for MIC in Kernel Registry" )

    using kregistry::KernelRegistry;
    using common::context::MIC;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD ;   // add it or delete it

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

#define LAMA_BLAS2_REGISTER(z, I, _)                                                        \
    KernelRegistry::set<BLASKernelTrait::gemv<ARITHMETIC_MIC_TYPE_##I> >( gemv, MIC, flag ); \

    BOOST_PP_REPEAT( ARITHMETIC_MIC_TYPE_CNT, LAMA_BLAS2_REGISTER, _ )

#undef LAMA_BLAS2_REGISTER
    // all other routines are not used in LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICBLAS2::RegisterGuard::RegisterGuard()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

MICBLAS2::RegisterGuard::~RegisterGuard()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

MICBLAS2::RegisterGuard MICBLAS2::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
