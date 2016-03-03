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
#include <scai/blaskernel/mic/MICBLAS3.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/mic/MICBLASWrapper.hpp>

// internal scai libraries
#include <scai/hmemo/mic/MICSyncToken.hpp>
#include <scai/hmemo/mic/MICContext.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/Container.hpp>

// external
#include <mkl_blas.h>

namespace scai
{

using tasking::MICSyncToken;
using hmemo::MICContext;

namespace blaskernel
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

template<typename ValueType>
void MICBLAS3::gemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const ValueType alpha,
    const ValueType* a,
    const IndexType lda,
    const ValueType* b,
    const IndexType ldb,
    const ValueType beta,
    ValueType* c,
    const IndexType ldc )
{
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
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

    SCAI_LOG_INFO( logger, "gemm<" << common::TypeTraits<ValueType>::id() << ">, ta = " << ta << ", tb = " << tb << ", a has shape " << m << " x " << n )

    int device = MICContext::getCurrentDevice();

    const ValueType* alphaPtr = &alpha;
    const ValueType* betaPtr = &beta;

#pragma offload target( mic : device ), in( ta, tb, m, n, k, alphaPtr[0:1], aPtr, lda, bPtr, ldb, betaPtr[0:1], cPtr, ldc )
    {
        const ValueType* a = static_cast<const ValueType*>( aPtr );
        const ValueType* b = static_cast<const ValueType*>( bPtr );
        ValueType* c = static_cast<ValueType*>( cPtr );

        MICBLASWrapper<ValueType>::gemm( ta, tb, m, n, k, *alphaPtr, a, lda, b, ldb, *betaPtr, c, ldc );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICBLAS3::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::MIC;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register BLAS3 OpenMP-routines for MIC at kernel registry [" << flag << "]" )

    KernelRegistry::set<BLASKernelTrait::gemm<ValueType> >( MICBLAS3::gemm, MIC, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICBLAS3::RegisterGuard::RegisterGuard()
{
    kregistry::mepr::Registrator<RegistratorV, ARITHMETIC_MIC_LIST>::call(
                                kregistry::KernelRegistry::KERNEL_ADD );
}

MICBLAS3::RegisterGuard::~RegisterGuard()
{
    kregistry::mepr::Registrator<RegistratorV, ARITHMETIC_MIC_LIST>::call(
                                kregistry::KernelRegistry::KERNEL_ERASE );
}

MICBLAS3::RegisterGuard MICBLAS3::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
