/**
 * @file MICBLAS3.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief MICBLAS3.cpp
 * @author Thomas Brandes
 * @date 05.07.2013
 */

// hpp
#include <scai/blaskernel/mic/MICBLAS3.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/mic/MICBLASWrapper.hpp>

// internal scai libraries
#include <scai/tasking/mic/MICSyncToken.hpp>
#include <scai/hmemo/mic/MICContext.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/TypeTraits.hpp>

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
    switch ( trans )
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

    switch ( order )
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
void MICBLAS3::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;
    KernelRegistry::set<BLASKernelTrait::gemm<ValueType> >( MICBLAS3::gemm, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICBLAS3::RegisterGuard::RegisterGuard()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

MICBLAS3::RegisterGuard::~RegisterGuard()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

MICBLAS3::RegisterGuard MICBLAS3::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
