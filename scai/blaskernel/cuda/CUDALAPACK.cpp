/**
 * @file CUDALAPACK.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief CUDALAPACK.cpp
 * @author lschubert
 * @date 06.07.2012
 */

// hpp
#include <scai/blaskernel/cuda/CUDALAPACK.hpp>

// local library
#include <scai/blaskernel/cuda/CUDABLAS1.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai library
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/cuda/CUDAError.hpp>

namespace scai
{

namespace blaskernel
{

/* ---------------------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( CUDALAPACK::logger, "CUDA.LAPACK" )

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void CUDALAPACK::laswp(
    const CBLAS_ORDER order,
    const IndexType n,
    ValueType* A_d,
    const IndexType lda,
    const IndexType k1,
    const IndexType k2,
    const IndexType* ipiv_h,
    const IndexType incx )
{
    IndexType info = 0;
    IndexType i = k1;

    if( order == CblasRowMajor )
    {
        IndexType feedback = 0;

        for( i = k1; i < k2 /*&& feedback == LAMA_STATUS_SUCCESS*/; ++i )
        {
            if( ipiv_h[i * incx] == i )
            {
                continue;
            }

            CUDABLAS1::swap( n, &A_d[ipiv_h[i * incx] * lda], incx, &A_d[i * lda], incx );
            SCAI_CHECK_CUDA_ERROR
        }

        info = -1 * (IndexType) feedback;
    }
    else if( order == CblasColMajor )
    {
        info = n + lda;
    }
    else
    {
        info = 1;
        COMMON_THROWEXCEPTION( "illegal order setting " << order )
    }

    if( info < 0 )
    {
        //TODO: throw exception
    }

//        return info;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDALAPACK::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::CUDA;

    SCAI_LOG_INFO( logger, "register LAPACK routines implemented by CuBLAS in KernelRegistry [" << flag << "]" )

    KernelRegistry::set<BLASKernelTrait::laswp<ValueType> >( CUDALAPACK::laswp, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDALAPACK::CUDALAPACK()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_CUDA_LIST>::call(
        kregistry::KernelRegistry::KERNEL_ADD );
}

CUDALAPACK::~CUDALAPACK()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_CUDA_LIST>::call(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

CUDALAPACK CUDALAPACK::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
