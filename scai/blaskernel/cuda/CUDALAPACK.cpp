/**
 * @file CUDALAPACK.cpp
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
 * @brief CUDALAPACK.cpp
 * @author lschubert
 * @date 06.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/blaskernel/cuda/CUDALAPACK.hpp>

// local library
#include <scai/blaskernel/cuda/CUDABLAS1.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai library
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/mepr/Container.hpp>

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
    using common::context::Host;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register LAPACK routines implemented by CuBLAS in KernelRegistry [" << flag << "]" )

    KernelRegistry::set<BLASKernelTrait::laswp<ValueType> >( CUDALAPACK::laswp, Host, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDALAPACK::CUDALAPACK()
{
    typedef common::mepr::ContainerV<RegistratorV, ARITHMETIC_CUDA> ValueTypes;

    kregistry::instantiate( kregistry::KernelRegistry::KERNEL_ADD, ValueTypes() );
}

CUDALAPACK::~CUDALAPACK()
{
    typedef common::mepr::ContainerV<RegistratorV, ARITHMETIC_CUDA> ValueTypes;

    kregistry::instantiate( kregistry::KernelRegistry::KERNEL_ERASE, ValueTypes() );
}

CUDALAPACK CUDALAPACK::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
