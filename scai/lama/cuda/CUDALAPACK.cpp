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
#include <scai/lama/cuda/CUDALAPACK.hpp>

// local library
#include <scai/lama/cuda/CUDABLAS1.hpp>

#include <scai/lama/BLASKernelTrait.hpp>

// internal scai library
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/cuda/CUDAError.hpp>

namespace scai
{

namespace lama
{

template<>
void CUDALAPACK::laswp(
    const CBLAS_ORDER order,
    const IndexType n,
    float* A_d,
    const IndexType lda,
    const IndexType k1,
    const IndexType k2,
    const IndexType* ipiv_h,
    const IndexType incx )
{
    int info = 0;
    int i = k1;

    if( order == CblasRowMajor )
    {
        int feedback = 0;

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
        BLASHelper::XERBLA_cpu( 0, info, "slaswp", "Illegal order setting." );
    }

    if( info < 0 )
    {
        //TODO: throw exception
    }

//        return info;
}

template<>
void CUDALAPACK::laswp(
    const CBLAS_ORDER order,
    const IndexType n,
    double* A_d,
    const IndexType lda,
    const IndexType k1,
    const IndexType k2,
    const IndexType* ipiv_h,
    const IndexType incx )
{
    int info = 0;
    int i = k1;

    if( order == CblasRowMajor )
    {
        int feedback = 0;

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
        BLASHelper::XERBLA_cpu( 0, info, "dlaswp", "Illegal order setting." );
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

void CUDALAPACK::registerKernels( bool deleteFlag )
{
    using kregistry::KernelRegistry;
    using common::context::CUDA;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD;

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

    KernelRegistry::set<BLASKernelTrait::laswp<float> >( laswp, CUDA, flag );
    KernelRegistry::set<BLASKernelTrait::laswp<double> >( laswp, CUDA, flag ); 

    // other routines are not used by LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDALAPACK::CUDALAPACK()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

CUDALAPACK::~CUDALAPACK()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

CUDALAPACK CUDALAPACK::guard;    // guard variable for registration

} /* end namespace lama */

} /* end namespace scai */
