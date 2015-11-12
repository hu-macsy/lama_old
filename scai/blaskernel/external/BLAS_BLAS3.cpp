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
#include <scai/blaskernel/external/BLAS_BLAS3.hpp>

// others
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/external/BLASWrapper.hpp>
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/ScalarType.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace scai
{

namespace blaskernel
{

using common::getScalarType;
using tasking::TaskSyncToken;

SCAI_LOG_DEF_LOGGER( BLAS_BLAS3::logger, "BLAS.BLAS3" )

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
    const IndexType ldc )
{
    SCAI_REGION( "BLAS.BLAS3.gemm" )

    SCAI_LOG_INFO( logger,
                   "gemm<" << getScalarType<ValueType>() << ">: " << "m = " << m << ", n = " << n << ", k = " << k 
                   << ", lda = " << lda << ", ldb = " << ldb << ", ldc = " << ldc << ", alpha = " << alpha << ", beta = " << beta )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported yet" )
    }

    BLASWrapper::gemm( order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void BLAS_BLAS3::registerKernels( bool deleteFlag )
{
    using kregistry::KernelRegistry;
    using common::context::Host;

    // using BLAS wrappers might be disabled explicitly by environment variable

    int level = 0;

    bool useBLAS = common::Settings::getEnvironment( level, "SCAI_USE_BLAS" );

    if( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS3 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS not set or 0)" )
        return;
    }
    else if( level > 2 )
    {
        SCAI_LOG_INFO( logger,
                       "BLAS3 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_INFO( logger, "set BLAS3 wrapper routines for Host Context in Interface" )

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_REPLACE;   // priority over OpenMPBLAS

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

    // Note: macro takes advantage of same name for routines and type definitions
    //       ( e.g. routine CUDABLAS1::sum<ValueType> is set for BLAS::BLAS1::sum variable

#define LAMA_BLAS3_REGISTER(z, I, _)                                                           \
    KernelRegistry::set<BLASKernelTrait::gemm<ARITHMETIC_HOST_TYPE_##I> >( gemm, Host, flag );  \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_BLAS3_REGISTER, _ )

#undef LAMA_BLAS3_REGISTER
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

BLAS_BLAS3::BLAS_BLAS3()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

BLAS_BLAS3::~BLAS_BLAS3()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

BLAS_BLAS3 BLAS_BLAS3::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
