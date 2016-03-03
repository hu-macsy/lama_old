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

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/external/BLASWrapper.hpp>
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/Container.hpp>

namespace scai
{

namespace blaskernel
{

using common::TypeTraits;
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
                   "gemm<" << TypeTraits<ValueType>::id() << ">: " << "m = " << m << ", n = " << n << ", k = " << k 
                   << ", lda = " << lda << ", ldb = " << ldb << ", ldc = " << ldc << ", alpha = " << alpha << ", beta = " << beta )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported yet" )
    }

    BLASTrait::BLASTrans ta = '-', tb = '-';

	switch( transA )
	{
		case CblasNoTrans:
			ta = 'N';
			break;
		case CblasTrans:
			ta = 'T';
			break;
		case CblasConjTrans:
			ta = 'C';
			break;
	}

	switch( transB )
	{
		case CblasNoTrans:
			tb = 'N';
			break;
		case CblasTrans:
			tb = 'T';
			break;
		case CblasConjTrans:
			tb = 'C';
			break;
	}

    if( order == CblasColMajor )
    {
    	BLASWrapper<ValueType>::gemm( ta, tb, static_cast<BLASTrait::BLASIndexType>( m ), static_cast<BLASTrait::BLASIndexType>( n ), static_cast<BLASTrait::BLASIndexType>( k ), alpha, A, static_cast<BLASTrait::BLASIndexType>( lda ), B, static_cast<BLASTrait::BLASIndexType>( ldb ), beta, C, static_cast<BLASTrait::BLASIndexType>( ldc ));
    }
    else if( order == CblasRowMajor )
    {
    	BLASWrapper<ValueType>::gemm( ta, tb, static_cast<BLASTrait::BLASIndexType>( n ), static_cast<BLASTrait::BLASIndexType>( m ), static_cast<BLASTrait::BLASIndexType>( k ), alpha, B, static_cast<BLASTrait::BLASIndexType>( ldb ), A, static_cast<BLASTrait::BLASIndexType>( lda ), beta, C, static_cast<BLASTrait::BLASIndexType>( ldc ));
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void BLAS_BLAS3::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::Host;
    using kregistry::KernelRegistry;

    bool useBLAS = false;
    int level = 0;

    useBLAS = common::Settings::getEnvironment( level, "SCAI_USE_BLAS" );

    if( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS3 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS not set or 0)" )
        return;
    }
    else if( level > 3 )
    {
        // only level 2 or level 3 wrappers might be used
        SCAI_LOG_INFO( logger,
                       "BLAS3 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_INFO( logger, "register BLAS3 wrapper routines for Host at kernel registry [" << flag << "]" )

    KernelRegistry::set<BLASKernelTrait::gemm<ValueType> >( BLAS_BLAS3::gemm, Host, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

BLAS_BLAS3::BLAS_BLAS3()
{
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_EXT_HOST_LIST>::call(
                            kregistry::KernelRegistry::KERNEL_REPLACE );
}

BLAS_BLAS3::~BLAS_BLAS3()
{
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_EXT_HOST_LIST>::call(
                            kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

BLAS_BLAS3 BLAS_BLAS3::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
