/**
 * @file BLAS_BLAS2.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief BLAS_BLAS2.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

// hpp
#include <scai/blaskernel/external/BLAS_BLAS2.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/external/BLASWrapper.hpp>
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai
{

using namespace tasking;    // for running kernels asynchronously

namespace blaskernel
{

using common::TypeTraits;

SCAI_LOG_DEF_LOGGER( BLAS_BLAS2::logger, "BLAS.BLAS2" )

/* ---------------------------------------------------------------------------------------*/
/*    gemv                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void BLAS_BLAS2::gemv(
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
    SCAI_REGION( "BLAS.BLAS2.gemv" )

    SCAI_LOG_INFO( logger,
                   "gemv<" << TypeTraits<ValueType>::id() << ">: " << "m = " << m << ", n = " << n << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

    if( m == 0 )
    {
        return; // empty X, Y, A
    }

    // N == 0: empty A, but deal with X, Y, we can handle this here

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported yet." )
    }

    // ToDo: error handling

//    BLASWrapper<ValueType>::gemv( order, transA, static_cast<BLASTrait::BLASIndexType>( m ), static_cast<BLASTrait::BLASIndexType>( n ), alpha, a, static_cast<BLASTrait::BLASIndexType>( lda ), x, static_cast<BLASTrait::BLASIndexType>( incX ), beta, y, static_cast<BLASTrait::BLASIndexType>( incY ));

    if( order == CblasColMajor )
    {
    	BLASTrait::BLASTrans ta = '-';

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

    	BLASWrapper<ValueType>::gemv( ta, static_cast<BLASTrait::BLASIndexType>( m ), static_cast<BLASTrait::BLASIndexType>( n ), alpha, a, static_cast<BLASTrait::BLASIndexType>( lda ), x, static_cast<BLASTrait::BLASIndexType>( incX ), beta, y, static_cast<BLASTrait::BLASIndexType>( incY ));
    }
    else if( order == CblasRowMajor )
    {
        BLASTrait::BLASTrans ta = '-';

		switch( transA )
		{
			case CblasNoTrans:
				ta = 'T';
				break;
			case CblasTrans:
				ta = 'N';
				break;
			case CblasConjTrans:
				ta = 'N';
				break;
		}

		if( common::isComplex( TypeTraits<ValueType>::stype ) && transA == CblasConjTrans )
		{
			COMMON_THROWEXCEPTION( "conj matrix vector multiply on complex numbers currently not supported" )
		}

		BLASWrapper<ValueType>::gemv( ta, static_cast<BLASTrait::BLASIndexType>( n ), static_cast<BLASTrait::BLASIndexType>( m ), alpha, a, static_cast<BLASTrait::BLASIndexType>( lda ), x, static_cast<BLASTrait::BLASIndexType>( incX ), beta, y, static_cast<BLASTrait::BLASIndexType>( incY ));
    }

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void BLAS_BLAS2::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::Host;

    bool useBLAS = false;
    int level = 0;

    useBLAS = common::Settings::getEnvironment( level, "SCAI_USE_BLAS" );

    if( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS2 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS not set or 0)" )
        return;
    }
    else if( level > 1 )
    {
        // only level 2 or level 3 wrappers might be used
        SCAI_LOG_INFO( logger,
                       "BLAS2 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_INFO( logger, "register BLAS2 wrapper routines for Host at kernel registry [" << flag << "]" )

    KernelRegistry::set<BLASKernelTrait::gemv<ValueType> >( BLAS_BLAS2::gemv, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

BLAS_BLAS2::BLAS_BLAS2()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_EXT_HOST_LIST>::call(
                            kregistry::KernelRegistry::KERNEL_REPLACE );
}

BLAS_BLAS2::~BLAS_BLAS2()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_EXT_HOST_LIST>::call(
                            kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

BLAS_BLAS2 BLAS_BLAS2::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
