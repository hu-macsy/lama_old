/**
 * @file BLAS_BLAS2.cpp
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
 * @brief BLAS_BLAS2.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
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
#include <scai/common/mepr/Container.hpp>

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

		if( common::scalar::isComplex( TypeTraits<ValueType>::stype ) && transA == CblasConjTrans )
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
void BLAS_BLAS2::Registrator<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::Host;
    using kregistry::KernelRegistry;

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

    KernelRegistry::set<BLASKernelTrait::gemv<ValueType> >( BLAS_BLAS2::gemv, Host, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

BLAS_BLAS2::BLAS_BLAS2()
{
    typedef common::mepr::Container<Registrator, ARITHMETIC_EXT_HOST> ValueTypes;

    common::mepr::instantiate( kregistry::KernelRegistry::KERNEL_REPLACE, ValueTypes() );
}

BLAS_BLAS2::~BLAS_BLAS2()
{
    typedef common::mepr::Container<Registrator, ARITHMETIC_EXT_HOST> ValueTypes;

    common::mepr::instantiate( kregistry::KernelRegistry::KERNEL_ERASE, ValueTypes() );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

BLAS_BLAS2 BLAS_BLAS2::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
