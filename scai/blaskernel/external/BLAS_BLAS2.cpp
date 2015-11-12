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
#include <scai/common/ScalarType.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace scai
{

using namespace tasking;    // for running kernels asynchronously

namespace blaskernel
{

using common::getScalarType;

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
                   "gemv<" << getScalarType<ValueType>() << ">: " << "m = " << m << ", n = " << n << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

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

    BLASWrapper::gemv( order, transA, static_cast<BLASWrapper::IndexType>( m ), static_cast<BLASWrapper::IndexType>( n ), alpha, a, static_cast<BLASWrapper::IndexType>( lda ), x, static_cast<BLASWrapper::IndexType>( incX ), beta, y, static_cast<BLASWrapper::IndexType>( incY ));

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void BLAS_BLAS2::registerKernels( bool deleteFlag )
{
    using kregistry::KernelRegistry;
    using common::context::Host;

    // using BLAS wrappers might be disabled explicitly by environment variable

    int level = 0;

    bool useBLAS = common::Settings::getEnvironment( level, "SCAI_USE_BLAS" );

    if( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS2 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS not set or 0)" )
        return;
    }
    else if( level > 2 )
    {
        SCAI_LOG_INFO( logger,
                       "BLAS2 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_INFO( logger, "set BLAS2 wrapper routines for Host Context in Interface" )

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_REPLACE;   // priority over OpenMPBLAS

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

#define LAMA_BLAS2_REGISTER(z, I, _)                                                        \
    KernelRegistry::set<BLASKernelTrait::gemv<ARITHMETIC_HOST_TYPE_##I> >( gemv, Host, flag ); \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_BLAS2_REGISTER, _ )

#undef LAMA_BLAS2_REGISTER

    // all other routines are not used in LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

BLAS_BLAS2::BLAS_BLAS2()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

BLAS_BLAS2::~BLAS_BLAS2()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

BLAS_BLAS2 BLAS_BLAS2::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
