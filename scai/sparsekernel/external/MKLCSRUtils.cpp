/**
 * @file MKLCSRUtils.cpp
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
 * @brief Implementation of CSR utilities with MKL
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/sparsekernel/external/MKLCSRUtils.hpp>

// local library
#include <scai/sparsekernel/openmp/OpenMPUtils.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
 
#include <scai/sparsekernel/external/MKLCSRTrait.hpp>
#include <scai/sparsekernel/external/MKLCSRWrapper.hpp>

#include <scai/sparsekernel/CSRKernelTrait.hpp>

// internal scai libraries

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/exception/UnsupportedException.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/preprocessor.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

// extern
#include <mkl_spblas.h>

namespace scai
{

using tasking::TaskSyncToken;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( MKLCSRUtils::logger, "MKL.CSRUtils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MKLCSRUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType /* nnz */,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_REGION( "MKL.scsrmv" )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << common::getScalarType<ValueType>() << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    typedef MKLCSRTrait::BLASIndexType BLASIndexType;
    typedef MKLCSRTrait::BLASTrans BLASTrans;
    typedef MKLCSRTrait::BLASMatrix BLASMatrix;

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if (common::TypeTraits<IndexType>::stype
                    != common::TypeTraits<BLASIndexType>::stype)
    {
        COMMON_THROWEXCEPTION("indextype mismatch");
    }

    if ( syncToken )
    {
        SCAI_UNSUPPORTED( "asynchronous execution not supported yet" )
        // ToDo: workaround required as boost::bind supports only up to 9 arguments
    }

    if ( y != result && beta != scai::common::constants::ZERO )
    {
        OpenMPUtils::set( result, y, numRows, common::reduction::COPY );
    }

    // performs y = alpha * A * x + beta * y


    BLASTrans transa = 'n';

    // General, - triangular, Non-Unit, C for zero-indexing

    BLASMatrix matdescra;

    matdescra[0] = 'g';
    matdescra[1] = ' ';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    // const_cast needed, MKL interface does not support it
    MKLCSRWrapper<ValueType>::csrmv( transa, numRows, numColumns, alpha, matdescra, csrValues,
                csrJA, csrIA, csrIA + 1, x, beta, result );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MKLCSRUtils::convertCSR2CSC(
    IndexType cscIA[],
    IndexType cscJA[],
    ValueType cscValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    IndexType numRows,
    IndexType numColumns,
    IndexType numValues )
{
    // Intel MKL supports csr to csc only for square matrices

    if( numRows == numColumns )
    {
        SCAI_REGION( "MKL.CSRUtils.convertCSR2CSC" )

        SCAI_LOG_INFO( logger, "convertCSR2CSC of matrix " << numRows << " x " << numColumns )

        IndexType job[] =
        { 0, 0, 0, 0, 0, 1 };

        MKLCSRWrapper<ValueType>::csr2csc( job, numRows, csrValues, csrJA, csrIA, cscValues, cscJA, cscIA );
    }
    else
    {
        OpenMPCSRUtils::convertCSR2CSC( cscIA, cscJA, cscValues, csrIA, csrJA, csrValues, numRows, numColumns,
                                        numValues );
    }
}


/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MKLCSRUtils::registerKernels( bool deleteFlag )
{
    bool useMKL = true;

    // using MKL for CSR might be disabled explicitly by environment variable

    common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" );

    if( !useMKL )
    {
        SCAI_LOG_INFO( logger, "MKL routines for Host Interface are disabled" )
        return;
    }

    // REGISTER1: give these routines priority in case of overriding

    SCAI_LOG_INFO( logger, "set CSR routines for MKL in Host Interface" )

    using kregistry::KernelRegistry;
    using common::context::Host;       // context for registration

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_REPLACE ;   // higher priority

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

#define LAMA_MKL_CSR_UTILS_REGISTER(z, I, _)                                                                   \
    KernelRegistry::set<CSRKernelTrait::normalGEMV<ARITHMETIC_HOST_TYPE_##I> >( normalGEMV, Host, flag );

    BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_MKL_CSR_UTILS_REGISTER, _ )

    // MKL conversion csr to csc has worse performance than our OpenMP Implementation
    // so we do not use it here.

//     KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<float> >( convertCSR2CSC, Host, flag );
//     KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<double> >( convertCSR2CSC, Host, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

MKLCSRUtils::MKLCSRUtils()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

MKLCSRUtils::~MKLCSRUtils()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

MKLCSRUtils MKLCSRUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
