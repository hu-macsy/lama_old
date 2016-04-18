/**
 * @file MICMKLCSRUtils.cpp
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
 * @brief Implementation of CSR utilities with MICMKL
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/sparsekernel/mic/MICMKLCSRUtils.hpp>

// local library
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/mic/MICMKLCSRWrapper.hpp>

// internal scai libraries
#include <scai/utilskernel/mic/MICUtils.hpp>
#include <scai/tracing.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tasking/mic/MICSyncToken.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Constants.hpp>

// external
#include <mkl_spblas.h>

namespace scai
{

using namespace hmemo;
using tasking::MICSyncToken;
using utilskernel::MICUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( MICMKLCSRUtils::logger, "MIC.MKLCSRUtils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICMKLCSRUtils::normalGEMV(
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
    // SCAI_REGION( "MIC.MKLscsrmv" )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << common::getScalarType<ValueType>() << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet" )
    }

    if( y != result && beta != 0 )
    {
        MICUtils::set( result, y, numRows, common::reduction::COPY );
    }

    // performs y = alpha * A * x + beta * y

    char transa = 'n';

    // General, - triangular, Non-Unit, C for zero-indexing

    char matdescra[6];

    matdescra[0] = 'g';
    matdescra[1] = ' ';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    // const_cast needed, MICMKL interface does not support it

    MICMKLCSRWrapper<ValueType>::csrmv( transa, numRows, numColumns, alpha, matdescra, csrValues, csrJA, csrIA, csrIA + 1, x, beta, result );

//    mkl_scsrmv( &transa, const_cast<IndexType*>( &numRows ), const_cast<IndexType*>( &numColumns ),
//                const_cast<float*>( &alpha ), matdescra, const_cast<float*>( csrValues ),
//                const_cast<IndexType*>( csrJA ), const_cast<IndexType*>( csrIA ), const_cast<IndexType*>( csrIA + 1 ),
//                const_cast<float*>( x ), const_cast<float*>( &beta ), result );
}

/* --------------------------------------------------------------------------- */

template<>
void MICMKLCSRUtils::normalGEMV(
    double result[],
    const double alpha,
    const double x[],
    const double beta,
    const double y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType /* nnz */,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const double csrValues[] )
{
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet" )
    }

    // SCAI_REGION( "MIC.MKLdcsrmv" )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<double>, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    if( y != result && beta != common::constants::ZERO )
    {
        MICUtils::set( result, y, numRows, common::reduction::COPY );
    }

    // performs y = alpha * A * x + beta * y

    size_t csrIAPtr = (size_t) csrIA;
    size_t csrJAPtr = (size_t) csrJA;
    size_t csrValuesPtr = (size_t) csrValues;
    size_t xPtr = (size_t) x;
    size_t resultPtr = (size_t) result;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( numRows, numColumns, alpha, xPtr, beta, resultPtr, \
                                       csrValuesPtr, csrJAPtr, csrIAPtr )
    {
        char transa = 'n';

        // General, - triangular, Non-Unit, C for zero-indexing

        char matdescra[6];

        matdescra[0] = 'g';
        matdescra[1] = ' ';
        matdescra[2] = 'n';
        matdescra[3] = 'c';

        IndexType* csrIA = (IndexType*) csrIAPtr;
        IndexType* csrJA = (IndexType*) csrJAPtr;
        double* csrValues = (double*) csrValuesPtr;
        double* x = (double*) xPtr;
        double* result = (double*) resultPtr;

        mkl_dcsrmv( &transa, const_cast<IndexType*>( &numRows ), const_cast<IndexType*>( &numColumns ),
                    const_cast<double*>( &alpha ), matdescra, csrValues, csrJA, csrIA, csrIA + 1, x,
                    const_cast<double*>( &beta ), result );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICMKLCSRUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    const common::context::ContextType ctx = common::context::MIC;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register CSRUtils MKL-routines for MIC at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICMKLCSRUtils::MICMKLCSRUtils()
{
    bool useMKL = true;

   // using MKL for CSR might be disabled explicitly by environment variable

   common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" );
   int level = 0;

   if( !useMKL || ( level <= 0 ) )
   {
       return;
   }

   const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;

   kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
}

MICMKLCSRUtils::~MICMKLCSRUtils()
{
    bool useMKL = true;

   // using MKL for CSR might be disabled explicitly by environment variable

   common::Settings::getEnvironment( useMKL, "SCAI_USE_MKL" );
   int level = 0;

   if( !useMKL || ( level <= 0 ) )
   {
       return;
   }

   const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;

   kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
}

MICMKLCSRUtils MICMKLCSRUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
