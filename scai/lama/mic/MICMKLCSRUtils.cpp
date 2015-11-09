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
#include <scai/lama/mic/MICMKLCSRUtils.hpp>

// local library
#include <scai/lama/mic/MICUtils.hpp>
#include <scai/lama/CSRKernelTrait.hpp>

// internal scai libraries
#include <scai/tracing.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/hmemo/mic/MICSyncToken.hpp>

#include <scai/common/Assert.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Constants.hpp>

// external
#include <mkl.h>

namespace scai
{

using namespace hmemo;
using tasking::MICSyncToken;

namespace lama
{

SCAI_LOG_DEF_LOGGER( MICMKLCSRUtils::logger, "MIC.MKLCSRUtils" )

/* --------------------------------------------------------------------------- */

template<>
void MICMKLCSRUtils::normalGEMV(
    float result[],
    const float alpha,
    const float x[],
    const float beta,
    const float y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType /* nnz */,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const float csrValues[] )
{
    // SCAI_REGION( "MIC.MKLscsrmv" )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<float>, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet" )
    }

    if( y != result && beta != 0 )
    {
        MICUtils::set( result, y, numRows );
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

    mkl_scsrmv( &transa, const_cast<IndexType*>( &numRows ), const_cast<IndexType*>( &numColumns ),
                const_cast<float*>( &alpha ), matdescra, const_cast<float*>( csrValues ),
                const_cast<IndexType*>( csrJA ), const_cast<IndexType*>( csrIA ), const_cast<IndexType*>( csrIA + 1 ),
                const_cast<float*>( x ), const_cast<float*>( &beta ), result );
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

    if( y != result && beta != scai::common::constants::ZERO )
    {
        MICUtils::set( result, y, numRows );
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

void MICMKLCSRUtils::registerKernels()
{
    bool useMKL = true;

    // using MKL for CSR might be disabled explicitly by environment variable

    common::Settings::getEnvironment( useMKL, "USE_MKL" );

    if( !useMKL )
    {
        return;
    }

    // REGISTER1: overwrites previous settings

    // LAMA_INTERFACE_REGISTER_T( CSRUtils, normalGEMV, float )

    SCAI_LOG_INFO( logger, "register some CSR kernels implemented by MKL for MIC in Kernel Registry" )

    using kregistry::KernelRegistry;
    using common::context::MIC;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD ;   // add it or delete it

    // ToDo : routine causes problems

    // KernelRegistry::set<CSRKernelTrait::normalGEMV<float> >( normalGEMV, MIC, flag );
    // KernelRegistry::set<CSRKernelTrait::normalGEMV<double> >( normalGEMV, MIC, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool MICMKLCSRUtils::registerInterface()
{
    registerKernels();
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool MICMKLCSRUtils::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
