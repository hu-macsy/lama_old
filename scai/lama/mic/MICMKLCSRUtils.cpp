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
 * @since 1.0.0
 */

// for dll_import
#include <scai/lama/mic/MICMKLCSRUtils.hpp>
#include <scai/lama/mic/MICUtils.hpp>

// others
#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>

// assert
#include <scai/lama/exception/LAMAAssert.hpp>
#include <scai/lama/Settings.hpp>

// trace
#include <scai/tracing.hpp>

#include <mkl.h>

namespace scai
{

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
    const float csrValues[],
    SyncToken* syncToken )
{
    SCAI_REGION( "MIC.MKLscsrmv" )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<float>, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
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
    const double csrValues[],
    SyncToken* syncToken )
{
    SCAI_REGION( "MIC.MKLdcsrmv" )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<double>, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    if( y != result && beta != 0 )
    {
        MICUtils::set( result, y, numRows );
    }

    // performs y = alpha * A * x + beta * y

    size_t csrIAPtr = (size_t) csrIA;
    size_t csrJAPtr = (size_t) csrJA;
    size_t csrValuesPtr = (size_t) csrValues;
    size_t xPtr = (size_t) x;
    size_t resultPtr = (size_t) result;

#pragma offload target( mic ), in( numRows, numColumns, alpha, xPtr, beta, resultPtr, \
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

void MICMKLCSRUtils::setInterface( CSRUtilsInterface& CSRUtils )
{
    bool useMKL = true;

    // using MKL for CSR might be disabled explicitly by environment variable

    Settings::getEnvironment( useMKL, "USE_MKL" );

    if( !useMKL )
    {
        return;
    }

    // REGISTER1: overwrites previous settings

    // LAMA_INTERFACE_REGISTER_T( CSRUtils, normalGEMV, float )

    LAMA_INTERFACE_REGISTER_T( CSRUtils, normalGEMV, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool MICMKLCSRUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::MIC );
    setInterface( interface.CSRUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool MICMKLCSRUtils::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
