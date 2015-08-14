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

// for dll_import
#include <scai/lama/openmp/MKLCSRUtils.hpp>
#include <scai/lama/openmp/OpenMPUtils.hpp>
#include <scai/lama/openmp/OpenMPCSRUtils.hpp>

// others
#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>

// assert
#include <scai/lama/exception/LAMAAssert.hpp>
#include <scai/lama/Settings.hpp>

// trace
#include <scai/tracing.hpp>
#include <scai/lama/macros/unused.hpp>

#include <mkl_spblas.h>

using tasking::SyncToken;

namespace lama
{

LAMA_LOG_DEF_LOGGER( MKLCSRUtils::logger, "MKL.CSRUtils" )

/* --------------------------------------------------------------------------- */

template<>
void MKLCSRUtils::normalGEMV(
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
    LAMA_REGION( "MKL.scsrmv" )

    LAMA_LOG_INFO( logger,
                   "normalGEMV<float>, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    if( y != result && beta != 0 )
    {
        OpenMPUtils::set( result, y, numRows );
    }

    // performs y = alpha * A * x + beta * y

    char transa = 'n';

    // General, - triangular, Non-Unit, C for zero-indexing

    char matdescra[6];

    matdescra[0] = 'g';
    matdescra[1] = ' ';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    // const_cast needed, MKL interface does not support it

    mkl_scsrmv( &transa, const_cast<IndexType*>( &numRows ), const_cast<IndexType*>( &numColumns ),
                const_cast<float*>( &alpha ), matdescra, const_cast<float*>( csrValues ),
                const_cast<IndexType*>( csrJA ), const_cast<IndexType*>( csrIA ), const_cast<IndexType*>( csrIA + 1 ),
                const_cast<float*>( x ), const_cast<float*>( &beta ), result );
}

/* --------------------------------------------------------------------------- */

template<>
void MKLCSRUtils::normalGEMV(
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
    LAMA_REGION( "MKL.dcsrmv" )

    LAMA_LOG_INFO( logger,
                   "normalGEMV<double>, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    if( y != result && beta != 0 )
    {
        OpenMPUtils::set( result, y, numRows );
    }

    // performs y = alpha * A * x + beta * y

    char transa = 'n';

    // General, - triangular, Non-Unit, C for zero-indexing

    char matdescra[6];

    matdescra[0] = 'g';
    matdescra[1] = ' ';
    matdescra[2] = 'n';
    matdescra[3] = 'c';

    // const_cast needed, MKL interface does not support it

    mkl_dcsrmv( &transa, const_cast<IndexType*>( &numRows ), const_cast<IndexType*>( &numColumns ),
                const_cast<double*>( &alpha ), matdescra, const_cast<double*>( csrValues ),
                const_cast<IndexType*>( csrJA ), const_cast<IndexType*>( csrIA ), const_cast<IndexType*>( csrIA + 1 ),
                const_cast<double*>( x ), const_cast<double*>( &beta ), result );
}

/* --------------------------------------------------------------------------- */

template<>
void MKLCSRUtils::convertCSR2CSC(
    IndexType cscIA[],
    IndexType cscJA[],
    float cscValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const float csrValues[],
    IndexType numRows,
    IndexType numColumns,
    IndexType numValues )
{
    // Intel MKL supports csr to csc only for square matrices

    if( numRows == numColumns )
    {
        LAMA_REGION( "MKL.CSRUtils.convertCSR2CSC" )

        LAMA_LOG_INFO( logger, "convertCSR2CSC of matrix " << numRows << " x " << numColumns )

        int job[] =
        { 0, 0, 0, 0, 0, 1 };

        int info = 0; // not used yet

        mkl_scsrcsc( job, &numRows, const_cast<float*>( csrValues ), const_cast<IndexType*>( csrJA ),
                     const_cast<IndexType*>( csrIA ), cscValues, cscJA, cscIA, &info );
    }
    else
    {
        OpenMPCSRUtils::convertCSR2CSC( cscIA, cscJA, cscValues, csrIA, csrJA, csrValues, numRows, numColumns,
                                        numValues );
    }
}

/* --------------------------------------------------------------------------- */

template<>
void MKLCSRUtils::convertCSR2CSC(
    IndexType cscIA[],
    IndexType cscJA[],
    double cscValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const double csrValues[],
    IndexType numRows,
    IndexType numColumns,
    IndexType numValues )
{
    // Intel MKL supports csr to csc only for square matrices

    if( numRows == numColumns )
    {
        LAMA_REGION( "MKL.CSRUtils.convertCSR2CSC" )

        LAMA_LOG_INFO( logger, "convertCSR2CSC of matrix " << numRows << " x " << numColumns )

        int job[] =
        { 0, 0, 0, 0, 0, 1 };

        int info = 0; // not used yet

        mkl_dcsrcsc( job, &numRows, const_cast<double*>( csrValues ), const_cast<IndexType*>( csrJA ),
                     const_cast<IndexType*>( csrIA ), cscValues, cscJA, cscIA, &info );
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

void MKLCSRUtils::setInterface( CSRUtilsInterface& CSRUtils )
{
    bool useMKL = true;

    // using MKL for CSR might be disabled explicitly by environment variable

    Settings::getEnvironment( useMKL, "LAMA_USE_MKL" );

    if( !useMKL )
    {
        LAMA_LOG_INFO( logger, "MKL routines for Host Interface are disabled" )
        return;
    }

    // REGISTER1: give these routines priority in case of overriding

    LAMA_LOG_INFO( logger, "set CSR routines for MKL in Host Interface" )

    LAMA_INTERFACE_REGISTER1_T( CSRUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER1_T( CSRUtils, normalGEMV, double )

    // MKL conversion csr to csc has worse performance than our OpenMP Implementation
    // so we do not use it here.

    // LAMA_INTERFACE_REGISTER1_T( CSRUtils, convertCSR2CSC, float )
    // LAMA_INTERFACE_REGISTER1_T( CSRUtils, convertCSR2CSC, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool MKLCSRUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( memory::context::Host );
    setInterface( interface.CSRUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool MKLCSRUtils::initialized = registerInterface();

} // namespace lama
