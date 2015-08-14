/**
 * @file CUSparseCSRUtils.cu
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Implementation of some CSR routines with CUSparse library 5.0
 * @author Thomas Brandes
 * @date 11.06.2013
 * @since 1.0.1
 */

#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>

#include <scai/lama/cuda/utils.cu.h>
#include <scai/lama/cuda/CUDAError.hpp>
#include <scai/lama/cuda/CUSparseCSRUtils.hpp>
#include <scai/memory/cuda/CUDAStreamSyncToken.hpp>
#include <scai/lama/Settings.hpp>

#include <cuda.h>
#include <cusparse_v2.h>

#include <scai/tracing.hpp>

#include <scai/lama/ContextFactory.hpp>

namespace lama
{

    SCAI_LOG_DEF_LOGGER( CUSparseCSRUtils::logger, "CUDA.CSRUtilsSparse" )

    /* --------------------------------------------------------------------------- */
    /*     cusparse handle is needed, set by CUDAContext                           */
    /* --------------------------------------------------------------------------- */

    extern cusparseHandle_t CUDAContext_cusparseHandle;

    /* --------------------------------------------------------------------------- */
    /*     Template specialization convertCSR2CSC<float>                           */
    /* --------------------------------------------------------------------------- */

    template<>
    void CUSparseCSRUtils::convertCSR2CSC(
                    IndexType cscIA[],
                    IndexType cscJA[],
                    float cscValues[],
                    const IndexType csrIA[],
                    const IndexType csrJA[],
                    const float csrValues[],
                    int numRows,
                    int numColumns,
                    int numValues )
    {
        SCAI_LOG_INFO( logger,
                        "convertCSR2CSC<float> -> cusparseScsr2csc" << ", matrix size = "
                        << numRows << " x " << numColumns << ", nnz = " << numValues )

        LAMA_CUSPARSE_CALL(
                        cusparseScsr2csc( CUDAContext_cusparseHandle,
                                        numRows, numColumns, numValues,
                                        csrValues, csrIA, csrJA,
                                        cscValues, cscJA, cscIA,
                                        CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO ),
                        "convertCSR2SCC<float>" )

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "convertCSR2CSC" )
    }

    /* --------------------------------------------------------------------------- */
    /*     Template specialization convertCSR2CSC<double>                          */
    /* --------------------------------------------------------------------------- */

    template<>
    void CUSparseCSRUtils::convertCSR2CSC(
                    IndexType cscIA[],
                    IndexType cscJA[],
                    double cscValues[],
                    const IndexType csrIA[],
                    const IndexType csrJA[],
                    const double csrValues[],
                    int numRows,
                    int numColumns,
                    int numValues )
    {
        SCAI_LOG_INFO( logger,
                        "convertCSR2CSC<double> -> cusparseDcsr2csc" << ", matrix size = "
                        << numRows << " x " << numColumns << ", nnz = " << numValues )

        LAMA_CUSPARSE_CALL(
                        cusparseDcsr2csc( CUDAContext_cusparseHandle,
                                        numRows, numColumns, numValues,
                                        csrValues, csrIA, csrJA,
                                        cscValues, cscJA, cscIA,
                                        CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO ),
                        "convertCSR2SCC<double>" )

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "convertCSR2CSC" )
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                             normalGEMV                                                             */
    /* ------------------------------------------------------------------------------------------------------------------ */

    template<>
    void CUSparseCSRUtils::normalGEMV(
                    float result[],
                    const float alpha,
                    const float x[],
                    const float beta,
                    const float y[],
                    const IndexType numRows,
                    const IndexType numColumns,
                    const IndexType nnz,
                    const IndexType csrIA[],
                    const IndexType csrJA[],
                    const float csrValues[],
                    SyncToken* syncToken )
    {
        SCAI_LOG_INFO( logger, "normalGEMV<float>" <<
                        " result[ " << numRows << "] = " << alpha << " * A(csr) * x + " << beta << " * y " )

        SCAI_LOG_DEBUG( logger, "x = " << x << ", y = " << y << ", result = " << result )

        LAMA_CHECK_CUDA_ACCESS

        cudaStream_t stream = 0; // default stream if no syncToken is given

        cusparseMatDescr_t descrCSR;

        LAMA_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

        cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
        cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

        if ( syncToken )
        {
            CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
            LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
            stream = cudaStreamSyncToken->getCUDAStream();
            LAMA_CUSPARSE_CALL( cusparseSetStream( CUDAContext_cusparseHandle, stream ),
                            "cusparseSetStream" )
        }

        if ( y != result && beta != 0 )
        {
            LAMA_CUDA_RT_CALL( cudaMemcpy( result, y, numRows * sizeof( float ), cudaMemcpyDeviceToDevice ),
                            "cudaMemcpy for result = y" )
        }

        // call result = alpha * op(A) * x + beta * result of cusparse
        // Note: alpha, beta are passed as pointers

        SCAI_LOG_INFO( logger, "Start cusparseScsrmv, stream = " << stream )

        LAMA_CUSPARSE_CALL( cusparseScsrmv( CUDAContext_cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        numRows, numColumns, nnz, &alpha, descrCSR,
                                        csrValues, csrIA, csrJA, x, &beta, result ),
                        "cusparseScsrmv" )

        if ( syncToken )
        {
            // set back stream for cusparse

            LAMA_CUSPARSE_CALL( cusparseSetStream( CUDAContext_cusparseHandle, 0 ),
                            "cusparseSetStream" )
        }
        else
        {
            LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cusparseXcsrgeamNnz" )
        }
    }

    template<>
    void CUSparseCSRUtils::normalGEMV(
                    double result[],
                    const double alpha,
                    const double x[],
                    const double beta,
                    const double y[],
                    const IndexType numRows,
                    const IndexType numColumns,
                    const IndexType nnz,
                    const IndexType csrIA[],
                    const IndexType csrJA[],
                    const double csrValues[],
                    SyncToken* syncToken )
    {
        SCAI_LOG_INFO( logger, "normalGEMV<double>" <<
                        " result[ " << numRows << "] = " << alpha << " * A(csr) * x + " << beta << " * y " )

        SCAI_LOG_DEBUG( logger, "x = " << x << ", y = " << y << ", result = " << result )

        LAMA_CHECK_CUDA_ACCESS

        cudaStream_t stream = 0; // default stream if no syncToken is given

        cusparseMatDescr_t descrCSR;

        LAMA_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

        cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
        cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

        if ( syncToken )
        {
            CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
            LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
            stream = cudaStreamSyncToken->getCUDAStream();
            LAMA_CUSPARSE_CALL( cusparseSetStream( CUDAContext_cusparseHandle, stream ),
                            "cusparseSetStream" )
        }

        if ( y != result && beta != 0 )
        {
            LAMA_CUDA_RT_CALL( cudaMemcpy( result, y, numRows * sizeof( double ), cudaMemcpyDeviceToDevice ),
                            "cudaMemcpy for result = y" )
        }

        // call result = alpha * op(A) * x + beta * result of cusparse
        // Note: alpha, beta are passed as pointers

        SCAI_LOG_INFO( logger, "Start cusparseDcsrmv, stream = " << stream )

        LAMA_CUSPARSE_CALL( cusparseDcsrmv( CUDAContext_cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        numRows, numColumns, nnz, &alpha, descrCSR,
                                        csrValues, csrIA, csrJA, x, &beta, result ),
                        "cusparseScsrmv" )

        if ( syncToken )
        {
            // set back stream for cusparse

            LAMA_CUSPARSE_CALL( cusparseSetStream( CUDAContext_cusparseHandle, 0 ),
                            "cusparseSetStream" )
        }
        else
        {
            LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cusparseXcsrgeamNnz" )
        }
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                             matrixAddSizes                                                         */
    /* ------------------------------------------------------------------------------------------------------------------ */

    IndexType CUSparseCSRUtils::matrixAddSizes(
                    IndexType cIA[],
                    const IndexType numRows,
                    const IndexType numColumns,
                    bool diagonalProperty,
                    const IndexType aIA[],
                    const IndexType aJA[],
                    const IndexType bIA[],
                    const IndexType bJA[] )
    {
        LAMA_REGION( "CUDA.CSR.matrixAddSizes" )

        SCAI_LOG_INFO(
                        logger,
                        "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )

        LAMA_CHECK_CUDA_ACCESS

        cusparseMatDescr_t descrCSR;

        LAMA_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

        cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
        cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

        int nnzA = 0; // aIA[ m ]
        int nnzB = 0;// bIA[ numColumns ]

        // we have not passed the values, so copy it from device to host

        cudaMemcpy( &nnzA, &aIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );
        cudaMemcpy( &nnzB, &bIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );

        int nnzC;

        LAMA_CUSPARSE_CALL(
                        cusparseXcsrgeamNnz( CUDAContext_cusparseHandle,
                                        numRows, numColumns,
                                        descrCSR, nnzA, aIA, aJA,
                                        descrCSR, nnzB, bIA, bJA,
                                        descrCSR, cIA, &nnzC ),
                        "cusparseXcsrgeamNnz" )

        // synchronization might be redundant due to the return value

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cusparseXcsrgeamNnz" )

        return nnzC;
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                             matrixMultiplySizes                                                    */
    /* ------------------------------------------------------------------------------------------------------------------ */

    IndexType CUSparseCSRUtils::matrixMultiplySizes(
                    IndexType cIA[],
                    const IndexType m,
                    const IndexType n,
                    const IndexType k,
                    bool diagonalProperty,
                    const IndexType aIA[],
                    const IndexType aJA[],
                    const IndexType bIA[],
                    const IndexType bJA[] )
    {
        LAMA_REGION( "CUDA.CSR.matrixMultiplySizes" )

        SCAI_LOG_INFO(
                        logger,
                        "matrixMutliplySizes for " << m << " x " << n << " matrix" << ", diagonalProperty = " << diagonalProperty )

        LAMA_CHECK_CUDA_ACCESS

        cusparseMatDescr_t descrCSR;

        LAMA_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

        cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
        cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

        int nnzA = 0; // aIA[ m ]
        int nnzB = 0;// bIA[ numColumns ]

        // we have not passed the values, so copy it

        cudaMemcpy( &nnzA, &aIA[m], sizeof( IndexType ), cudaMemcpyDeviceToHost );
        cudaMemcpy( &nnzB, &bIA[k], sizeof( IndexType ), cudaMemcpyDeviceToHost );

        int nnzC;

        SCAI_LOG_DEBUG( logger, "multSizes, A is " << m << " x " << k << ", nnz = " << nnzA
                        << ", B is " << k << " x " << n << ", nnz = " << nnzB
                        << ", C = " << m << " x " << n )

        LAMA_CUSPARSE_CALL(
                        cusparseXcsrgemmNnz( CUDAContext_cusparseHandle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        m, n, k,
                                        descrCSR, nnzA, aIA, aJA,
                                        descrCSR, nnzB, bIA, bJA,
                                        descrCSR, cIA, &nnzC ),
                        "cusparseXcsrgemmNnz" )

        // synchronization might be redundant due to the return value

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "convertCSR2CSC" )

        SCAI_LOG_DEBUG( logger, "matrixMultiplySizes, nnzC = " << nnzC )

        return nnzC;
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                             matrixAdd                                                              */
    /* ------------------------------------------------------------------------------------------------------------------ */

    template<>
    void CUSparseCSRUtils::matrixAdd(
                    IndexType cJA[],
                    float cValues[],
                    const IndexType cIA[],
                    const IndexType numRows,
                    const IndexType numColumns,
                    bool diagonalProperty,
                    const float alpha,
                    const IndexType aIA[],
                    const IndexType aJA[],
                    const float aValues[],
                    const float beta,
                    const IndexType bIA[],
                    const IndexType bJA[],
                    const float bValues[] )
    {
        LAMA_REGION( "CUDA.CSR.matrixAdd" )

        SCAI_LOG_INFO( logger, "matrixAdd for " << numRows << "x" << numColumns << " matrix" )

        LAMA_CHECK_CUDA_ACCESS

        cusparseMatDescr_t descrCSR;

        LAMA_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

        cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
        cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

        int nnzA = 0; // aIA[ m ]
        int nnzB = 0;// bIA[ numColumns ]

        // we have not passed the values, so copy it

        cudaMemcpy( &nnzA, &aIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );
        cudaMemcpy( &nnzB, &bIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );

        // cIA requires const_cast, but will not be modified

        LAMA_CUSPARSE_CALL(
                        cusparseScsrgeam( CUDAContext_cusparseHandle,
                                        numRows, numColumns,
                                        &alpha, descrCSR, nnzA, aValues, aIA, aJA,
                                        &beta, descrCSR, nnzB, bValues, bIA, bJA,
                                        descrCSR, cValues, const_cast<IndexType*>( cIA ), cJA ),
                        "cusparseScsrgeam" )

        // synchronization might be redundant

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cusparseScsrgeam" )
    }

    template<>
    void CUSparseCSRUtils::matrixAdd(
                    IndexType cJA[],
                    double cValues[],
                    const IndexType cIA[],
                    const IndexType numRows,
                    const IndexType numColumns,
                    bool diagonalProperty,
                    const double alpha,
                    const IndexType aIA[],
                    const IndexType aJA[],
                    const double aValues[],
                    const double beta,
                    const IndexType bIA[],
                    const IndexType bJA[],
                    const double bValues[] )
    {
        LAMA_REGION( "CUDA.CSR.matrixAdd" )

        SCAI_LOG_INFO( logger, "matrixAdd for " << numRows << "x" << numColumns << " matrix" )

        LAMA_CHECK_CUDA_ACCESS

        cusparseMatDescr_t descrCSR;

        LAMA_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

        cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
        cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

        int nnzA = 0; // aIA[ numRows ]
        int nnzB = 0;// bIA[ numColumns ]

        // we have not passed the number of non-zero values for A, B, so copy it

        cudaMemcpy( &nnzA, &aIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );
        cudaMemcpy( &nnzB, &bIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );

        // cIA requires const_cast, but will not be modified

        LAMA_CUSPARSE_CALL(
                        cusparseDcsrgeam( CUDAContext_cusparseHandle,
                                        numRows, numColumns,
                                        &alpha, descrCSR, nnzA, aValues, aIA, aJA,
                                        &beta, descrCSR, nnzB, bValues, bIA, bJA,
                                        descrCSR, cValues, const_cast<IndexType*>( cIA ), cJA ),
                        "cusparseDcsrgeam" )

        // synchronization might be redundant

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cusparseDcsrgeam" )
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                             matrixMultiply                                                         */
    /* ------------------------------------------------------------------------------------------------------------------ */

    template<>
    void CUSparseCSRUtils::matrixMultiply(
                    const IndexType cIA[],
                    IndexType cJA[],
                    float cValues[],
                    const IndexType m,
                    const IndexType n,
                    const IndexType k,
                    const float alpha,
                    bool diagonalProperty,
                    const IndexType aIA[],
                    const IndexType aJA[],
                    const float aValues[],
                    const IndexType bIA[],
                    const IndexType bJA[],
                    const float bValues[] )
    {
        LAMA_REGION( "CUDA.CSR.matrixMultiply" )

        SCAI_LOG_INFO( logger, "matrixMultiply, result is " << m << "x" << n << " CSR storage" )

        LAMA_CHECK_CUDA_ACCESS

        cusparseMatDescr_t descrCSR;

        LAMA_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

        cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
        cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

        int nnzA = 0; // aIA[ m ]
        int nnzB = 0;// bIA[ numColumns ]

        // we have not passed the number of non-zero values for A, B, so copy it

        cudaMemcpy( &nnzA, &aIA[m], sizeof( IndexType ), cudaMemcpyDeviceToHost );
        cudaMemcpy( &nnzB, &bIA[k], sizeof( IndexType ), cudaMemcpyDeviceToHost );

        LAMA_ASSERT_EQUAL_ERROR( static_cast<float>( 1 ), alpha );

        LAMA_CUSPARSE_CALL(
                        cusparseScsrgemm( CUDAContext_cusparseHandle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        m, n, k,
                                        descrCSR, nnzA, aValues, aIA, aJA,
                                        descrCSR, nnzB, bValues, bIA, bJA,
                                        descrCSR, cValues, cIA, cJA ),
                        "cusparseScsrgemm" )

        // synchronization might be redundant d

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "csrSparseMatmulS" )
    }

    template<>
    void CUSparseCSRUtils::matrixMultiply(
                    const IndexType cIA[],
                    IndexType cJA[],
                    double cValues[],
                    const IndexType m,
                    const IndexType n,
                    const IndexType k,
                    const double alpha,
                    bool diagonalProperty,
                    const IndexType aIA[],
                    const IndexType aJA[],
                    const double aValues[],
                    const IndexType bIA[],
                    const IndexType bJA[],
                    const double bValues[] )
    {
        LAMA_REGION( "CUDA.CSR.matrixMultiply" )

        SCAI_LOG_INFO( logger, "matrixMultiply, result is " << m << "x" << n << " CSR storage" )

        LAMA_CHECK_CUDA_ACCESS

        cusparseMatDescr_t descrCSR;

        LAMA_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

        cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
        cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

        int nnzA = 0; // aIA[ m ]
        int nnzB = 0;// bIA[ n ]

        // we have not passed the number of non-zero values for A, B, so copy it

        cudaMemcpy( &nnzA, &aIA[m], sizeof( IndexType ), cudaMemcpyDeviceToHost );
        cudaMemcpy( &nnzB, &bIA[k], sizeof( IndexType ), cudaMemcpyDeviceToHost );

        LAMA_ASSERT_EQUAL_ERROR( static_cast<double>( 1 ), alpha );

        LAMA_CUSPARSE_CALL(
                        cusparseDcsrgemm( CUDAContext_cusparseHandle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        m, n, k,
                                        descrCSR, nnzA, aValues, aIA, aJA,
                                        descrCSR, nnzB, bValues, bIA, bJA,
                                        descrCSR, cValues, cIA, cJA ),
                        "cusparseDcsrgemm" )

        // synchronization might be redundant d

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "csrSparseMatmulD" )
    }

    /* ------------------------------------------------------------------------------------------------------------------ */

    /* --------------------------------------------------------------------------- */
    /*     Template instantiations via registration routine                        */
    /* --------------------------------------------------------------------------- */

    void CUSparseCSRUtils::setInterface( CSRUtilsInterface& CSRUtils )
    {
        SCAI_LOG_INFO( logger, "set CSR routines for CUSparse in Interface" )

        bool useCUSparse = true;

        // using CUSparse for CSR might be disabled explicitly by environment variable

        Settings::getEnvironment( useCUSparse, "LAMA_USE_CUSPARSE" );

        if ( !useCUSparse )
        {
            return;
        }

        // REGISTER1: overwrites previous settings

        LAMA_INTERFACE_REGISTER1_T( CSRUtils, normalGEMV, float )
        LAMA_INTERFACE_REGISTER1_T( CSRUtils, normalGEMV, double )

        LAMA_INTERFACE_REGISTER1_T( CSRUtils, convertCSR2CSC, float )
        LAMA_INTERFACE_REGISTER1_T( CSRUtils, convertCSR2CSC, double )

        LAMA_INTERFACE_REGISTER1( CSRUtils, matrixAddSizes )
        LAMA_INTERFACE_REGISTER1( CSRUtils, matrixMultiplySizes )

        LAMA_INTERFACE_REGISTER1_T( CSRUtils, matrixAdd, float )
        LAMA_INTERFACE_REGISTER1_T( CSRUtils, matrixAdd, double )

        LAMA_INTERFACE_REGISTER1_T( CSRUtils, matrixMultiply, float )
        LAMA_INTERFACE_REGISTER1_T( CSRUtils, matrixMultiply, double )
    }

    /* --------------------------------------------------------------------------- */
    /*    Static registration of the Utils routines                                */
    /* --------------------------------------------------------------------------- */

    bool CUSparseCSRUtils::registerInterface()
    {
        LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
        setInterface( interface.CSRUtils );
        return true;
    }

    /* --------------------------------------------------------------------------- */
    /*    Static initialiazion at program start                                    */
    /* --------------------------------------------------------------------------- */

    bool CUSparseCSRUtils::initialized = registerInterface();

} // namespace lama
