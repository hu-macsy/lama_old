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

// hpp
#include <scai/lama/cuda/CUSparseCSRUtils.hpp>

// local library
#include <scai/lama/cuda/CUSPARSEWrapper.hpp>
#include <scai/lama/UtilKernelTrait.hpp>
#include <scai/lama/CSRKernelTrait.hpp>

// internal scai libraries
#include <scai/hmemo/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Constants.hpp>

// CUDA
#include <cuda.h>
#include <cusparse_v2.h>

namespace scai
{

using tasking::CUDAStreamSyncToken;

/* --------------------------------------------------------------------------- */
/*     cusparse handle is needed, set by CUDAContext                           */
/* --------------------------------------------------------------------------- */

extern cusparseHandle_t CUDAContext_cusparseHandle;

namespace lama
{

SCAI_LOG_DEF_LOGGER( CUSparseCSRUtils::logger, "CUDA.CSRUtilsSparse" )

/* --------------------------------------------------------------------------- */
/*     Template specialization convertCSR2CSC<float>                           */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUSparseCSRUtils::convertCSR2CSC(
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
    SCAI_LOG_INFO( logger,
                    "convertCSR2CSC<" << common::getScalarType<ValueType>() << "> -> cusparseScsr2csc" << ", matrix size = "
                    << numRows << " x " << numColumns << ", nnz = " << numValues )

    typedef CUSPARSETrait::BLASIndexType BLASIndexType;

    if (common::TypeTraits<IndexType>::stype
                    != common::TypeTraits<BLASIndexType>::stype)
    {
        COMMON_THROWEXCEPTION("indextype mismatch");
    }

    CUSPARSEWrapper<ValueType>::csr2csc( CUDAContext_cusparseHandle,
                    numRows, numColumns, numValues,
                    csrValues, csrIA, csrJA,
                    cscValues, cscJA, cscIA,
                    CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "convertCSR2CSC" )
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             normalGEMV                                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUSparseCSRUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType nnz,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_LOG_INFO( logger, "normalGEMV<" << common::getScalarType<ValueType>() << ">" <<
                   " result[ " << numRows << "] = " << alpha << " * A(csr) * x + " << beta << " * y " )

    SCAI_LOG_DEBUG( logger, "x = " << x << ", y = " << y << ", result = " << result )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    cusparseMatDescr_t descrCSR;

    SCAI_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

    cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
    cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    if ( y != result && beta != 0.0f )
    {
        SCAI_CUDA_RT_CALL( cudaMemcpy( result, y, numRows * sizeof( ValueType ), cudaMemcpyDeviceToDevice ),
                           "cudaMemcpy for result = y" )
    }

    // call result = alpha * op(A) * x + beta * result of cusparse
    // Note: alpha, beta are passed as pointers

    SCAI_LOG_INFO( logger, "Start cusparseXcsrmv, stream = " << stream )

    CUSPARSEWrapper<ValueType>::csrmv( CUDAContext_cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        numRows, numColumns, nnz, &alpha, descrCSR,
                                        csrValues, csrIA, csrJA, x, &beta, result );

    if ( syncToken )
    {
        // set back stream for cusparse

        SCAI_CUSPARSE_CALL( cusparseSetStream( CUDAContext_cusparseHandle, 0 ),
                            "cusparseSetStream" )
    }
    else
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cusparseXcsrmv" )
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
    SCAI_REGION( "CUDA.CSR.matrixAddSizes" )

    SCAI_LOG_INFO(
        logger,
        "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )

    SCAI_CHECK_CUDA_ACCESS

    cusparseMatDescr_t descrCSR;

    SCAI_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

    cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
    cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

    int nnzA = 0; // aIA[ m ]
    int nnzB = 0;// bIA[ numColumns ]

    // we have not passed the values, so copy it from device to host

    cudaMemcpy( &nnzA, &aIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );
    cudaMemcpy( &nnzB, &bIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );

    int nnzC;

    SCAI_CUSPARSE_CALL(
        cusparseXcsrgeamNnz( CUDAContext_cusparseHandle,
                             numRows, numColumns,
                             descrCSR, nnzA, aIA, aJA,
                             descrCSR, nnzB, bIA, bJA,
                             descrCSR, cIA, &nnzC ),
        "cusparseXcsrgeamNnz" )

    // synchronization might be redundant due to the return value

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cusparseXcsrgeamNnz" )

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
    SCAI_REGION( "CUDA.CSR.matrixMultiplySizes" )

    SCAI_LOG_INFO(
        logger,
        "matrixMutliplySizes for " << m << " x " << n << " matrix" << ", diagonalProperty = " << diagonalProperty )

    SCAI_CHECK_CUDA_ACCESS

    cusparseMatDescr_t descrCSR;

    SCAI_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

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

    SCAI_CUSPARSE_CALL(
        cusparseXcsrgemmNnz( CUDAContext_cusparseHandle,
                             CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                             m, n, k,
                             descrCSR, nnzA, aIA, aJA,
                             descrCSR, nnzB, bIA, bJA,
                             descrCSR, cIA, &nnzC ),
        "cusparseXcsrgemmNnz" )

    // synchronization might be redundant due to the return value

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "convertCSR2CSC" )

    SCAI_LOG_DEBUG( logger, "matrixMultiplySizes, nnzC = " << nnzC )

    return nnzC;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             matrixAdd                                                              */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUSparseCSRUtils::matrixAdd(
    IndexType cJA[],
    ValueType cValues[],
    const IndexType cIA[],
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const ValueType alpha,
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const ValueType beta,
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[] )
{
    SCAI_REGION( "CUDA.CSR.matrixAdd" )

    SCAI_LOG_INFO( logger, "matrixAdd for " << numRows << "x" << numColumns << " matrix" )

    SCAI_CHECK_CUDA_ACCESS

    cusparseMatDescr_t descrCSR;

    SCAI_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

    cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
    cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

    int nnzA = 0; // aIA[ m ]
    int nnzB = 0;// bIA[ numColumns ]

    // we have not passed the values, so copy it

    cudaMemcpy( &nnzA, &aIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );
    cudaMemcpy( &nnzB, &bIA[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );

    // cIA requires const_cast, but will not be modified

    CUSPARSEWrapper<ValueType>::csrgeam( CUDAContext_cusparseHandle,
                          numRows, numColumns,
                          &alpha, descrCSR, nnzA, aValues, aIA, aJA,
                          &beta, descrCSR, nnzB, bValues, bIA, bJA,
                          descrCSR, cValues, const_cast<IndexType*>( cIA ), cJA );

    // synchronization might be redundant

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cusparseXcsrgeam" )
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             matrixMultiply                                                         */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUSparseCSRUtils::matrixMultiply(
    const IndexType cIA[],
    IndexType cJA[],
    ValueType cValues[],
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const ValueType alpha,
    bool diagonalProperty,
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[] )
{
    SCAI_REGION( "CUDA.CSR.matrixMultiply" )

    SCAI_LOG_INFO( logger, "matrixMultiply, result is " << m << "x" << n << " CSR storage" )

    SCAI_CHECK_CUDA_ACCESS

    cusparseMatDescr_t descrCSR;

    SCAI_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )

    cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
    cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

    int nnzA = 0; // aIA[ m ]
    int nnzB = 0;// bIA[ numColumns ]

    // we have not passed the number of non-zero values for A, B, so copy it

    cudaMemcpy( &nnzA, &aIA[m], sizeof( IndexType ), cudaMemcpyDeviceToHost );
    cudaMemcpy( &nnzB, &bIA[k], sizeof( IndexType ), cudaMemcpyDeviceToHost );

    if ( alpha != common::constants::ONE )
    {
        COMMON_THROWEXCEPTION( "cusparseMatrixMultiply only supports alpha = 1, but alpha = " << alpha )
    }

   CUSPARSEWrapper<ValueType>::csrgemm( CUDAContext_cusparseHandle,
                          CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                          m, n, k,
                          descrCSR, nnzA, aValues, aIA, aJA,
                          descrCSR, nnzB, bValues, bIA, bJA,
                          descrCSR, cValues, cIA, cJA );

    // synchronization might be redundant d

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "csrSparseMatmulX" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUSparseCSRUtils::registerKernels( bool deleteFlag )
{
    SCAI_LOG_INFO( logger, "set CSR routines for CUSparse in Interface" )

    bool useCUSparse = true;

    // using CUSparse for CSR might be disabled explicitly by environment variable

    common::Settings::getEnvironment( useCUSparse, "SCAI_CUDA_USE_CUSPARSE" );

    if ( !useCUSparse )
    {
        return;
    }

    // REGISTER1: overwrites previous settings

    using kregistry::KernelRegistry;
    using common::context::CUDA;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_REPLACE;   // priority over OpenMPBLAS

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

    KernelRegistry::set<CSRKernelTrait::matrixAddSizes>( matrixAddSizes, CUDA, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiplySizes>( matrixMultiplySizes, CUDA, flag );


#define LAMA_CUSPARSE_CSR_REGISTER(z, I, _)                                                                         \
    KernelRegistry::set<CSRKernelTrait::normalGEMV<ARITHMETIC_CUDA_TYPE_##I> >( normalGEMV, CUDA, flag );           \
    KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<ARITHMETIC_CUDA_TYPE_##I> >( convertCSR2CSC, CUDA, flag );   \
    KernelRegistry::set<CSRKernelTrait::matrixAdd<ARITHMETIC_CUDA_TYPE_##I> >( matrixAdd, CUDA, flag );             \
    KernelRegistry::set<CSRKernelTrait::matrixMultiply<ARITHMETIC_CUDA_TYPE_##I> >( matrixMultiply, CUDA, flag );

    // loop over all supported CUDA types

    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT, LAMA_CUSPARSE_CSR_REGISTER, _ )

#undef LAMA_CUSPARSE_CSR_REGISTER
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUSparseCSRUtils::CUSparseCSRUtils()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

CUSparseCSRUtils::~CUSparseCSRUtils()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

CUSparseCSRUtils CUSparseCSRUtils::guard;    // guard variable for registration

} /* end namespace lama */

} /* end namespace scai */
