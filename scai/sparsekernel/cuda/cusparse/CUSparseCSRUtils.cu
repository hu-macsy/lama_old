/**
 * @file CUSparseCSRUtils.cu
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Implementation of some CSR routines with CUSparse library 5.0
 * @author Thomas Brandes
 * @date 11.06.2013
 */

// hpp
#include <scai/sparsekernel/cuda/cusparse/CUSparseCSRUtils.hpp>

// local library
#include <scai/sparsekernel/cuda/cusparse/CUSPARSEWrapper.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/cuda/CUDACSRUtils.hpp>
#include <scai/utilskernel/cuda/CUDASparseUtils.hpp>
#include <scai/utilskernel/cuda/CUDAUtils.hpp>

// internal scai libraries
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Constants.hpp>

// CUDA
#include <cuda.h>
#include <cusparse_v2.h>

namespace scai
{

using tasking::CUDAStreamSyncToken;

namespace sparsekernel
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
    if ( numValues == 0 )
    {
        // handle this case on its own, caused undeterministic problems

        utilskernel::CUDAUtils::setVal<IndexType>( cscIA, numColumns + 1, 0, common::BinaryOp::COPY );

        return;
    }

    SCAI_REGION( "CUSparse.CSR.convert2CSC" )

    SCAI_LOG_INFO( logger,
                   "convertCSR2CSC<" << common::getScalarType<ValueType>() << "> -> cusparseScsr2csc" << ", matrix size = "
                   << numRows << " x " << numColumns << ", nnz = " << numValues )

    typedef CUSPARSETrait::BLASIndexType BLASIndexType;

    if ( common::TypeTraits<IndexType>::stype
            != common::TypeTraits<BLASIndexType>::stype )
    {
        COMMON_THROWEXCEPTION( "indextype mismatch" );
    }

    // note: SCAI_CHECK_CUDA_ACCESS not required due to getCurrentCUDACtx
    cusparseHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuSparseHandle();
    CUSPARSEWrapper<ValueType>::csr2csc( handle,
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
    const ValueType csrValues[],
    const common::MatrixOp op )
{
    SCAI_REGION( "CUSparse.CSR.normalGEMV" )

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

    cusparseHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuSparseHandle();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
        SCAI_CUSPARSE_CALL( cusparseSetStream( handle, stream ), "cusparseSetStream" )
    }

    if ( y != result && beta != ValueType( 0 ) )
    {
        const IndexType n = common::isTranspose( op ) ? numColumns : numRows;

        SCAI_CUDA_RT_CALL( cudaMemcpy( result, y, n * sizeof( ValueType ), cudaMemcpyDeviceToDevice ),
                           "cudaMemcpy for result = y" )
    }

    // call result = alpha * op(A) * x + beta * result of cusparse
    // Note: alpha, beta are passed as pointers
    SCAI_LOG_INFO( logger, "Start cusparseXcsrmv, stream = " << stream )

    CUSPARSEWrapper<ValueType>::csrmv( handle,
                                       common::isTranspose( op ) ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE,
                                       numRows, numColumns, nnz, &alpha, descrCSR,
                                       csrValues, csrIA, csrJA, x, &beta, result );

    if ( syncToken )
    {
        // set back stream for cusparse
        SCAI_CUSPARSE_CALL( cusparseSetStream( handle, 0 ), "cusparseSetStream" )
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
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType bIA[],
    const IndexType bJA[] )
{
    SCAI_REGION( "CUDA.cuCSR.matrixAddSizes" )
    SCAI_LOG_INFO(
        logger,
        "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" )

    // verify that the matrices are sorted, otherwise the kernel crashes
  
    bool hasSortedRows = CUDACSRUtils::hasSortedRows( aIA, aJA, numRows, numColumns, 1 );
    SCAI_ASSERT( hasSortedRows, "matrix A is not sorted, required for cuSparse call" )

    hasSortedRows = CUDACSRUtils::hasSortedRows( bIA, bJA, numRows, numColumns, 1 );
    SCAI_ASSERT( hasSortedRows, "matrix B is not sorted, required for cuSparse call" )

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
    cusparseHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuSparseHandle();
    SCAI_CUSPARSE_CALL(
        cusparseXcsrgeamNnz( handle,
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
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType bIA[],
    const IndexType bJA[] )
{
    SCAI_REGION( "CUDA.CSR.matrixMultiplySizes" )
    SCAI_LOG_INFO(
        logger,
        "matrixMutliplySizes for " << m << " x " << n << " matrix" )

    bool hasSortedRows = CUDACSRUtils::hasSortedRows( aIA, aJA, m, k, 1 );
    SCAI_ASSERT( hasSortedRows, "csr data of first matrix is not sorted, required for cuSparse call of matrixMultiply" )

    hasSortedRows = CUDACSRUtils::hasSortedRows( bIA, bJA, k, n, 1 );
    SCAI_ASSERT( hasSortedRows, "csr data of second matrix is not sorted, required for cuSparse call of matrixMultiply" )

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
    cusparseHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuSparseHandle();
    SCAI_CUSPARSE_CALL(
        cusparseXcsrgemmNnz( handle,
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
    const ValueType alpha,
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const ValueType beta,
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[] )
{
    SCAI_REGION( "CUDA.cuCSR.matrixAdd" )
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
    cusparseHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuSparseHandle();
    CUSPARSEWrapper<ValueType>::csrgeam( handle,
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

    if ( alpha != common::Constants::ONE )
    {
        COMMON_THROWEXCEPTION( "cusparseMatrixMultiply only supports alpha = 1, but alpha = " << alpha )
    }

    cusparseHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuSparseHandle();
    CUSPARSEWrapper<ValueType>::csrgemm( handle,
                                         CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                         m, n, k,
                                         descrCSR, nnzA, aValues, aIA, aJA,
                                         descrCSR, nnzB, bValues, bIA, bJA,
                                         descrCSR, cValues, cIA, cJA );
    // synchronization might be redundant d
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "csrSparseMatmulX" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUSparseCSRUtils::sortRows(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues )
{
    SCAI_LOG_INFO( logger, "sort elements in each of " << numRows << " rows" )

    SCAI_REGION( "cuSparse.CSR.sortRow" )

    cusparseHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuSparseHandle();

    cusparseMatDescr_t descrCSR;
    SCAI_CUSPARSE_CALL( cusparseCreateMatDescr( &descrCSR ), "cusparseCreateMatDescr" )
    cusparseSetMatType( descrCSR, CUSPARSE_MATRIX_TYPE_GENERAL );
    cusparseSetMatIndexBase( descrCSR, CUSPARSE_INDEX_BASE_ZERO );

    size_t pBufferSizeInBytes = 0;
    void* pBuffer = NULL;
    int* perm = NULL;

    // step 1: allocate buffer
    cusparseXcsrsort_bufferSizeExt( handle, numRows, numColumns, numValues, csrIA, csrJA, &pBufferSizeInBytes );
    SCAI_CUDA_RT_CALL( cudaMalloc( &pBuffer, sizeof( char )* pBufferSizeInBytes ), "alloc sort buffer" )

    size_t sizeOfValuesInBytes = sizeof( ValueType )* numValues;

    SCAI_LOG_DEBUG( logger, "allocated pBuffer for sort, size = " << pBufferSizeInBytes << " bytes, #rows = " << numRows 
                    << ", nnz = " << numValues << ", sizeof( values ) = " << sizeOfValuesInBytes )

    // step 2: setup permutation vector P to identity
    SCAI_CUDA_RT_CALL( cudaMalloc( ( void** )&perm, sizeof( int ) * numValues ), "alloc perm" )
    cusparseCreateIdentityPermutation( handle, numValues, perm );

    // step 3: sort CSR format
    cusparseXcsrsort( handle, numRows, numColumns, numValues, descrCSR, csrIA, csrJA, perm, pBuffer );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cusparsXcsrsort" )

    // step 4: gather sorted csrVal
    ValueType* tmp = NULL;
    SCAI_CUDA_RT_CALL( cudaMalloc( ( void** ) &tmp, sizeOfValuesInBytes ), "alloc tmp for sort" )
    SCAI_CUDA_RT_CALL( cudaMemcpy( tmp, csrValues, sizeOfValuesInBytes, cudaMemcpyDeviceToDevice ), "memcopy" )

    utilskernel::CUDASparseUtils::setGather( csrValues, tmp, perm, common::BinaryOp::COPY, numValues );

    SCAI_CUDA_RT_CALL( cudaFree( tmp ), "freeTmp" )
    SCAI_CUDA_RT_CALL( cudaFree( perm ), "freePerm" )
    SCAI_CUDA_RT_CALL( cudaFree( pBuffer ), "freeBuffer" )
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUSparseCSRUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_INFO( logger, "register CUSparseCSRUtils CUSparse-routines for CUDA at kernel registry [" << flag << "]" )
    KernelRegistry::set<CSRKernelTrait::matrixAddSizes>( matrixAddSizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiplySizes>( matrixMultiplySizes, ctx, flag );
}

template<typename ValueType>
void CUSparseCSRUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_INFO( logger, "register CUSparseCSRUtils CUSparse-routines for CUDA at kernel registry [" << flag
                   << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<ValueType> >( convertCSR2CSC, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAdd<ValueType> >( matrixAdd, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiply<ValueType> >( matrixMultiply, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sortRows<ValueType> >( sortRows, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUSparseCSRUtils::CUSparseCSRUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;
    bool useCUSparse = true;
    common::Settings::getEnvironment( useCUSparse, "SCAI_CUDA_USE_CUSPARSE" );

    // replace the own CUDA kernels as cuSPARSE library might be more efficient

    if ( useCUSparse )
    {
        Registrator::registerKernels( flag );
        kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    }
}

CUSparseCSRUtils::~CUSparseCSRUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    bool useCUSparse = true;
    common::Settings::getEnvironment( useCUSparse, "SCAI_CUDA_USE_CUSPARSE" );

    if ( useCUSparse )
    {
        Registrator::registerKernels( flag );
        kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    }
}

CUSparseCSRUtils CUSparseCSRUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
