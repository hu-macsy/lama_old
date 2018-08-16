/**
 * @file CUSPARSEWrapper.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief ToDo: Missing description in ./sparsekernel/cuda/CUSPARSEWrapper.hpp
 * @author eschricker
 * @date 24.08.2015
 */

#pragma once

// internal scai libraries
#include <scai/sparsekernel/cuda/cusparse/CUSPARSETrait.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/TypeTraits.hpp>

// CUDA
#ifdef SCAI_COMPLEX_SUPPORTED
#include <cuComplex.h>
#endif

#include <cublas_v2.h>

namespace scai
{

namespace sparsekernel
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT CUSPARSEWrapper;

#define CUSPARSEWRAPPER_DEF( ValueType, CUSPARSEValueType, prefix )                                                             \
    template<>                                                                                                                  \
    class COMMON_DLL_IMPORTEXPORT CUSPARSEWrapper<ValueType>                                                                    \
    {                                                                                                                           \
    public:                                                                                                                     \
        typedef CUSPARSETrait::BLASIndexType BLASIndexType;                                                                     \
        typedef CUSPARSETrait::BLASTrans BLASTrans;                                                                             \
        typedef CUSPARSETrait::BLASHandle BLASHandle;                                                                           \
        typedef CUSPARSETrait::BLASMatrix BLASMatrix;                                                                           \
        typedef CUSPARSETrait::BLASOperationType BLASOperationType;                                                             \
        typedef CUSPARSETrait::BLASIndexBase BLASIndexBase;                                                                     \
        \
        static void csr2csc(                                                                                                    \
                BLASHandle handle,                                                                                                  \
                BLASIndexType m,                                                                                                    \
                BLASIndexType n,                                                                                                    \
                BLASIndexType nnz,                                                                                                  \
                const ValueType *csrVal,                                                                                            \
                const BLASIndexType *csrRowPtr,                                                                                     \
                const BLASIndexType *csrColInd,                                                                                     \
                ValueType *cscVal,                                                                                                  \
                BLASIndexType *cscRowInd,                                                                                           \
                BLASIndexType *cscColPtr,                                                                                           \
                cusparseAction_t copyValues,                                                                                        \
                cusparseIndexBase_t idxBase )                                                                                       \
        {                                                                                                                       \
            SCAI_CUSPARSE_CALL(                                                                                                 \
                    CUSPARSE_BLAS_NAME( csr2csc, prefix )( handle, m, n, nnz,                                           \
                            reinterpret_cast<const CUSPARSEValueType*>( csrVal ),        \
                            csrRowPtr, csrColInd,                                        \
                            reinterpret_cast<CUSPARSEValueType*>( cscVal ), cscRowInd,   \
                            cscColPtr, copyValues, idxBase ),                            \
                    "CUSPARSEWrapper::csr2csc" )                                                                        \
        }                                                                                                                       \
        \
        static void csrmv(                                                                                                      \
                BLASHandle handle,                                                                                                  \
                BLASTrans transA,                                                                                                   \
                BLASIndexType m,                                                                                                    \
                BLASIndexType n,                                                                                                    \
                BLASIndexType nnz,                                                                                                  \
                const ValueType *alpha,                                                                                             \
                const BLASMatrix descrA,                                                                                            \
                const ValueType *csrValA,                                                                                           \
                const BLASIndexType *csrRowPtrA,                                                                                    \
                const BLASIndexType *csrColIndA,                                                                                    \
                const ValueType *x,                                                                                                 \
                const ValueType *beta,                                                                                              \
                ValueType *y )                                                                                                      \
        {                                                                                                                       \
            SCAI_CUSPARSE_CALL(                                                                                                 \
                    CUSPARSE_BLAS_NAME( csrmv, prefix )( handle, transA, m, n, nnz,                                     \
                            reinterpret_cast<const CUSPARSEValueType*>( alpha ), descrA,   \
                            reinterpret_cast<const CUSPARSEValueType*>( csrValA ),         \
                            csrRowPtrA, csrColIndA,                                        \
                            reinterpret_cast<const CUSPARSEValueType*>( x ),               \
                            reinterpret_cast<const CUSPARSEValueType*>( beta ),            \
                            reinterpret_cast<CUSPARSEValueType*>( y ) ),                   \
                    "CUSPARSEWrapper::csrmv" )                                                                          \
        }                                                                                                                       \
        \
        static void csrgeam(                                                                                                    \
                BLASHandle handle,                                                                                                  \
                BLASIndexType m,                                                                                                    \
                BLASIndexType n,                                                                                                    \
                const ValueType *alpha,                                                                                             \
                const BLASMatrix descrA,                                                                                            \
                BLASIndexType nnzA,                                                                                                 \
                const ValueType *csrValA,                                                                                           \
                const BLASIndexType *csrRowPtrA,                                                                                    \
                const BLASIndexType *csrColIndA,                                                                                    \
                const ValueType *beta,                                                                                              \
                const BLASMatrix descrB,                                                                                            \
                BLASIndexType nnzB,                                                                                                 \
                const ValueType *csrValB,                                                                                           \
                const BLASIndexType *csrRowPtrB,                                                                                    \
                const BLASIndexType *csrColIndB,                                                                                    \
                const BLASMatrix descrC,                                                                                            \
                ValueType *csrValC,                                                                                                 \
                BLASIndexType *csrRowPtrC,                                                                                          \
                BLASIndexType *csrColIndC )                                                                                         \
        {                                                                                                                       \
            SCAI_CUSPARSE_CALL(                                                                                                 \
                    CUSPARSE_BLAS_NAME( csrgeam, prefix )( handle, m, n,                                                \
                            reinterpret_cast<const CUSPARSEValueType*>( alpha ),         \
                            descrA, nnzA,                                                \
                            reinterpret_cast<const CUSPARSEValueType*>( csrValA ),       \
                            csrRowPtrA, csrColIndA,                                      \
                            reinterpret_cast<const CUSPARSEValueType*>( beta ),          \
                            descrB, nnzB,                                                \
                            reinterpret_cast<const CUSPARSEValueType*>( csrValB ),       \
                            csrRowPtrB, csrColIndB, descrC,                              \
                            reinterpret_cast<CUSPARSEValueType*>( csrValC ),             \
                            csrRowPtrC, csrColIndC ),                                    \
                    "CUSPARSEWrapepr::csrgeam" )                                                                        \
        }                                                                                                                       \
        \
        static void csrgemm(                                                                                                    \
                BLASHandle handle,                                                                                                  \
                BLASTrans transA,                                                                                                   \
                BLASTrans transB,                                                                                                   \
                IndexType m,                                                                                                        \
                IndexType n,                                                                                                        \
                IndexType k,                                                                                                        \
                const BLASMatrix descrA,                                                                                            \
                const IndexType nnzA,                                                                                               \
                const ValueType *csrValA,                                                                                           \
                const IndexType *csrRowPtrA,                                                                                        \
                const IndexType *csrColIndA,                                                                                        \
                const BLASMatrix descrB,                                                                                            \
                const IndexType nnzB,                                                                                               \
                const ValueType *csrValB,                                                                                           \
                const IndexType *csrRowPtrB,                                                                                        \
                const IndexType *csrColIndB,                                                                                        \
                const BLASMatrix descrC,                                                                                            \
                ValueType *csrValC,                                                                                                 \
                const IndexType *csrRowPtrC,                                                                                        \
                IndexType *csrColIndC )                                                                                             \
        {                                                                                                                       \
            SCAI_CUSPARSE_CALL(                                                                                                 \
                    CUSPARSE_BLAS_NAME( csrgemm, prefix )( handle, transA, transB, m, n, k, descrA, nnzA,               \
                            reinterpret_cast<const CUSPARSEValueType*>( csrValA ),       \
                            csrRowPtrA, csrColIndA, descrB, nnzB,                        \
                            reinterpret_cast<const CUSPARSEValueType*>( csrValB ),       \
                            csrRowPtrB, csrColIndB, descrC,                              \
                            reinterpret_cast<CUSPARSEValueType*>( csrValC ),             \
                            csrRowPtrC, csrColIndC ),                                    \
                    "CUSPARSEWrapper::csrgemm" )                                                                        \
        }                                                                                                                       \
        \
    };

CUSPARSEWRAPPER_DEF( float, float, S )
CUSPARSEWRAPPER_DEF( double, double, D )

#ifdef SCAI_COMPLEX_SUPPORTED
CUSPARSEWRAPPER_DEF( ComplexFloat, cuFloatComplex, C )
CUSPARSEWRAPPER_DEF( ComplexDouble, cuDoubleComplex, Z )
#endif

#undef CUSPARSEWRAPPER_DEF

} /* end namespace sparsekernel */

} /* end namespace scai */

