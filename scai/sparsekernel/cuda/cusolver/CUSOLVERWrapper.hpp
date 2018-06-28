/**
 * @file sparsekernel/cuda/cusolver/CUSOLVERWrapper.hpp
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
 * @brief Wrapping for cusolver calls
 * @author Lauretta Schubert
 * @date 18.07.2016
 */

#pragma once

#include <scai/sparsekernel/cuda/cusolver/CUSOLVERTrait.hpp>

#include <scai/common/cuda/CUDAError.hpp>

// CUDA
#ifdef SCAI_COMPLEX_SUPPORTED
#include <cuComplex.h>
#endif

#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#elif ( CUDART_VERSION >= 7050 )

#include <cusolverSp.h>

namespace scai
{

namespace sparsekernel
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT CUSOLVERWrapper;

#define CUSOLVERWRAPPER_DEF( ValueType, CUSOLVERValueType, prefix )                                                             \
    template<>                                                                                                                  \
    class COMMON_DLL_IMPORTEXPORT CUSOLVERWrapper<ValueType>                                                                    \
    {                                                                                                                           \
    public:                                                                                                                     \
        typedef CUSOLVERTrait::SOLVERDnHandle SOLVERDnHandle;                                                                   \
        typedef CUSOLVERTrait::SOLVERSpHandle SOLVERSpHandle;                                                                   \
        typedef CUSOLVERTrait::BLASMatrix BLASMatrix;                                                                           \
        \
        static void csrQR(                                                                                                      \
                SOLVERSpHandle handle,                                                                                          \
                IndexType n,                                                                                                    \
                IndexType nnz,                                                                                                  \
                const BLASMatrix descrA,                                                                                        \
                const ValueType *csrValA,                                                                                       \
                const IndexType *csrRowPtrA,                                                                                    \
                const IndexType *csrColIndA,                                                                                    \
                const ValueType *b,                                                                                             \
                ValueType tol,                                                                                                  \
                IndexType reorder,                                                                                              \
                ValueType *x,                                                                                                   \
                IndexType *singularity )                                                                                        \
        {                                                                                                                       \
            SCAI_CUSOLVER_CALL(                                                                                                 \
                    CUSOLVER_SP_NAME( csrlsvqr, prefix )( handle, n, nnz, descrA,                                               \
                            reinterpret_cast<const CUSOLVERValueType*>( csrValA ),                                              \
                            csrRowPtrA, csrColIndA,                                                                             \
                            reinterpret_cast<const CUSOLVERValueType*>( b ),                                                    \
                            tol, reorder,                                                                                       \
                            reinterpret_cast<CUSOLVERValueType*>( x ),                                                          \
                            singularity ),                                                                                      \
                    "CUSOLVERWrapper::csrlsvlu" )                                                                               \
        }                                                                                                                       \
        \
    };

CUSOLVERWRAPPER_DEF( float, float, S )
CUSOLVERWRAPPER_DEF( double, double, D )

#ifdef SCAI_COMPLEX_SUPPORTED
CUSOLVERWRAPPER_DEF( ComplexFloat, cuFloatComplex, C )
CUSOLVERWRAPPER_DEF( ComplexDouble, cuDoubleComplex, Z )
#endif

#undef CUSOLVERWRAPPER_DEF

} /* end namespace sparsekernel */

} /* end namespace scai */

#endif
