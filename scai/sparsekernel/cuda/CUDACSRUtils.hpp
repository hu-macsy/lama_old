/**
 * @file CUDACSRUtils.hpp
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
 * @brief General conversion routines for CSR sparse matrices implemented in CUDA.
 * @author Thomas Brandes
 * @date 03.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// scai internal libraries
#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/MatrixOp.hpp>

#include <cuda_runtime_api.h>

#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#endif

namespace scai
{

namespace sparsekernel
{

/** Static class that provides CUDA implementaions for the routines
 *  needed for operations on CSR storage as specified in CSRKernelTrait.
 *  Routines will be registered at KernelRegistry during the static initialization.
 */
class COMMON_DLL_IMPORTEXPORT CUDACSRUtils
{
public:
    /** Helper routine for conversion CSR to CSR format.  */

    static IndexType sizes2offsets( IndexType sizes[], const IndexType numRows );

    static void offsets2sizes( IndexType sizes[], const IndexType offsets[], const IndexType n );

    /** Implementation for CSRKernelTrait::getValuePosCol */

    static IndexType getValuePosCol(
        IndexType row[],
        IndexType pos[],
        const IndexType j,
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType csrJA[],
        const IndexType numValues );

    static bool hasDiagonalProperty( 
        const IndexType numDiagonals, 
        const IndexType csrIA[], 
        const IndexType csrJA[],
        const bool );

    /** Matrix transpose for CSR matrices on CUDA device. */

    template<typename ValueType>
    static void convertCSR2CSC(
        IndexType cscIA[],
        IndexType cscJA[],
        ValueType cscValues[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        IndexType numRows,
        IndexType numColumns,
        IndexType numValues );

    /** Implementation for CSRKernelTrait::Mult::scaleRows  */

    template<typename ValueType>
    static void scaleRows(
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const ValueType values[] );

    /** Implementation for CSRKernelTrait::Mult::normalGEMV  */

    template<typename ValueType>
    static void normalGEMV(
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
        const common::MatrixOp op );

    /** Implementation for CSRKernelTrait::Mult::sparseGEMV  */

    template<typename ValueType>
    static void sparseGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        const common::MatrixOp op );

    /** Implementation for CSRKernelTrait::Solver::jacobi  */

    template<typename ValueType>
    static void jacobi(
        ValueType* const solution,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        const ValueType rhs[],
        const ValueType oldSolution[],
        const ValueType omega,
        const IndexType numRows );

    /** Implementation for CSRKernelTrait::Solver::jacobiHalo
     */

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const ValueType localDiagonal[],
        const IndexType haloIA[],
        const IndexType haloJA[],
        const ValueType haloValues[],
        const IndexType haloRowIndexes[],
        const ValueType oldSolution[],
        const ValueType omega,
        const IndexType numNonEmptyRows );

    /** Implementation for CSRKernelTrait::Offsets::matrixAddSizes  */

    static IndexType matrixAddSizes(
        IndexType cSizes[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType aIA[],
        const IndexType aJA[],
        const IndexType bIA[],
        const IndexType bJA[] );

    /** Implementation for CSRKernelTrait::Offsets::matrixMultiplySizes  */

    static IndexType matrixMultiplySizes(
        IndexType cSizes[],
        const IndexType m,
        const IndexType n,
        const IndexType k,
        const IndexType aIA[],
        const IndexType aJA[],
        const IndexType bIA[],
        const IndexType bJA[] );

    /** Implementation for CSRKernelTrait::Mult::matrixAdd */

    template<typename ValueType>
    static void matrixAdd(
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
        const ValueType bValues[] );

    /** Implementation for CSRKernelTrait::Mult::matrixMultiply */

    template<typename ValueType>
    static void matrixMultiply(
        const IndexType cIa[],
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
        const ValueType bValues[] );

    /** CUDA Implementation for CSRKernelTrait::setDiagonalFirst  */

    template<typename ValueType>
    static IndexType setDiagonalFirst(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType numDiagonals,
        const IndexType csrIA[] );

    /** CUDA implementation for CSRKernelTrait::sortRows */

    template<typename ValueType>
    static void sortRows(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues );

    /** CUDA implementation for CSRKernelTrait::hasSortedRows */

    static bool hasSortedRows(
        const IndexType csrIA[],
        const IndexType csrJA[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType nnz );

    /** CUDA Implementation for CSRKernelTrait::countNonZeros */

    template<typename ValueType>
    static void countNonZeros(
        IndexType sizes[],
        const IndexType ia[],
        const IndexType ja[],
        const ValueType values[],
        const IndexType numRows,
        const RealType<ValueType> eps );

    /** CUDA Implementation for CSRKernelTrait::compress */

    template<typename ValueType>
    static void compress(
        IndexType newJA[],
        ValueType newValues[],
        const IndexType newIA[],
        const IndexType ia[],
        const IndexType ja[],
        const ValueType values[],
        const IndexType numRows,
        const RealType<ValueType> eps );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static unsigned int lastHashTableSize;// local variable to handhover hash table size for multiply

    /** Struct for registration of methods without template arguments */

    struct Registrator
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Struct for registration of methods with one template argument.
     *
     *  Registration function is wrapped in struct/class that can be used as template
     *  argument for metaprogramming classes to expand for each supported type
     */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUDACSRUtils();

    /** Destructor for unregistration. */

    ~CUDACSRUtils();

    /** Static variable for registration at static initialization. */

    static CUDACSRUtils guard;
};

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
