/**
 * @file OpenMPELLUtils.hpp
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
 * @brief General conversion routines for ELL sparse matrices.
 * @author Thomas Brandes
 * @date 03.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries

#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/TypeTraits.hpp>

#include <utility>

namespace scai
{

namespace sparsekernel
{

/** This class provides routines to converse ELL storage data to CSR storage data and vice versa.
 *
 *  All routines work on already allocated data and utilize OpenMP for their parallelization.
 */

class COMMON_DLL_IMPORTEXPORT OpenMPELLUtils
{

public:

    /** Sort the entries of one row in ascending order */

    template<typename ValueType>
    static void sortRowElements(
        IndexType ellJA[],
        ValueType ellValues[],
        const IndexType ellIA[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const bool diagonalFlag );

private:

    /** Addressing function for the arrays ellJA[numRows*numValuesPerRow] and ellValues: column-major order */

    static inline IndexType ellindex(
        const IndexType i,
        const IndexType jj,
        const IndexType numRows,
        const IndexType /* numValuesPerRow */ )
    {
        return jj * numRows + i; // column major-order
        // return i * numValuesPerRow + jj;    // row major-order
    }

    /** Returns the maximal absolute value of the ELL storage. */

    template<typename ValueType>
    static ValueType absMaxVal(
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const ValueType ellValues[] );

    static void check(
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType numColumns,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const char* msg );

    /** Returns one row of the matrix */

    template<typename ValueType>
    static void getRow(
        ValueType row[],
        const IndexType i,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Host implementation for ELLKernelTrait::getValuePos */

    static IndexType getValuePos(
        const IndexType i,
        const IndexType j,
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[] );

    /** Implementation for ELLKernelTrait::getDiagonalPositions */

    static IndexType getDiagonalPositions(
        IndexType diagonalPositions[],
        const IndexType numDiagonals,
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[] );

    /** Implementation for ELLKernelTrait::getColumnPositions */

    static IndexType getColumnPositions(
        IndexType row[],
        IndexType pos[],
        const IndexType j,
        const IndexType ellIA[],
        const IndexType numRows,
        const IndexType ellJA[],
        const IndexType numValuesPerRow );

    /** Implementation for ELLKernelTrait::setRows */

    template<typename ValueType>
    static void setRows(
        ValueType ellValues[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const ValueType values[],
        const common::BinaryOp op );

    /** Implementation for ELLKernelTrait::setColumns */

    template<typename ValueType>
    static void setColumns(
        ValueType ellValues[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType values[],
        const common::BinaryOp op );

    /** Implementation for ELLKernelTrait::compressIA */

    template<typename ValueType>
    static void compressIA(
        IndexType newIA[],
        const IndexType ellIA[],
        const IndexType ellJA[],
        const ValueType ellValues[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const RealType<ValueType> eps );

    /** Implementation for ELLKernelTrait::compressValues */

    template<typename ValueType>
    static void compressValues(
        IndexType newJA[],
        ValueType newValues[],
        const IndexType newNumValuesPerRow,
        const IndexType ellIA[],
        const IndexType ellJA[],
        const ValueType ellValues[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const RealType<ValueType> eps );

    /** Implementation for ELLKernelTrait::getCSRValues */

    template<typename ValueType>
    static void getCSRValues(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    template<typename ValueType>
    static void fillELLValues(
        IndexType ellJA[],
        ValueType ellValues[],
        const IndexType ellSizes[],
        const IndexType numRows,
        const IndexType numValuesPerRow );

    /** Implementation for ELLKernelTrait::setCSRValues */

    template<typename ValueType>
    static void setCSRValues(
        IndexType ellJA[],
        ValueType ellValues[],
        const IndexType ellSizes[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[] );

    /** Implementation for ELLKernelTrait::matrixMultiplySizes */

    static void matrixMultiplySizes(
        IndexType cSizes[],
        const IndexType m,
        const IndexType n,
        const IndexType k,
        const bool diagonalProperty,
        const IndexType aSizes[],
        const IndexType aJA[],
        const IndexType aNumValuesPerRow,
        const IndexType bSizes[],
        const IndexType bJA[],
        const IndexType bNumValuesPerRow );

    /** Implementation for ELLKernelTrait::matrixMultiply */

    template<typename ValueType>
    static void matrixMultiply(
        IndexType cJA[],
        ValueType cValues[],
        const IndexType cSizes[],
        const IndexType cNumValuesPerRow,
        const IndexType m,
        const IndexType n,
        const IndexType k,
        const bool diagonalProperty,
        const ValueType alpha,
        const IndexType aSizes[],
        const IndexType aJA[],
        const ValueType aValues[],
        const IndexType aNumValuesPerRow,
        const IndexType bSizes[],
        const IndexType bJA[],
        const ValueType bValues[],
        const IndexType bNumValuesPerRow );

    /** Implementation for ELLKernelTrait::matrixAddSizes */

    static void matrixAddSizes(
        IndexType csizes[],
        const IndexType m,
        const IndexType n,
        const bool diagonalProperty,
        const IndexType aSizes[],
        const IndexType aJA[],
        const IndexType aNumValuesPerRow,
        const IndexType bSizes[],
        const IndexType bJA[],
        const IndexType bNumValuesPerRow );

    /** Implementation for ELLKernelTrait::matrixAdd */

    template<typename ValueType>
    static void matrixAdd(
        IndexType cJA[],
        ValueType cValues[],
        const IndexType cSizes[],
        const IndexType cNumValuesPerRow,
        const IndexType m,
        const IndexType n,
        const bool diagonalProperty,
        const ValueType alpha,
        const IndexType aSizes[],
        const IndexType aJA[],
        const ValueType aValues[],
        const IndexType aNumValuesPerRow,
        const ValueType beta,
        const IndexType bSizes[],
        const IndexType bJA[],
        const ValueType bValues[],
        const IndexType bNumValuesPerRow );

    /** Implementation for ELLKernelTrait::jacobi */

    template<typename ValueType>
    static void jacobi(
        ValueType solution[],
        const IndexType numRows,
        const IndexType ellNumValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[],
        const ValueType oldSolution[],
        const ValueType rhs[],
        const ValueType omega );

    /** Implementation for ELLKernelTrait::jacobiHalo */

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const IndexType numRows,
        const ValueType diagonal[],
        const IndexType ellNumValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[],
        const IndexType rowIndexes[],
        const IndexType numNonEmptyRows,
        const ValueType oldSolution[],
        const ValueType omega );

    /** Implementation for ELLKernelTrait::normalGEMV  */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numNonZerosPerRow,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        const common::MatrixOp op );

    /** Implementation for ELLKernelTrait::sparseGEMV  */

    template<typename ValueType>
    static void sparseGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numRows,
        const IndexType numNonZerosPerRow,
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        const common::MatrixOp op );

private:

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

    OpenMPELLUtils();

    /** Destructor for unregistration. */

    ~OpenMPELLUtils();

    /** Static variable for registration at static initialization. */

    static OpenMPELLUtils guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
