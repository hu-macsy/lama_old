/**
 * @file OpenMPCSRUtils.hpp
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
 * @brief Implementation of CSR utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/UnaryOp.hpp>
#include <scai/common/MatrixOp.hpp>

namespace scai
{

namespace sparsekernel
{

/** This class provides routines on compressed sparse row data
 */

class COMMON_DLL_IMPORTEXPORT OpenMPCSRUtils
{
public:

    /** Implementation for CSRKernelTrait::getValuePos */

    static IndexType getValuePos(
        const IndexType i,
        const IndexType j,
        const IndexType csrIA[],
        const IndexType csrJA[] );

    /** Implementation for CSRKernelTrait::getColumnPositions */

    static IndexType getColumnPositions(
        IndexType row[],
        IndexType pos[],
        const IndexType j,
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType csrJA[],
        const IndexType numValues );

    /** Implementation for CSRKernelTrait::nonEmptyRows */

    static IndexType nonEmptyRows( 
        IndexType rowIndexes[],
        const IndexType csrIA[], 
        const IndexType numRows );

    /** This function build an offset array from a counter array.
     *
     *  @param[in,out] array contains counter values and later the offsets
     *  @param[in]     numValues is the number of values, array must contain one additional value
     *  @returns       the total number of values
     *
     *  \code
     *    array  :    3    7   8   4   2
     *    array  :    0   10  15  12  16    -> returns 18
     *  \endcode
     *
     *  CSR sparse representation of matrices store the sum of all values at an additional
     *  position at the end of the array.
     */
    static IndexType scan( IndexType array[], const IndexType numValues );

    /** Implementation for CSRKernelTrait::sizes2offsets */

    static IndexType sizes2offsets( IndexType sizes[], const IndexType numRows );

    /** Implementation for CSRKernelTrait::offsets2sizes */

    static void offsets2sizes( IndexType sizes[], const IndexType offsets[], const IndexType n );

    /** Implementation for CSRKernelTrait::gatherSizes */

    static void gatherSizes(
        IndexType sizes[],
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType rowIndexes[],
        const IndexType nIndexes );

    /** Implementation for CSRKernelTrait::validOffsets  */

    static bool validOffsets( const IndexType array[], const IndexType n, const IndexType total );

    /** Implementation for CSRKernelTrait::hasDiagonalProperty using OpenMP parallelization. */

    static bool hasDiagonalProperty( 
        const IndexType numDiagonals, 
        const IndexType csrIA[], 
        const IndexType csrJA[], 
        const bool isSorted );

    /** Host implementation for CSRKernelTrait::sortRowElements using OpenMP parallelization. */

    template<typename ValueType>
    static void sortRows(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType nnz );

    /** Host implementation for CSRKernelTrait::hasSortedRows using OpenMP parallelization. */

    static bool hasSortedRows(
        const IndexType csrIA[],
        const IndexType csrJA[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType nnz );

    static IndexType getPosDiagonal(
        IndexType pos[],
        const IndexType numDiagonals,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const bool isSorted );

    /** Implementation for CSRKernelTrait::shiftDiagonal  */

    template<typename ValueType>
    static IndexType shiftDiagonal(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType numDiagonals,
        const IndexType csrIA[] );

    /** Implementation for CSRKernelTrait::convertCSR2CSC  */

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

    /** Implementation for CSRKernelTrait::setRows  */

    template<typename ValueType>
    static void setRows(
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const ValueType values[], 
        common::BinaryOp op );


    template<typename ValueType>
    static void reduce(
        ValueType result[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        const IndexType numRows,
        const IndexType dim,
        const common::BinaryOp reduceOp,
        const common::UnaryOp elemOp );

    /** Implementation for CSRKernelTrait::normalGEMV  */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[], 
        const common::MatrixOp op );

    /** Implementation for CSRKernelTrait::sparseGEMV  */

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

    /** Implementation for CSRKernelTrait::Mult::gemm  */

    template<typename ValueType>
    static void gemm(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType m,
        const IndexType n,
        const IndexType p,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[], 
        const common::MatrixOp op );

    /** Implementation for CSRKernelTrait::jacobi */

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

    /** Implementation for CSRKernelTrait::jacobiHalo  */

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const ValueType diagonal[],
        const IndexType haloIA[],
        const IndexType haloJA[],
        const ValueType haloValues[],
        const IndexType haloRowIndexes[],
        const ValueType oldSolution[],
        const ValueType omega,
        const IndexType numNonEmptyRows );

    /** Implementation for CSRKernelTrait::decomposition */

    template<typename ValueType>
    static void decomposition(
        ValueType* const solution,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        const ValueType rhs[],
        const IndexType numRows,
        const IndexType nnz,
        const bool isSymmetic );

    /** Implementation for CSRKernelTrait::Offsets::matrixAddSizes  */

    static IndexType matrixAddSizes(
        IndexType cSizes[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType aIA[],
        const IndexType aJA[],
        const IndexType bIA[],
        const IndexType bJA[] );

    /** Implementation for CSRKernelTrait::Offsets::binaryOpSizes  */

    static IndexType binaryOpSizes(
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

    /** binary operation for two CSR storages */

    template<typename ValueType>
    static void binaryOp(
        IndexType cJA[],
        ValueType cValues[],
        const IndexType cIA[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType aIA[],
        const IndexType aJA[],
        const ValueType aValues[],
        const IndexType bIA[],
        const IndexType bJA[],
        const ValueType bValues[],
        common::BinaryOp op );

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

    /** Implementation for CSRKernelTrait::absMaxDiffVal */

    template<typename ValueType>
    static RealType<ValueType> absMaxDiffVal(
        IndexType numRows,
        bool sortedRows,
        const IndexType csrIA1[],
        const IndexType csrJA1[],
        const ValueType csrValues1[],
        const IndexType csrIA2[],
        const IndexType csrJA2[],
        const ValueType csrValues2[] );

    /** Implementation for CSRKernelTrait::countNonZeros */

    template<typename ValueType>
    static void countNonZeros(
        IndexType sizes[],
        const IndexType ia[],
        const IndexType ja[],
        const ValueType values[],
        const IndexType numRows,
        const RealType<ValueType> eps );

    /** Implementation for CSRKernelTrait::compress */

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

    template<typename ValueType>
    static void getDiagonal(
        ValueType diagonal[],
        const IndexType numDiagonals,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        const bool isSorted );

    template<typename ValueType>
    static bool setDiagonalV(
        ValueType csrValues[],
        const ValueType diagonal[],
        const IndexType numDiagonals,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const bool isSorted );

    template<typename ValueType>
    static bool setDiagonal(
        ValueType csrValues[],
        const ValueType diagonal,
        const IndexType numDiagonals,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const bool isSorted );

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

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

    OpenMPCSRUtils();

    /** Destructor for unregistration. */

    ~OpenMPCSRUtils();

    /** Static variable for registration at static initialization. */

    static OpenMPCSRUtils guard;

    static IndexType scanSerial( IndexType array[], const IndexType numValues );

    static IndexType scanParallel( PartitionId numThreads, IndexType array[], const IndexType numValues );

    template<typename ValueType>
    static RealType<ValueType> absMaxDiffRowSorted(
        const IndexType n1, const IndexType csrJA1[], const ValueType csrValues1[],
        const IndexType n2, const IndexType csrJA2[], const ValueType csrValues2[] );

    template<typename ValueType>
    static RealType<ValueType> absMaxDiffRowUnsorted(
        const IndexType n1, const IndexType csrJA1[], const ValueType csrValues1[],
        const IndexType n2, const IndexType csrJA2[], const ValueType csrValues2[] );
};

} /* end namespace sparsekernel */

} /* end namespace scai */
