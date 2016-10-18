/**
 * @file OpenMPCSRUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
#include <scai/utilskernel/BinaryOp.hpp>

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

    /** Implementation for CSRKernelTrait::getValuePosCol */

    static IndexType getValuePosCol( 
        IndexType row[], 
        IndexType pos[],
        const IndexType j, 
        const IndexType csrIA[], 
        const IndexType numRows,
        const IndexType csrJA[],
        const IndexType numValues );

    /** This method computes the total number of non-zero rows by the offset array
     *
     */

    static IndexType countNonEmptyRowsByOffsets( const IndexType offsets[], const IndexType numRows );

    /** Build a vector of indexes for non-empty rows. */

    static void setNonEmptyRowsByOffsets(
        IndexType rowIndexes[],
        const IndexType numNonEmptyRows,
        const IndexType offsets[],
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

    /** Implementation for CSRKernelTrait::Offsets::sizes2offsets */

    static IndexType sizes2offsets( IndexType sizes[], const IndexType numRows );

    /** Implementation for CSRKernelTrait::Offsets::offsets2sizes */

    static void offsets2sizes( IndexType sizes[], const IndexType offsets[], const IndexType n );

    /** offset2sizes for indexed rows */

    static void offsets2sizesGather(
        IndexType sizes[],
        const IndexType offsets[],
        const IndexType rowIndexes[],
        const IndexType numRows );

    /** Implementation for CSRKernelTrait::Offsets::validOffsets  */

    static bool validOffsets( const IndexType array[], const IndexType n, const IndexType total );

    /** This function checks if csr data has diagonal property.
     *
     *  @param[in] numDiagonals is min( numRows, numColumns )
     *  @param[in] csrIA is the csr offset array
     *  @param[in] csrJA is array with column indexes
     */

    static bool hasDiagonalProperty( const IndexType numDiagonals, const IndexType csrIA[], const IndexType csrJA[] );

    /** Host implementation for CSRKernelTrait::sortRowElements using OpenMP parallelization. */

    template<typename ValueType>
    static void sortRowElements(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const bool diagonalFlag );

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

    /** Implementation for CSRKernelTrait::scaleRows  */

    template<typename ValueType1, typename ValueType2>
    static void scaleRows(
        ValueType1 csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const ValueType2 values[] );

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
        const IndexType nnz,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[] );

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
        const ValueType csrValues[] );

    /** Implementation for CSRKernelTrait::normalGEVM */

    template<typename ValueType>
    static void normalGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[] );

    /** Implementation for CSRKernelTrait::sparseGEVM  */

    template<typename ValueType>
    static void sparseGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numColumns,
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[] );

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
        const ValueType csrValues[] );

    /** Implementation for CSRKernelTrait::Jacobi::jacobi(Async/Halo) */

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

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const IndexType localIA[],
        const ValueType localValues[],
        const IndexType haloIA[],
        const IndexType haloJA[],
        const ValueType haloValues[],
        const IndexType haloRowIndexes[],
        const ValueType oldSolution[],
        const ValueType omega,
        const IndexType numNonEmptyRows );

    /** Implementation for CSRKernelTrait::Jacobi::jacobiHaloWithDiag
     *  @since 1.1.0
     */

    template<typename ValueType>
    static void jacobiHaloWithDiag(
        ValueType solution[],
        const ValueType localDiagValues[],
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
        bool diagonalProperty,
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
        bool diagonalProperty,
        const IndexType aIA[],
        const IndexType aJA[],
        const IndexType bIA[],
        const IndexType bJA[] );

    /** Implementation for CSRKernelTrait::Offsets::matrixMultiplyJA  */

    static void matrixMultiplyJA(
        IndexType cJA[],
        const IndexType cIA[],
        const IndexType numRows,
        const IndexType numColumns,
        bool diagonalProperty,
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
        bool diagonalProperty,
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
        bool diagonalProperty,
        const IndexType aIA[],
        const IndexType aJA[],
        const ValueType aValues[],
        const IndexType bIA[],
        const IndexType bJA[],
        const ValueType bValues[] );

    /** Implementation for CSRKernelTrait::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal(
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
        const ValueType eps,
        const bool diagonalFlag );

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
        const ValueType eps,
        const bool diagonalFlag );

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

    /** Struct for registration of methods with two template arguments.
     *
     *  Registration function is wrapped in struct/class that can be used as template 
     *  argument for metaprogramming classes to expand for all supported types.
     */

    template<typename ValueType, typename OtherValueType>
    struct RegistratorVO
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
    static ValueType absMaxDiffRowSorted(
        const IndexType n1, const IndexType csrJA1[], const ValueType csrValues1[],
        const IndexType n2, const IndexType csrJA2[], const ValueType csrValues2[] );

    template<typename ValueType>
    static ValueType absMaxDiffRowUnsorted(
        const IndexType n1, const IndexType csrJA1[], const ValueType csrValues1[],
        const IndexType n2, const IndexType csrJA2[], const ValueType csrValues2[] );

    /** Help routine to use normalGEMV with maximal 9 args. */

    template<typename ValueType>
    static void normalGEMV_s(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[] );

    template<typename ValueType>
    static void normalGEVM_s(
        ValueType result[],
        std::pair<ValueType, const ValueType*> ax,
        std::pair<ValueType, const ValueType*> by,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[] );
};

} /* end namespace sparsekernel */

} /* end namespace scai */
