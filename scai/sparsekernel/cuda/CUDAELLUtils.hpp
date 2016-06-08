/**
 * @file CUDAELLUtils.hpp
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
 * @brief General conversion routines for ELL sparse matrices.
 * @author Thomas Brandes
 * @date 03.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai library
#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace sparsekernel
{

/** This class provides routines to converse ELL storage data to CSR storage data and vice versa.
 *
 *  All routines work on already allocated data and utilize CUDA for their parallelization.
 */

class COMMON_DLL_IMPORTEXPORT CUDAELLUtils
{
public:

    /** Addressing function for the arrays ia and ja: column-wise */

    static inline IndexType ellindex( const IndexType i, const IndexType jj, const IndexType numRows )
    {
        return jj * numRows + i;
    }

    /** This method computes the total number of non-zero rows by the size array  */

    static IndexType countNonEmptyRowsBySizes( const IndexType sizes[], const IndexType numRows );

    /** check diagonal property. ELL format with diagonal property: diagonal is just the first column in mValues */

    static bool hasDiagonalProperty( const IndexType numDiagonals, const IndexType ellJA[] );

    /** Build a vector of indexes for non-empty rows. */

    static void setNonEmptyRowsBySizes(
        IndexType rowIndexes[],
        const IndexType numNonEmptyRows,
        const IndexType sizes[],
        const IndexType numRows );

    static void check(
        const IndexType mNumRows,
        const IndexType mNumValuesPerRow,
        const IndexType mNumColumns,
        const IndexType* ia,
        const IndexType* ja,
        const char* msg );

    /** Returns one row of the matrix */

    template<typename ValueType, typename OtherValueType>
    static void getRow(
        OtherValueType row[],
        const IndexType i,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Returns one value of the matrix */

    template<typename ValueType>
    static ValueType getValue(
        const IndexType i,
        const IndexType j,
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Scales matrix using an vector */

    template<typename ValueType, typename OtherValueType>
    static void scaleValue(
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        ValueType ellValues[],
        const OtherValueType values[] );

    /** Implementation for ELLKernelTrait::Conversions::getCSRValues */

    template<typename ELLValueType, typename CSRValueType>
    static void getCSRValues(
        IndexType csrJA[],
        CSRValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ELLValueType ellValues[] );

    template<typename ValueType>
    static void fillELLValues(
        IndexType ellJA[],
        ValueType ellValues[],
        const IndexType ellSizes[],
        const IndexType numRows,
        const IndexType numValuesPerRow );

    /** Helper routine for conversion CSR to ELL format.  */

    template<typename ELLValueType, typename CSRValueType>
    static void setCSRValues(
        IndexType ellJA[],
        ELLValueType ellValues[],
        const IndexType ellSizes[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const CSRValueType csrValues[] );

    /** Implementation for ELLKernelTrait::normalGEMV on CUDA devices. */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellIA[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Implementation for ELLKernelTrait::normalGEVM  */

    template<typename ValueType>
    static void normalGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Implementation for ELLKernelTrait::sparseGEMV  */

    template<typename ValueType>
    static void sparseGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Implementation for ELLKernelTrait::sparseGEVM  */

    template<typename ValueType>
    static void sparseGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numNonZerosPerRow,
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType ellIA[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Implementation for ELLKernelTrait::jacobi  */

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

    /** Implementation for ELLKernelTrait::Solver::jacobiHalo  */

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

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Routine that registers all methods at the kernel registry. */

    SCAI_KREGISTRY_DECL_REGISTRATOR( Registrator )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorVO, template<typename ValueType, typename OtherValueType> )

    /** Constructor for registration. */

    CUDAELLUtils();

    /** Destructor for unregistration. */

    ~CUDAELLUtils();

    /** Static variable for registration at static initialization. */

    static CUDAELLUtils guard;
};

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
