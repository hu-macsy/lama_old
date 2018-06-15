/**
 * @file DIAKernelTrait.hpp
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
 * @brief Struct with traits for all DIA storage methods provided as kernels.
 * @author Thomas Brandes
 * @date 23.10.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/MatrixOp.hpp>

namespace scai
{

namespace sparsekernel
{

/** Traits for LAMA kernels used in DIA storage.  */

struct DIAKernelTrait
{
    struct getValuePos
    {
        /** Returns position of element (i,j) in values array
         *
         *  @param[in] i is the row of the element
         *  @param[in] j is the column of the element
         *  @param[in] diaOffsets diagonal offsets, size is numDiagonals
         *  @param[in] numDiagonals number of stored diagonals
         *  @returns  offset of element in values array, invalidIndex if not found
         */

        typedef IndexType ( *FuncType ) (
            const IndexType i,
            const IndexType j,
            const IndexType numRows,
            const IndexType diaOffsets[],
            const IndexType numDiagonals );

        static const char* getId()
        {
            return "DIA.getValuePos";
        }
    };

    template<typename ValueType>
    struct getCSRSizes
    {
        /** Type definition of function pointer for counting sparse values in DIA storage
         *
         *  @param[out] csrSizes array with number of non-zero entries in each row
         *  @param[in] numRows is the number of rows
         *  @param[in] numColumns is the number of columns
         *  @param[in] numDiagonals number of diagonals used in the DIA format
         *  @param[in] diaOffsets diagonal offsets, size is numDiagonals
         *  @param[in] diaValues are stored values of the diagonals
         *
         *  Note: the diagonals might contain zero entries so the number of non-zero
         *        elements might be less than number of stored elements
         *
         *  - csrSizes must have been allocated with at least numRows entries
         *  - diaOffsets has at least numDiagonals entries
         *  - diaValues has numDiagonals x max(numRows, numColumns) entries
         */

        typedef void ( *FuncType )(
            IndexType csrSizes[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const ValueType diaValues[] );

        static const char* getId()
        {
            return "DIA.getCSRSizes";
        }
    };

    template<typename ValueType>
    struct getCSRValues
    {
        /** Type definition of function pointer for conversion of DIA storage data to CSR data.
         *
         *  @param[out] csrJA will contain the column indexes
         *  @param[out] csrValues will contain the matrix elements
         *  @param[in] csrIA is the array with the offsets (must already be available before)
         *  @param[in] numRows is the number of rows
         *  @param[in] numColumns is the number of columns
         *  @param[in] numDiagonals number of diagonals used in the DIA format
         *  @param[in] diaOffsets diagonal offsets, size is numDiagonals
         *  @param[in] diaValues are stored values of the diagonals
         *
         *   - csrIA has numRows + 1 entries
         *   - csrJA and csrValues must have at least numValues entries, numValues = csrIA[numRows]
         */
        typedef void ( *FuncType ) (
            IndexType csrJA[],
            ValueType csrValues[],
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const ValueType diaValues[] );

        static const char* getId()
        {
            return "DIA.getCSRValues";
        }
    };

    template<typename ValueType>
    struct normalGEMV
    {
        /** result = alpha * DIA-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numValues is the number of diagonals in DIA storage
         *  @param diaOffsets, diaValues are arrays of DIA storage
         *  @param op specifies an operation implicitly applied to the matrix storage
         */

        typedef void ( *FuncType ) (
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const ValueType beta,
            const ValueType y[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const ValueType diaValues[],
            const common::MatrixOp op );

        static const char* getId()
        {
            return "DIA.normalGEMV";
        }
    };

    /** Structure with type definitions for solver routines */

    template<typename ValueType>
    struct jacobi
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         */
        typedef void ( *FuncType ) (
            ValueType solution[],
            const IndexType n,
            const IndexType numDiagonals,
            const IndexType diaOffset[],
            const ValueType diaValues[],
            const ValueType oldSolution[],
            const ValueType rhs[],
            const ValueType omega );

        static const char* getId()
        {
            return "DIA.jacobi";
        }
    };

    template<typename ValueType>
    struct jacobiHalo
    {
        /** Compute one iteration step in Jacobi method for halo
         *
         *  \code
         *      solution -= omega * ( dia_halo * oldSolution ) ./ diagonal 
         *  \endcode
         *
         */
        typedef void ( *FuncType ) (
            ValueType solution[],
            const ValueType diagonal[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffset[],
            const ValueType diaValues[],
            const ValueType oldSolution[],
            const ValueType omega );

        static const char* getId()
        {
            return "DIA.jacobiHalo";
        }
    };

    /** Structure with type definitions for reduction routines */

    template<typename ValueType>
    struct absMaxVal
    {
        /** This method returns the maximal absolute value of a DIA storage. */

        typedef ValueType ( *FuncType ) (
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const ValueType diaValues[]
        );

        static const char* getId()
        {
            return "DIA.absMaxVal";
        }
    };
};

} /* end namespace sparsekernel */

} /* end namespace scai */
