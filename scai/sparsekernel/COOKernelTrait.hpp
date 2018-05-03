/**
 * @file COOKernelTrait.hpp
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
 * @brief Struct with traits for all COO storage methods provided as kernels.
 * @author Thomas Brandes
 * @date 22.10.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/MatrixOp.hpp>

namespace scai
{

namespace sparsekernel
{

/** Kernel traits for functions to be used in COO storage.
 *
 *  Note: routines to build CSR data from COO data are not required any more
 *        as this is done by bucket sort
 */

struct COOKernelTrait
{
    struct getValuePos
    {
        /** Returns position of element (i,j) in ja/values array
         *
         *  @param[in] i is the row of the element
         *  @param[in] j is the column of the element
         *  @param[in] cooIA is the COO ia array
         *  @param[in] cooJA is the COO ja array
         *  @returns  offset of element in values array, invalidIndex if not found
         */

        typedef IndexType ( *FuncType ) (
            const IndexType i,
            const IndexType j,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.getValuePos";
        }
    };

    struct getValuePosRow
    {
        /** This method returns for a certain row of the COO matrix all
         *  col indexes for which elements exist and the corresponding positions
         *  in the cooIA/cooJA/cooValues array
         *
         *  @param[out] col indexes of cols that have an entry for row i
         *  @param[out] pos positions of entries with row = i in cooXXX array
         *  @param[in] i is the row of which positions are required
         *  @param[in] cooIA is the COO array with row indexes
         *  @param[in] numColumns is the number of columns
         *  @param[in] cooJA is the COO array with col indexes
         *  @param[in] numValues is the number of non-zero values
         *  @returns  number of entries with row index = i
         */
        typedef IndexType ( *FuncType ) (
            IndexType col[],
            IndexType pos[],
            const IndexType i,
            const IndexType cooIA[],
            const IndexType numColumns,
            const IndexType cooJA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.getValuePosRow";
        }
    };

    struct getValuePosCol
    {
        /** This method returns for a certain column of the COO matrix all
         *  row indexes for which elements exist and the corresponding positions
         *  in the cooIA/cooJA/cooValues array
         *
         *  @param[out] row indexes of rows that have an entry for column j
         *  @param[out] pos positions of entries with col = j in cooJA,
         *  @param[in] j is the column of which positions are required
         *  @param[in] cooIA is the COO array with row indexes
         *  @param[in] numRows is the number of rows
         *  @param[in] cooJA is the COO array with col indexes
         *  @param[in] numValues is the number of non-zero values
         *  @returns  number of entries with col index = j
         */
        typedef IndexType ( *FuncType ) (
            IndexType row[],
            IndexType pos[],
            const IndexType j,
            const IndexType cooIA[],
            const IndexType numRows,
            const IndexType cooJA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.getValuePosCol";
        }
    };

    struct hasDiagonalProperty
    {
        /** Routine checks for diagonal property, first n entries are the diagonal elements.
         *
         *  @param[in] cooIA row indexes
         *  @param[in] cooJA column indexes
         *  @param[in] n number of diagonal elements
         *  @return true if first n entries stand for the diagonal elements
         *
         *  Attention: do not call this routine if n > numValues (size of cooIA, cooJA) where
         *             diagonal property is already false
         */

        typedef bool ( *FuncType )(
            const IndexType cooIA[],
            const IndexType cooJA[],
            const IndexType n );

        static const char* getId()
        {
            return "COO.hasDiagonalProperty";
        }
    };

    struct offsets2ia
    {
        /** Routine for conversion of CSR offset array to COO ia array
         *
         *  @param[out] cooIA is the array with all row indexes
         *  @param[in] numValues number of non-zero values, size of cooIA
         *  @param[in] csrIA is the CSR row offset array, size is numRows+1
         *  @param[in] numRows number of rows
         *  @param[in] numDiagonals is number of diagonals
         *
         *  The diagonal values will be stored at the beginning of the array cooIA.
         */

        typedef void ( *FuncType )(
            IndexType cooIA[],
            const IndexType numValues,
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType numDiagonals );

        static const char* getId()
        {
            return "COO.offsets2ia";
        }
    };

    template<typename COOValueType, typename OtherValueType>
    struct scaleRows
    {

        /** This operation multiplies each row with an own value.
         *
         *  @param[in,out] cooValues matrix data that is scaled
         *  @param[in]     cooIA are the row indexes
         *  @param[in]     rowValues array with scale factor for each row
         *  @param[in]     numValues number of entries in cooValues and cooIA
         *
         *  This routine supports different precision for matrix values and scale values.
         */

        typedef void ( *FuncType ) (
            COOValueType cooValues[],
            const OtherValueType rowValues[],
            const IndexType cooIA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.scaleRows";
        }
    };

    template<typename COOValueType, typename CSRValueType>
    struct setCSRData
    {

        /** Conversion of CSR data (ja, values) to COO data.
         *
         *  @param[out] cooValues is new COO data with diagonal elements at first
         *  @param[in] csrValues is given CSR data with diagonal elements first in each row
         *  @param[in] numValues is size of arrays cooValues and csrValues
         *  @param[in] csrIA is CSR offset array
         *  @param[in] numRows is number of rows
         *  @param[in] numDiagonals is number of diagonal elements to be stored at the beginning
         *  @param[in] csrDiagonalProperty is true if CSR data has diagonal property
         *
         *  Note: Diagonal elements must be first in each row for CSR data, no resort done for this
         *  Note: For numDiagonals == 0, this routine can be replaced with Utils::set.
         */

        typedef void ( *FuncType ) (
            COOValueType cooValues[],
            const CSRValueType csrValues[],
            const IndexType numValues,
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType numDiagonals );

        static const char* getId()
        {
            return "COO.setCSRData";
        }
    };

    template<typename ValueType>
    struct normalGEMV
    {
        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for vectors result and b and rows of matrix
         *  @param numColumns is number of elements for columns of matrix, size of result, y
         *  @param cooIA, cooJA, cooValues are arrays of COO storage
         *  @param numValues is the size of the coo arrays
         */

        typedef void ( *FuncType ) (
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const ValueType beta,
            const ValueType y[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType nnz,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const ValueType cooValues[],
            const common::MatrixOp op );

        static const char* getId()
        {
            return "COO.normalGEMV";
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
            ValueType* const solution,
            const IndexType cooNumValues,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const ValueType cooValues[],
            const ValueType oldSolution[],
            const ValueType rhs[],
            const ValueType omega,
            const IndexType numRows );

        static const char* getId()
        {
            return "COO.jacobi";
        }
    };
};

} /* end namespace sparsekernel */

} /* end namespace scai */
