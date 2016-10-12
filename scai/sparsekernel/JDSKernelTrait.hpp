/**
 * @file JDSKernelTrait.hpp
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
 * @brief Struct with traits for all JDS storage methods provided as kernels.
 * @author Thomas Brandes
 * @date 23.10.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/utilskernel/ReductionOp.hpp>

namespace scai
{

namespace sparsekernel
{

/** Kernel traits for functions to be used in JDS storage.  */

struct JDSKernelTrait
{
    template<typename ValueType>
    struct jacobi
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         */
        typedef void ( *FuncType )(
            ValueType* const solution,
            const IndexType numRows,
            const IndexType jdsPerm[],
            const IndexType jdsIlg[],
            const IndexType jdsNumDiagonals,
            const IndexType jdsDlg[],
            const IndexType jdsJA[],
            const ValueType jdsValues[],
            const ValueType oldSolution[],
            const ValueType rhs[],
            const ValueType omega );

        static const char* getId()
        {
            return "JDS.jacobi";
        }
    };

    /** Structure with type definitions for solver routines */
    template<typename ValueType>
    struct jacobiHalo
    {
        /** Method to compute one iteration step in Jacobi method with halo.  */

        typedef void ( *FuncType )(
            ValueType solution[],
            const IndexType numRows,
            const ValueType invDiagonal[],
            const IndexType numDiagonals,
            const IndexType jdsHaloPerm[],
            const IndexType jdsHaloIlg[],
            const IndexType jdsHaloDlg[],
            const IndexType jdsHaloJA[],
            const ValueType jdsHaloValues[],
            const ValueType oldSolution[],
            const ValueType omega );

        static const char* getId()
        {
            return "JDS.jacobiHalo";
        }
    };

    struct sortRows
    {
        /** Stable sorting of values in array in descending order.
         *
         *  @param[in,out] array are the values to be sorted
         *  @param[in,out] perm, where perm[i] has the value of the original position
         *  @param[in]    n is the number of values to be sorted
         *
         *  \code
         *           array =   1  4   1  8  5  7
         *           perm  =   0  1   2  3  4  5
         +
         *           array =   8  7   5  4  1  1
         *           perm  =   3  5   4  1  0  2
         *  \endcode
         */

        typedef void ( *FuncType ) ( IndexType array[],
                                     IndexType perm[],
                                     const IndexType n );

        static const char* getId()
        {
            return "JDS.sortRows";
        }
    };

    struct setInversePerm
    {
        /** Compute the inverse permutation for a given permutation.
         *
         *  inversePerm [ perm [i] ] == i , 0 <= i < n
         *
         *  @param[out] inversePerm, size = n, will contain the inverse permutation
         *  @param[in] perm, size = n, is input permuation of 0, ..., n-1
         *  @param[in] n specifies the size of perm and inversePerm
         *
         *  /code
         *       perm      2  5  1  4  6  3  0
         *     inperm      6  2  0  5  3  1  4
         *  /endcode
         */

        typedef void ( *FuncType ) ( IndexType inversePerm[],
                                     const IndexType perm[],
                                     const IndexType n );

        static const char* getId()
        {
            return "JDS.setInversePerm";
        }
    };

    struct ilg2dlg
    {

        /** Compute dlg array from ilg array.
         *
         *  @param[out] dlg is the array with sizes of the columns
         *  @param[in]  numDiagonals is size of dlg
         *  @param[in]  ilg is the array with sizes of the rows
         *  @param[in]  numRows is the number of rows, size of ilg
         *
         *  The values in ilg must be descending. The same will be true
         *  for the output array dlg.
         *
         *  /code
         *       ilg       4  3  2  2  1   dlg
         *       5         x  x  x  x  x
         *       4         x  x  x  x
         *       2         x  x
         *       1         x
         *  /endcode
         */

        typedef IndexType ( *FuncType ) ( IndexType dlg[], const IndexType numDiagonals,
                                          const IndexType ilg[], const IndexType numRows );

        static const char* getId()
        {
            return "JDS.ilg2dlg";
        }
    };

    template<typename JDSValueType, typename CSRValueType>
    struct getCSRValues
    {
        /** Conversion of JDS storage data to CSR data
         *
         *  @param[out]  csrJA is array with column indexes
         *  @param[out]  csrValues is array with non-zero values
         *  @param[in]   csrIA is offset array (must be computed before)
         *  @param[in]   numRows number of rows in matrix
         *  @param[in]   jdsPerm with jdsPerm[ii] is original index of row i
         *  @param[in]   jdsILG with size of entries in row i
         *  @param[in]   jdsDLG distances of columns
         *  @param[in]   jdsJA column indexes
         *  @param[in]   jdsValues matrix values
         */

        typedef void ( *FuncType ) ( IndexType csrJA[],
                                     CSRValueType csrValues[],
                                     const IndexType csrIA[],
                                     const IndexType numRows,
                                     const IndexType jdsPerm[],
                                     const IndexType jdsILG[],
                                     const IndexType jdsDLG[],
                                     const IndexType jdsJA[],
                                     const JDSValueType jdsValues[] );

        static const char* getId()
        {
            return "JDS.getCSRValues";
        }
    };

    template<typename JDSValueType, typename CSRValueType>
    struct setCSRValues
    {

        /** Conversion of CSR storage data to JDS data
         *
         *  @param[out]  jdsJA column indexes
         *  @param[out]  jdsValues matrix values
         *  @param[in]   numRows number of rows in matrix
         *  @param[in]   jdsPerm with jdsPerm[ii] is original index of row i
         *  @param[in]   jdsILG with size of entries in row i
         *  @param[in]   numDiagonals size of array jdsDLG
         *  @param[in]   jdsDLG distances of columns
         *  @param[in]   csrIA is offset array (must be computed before)
         *  @param[in]   csrJA is array with column indexes
         *  @param[in]   csrValues is array with non-zero values
         */

        typedef void( *FuncType ) ( IndexType jdsJA[],
                                    JDSValueType jdsValues[],
                                    const IndexType numRows,
                                    const IndexType jdsPerm[],
                                    const IndexType jdsILG[],
                                    const IndexType numDiagonals,
                                    const IndexType jdsDLG[],
                                    const IndexType csrIA[],
                                    const IndexType csrJA[],
                                    const CSRValueType csrValues[] );

        static const char* getId()
        {
            return "JDS.setCSRValues";
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
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numDiagonals are number of diagonals, is size of jdsDLG
         *  @param jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues are arrays of JDS storage
         */

        typedef void ( *FuncType ) ( ValueType result[],
                                     const ValueType alpha,
                                     const ValueType x[],
                                     const ValueType beta,
                                     const ValueType y[],
                                     const IndexType numRows,
                                     const IndexType jdsPerm[],
                                     const IndexType jdsILG[],
                                     const IndexType ndlg,
                                     const IndexType jdsDLG[],
                                     const IndexType jdsJA[],
                                     const ValueType jdsValues[] );

        static const char* getId()
        {
            return "JDS.normalGEMV";
        }
    };

    template<typename ValueType>
    struct normalGEVM
    {
        typedef void ( *FuncType ) ( ValueType result[],
                                     const ValueType alpha,
                                     const ValueType x[],
                                     const ValueType beta,
                                     const ValueType y[],
                                     const IndexType numColumns,
                                     const IndexType jdsPerm[],
                                     const IndexType jdsILG[],
                                     const IndexType ndlg,
                                     const IndexType jdsDLG[],
                                     const IndexType jdsJA[],
                                     const ValueType jdsValues[] );

        static const char* getId()
        {
            return "JDS.normalGEVM";
        }
    };

    template<typename ValueType, typename OtherValueType>
    struct getRow
    {
        typedef void ( *FuncType ) ( OtherValueType row[],
                                     const IndexType i,
                                     const IndexType numColumns,
                                     const IndexType numRows,
                                     const IndexType perm[],
                                     const IndexType ilg[],
                                     const IndexType dlg[],
                                     const IndexType ja[],
                                     const ValueType values[] );

        static const char* getId()
        {
            return "JDS.getRow";
        }
    };

    struct getValuePos
    {
        typedef IndexType ( *FuncType ) ( const IndexType i,
                                          const IndexType j,
                                          const IndexType numRows,
                                          const IndexType ilg[],
                                          const IndexType dlg[],
                                          const IndexType perm[],
                                          const IndexType ja[] );

        static const char* getId()
        {
            return "JDS.getValuePos";
        }
    };

    struct getValuePosCol
    {
        /** This method returns for a certain column of the JDS matrix all
         *  row indexes for which elements exist and the corresponding positions
         *  in the jdsJA/jdsValues array
         *
         *  @param[out] row indexes of rows that have an entry for column j
         *  @param[out] pos positions of entries with col = j in csrJA, 
         *  @param[in] j is the column of which positions are required
         *  @param[in] numRows is the number of rows
         *  @param[in] ilg
         *  @param[in] dlg
         *  @param[in] perm
         *  @param[in] ja
         *  @returns  number of entries with col index = j
         */
        typedef IndexType ( *FuncType ) (
            IndexType row[],
            IndexType pos[],
            const IndexType j,
            const IndexType numRows,
            const IndexType ilg[],
            const IndexType dlg[],
            const IndexType perm[],
            const IndexType ja[] );

        static const char* getId()
        {
            return "JDS.getValuePosCol";
        }
    };

    template<typename ValueType, typename OtherValueType>
    struct scaleValue
    {
        typedef void ( *FuncType ) ( const IndexType numRows,
                                     const IndexType perm[],
                                     const IndexType ilg[],
                                     const IndexType dlg[],
                                     ValueType mValues[],
                                     const OtherValueType values[] );

        static const char* getId()
        {
            return "JDS.scaleValue";
        }
    };

    struct checkDiagonalProperty
    {
        typedef bool ( *FuncType ) ( const IndexType numDiagonals,
                                     const IndexType numRows,
                                     const IndexType numColumns,
                                     const IndexType perm[],
                                     const IndexType ja[],
                                     const IndexType dlg[] );

        static const char* getId()
        {
            return "JDS.checkDiagonalProperty";
        }
    };
};

} /* end namespace sparsekernel */

} /* end namespace scai */
