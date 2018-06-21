/**
 * @file JDSKernelTrait.hpp
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
 * @brief Struct with traits for all JDS storage methods provided as kernels.
 * @author Thomas Brandes
 * @date 23.10.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/MatrixOp.hpp>

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

    /** Structure with trait for jacobi iterate on halo for JDS storage */

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

    template<typename ValueType>
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
                                     ValueType csrValues[],
                                     const IndexType csrIA[],
                                     const IndexType numRows,
                                     const IndexType jdsPerm[],
                                     const IndexType jdsILG[],
                                     const IndexType jdsDLG[],
                                     const IndexType jdsJA[],
                                     const ValueType jdsValues[] );

        static const char* getId()
        {
            return "JDS.getCSRValues";
        }
    };

    template<typename ValueType>
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
                                    ValueType jdsValues[],
                                    const IndexType numRows,
                                    const IndexType jdsPerm[],
                                    const IndexType jdsILG[],
                                    const IndexType numDiagonals,
                                    const IndexType jdsDLG[],
                                    const IndexType csrIA[],
                                    const IndexType csrJA[],
                                    const ValueType csrValues[] );

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
         *  @param numColumns is number of columns for the matrix
         *  @param numDiagonals are number of diagonals, is size of jdsDLG
         *  @param jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues are arrays of JDS storage
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
            const IndexType jdsPerm[],
            const IndexType jdsILG[],
            const IndexType ndlg,
            const IndexType jdsDLG[],
            const IndexType jdsJA[],
            const ValueType jdsValues[],
            const common::MatrixOp op );

        static const char* getId()
        {
            return "JDS.normalGEMV";
        }
    };

    template<typename ValueType>
    struct getRow
    {
        typedef void ( *FuncType ) ( ValueType row[],
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

    template<typename ValueType>
    struct setRow
    {
        typedef void ( *FuncType ) ( ValueType values[],
                                     const IndexType i,
                                     const IndexType numColumns,
                                     const IndexType numRows,
                                     const IndexType perm[],
                                     const IndexType ilg[],
                                     const IndexType dlg[],
                                     const IndexType ja[],
                                     const ValueType row[],
                                     const common::BinaryOp op );

        static const char* getId()
        {
            return "JDS.setRow";
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

    struct getDiagonalPositions
    {
        /** Get all positions for the diagonal entries of a JDS storage. 
         */
        typedef IndexType ( *FuncType ) ( 
            IndexType diagonalPositions[],
            const IndexType numDiagonals,
            const IndexType numRows,
            const IndexType jdsDLG[],
            const IndexType jdsILG[],
            const IndexType jdsPerm[],
            const IndexType jdsJA[] );

        static const char* getId()
        {
            return "JDS.getDiagonalPositions";
        }
    };

    struct getRowPositions
    {
        /** This method returns for a certain row of the JDS matrix all
         *  corresponding positions in the jdsJA/jdsValues array belonging to the row
         *
         *  @param[out] pos positions of entries with row = i 
         *  @param[in] i is the row of which positions are required
         *  @param[in] numRows is the number of rows
         *  @param[in] ilg
         *  @param[in] dlg
         *  @param[in] perm
         *  @returns  number of non-zero entries in row i
         */
        typedef IndexType ( *FuncType ) (
            IndexType pos[],
            const IndexType i,
            const IndexType numRows,
            const IndexType ilg[],
            const IndexType dlg[],
            const IndexType perm[] );

        static const char* getId()
        {
            return "JDS.getRowPositions";
        }
    };

    struct getColumnPositions
    {
        /** This method returns for a certain column of the JDS matrix all
         *  row indexes for which elements exist and the corresponding positions
         *  in the jdsJA/jdsValues array
         *
         *  @param[out] row indexes of rows that have an entry for column j
         *  @param[out] pos positions of entries with jdsJA[pos[i]] = j
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
            return "JDS.getColumnPositions";
        }
    };

    template<typename ValueType>
    struct setRows
    {
        /** This method scales each row of the matrix with a certain value
         *
         * @param[in]     numRows    is the number of rows
         * @param[in]     perm       perm[i] is the original position of row i
         * @param[in]     ilg        ilg[i] number of entries in row i
         * @param[in]     dlg        diagonals
         * @param[in,out] jdsValues  array containing all non-zero values, are scaled row-wise
         * @param[in]     rowValues  rowvalues[i] is applied to  row i
         * @param[in]     op         binary operation specifies how to update
         */
        typedef void ( *FuncType ) ( ValueType jdsValues[],
                                     const IndexType numRows,
                                     const IndexType perm[],
                                     const IndexType ilg[],
                                     const IndexType dlg[],
                                     const ValueType rowValues[],
                                     const common::BinaryOp op );

        static const char* getId()
        {
            return "JDS.setRows";
        }
    };
};

} /* end namespace sparsekernel */

} /* end namespace scai */
