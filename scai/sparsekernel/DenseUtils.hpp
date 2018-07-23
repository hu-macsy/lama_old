/**
 * @file DenseUtils.hpp
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
 * @brief Utility functions for dense storage data
 * @author Thomas Brandes
 * @date 28.05.2018
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/MatrixOp.hpp>

namespace scai
{

namespace sparsekernel
{

class COMMON_DLL_IMPORTEXPORT DenseUtils
{
public:

    /** 
     *  @brief Conversion of dense matrix storage data to the CSR storage format.
     *
     *  \param[in] denseValues is an array of size numRows * numColumns, values are stored row-wise
     *  \param[in] numRows, numColumns are the dimension of the storage
     *  \param[out] csrIA will be the offset array for the row sizes
     *  \param[out] csrJA contains the columns of the nonz-zero entries
     *  \param[out] csrValues contains the non-zero values
     *  \param[in] prefLoc is the context where conversion should be executed.
     */
    template<typename ValueType>
    static void convertDense2CSR( 
        hmemo::HArray<IndexType>& csrIA, 
        hmemo::HArray<IndexType>& csrJA, 
        hmemo::HArray<ValueType>& csrValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<ValueType>& denseValues,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief Determine for each row the number of non-zero entries
     *
     *  \param[in] numRows, numColumns are the dimension of the storage
     *  \param[in] denseValues is an array of size numRows * numColumns, values are stored row-wise
     *  \param[out] rowSizes will contain for each row the number of non-zero elements
     *  \param[in] prefLoc is the context where conversion should be executed.
     */
    template<typename ValueType>
    static void getSparseRowSizes( 
        hmemo::HArray<IndexType>& rowSizes,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<ValueType>& denseValues,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief Convert CSR storage data to a dense storage
     *
     *  \param[out] denseValues is an array of size numRows * numColumns, values are stored row-wise
     *  \param[in] numRows, numColumns are the dimensions of the storage
     *  \param[in] csrIA, csrJA, csrValues are the CSR array to be converted
     *  \param[in] prefLoc is the context where conversion should be executed.
     */
    template<typename ValueType>
    static void convertCSR2Dense(
        hmemo::HArray<ValueType>& denseValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Jacobi iteration step with dense storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     */
    template<typename ValueType>
    static void jacobi(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const IndexType n,
        const hmemo::HArray<ValueType>& denseValues,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Jacobi halo iteration step with dense storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     */
    template<typename ValueType>
    static void jacobiHalo(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& diagonal,
        const hmemo::HArray<ValueType>& oldSolution,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<ValueType>& denseValues,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief Invert a dense storage
     * 
     *  Note: throws exception if matrix cannot be inverted
     */
    template<typename ValueType>
    static void invert( 
        hmemo::HArray<ValueType>& a,
        const IndexType n,
        hmemo::ContextPtr prefLoc );

    /**
     * @brief Set/update each element in a 2D dense array
     */
    template<typename ValueType>
    static void setScalar( 
        hmemo::HArray<ValueType>& denseValues,
        const IndexType numRows,
        const IndexType numColumns,
        const ValueType scalar,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc );

    /**
     * @brief Set/update each row in a 2D dense array with an individual value
     *
     *  \code
     *      for all i = 0, ..., n-1; j = 0, ..., m-1
     *      denseValues( i, j ) = denseValues( i, j ) <op> rowValues( i )
     *  \endcode
     */
    template<typename ValueType>
    static void setRows( 
        hmemo::HArray<ValueType>& denseValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<ValueType>& rowValues,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc );

    /**
     * @brief Set/update each column in a 2D dense array with an individual value
     *
     *  \code
     *      for all i = 0, ..., n-1; j = 0, ..., m-1
     *      denseValues( i, j ) = denseValues( i, j ) <op> columnValues( j )
     *  \endcode
     */
    template<typename ValueType>
    static void setColumns(
        hmemo::HArray<ValueType>& denseValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<ValueType>& columnValues,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc );
   
    /**
     * @brief Count the number of non-zero elements in a 2D dense array.
     */
    template<typename ValueType>
    static IndexType getNumValues(
        const hmemo::HArray<ValueType>& denseValues,
        const IndexType numRows,
        const IndexType numColumns,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Dense matrix-vector multiplication: result = alpha * denseMatrix * x + beta * y
     */
    template<typename ValueType>
    static tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<ValueType>& denseValues,
        const common::MatrixOp op,
        bool async,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Dense matrix-vector multiplication: result = alpha * denseMatrix * x
     */
    template<typename ValueType>
    static tasking::SyncToken* gemv0(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<ValueType>& denseValues,
        const common::MatrixOp op,
        bool async,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Dense matrix multiplication: c = alpha * a * b + beta * c
     *
     *  @param[in,out] c is matrix storage of size m x n 
     *  @param[in] a is matrix storage of size m x k
     *  @param[in] b is matrix storage of size k x n
     *  @param[in] alpha scaling factor for matrix product
     *  @param[in] opB matrix operation for b (normal, transpose)
     *  @param[in] opA matrix operation for a (normal, transpose)
     *  @param[in] beta scaling factor for matrix c
     *  @param[in] m, n, k stand for the sizes of the matrices
     *  @param[in] prefLoc specifies the context where the operation should be done
     */
    template<typename ValueType>
    static void gemm(
        hmemo::HArray<ValueType>& c,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& a,
        const common::MatrixOp opA,
        const hmemo::HArray<ValueType>& b,
        const common::MatrixOp opB,
        const ValueType beta,
        const IndexType m,
        const IndexType n,
        const IndexType k,
        hmemo::ContextPtr prefLoc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
