/**
 * @file MatrixCreator.hpp
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
 * @brief Class with static methods for the creation of distribted matrices.
 * @author Thomas Brandes
 * @date 01.02.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

// std
#include <memory>

namespace scai
{

namespace lama
{

/**
 * @brief This class provides some static methods to build or create sparse matrices for
 *        certain problem classes.
 */

class COMMON_DLL_IMPORTEXPORT MatrixCreator
{
public:

    /** Fill an allocated matrix in CSR format with random values
     *
     *  param[inout] matrix size and distribution of matrix remains unchanged
     *  param[in] density specifies the density of sparse entries (0.0 is empty, 1.0 is full )
     *
     *  The non-zero values of the matrix will be in the range [0,1] and might be adapted by
     *  the scale method of matrices.
     */

    static void fillRandom( _Matrix& matrix, float density );

    /** Typed version of fill random, called via meta-programming */

    template<typename ValueType>
    static void fillRandomImpl( Matrix<ValueType>& matrix, float density );

    /** Builds a block distributed matrix with random values from scratch.
     *
     *  param[out] matrix will be new defined as required
     *  param[in] size is the size of the square matrix (number of rows and number of columns)
     *  param[in] density specifies the density of sparse entries (0.0 is empty, 1.0 is full )
     */
    static void buildRandom( _Matrix& matrix, const IndexType size, const float density );

    /** Build a sparse matrix representing the discretization of the Laplacian operator
     *  on a one-dimensional structured grid.
     *
     *  @param[out] matrix in CSR format
     *  @param[in] stencilType must be 3
     *  @param[in] dim is the grid size
     */

    static void buildPoisson1D( _Matrix& matrix, const IndexType stencilType, const IndexType dim );

    /** Build a sparse matrix representing the discretization of the Laplacian operator
     *  on a two-dimensional structured grid.
     *
     *  @param[out] matrix in any format
     *  @param[in] stencilType must be 5, or 9
     *  @param[in] dim1, dim2 are the sizes of the two-dimensional grid
     */

    static void buildPoisson2D(
        _Matrix& matrix,
        const IndexType stencilType,
        const IndexType dim1,
        const IndexType dim2 );

    /** Build a sparse matrix representing the discretization of the Laplacian operator
     *  on a three-dimensional structured grid.
     *
     *  @param[out] matrix in CSR format
     *  @param[in] stencilType must be 7, 19, or 27
     *  @param[in] dim1, dim2, dim3 are the sizes of the three-dimensional grid
     */

    static void buildPoisson3D(
        _Matrix& matrix,
        const IndexType stencilType,
        const IndexType dim1,
        const IndexType dim2,
        const IndexType dim3 );

    /** Build a sparse matrix representing the discretization of the Laplacian operator
     *  on a structured grid in up to three dimensions.
     *
     *  @param[out] matrix in CSR format
     *  @param[in] dimension must be 1, 2, or 3
     *  @param[in] stencilType type of stencil, must be supported for corresponding dimension
     *  @param[in] dimX, dimY, dimZ are the sizes of the structured grid
     *
     *  Note: dimZ is ignored for dimension < 3, dimY is ignored for dimension < 2.
     */
    static void buildPoisson(
        _Matrix& matrix,
        const IndexType dimension,
        const IndexType stencilType,
        const IndexType dimX,
        const IndexType dimY,
        const IndexType dimZ );

    /** Query routine to ask for a supported stencil type
     *
     *  @param[in] dimension is dimension of the structured grid
     *  @param[in] stencilType number of points for the stencil
     */

    static bool supportedStencilType( const IndexType dimension, const IndexType stencilType );

    /** Build a (distributed) sparse matrix from a sparse storage by  replicating it in the diagonals.
     *
     *  @param[out] matrix is the generated distributed matrix
     *  @param[in] storage is the matrix storage that is replicated
     *  @param[in] nRepeat number of repitions in the diagonals
     *
     *  \code
     *     storage :   3  0  1  0
     *                 0  1  0  -1
     *
     *     buildReplicatedDiag( matrix, storage, 2 )
     *
     *     matrix  :   3  0  1  0   0  0  0  0
     *                 0  1  0  -1  0  0  0  0
     *                 0  0  0   0  3  0  1  0
     *                 0  0  0   0  0  1  0  -1
     *  \endcode
     */

    template<typename ValueType>
    static void buildReplicatedDiag( SparseMatrix<ValueType>& matrix,
                                     const MatrixStorage<ValueType>& storage,
                                     const IndexType nRepeat );

    /** Build a (distributed) sparse matrix from a (replicated) sparse storage by
     *  replicating it m times in the rows and n times in the columns
     *
     *  @param[out] matrix is the generated distributed matrix
     *  @param[in] storage is the matrix storage that is replicated
     *  @param[in] nRepeatRow number of repitions in the row
     *  @param[in] nRepeatCol number of repitions in the columns
     *
     *  \code
     *     storage :   3  0  1  0
     *                 0  1  0  -1
     *
     *     buildReplicated( matrix, storage, 3, 2 )
     *
     *     matrix  :   3  0  1  0   0  0  0  0
     *                 0  1  0  -1  0  0  0  0
     *                 3  0  1   0  3  0  1  0
     *                 0  1  0  -1  0  1  0  -1
     *                 3  0  1   0  3  0  1  0
     *                 0  1  0  -1  0  1  0  -1
     *  \endcode
     */

    template<typename ValueType>
    static void buildReplicated( SparseMatrix<ValueType>& matrix,
                                 const MatrixStorage<ValueType>& storage,
                                 const IndexType nRepeatRow,
                                 const IndexType nRepeatCol );

private:

    /** Help routine that generates a random sparsity pattern */

    static void randomCSRPattern( std::vector<IndexType>& csrIA, std::vector<IndexType>& csrJA,
                                  const IndexType numRows, const IndexType numColums, float density );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */
