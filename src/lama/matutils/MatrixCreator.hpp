/**
 * @file MatrixCreator.hpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Class with static methods for the creation of distribted matrices.
 * @author Thomas Brandes
 * @date 01.02.2012
 * $Id$
 */

#ifndef LAMA_MATRIX_CREATOR_HPP_
#define LAMA_MATRIX_CREATOR_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/matrix/CSRSparseMatrix.hpp>

#include <memory>

/** This class creates 'distributed' marices for poisson solvers.
 *
 *  \code
 *  std::auto_ptr<CSRSparseMatrix<double> > MatrixCreator<double>::createPoisson2D( 9, 4, 4 )  );
 *  \endcode
 *
 *  The matrix A will have a general distribution and is in CSR format.
 */

namespace lama
{

/**
 * @brief This class provides some static methods to build or create sparse matrices for
 *        certain problem classes.
 *
 * @tparam T is the value type of the matrix values.
 */

template<typename T>
class LAMA_DLL_IMPORTEXPORT MatrixCreator
{
public:

    /** Fill a matrix in CSR format with random values
     *
     *  param[inout] matrix size and distribution of matrix remains unchanged
     *  param[in] density specifies the density of sparse entries (0.0 is empty, 1.0 is full )
     *
     */

    static void fillRandom( Matrix& matrix, double density );

    /** Builds a block distributed CSR matrix with random values.
     *
     *  param[out] matrix will be new defined as required
     *  param[in] size is the size of the square matrix (number of rows and number of columns)
     *  param[in] density specifies the density of sparse entries (0.0 is empty, 1.0 is full )
     */

    static void buildRandom( CSRSparseMatrix<T>& matrix, const IndexType size, const double density );

    /** Build a sparse matrix representing the discretization of the Laplacian operator
     *  on a one-dimensional structured grid.
     *
     *  @param[out] matrix in CSR format
     *  @param[in] stencilType must be 3
     *  @param[in] dim is the grid size
     */

    static void buildPoisson1D( CSRSparseMatrix<T>& matrix, const IndexType stencilType, const IndexType dim );

    /** Build a sparse matrix representing the discretization of the Laplacian operator
     *  on a two-dimensional structured grid.
     *
     *  @param[out] matrix in CSR format
     *  @param[in] stencilType must be 5, or 9
     *  @param[in] dim1, dim2 are the sizes of the two-dimensional grid
     */

    static void buildPoisson2D(
        CSRSparseMatrix<T>& matrix,
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
        CSRSparseMatrix<T>& matrix,
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
        CSRSparseMatrix<T>& matrix,
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

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace lama

#endif  // LAMA_MATRIX_CREATOR_HPP_
