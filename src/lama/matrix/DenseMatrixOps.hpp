/**
 * @file DenseMatrixOps.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Definition of matrix operations for dense matrices.
 * @author Thomas Brandes
 * @date 04.01.2012
 * $Id$
 */
#ifndef LAMA_DENSE_MATRIX_OPS_HPP_
#define LAMA_DENSE_MATRIX_OPS_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/matrix/DenseMatrix.hpp>

namespace lama
{

/** Static class that offers some matrix operations for dense
 *  matrices.
 *
 *  @todo invert as method for DenseStorage, MatrixStorage
 *  @todo invert as method for DenseMatrix, Matrix
 */

class LAMA_DLL_IMPORTEXPORT DenseMatrixOps
{
public:

    /** Compute inverse of a square replicated matrix in-place. */

    template<typename T>
    static void invertReplicated( DenseMatrix<T>& matrix );

    /** Compute inverse of a square block-cyclic distributed matrix in-place. */

    template<typename T>
    static void invertCyclic( DenseMatrix<T>& matrix );

    /** Compute inverse of a square distributed or replicated matrix in-place. */

    template<typename T>
    static void invert( DenseMatrix<T>& matrix );

private:

    LAMA_LOG_DECL_STATIC_LOGGER(logger);

    static int iceil(int inum, int idenom)
    {
        return (inum+idenom-1)/idenom;
    };
};

} // namespace lama

#endif // LAMA_DENSE_MATRIX_OPS_HPP_
