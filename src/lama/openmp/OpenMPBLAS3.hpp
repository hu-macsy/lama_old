/**
 * @file blas/OpenMPBLAS3.hpp
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
 * @brief OpenMPBLAS3.hpp
 * @author Eric Schricker
 * @date 05.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_OPENMPBLAS3_HPP_
#define LAMA_OPENMPBLAS3_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/SyncToken.hpp>

#include <lama/openmp/BLASHelper.hpp>

namespace lama
{

/** Implementations of methods for lama::BLAS3Interface with OpenMP.
 *
 *  @todo Move all method documentations to LAMAInterface and make references here
 *  @todo Add information here about use of native BLAS1 libraries
 */

class LAMA_DLL_IMPORTEXPORT OpenMPBLAS3
{
public:

    /** OpenMP implementation for BLAS3Interface<T>::gemm */

    template<typename T>
    static void gemm(
        const CBLAS_ORDER order,
        const CBLAS_TRANSPOSE TransA,
        const CBLAS_TRANSPOSE TransB,
        const IndexType M,
        const IndexType N,
        const IndexType K,
        const T alpha,
        const T* A,
        const IndexType lda,
        const T* B,
        const IndexType ldb,
        const T beta,
        T* C,
        const IndexType ldc,
        SyncToken* syncToken );

    /**
     * @brief symm performs one of the matrix-matrix operations
     *
     * C = alpha * A * B + beta * C or C = alpha * B * A + beta * C,
     *
     * where alpha and beta are scalars,
     * A is a symmetric matrix consisting of single‐precision elements
     * and is stored in either lower or upper storage mode.
     * B and C are m×n matrices consisting of elements.
     *
     * @param[in] order   TODO[doxy] Complete Description.
     * @param[in] side    specifies whether the symmetric matrix A appears on the
     *                    left-hand side or right-hand side of Matrix B.
     *                    If side == 'L' or 'l', C = alpha * A * B + beta * C.
     *                    If side == 'R' or 'R', C = alpha * B * A + beta * C.
     * @param[in] uplo    specifies whether the symmetric matrix A is stored
     *                    in upper or lower storage mode.
     *                    If uplo == 'U' or 'u',
     *                    only the upper triangular part of the symmetric matrix is referenced,
     *                    and the elements of the strictly lower triangular part are inferred
     *                    from those in the upper triangular part.
     *                    If uplo == 'L' or 'l',
     *                    only the lower triangular part of the symmetric matrix is referenced,
     *                    and the elements of the strictly upper triangular part are inferred
     *                    from those in the lower triangular part.
     * @param[in] m       specifies the number of rows of matrix C,
     *                    and the number of rows of matrix B.
     *                    It also specifies the dimensions of symmetric matrix A when side == 'L' or 'l';
     *                    m must be at least zero.
     * @param[in] n       specifies the number of columns of matrix C,
     *                    and the number of columns of matrix B.
     *                    It also specifies the dimensions of symmetric matrix A when side == 'R' or 'r'';
     *                    n must be at least zero.
     * @param[in] alpha   scalar multiplier applied to A * B or B * A
     * @param[in] A       array of dimensions (lda, ka),
     *                    where ka is m when side == 'L' or 'l' and is n otherwise.
     *                    If side == 'L' or 'l',
     *                    the leading m×m part of array A must
     *                    contain the symmetric matrix such that when uplo == 'U' or 'u',
     *                    the leading m×m part stores the upper triangular part of the symmetric matrix,
     *                    and the strictly lower triangular part of A is not referenced;
     *                    and when uplo == 'L' or 'l',
     *                    the leading m×m part stores the lower triangular part of the symmetric matrix,
     *                    and the strictly upper triangular part is not referenced.
     *                    If side == 'R' or 'r',
     *                    the leading n×n part of array A must contain the symmetric matrix such that
     *                    when uplo == 'U' or 'u', the leading n×n part stores the
     *                    upper triangular part of the symmetric matrix,
     *                    and the strictly lower triangular part of A is not referenced;
     *                    and when uplo == 'L' or 'l',
     *                    the leading n×n part stores the lower triangular
     *                    part of the symmetric matrix, and the strictly upper triangular part is
     *                    not referenced.
     * @param[in] lda     leading dimension of A.
     *                    When side == 'L' or 'l', it must be at least max(1, m)
     *                    and at least max(1, n) otherwise.
     * @param[in] B       array of dimensions (ldb, n).
     *                    On entry, the leading m×n part of the array contains the matrix B.
     * @param[in] ldb     leading dimension of B;
     *                    ldb must be at least max(1,m).
     * @param[in] beta    scalar multiplier applied to C.
     *                    If beta is zero, C does not have to be a valid input
     * @param[in] C       array of dimensions (ldc,n);
     * @param[in] ldc     leading dimension of C.
     *                    ldc must be at least max(1,m)
     * TODO[doxy] Is the following description correct?
     * @param[out] syncToken  updated according to C = alpha * A * B + beta * C
     *                                          or C = alpha * B * A + beta * C
     */

    /** Routine that sets functions pointers belonging to BLAS1 in a BLASInterface.
     *
     *  @param[inout] BLASInterface struct to register all routines implemented in CUDA
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void setInterface( struct BLASInterface& BLAS );

private:

    static bool initialized;

    static bool registerInterface();

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} /* namespace lama */

#endif // LAMA_OPENMPBLAS3_HPP_
