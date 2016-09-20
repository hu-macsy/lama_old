/**
 * @file BLAS_BLAS3.hpp
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
 * @brief BLAS_BLAS3.hpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>
#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

/** Implementations of blas3 methods for BLASKernelTrait by using BLAS library. */

class COMMON_DLL_IMPORTEXPORT BLAS_BLAS3
{
public:

    /** OpenMP implementation for BLAS3Interface<ValueType>::gemm */

    template<typename ValueType>
    static void gemm(
        const CBLAS_ORDER order,
        const CBLAS_TRANSPOSE TransA,
        const CBLAS_TRANSPOSE TransB,
        const IndexType M,
        const IndexType N,
        const IndexType K,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* B,
        const IndexType ldb,
        const ValueType beta,
        ValueType* C,
        const IndexType ldc );

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
     */
    template<typename ValueType>
    static void symm(
        const CBLAS_ORDER order,
        const CBLAS_SIDE side,
        const CBLAS_UPLO uplo,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* B,
        const IndexType ldb,
        const ValueType beta,
        ValueType* C,
        const IndexType ldc );

    /**
     * @brief trmm performs one of the matrix-matrix operations
     *
     * B = alpha * op(A) * B or
     * B = alpha * B * op(A)
     * where op(A) = A or op(A) = AT
     *
     * alpha is a scalar,
     * B is an m×n matrix consisting of double‐precision elements,
     * A is a unit or non‐unit, upper or lower triangular matrix consisting of double‐precision elements.
     * matrix A and B are stored in column‐major format,
     * lda and ldb are the leading dimensions of the two‐dimensional arrays
     * that contain A and B, respectively.
     *
     * @param[in] order   TODO[doxy] Complete Description.
     * @param[in] side    specifies whether op(A) multiplies B from the left or right.
     *                    If side == 'L' or 'l', C = alpha * op(A) * B
     *                    If side == 'R' or 'r', C = alpha * B * op(A)
     * @param[in] uplo    specifies whether the  matrix A is an upper or lower triangular matrix.
     *                    If uplo == 'U' or 'u',
     *                    A is an upper triangular matrix.
     *                    If uplo == 'L' or 'l',
     *                    A is a lower triangular matrix.
     * @param[in] transA  specifies the form of op(A) to be used in the matrix multiplication.
     *                    If transa == 'N' or 'n', op(A) = A.
     *                    If transa == 'T','t','C','c', op(A) = AT.
     * @param[in] diag    specifies whether or not A is a unit triangular matrix.
     *                    If diag == 'U' or 'u',
     *                    A is assumed to be unit triangular.
     *                    If diag == 'N' or 'n',
     *                    A is not assumed to be unit triangular.
     * @param[in] m       the number of rows of matrix B;
     *                    m must be at least zero.
     * @param[in] n       the number of columns of matrix B;
     *                    n must be at least zero.
     * @param[in] alpha   scalar multiplier applied to op(A) * B or B * op(A),
     *                    respectively.
     *                    If alpha is zero,
     *                    no accesses are made to matrix A and no read accesses are made to matrix B.
     * @param[in] A       array of dimensions (lda, k).
     *                    If side == 'L' or 'l', k = m.
     *                    If side == 'R' or 'r', k = n.
     *                    If uplo == 'U' or 'u',
     *                    the leading k×k upper triangular part of the array A
     *                    must contain the upper triangular matrix,
     *                    and the strictly lower triangular part of A is not referenced.
     *                    If uplo == 'L' or 'l',
     *                    the leading k×k lower triangular part of the array A
     *                    must contain the lower triangular matrix,
     *                    and the strictly upper triangular part of A is not referenced.
     *                    When diag == 'U' or 'u',
     *                    the diagonal elements of A are not referenced and are assumed to be unity.
     * @param[in] lda     leading dimension of A.
     *                    When side == 'L' or 'l', it must be at least max(1, m)
     *                    and at least max(1, n) otherwise.
     * @param[in] B       array of dimensions (ldb, n).
     *                    On entry, the leading m×n part of the array contains the matrix B.
     *                    It is overwritten with the transformed matrix on exit
     * @param[in] ldb     leading dimension of B;
     *                    ldb must be at least max(1,m).
     */
    template<typename ValueType>
    static void trmm(
        const CBLAS_ORDER order,
        const CBLAS_SIDE side,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE transA,
        const CBLAS_DIAG diag,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        ValueType* B,
        const IndexType ldb );

    /**
     * @brief trsm solves one of the matrix equations
     *
     * op(A) * X = alpha * B or
     * X * op(A) = alpha * B
     * where op(A) = A or op(A) = AT
     *
     * alpha is a scalar,
     * X and B are m×n matrices that consist of elements.
     * A is a unit or non‐unit, upper or lower, triangular matrix.
     * The result matrix X overwrites input matrix B;
     * that is, on exit the result is stored in B.
     * matrix A and B are stored in column‐major format,
     * lda and ldb are the leading dimensions of the two‐dimensional arrays
     * that contain A and B, respectively.
     *
     * @param[in] order   TODO[doxy] Complete Description.
     * @param[in] side    specifies whether op(A) appears on the left or right of X:
     *                    If side == 'L' or 'l', op(A) * X = alpha * B
     *                    If side == 'R' or 'R', X * op(A) = alpha * B
     * @param[in] uplo    specifies whether the  matrix A is an upper or lower triangular matrix.
     *                    If uplo == 'U' or 'u',
     *                    A is an upper triangular matrix.
     *                    If uplo == 'L' or 'l',
     *                    A is a lower triangular matrix.
     * @param[in] transA  specifies the form of op(A) to be used in the matrix multiplication.
     *                    If transa == 'N' or 'n', op(A) = A.
     *                    If transa == 'T','t','C','c', op(A) = AT.
     * @param[in] diag    specifies whether or not A is a unit triangular matrix.
     *                    If diag == 'U' or 'u',
     *                    A is assumed to be unit triangular.
     *                    If diag == 'N' or 'n',
     *                    A is not assumed to be unit triangular.
     * @param[in] m       the number of rows of matrix B;
     *                    m must be at least zero.
     * @param[in] n       the number of columns of matrix B;
     *                    n must be at least zero.
     * @param[in] alpha   scalar multiplier applied to B.
     *                    When alpha is zero,
     *                    A is not referenced and B does not have to be a valid input.
     * @param[in] A       array of dimensions (lda, k),
     *                    where k is m when side == 'L' or 'l'
     *                    and is n when side == 'R' or 'r'.
     *                    If uplo == 'U' or 'u',
     *                    the leading k×k upper triangular part of the array A
     *                    must contain the upper triangular matrix,
     *                    and the strictly lower triangular matrix of A is not referenced.
     *                    When uplo == 'L' or 'l',
     *                    the leading k×k lower triangular part of the array A
     *                    must contain the lower triangular matrix,
     *                    and the strictly upper triangular part of A is not referenced.
     *                    If diag == 'U' or 'u',
     *                    the diagonal elements of A are not referenced and are assumed to be unity.
     * @param[in] lda     leading dimension of the two-dimensional array containing A.
     *                    When side == 'L' or 'l', lda must be at least max(1, m).
     *                    When side == 'R' or 'r', lda must be at least max(1, n).
     * @param[in] B       array of dimensions (ldb, n);
     *                    ldb must be at least max(1, m). T
     *                    he leading m×n part of the array B must contain
     *                    the righthand side matrix B.
     *                    On exit B is overwritten by the solution matrix X.
     * @param[in] ldb     leading dimension of the two-dimensional array containing B;
     *                    ldb must be at least max(1, m).
     */
    template<typename ValueType>
    static void trsm(
        const CBLAS_ORDER order,
        const CBLAS_SIDE side,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE transA,
        const CBLAS_DIAG diag,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        ValueType* B,
        const IndexType ldb );

    /**
     * @brief syrk performs one of the symmetric rank k operations
     *
     * C = alpha * A * AT + beta * C or C = alpha * AT * A + beta * C
     *
     * where alpha and beta are scalars.
     * C is an n×n symmetric matrix consisting of elements
     * and is stored in either lower or upper storage mode.
     * A is a matrix consisting of elements with dimensions of n×k
     * in the first case and k×n in the second case.
     *
     * @param[in] order   TODO[doxy] Complete Description.
     * @param[in] uplo    specifies whether the symmetric matrix C is
     *                    stored in upper or lower storage mode.
     *                    If uplo == 'U' or 'u',
     *                    only the upper triangular part of the symmetric matrix is referenced,
     *                    and the elements of the strictly lower triangular part are inferred
     *                    from those in the upper triangular part.
     *                    If uplo == 'L' or 'l',
     *                    only the lower triangular part of the symmetric matrix is referenced,
     *                    and the elements of the strictly upper triangular part are inferred
     *                    from those in the lower triangular part.
     * @param[in] trans   specifies the operation to be performed.
     *                    If trans == 'N', or 'n'
     *                    C = alpha * A * AT + beta * C.
     *                    If trans == 'T','t','C','c'
     *                    C = alpha * AT * A + beta * C.
     * @param[in] n       specifies the number of rows and the number columns of matrix C.
     *                    If trans == 'N' or 'n',
     *                    n specifies the number of rows of matrix A.
     *                    If trans == 'T', 't', 'C', or 'c',
     *                    n specifies the number of columns of matrix A; n must be at least zero.
     * @param[in] k       If trans == 'N' or 'n',
     *                    k specifies the number of columns of matrix A.
     *                    If trans == 'T', 't', 'C', or 'c',
     *                    k specifies the number of rows of matrix A;
     *                    k must be at least zero.
     * @param[in] alpha   scalar multiplier applied to A * AT or AT * A
     * @param[in] A       array of dimensions (lda, ka),
     *                    where ka is k when trans == 'N' or 'n'
     *                    and is n otherwise.
     *                    When trans == 'N' or 'n',
     *                    the leading n×k part of array A contains the matrix A;
     *                    otherwise,
     *                    the leading k×n part of the array contains the matrix A.
     * @param[in] lda     leading dimension of A.
     *                    When trans == 'N' or 'n', lda must be at least max(1, n).
     *                    Otherwise lda must be at least max(1, k).
     * @param[in] beta    scalar multiplier applied to C.
     *                    If beta is zero, C is not read.
     * @param[in] C       array of dimensions (ldc, n).
     *                    If uplo == 'U' or 'u',
     *                    the leading n×n triangular part of the array C must contain
     *                    the upper triangular part of the symmetric matrix C,
     *                    and the strictly lower triangular part of C is not referenced.
     *                    On exit, the upper triangular part of C is overwritten
     *                    by the upper triangular part of the updated matrix.
     *                    If uplo == 'L' or 'l',
     *                    the leading n×n triangular part of the array C must contain
     *                    the lower triangular part of the symmetric matrix C,
     *                    and the strictly upper triangular part of C is not referenced.
     *                    On exit, the lower triangular part of C is overwritten
     *                    by the lower triangular part of the updated matrix.
     * @param[in] ldc     leading dimension of C.
     *                    ldc must be at least max(1,n)
     */
    template<typename ValueType>
    static void syrk(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const IndexType n,
        const IndexType k,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType beta,
        ValueType* C,
        const IndexType ldc );

    /**
     * @brief syrk2 performs one of the symmetric rank 2k operations
     *
     * C = alpha * A * BT + alpha * B * AT + beta * C or
     * C = alpha * AT * B + alpha * BT * A + beta * C
     *
     * where alpha and beta are sscalars.
     * C is an n×n symmetric matrix consisting of elements
     * and is stored in either lower or upper storage mode.
     * A and B are matrices consisting of elements with dimension
     * of n×k in the first case and k×n in the second case.
     *
     * @param[in] order   TODO[doxy] Complete Description.
     * @param[in] uplo    specifies whether the symmetric matrix C is
     *                    stored in upper or lower storage mode.
     *                    If uplo == 'U' or 'u',
     *                    only the upper triangular part of the symmetric matrix is referenced,
     *                    and the elements of the strictly lower triangular part are inferred
     *                    from those in the upper triangular part.
     *                    If uplo == 'L' or 'l',
     *                    only the lower triangular part of the symmetric matrix is referenced,
     *                    and the elements of the strictly upper triangular part are inferred
     *                    from those in the lower triangular part.
     * @param[in] trans   specifies the operation to be performed.
     *                    If trans = 'N' or 'n'
     *                    C = alpha * A * BT + alpha * B * AT + beta * C
     *                    If trans == 'T','t','C','c'
     *                    alpha * AT * B + alpha * BT * A + beta * C
     * @param[in] n       specifies the number of rows and the number columns of matrix C.
     *                    If trans == 'N', or 'n'
     *                    n specifies the number of rows of matrix A.
     *                    If trans == 'T','t','C','c'
     *                    n specifies the number of columns of matrix A;
     *                    n must be at least zero.
     * @param[in] k       If trans == 'N', or 'n'
     *                    k specifies the number of columns of matrix A.
     *                    If trans == 'T','t','C','c'
     *                    k specifies the number of rows of matrix A;
     *                    k must be at least zero.
     * @param[in] alpha   scalar multiplier applied to A * AT or AT * A
     * @param[in] A       array of dimensions (lda, ka),
     *                    where ka is k when trans == 'N' or 'n'
     *                    and is n otherwise.
     *                    When trans == 'N' or 'n',
     *                    the leading n×k part of array A must contain the matrix A,
     *                    otherwise,
     *                    the leading k×n part of the array must contain the matrix A.
     * @param[in] lda     leading dimension of A.
     *                    When trans == 'N' or 'n', lda must be at least max(1, n).
     *                    Otherwise lda must be at least max(1, k).
     * @param[in] B       TODO[doxy] Complete Description.
     * @param[in] ldb     TODO[doxy] Complete Description.
     * @param[in] beta    scalar multiplier applied to C.
     *                    If beta is zero, C does not have to be a valid input.
     * @param[in] C       array of dimensions (ldc, n).
     *                    If uplo == 'U' or 'u',
     *                    the leading n×n triangular part of the array C must contain
     *                    the upper triangular part of the symmetric matrix C,
     *                    and the strictly lower triangular part of C is not referenced.
     *                    On exit, the upper triangular part of C is overwritten
     *                    by the upper triangular part of the updated matrix.
     *                    If uplo == 'L' or 'l',
     *                    the leading n×n triangular part of the array C must contain
     *                    the lower triangular part of the symmetric matrix C,
     *                    and the strictly upper triangular part of C is not referenced.
     *                    On exit, the lower triangular part of C is overwritten
     *                    by the lower triangular part of the updated matrix.
     * @param[in] ldc     leading dimension of C.
     *                    ldc must be at least max(1,n)
     */
    template<typename ValueType>
    static void syrk2(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const IndexType n,
        const IndexType k,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* B,
        const IndexType ldb,
        const ValueType beta,
        ValueType* C,
        const IndexType ldc );

private:

    /** Routine that registers all methods at the kernel registry. */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    BLAS_BLAS3();

    /** Destructor for unregistration. */

    ~BLAS_BLAS3();

    /** Static variable for registration at static initialization. */

    static BLAS_BLAS3 guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace blaskernel */

} /* end namespace scai */
