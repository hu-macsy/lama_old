/**
 * @file TestSparseMatrices.hpp
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
 * @brief TestSparseMatrices.hpp
 * @author Jiri Kraus
 * @date 06.04.2011
 */

#pragma once

#include <memory>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/common/unique_ptr.hpp>

/**
 * @brief The class TestSparseMatrices summaries test sparse test matrices that
 *        are needed by the unit tests.
 *
 * The class TestSparseMatrices summaries test sparse test matrices that are
 * needed by the unit tests. TestSparseMatrices consists only of static methods
 * that return copies to Sparsematrix. Some of these matrix are related to
 * each other. E.g. n4m4InverseTestMatrix1 is the inverse of n4m4TestMatrix1.
 */
class TestSparseMatrices
{
public:

    /**
     * @brief Creates a 4x4 matrix (TestMatrix1).
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4TestMatrix1();

    /**
     * @brief Creates the inverse of TestMatrix1.
     *
     * @see n4m4TestMatrix1
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4InverseTestMatrix1();

    /**
     * @brief Creates a diagonal 4x4 matrix.
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4DiagonalMatrix();

    /**
     * @brief Creates a symmetric 4x4 matrix.
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4SymmetricMatrix();

    /**
     * @brief Creates 4x6 matrix.
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m6NoneSquareMatrix();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixARes = n4m4MatrixA1 * n4m4MatrixA2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4MatrixA1();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixARes = n4m4MatrixA1 * n4m4MatrixA2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4MatrixA2();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixARes = n4m4MatrixA1 * n4m4MatrixA2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4MatrixARes();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixBRes = n4m4MatrixB1 * n4m4MatrixB2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4MatrixB1();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixBRes = n4m4MatrixB1 * n4m4MatrixB2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4MatrixB2();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixBRes = n4m4MatrixB1 * n4m4MatrixB2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4MatrixBRes();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixCRes = n4m4MatrixC1 * n4m4MatrixC1
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4MatrixC1();
    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixCRes = n4m4MatrixC1 * n4m4MatrixC1
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4MatrixCRes();

    /**
     * @brief Creates a 6x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m6MatrixDRes = n6m4MatrixD1 * n4m6MatrixD2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n6m4MatrixD1();

    /**
     * @brief Creates a 4x6 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m6MatrixDRes = n6m4MatrixD1 * n4m6MatrixD2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m6MatrixD2();

    /**
     * @brief Creates a 6x6 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m6MatrixDRes = n6m4MatrixD1 * n4m6MatrixD2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n6m6MatrixDRes();

    /**
     * @brief Creates a 6x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m3MatrixERes = n6m4MatrixE1 * n4m3MatrixE2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n6m4MatrixE1();

    /**
     * @brief Creates a 4x3 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m3MatrixERes = n6m4MatrixE1 * n4m3MatrixE2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m3MatrixE2();

    /**
     * @brief Creates a 6x3 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m3MatrixERes = n6m4MatrixE1 * n4m3MatrixE2
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n6m3MatrixERes();

    /**
     * @brief Creates a 8x4 matrix.
     *
     * Creates a 8x4 matrix which is an direct interpolation for the matrix
     * returned by n8m8Laplace1D().
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n8m4Interpol();

    /**
     * @brief Creates a 4x8 matrix.
     *
     * Creates a 4x8 matrix which is an transpose of the matrix returned by
     * n8m4Interpol().
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m8InterpolTranspose();

    /**
     * @brief Creates a 8x8 matrix.
     *
     * Creates a 8x8 matrix which is the discretization of a 1D laplace equation.
     *
     * @return the matrix
     */
    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n8m8Laplace1D();

    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n8m4GalerkinTemp();

    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4Galerkin();

    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> nnIdentityMatrix( IndexType n );

    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n4m4IdentityMatrix();

    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n6m6TestMatrix();


    template<typename ValueType>
    static scai::lama::CSRSparseMatrix<ValueType> n6m6Full();
//    template<typename ValueType>
//    static scai::lama::CSRSparseMatrix<ValueType> n8m4MatrixF1();
//
//    template<typename ValueType>
//    static scai::lama::CSRSparseMatrix<ValueType> n4m6MatrixF2();
//
//    template<typename ValueType>
//    static scai::lama::CSRSparseMatrix<ValueType> n8m6MatrixFRes();
//
};

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4TestMatrix1()
{
    const IndexType n = 4;
    ValueType randomValues[] =
    {
        0.436213f, 0.683202f, 0.531013f, 0.422152f, 0.4632f, 0.168648f, 0.967549f, 0.498486f, 0.126115f, 0.708545f,
        0.131853f, 0.820422f, 0.992481f, 0.202542f, 0.47369f, 0.947076f
    };
    scai::lama::CSRSparseMatrix<ValueType> randomMatrix;
    randomMatrix.setRawDenseData( n, n, randomValues );
    return randomMatrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4InverseTestMatrix1()
{
    const IndexType n = 4;
    ValueType randomValuesInverse[] =
    {
        1.26932f, -1.06467f, -1.28787f, 1.11023f, 1.72369f, -0.749837f, 0.0459759f, -0.41348f, -0.0443016f, 1.36936f,
        -0.110096f, -0.605633f, -1.67664f, 0.591171f, 1.39484f, 0.283764f
    };
    scai::lama::CSRSparseMatrix<ValueType> randomMatrixInverse;
    randomMatrixInverse.setRawDenseData( n, n, randomValuesInverse );
    return randomMatrixInverse;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4IdentityMatrix()
{
    const IndexType n = 4;
    ValueType identityValues[] =
    { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };
    scai::lama::CSRSparseMatrix<ValueType> n4m4IdentM;
    n4m4IdentM.setRawDenseData( n, n, identityValues );
    return n4m4IdentM;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::nnIdentityMatrix( IndexType n )
{
    scai::common::scoped_array<ValueType> identityValues( new ValueType[n * n] );

    for ( int i = 1; i <= n * n; ++i )
    {
        identityValues[i - 1] = 0.0;
    }

    for ( int i = 0; i < n; ++i )
    {
        identityValues[i * n + i] = 1.0;
    }

    scai::lama::CSRSparseMatrix<ValueType> nnIdentM;
    nnIdentM.setRawDenseData( n, n, identityValues.get() );
    return nnIdentM;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4DiagonalMatrix()
{
    const IndexType n = 4;
    ValueType values[] =
    { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.0f };
    scai::lama::CSRSparseMatrix<ValueType> m;
    m.setRawDenseData( n, n, values );
    return m;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n6m6TestMatrix()
{
    const IndexType n = 6;
    ValueType values[] =
    {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    scai::lama::CSRSparseMatrix<ValueType> sparseM;
    sparseM.setRawDenseData( n, n, values );
    return sparseM;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4SymmetricMatrix()
{
    const IndexType n = 4;
    ValueType values[] =
    { 1.0f, 0.0f, 0.0f, 5.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f, 5.0f, 0.0f, 0.0f, 4.0f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, n, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m6NoneSquareMatrix()
{
    const IndexType n = 4;
    const IndexType m = 6;
    ValueType values[] =
    {
        0.436213f, 0.0f, 0.0f, 0.422152f, 0.4632f, 0.168648f, 0.0f, 0.0f, 0.0f, 0.708545f, 0.131853f, 0.820422f, 0.1f,
        0.0f, 0.131853f, 0.820422f, 0.0f, 0.1f, 0.0f, 0.820422f, 0.0f, 0.202542f, 0.0f, 0.947076f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4MatrixA1()
{
    const IndexType n = 4;
    ValueType values[] =
    { 0.6f, 0.0f, 0.0f, 0.4f, 0.7f, 0.4f, 0.0f, 0.0f, 0.0f, 0.0f, 0.9f, 0.4f, 0.2f, 0.5f, 0.0f, 0.3f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, n, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4MatrixA2()
{
    const IndexType n = 4;
    ValueType values[] =
    { 0.1f, 0.4f, 0.0f, 0.0f, 0.0f, 0.1f, 0.0f, 0.0f, 0.9f, 0.0f, 0.9f, 0.0f, 0.0f, 0.0f, 0.2f, 0.1f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, n, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4MatrixARes()
{
    const IndexType n = 4;
    ValueType values[] =
    { 0.06f, 0.24f, 0.08f, 0.04f, 0.07f, 0.32f, 0.0f, 0.0f, 0.81f, 0.0f, 0.89f, 0.04f, 0.02f, 0.13f, 0.06f, 0.03f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, n, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4MatrixB1()
{
    const IndexType n = 4;
    ValueType values[] =
    { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.0f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, n, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4MatrixB2()
{
    const IndexType n = 4;
    ValueType values[] =
    { 4.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, n, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4MatrixBRes()
{
    const IndexType n = 4;
    ValueType values[] =
    { 4.0f, 0.0f, 0.0f, 0.0f, 0.0f, 6.0f, 0.0f, 0.0f, 0.0f, 0.0f, 6.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.0f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, n, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4MatrixC1()
{
    const IndexType n = 4;
    ValueType values[] =
    { 2.0f, -1.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, -1.0f, 2.0f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, n, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4MatrixCRes()
{
    const IndexType n = 4;
    ValueType values[] =
    { 5.0f, -4.0f, 1.0f, 0.0f, -4.0f, 6.0f, -4.0f, 1.0f, 1.0f, -4.0f, 6.0f, -4.0f, 0.0f, 1.0f, -4.0f, 5.0f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, n, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n6m4MatrixD1()
{
    const IndexType n = 6;
    const IndexType m = 4;
    ValueType values[] =
    {
        6.0f, 0.0f, 0.0f, 4.0f, 7.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f, 9.0f, 4.0f, 2.0f, 5.0f, 0.0f, 3.0f, 2.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f, 2.0f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m6MatrixD2()
{
    const IndexType n = 4;
    const IndexType m = 6;
    ValueType values[] =
    {
        1.0f, 4.0f, 0.0f, 0.0f, 4.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 9.0f, 1.0f, 9.0f, 0.0f, 9.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 2.0f, 1.0f, 0.0f, 3.0f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n6m6MatrixDRes()
{
    const IndexType n = 6;
    const IndexType m = 6;
    ValueType values[] =
    {
        6.0f, 24.0f, 8.0f, 4.0f, 24.0f, 12.0f, 7.0f, 32.0f, 0.0f, 0.0f, 64.0f, 4.0f, 81.0f, 0.0f, 89.0f, 4.0f, 0.0f,
        12.0f, 2.0f, 13.0f, 6.0f, 3.0f, 53.0f, 14.0f, 2.0f, 8.0f, 2.0f, 1.0f, 8.0f, 3.0f, 0.0f, 1.0f, 4.0f, 2.0f, 9.0f,
        7.0f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n6m4MatrixE1()
{
    const IndexType n = 6;
    const IndexType m = 4;
    ValueType values[] =
    {
        6.0f, 0.0f, 0.0f, 4.0f, 7.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f, 9.0f, 4.0f, 2.0f, 5.0f, 0.0f, 3.0f, 2.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f, 2.0f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m3MatrixE2()
{
    const IndexType n = 4;
    const IndexType m = 3;
    ValueType values[] =
    { 1.0f, 4.0f, 0.0f, 0.0f, 1.0f, 0.0f, 9.0f, 0.0f, 9.0f, 0.0f, 0.0f, 2.0f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n6m3MatrixERes()
{
    IndexType n = 6;
    IndexType m = 3;
    ValueType values[] =
    { 6.0f, 24.0f, 8.0f, 7.0f, 32.0f, 0.0f, 81.0f, 0.0f, 89.0f, 2.0f, 13.0f, 6.0f, 2.0f, 8.0f, 2.0f, 0.0f, 1.0f, 4.0f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n8m4Interpol()
{
    IndexType n = 8;
    IndexType m = 4;
    ValueType values[] =
    {
        0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.5f, 0.5f,
        0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m8InterpolTranspose()
{
    IndexType n = 4;
    IndexType m = 8;
    ValueType values[] =
    {
        0.5f, 1.0f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.5f, 1.0f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.5f, 1.0f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.5f, 1.0f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n8m8Laplace1D()
{
    IndexType n = 8;
    IndexType m = 8;
    ValueType values[] =
    {
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n8m4GalerkinTemp()
{
    IndexType n = 8;
    IndexType m = 4;
    ValueType values[] =
    {
        0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -0.5f, 1.0f, -0.5f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, -0.5f, 1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -0.5f, 1.5f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n4m4Galerkin()
{
    IndexType n = 4;
    IndexType m = 4;
    ValueType values[] =
    { 1.0f, -0.5f, 0.0f, 0.0f, -0.5f, 1.0f, -0.5f, 0.0f, 0.0f, -0.5f, 1.0f, -0.5f, 0.0f, 0.0f, -0.5f, 1.5f };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}

template<typename ValueType>
scai::lama::CSRSparseMatrix<ValueType> TestSparseMatrices::n6m6Full()
{
    IndexType n = 6;
    IndexType m = 6;
    ValueType values[] =
    {
        1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 16.0f, 17.0f,
        18.0f, 19.0f, 20.0f, 21.0f, 22.0f, 23.0f, 24.0f, 25.0f, 26.0f, 27.0f, 28.0f, 29.0f, 30.0f, 31.0f, 32.0f, 33.0f, 34.0f,
        35.0f, 36.0f
    };
    scai::lama::CSRSparseMatrix<ValueType> matrix;
    matrix.setRawDenseData( n, m, values );
    return matrix;
}