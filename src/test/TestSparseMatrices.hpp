/**
 * @file TestSparseMatrices.hpp
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
 * @brief TestSparseMatrices.hpp
 * @author Jiri Kraus
 * @date 06.04.2011
 * @since 1.0.0
 */

#ifndef LAMA_TEST_SPARSE_MATRICES_HPP_
#define LAMA_TEST_SPARSE_MATRICES_HPP_

#include <memory>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <boost/scoped_array.hpp>

using namespace lama;

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
    template<typename T>
    static CSRSparseMatrix<T> n4m4TestMatrix1();

    /**
     * @brief Creates the inverse of TestMatrix1.
     *
     * @see n4m4TestMatrix1
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4InverseTestMatrix1();

    /**
     * @brief Creates a diagonal 4x4 matrix.
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4DiagonalMatrix();

    /**
     * @brief Creates a symmetric 4x4 matrix.
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4SymmetricMatrix();

    /**
     * @brief Creates 4x6 matrix.
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m6NoneSquareMatrix();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixARes = n4m4MatrixA1 * n4m4MatrixA2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4MatrixA1();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixARes = n4m4MatrixA1 * n4m4MatrixA2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4MatrixA2();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixARes = n4m4MatrixA1 * n4m4MatrixA2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4MatrixARes();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixBRes = n4m4MatrixB1 * n4m4MatrixB2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4MatrixB1();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixBRes = n4m4MatrixB1 * n4m4MatrixB2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4MatrixB2();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixBRes = n4m4MatrixB1 * n4m4MatrixB2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4MatrixBRes();

    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixCRes = n4m4MatrixC1 * n4m4MatrixC1
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4MatrixC1();
    /**
     * @brief Creates a 4x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n4m4MatrixCRes = n4m4MatrixC1 * n4m4MatrixC1
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m4MatrixCRes();

    /**
     * @brief Creates a 6x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m6MatrixDRes = n6m4MatrixD1 * n4m6MatrixD2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n6m4MatrixD1();

    /**
     * @brief Creates a 4x6 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m6MatrixDRes = n6m4MatrixD1 * n4m6MatrixD2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m6MatrixD2();

    /**
     * @brief Creates a 6x6 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m6MatrixDRes = n6m4MatrixD1 * n4m6MatrixD2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n6m6MatrixDRes();

    /**
     * @brief Creates a 6x4 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m3MatrixERes = n6m4MatrixE1 * n4m3MatrixE2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n6m4MatrixE1();

    /**
     * @brief Creates a 4x3 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m3MatrixERes = n6m4MatrixE1 * n4m3MatrixE2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m3MatrixE2();

    /**
     * @brief Creates a 6x3 matrix.
     *
     * This matrix is used for matrix multiplication tests.
     * The following matrices are needed to perform the test:
     * n6m3MatrixERes = n6m4MatrixE1 * n4m3MatrixE2
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n6m3MatrixERes();

    /**
     * @brief Creates a 8x4 matrix.
     *
     * Creates a 8x4 matrix which is an direct interpolation for the matrix
     * returned by n8m8Laplace1D().
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n8m4Interpol();

    /**
     * @brief Creates a 4x8 matrix.
     *
     * Creates a 4x8 matrix which is an transpose of the matrix returned by
     * n8m4Interpol().
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n4m8InterpolTranspose();

    /**
     * @brief Creates a 8x8 matrix.
     *
     * Creates a 8x8 matrix which is the discretization of a 1D laplace equation.
     *
     * @return the matrix
     */
    template<typename T>
    static CSRSparseMatrix<T> n8m8Laplace1D();

    template<typename T>
    static CSRSparseMatrix<T> n8m4GalerkinTemp();

    template<typename T>
    static CSRSparseMatrix<T> n4m4Galerkin();

    template<typename T>
    static CSRSparseMatrix<T> nnIdentityMatrix( IndexType n );

    template<typename T>
    static CSRSparseMatrix<T> n4m4IdentityMatrix();

    template<typename T>
    static CSRSparseMatrix<T> n6m6TestMatrix();

//    template<typename T>
//    static CSRSparseMatrix<T> n8m4MatrixF1();
//
//    template<typename T>
//    static CSRSparseMatrix<T> n4m6MatrixF2();
//
//    template<typename T>
//    static CSRSparseMatrix<T> n8m6MatrixFRes();
//
};

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4TestMatrix1()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType randomValues[] =
    {   0.436213f, 0.683202f, 0.531013f, 0.422152f, 0.4632f, 0.168648f, 0.967549f, 0.498486f, 0.126115f, 0.708545f,
        0.131853f, 0.820422f, 0.992481f, 0.202542f, 0.47369f, 0.947076f
    };

    CSRSparseMatrix<T> randomMatrix;
    randomMatrix.setRawDenseData( n, n, randomValues );

    return randomMatrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4InverseTestMatrix1()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType randomValuesInverse[] =
    {   1.26932f, -1.06467f, -1.28787f, 1.11023f, 1.72369f, -0.749837f, 0.0459759f, -0.41348f, -0.0443016f, 1.36936f,
        -0.110096f, -0.605633f, -1.67664f, 0.591171f, 1.39484f, 0.283764f
    };

    CSRSparseMatrix<T> randomMatrixInverse;
    randomMatrixInverse.setRawDenseData( n, n, randomValuesInverse );

    return randomMatrixInverse;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4IdentityMatrix()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType identityValues[] =
    { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };

    CSRSparseMatrix<T> n4m4IdentM;
    n4m4IdentM.setRawDenseData( n, n, identityValues );

    return n4m4IdentM;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::nnIdentityMatrix( IndexType n )
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    boost::scoped_array<ValueType> identityValues( new ValueType[n * n] );

    for ( int i = 1; i <= n * n; ++i )
    {
        identityValues[i - 1] = 0.0;
    }

    for ( int i = 0; i < n; ++i )
    {
        identityValues[i * n + i] = 1.0;
    }

    CSRSparseMatrix<T> nnIdentM;
    nnIdentM.setRawDenseData( n, n, identityValues.get() );

    return nnIdentM;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4DiagonalMatrix()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.0f };

    CSRSparseMatrix<T> m;
    m.setRawDenseData( n, n, values );

    return m;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n6m6TestMatrix()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 6;

    ValueType values[] =
    {   0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };

    CSRSparseMatrix<T> sparseM;
    sparseM.setRawDenseData( n, n, values );

    return sparseM;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4SymmetricMatrix()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 1.0f, 0.0f, 0.0f, 5.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f, 5.0f, 0.0f, 0.0f, 4.0f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, n, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m6NoneSquareMatrix()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;
    const IndexType m = 6;

    ValueType values[] =
    {   0.436213f, 0.0f, 0.0f, 0.422152f, 0.4632f, 0.168648f, 0.0f, 0.0f, 0.0f, 0.708545f, 0.131853f, 0.820422f, 0.1f,
        0.0f, 0.131853f, 0.820422f, 0.0f, 0.1f, 0.0f, 0.820422f, 0.0f, 0.202542f, 0.0f, 0.947076f
    };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4MatrixA1()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 0.6f, 0.0f, 0.0f, 0.4f, 0.7f, 0.4f, 0.0f, 0.0f, 0.0f, 0.0f, 0.9f, 0.4f, 0.2f, 0.5f, 0.0f, 0.3f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, n, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4MatrixA2()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 0.1f, 0.4f, 0.0f, 0.0f, 0.0f, 0.1f, 0.0f, 0.0f, 0.9f, 0.0f, 0.9f, 0.0f, 0.0f, 0.0f, 0.2f, 0.1f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, n, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4MatrixARes()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 0.06f, 0.24f, 0.08f, 0.04f, 0.07f, 0.32f, 0.0f, 0.0f, 0.81f, 0.0f, 0.89f, 0.04f, 0.02f, 0.13f, 0.06f, 0.03f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, n, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4MatrixB1()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.0f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, n, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4MatrixB2()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 4.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f, 0.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, n, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4MatrixBRes()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 4.0f, 0.0f, 0.0f, 0.0f, 0.0f, 6.0f, 0.0f, 0.0f, 0.0f, 0.0f, 6.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.0f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, n, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4MatrixC1()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 2.0f, -1.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, -1.0f, 2.0f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, n, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4MatrixCRes()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;

    ValueType values[] =
    { 5.0f, -4.0f, 1.0f, 0.0f, -4.0f, 6.0f, -4.0f, 1.0f, 1.0f, -4.0f, 6.0f, -4.0f, 0.0f, 1.0f, -4.0f, 5.0f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, n, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n6m4MatrixD1()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 6;
    const IndexType m = 4;

    ValueType values[] =
    {   6.0f, 0.0f, 0.0f, 4.0f, 7.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f, 9.0f, 4.0f, 2.0f, 5.0f, 0.0f, 3.0f, 2.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f, 2.0f
    };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m6MatrixD2()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;
    const IndexType m = 6;

    ValueType values[] =
    {   1.0f, 4.0f, 0.0f, 0.0f, 4.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 9.0f, 1.0f, 9.0f, 0.0f, 9.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 2.0f, 1.0f, 0.0f, 3.0f
    };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n6m6MatrixDRes()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 6;
    const IndexType m = 6;

    ValueType values[] =
    {   6.0f, 24.0f, 8.0f, 4.0f, 24.0f, 12.0f, 7.0f, 32.0f, 0.0f, 0.0f, 64.0f, 4.0f, 81.0f, 0.0f, 89.0f, 4.0f, 0.0f,
        12.0f, 2.0f, 13.0f, 6.0f, 3.0f, 53.0f, 14.0f, 2.0f, 8.0f, 2.0f, 1.0f, 8.0f, 3.0f, 0.0f, 1.0f, 4.0f, 2.0f, 9.0f,
        7.0f
    };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n6m4MatrixE1()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 6;
    const IndexType m = 4;

    ValueType values[] =
    {   6.0f, 0.0f, 0.0f, 4.0f, 7.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f, 9.0f, 4.0f, 2.0f, 5.0f, 0.0f, 3.0f, 2.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f, 2.0f
    };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m3MatrixE2()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    const IndexType n = 4;
    const IndexType m = 3;

    ValueType values[] =
    { 1.0f, 4.0f, 0.0f, 0.0f, 1.0f, 0.0f, 9.0f, 0.0f, 9.0f, 0.0f, 0.0f, 2.0f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n6m3MatrixERes()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    IndexType n = 6;
    IndexType m = 3;

    ValueType values[] =
    { 6.0f, 24.0f, 8.0f, 7.0f, 32.0f, 0.0f, 81.0f, 0.0f, 89.0f, 2.0f, 13.0f, 6.0f, 2.0f, 8.0f, 2.0f, 0.0f, 1.0f, 4.0f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n8m4Interpol()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    IndexType n = 8;
    IndexType m = 4;

    ValueType values[] =
    {   0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.5f, 0.5f,
        0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f
    };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;

}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m8InterpolTranspose()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    IndexType n = 4;
    IndexType m = 8;

    ValueType values[] =
    {   0.5f, 1.0f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.5f, 1.0f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.5f, 1.0f, 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.5f, 1.0f
    };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;

}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n8m8Laplace1D()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    IndexType n = 8;
    IndexType m = 8;

    ValueType values[] =
    {   2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f
    };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;

}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n8m4GalerkinTemp()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    IndexType n = 8;
    IndexType m = 4;

    ValueType values[] =
    {   0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -0.5f, 1.0f, -0.5f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, -0.5f, 1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -0.5f, 1.5f
    };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;
}

template<typename T>
CSRSparseMatrix<T> TestSparseMatrices::n4m4Galerkin()
{
    typedef typename CSRSparseMatrix<T>::ValueType ValueType;

    IndexType n = 4;
    IndexType m = 4;

    ValueType values[] =
    { 1.0f, -0.5f, 0.0f, 0.0f, -0.5f, 1.0f, -0.5f, 0.0f, 0.0f, -0.5f, 1.0f, -0.5f, 0.0f, 0.0f, -0.5f, 1.5f };

    CSRSparseMatrix<T> matrix;
    matrix.setRawDenseData( n, m, values );

    return matrix;

}

#endif // LAMA_TEST_SPARSE_MATRICES_HPP_
