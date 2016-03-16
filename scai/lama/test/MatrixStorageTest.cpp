/**
 * @file MatrixStorageTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Contains the implementation of the class MatrixStorageTest
 * @author Thomas Brandes
 * @date 02.03.2012
 */

#include <scai/lama/test/MatrixStorageTest.hpp>
#include <scai/lama/test/SameMatrixHelper.hpp>

#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/ELLStorage.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/lama/test/TestMacros.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/SCAITypes.hpp>

#include <iostream>
#include <fstream>

using namespace scai;
using namespace lama;
using namespace utilskernel;
using namespace hmemo;

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, MatrixStorageTest<ValueType>::logger, "Test.MatrixStorageTest" )

/* ========================================================================= */
/* Implementations of some helpfunction                                      */
/* ========================================================================= */

template<typename ValueType>
void getMatrix_7_4 ( IndexType& numRows,
                     IndexType& numColumns,
                     LArray<IndexType>& matrixRowSizes,
                     LArray<IndexType>& matrixJA,
                     LArray<ValueType>& matrixValues,
                     LArray<ValueType>& denseValues )
{
    numRows = 7;
    numColumns = 4;
    const IndexType numValues = 12;
    // Attention: ia array is not an offset array, it contains number of values in each row
    IndexType ia[] = { 2, 1, 2, 3, 2, 0, 2 };
    ValueType values[] = { 6.0f, 4.0f, 7.0f, -9.3f, 4.0f, 2.0f, 5.0f, 3.0f, 2.0f, 1.0f, 1.0f, 2.0f };
    IndexType ja[] = { 0, 3, 0, 2, 3, 0, 1, 3, 0, 3, 1, 3 };
    ValueType Resultmatrix[] =
    {
        6, 0, 0, 4,
        7, 0, 0, 0,
        0, 0, -9.3f , 4,
        2, 5, 0, 3,
        2, 0, 0, 1,
        0, 0, 0, 0,
        0, 1, 0, 2
    };
    matrixRowSizes = LArray<IndexType>( numRows, ia );
    matrixJA = LArray<IndexType>( numValues, ja );
    matrixValues = LArray<ValueType>( numValues, values );
    denseValues = LArray<ValueType>( numRows * numColumns, Resultmatrix );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorageTest<ValueType>::setDenseData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    static ValueType values[] = { 6, 0, 0, 4, 7, 0, 0, 0, 0, 0, -9.3f, 4, 2, 5, 0, 3 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorageTest<ValueType>::setDenseDataNotSquare( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 6;
    const IndexType numColumns = 4;
    static ValueType values[] = { 6, 0, 0, 4, 7, 0, 0, 0, 0, 0, -9.3f, 4, 2, 5, 0, 3, 0, 1, 2, 0, 8, 0, 0, 1 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorageTest<ValueType>::setDenseDataSymmetric( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    static ValueType values[] = { 1, 2, 0, 5, 2, 1, 3, 0, 0, 3, 1, 4, 5, 0, 4, 1 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorageTest<ValueType>::setDenseLocal( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    static ValueType values[] = { 10, 0, 0, 4, 3, 10, 0, 0, 0, 0, -9.3f, 4, 1, 5, 0, 13 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorageTest<ValueType>::setDenseHalo( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 3;
    static ValueType values[] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorageTest<ValueType>::setDenseRandom( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    static ValueType values[] =
    {
        0.436213f, 0.683202f, 0.531013f, 0.422152f,
        0.4632f,   0.168648f, 0.967549f, 0.498486f,
        0.126115f, 0.708545f, 0.131853f, 0.820422f,
        0.992481f, 0.202542f, 0.47369f,  0.947076f
    };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorageTest<ValueType>::setDenseRandomInverse( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    static ValueType values[] =
    {
        1.26932f, -1.06467f, -1.28787f, 1.11023f, 1.72369f, -0.749837f, 0.0459759f, -0.41348f, -0.0443016f, 1.36936f,
        -0.110096f, -0.605633f, -1.67664f, 0.591171f, 1.39484f, 0.283764f
    };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ========================================================================= */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, emptyTest )
{
    SCAI_LOG_INFO( logger, "emptyTest" );
    mMatrixStorage.clear();
    // verify that empty matrix has diagonal property
    BOOST_CHECK( mMatrixStorage.hasDiagonalProperty() );
    mMatrixStorage.allocate( 1, 1 );

    if ( mMatrixStorage.getFormat() == Format::DENSE )
    {
        // only dense matrix keeps its diagonal property
        BOOST_CHECK( mMatrixStorage.hasDiagonalProperty() );
    }
    else
    {
        BOOST_CHECK( ! mMatrixStorage.hasDiagonalProperty() );
    }

    mMatrixStorage.purge();
   BOOST_CHECK( mMatrixStorage.hasDiagonalProperty() );
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ========================================================================= */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, setIdentityTest )

const IndexType n = 15;

SCAI_LOG_INFO( logger, "setIdentity, uses n = " << n );

mMatrixStorage.clear();

mMatrixStorage.setIdentity( n );

BOOST_REQUIRE_EQUAL( n, mMatrixStorage.getNumRows() );
BOOST_REQUIRE_EQUAL( n, mMatrixStorage.getNumColumns() );
BOOST_REQUIRE_EQUAL( n, mMatrixStorage.getNumValues() );

LArray<double> row;

for ( IndexType i = 0; i < n; ++i )
{
    mMatrixStorage.getRow( row, i );
    ReadAccess<double> rRow( row );

    for ( IndexType j = 0; j < n; ++j )
    {
        if ( i == j )
        {
            BOOST_CHECK_EQUAL( 1.0, rRow[j] );
        }
        else
        {
            BOOST_CHECK_EQUAL( 0.0, rRow[j] );
        }
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, setCSRDataTest )
{
    SCAI_LOG_INFO( logger, "setCSRDataTest" )
    mMatrixStorage.clear();
    IndexType numRows;
    IndexType numColumns;
    LArray<IndexType> matrixRowSizes;
    LArray<IndexType> matrixJA;
    LArray<ValueType> matrixValues;
    LArray<ValueType> matrixDense;
    getMatrix_7_4 ( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues, matrixDense );
    IndexType numValues = matrixJA.size();
    _MatrixStorage::sizes2offsets( matrixRowSizes );
// Now we can use matrixRowSizes as IA array
    mMatrixStorage.setCSRData( numRows, numColumns, numValues, matrixRowSizes, matrixJA, matrixValues );
    SCAI_LOG_INFO( logger, "set CSR data (" << numRows << " x " << numColumns
                   << ", nnz = " << numValues << ") : matrix = " << mMatrixStorage )
    mMatrixStorage.prefetch();
    ReadAccess<ValueType> results( matrixDense );

    for ( IndexType i = 0; i < numRows; i++ )
    {
        for ( IndexType j = 0; j < numColumns; j++ )
        {
            BOOST_CHECK_EQUAL( mMatrixStorage.getValue( i, j ) , results[ i * numColumns + j ] );
        }
    }

    mMatrixStorage.wait();
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, buildCSRDataTest )
{
    mMatrixStorage.clear();
    IndexType numRows;
    IndexType numColumns;
    LArray<IndexType> matrixRowSizes;
    LArray<IndexType> matrixJA;
    LArray<ValueType> matrixValues;
    LArray<ValueType> matrixDense;
    getMatrix_7_4 ( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues, matrixDense );
    IndexType numValues = matrixJA.size();
// IA array not available yet, we have only sizes
    _MatrixStorage::sizes2offsets( matrixRowSizes );
// Now we can use matrixRowSizes as IA array
    mMatrixStorage.setCSRData( numRows, numColumns, numValues, matrixRowSizes, matrixJA, matrixValues );
// make sure that storage has not diagonal property, otherwise it will build wrong CSR data
    BOOST_REQUIRE_EQUAL( mMatrixStorage.hasDiagonalProperty(), false );
// make sure that we have all values stored
    BOOST_CHECK_EQUAL( numValues, mMatrixStorage.getNumValues() );
    SCAI_LOG_INFO( logger, "set CSR data (" << numRows << " x " << numColumns
                   << ", nnz = " << numValues << ") : matrix = " << mMatrixStorage )
    LArray<IndexType> csrIA;
    LArray<IndexType> csrJA;
    LArray<ValueType> csrValues;
    mMatrixStorage.buildCSRData( csrIA, csrJA, csrValues );
    BOOST_REQUIRE_EQUAL( numRows + 1, csrIA.size() );
    BOOST_REQUIRE_EQUAL( numValues, csrJA.size() );
    BOOST_REQUIRE_EQUAL( numValues, csrValues.size() );

    // check the IA array ( csrIA are offsets, matrixRowSizes was converted to offsets

    {
        ReadAccess<IndexType> ia1( matrixRowSizes );
        ReadAccess<IndexType> ia2( csrIA );

        for ( IndexType i = 0; i <= numRows; ++i )
        {
            BOOST_CHECK_EQUAL( ia1[i], ia2[i] );
        }
    }

    ReadAccess<IndexType> ja1( matrixJA );
    ReadAccess<IndexType> ja2( csrJA );
    ReadAccess<ValueType> values1( matrixValues );
    ReadAccess<ValueType> values2( csrValues );

    for ( IndexType i = 0; i < numValues; i++ )
    {
        BOOST_CHECK_EQUAL( ja1[i], ja2[i] );
        BOOST_CHECK_EQUAL( values1[i], values2[i] );
    }
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, purgeTest )

SCAI_LOG_INFO( logger, "PurgeTest" )

setDenseData( mMatrixStorage );

size_t usedBytes = mMatrixStorage.getMemoryUsage();

mMatrixStorage.purge();

// Note: clear() also decreases the usage, but will not free memory

size_t usedBytesPurged = mMatrixStorage.getMemoryUsage();

// make sure that really some data has been freed

BOOST_CHECK ( usedBytesPurged < usedBytes );

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, diagonalTest )

setDenseData( mMatrixStorage );

SCAI_LOG_INFO( logger, "diagonalTest: matrix = " << mMatrixStorage )

LArray<ValueType> diag;
mMatrixStorage.getDiagonal( diag );

BOOST_CHECK_EQUAL( diag.size(), mMatrixStorage.getNumRows() ); // square matrix

SCAI_LOG_INFO( logger, "diagonalTest: get diagonal = " << diag )

{
    ReadAccess<ValueType> rDiag ( diag );

    for ( IndexType i = 0; i < diag.size(); ++i )
    {
        Scalar s = mMatrixStorage.getValue( i, i );
        BOOST_CHECK_EQUAL( s.getValue<ValueType>(), rDiag[i] );
    }
}

mMatrixStorage.setDiagonal( 0.0 );

for ( IndexType i = 0; i < diag.size(); ++i )
{
    Scalar s = mMatrixStorage.getValue( i, i );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 0.0 );
}

mMatrixStorage.setDiagonalV( diag );

{
    ReadAccess<ValueType> rDiag ( diag );

    for ( IndexType i = 0; i < diag.size(); ++i )
    {
        Scalar s = mMatrixStorage.getValue( i, i );
        BOOST_CHECK_EQUAL( s.getValue<ValueType>(), rDiag[i] );
    }
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, scaleTest )

setDenseData( mMatrixStorage );

DenseStorage<ValueType> tmp;

setDenseData( tmp );

mMatrixStorage.scale( 2.0 );

for ( IndexType i = 0; i < mMatrixStorage.getNumRows(); ++i )
{
    for ( IndexType j = 0; j < mMatrixStorage.getNumColumns(); ++j )
    {
        SCAI_CHECK_CLOSE( 2.0 * tmp.getValue( i, j ),
                          mMatrixStorage.getValue( i, j ), 1 );
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, normTest )

setDenseData( mMatrixStorage );

SCAI_LOG_INFO( logger, "maxNormTest, mMatrixStorage = " << mMatrixStorage )

ValueType maxNorm = mMatrixStorage.maxNorm();

ValueType expected = 9.3f; // maximal absolute value

SCAI_CHECK_CLOSE( maxNorm, expected, 1 );

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, vectorMultTest )

setDenseData( mMatrixStorage );

DenseStorage<ValueType> orig; // used for comparing results
setDenseData( orig );

SCAI_LOG_INFO( logger, "Test matrixTimesVector" )
{
    const ValueType xVal = 1.0;
    const ValueType yVal = 1.0;
    const ValueType alpha = 1.0;
    const ValueType beta = 2.0;
    LArray<ValueType> x( mMatrixStorage.getNumColumns(), xVal );
    LArray<ValueType> y( mMatrixStorage.getNumRows(), yVal );
    // due to use of LArrayView we have to give result the correct size
    LArray<ValueType> result ( mMatrixStorage.getNumRows() );
    mMatrixStorage.matrixTimesVector( result, alpha, x, beta, y );
    BOOST_CHECK_EQUAL( result.size(), mMatrixStorage.getNumRows() );
    ReadAccess<ValueType> values( orig.getData() );
    ReadAccess<ValueType> res( result );
    int ncol = mMatrixStorage.getNumColumns();
    int nrow = mMatrixStorage.getNumRows();

    for ( IndexType i = 0; i < nrow; ++i )
    {
        ValueType sum = 0.0;

        for ( IndexType j = 0; j < ncol; ++j )
        {
            sum += values[ i * ncol + j] * alpha * xVal;
        }

        sum += beta * yVal;
        SCAI_CHECK_CLOSE( sum, res[i], 0.1 );
    }
}

SCAI_LOG_INFO( logger, "Test " << mMatrixStorage.getTypeName() << "::matrixTimesVectorAsync" )
{
    const ValueType xVal = static_cast<ValueType>( 1.5 );
    const ValueType yVal = static_cast<ValueType>( 2.0 );
    const ValueType alpha = static_cast<ValueType>( 0.5f );
    const ValueType beta = static_cast<ValueType>( 1.0f );
    LArray<ValueType> x( mMatrixStorage.getNumColumns(), xVal );
    LArray<ValueType> y( mMatrixStorage.getNumRows(), yVal );
    LArray<ValueType> result( mMatrixStorage.getNumRows() );
    // asynchronous execution, only checks correct calling
    {
        common::unique_ptr<tasking::SyncToken> token ( mMatrixStorage.matrixTimesVectorAsync( result, alpha, x, beta, y ) );
        // unique pointer implies delete for token so 
        // destructor of the token does the synchronization at the end of this scope
    }
    BOOST_CHECK_EQUAL( result.size(), mMatrixStorage.getNumRows() );
    ReadAccess<ValueType> values( orig.getData() );
    ReadAccess<ValueType> res( result );
    int ncol = mMatrixStorage.getNumColumns();
    int nrow = mMatrixStorage.getNumRows();

    for ( IndexType i = 0; i < nrow; ++i )
    {
        ValueType sum = 0.0;

        for ( IndexType j = 0; j < ncol; ++j )
        {
            sum += values[ i * ncol + j ] * alpha * xVal;
        }

        sum += beta * yVal;
        SCAI_CHECK_CLOSE( sum, res[i], 0.1f );
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, vectorTimesMatrixTest )

setDenseData( mMatrixStorage );

DenseStorage<ValueType> orig; // used for comparing results
setDenseData( orig );

SCAI_LOG_INFO( logger, "Test vectorTimesMatrix" )
{
    const ValueType xVal = 1.0;
    const ValueType yVal = 1.0;
    const ValueType alpha = 1.5;
    const ValueType beta = 2.3;
    IndexType m = mMatrixStorage.getNumColumns();
    IndexType n = mMatrixStorage.getNumRows();
    LArray<ValueType> x( n, xVal );
    LArray<ValueType> y( m, yVal );
    // due to use of LArrayView we have to give result the correct size
    LArray<ValueType> result ( m );
    mMatrixStorage.vectorTimesMatrix( result, alpha, x, beta, y );
    BOOST_CHECK_EQUAL( result.size(), m );
    ReadAccess<ValueType> values( orig.getData() );
    ReadAccess<ValueType> res( result );

    for ( IndexType j = 0; j < m; ++j )
    {
        ValueType sum = 0.0;

        for ( IndexType i = 0; i < n; ++i )
        {
            sum += values[ i * m + j ] * alpha * xVal;
        }

        sum += beta * yVal;
        SCAI_CHECK_CLOSE( sum, res[j], 0.1f );
    }
}

SCAI_LOG_INFO( logger, "Test vectorTimesMatrixAsync" )
{
    const ValueType xVal = 1.0;
    const ValueType yVal = 1.0;
    const ValueType alpha = 1.0;
    const ValueType beta = 2.0;
    IndexType m = mMatrixStorage.getNumColumns();
    IndexType n = mMatrixStorage.getNumRows();
    LArray<ValueType> x( n, xVal );
    LArray<ValueType> y( m, yVal );
    // due to use of LArrayView we have to give result the correct size
    LArray<ValueType> result ( m );
    // asynchronous execution, only checks correct calling
    {
        common::unique_ptr<tasking::SyncToken> token ( mMatrixStorage.vectorTimesMatrixAsync( result, alpha, x, beta, y ) );
        // free of token at end of this scope does the synchronization
    }
    SCAI_LOG_TRACE( logger, "vectorTimesMatrixAsync synchronized" )
    BOOST_CHECK_EQUAL( result.size(), m );
    ReadAccess<ValueType> values( orig.getData() );
    ReadAccess<ValueType> res( result );

    for ( IndexType j = 0; j < m; ++j )
    {
        ValueType sum = 0.0;

        for ( IndexType i = 0; i < n; ++i )
        {
            sum += values[ i * m + j ] * alpha * xVal;
        }

        sum += beta * yVal;
        SCAI_CHECK_CLOSE( sum, res[j], 0.1f );
    }
}

SCAI_LOG_INFO( logger, "Test vectorTimesMatrixAsync 2" )
{
    IndexType n = mMatrixStorage.getNumRows();
    ValueType* xValues = ( ValueType* ) malloc( n * sizeof( ValueType ) );

    for ( IndexType i = 0; i < n; ++i )
    {
        xValues[i] = static_cast<ValueType>( i + 1 );
    }

    IndexType m = mMatrixStorage.getNumColumns();
    ValueType* yValues = ( ValueType* ) malloc( m * sizeof( ValueType ) );

    for ( IndexType i = 0; i < m; ++i )
    {
        yValues[i] = static_cast<ValueType>( m - i );
    }

    const ValueType alpha = 1.0;
    const ValueType beta = 2.0;
    LArray<ValueType> x( n, xValues );
    LArray<ValueType> y( m, yValues );
    // due to use of LArrayView we have to give result the correct size
    LArray<ValueType> result ( m );
    // asynchronous execution, only checks correct calling
    {
        common::unique_ptr<tasking::SyncToken> token ( mMatrixStorage.vectorTimesMatrixAsync( result, alpha, x, beta, y ) );
        // free of token at end of this scope does the synchronization
    }
    SCAI_LOG_TRACE( logger, "vectorTimesMatrixAsync synchronized" )
    BOOST_CHECK_EQUAL( result.size(), m );
    ReadAccess<ValueType> values( orig.getData() );
    ReadAccess<ValueType> res( result );

    for ( IndexType j = 0; j < m; ++j )
    {
        ValueType sum = 0.0;

        for ( IndexType i = 0; i < n; ++i )
        {
            sum += values[ i * m + j ] * alpha * xValues[i];
        }

        sum += beta * yValues[j];
        SCAI_CHECK_CLOSE( sum, res[j], 0.1f );
    }
}

setDenseDataNotSquare( mMatrixStorage );
setDenseDataNotSquare( orig );

SCAI_LOG_INFO( logger, "Test vectorTimesMatrixAsync 3" )
{
    IndexType n = mMatrixStorage.getNumRows();
    ValueType* xValues = ( ValueType* ) malloc( n * sizeof( ValueType ) );

    for ( IndexType i = 0; i < n; ++i )
    {
        xValues[i] = static_cast<ValueType>( i + 1 );
    }

    IndexType m = mMatrixStorage.getNumColumns();
    ValueType* yValues = ( ValueType* ) malloc( m * sizeof( ValueType ) );

    for ( IndexType i = 0; i < m; ++i )
    {
        yValues[i] = static_cast<ValueType>( m - i );
    }

    const ValueType alpha = 1.0;
    const ValueType beta = 2.0;
    LArray<ValueType> x( n, xValues );
    LArray<ValueType> y( m, yValues );
    // due to use of LArrayView we have to give result the correct size
    LArray<ValueType> result ( m );
    // asynchronous execution, only checks correct calling
    {
        common::unique_ptr<tasking::SyncToken> token ( mMatrixStorage.vectorTimesMatrixAsync( result, alpha, x, beta, y ) );
        // free of token at end of this scope does the synchronization
    }
    SCAI_LOG_TRACE( logger, "vectorTimesMatrixAsync synchronized" )
    BOOST_CHECK_EQUAL( result.size(), m );
    ReadAccess<ValueType> values( orig.getData() );
    ReadAccess<ValueType> res( result );

    for ( IndexType j = 0; j < m; ++j )
    {
        ValueType sum = 0.0;

        for ( IndexType i = 0; i < n; ++i )
        {
            sum += values[ i * m + j ] * alpha * xValues[i];
        }

        sum += beta * yValues[j];
        SCAI_CHECK_CLOSE( sum, res[j], 0.1f );
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, numericalTest )

std::string prefix = scai::test::Configuration::getPath();
CSRSparseMatrix<ValueType> symA( prefix + "/" + "nos6.mtx" );
mMatrixStorage = symA.getLocalStorage();

DenseStorage<ValueType> orig; // used for comparing results
orig = mMatrixStorage;

DenseVector<ValueType> vec( prefix + "/" + "p2" );

bool symmetry = symA.checkSymmetry();
BOOST_CHECK_EQUAL( symmetry, true );

IndexType m = mMatrixStorage.getNumColumns();

const ValueType alpha = 1.0;
const ValueType beta = 0.0;

LArray<ValueType> x( vec.getLocalValues() );
LArray<ValueType> y( m, 0.0 );
// due to use of LArrayView we have to give result the correct size
LArray<ValueType> result ( m );
LArray<ValueType> result2 ( m );

DenseStorage<ValueType> denseStorage ( symA.getLocalStorage() ); //mMatrixStorage );
LArray<ValueType> denseResult ( m );
LArray<ValueType> denseResult2 ( m );

// asynchronous execution, only checks correct calling

{
    common::unique_ptr<tasking::SyncToken> token ( mMatrixStorage.vectorTimesMatrixAsync( result, alpha, x, beta, y ) );
    common::unique_ptr<tasking::SyncToken> token2 ( mMatrixStorage.matrixTimesVectorAsync( result2, alpha, x, beta, y ) );
    common::unique_ptr<tasking::SyncToken> token3 ( denseStorage.vectorTimesMatrixAsync( denseResult, alpha, x, beta, y ) );
    common::unique_ptr<tasking::SyncToken> token4 ( denseStorage.matrixTimesVectorAsync( denseResult2, alpha, x, beta, y ) );
    // free of token at end of this scope does the synchronization
}
SCAI_LOG_TRACE( logger, "vectorTimesMatrixAsync and matrixTimesVectorAsync synchronized" )

BOOST_CHECK_EQUAL( result.size(), mMatrixStorage.getNumColumns() );

ReadAccess<ValueType> xValues( x );
ReadAccess<ValueType> values( orig.getData() );
ReadAccess<ValueType> res( result );
ReadAccess<ValueType> res2( result2 );
ReadAccess<ValueType> resD( denseResult );
ReadAccess<ValueType> res2D( denseResult2 );

int ncol = mMatrixStorage.getNumColumns();
int nrow = mMatrixStorage.getNumRows();

for ( IndexType j = 0; j < ncol; ++j )
{
    ValueType sum = 0.0;

    if ( symA.hasDiagonalProperty() )
    {
        sum += values[ j * ncol ] * xValues[j];

        for ( IndexType i = 0; i < nrow; ++i )
        {
            if ( i != j )
            {
                return;
            }

            sum += values[ i * ncol + j ] * xValues[i];
        }
    }
    else
    {
        for ( IndexType i = 0; i < nrow; ++i )
        {
            sum += values[ i * ncol + j ] * xValues[i];
        }
    }

    sum *= alpha;
    SCAI_CHECK_CLOSE( sum, res[j], 0.1f );
    SCAI_CHECK_CLOSE( res[j], res2[j], 0.1f );
    SCAI_CHECK_CLOSE( resD[j], res2D[j], 0.1f );
    SCAI_CHECK_CLOSE( resD[j], res[j], 0.1f );
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, matrixMultTest )

SCAI_LOG_INFO( logger, "Test " << typeid( ValueType ).name() << "::matrixTimesMatrix" )

const IndexType n = 3;
const ValueType aValue = 2.0;
const ValueType bValue = 3.0;

LArray<ValueType> aDiagonal( n, aValue );
LArray<ValueType> bDiagonal( n, bValue );

//TODO: USE TEMPLATED STORAGE!!!!!!
ELLStorage<ValueType> a;
ELLStorage<ValueType> b;
a.setIdentity( n );
b.setIdentity( n );
a.setDiagonal( aValue );
b.setDiagonal( bValue );

mMatrixStorage.matrixTimesMatrix( 1.0, a, b, 0.0, mMatrixStorage );

BOOST_CHECK_EQUAL( n, mMatrixStorage.getNumRows() );
BOOST_CHECK_EQUAL( n, mMatrixStorage.getNumColumns() );

for ( IndexType i = 0; i < n; ++i )
{
    for ( IndexType j = 0; j < n; ++j )
    {
        if ( i == j )
        {
            SCAI_CHECK_CLOSE( mMatrixStorage.getValue( i, j ), aValue * bValue, 1 );
        }
        else
        {
            SCAI_CHECK_CLOSE( mMatrixStorage.getValue( i, j ), 0.0, 1 );
        }
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, matrixAddTest )

SCAI_LOG_INFO( logger, "Test " << typeid( ValueType ).name() << "::matrixTimesMatrix" );

IndexType n = 3;
const ValueType aValue = 2.0;
const ValueType bValue = 3.0;
const ValueType cValue = -3.0;

LArray<ValueType> aDiagonal( n, aValue );
LArray<ValueType> bDiagonal( n, bValue );
LArray<ValueType> cDiagonal( n, cValue );

//TODO: USE TEMPLATED STORAGE!!!!!!
ELLStorage<ValueType> a;
ELLStorage<ValueType> b;
ELLStorage<ValueType> c;
a.setIdentity( n );
b.setIdentity( n );
c.setIdentity( n );
a.setDiagonal( aValue );
b.setDiagonal( bValue );
c.setDiagonal( cValue );

mMatrixStorage.matrixTimesMatrix( 1.0, a, b, 1.0, c );

BOOST_CHECK_EQUAL( n, mMatrixStorage.getNumRows() );
BOOST_CHECK_EQUAL( n, mMatrixStorage.getNumColumns() );

for ( IndexType i = 0; i < n; ++i )
{
    for ( IndexType j = 0; j < n; ++j )
    {
        if ( i == j )
        {
            SCAI_CHECK_CLOSE( mMatrixStorage.getValue( i, j ), aValue * bValue + cValue, 1 );
        }
        else
        {
            SCAI_CHECK_CLOSE( mMatrixStorage.getValue( i, j ), 0.0, 1 );
        }
    }
}

setDenseData( mMatrixStorage );

n = mMatrixStorage.getNumRows();

// setDenseData must return square matrix

BOOST_CHECK_EQUAL( n, mMatrixStorage.getNumColumns() );

b.setIdentity( n );

// Note: we have two aliases

mMatrixStorage.matrixTimesMatrix( 1.0, mMatrixStorage, b, -1.0, mMatrixStorage );

for ( IndexType i = 0; i < n; ++i )
{
    for ( IndexType j = 0; j < n; ++j )
    {
        SCAI_CHECK_CLOSE( mMatrixStorage.getValue( i, j ), 0.0, 1 );
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, matrixMultTest1 )

SCAI_LOG_INFO( logger, "matrixMultTest1" )

//TODO: USE TEMPLATED STORAGE!!!!!!
ELLStorage<ValueType> a;
ELLStorage<ValueType> b;
setDenseData( a );
setDenseData( b );

// both are square

const IndexType n = b.getNumColumns();

mMatrixStorage.matrixTimesMatrix( 1.0, a, b, 0.0, mMatrixStorage );

LArray<ValueType> x( n, 1.0 );
LArray<ValueType> dummy( n, 0.0 );
LArray<ValueType> y1( n );
LArray<ValueType> y2( n );
LArray<ValueType> tmp( n );

// check: storage * x == a * b * x

a.matrixTimesVector( tmp, 1.0, x, 0.0, dummy );
b.matrixTimesVector( y1, 1.0, tmp, 0.0, dummy );
mMatrixStorage.matrixTimesVector( y2, 1.0, x, 0.0, dummy );

BOOST_CHECK_EQUAL( y1.size(), y2.size() );

ReadAccess<ValueType> rY1( y1 );
ReadAccess<ValueType> rY2( y2 );

for ( IndexType i = 0; i < n; ++i )
{
    SCAI_CHECK_CLOSE( rY1[i], rY2[i], 1 );
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorageTest<ValueType>::jacobiTest( const ValueType omega )
{
    SCAI_LOG_INFO( logger, "jacobiTest: matrix = " << mMatrixStorage << ", omega = " << omega )
    setDenseLocal( mMatrixStorage );
    CSRStorage<ValueType> tmp;
    setDenseLocal( tmp );
    const IndexType n = mMatrixStorage.getNumRows();
    LArray<ValueType> oldSolution( n, 1.0 );
    LArray<ValueType> rhs( n, 2.0 );
    LArray<ValueType> solution1( n );
    LArray<ValueType> solution2( n );
    SCAI_LOG_INFO( logger, "initialized: solution1 = " << solution1 )
    // jacobi:  solution1 = omega * ( rhs - B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
    mMatrixStorage.jacobiIterate( solution1, oldSolution, rhs, omega );
    LArray<ValueType> diagonal;
    tmp.getDiagonal( diagonal );
    {
        // invert the diagonal
        WriteAccess<ValueType> wDiagonal( diagonal );

        for ( IndexType i = 0; i < n; i++ )
        {
            wDiagonal[i] = static_cast<ValueType>( 1.0 ) / wDiagonal[i];
        }
    }
    tmp.setDiagonal( 0.0 );
    // intermediate step: solution2 = rhs - tmp * oldSolution
    tmp.matrixTimesVector( solution2, -1.0, oldSolution, 1.0, rhs );
    {
        // invert the diagonal
        WriteAccess<ValueType> wSolution( solution2 );
        ReadAccess<ValueType> rDiagonal( diagonal );
        ReadAccess<ValueType> rOldSolution( oldSolution );

        for ( IndexType i = 0; i < n; i++ )
        {
            wSolution[i] = static_cast<ValueType>( omega * wSolution[i] * rDiagonal[i]
                                                   + ( 1.0 - omega ) * rOldSolution[i] );
        }
    }
    // by this way solution1 and solution2 must be equal
    {
        ReadAccess<ValueType> rSolution1( solution1 );
        ReadAccess<ValueType> rSolution2( solution2 );

        for ( IndexType i = 0; i < n; ++i )
        {
            SCAI_CHECK_CLOSE( rSolution1[i], rSolution2[i], 1 );
        }
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, jacobiTest )
jacobiTest( 1.0f );
jacobiTest( 0.8f );
jacobiTest( 0.5f );
LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, jacobiHaloTest )

const ValueType omega = static_cast<ValueType> ( 0.35 );

SCAI_LOG_INFO( logger, "jacobiHaloTest: matrix = " << mMatrixStorage
               << ", omega = " << omega )

setDenseHalo( mMatrixStorage );

common::shared_ptr<MatrixStorage<ValueType> > local( MatrixStorage<ValueType>::create( mMatrixStorage.getCreateValue() ) );

setDenseLocal( *local );

CSRStorage<ValueType> tmp;
setDenseHalo( tmp );

const IndexType numRows = mMatrixStorage.getNumRows();
const IndexType numCols = mMatrixStorage.getNumColumns();

LArray<ValueType> oldSolution ( numCols, 1.0 );

LArray<ValueType> solution1( numRows, 2.0 );
LArray<ValueType> solution2( numRows, 2.0 );

SCAI_LOG_INFO( logger, "initialized: solution1 = " << solution1 )

// jacobi:  solution1 -=  omega * ( B(halo) * oldSolution) * dinv

mMatrixStorage.jacobiIterateHalo( solution1, *local, oldSolution, omega );

LArray<ValueType> diagonal;

local->getDiagonal( diagonal );

{
    // invert the diagonal
    WriteAccess<ValueType> wDiagonal( diagonal );

    for ( IndexType i = 0; i < numRows; i++ )
    {
        wDiagonal[i] = static_cast<ValueType>( 1.0 ) / wDiagonal[i];
    }
}

// scale the halo with the diagonal

tmp.scaleRows( diagonal );

// solution2 = - omega * tmp * oldSolution + solution2

tmp.matrixTimesVector( solution2, -omega, oldSolution, 1.0, solution2 );

// by this way solution1 and solution2 must be equal

{
    ReadAccess<ValueType> rSolution1( solution1 );
    ReadAccess<ValueType> rSolution2( solution2 );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        SCAI_CHECK_CLOSE( rSolution1[i], rSolution2[i], 1 );
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, inverseTest )

SCAI_LOG_INFO( logger, "inverseTest for " << mMatrixStorage.getTypeName() )

setDenseRandom( mMatrixStorage );
common::shared_ptr<MatrixStorage<ValueType> > identity( MatrixStorage<ValueType>::create( mMatrixStorage.getCreateValue() ) );
identity->setIdentity( mMatrixStorage.getNumRows() );

// 1. check: use the inverse of the identity, which is expected to be the identity again
{
    common::shared_ptr<MatrixStorage<ValueType> > identity2( MatrixStorage<ValueType>::create( mMatrixStorage.getCreateValue() ) );
    identity2->setIdentity( mMatrixStorage.getNumRows() );
    identity2->invert( mMatrixStorage );
    // testSameMatrixStorage( *identity, *identity2, Scalar( common::TypeTraits<ValueType>::small() ) );
}

// 2. check: check invert on matrix ( multiply inverted matrix with original matrix and check whether the result is the
// identity
{
    common::shared_ptr<MatrixStorage<ValueType> > inverse( MatrixStorage<ValueType>::create( mMatrixStorage.getCreateValue() ) );
    common::shared_ptr<MatrixStorage<ValueType> > result( MatrixStorage<ValueType>::create( mMatrixStorage.getCreateValue() ) );
    inverse->invert( mMatrixStorage );
    result->matrixTimesMatrix( static_cast<ValueType>(1.0), mMatrixStorage, *inverse, static_cast<ValueType>(0.0), *inverse );
    testSameMatrixStorage( *identity, *result, Scalar( common::TypeTraits<ValueType>::small() ) );
}

// 3. check: check invert on matrix - in-place ( multiply inverted matrix with original matrix and check whether the
// result is the identity
{
    common::shared_ptr<MatrixStorage<ValueType> > inverse( mMatrixStorage.copy() );
    common::shared_ptr<MatrixStorage<ValueType> > result( MatrixStorage<ValueType>::create( mMatrixStorage.getCreateValue() ) );
    inverse->invert( mMatrixStorage );
    result->matrixTimesMatrix( ValueType( 1 ), mMatrixStorage, *inverse, ValueType( 0 ), *inverse);
    testSameMatrixStorage( *identity, *result, Scalar( common::TypeTraits<ValueType>::small() ) );
}

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, symmetryTest )

std::string prefix = scai::test::Configuration::getPath();
CSRSparseMatrix<ValueType> sym( prefix + "/" + "nos6.mtx" );
mMatrixStorage = sym.getLocalStorage();

DenseStorage<ValueType> orig; // used for comparing results
orig = mMatrixStorage;

bool symmetry = sym.checkSymmetry();
BOOST_CHECK_EQUAL( symmetry, true );

CSRSparseMatrix<ValueType> asym( prefix + "/" + "impcol_b.mtx" );
mMatrixStorage = asym.getLocalStorage();

orig = mMatrixStorage;

symmetry = asym.checkSymmetry();
BOOST_CHECK_EQUAL( symmetry, false );

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( MatrixStorageTest, ValueType, writeAtTest )

LAMA_WRITEAT_TEST( mMatrixStorage )

LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ------------------------------------------------------------------------- */

template<typename StorageType>                                                                                     \
void MatrixStorageTest<StorageType>::runTests()
{
    emptyTest();
    purgeTest();
    setIdentityTest();
    setCSRDataTest();
    buildCSRDataTest();
    diagonalTest();
    scaleTest();
    normTest();
    vectorMultTest();
    vectorTimesMatrixTest();
    numericalTest();
    jacobiTest();
    jacobiHaloTest();
    matrixMultTest();
    matrixAddTest();
    matrixMultTest1();
    inverseTest();
    symmetryTest();
    writeAtTest();
}

/* ------------------------------------------------------------------------- */

SCAI_COMMON_INST_CLASS( MatrixStorageTest, ARITHMETIC_HOST_CNT, ARITHMETIC_HOST )
