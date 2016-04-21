
#pragma once

#include <scai/common/SCAITypes.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/lama/storage/MatrixStorage.hpp>

template<typename ValueType>
void getMatrix_7_4 ( IndexType& numRows,
                     IndexType& numColumns,
                     scai::utilskernel::LArray<IndexType>& matrixRowSizes,
                     scai::utilskernel::LArray<IndexType>& matrixJA,
                     scai::utilskernel::LArray<ValueType>& matrixValues,
                     scai::utilskernel::LArray<ValueType>& denseValues )
{
    // Attention: ia array is not an offset array, it contains number of values in each row

    IndexType ia[] = { 2, 1, 2, 3, 2, 0, 2 };
    ValueType values[] = { 6.0f, 4.0f, 7.0f, -9.3f, 4.0f, 2.0f, 5.0f, 3.0f, 2.0f, 1.0f, 1.0f, 2.0f };
    IndexType ja[] = { 0, 3, 0, 2, 3, 0, 1, 3, 0, 3, 1, 3 };

    numRows = sizeof( ia ) / sizeof( IndexType );

    size_t numValues = sizeof( ja ) / sizeof( IndexType );

    BOOST_REQUIRE_EQUAL( (size_t ) numValues, sizeof( values ) / sizeof( ValueType ) );

    ValueType resultMatrix[] =
    {
        6, 0, 0, 4,
        7, 0, 0, 0,
        0, 0, -9.3f , 4,
        2, 5, 0, 3,
        2, 0, 0, 1,
        0, 0, 0, 0,
        0, 1, 0, 2
    };

    matrixRowSizes.init( ia, numRows );
    matrixJA.init( ja, numValues );
    matrixValues.init( values, numValues );

    numColumns = matrixJA.max() + 1;

    BOOST_REQUIRE_EQUAL( static_cast<size_t>( numRows * numColumns ), sizeof( resultMatrix ) / sizeof( ValueType ) );

    denseValues.init( resultMatrix, numRows * numColumns );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void setDenseData( scai::lama::MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 8;
    const IndexType numColumns = 2;
    static ValueType values[] = { 6, 0, 0, 4, 7, 0, 0, 0, 0, 0, -9.3f, 4, 2, 5, 0, 3 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void setDenseSquareData( scai::lama::MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    // just make sure that number of entries in values matches the matrix size
    static ValueType values[] = { 10, 0, 0, 4, 3, 10, 0, 0, 0, 0, -9.3f, 4, 1, 5, 0, 13 };
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    ValueType eps = static_cast<ValueType>( 1E-5 );
    storage.setRawDenseData( numRows, numColumns, values, eps );

    // Note: diagonal does not contain any zeros
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void setDenseRandom( scai::lama::MatrixStorage<ValueType>& storage )
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

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void setDenseHalo( scai::lama::MatrixStorage<ValueType>& storage )
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

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void setSymDenseData( scai::lama::MatrixStorage<ValueType>& storage )
{
    /* Matrix:     1  2  0  5
     *             2  1  3  0
     *             0  3  1  4
     *             5  0  4  2
     */

    const IndexType numRows = 4;
    const IndexType numColumns = 4;

    static ValueType values[] = { 1, 2, 0, 5, 2, 1, 3, 0, 0, 3, 1, 4, 5, 0, 4, 2 };

    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );

    ValueType eps = static_cast<ValueType>( 1E-5 );

    storage.setRawDenseData( numRows, numColumns, values, eps );
}

