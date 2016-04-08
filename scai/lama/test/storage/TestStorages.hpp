
#include <scai/common/SCAITypes.hpp>
#include <scai/utilskernel/LArray.hpp>

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
