/**
 * @file TestData1.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Contains test data of one matrix for different formats
 * @author Thomas Brandes
 * @date 20.12.2016
 */

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

namespace scai
{

namespace data1
{

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getDenseTestData(
    IndexType& numRows,
    IndexType& numColumns,
    hmemo::HArray<ValueType>& denseValues )
{
    /*   Matrix:       6  0  0  4         0  3  -    6  4  -
                       7  0  0  0         0  -  -    7  -  -
                       0  0  9  4         2  3  -    9  4  -
                       2  5  0  3         0  1  3    2  5  3
                       2  0  0  1         0  3  -    2  1  -
                       0  0  0  0         -  -  -    -  -  -
                       0  1  0  2         1  3  -    1  2  -
     */

    const ValueType values[]  = { 6, 0, 0, 4,
                                  7, 0, 0, 0,
                                  0, 0, 9, 4,
                                  2, 5, 0, 3,
                                  2, 0, 0, 1,
                                  0, 0, 0, 0,
                                  0, 1, 0, 2
                                };

    numRows    = 7;
    numColumns = 4;

    IndexType numValues  = sizeof( values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numValues, numRows * numColumns, "size mismatch" )

    denseValues.setRawData( numValues, values );
}

/* ------------------------------------------------------------------------------------- */

/** Return  a vector with indexes of all rows that are not empty */

static inline void getRowIndexes( hmemo::HArray<IndexType>& rowIndexes )
{
    const IndexType nonzero_rows[]  = { 0, 1, 2, 3, 4, 6 };   // only row 5 is empty

    const IndexType numNonEmptyRows = sizeof( nonzero_rows ) / sizeof( IndexType );

    rowIndexes.setRawData( numNonEmptyRows, nonzero_rows );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getCSRTestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    hmemo::HArray<IndexType>& csrIA,
    hmemo::HArray<IndexType>& csrJA,
    hmemo::HArray<ValueType>& csrValues )
{
    /*   Matrix:       6  0  0  4         0  3  -    6  4  -
                       7  0  0  0         0  -  -    7  -  -
                       0  0  9  4         2  3  -    9  4  -
                       2  5  0  3         0  1  3    2  5  3
                       2  0  0  1         0  3  -    2  1  -
                       0  0  0  0         -  -  -    -  -  -
                       0  1  0  2         1  3  -    1  2  -
     */

    const IndexType ia_values[]  = { 0,    2, 3,    5,       8,    10, 10,   12 };
    const IndexType ja_values[]  = { 0, 3, 0, 2, 3, 0, 1, 3, 0, 3,     1, 3 };
    const ValueType nz_values[]  = { 6, 4, 7, 9, 4, 2, 5, 3, 2, 1,     1, 2 };

    numRows    = sizeof( ia_values ) / sizeof( IndexType ) - 1;
    numColumns = 4;
    numValues  = sizeof( ja_values ) / sizeof( IndexType );

    SCAI_ASSERT_EQ_ERROR( numValues, ia_values[numRows], "size mismatch" )

    csrIA.setRawData( numRows + 1, ia_values );
    csrJA.setRawData( numValues, ja_values );
    csrValues.setRawData( numValues, nz_values );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getCSCTestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    hmemo::HArray<IndexType>& cscIA,
    hmemo::HArray<IndexType>& cscJA,
    hmemo::HArray<ValueType>& cscValues )
{
    /*   Matrix:       6  0  0  4         0    3   2   0        6  5  9  4
                       7  0  0  0         1    6       2        7  1     4
                       0  0  9  4         3            3        2        3
                       2  5  0  3         4            4        2        1
                       2  0  0  1                      6                 2
                       0  0  0  0
                       0  1  0  2
     */

    const IndexType ja_values[]  = { 0,          4,    6, 7,             12 };
    const IndexType ia_values[]  = { 0, 1, 3, 4, 3, 6, 2, 0, 2, 3, 4, 6 };
    const ValueType nz_values[]  = { 6, 7, 2, 2, 5, 1, 9, 4, 4, 3, 1, 2 };

    numRows    = 7;
    numColumns = sizeof( ja_values ) / sizeof( IndexType ) - 1;
    numValues  = sizeof( ia_values ) / sizeof( IndexType );

    SCAI_ASSERT_EQ_ERROR( numValues, ja_values[numColumns], "size mismatch" )

    cscJA.setRawData( numColumns + 1, ja_values );
    cscIA.setRawData( numValues, ia_values );
    cscValues.setRawData( numValues, nz_values );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getELLTestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValuesPerRow,
    hmemo::HArray<IndexType>& ellIA,
    hmemo::HArray<IndexType>& ellJA,
    hmemo::HArray<ValueType>& ellValues )
{
    /*   Matrix:       6  0  0  4         0  3  -    6  4  -
                       7  0  0  0         0  -  -    7  -  -
                       0  0  9  4         2  3  -    9  4  -
                       2  5  0  3         0  1  3    2  5  3
                       2  0  0  1         0  3  -    2  1  -
                       0  0  0  0         -  -  -    -  -  -
                       0  1  0  2         1  3  -    1  2  -
     */

    const IndexType x = 0;
    const ValueType v = 0;

    const IndexType ia_sizes[]   = { 2, 1, 2, 3, 2, 0, 2 };
    const IndexType ja_values[]  = { 0, 0, 2, 0, 0, x, 1, 3, x, 3, 1, 3, x, 3, x, x, x, 3, x, x, x };
    const ValueType nz_values[]  = { 6, 7, 9, 2, 2, v, 1, 4, v, 4, 5, 1, v, 2, v, v, v, 3, v, v, v };

    numRows         = sizeof( ia_sizes ) / sizeof( IndexType );
    numColumns      = 4;

    SCAI_ASSERT_EQ_ERROR( sizeof( ja_values ) / sizeof( IndexType ),
                          sizeof( nz_values ) / sizeof( ValueType ),
                          "ja and values of ELL format must have same number of entries" )

    IndexType numValues  = sizeof( ja_values ) / sizeof( IndexType );

    numValuesPerRow = numValues / numRows;

    SCAI_ASSERT_EQ_ERROR( numValues, numValuesPerRow * numRows, "size mismatch" )

    ellIA.setRawData( numRows, ia_sizes );
    ellJA.setRawData( numValues, ja_values );
    ellValues.setRawData( numValues, nz_values );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getJDSTestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numDiagonals,
    hmemo::HArray<IndexType>& jdsPerm,
    hmemo::HArray<IndexType>& jdsILG,
    hmemo::HArray<IndexType>& jdsDLG,
    hmemo::HArray<IndexType>& jdsJA,
    hmemo::HArray<ValueType>& jdsValues )
{
    /*   Matrix:       6  0  0  4         0  3  -    6  4  -
                       7  0  0  0         0  -  -    7  -  -
                       0  0  9  4         2  3  -    9  4  -
                       2  5  0  3         0  1  3    2  5  3
                       2  0  0  1         0  3  -    2  1  -
                       0  0  0  0         -  -  -    -  -  -
                       0  1  0  2         1  3  -    1  2  -

                                   3      0  1  3    2  5  3
                                   0      0  3  -    6  4  -
                                   2      2  3  -    9  4  -
                                   4      0  3  -    2  1  -
                                   6      1  3  -    1  2  -
                                   1      0  -  -    7  -  -
                                   5      -  -  -    -  -  -
     */

    const IndexType perm_values[]  = { 3, 0, 2, 4, 6, 1, 5 };
    const IndexType ilg_values[]   = { 3, 2, 2, 2, 2, 1, 0 };
    const IndexType dlg_values[]   = { 6, 5, 1 };
    const IndexType ja_values[]    = { 0, 0, 2, 0, 1, 0, 1, 3, 3, 3, 3, 3 };
    const ValueType nz_values[]    = { 2, 6, 9, 2, 1, 7, 5, 4, 4, 1, 2, 3 };

    numRows         = sizeof( perm_values ) / sizeof( IndexType );
    numColumns      = 4;
    numDiagonals    = sizeof( dlg_values ) / sizeof( IndexType );

    IndexType numValues  = sizeof( ja_values ) / sizeof( IndexType );

    IndexType ilg_n = sizeof( ilg_values ) / sizeof( IndexType );
    SCAI_ASSERT_EQ_ERROR( ilg_n, numRows, "ilg and perm must have same number of entries" );

    IndexType nz_n = sizeof( nz_values ) / sizeof( ValueType );
    SCAI_ASSERT_EQ_ERROR( nz_n, numValues, "ja and values must have same number of entries" );

    jdsPerm.setRawData( numRows, perm_values );
    jdsILG.setRawData( numRows, ilg_values );
    jdsDLG.setRawData( numDiagonals, dlg_values );
    jdsJA.setRawData( numValues, ja_values );
    jdsValues.setRawData( numValues, nz_values );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getDIATestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numDiagonals,
    hmemo::HArray<IndexType>& diaOffsets,
    hmemo::HArray<ValueType>& diaValues )
{
    /*                -5 -4 -3 -2 -1  0  1  2  3

        Matrix:     x  x  x  x  x  x  6  0  0  4
                       x  x  x  x  x  7  0  0  0  x
                          x  x  x  x  0  0  9  4  x  x
                             x  x  x  2  5  0  3  x  x  x
                                x  x  2  0  0  1  x  x  x  x
                                   x  0  0  0  0  x  x  x  x  x
                                      0  1  0  2  x  x  x  x  x  x  */

    const IndexType diag_offsets[] = { static_cast<IndexType>( -5 ),
                                       static_cast<IndexType>( -4 ),
                                       static_cast<IndexType>( -3 ),
                                       static_cast<IndexType>( -2 ),
                                       static_cast<IndexType>( -1 ),
                                       0, 1, 3
                                     };

    const ValueType x = 0;  // just a stupid value as these entries should never be used

    const ValueType diag_values[]  = { x, x, x, x, x, 0, 1,
                                       x, x, x, x, 2, 0, 0,
                                       x, x, x, 2, 0, 0, 2,
                                       x, x, 0, 5, 0, 0, x,
                                       x, 7, 0, 0, 1, x, x,
                                       6, 0, 9, 3, x, x, x,
                                       0, 0, 4, x, x, x, x,
                                       4, x, x, x, x, x, x
                                     };

    numRows      = 7;
    numColumns   = 4;
    numDiagonals = sizeof( diag_offsets ) / sizeof( IndexType );

    const IndexType diag_nvalues = sizeof( diag_values ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( diag_nvalues, numRows * numDiagonals );

    diaOffsets.setRawData( numDiagonals, diag_offsets );
    diaValues.setRawData( diag_nvalues, diag_values );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getCOOTestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    hmemo::HArray<IndexType>& cooIA,
    hmemo::HArray<IndexType>& cooJA,
    hmemo::HArray<ValueType>& cooValues )
{
    /*   Matrix:       6  0  0  4
                       7  0  0  0
                       0  0  9  4
                       2  5  0  3
                       2  0  0  1
                       0  0  0  0
                       0  1  0  2
     */

    const IndexType ia_values[]  = { 0, 0, 1, 2, 2, 3, 3, 3, 4, 4, 6, 6 };
    const IndexType ja_values[]  = { 0, 3, 0, 2, 3, 0, 1, 3, 0, 3, 1, 3 };
    const ValueType nz_values[]  = { 6, 4, 7, 9, 4, 2, 5, 3, 2, 1, 1, 2 };

    numRows     = 7;
    numColumns  = 4;
    numValues   = sizeof( nz_values ) / sizeof( ValueType );

    cooIA.setRawData( numValues, ia_values );
    cooJA.setRawData( numValues, ja_values );
    cooValues.setRawData( numValues, nz_values );
}

/* ------------------------------------------------------------------------------------- */

/** This method computes the result of a matrix-vector multiplication via dense data by hand */

template<typename ValueType>
static hmemo::HArray<ValueType> getGEMVNormalResult(
    const ValueType alpha,
    const hmemo::HArray<ValueType>& x,
    const ValueType beta,
    const hmemo::HArray<ValueType>& y )
{
    IndexType numRows;
    IndexType numColumns;

    hmemo::HArray<ValueType> denseValues;

    getDenseTestData( numRows, numColumns, denseValues );

    SCAI_ASSERT_EQ_ERROR( x.size(), numColumns, "size mismatch for x" )
    SCAI_ASSERT_EQ_ERROR( y.size(), numRows, "size mismatch for y" )

    hmemo::HArray<ValueType> res;

    {
        using hmemo::ReadAccess;
        using hmemo::WriteOnlyAccess;

        ReadAccess<ValueType> rX( x );
        ReadAccess<ValueType> rY( y );
        WriteOnlyAccess<ValueType> wRes( res, numRows );
        ReadAccess<ValueType> rDense( denseValues );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType v = 0;

            for ( IndexType j = 0; j < numColumns; ++j )
            {
                v += rDense[ i * numColumns + j ] * rX[j];
            }

            wRes[i] = alpha * v + beta * rY[i];
        }
    }

    return res;
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getReduceResult( hmemo::HArray<ValueType>& res,
                             IndexType dim )
{
    using hmemo::ReadAccess;
    using hmemo::WriteOnlyAccess;

    IndexType numRows;
    IndexType numColumns;

    hmemo::HArray<ValueType> denseValues;

    getDenseTestData( numRows, numColumns, denseValues );

    if ( dim == 0 )
    {
        WriteOnlyAccess<ValueType> wRes( res, numRows );
        ReadAccess<ValueType> rDense( denseValues );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType v = 0;

            for ( IndexType j = 0; j < numColumns; ++j )
            {
                ValueType elem = rDense[ i * numColumns + j ];
                v += elem;
            }

            wRes[i] = v;
        }
    }
    else if ( dim == 1 )
    {
        WriteOnlyAccess<ValueType> wRes( res, numColumns );
        ReadAccess<ValueType> rDense( denseValues );

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            ValueType v = 0;

            for ( IndexType i = 0; i < numRows; ++i )
            {
                ValueType elem = rDense[ i * numColumns + j ];
                v += elem;
            }

            wRes[j] = v;
        }
    }
}

/* ------------------------------------------------------------------------------------- */

/** This method computes the result of a vector-matrix multiplication via dense data by hand */

template<typename ValueType>
static hmemo::HArray<ValueType> getGEMVTransposeResult( 
    const ValueType alpha,
    const hmemo::HArray<ValueType>& x,
    const ValueType beta,
    const hmemo::HArray<ValueType>& y )
{
    IndexType numRows;
    IndexType numColumns;

    hmemo::HArray<ValueType> denseValues;

    getDenseTestData( numRows, numColumns, denseValues );

    SCAI_ASSERT_EQ_ERROR( x.size(), numRows, "size mismatch for x" )
    SCAI_ASSERT_EQ_ERROR( y.size(), numColumns, "size mismatch for y" )

    hmemo::HArray<ValueType> res;

    {
        using hmemo::ReadAccess;
        using hmemo::WriteOnlyAccess;

        ReadAccess<ValueType> rX( x );
        ReadAccess<ValueType> rY( y );
        WriteOnlyAccess<ValueType> wRes( res, numColumns );
        ReadAccess<ValueType> rDense( denseValues );

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            ValueType v = 0;

            for ( IndexType i = 0; i < numRows; ++i )
            {
                v += rX[i] * rDense[ i * numColumns + j ];
            }

            wRes[j] = alpha * v + beta * rY[j];
        }
    }

    return res;
}

/* ------------------------------------------------------------------------------------- */

/** This method computes the result of a jacobi step by hand */

template<typename ValueType>
static void getJacobiHaloResult( hmemo::HArray<ValueType>& solution,
                                 const hmemo::HArray<ValueType>& oldSolution,
                                 const hmemo::HArray<ValueType>& diag,
                                 const ValueType omega )
{
    IndexType numRows;
    IndexType numColumns;

    hmemo::HArray<ValueType> denseValues;

    getDenseTestData( numRows, numColumns, denseValues );

    SCAI_ASSERT_EQ_ERROR( solution.size(), numRows, "size mismatch for solution" )
    SCAI_ASSERT_EQ_ERROR( oldSolution.size(), numRows, "size mismatch for oldSolution" )
    SCAI_ASSERT_EQ_ERROR( diag.size(), numRows, "size mismatch for diagonal" )

    {
        using hmemo::ReadAccess;
        using hmemo::WriteAccess;

        ReadAccess<ValueType> rOld( oldSolution );
        ReadAccess<ValueType> rDiag( diag );
        WriteAccess<ValueType> wSol( solution );
        ReadAccess<ValueType> rDense( denseValues );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType v    = 0;

            for ( IndexType j = 0; j < numColumns; ++j )
            {
                v += rDense[ i * numColumns + j ] * rOld[j];
            }

            wSol[i] -= omega * v / rDiag[i];
        }
    }
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static ValueType getMaxVal()
{
    return ValueType( 9 );
}

/* ------------------------------------------------------------------------------------- */

}

}
