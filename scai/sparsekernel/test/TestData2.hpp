/**
 * @file TestData2.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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

namespace data2
{

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getDenseTestData(
    IndexType& numRows,
    IndexType& numColumns,
    hmemo::HArray<ValueType>& denseValues )
{
    /*   Matrix:       6  1  0  4
                      -2  8  3  0
                       0  1  9  4
                      -1  0  1  3
     */

    denseValues = { 6, 1, 0, 4,
                   -2, 8, 3, 0,
                    0, 1, 9, 4,
                   -1, 0, 1, 3   };

    numRows    = 4;
    numColumns = 4;

    IndexType numValues  = denseValues.size();

    SCAI_ASSERT_EQ_ERROR( numValues, numRows * numColumns, "size mismatch" )
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
    /*   Matrix:       6  1  0  4        0  1  3    6  1  4
                      -2  8  3  0        1  0  2    8 -2  3
                       0  1  9  4        2  1  3    9  1  4
                      -1  0  1  3        3  0  2    3 -1  1

         Note: CSR data is sorted here, but is not mandatory
     */

    csrIA     = { 0,          3,          6,       9,       12 };
    csrJA     = { 0,  1,  3,  0,  1,  2,  1, 2, 3, 0,  2, 3 };
    csrValues = { 6,  1,  4,  -2, 8,  3,  1, 9, 4, -1, 1, 3 };

    numRows    = csrIA.size() - 1;
    numColumns = 4;
    numValues  = csrIA[numRows];

    SCAI_ASSERT_EQ_ERROR( numValues, csrJA.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( numValues, csrValues.size(), "size mismatch" )
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
    /*   Matrix:       6  1  0  4        0  1  3    6  1  4
                      -2  8  3  0        0  1  2    -2 8  3
                       0  1  9  4        1  2  3    1  9  4
                      -1  0  1  3        0  2  3   -1  1  3
     */

    ellIA     = { 3, 3, 3, 3 };
    ellJA     = { 0, 0, 1, 0, 1,  1, 2,  2, 3, 2, 3, 3 };
    ellValues = { 6, -2, 1, -1, 1, 8, 9, 1, 4, 3, 4, 3 };

    numRows         = ellIA.size();
    numColumns      = 4;

    SCAI_ASSERT_EQ_ERROR( ellJA.size(), ellValues.size(), "serious size mismatch for ELL data" )

    IndexType ellSize = ellJA.size();

    numValuesPerRow = ellSize / numRows;

    SCAI_ASSERT_EQ_ERROR( ellSize, numValuesPerRow * numRows, "size of ELL arrays not multiple of numRows = " << numRows )
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
    /*   Matrix:       6  1  0  4        0  1  3    6  1  4
                      -2  8  3  0        1  0  2    8 -2  3
                       0  1  9  4        1  2  3    1  9  4
                      -1  0  1  3        3  0  2    3 -1  1
     */

    const IndexType perm_values[] = { 0, 1, 2, 3 };
    const IndexType ilg_values[]  = { 3, 3, 3, 3 };
    const IndexType dlg_values[]  = { 4, 4, 4 };
    const IndexType ja_values[]   = { 0, 1, 1, 3, 1,  0, 2,  0, 3, 2, 3, 2 };
    const ValueType nz_values[]   = { 6, 8, 1, 3, 1, -2, 9, -1, 4, 3, 4, 1 };

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
    /*   Matrix:
                                            -3    -1  0  1     3
                      6  1  0  4             x  x  x  6  1  0  4
                     -2  8  3  0             x  x -2  8  3  0  x
                      0  1  9  4             x  0  1  9  4  x  x
                     -1  0  1  3            -1  0  1  3  x  x  x
    */

    const IndexType diag_offsets[] = { 0,
                                       static_cast<IndexType>( -3 ),
                                       static_cast<IndexType>( -1 ),
                                       1, 3
                                     };

    const ValueType x = 0;  // just a stupid value as these entries should never be used

    const ValueType diag_values[]  = { 6, 8, 9, 3,
                                       x, x, x, -1,
                                       x, -2, 1, 1,
                                       1, 3, 4, x,
                                       4, x, x, x
                                     };

    numRows      = 4;
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
    /*   Matrix:       6  1  -  4
                      -2  8  3  -
                       -  1  9  4
                      -1  -  1  3

         Note: COO data is sorted with rows as first key and columns as second key.
     */

    cooIA     = { 0, 0, 0,  1, 1, 1, 2, 2, 2,  3, 3, 3 };
    cooJA     = { 0, 1, 3,  0, 1, 2, 1, 2, 3,  0, 2, 3 };
    cooValues = { 6, 1, 4, -2, 8, 3, 1, 9, 4, -1, 1, 3 };

    SCAI_ASSERT_EQ_ERROR( cooIA.size(), cooJA.size(), "serious size mismatch on COO data" )
    SCAI_ASSERT_EQ_ERROR( cooIA.size(), cooValues.size(), "serious size mismatch on COO data" )

    numRows     = 4;
    numColumns  = 4;
    numValues   = cooIA.size();
}

/* ------------------------------------------------------------------------------------- */

/** This method computes the result of a jacobi step by hand */

template<typename ValueType>
static void getJacobiResult( hmemo::HArray<ValueType>& solution,
                             const hmemo::HArray<ValueType>& oldSolution,
                             const ValueType omega,
                             const hmemo::HArray<ValueType>& rhs )
{
    IndexType numRows;
    IndexType numColumns;

    hmemo::HArray<ValueType> denseValues;

    getDenseTestData( numRows, numColumns, denseValues );

    SCAI_ASSERT_EQ_ERROR( oldSolution.size(), numRows, "size mismatch for x" )
    SCAI_ASSERT_EQ_ERROR( rhs.size(), numRows, "size mismatch for x" )

    {
        hmemo::ReadAccess<ValueType> rOld( oldSolution );
        hmemo::ReadAccess<ValueType> rRhs( rhs );
        hmemo::WriteOnlyAccess<ValueType> wSol( solution, numRows );
        hmemo::ReadAccess<ValueType> rDense( denseValues );

        ValueType one = 1;

        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType v    = 0;
            ValueType diag = one;

            for ( IndexType j = 0; j < numColumns; ++j )
            {
                if ( i == j )
                {
                    diag = rDense[ i * numColumns + j ];
                }
                else
                {
                    v += rDense[ i * numColumns + j ] * rOld[j];
                }
            }

            wSol[i] = omega * ( rRhs[i] - v ) / diag + ( one - omega ) * rOld[i];
        }
    }
}

/* ------------------------------------------------------------------------------------- */

}

}
