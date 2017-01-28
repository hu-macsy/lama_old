/**
 * @file TestData2.hpp
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

using scai::hmemo::HArray;

namespace data2
{

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getDenseTestData(
    IndexType& numRows,
    IndexType& numColumns,
    HArray<ValueType>& denseValues )
{
    /*   Matrix:       6  1  0  4
                      -2  8  3  0
                       0  1  9  4
                      -1  0  1  3
     */

    const ValueType values[]  = { 6, 1, 0, 4,
                                  -2, 8, 3, 0,
                                  0, 1, 9, 4,
                                  -1, 0, 1, 3
                                };

    numRows    = 4;
    numColumns = 4;

    IndexType numValues  = sizeof( values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numValues, numRows * numColumns, "size mismatch" )

    denseValues.init( values, numValues );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getCSRTestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues )
{
    /*   Matrix:       6  1  0  4        0  1  3    6  1  4
                      -2  8  3  0        1  0  2    8 -2  3
                       0  1  9  4        2  1  3    9  1  4
                      -1  0  1  3        3  0  2    3 -1  1

         Note: diagonals are always the first entry
     */

    const IndexType ia_values[]  = { 0,          3,          6,       9,       12 };
    const IndexType ja_values[]  = { 0,  1,  3,  1,  0,  2,  2, 1, 3, 3,  0, 2 };
    const ValueType nz_values[]  = { 6,  1,  4,  8, -2,  3,  9, 1, 4, 3, -1, 1 };

    numRows    = sizeof( ia_values ) / sizeof( IndexType ) - 1;
    numColumns = 4;
    numValues  = sizeof( ja_values ) / sizeof( IndexType );

    SCAI_ASSERT_EQ_ERROR( numValues, ia_values[numRows], "size mismatch" )

    csrIA.init( ia_values, numRows + 1 );
    csrJA.init( ja_values, numValues );
    csrValues.init( nz_values, numValues );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getELLTestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValuesPerRow,
    HArray<IndexType>& ellIA,
    HArray<IndexType>& ellJA,
    HArray<ValueType>& ellValues )
{
    /*   Matrix:       6  1  0  4        0  1  3    6  1  4
                      -2  8  3  0        1  0  2    8 -2  3
                       0  1  9  4        2  1  3    9  1  4
                      -1  0  1  3        3  0  2    3 -1  1
     */

    const IndexType ia_sizes[]   = { 3, 3, 3, 3 };
    const IndexType ja_values[]  = { 0, 1, 2, 3, 1,  0, 1,  0, 3, 2, 3, 2 };
    const ValueType nz_values[]  = { 6, 8, 9, 3, 1, -2, 1, -1, 4, 3, 4, 1 };

    numRows         = sizeof( ia_sizes ) / sizeof( IndexType );
    numColumns      = 4;

    IndexType numValues  = sizeof( ja_values ) / sizeof( IndexType );

    numValuesPerRow = numValues / numRows;

    SCAI_ASSERT_EQ_ERROR( numValues, numValuesPerRow * numRows, "size mismatch" )

    ellIA.init( ia_sizes, numRows );
    ellJA.init( ja_values, numValues );
    ellValues.init( nz_values, numValues );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getJDSTestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numDiagonals,
    HArray<IndexType>& jdsPerm,
    HArray<IndexType>& jdsILG,
    HArray<IndexType>& jdsDLG,
    HArray<IndexType>& jdsJA,
    HArray<ValueType>& jdsValues )
{
    /*   Matrix:       6  1  0  4        0  1  3    6  1  4
                      -2  8  3  0        1  0  2    8 -2  3
                       0  1  9  4        2  1  3    9  1  4
                      -1  0  1  3        3  0  2    3 -1  1
     */

    const IndexType perm_values[] = { 0, 1, 2, 3 };
    const IndexType ilg_values[]  = { 3, 3, 3, 3 };
    const IndexType dlg_values[]  = { 4, 4, 4 };
    const IndexType ja_values[]   = { 0, 1, 2, 3, 1,  0, 1,  0, 3, 2, 3, 2 };
    const ValueType nz_values[]   = { 6, 8, 9, 3, 1, -2, 1, -1, 4, 3, 4, 1 };

    numRows         = sizeof( perm_values ) / sizeof( IndexType );
    numColumns      = 4;
    numDiagonals    = sizeof( dlg_values ) / sizeof( IndexType );

    IndexType numValues  = sizeof( ja_values ) / sizeof( IndexType );

    IndexType ilg_n = sizeof( ilg_values ) / sizeof( IndexType );
    SCAI_ASSERT_EQ_ERROR( ilg_n, numRows, "ilg and perm must have same number of entries" );

    IndexType nz_n = sizeof( nz_values ) / sizeof( ValueType );
    SCAI_ASSERT_EQ_ERROR( nz_n, numValues, "ja and values must have same number of entries" );

    jdsPerm.init( perm_values, numRows );
    jdsILG.init( ilg_values, numRows );
    jdsDLG.init( dlg_values, numDiagonals );
    jdsJA.init( ja_values, numValues );
    jdsValues.init( nz_values, numValues );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getDIATestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numDiagonals,
    HArray<IndexType>& diaOffsets,
    HArray<ValueType>& diaValues )
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

    diaOffsets.init( diag_offsets, numDiagonals );
    diaValues.init( diag_values, diag_nvalues );
}

/* ------------------------------------------------------------------------------------- */

template<typename ValueType>
static void getCOOTestData(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValues,
    HArray<IndexType>& cooIA,
    HArray<IndexType>& cooJA,
    HArray<ValueType>& cooValues )
{
    /*   Matrix:       6  1  -  4
                      -2  8  3  -
                       -  1  9  4
                      -1  -  1  3
     */

    const IndexType ia_values[]  = { 0, 1, 2, 3, 0, 0,  1, 1, 2,  2, 3, 3 };
    const IndexType ja_values[]  = { 0, 1, 2, 3, 1, 3,  0, 2, 1,  3, 0, 2 };
    const ValueType nz_values[]  = { 6, 8, 9, 3, 1, 4, -2, 3, 1,  4, -1, 1 };

    numRows     = 4;
    numColumns  = 4;
    numValues   = sizeof( nz_values ) / sizeof( ValueType );

    cooIA.init( ia_values, numValues );
    cooJA.init( ja_values, numValues );
    cooValues.init( nz_values, numValues );
}

/* ------------------------------------------------------------------------------------- */

/** This method computes the result of a jacobi step by hand */

template<typename ValueType>
static void getJacobiResult( HArray<ValueType>& solution,
                             const HArray<ValueType>& oldSolution,
                             const ValueType omega,
                             const HArray<ValueType>& rhs )
{
    IndexType numRows;
    IndexType numColumns;

    HArray<ValueType> denseValues;

    getDenseTestData( numRows, numColumns, denseValues );

    SCAI_ASSERT_EQ_ERROR( oldSolution.size(), numRows, "size mismatch for x" )
    SCAI_ASSERT_EQ_ERROR( rhs.size(), numRows, "size mismatch for x" )

    {
        using scai::hmemo::ReadAccess;
        using scai::hmemo::WriteOnlyAccess;

        ReadAccess<ValueType> rOld( oldSolution );
        ReadAccess<ValueType> rRhs( rhs );
        WriteOnlyAccess<ValueType> wSol( solution, numRows );
        ReadAccess<ValueType> rDense( denseValues );

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
