/**
 * @file resortMatrix.cpp
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
 * @brief ToDo: Missing description in ./partitioning/examples/resortMatrix.cpp
 * @author Thomas Brandes
 * @date 31.08.2017
 */

#include <scai/lama.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/lama/io/FileIO.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <memory>

using namespace scai;
using namespace hmemo;
using namespace utilskernel;
using namespace lama;

typedef float real_t;

typedef DefaultReal ValueType;

int main( int narg, const char* argv[] )
{
    if ( narg < 5 )
    {
        std::cout << "Please call: " << argv[0] << " <inMatrixFileName> <outMatrixFileName> <distRowFile> <distColFile>" << std::endl;
        return -1;
    }

    std::string inFileName = argv[1];
    std::string outFileName = argv[2];
    std::string rowDistFileName = argv[3];
    std::string colDistFileName = argv[4];

    // read in the CSR storage

    CSRStorage<ValueType> storage;

    storage.readFromFile( inFileName );

    std::cout << "Read matrix " << storage << std::endl;

    IndexType numRows    = storage.getNumRows();
    IndexType numColumns = storage.getNumColumns();
    IndexType numValues  = storage.getNumValues();

    HArray<IndexType> rowDist;

    FileIO::read( rowDist, rowDistFileName );

    SCAI_ASSERT_EQ_ERROR( rowDist.size(), numRows, "array length of distRow does not match #rows in input matrix" )

    HArray<IndexType> colDist;

    FileIO::read( colDist, colDistFileName );

    SCAI_ASSERT_EQ_ERROR( colDist.size(), numColumns, "array length of distRow does not match #rows in input matrix" )

    IndexType np_row = HArrayUtils::max( rowDist ) + 1;
    IndexType np_col = HArrayUtils::max( colDist ) + 1;

    std::cout << "Row dist for " << np_row << " processors, col dist for " << np_col << " processors" << std::endl;

    // Now make a bucket sort of the indexes according to the distributions

    HArray<IndexType> rowPermutation;
    HArray<IndexType> colPermutation;
    HArray<IndexType> invColPermutation;

    HArray<IndexType> rowOffsets;   // offset array required as temporary, will not be used later
    HArray<IndexType> colOffsets;   // offset array required as temporary, will not be used later

    HArrayUtils::bucketSort( rowOffsets, rowPermutation, rowDist, np_row );
    HArrayUtils::bucketSort( colOffsets, colPermutation, colDist, np_col );

    std::cout << "Bucketsort rows, offsets = " << rowOffsets << ", perm = " << rowPermutation << std::endl;
    std::cout << "Bucketsort cols, offsets = " << colOffsets << ", perm = " << colPermutation << std::endl;

    // Column j in new matrix is colPermutation[j] of old matrix
    // inverse perm needed to know new column for an old column index

    HArrayUtils::inversePerm( invColPermutation, colPermutation );

    // Now build a new CSR storage

    HArray<IndexType> ia;
    HArray<IndexType> ja;
    HArray<ValueType> values;

    {
        ReadAccess<IndexType> oldIA( storage.getIA() );
        ReadAccess<IndexType> oldJA( storage.getJA() );
        ReadAccess<ValueType> oldValues( storage.getValues() );

        WriteOnlyAccess<IndexType> newIA( ia, numRows + 1 );
        WriteOnlyAccess<IndexType> newJA( ja, numValues );
        WriteOnlyAccess<ValueType> newValues( values, numValues );

        IndexType offset = 0;

        for ( IndexType i = 0; i < numRows; ++i )
        {
            newIA[i] = offset;

            IndexType oldI = rowPermutation[i];
 
            SCAI_ASSERT_VALID_INDEX_ERROR( oldI, numRows, "wrong permutation" )

            for ( IndexType jj = oldIA[ oldI ]; jj < oldIA[ oldI + 1 ]; ++jj )
            {
                newValues[ offset ] = oldValues[ jj ];

                IndexType oldJ = oldJA[jj];
                IndexType j    = invColPermutation[ oldJ];
                
                newJA[ offset ] = j;

                offset++;
            }
        }

        newIA[numRows] = offset;

        SCAI_ASSERT_EQ_ERROR( offset, numValues, "serious wrong counting" )
    }

    std::cout << "Build new csr data, ia = " << ia << ", ja = " << ja << ", values = " << values << std::endl;

    storage.swap( ia, ja, values );

    std::cout << "Write new matrix to file " << outFileName << std::endl;

    storage.writeToFile( outFileName );
}
