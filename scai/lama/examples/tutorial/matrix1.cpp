/**
 * @file matrix1.cpp
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
 * @brief matrix1.cpp is an example to show how getRow works
 * @author Thomas Brandes
 * @date 03.02.2016
 */

#include <scai/lama.hpp>

// _Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <iostream>
#include <memory>
#include <stdlib.h>

using namespace scai;
using namespace hmemo;
using namespace dmemo;
using namespace lama;

/** Take default real type for this example. */

typedef DefaultReal ValueType;

//
// EXAMPLE multiplication of a dense vector with a sparse matrix in CSR format.
//

static inline ValueType mv( const IndexType i, const IndexType j )
{
    return ValueType( 3 * i - j );
}

int main()
{
    IndexType perm [] = { 5, 2, 1, 0, 3, 4 };
    const IndexType irow = 3;
    const IndexType N = sizeof( perm ) / sizeof( IndexType );

    HArray<ValueType> values( N * N );

    {
        WriteOnlyAccess<ValueType> writeValues( values, N * N );

        for ( IndexType i = 0; i < N; ++i )
        {
            for ( IndexType j = 0; j < N; ++j )
            {
                writeValues[ i * N + j ] = mv( i, j );
            }
        }
    }

    DenseStorage<ValueType> denseData( N, N, std::move( values ) );

    //      0  -1 -2 -3 -4 -5
    //      3   2  1  0 -1 -2
    //      6   5  4  3  2  1
    //      9   8  7  6  5  4
    //     12  11 10  9  8  7
    //     15  14 13 12 11 10

    // convert storage -> matrix, is replicated on each processor

    auto a = convert<CSRSparseMatrix<ValueType>>( denseData );

    a.assign( denseData );   // replicated matrix

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::vector<IndexType> myGlobalIndexes;

    for ( IndexType i = 0; i < N; ++i )
    {
        if ( ( perm[i] % comm->getSize() ) == comm->getRank() )
        {
            myGlobalIndexes.push_back( i );
        }
    }

    std::cout << *comm << ": have " << myGlobalIndexes.size() << " indexes" << std::endl;

    hmemo::HArrayRef<IndexType> indexes( myGlobalIndexes );

    auto dist = std::make_shared<GeneralDistribution>( N, indexes, comm );

    a.redistribute( dist, dist );

    std::cout << "Communicator = " << *comm << std::endl;

    DenseVector<ValueType> row;   // will be initialized 

    a.getRow( row, irow );

    std::cout << "a( " << irow << ", : ) = " << row << std::endl;

    ReadAccess<ValueType> rowRead( row.getLocalValues() );

    int errors = 0;
    std::cout << "Values = ";

    for ( IndexType j = 0; j < rowRead.size(); ++ j )
    {
        std::cout << " " << rowRead[j];

        if ( rowRead[j] != mv( irow, dist->local2global( j ) ) )
        {
            std::cout << " Error";
            errors++;
        }
    }

    std::cout << std::endl;
    std::cout << "Errors = " << errors << std::endl;
    a.writeToFile( "MatrixA.mtx" );
}
