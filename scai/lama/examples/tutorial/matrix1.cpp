/**
 * @file matrix1.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief matrix1.cpp is an example to show how getRow works
 * @author Thomas Brandes
 * @date 03.02.2016
 */

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <scai/common/unique_ptr.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace hmemo;
using namespace dmemo;
using namespace lama;

//
// EXAMPLE multiplication of a dense vector with a sparse matrix in CSR format.
//

static inline float mv( const IndexType i, const IndexType j )
{
    return float( 3 * i - j );
}

int main()
{
    IndexType perm [] = { 5, 2, 1, 0, 3, 4 };

    const IndexType irow = 3;

    const int N = sizeof( perm ) / sizeof( IndexType );
 
    CSRSparseMatrix<float> a;

    common::scoped_array<float> values( new float[ N * N ] );

    for ( IndexType i = 0; i < N; ++i )  
    {
        for ( IndexType j = 0; j < N; ++j )  
        {
            values[ i * N + j ] = mv( i, j );
        }
    }

    DistributionPtr rep( new NoDistribution( N ) );

    a.setRawDenseData( rep, rep, values.get() );

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

    DistributionPtr dist( new GeneralDistribution( N, myGlobalIndexes, comm ) );

    a.redistribute( dist, dist );

    std::cout << "Communicator = " << *comm << std::endl;

    DenseVector<double> row( dist );     // any type, any distribution

    a.getRow( row, irow );

    std::cout << "a( " << irow << ", : ) = " << row << std::endl;

    ReadAccess<double> rowRead( row.getLocalValues() );

    int errors = 0;

    std::cout << "Values = ";

    for ( IndexType j = 0; j < rowRead.size(); ++ j )
    {
        std::cout << " " << rowRead[j];

        if ( rowRead[j] != double( mv( irow, j ) ) )
        {
            std::cout << " Error";
            errors++;
        }
    }
    std::cout << std::endl;

    std::cout << "Errors = " << errors << std::endl;

    a.writeToFile( "MatrixA", File::MATRIX_MARKET );
}
