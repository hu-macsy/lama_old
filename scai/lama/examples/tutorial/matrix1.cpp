/**
 * @file matrix1.cpp
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
 * @brief matrix1.cpp is an example to show how getRow works
 * @author Thomas Brandes
 * @date 03.02.2016
 */

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/distribution/BlockDistribution.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <scai/common/unique_ptr.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace hmemo;
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

    CommunicatorPtr comm = Communicator::getCommunicator();

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

    DenseVector<double> row;     // any type
    row.resize( dist );          // any distribution   

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
