/**
 * @file gridFillVector.hpp
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
 * @brief Example program of filling a distributed vector
 * @author Thomas Brandes
 * @date 29.06.2018
 */

#include <scai/lama.hpp>

#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GridDistribution.hpp>

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridWriteAccess.hpp>

using namespace scai;
using namespace lama;

using hmemo::HArray;

const IndexType N1 = 100;
const IndexType N2 = 100;
const IndexType N3 = 100;

typedef double ValueType;

/** This method computes the entry for the vector of size N at (global) position i ) */

ValueType getValue( 
    const IndexType i1, const IndexType N1,
    const IndexType i2, const IndexType N2, 
    const IndexType i3, const IndexType N3 )
{
    return   sqrt( ValueType( i1 ) / ValueType( N1 ) )
           + sin( ValueType( i2 ) / ValueType( N2 ) )
           + cos( ValueType( i3 ) / ValueType( N3 ) );
}

/** This method computes the full vector as a single array */

void readData( HArray<ValueType>& data, const IndexType N1, const IndexType N2, const IndexType N3 )
{
    auto write = hmemo::hostWriteOnlyAccess( data, N1 * N2 * N3 );

    for ( IndexType i1 = 0; i1 < N1; ++i1 )
    {
        for ( IndexType i2 = 0; i2 < N2; ++i2 )
        {
             for ( IndexType i3 = 0; i3 < N3; ++i3 )
             {
                 write[i1 * N2 * N3 + i2 * N3 + i3 ] = getValue( i1, N1, i2, N2, i3, N3 );
             }
        }
    }
}

/** Main program runs different versions of dense vector setup */

int main( )
{
    // Take the default communicator (usually MPI_WORLD)

    auto comm = dmemo::Communicator::getCommunicatorPtr();

    // Set up the grid distribution

    common::Grid3D grid( N1, N2, N3 );

    auto dist = std::make_shared<dmemo::GridDistribution>( grid, comm );

    DenseVector<ValueType> v;

    // Solution 1 : replicated fill up

    {
        SCAI_REGION( "main.repFillUp" )
        
        HArray<ValueType> data;

        readData( data, N1, N2, N3 );   // all processors read the data 

        v = distribute<DenseVector<ValueType>>( data, dist );

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 1: " << v.sum() << std::endl;

    // Solution 2 : host process fills up

    {
        SCAI_REGION( "main.hostFillUp" );

        PartitionId host = 0;

        auto singleDist = std::make_shared<dmemo::SingleDistribution>( N1 * N2 * N3, comm, host );

        HArray<ValueType> data;

        if ( comm->getRank() == host )
        {
            readData( data, N1, N2, N3 );
        }
 
        v = DenseVector<ValueType>( singleDist, data );

        v.redistribute( dist );

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 2: " << v.sum() << std::endl;

    // Solution 3 : elementwise fill up

    {
        SCAI_REGION( "main.singleFillUp" )

        v = DenseVector<ValueType>( dist, ValueType( 0 ) );
        
        for ( IndexType i1 = 0; i1 < N1; ++i1 )
        {
            for ( IndexType i2 = 0; i2 < N2; ++i2 )
            {
                 for ( IndexType i3 = 0; i3 < N3; ++i3 )
                 {
                     v[i1 * N2 * N3 + i2 * N3 + i3 ] = getValue( i1, N1, i2, N2, i3, N3 );
                 }
            }
        }

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 3: " << v.sum() << std::endl;

    // Solution 4 : assembly set up

    {
        SCAI_REGION( "main.assembly" );

        VectorAssembly<ValueType> assembly;

        IndexType rank = comm->getRank();
        IndexType size = comm->getSize();

        // It does not matter which processor assembles which value

        for ( IndexType i1 = rank; i1 < N1; i1++ )
        {
            for ( IndexType i2 = 0; i2 < N2; ++i2 )
            {
                 for ( IndexType i3 = 0; i3 < N3; ++i3 )
                 {
                     ValueType v = getValue( i1, N1, i2, N2, i3, N3 );
                     IndexType globalIndex = i1 * N2 * N3 + i2 * N3 + i3 ;
                     assembly.push( globalIndex, v );
                 }
            }
        }

        v = DenseVector<ValueType>( dist, ValueType( 0 ) );
        v.fillFromAssembly( assembly );

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 4: " << v.sum() << std::endl;

    // Solution 5 : local fill up

    {
        SCAI_REGION( "main.localFillUp" )

        HArray<IndexType> myGlobalIndexes;

        dist->getOwnedIndexes( myGlobalIndexes );

        const IndexType localN = myGlobalIndexes.size();

        HArray<ValueType> localData;

        {
            auto readIndexes = hmemo::hostReadAccess( myGlobalIndexes );
            auto localWrite = hmemo::hostWriteOnlyAccess( localData, localN );

            for ( IndexType i = 0; i < localN; ++i )
            {
                IndexType pos[3];

                // compute the 'global' grid position

                grid.gridPos( pos, readIndexes[i] );

                localWrite[i] = getValue( pos[0], N1, pos[1], N2, pos[2], N3 );
            }
        }
 
        v = DenseVector<ValueType>( dist, localData );

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 5: " << v.sum() << std::endl;

    // Solution 6 : gridWriteAccess

    {
        SCAI_REGION( "main.localGridFillUp" )

        GridVector<ValueType> gv( dist );

        {
            const IndexType* lb = dist->localLB();   // offsets for local->global

            GridWriteAccess<ValueType> localWrite( gv );

            const IndexType localN1 = localWrite.size( 0 );
            const IndexType localN2 = localWrite.size( 1 );
            const IndexType localN3 = localWrite.size( 2 );

            for ( IndexType i1 = 0; i1 < localN1; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < localN2; ++i2 )
                {
                     for ( IndexType i3 = 0; i3 < localN3; ++i3 )
                     {
                         localWrite( i1, i2, i3 ) = getValue( i1 + lb[0], N1, i2 + lb[1], N2, i3 + lb[2], N3 );
                     }
                }
            }
        }

        v = gv;

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 6: " << v.sum() << std::endl;
}
