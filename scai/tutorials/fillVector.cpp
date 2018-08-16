/**
 * @file fillVector.hpp
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

using namespace scai;
using namespace lama;

using hmemo::HArray;

const IndexType N = 1000000;

typedef double ValueType;

/** This method computes the entry for the vector of size N at (global) position i ) */

ValueType getValue( const IndexType i, const IndexType N )
{
    return sqrt( ValueType( i ) / ValueType( N ) );
}

/** This method computes the full vector as a single array */

void readData( HArray<ValueType>& data, const IndexType N )
{
    auto write = hmemo::hostWriteOnlyAccess( data, N );

    for ( IndexType i = 0; i < N; ++i )
    {
        write[i] = getValue( i, N );
    }
}

/** Main program runs different versions of dense vector setup */

int main( )
{
    // Take the default communicator (usually MPI_WORLD)

    auto comm = dmemo::Communicator::getCommunicatorPtr();

    // For the following set ups this might be any distribution

    auto dist = std::make_shared<dmemo::BlockDistribution>( N, comm );

    DenseVector<ValueType> v;

    // Solution 1 : replicated fill up
    //     + most efficient, no communication
    //     - full data required on each processor

    {
        SCAI_REGION( "main.repFillUp" )
        
        HArray<ValueType> data;

        readData( data, N );   // all processors read the data 

        v = distribute<DenseVector<ValueType>>( data, dist );

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 1: " << v.sum() << std::endl;

    // Solution 2 : host process fills up
    //      +  can be used if only host processor can compute the data
    //      -  communication required to distribute the data

    {
        SCAI_REGION( "main.hostFillUp" );

        PartitionId host = 0;

        // define a distribution where all data is owned by the first processor

        auto singleDist = std::make_shared<dmemo::SingleDistribution>( N, comm, host );

        HArray<ValueType> data;

        if ( comm->getRank() == host )
        {
            readData( data, N );
        }
 
        // Define a vector where only the first processor owns all data

        v = DenseVector<ValueType>( singleDist, std::move( data ) );

        v.redistribute( dist );

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 2: " << v.sum() << std::endl;

    // Solution 3 : elementwise fill up
    //      +  most convenient use, only owner process sets data
    //      -  elementwise access is rather expensive

    {
        SCAI_REGION( "main.singleFillUp" )

        v = DenseVector<ValueType>( dist, ValueType( 0 ) );
        
        for ( IndexType i = 0; i < N; ++i )
        {
            v[i] = getValue( i, N );
            // same as v.setValue( i, getValue( i, N ) );
        }

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 3: " << v.sum() << std::endl;

    // Solution 4 : assembly set up
    //      +  each processor assembles some data, no full allocation required
    //      +  the more data is assembled locally the less communication is required
    //      -  redistribution of assembled data required 

    {
        SCAI_REGION( "main.assembly" );

        VectorAssembly<ValueType> assembly;

        IndexType rank = comm->getRank();
        IndexType size = comm->getSize();

        // It does not matter which processor assembles which value

        for ( IndexType i = rank; i < N; i += size )
        {
            assembly.push( i, getValue( i, N ) );
        }

        v = DenseVector<ValueType>( dist, ValueType( 0 ) );
        v.fillFromAssembly( assembly );

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 4: " << v.sum() << std::endl;

    // Solution 5 : local fill up
    //      +  each processor sets only its owned values
    //      -  array of owned global indexes is needed

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
                localWrite[i] = getValue( readIndexes[i], N );
            }
        }
 
        v = DenseVector<ValueType>( dist, localData );

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 5: " << v.sum() << std::endl;

    // Solution 6 : local fill up, take advantage of block distribution
    //      +  each processor sets only its owned values
    //      +  owned global indexes are is just a range of values
    //      -  works only for block distribution

    {
        SCAI_REGION( "main.localBlockFillUp" )

        const IndexType localN = dist->getLocalSize();

        HArray<ValueType> localData;

        {
            auto localWrite = hmemo::hostWriteOnlyAccess( localData, localN );

            IndexType lb = dist->lb();

            for ( IndexType i = 0; i < localN; ++i )
            {
                localWrite[i] = getValue( lb + i, N );
            }
        }

        v = DenseVector<ValueType>( dist, localData );

        comm->synchronize();  // just for timing
    }

    std::cout << "Result 6: " << v.sum() << std::endl;
}
