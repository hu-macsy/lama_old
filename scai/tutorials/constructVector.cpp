/**
 * @file constructVector.hpp
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
 * @brief Example program of constructing dense / sparse vectors
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
    auto ctx  = hmemo::Context::getContextPtr();

    // For the following set ups this might be any distribution

    auto dist = std::make_shared<dmemo::BlockDistribution>( N, comm );

    ValueType initVal = 1;

    DenseVector<ValueType> denseVector1;
    DenseVector<ValueType> denseVector2( ctx );
    // not available: DenseVector<ValueType> denseVector3( N );
    DenseVector<ValueType> denseVector4( N, initVal );
   
    HArray<ValueType> globalData( N, initVal );
    DenseVector<ValueType> denseVector5( globalData );

    DenseVector<ValueType> denseVector6( dist, initVal );

    HArray<ValueType> localData( dist->getLocalSize(), initVal );
    DenseVector<ValueType> denseVector7( dist, std::move( localData ) );
    DenseVector<ValueType> denseVector8( denseVector7 );

    auto dv1 = fill<DenseVector<ValueType>>( dist, initVal );
    dv1.writeToFile( "data.mtx" );
    auto dv2 = eval<DenseVector<ValueType>>( 2 * dv1 );
    auto dv3 = distribute<DenseVector<ValueType>>( globalData, dist );
    auto dv4 = read<DenseVector<ValueType>>( "data.mtx" );
    auto dv5 = linearDenseVector<ValueType>( ValueType( 1 ) , ValueType( 0.1 ), N );
}
