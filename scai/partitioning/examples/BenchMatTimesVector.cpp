/**
 * @file partitioning/examples/BenchMatTimesVector.cpp
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
 * @brief Benchmark of matrix-vector multiplication with arbitrary row/col distribution
 * @author Thomas Brandes
 * @date 03.06.2013
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

// _Matrix & vector related includes

#include <scai/lama/matrix/all.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/lama/io/PartitionIO.hpp>
#include <scai/partitioning/Partitioning.hpp>

#include <scai/tracing.hpp>

#include <memory>

using namespace scai;
using namespace lama;
using namespace dmemo;
using namespace partitioning;
using namespace std;
using scai::common::Walltime;

typedef DefaultReal ValueType;

static void bench( CSRSparseMatrix<ValueType>& mat )
{
    auto y = fillDenseVector<ValueType>( mat.getColDistributionPtr(), 1 );
    auto x = fillDenseVector<ValueType>( mat.getRowDistributionPtr(), 0 );

    x = x + mat * y;  // warm up
    x = x + mat * y;  // warm up

    mat.setCommunicationKind( SyncKind::SYNCHRONOUS );

    CommunicatorPtr comm = mat.getRowDistribution().getCommunicatorPtr();

    double time = Walltime::get();

    // measure 3 iterations
 
    for ( IndexType i = 0; i < 3; ++i )
    {
        x = x + mat * y;
    }

    time = Walltime::get() - time;
    time = time / 3;

    cout << *comm << "matrixTimesVector: " << time  << " s for 3 iterations" << std::endl;

    IndexType niter = static_cast<IndexType>( 5.0 / time ) + 1;

    comm->bcast( &niter, 1, 0 );

    // measure niter iterations
 
    time = Walltime::get();

    for ( IndexType i = 0; i < niter; ++i )
    {
        x = x + mat * y;
    }

    time = Walltime::get() - time;
    time = time / niter * 1000.0;    // time per iteration, scaled to ms

    cout << *comm << "matrixTimesVector: " << time  << " ms / Iteration, done " << niter << std::endl;
}

int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.Bench.main" )

    common::Settings::parseArgs( argc, argv );

    if ( argc < 3 )
    {
        std::cout << "Please call: " << argv[0] << " matrixFileName [ FILE <rowdist> <coldist> | BLOCK | CYCLIC | METIS | .... ] " << std::endl;
        return -1;
    }

    CSRSparseMatrix<ValueType> m;

    m.readFromFile( argv[1] );

    if ( strcmp( argv[2], "FILE" ) == 0 )
    {
        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        DistributionPtr rowDist = PartitionIO::readDistribution( argv[3], comm );
        DistributionPtr colDist = PartitionIO::readDistribution( argv[4], comm );

        m.redistribute( rowDist, colDist );
    }
    else if ( !Partitioning::canCreate( argv[2] ) )
    {
        std::cout << "ERROR: Partitioning kind = " << argv[2] << " not supported" << std::endl;

        vector<string> values;  // string is create type for the factory

        Partitioning::getCreateValues( values );

        std::cout << "Supported partitiong types are:";

        for ( size_t i = 0; i < values.size(); ++i )
        {
            std::cout << " " << values[i];
        }

        std::cout << endl;

        return -1;
    }
    else
    {
        PartitioningPtr thePartitioning( Partitioning::create( argv[2] ) );

        thePartitioning->rectangularRedistribute( m, 1.0 );
    }

    CommunicatorPtr comm = m.getRowDistribution().getCommunicatorPtr();

    std::cout << *comm << ": matrix = " << m << std::endl;

    bench( m );
}

