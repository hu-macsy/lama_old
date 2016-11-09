/**
 * @file sort.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Benchmark program to measure runtime of (parallel) sort
 * @author Thomas Brandes
 * @date 09.11.2016
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/OpenMP.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;
using namespace dmemo;
using namespace std;
using scai::common::Walltime;

#define HOST_PRINT( rank, msg )   \
{                                 \
    if ( rank == 0 )              \
    {                             \
        cout << msg << endl;      \
    }                             \
}                                 \

/** Generic routine for benchmarking sort routine.
 *
 *  @tparam    ValueType type of values to be sorted
 *  @param[in] N         number of values to sort
 */
template<typename ValueType>
void bench( const IndexType N )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    PartitionId rank = comm->getRank();

    DistributionPtr blockDist( new BlockDistribution( N, comm ) );

    // generate random numbers

    DenseVector<ValueType> X;
    DenseVector<IndexType> perm;

    float fillRate = 1.0f;

    srand( 131 + comm->getRank() );

    X.setRandom( blockDist, fillRate );

    DenseVector<ValueType> Xsave( X );  // save for comparison

    bool debug = false;

    bool ascending = true;

    double tmpTime = Walltime::get();

    X.sort( ascending );

    tmpTime = Walltime::get() - tmpTime;

    HOST_PRINT( rank, "Sort time (  NO  perm) : " << tmpTime << " seconds" )

    X = Xsave;  // restore old values 

    tmpTime = Walltime::get();

    X.sort( perm, ascending );

    tmpTime = Walltime::get() - tmpTime;

    HOST_PRINT( rank, "Sort time ( with perm) : " << tmpTime << " seconds" )

    cout << *comm << ": sorted vector X = " << X << endl;

    // The following code might be used for debugging

    if ( debug )
    {
        const utilskernel::LArray<ValueType>& localValues = X.getLocalValues();
        const utilskernel::LArray<IndexType>& permValues = perm.getLocalValues();

        for ( IndexType i = 0; i < X.getDistribution().getLocalSize(); ++i )
        {
            cout << "X[local:" << i << "] = " << localValues[i] << endl;
            cout << "perm[local:" << i << "] = " << permValues[i] << endl;
        }
    }

    // check for sorted values

    bool isSorted = X.isSorted( ascending );

    SCAI_ASSERT( isSorted, "Vector not sorted correctly" )

    // check for correct permutation by testing X == Xsave[perm]

    HOST_PRINT( rank, "X.isSorted( " << ascending << " ) = " << isSorted )

    // check the sorted values

    DenseVector<ValueType> Xcomp;

    tmpTime = Walltime::get();

    Xcomp.gather( Xsave, perm );
    Xcomp -= X;
    Scalar maxDiff = Xcomp.maxNorm();

    tmpTime = Walltime::get() - tmpTime;

    HOST_PRINT( rank, "Gather/compare time: " << tmpTime << " seconds" )

    SCAI_ASSERT_EQUAL( 0, maxDiff, "wrong sort result" )
}

/* ----------------------------------------------------------------------------- */

int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.sort" )

    common::Settings::parseArgs( argc, argv );

    int nThreads;

    if ( scai::common::Settings::getEnvironment( nThreads, "SCAI_NUM_THREADS" ) )
    {
        omp_set_num_threads( nThreads );
    }

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    PartitionId rank = comm->getRank();

    if ( argc != 2 )
    {
        HOST_PRINT( rank, "Usage: " << argv[0] << " <n>" )
        HOST_PRINT( rank, " - n is the size of vector to sort" )
        return -1;
    }

    istringstream val( argv[1] );

    IndexType n = 0;

    val >> n;

    if ( val.fail() )
    {
        HOST_PRINT( rank, argv[1] << ": is not a legal size value" )
        return -1;
    }

    bench<double>( n );

    return 0;
}
