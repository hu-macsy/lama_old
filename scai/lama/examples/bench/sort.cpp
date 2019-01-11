/**
 * @file sort.cpp
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
 * @brief Benchmark program to measure runtime of (parallel) sort
 * @author Thomas Brandes
 * @date 09.11.2016
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

// _Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/mepr/TypeList.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;
using namespace dmemo;
using namespace std;

using common::Walltime;

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
static void bench( const IndexType N )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    PartitionId rank = comm->getRank();

    DistributionPtr blockDist( new BlockDistribution( N, comm ) );

    // generate random numbers

    DenseVector<ValueType> X;
    DenseVector<IndexType> perm;

    srand( 131 + comm->getRank() );

    X.setRandom( blockDist, 1 );

    DenseVector<ValueType> Xsave( X );  // save for comparison

    bool debug = false;

    common::Settings::getEnvironment( debug, "SCAI_DEBUG" );

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
        const hmemo::HArray<ValueType>& localValues = X.getLocalValues();
        const hmemo::HArray<IndexType>& permValues = perm.getLocalValues();

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

    Xsave.gatherFrom( Xcomp, perm );
    Xcomp -= X;
    RealType<ValueType> maxDiff = Xcomp.maxNorm();

    tmpTime = Walltime::get() - tmpTime;

    SCAI_ASSERT_EQUAL( 0, maxDiff, "wrong sort result (gather test)" )

    HOST_PRINT( rank, "Gather/compare time: " << tmpTime << " seconds" )

    // equivalent check by testing Xcmp[perm] = X, Xcmp == Xsave

    tmpTime = Walltime::get();

    // in contrary to gather the scatter method requires an allocated array

    Xcomp.allocate( Xsave.getDistributionPtr() );
    Xcomp.scatter( perm, true, X );

    if ( debug )
    {
        const hmemo::HArray<ValueType>& localValues1 = Xcomp.getLocalValues();
        const hmemo::HArray<ValueType>& localValues2 = Xsave.getLocalValues();

        SCAI_ASSERT_EQUAL( localValues1.size(), localValues2.size(), "serious mismatch" )

        for ( IndexType i = 0; i < localValues1.size(); ++i )
        {
            cout << "Xcomp[local:" << i << "] = " << localValues1[i] << endl;
            cout << "Xsave[local:" << i << "] = " << localValues2[i] << endl;
        }
    }

    Xcomp -= Xsave;
    maxDiff = Xcomp.maxNorm();

    tmpTime = Walltime::get() - tmpTime;

    // check for sorted values

    SCAI_ASSERT_EQUAL( 0, maxDiff, "wrong sort result (scatter test)" )

    HOST_PRINT( rank, "Scatter/compare time: " << tmpTime << " seconds" )
}

/* ----------------------------------------------------------------------------- */

template<typename TList> struct Calling;

// termination call

template<> struct Calling<common::mepr::NullType>
{
    static bool callBench( const common::ScalarType, const IndexType )
    {
        return false;
    }
};

// call bench for header T and recursive call for tail of list

template<typename HeadType, typename TailTypes>
struct Calling<common::mepr::TypeList<HeadType, TailTypes> >
{
    static bool callBench( const common::ScalarType stype, const IndexType n )
    {
        if ( common::TypeTraits<HeadType>::stype == stype )
        {
            bench<HeadType>( n );
            return true;
        }
        else
        {
            return Calling<TailTypes>::callBench( stype, n );
        }
    }
};

/* ----------------------------------------------------------------------------- */

int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.sort" )

    common::Settings::parseArgs( argc, argv );

    // evaluate SCAI_NUM_THREADS

    int nThreads;

    if ( common::Settings::getEnvironment( nThreads, "SCAI_NUM_THREADS" ) )
    {
        omp_set_num_threads( nThreads );
    }

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    PartitionId rank = comm->getRank();

    // evaluate SCAI_TYPE

    string typeString;

    common::ScalarType dataType = common::ScalarType::DOUBLE;

    if ( common::Settings::getEnvironment( typeString, "SCAI_TYPE" ) )
    {
        dataType = common::str2ScalarType( typeString.c_str() );

        if ( dataType == common::ScalarType::UNKNOWN )
        {
            HOST_PRINT( rank, "SCAI_TYPE=" << typeString << ": is not a known data type" )
            return -1;
        }
    }

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

    if ( Calling<SCAI_ARRAY_TYPES_HOST_LIST>::callBench( dataType, n ) )
    {
        HOST_PRINT( rank, "sort bench for dataType " << dataType << " completed" )
        return 0;
    }
    else
    {
        HOST_PRINT( rank, "sort for dataType " << dataType << " unsupported" )
        return -1;
    }
}
