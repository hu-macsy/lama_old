/**
 * @file scan.cpp
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
 * @brief Benchmark program to measure runtime of (parallel) scan on DenseVector
 * @author Thomas Brandes
 * @date 09.11.2016
 */

#include <iostream>
#include <iomanip>
#include <numeric>

#include <scai/lama.hpp>

// _Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

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
     
/** Own implementation of scan provided by Moritz von Looz  */

template<typename ValueType>
DenseVector<ValueType> computeGlobalPrefixSum( const DenseVector<ValueType>& input )
{
    scai::dmemo::CommunicatorPtr comm = input.getDistributionPtr()->getTargetCommunicatorPtr();

    const IndexType p = comm->getSize();

    //first, check that the input is some block distribution
    const IndexType localN = input.getDistributionPtr()->getBlockDistributionSize();

    if ( localN == invalidIndex )
    {
        throw std::logic_error( "Global prefix sum only implemented for block distribution." );
    }

    //get local prefix sum
    scai::hmemo::ReadAccess<ValueType> localValues( input.getLocalValues() );
    std::vector<ValueType> localPrefixSum( localN );
    std::partial_sum( localValues.get(), localValues.get() + localN, localPrefixSum.begin() );

    ValueType localSum[1] = {localPrefixSum[localN - 1]};

    //communicate local sums
    ValueType allOffsets[p];
    comm->gather( allOffsets, 1, 0, localSum );

    //compute prefix sum of offsets.
    std::vector<ValueType> offsetPrefixSum( p + 1, 0 );

    if ( comm->getRank() == 0 )
    {
        //shift begin of output by one, since the first offset is 0
        std::partial_sum( allOffsets, allOffsets + p, offsetPrefixSum.begin() + 1 );
    }

    //remove last value, since it would be the offset for the p+1th processor
    offsetPrefixSum.resize( p );

    //communicate offsets
    ValueType myOffset[1];
    comm->scatter( myOffset, 1, 0, offsetPrefixSum.data() );

    // get results by adding local sums and offsets

    HArray<ValueType> resultValues;

    scai::hmemo::WriteOnlyAccess<ValueType> wResult( resultValues, localN );

    for ( IndexType i = 0; i < localN; i++ )
    {
        wResult[i] = localPrefixSum[i] + myOffset[0];
    }

    return DenseVector<ValueType>( input.getDistributionPtr(), std::move( resultValues ) );
}

/** Generic routine for benchmarking scan routine.
 *
 *  @tparam    ValueType type of values to be scanned
 *  @param[in] N         number of values to scan
 */
template<typename ValueType>
static void bench( const IndexType N )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    PartitionId rank = comm->getRank();

    DistributionPtr blockDist( new BlockDistribution( N, comm ) );

    // generate random numbers

    DenseVector<ValueType> X;

    srand( 131 + comm->getRank() );

    X.setRandom( blockDist, 1 );

    DenseVector<ValueType> Y;  // save for comparison

    bool debug = false;

    common::Settings::getEnvironment( debug, "SCAI_DEBUG" );

    double tmpTime = Walltime::get();

    Y.scan( X );

    double scanTime1 = Walltime::get() - tmpTime;

    HOST_PRINT( rank, "Scan time 1 : " << scanTime1 << " seconds" )

    tmpTime = Walltime::get();

    DenseVector<ValueType> Z = computeGlobalPrefixSum( X );

    double scanTime2 = Walltime::get() - tmpTime;

    HOST_PRINT( rank, "Scan time 2 : " << scanTime2 << " seconds" )

    Y -= Z;
    auto diff = Y.maxNorm();

    HOST_PRINT( rank, "diff between the two solutions: " << diff )
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
        HOST_PRINT( rank, " - n is the size of vector to scan" )
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
        HOST_PRINT( rank, "scan bench for dataType " << dataType << " completed" )
        return 0;
    }
    else
    {
        HOST_PRINT( rank, "scan for dataType " << dataType << " unsupported" )
        return -1;
    }
}
