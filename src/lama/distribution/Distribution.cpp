/**
 * @file Distribution.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Implementation of methods for base class Distribution.
 * @author Jiri Kraus
 * @date 22.02.2011
 * $Id$
 */

// hpp
#include <lama/distribution/Distribution.hpp>

// others
#include <lama/CommunicatorFactory.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// boost
#include <boost/scoped_array.hpp>

namespace lama
{

/* ------  Static class variables --------------------------------------- */

LAMA_LOG_DEF_LOGGER( Distribution::logger, "Distribution" )

/* ------  Constructor  ------------------------------------------------- */

Distribution::Distribution( const IndexType globalSize )
    : mGlobalSize( globalSize ), mCommunicator( CommunicatorFactory::get( "none" ) )
{
    LAMA_LOG_INFO( logger, "Distribution(" << mGlobalSize << ") onto NoCommunicator" )
}

/* ------  Constructor  ------------------------------------------------- */

Distribution::Distribution( const IndexType globalSize, const CommunicatorPtr communicator )
    : mGlobalSize( globalSize ), mCommunicator( communicator )
{
    if ( !mCommunicator )
    {
        LAMA_THROWEXCEPTION( "Distribution without a Communicator is not allowed" )
    }

    LAMA_LOG_INFO( logger, "Distribution( size = " << globalSize << ", comm = " << *mCommunicator << " )" )
}

/* -------  Destructor  ------------------------------------------------- */

Distribution::~Distribution()
{
    LAMA_LOG_INFO( logger, "~Distribution, global size " << mGlobalSize << ", communicator = " << mCommunicator )
}

/* ---------------------------------------------------------------------- */

bool Distribution::operator==( const Distribution& other ) const
{
    // two distributions are the same if they are both replicated
    LAMA_LOG_TRACE( logger, "check " << *this << " == " << other )
    bool isSame = false;

    if ( this == &other )
    {
        isSame = true;
        LAMA_LOG_DEBUG( logger, *this << " == " << other << ": pointer equal" )
    }
    else if ( isReplicated() && other.isReplicated() )
    {
        isSame = getGlobalSize() == other.getGlobalSize();
        LAMA_LOG_DEBUG( logger, *this << " == " << other << ": both are replicated, same size" )
    }
    else
    {
        isSame = isEqual( other );
        LAMA_LOG_DEBUG( logger, *this << " == " << other << ": " << isSame )
    }

    return isSame;
}

/* ---------------------------------------------------------------------- */

bool Distribution::operator!=( const Distribution& other ) const
{
    return !( *this == other );
}

/* ---------------------------------------------------------------------- */

const Communicator& Distribution::getCommunicator() const
{
    return *mCommunicator;
}

/* ---------------------------------------------------------------------- */

CommunicatorPtr Distribution::getCommunicatorPtr() const
{
    return mCommunicator;
}

/* ---------------------------------------------------------------------- */

PartitionId Distribution::getNumPartitions() const
{
    // mCommunicator is never NULL, but just in case
    LAMA_ASSERT( mCommunicator, "Distribution without a Communicator is not allowed" )
    return mCommunicator->getSize();
}

/* ---------------------------------------------------------------------- */

void Distribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "Distribution";
}

/* ---------------------------------------------------------------------- */

void Distribution::computeOwners(
    const std::vector<IndexType>& requiredIndexes,
    std::vector<PartitionId>& owners ) const
{
    LAMA_LOG_INFO( logger, "compute owners via communicator (default)" )
    // use communicator to compute ownership on each processor
    mCommunicator->computeOwners( requiredIndexes, *this, owners );
}

/* ---------------------------------------------------------------------- */

std::ostream& operator<<( std::ostream& stream, Distribution const& dist )
{
    dist.writeAt( stream );
    return stream;
}

/* ---------------------------------------------------------------------- */

template<typename T1,typename T2>
void Distribution::replicate( T1* allValues, const T2* localValues ) const
{
    const Communicator& comm = getCommunicator();
    IndexType currentSize = getLocalSize();
    // Implemenation via cyclic shifting of the vector data and distribution
    IndexType maxLocalSize = comm.max( currentSize );
    LAMA_LOG_DEBUG( logger, comm << ": local size = " << currentSize << ", maximal local size  = " << maxLocalSize )
    // Only allocate the needed size of the Arrays
    boost::scoped_array<T1> valuesSend( new T1[maxLocalSize] );
    boost::scoped_array<T1> valuesReceive( new T1[maxLocalSize] );
    boost::scoped_array<IndexType> indexesSend( new IndexType[maxLocalSize] );
    boost::scoped_array<IndexType> indexesReceive( new IndexType[maxLocalSize] );
    IndexType countValues = 0; // count values set in the global vector

    // set my owned indexes and my values

    for ( IndexType i = 0; i < currentSize; i++ )
    {
        IndexType globalIndex = local2global( i );
        indexesSend[i] = globalIndex;
        valuesSend[i] = static_cast<T1>( localValues[i] ); // type conversion here
        allValues[globalIndex] = static_cast<T1>( localValues[i] ); // type conversion here
    }

    countValues += currentSize;
    // now nproc - 1 steps for cyclic shifting
    PartitionId np = comm.getSize(); // number partitions

    for ( PartitionId ip = 0; ip < np - 1; ++ip )
    {
        IndexType newSize1 = 1;
        newSize1 = comm.shift( indexesReceive.get(), maxLocalSize, indexesSend.get(), currentSize, 1 );
        IndexType newSize2 = -1;
        newSize2 = comm.shift( valuesReceive.get(), maxLocalSize, valuesSend.get(), currentSize, 1 );
        LAMA_ASSERT_EQUAL_DEBUG( newSize1, newSize2 )
        currentSize = newSize1;

        // sort in the received values

        for ( IndexType i = 0; i < currentSize; i++ )
        {
            const IndexType globalIndex = indexesReceive[i];
            allValues[globalIndex] = valuesReceive[i]; // implicit type conversion done here
        }

        countValues += currentSize;
        indexesReceive.swap( indexesSend );
        valuesReceive.swap( valuesSend );
    }

    // # globalSize values must have been counted
    LAMA_ASSERT_EQUAL_DEBUG( countValues, getGlobalSize() )
}

/* ---------------------------------------------------------------------- */

template<typename T1,typename T2>
void Distribution::replicateN( T1* allValues, const T2* localValues, const IndexType n ) const
{
    const Communicator& comm = getCommunicator();
    IndexType currentSize = getLocalSize();
    // Implemenation via cyclic shifting of the vector data and distribution
    IndexType maxLocalSize = comm.max( currentSize );
    LAMA_LOG_DEBUG( logger, comm << ": local size = " << currentSize << ", maximal local size  = " << maxLocalSize )
    // Only allocate the needed size of the Arrays
    boost::scoped_array<T1> valuesSend( new T1[maxLocalSize * n] );
    boost::scoped_array<T1> valuesReceive( new T1[maxLocalSize * n] );
    boost::scoped_array<IndexType> indexesSend( new IndexType[maxLocalSize] );
    boost::scoped_array<IndexType> indexesReceive( new IndexType[maxLocalSize] );
    IndexType countLines = 0; // count values set in the global vector

    // set my owned indexes and my values

    for ( IndexType i = 0; i < currentSize; i++ )
    {
        IndexType globalIndex = local2global( i );
        indexesSend[i] = globalIndex;

        for ( IndexType j = 0; j < n; j++ )
        {
            valuesSend[i * n + j] = localValues[i * n + j]; // type conversion here
            allValues[globalIndex * n + j] = localValues[i * n + j]; // type conversion here
        }
    }

    countLines += currentSize;
    // now nproc - 1 steps for cyclic shifting
    PartitionId np = comm.getSize(); // number partitions

    for ( PartitionId ip = 0; ip < np - 1; ++ip )
    {
        IndexType newSize1 = -1;
        newSize1 = comm.shift( indexesReceive.get(), maxLocalSize, indexesSend.get(), currentSize, 1 );
        IndexType newSize2 = -1;
        newSize2 = comm.shift( valuesReceive.get(), n * maxLocalSize, valuesSend.get(), n * currentSize, 1 );
        LAMA_ASSERT_EQUAL_DEBUG( newSize1 * n, newSize2 )
        currentSize = newSize1;

        // sort in the received values

        for ( IndexType i = 0; i < currentSize; i++ )
        {
            const IndexType globalIndex = indexesReceive[i];

            for ( IndexType j = 0; j < n; j++ )
            {
                allValues[globalIndex * n + j] = valuesReceive[i * n + j];
            }
        }

        countLines += currentSize;
        indexesReceive.swap( indexesSend );
        valuesReceive.swap( valuesSend );
    }

    // # globalSize values must have been counted
    LAMA_ASSERT_EQUAL_DEBUG( countLines, getGlobalSize() )
}

/* ---------------------------------------------------------------------- */

template<typename T>
static IndexType fillGlobal(
    T* allValues,
    const IndexType* allOffsets,
    const IndexType* indexes,
    const IndexType numIndexes,
    const T* values )
{
    IndexType counter = 0; // traverses values, counts the number of filled values

    for ( IndexType i = 0; i < numIndexes; i++ )
    {
        const IndexType allIndex = indexes[i];
        const IndexType offset = allOffsets[allIndex];
        const IndexType size = allOffsets[allIndex + 1] - offset;

        /* No logger available here:

         std::cout << "fill index " << i << " of " << numIndexes << ": global = "
         << allIndex << ", size = " << size << ", offset = " << offset  << std::endl;
         */

        for ( IndexType j = 0; j < size; j++ )
        {
            allValues[offset + j] = values[counter + j];
        }

        counter += size;
    }

    return counter; // number of filled elements
}

/* ---------------------------------------------------------------------- */

template<typename T>
void Distribution::replicateRagged( T* allValues, const T* localValues, const IndexType* allOffsets ) const
{
    IndexType currentElemSize = getLocalSize();
    const Communicator& comm = getCommunicator();
    IndexType maxLocalElemSize = comm.max( currentElemSize );
    // we need the offsets for allValues to store directly the values
    IndexType maxLocalSize = comm.max( currentElemSize );
    boost::scoped_array<IndexType> indexesSend( new IndexType[maxLocalSize] );
    boost::scoped_array<IndexType> indexesReceive( new IndexType[maxLocalElemSize] );

    // pack my indexes, work as global indexes for allValues

    for ( IndexType i = 0; i < currentElemSize; i++ )
    {
        indexesSend[i] = local2global( i );
    }

    // fill my local values in all values
    IndexType currentDataSize = fillGlobal( allValues, allOffsets, indexesSend.get(), currentElemSize, localValues );
    LAMA_LOG_DEBUG( logger,
                    comm << ": filled my local data: " << currentElemSize << " buckets with " << currentDataSize << " values" )
    // get maximal size of data values and allocate send buffer
    IndexType maxLocalDataSize = comm.max( currentDataSize );
    // Implemenation via cyclic shifting of the vector data and distribution
    LAMA_LOG_DEBUG( logger, "maximal data size for exchange = " << maxLocalDataSize )
    boost::scoped_array<T> valuesSend( new T[maxLocalDataSize] );
    boost::scoped_array<T> valuesReceive( new T[maxLocalDataSize] );

    // fill my local values in send buffer

    for ( IndexType i = 0; i < currentDataSize; i++ )
    {
        valuesSend[i] = localValues[i];
    }

    // Define counters for exchanged values
    IndexType countElemValues = currentElemSize;
    IndexType countDataValues = currentDataSize;
    // now nproc - 1 steps for cyclic shifting
    PartitionId np = comm.getSize(); // number partitions

    for ( PartitionId ip = 0; ip < np - 1; ++ip )
    {
        int direction = 1;
        IndexType newSize1 = comm.shift( indexesReceive.get(), maxLocalSize, indexesSend.get(), currentElemSize,
                                         direction );
        IndexType newSize2 = comm.shift( valuesReceive.get(), maxLocalDataSize, valuesSend.get(), currentDataSize,
                                         direction );
        // fill the received values in allValues
        IndexType size = -1;
        size = fillGlobal( allValues, allOffsets, indexesReceive.get(), newSize1, valuesReceive.get() );
        LAMA_LOG_DEBUG( logger,
                        comm << ": filled received data: " << newSize1 << " buckets with " << size << " values" )
        LAMA_ASSERT_EQUAL_DEBUG( size, newSize2 )
        currentElemSize = newSize1;
        currentDataSize = newSize2;
        countElemValues += currentElemSize;
        countDataValues += currentDataSize;
        // swap the send and receive buffer
        indexesReceive.swap( indexesSend );
        valuesReceive.swap( valuesSend );
    }

    // verify that all values are available
    LAMA_ASSERT_EQUAL_DEBUG( countElemValues, getGlobalSize() )
    LAMA_ASSERT_EQUAL_DEBUG( countDataValues, allOffsets[getGlobalSize()] )
}

/* ---------------------------------------------------------------------- */

// Instantiation of all relevant replicate routines
template LAMA_DLL_IMPORTEXPORT void Distribution::replicate( double* allValues, const double* localValues ) const;
template LAMA_DLL_IMPORTEXPORT void Distribution::replicate( double* allValues, const float* localValues ) const;
template LAMA_DLL_IMPORTEXPORT void Distribution::replicate( float* allValues, const double* localValues ) const;
template LAMA_DLL_IMPORTEXPORT void Distribution::replicate( float* allValues, const float* localValues ) const;
template LAMA_DLL_IMPORTEXPORT void Distribution::replicate( IndexType* allValues, const IndexType* localValues ) const;

template LAMA_DLL_IMPORTEXPORT void Distribution::replicateN(
    double* allValues,
    const double* localValues,
    const IndexType n ) const;
template LAMA_DLL_IMPORTEXPORT void Distribution::replicateN(
    float* allValues,
    const float* localValues,
    const IndexType n ) const;

template LAMA_DLL_IMPORTEXPORT void Distribution::replicateRagged(
    double* allValues,
    const double* localValues,
    const IndexType* allOffsets ) const;
template LAMA_DLL_IMPORTEXPORT void Distribution::replicateRagged(
    float* allValues,
    const float* localValues,
    const IndexType* allOffsets ) const;
template LAMA_DLL_IMPORTEXPORT void Distribution::replicateRagged(
    IndexType* allValues,
    const IndexType* localValues,
    const IndexType* allOffsets ) const;

/* ---------------------------------------------------------------------- */

}
