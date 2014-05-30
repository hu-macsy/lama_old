/**
 * @file Distribution.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @since 1.0.0
 */

// hpp
#include <lama/distribution/Distribution.hpp>

// others
#include <lama/CommunicatorFactory.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// tracing
#include <lama/tracing.hpp>

// boost
#include <boost/scoped_array.hpp>
#include <boost/preprocessor.hpp>

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
    LAMA_LOG_DEBUG( logger, "~Distribution, global size " << mGlobalSize << ", communicator = " << mCommunicator )
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
    IndexType n = requiredIndexes.size();
    owners.resize( n );
    mCommunicator->computeOwners( &owners[0], *this, &requiredIndexes[0], n );
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
    LAMA_REGION( "Distribution.replicate" )

    const Communicator& comm = getCommunicator();
    IndexType currentSize = getLocalSize();
    // Implemenation via cyclic shifting of the vector data and distribution
    IndexType maxLocalSize = comm.max( currentSize );

    LAMA_LOG_INFO( logger, comm << ": replicate localValues<" << Scalar::getType<T2>()
                          << ">[ " << currentSize << ", max = " << maxLocalSize  << " ] "
                          << " to allValues<" << Scalar::getType<T1>() 
                          << ">[ " << getGlobalSize() << " ]" )

    // Only allocate the needed size of the Arrays

    LAMAArray<T1> valuesSend;
    LAMAArray<T1> valuesReceive;
    LAMAArray<IndexType> indexesSend;
    LAMAArray<IndexType> indexesReceive;

    IndexType countValues = 0; // count values set in the global vector

    // set my owned indexes and my values

    ContextPtr commContext = comm.getCommunicationContext();

    // capacity of send arrays should also be sufficient for receiving data

    indexesSend.reserve( commContext, maxLocalSize );
    valuesSend.reserve( commContext, maxLocalSize );

    LAMA_LOG_INFO( logger, "replicate on this communication context: " << *commContext )

    { 
        WriteOnlyAccess<IndexType> wIndexesSend( indexesSend, commContext, currentSize );
        WriteOnlyAccess<T1> wValuesSend( valuesSend, commContext, currentSize );

        // current workaround as commContext works like HostContext

        IndexType* pIndexesSend = wIndexesSend.get();
        T1* pValuesSend = wValuesSend.get();

        for ( IndexType i = 0; i < currentSize; i++ )
        {
            IndexType globalIndex = local2global( i );
            LAMA_ASSERT_DEBUG( globalIndex < getGlobalSize(), *this << ": global index " << globalIndex << " illegal" )
            pIndexesSend[i] = globalIndex;
            pValuesSend[i] = static_cast<T1>( localValues[i] ); // type conversion here
            allValues[globalIndex] = static_cast<T1>( localValues[i] ); // type conversion here
        }
    }

    // capacity of receive arrays should be sufficient for receiving data to avoid reallocations

    indexesReceive.reserve( commContext, maxLocalSize );
    valuesReceive.reserve( commContext, maxLocalSize );

    countValues += currentSize;
    // now nproc - 1 steps for cyclic shifting
    PartitionId np = comm.getSize(); // number partitions

    for ( PartitionId ip = 0; ip < np - 1; ++ip )
    {
        comm.shiftArray( indexesReceive, indexesSend, 1);
        comm.shiftArray( valuesReceive, valuesSend, 1);

        LAMA_ASSERT_EQUAL_DEBUG( valuesReceive.size(), indexesReceive.size() )

        currentSize = valuesReceive.size();

        // sort in the received values

        { 
            ReadAccess<IndexType> readIndexesReceive( indexesReceive, commContext );
            ReadAccess<T1> readValuesReceive( valuesReceive, commContext );

            // current workaround as commContext works like HostContext

            const IndexType* rIndexesReceive = readIndexesReceive.get();
            const T1* rValuesReceive = readValuesReceive.get();
        
            for ( IndexType i = 0; i < currentSize; i++ )
            {
                const IndexType globalIndex = rIndexesReceive[i];
                LAMA_ASSERT_DEBUG( globalIndex < getGlobalSize(), *this << ": global index " << globalIndex << " illegal" )
                allValues[globalIndex] = rValuesReceive[i]; // implicit type conversion done here
            }
        }

        countValues += currentSize;
        indexesReceive.swap( indexesSend );
        valuesReceive.swap( valuesSend );
    }

    LAMA_LOG_INFO( logger, "replicated by " << ( np - 1 ) << " array shifts" )

    // # globalSize values must have been counted
    LAMA_ASSERT_EQUAL_DEBUG( countValues, getGlobalSize() )
}

/* ---------------------------------------------------------------------- */

template<typename T1,typename T2>
void Distribution::replicateN( T1* allValues, const T2* localValues, const IndexType n ) const
{
    LAMA_REGION( "Distribution.replicateN" )

    const Communicator& comm = getCommunicator();

    // Implemenation via cyclic shifting of the vector data and distribution

    // maximal number of elements needed to allocate sufficient receive buffer but also to avoid reallocations

    IndexType currentSize = getLocalSize();
    IndexType maxLocalSize = comm.max( currentSize );

    LAMA_LOG_INFO( logger, comm << ": replicateN, n = " << n << ", localValues<" << Scalar::getType<T2>()
                          << ">[ " << currentSize << ", max = " << maxLocalSize  << " ] "
                          << " to allValues<" << Scalar::getType<T1>() 
                          << ">[ " << getGlobalSize() << " ]" )

    // Only allocate the needed size of the Arrays

    LAMAArray<T1> valuesSend;
    LAMAArray<T1> valuesReceive;
    LAMAArray<IndexType> indexesSend;
    LAMAArray<IndexType> indexesReceive;

    // set my owned indexes and my values

    ContextPtr commContext = comm.getCommunicationContext();

    // capacity of send arrays should also be sufficient for receiving data

    indexesSend.reserve( commContext, maxLocalSize );
    valuesSend.reserve( commContext, maxLocalSize * n );

    LAMA_LOG_INFO( logger, "replicate on this communication context: " << *commContext )

    {
        WriteOnlyAccess<IndexType> wIndexesSend( indexesSend, commContext, currentSize );
        WriteOnlyAccess<T1> wValuesSend( valuesSend, commContext, currentSize * n );

        // current workaround as commContext works like HostContext

        IndexType* pIndexesSend = wIndexesSend.get();
        T1* pValuesSend = wValuesSend.get();

        for ( IndexType i = 0; i < currentSize; i++ )
        {
            IndexType globalIndex = local2global( i );
            LAMA_ASSERT_DEBUG( globalIndex < getGlobalSize(), *this << ": global index " << globalIndex << " illegal" )
            pIndexesSend[i] = globalIndex;
            for ( IndexType j = 0; j < n; j++ )
            {
                pValuesSend[i * n + j] = static_cast<T1>( localValues[i * n + j] );         // type conversion here
                allValues[globalIndex * n + j] = static_cast<T1>( localValues[i * n + j] ); // type conversion here
            }
            pValuesSend[i] = static_cast<T1>( localValues[i] ); // type conversion here
            allValues[globalIndex] = static_cast<T1>( localValues[i] ); // type conversion here
        }
    }

    IndexType countLines = currentSize;   // count values set in the global vector, currentSize has been done

    // capacity of receive arrays should be sufficient for receiving data to avoid reallocations

    indexesReceive.reserve( commContext, maxLocalSize );
    valuesReceive.reserve( commContext, maxLocalSize * n );

    // now nproc - 1 steps for cyclic shifting

    PartitionId np = comm.getSize(); // number partitions

    for ( PartitionId ip = 0; ip < np - 1; ++ip )
    {
        comm.shiftArray( indexesReceive, indexesSend, 1 );
        currentSize = indexesReceive.size();               // next numbe of lines to deal with

        comm.shiftArray( valuesReceive, valuesSend, 1 );

        LAMA_ASSERT_EQUAL_DEBUG( currentSize * n, valuesReceive.size() )

        // sort in the received values

        {
            ReadAccess<IndexType> readIndexesReceive( indexesReceive, commContext );
            ReadAccess<T1> readValuesReceive( valuesReceive, commContext );

            // current workaround as commContext works like HostContext

            const IndexType* rIndexesReceive = readIndexesReceive.get();
            const T1* rValuesReceive = readValuesReceive.get();

            for ( IndexType i = 0; i < currentSize; i++ )
            {
                const IndexType globalIndex = rIndexesReceive[i];
                LAMA_ASSERT_DEBUG( globalIndex < getGlobalSize(), *this << ": global index " << globalIndex << " illegal" )

                for ( IndexType j = 0; j < n; j++ )
                {
                    allValues[globalIndex * n + j] = rValuesReceive[i * n + j];
                }
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
void Distribution::replicateRagged( T allValues[], const T localValues[], const IndexType allOffsets[] ) const
{
    IndexType currentElemSize = getLocalSize();
    const Communicator& comm = getCommunicator();

    // we need the offsets for allValues to store directly the values
    IndexType maxLocalSize = comm.max( currentElemSize );

    LAMAArray<IndexType> indexesSend;
    LAMAArray<IndexType> indexesReceive;

    // boost::scoped_array<IndexType> indexesSend( new IndexType[maxLocalSize] );
    // boost::scoped_array<IndexType> indexesReceive( new IndexType[maxLocalElemSize] );

    ContextPtr commContext = comm.getCommunicationContext();

    indexesReceive.reserve( commContext, maxLocalSize );
    indexesSend.reserve( commContext, maxLocalSize );

    IndexType currentDataSize = 0;

    {
        // due to currentElemSize <= maxLocalSize

        WriteOnlyAccess<IndexType> wIndexesSend( indexesSend, commContext, currentElemSize );

        IndexType* pIndexesSend = wIndexesSend.get();  // is a HostContext

        // pack my indexes, work as global indexes for allValues

        for ( IndexType i = 0; i < currentElemSize; i++ )
        {
            pIndexesSend[i] = local2global( i );
        }

        // fill my local values in all values
        currentDataSize = fillGlobal( allValues, allOffsets, pIndexesSend, currentElemSize, localValues );
    }

    LAMA_LOG_DEBUG( logger,
                    comm << ": filled my local data: " << currentElemSize << " buckets with " << currentDataSize << " values" )
    // get maximal size of data values and allocate send buffer
    IndexType maxLocalDataSize = comm.max( currentDataSize );
    // Implemenation via cyclic shifting of the vector data and distribution
    LAMA_LOG_DEBUG( logger, "maximal data size for exchange = " << maxLocalDataSize )

    LAMA_LOG_INFO( logger, comm << ": replicateRagged<" << Scalar::getType<T>()
                          << ">, localValues[ " << currentDataSize << ", max = " << maxLocalDataSize  << " ] "
                          << " to allValues [ allOffsets[ " << getGlobalSize() << " ] ]" )

    LAMAArray<T> valuesSend;
    LAMAArray<T> valuesReceive;

    valuesSend.reserve( commContext, maxLocalDataSize );
    valuesReceive.reserve( commContext, maxLocalDataSize );

    {
        WriteOnlyAccess<T> wValuesSend( valuesSend, commContext, currentDataSize );
 
        T* pValuesSend = wValuesSend.get();   

        // fill my local values in send buffer

        for ( IndexType i = 0; i < currentDataSize; i++ )
        {
            pValuesSend[i] = localValues[i];
        }
    }

    // Define counters for exchanged values
    IndexType countElemValues = currentElemSize;
    IndexType countDataValues = currentDataSize;
    // now nproc - 1 steps for cyclic shifting
    PartitionId np = comm.getSize(); // number partitions

    for ( PartitionId ip = 0; ip < np - 1; ++ip )
    {
        comm.shiftArray( indexesReceive, indexesSend, 1);
        comm.shiftArray( valuesReceive, valuesSend, 1);

        IndexType newSize1 = indexesReceive.size();
        IndexType newSize2 = valuesReceive.size();

        // copy the received values in allValues at the right place

        {
            ReadAccess<IndexType> rIndexesReceive( indexesReceive, commContext );
            ReadAccess<T> rValuesReceive( valuesReceive, commContext );
 
            IndexType size = -1;
            size = fillGlobal( allValues, allOffsets, rIndexesReceive.get(), newSize1, rValuesReceive.get() );
            LAMA_LOG_DEBUG( logger,
                            comm << ": filled received data: " << newSize1 << " buckets with " << size << " values" )
            LAMA_ASSERT_EQUAL_DEBUG( size, newSize2 )
        }

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

// Macro to instantiate for type pair ARRAY_TYPE##I, ARRAY_TYPE##J

#define LAMA_DISTRIBUTE2_INSTANTIATE(z, J, TYPE)                           \
                                                                           \
template LAMA_DLL_IMPORTEXPORT void Distribution::replicate(               \
    TYPE allValues[],                                                      \
    const ARRAY_TYPE##J localValues[] ) const;                             \


#define LAMA_DISTRIBUTE_INSTANTIATE(z, I, _)                               \
template LAMA_DLL_IMPORTEXPORT void Distribution::replicateRagged(         \
    ARRAY_TYPE##I allValues[],                                             \
    const ARRAY_TYPE##I localValues[],                                     \
    const IndexType allOffsets[] ) const;                                  \
                                                                           \
template LAMA_DLL_IMPORTEXPORT void Distribution::replicateN(              \
    ARRAY_TYPE##I allValues[],                                             \
    const ARRAY_TYPE##I localValues[],                                     \
    const IndexType n ) const;                                             \
                                                                           \
    BOOST_PP_REPEAT( ARRAY_TYPE_CNT,                                       \
                     LAMA_DISTRIBUTE2_INSTANTIATE,                         \
                     ARRAY_TYPE##I )                                       \

// template instantiation for the supported data types

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_DISTRIBUTE_INSTANTIATE, _ )

#undef LAMA_DISTRIBUTE_INSTANTIATE

/* ---------------------------------------------------------------------- */

}
